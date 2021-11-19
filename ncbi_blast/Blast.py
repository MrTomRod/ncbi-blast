import os
import re
import html
from io import StringIO
from subprocess import run, PIPE
from tempfile import NamedTemporaryFile

from Bio import SeqIO


def is_installed(program):
    """
    Test if a program is installed.

    :param program: path to executable or command
    :return: if program executable: program; else None
    """

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:  # check if path to program is valid
        return is_exe(program)
    else:  # check if program is in PATH
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
        return False


class Blast:
    def __init__(self, blast_path: str = None, outfmt: int = 6, blast_columns: [str] = None, verbose: bool = True):
        assert outfmt in range(0, 19), 'outfmt must be between 0 and 18'
        if blast_columns is not None:
            assert outfmt in {6, 7, 10}, f'blast_columns can only be specified with outfmt 6, 7 and 10'

        self.verbose = verbose
        self._outfmt = outfmt
        self._blast_columns = blast_columns
        self.executables = dict(
            makeblastdb='makeblastdb',
            blastn='blastn',
            blastp='blastp',
            blastx='blastx',
            tblastn='tblastn',
            tblastx='tblastx'
        )

        if blast_path is not None:
            self.executables = {tool: os.path.join(blast_path, tool) for tool in self.executables}

        for tool, executable in self.executables.items():
            assert is_installed(executable), F'{tool} is not executable ({executable})'

    def __str__(self) -> str:
        return f'Blast:outfmt={self.outfmt};blast_columns={self._blast_columns}'

    @property
    def outfmt(self) -> str:
        if self._blast_columns is None:
            return str(self._outfmt)
        else:
            return f'"{self._outfmt} {" ".join(self._blast_columns)}"'

    def version(self):
        command = [self.executables["blastp"], '-version']
        subprocess = run(command, stdout=PIPE, stderr=PIPE, encoding='ascii')
        assert subprocess.returncode == 0, self.error_message('blastp', command, subprocess)
        return subprocess.stdout.split('\n')[0].split(' ')[1]

    def mkblastdb(self, file, dbtype, taxid=None, title=None, overwrite=True):
        assert os.path.isfile(file), F'file does not exist: {file}'
        assert ' ' not in file, F'file paths may not contain blanks: {file}'
        assert dbtype in ['prot', 'nucl'], 'dbtype must be either prot or nucl'
        if dbtype == 'prot': file_extensions = ['.pdb', '.phr', '.pin', '.pot', '.psq', '.ptf', '.pto']
        if dbtype == 'nucl': file_extensions = ['.ndb', '.nhr', '.nin', '.not', '.nsq', '.ntf', '.nto']

        if not overwrite:
            n_files_missing = sum([not os.path.isfile(file + ext) for ext in file_extensions])
            if n_files_missing == 0: return

        command = [self.executables["makeblastdb"], "-in", file, "-dbtype", dbtype]

        if title:
            command.extend(['-title', title])
        if taxid:
            command.extend(['-taxid', str(taxid)])
        subprocess = run(command, stdout=PIPE, stderr=PIPE, encoding='ascii')

        error_message = self.error_message('makeblastdb', ' '.join(command), subprocess)
        assert subprocess.stderr == '', error_message
        assert subprocess.returncode == 0, error_message

    def blastp(self, fasta_string, db, **kwargs):
        if not self.is_protein(fasta_string):
            raise TypeError('fasta_string is not a valid protein sequence')
        return self.blast(fasta_string, db, mode='blastp', **kwargs)

    def blastn(self, fasta_string, db, **kwargs):
        if not self.is_dna(fasta_string):
            raise TypeError('fasta_string is not a valid DNA sequence')
        return self.blast(fasta_string, db, mode='blastn', **kwargs)

    def blastx(self, fasta_string, db, **kwargs):
        if not self.is_dna(fasta_string):
            raise TypeError('fasta_string is not a valid DNA sequence')

        return self.blast(fasta_string, db, mode='blastx', **kwargs)

    def tblastn(self, fasta_string, db, **kwargs):
        if not self.is_protein(fasta_string):
            raise TypeError('fasta_string is not a valid protein sequence')
        return self.blast(fasta_string, db, mode='tblastn', **kwargs)

    def blast(self, fasta_string, db, mode='blastp', **kwargs):
        assert mode in ['blastp', 'blastn', 'blastx', 'tblastn'], F"mode must be 'blastp', 'blastn', 'blastx' or 'tblastn', is: {mode}"
        assert type(db) in [str, list], F'db must be a string (path to db) or a list of strings (multiple db paths)'

        # blast database
        if type(db) is str:
            db = [db]

        for file in db:
            assert os.path.isfile(file), F'file does not exist: {file}'
            assert ' ' not in file, F'file paths may not contain blanks: {file}'
        db = F'"{" ".join(db)}"'

        # kwargs
        kwargs = self.kwargs_as_list(self.clean_kwargs({f'-{k}': str(a) for k, a in kwargs.items()}))

        # write fasta to temporary file
        temp_file = NamedTemporaryFile(mode='w', delete=True)
        with temp_file as f:
            f.write(fasta_string)
            f.flush()

            # execute command
            command = [self.executables[mode], '-query', temp_file.name, '-db', db, '-outfmt', self.outfmt] + kwargs
            subprocess = run(' '.join(command), stdout=PIPE, stderr=PIPE, encoding='ascii', shell=True)

        temp_file.close()  # delete temporary file

        assert subprocess.returncode == 0, self.error_message(mode, ' '.join(command), subprocess)

        return subprocess.stdout.rstrip()

    def is_protein_and_not_dna(self, fasta_string):
        if self.is_dna(fasta_string):
            return False
        return self.is_protein(fasta_string)

    def is_protein(self, fasta_string):
        """
        Note: A DNA sequence is also a valid protein sequence.
        """
        parsed_seq = [i for i in SeqIO.parse(StringIO(fasta_string), 'fasta')]
        for seq in parsed_seq:
            if not self.__verify_alphabet(seq.seq, 'ACDEFGHIKLMNPQRSTVWY'):
                return False
        return True

    def is_dna(self, fasta_string):
        parsed_seq = [i for i in SeqIO.parse(StringIO(fasta_string), 'fasta')]
        for seq in parsed_seq:
            if not self.__verify_alphabet(seq.seq, 'GATCRYWSMKHBVDN'):
                return False
        return True

    @staticmethod
    def __verify_alphabet(sequence, letters):
        for letter in sequence:
            if letter not in letters:
                return False
        return True

    @staticmethod
    def parse_kwarg_string(kwarg_string: str) -> dict[str:str]:
        if kwarg_string == '':
            return {}

        kwargs = re.split(r'[ =]', kwarg_string)
        assert len(kwargs) % 2 == 0, f'One argument is missing a value! {kwarg_string}'
        kwargs = dict(zip(
            kwargs[::2],  # even elements (argument)
            kwargs[1::2]  # odd elements (value)
        ))

        return Blast.clean_kwargs(kwargs)

    @staticmethod
    def clean_kwargs(kwargs: dict[str:str]) -> dict[str:str]:
        """
        Ensure the arguments contain no dangerous characters
        """
        # blast keywords start with '-', followed by lower case letters and '_'
        re_valid_kw = re.compile(pattern=r'\-[a-z_]+')
        # blast arguments may contain letters, numbers, '-' and '.'
        re_valid_arg = re.compile(pattern=r'[a-zA-Z0-9\-\.]+')

        for kw, arg in kwargs.items():
            assert bool(re_valid_kw.fullmatch(kw)), f'Invalid parameter: {html.escape(kw)}'
            assert bool(re_valid_arg.fullmatch(arg)), f'Invalid parameter: {html.escape(arg)}'

        return kwargs

    @staticmethod
    def kwargs_as_list(kwargs: dict[str:str]) -> [str]:
        return [kw_or_arg for kw_arg_tuple in kwargs.items() for kw_or_arg in kw_arg_tuple]

    def error_message(self, tool: str, command: [str], subprocess) -> str:
        if self.verbose:
            return f'{tool} command failed: "{" ".join(command)},\n stdout: {subprocess.stdout}",\n stderr: {subprocess.stderr}'
        else:
            return subprocess.stderr
