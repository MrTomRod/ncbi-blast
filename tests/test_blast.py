from unittest import TestCase
from ncbi_blast.Blast import Blast
import os


class TestMakeBlastDBs(TestCase):
    def setUp(self) -> None:
        self.blast = Blast()

        self.prot_file_extensions = ['.pdb', '.phr', '.pin', '.pot', '.psq', '.ptf', '.pto']
        self.nucl_file_extensions = ['.ndb', '.nhr', '.nin', '.not', '.nsq', '.ntf', '.nto']

        self.test_fna = '../data/fna/assembly.fna'
        self.test_faa = '../data/faa/prot_prot.faa'
        self.test_ffn = '../data/ffn/prot_nucl_1.ffn'

        for file in [self.test_fna, self.test_faa, self.test_ffn]:
            assert os.path.isfile(file)

    def tearDown(self) -> None:
        pass

    def test_fna(self):
        file = self.test_fna
        for ext in self.nucl_file_extensions:
            _file = file + ext
            if os.path.isfile(_file): os.remove(_file)

        self.blast.mkblastdb(file, 'nucl', taxid=None, title=None, overwrite=True)
        print(self.blast.blast_db_info(file))

        for ext in self.nucl_file_extensions:
            _file = file + ext
            self.assertTrue(os.path.isfile(_file))
            os.remove(_file)

    def test_ffn(self):
        file = self.test_ffn
        for ext in self.nucl_file_extensions:
            _file = file + ext
            if os.path.isfile(_file): os.remove(_file)

        self.blast.mkblastdb(file, 'nucl', taxid=None, title=None, overwrite=True)
        print(self.blast.blast_db_info(file))

        for ext in self.nucl_file_extensions:
            _file = file + ext
            self.assertTrue(os.path.isfile(_file))
            os.remove(_file)

    def test_faa(self):
        file = self.test_faa
        for ext in self.prot_file_extensions:
            _file = file + ext
            if os.path.isfile(_file): os.remove(_file)

        self.blast.mkblastdb(file, 'prot', taxid=None, title=None, overwrite=True)
        print(self.blast.blast_db_info(file))

        for ext in self.prot_file_extensions:
            _file = file + ext
            self.assertTrue(os.path.isfile(_file))
            os.remove(_file)

    def test_faa_fail(self):
        file = self.test_faa
        with self.assertRaises(AssertionError):
            self.blast.mkblastdb(file, 'nucl', taxid=None, title=None, overwrite=True)

        for ext in self.nucl_file_extensions:
            _file = file + ext
            os.remove(_file)


class TestBlast(TestCase):
    def setUp(self) -> None:
        self.blast = Blast()

        self.prot_file_extensions = ['.pdb', '.phr', '.pin', '.pot', '.psq', '.ptf', '.pto']
        self.nucl_file_extensions = ['.ndb', '.nhr', '.nin', '.not', '.nsq', '.ntf', '.nto']

        self.test_fna = '../data/fna/assembly.fna'
        self.blast.mkblastdb(self.test_fna, 'nucl', taxid=None, title=None, overwrite=True)

        self.test_faa = '../data/faa/prot_prot.faa'
        self.blast.mkblastdb(self.test_faa, 'prot', taxid=None, title=None, overwrite=True)

        self.test_ffn_1 = '../data/ffn/prot_nucl_1.ffn'
        self.blast.mkblastdb(self.test_ffn_1, 'nucl', taxid=None, title=None, overwrite=True)

        self.test_ffn_2 = '../data/ffn/prot_nucl_2.ffn'
        self.blast.mkblastdb(self.test_ffn_2, 'nucl', taxid=None, title=None, overwrite=True)

        self.single_prot_query = open('../data/query_single_prot.fasta').read()
        self.single_nucl_query = open('../data/query_single_nucl.fasta').read()
        self.multi_prot_query = open('../data/query_multiple_prot.fasta').read()
        self.multi_nucl_query = open('../data/query_multiple_nucl.fasta').read()

        for file in [self.test_fna, self.test_faa, self.test_ffn_1, self.test_ffn_2]:
            assert os.path.isfile(file)

    def tearDown(self) -> None:
        all_extensions = self.nucl_file_extensions
        all_extensions.extend(self.prot_file_extensions)
        for file in [self.test_fna, self.test_faa, self.test_ffn_1]:
            for ext in all_extensions:
                _file = file + ext
                if os.path.isfile(_file): os.remove(_file)

    def test_single_get_alphabet_prot(self):
        self.assertTrue(self.blast.is_protein(self.single_prot_query))
        self.assertFalse(self.blast.is_dna(self.single_prot_query))

    def test_multi_get_alphabet_prot(self):
        self.assertTrue(self.blast.is_protein(self.multi_prot_query))
        self.assertFalse(self.blast.is_dna(self.multi_prot_query))

    def test_single_get_alphabet_nucl(self):
        self.assertTrue(self.blast.is_dna(self.single_nucl_query))
        self.assertTrue(self.blast.is_protein(self.single_nucl_query))
        self.assertFalse(self.blast.is_protein_and_not_dna(self.single_nucl_query))
        self.assertTrue(self.blast.is_protein('>abc\nATCATCEDAGD'))
        self.assertFalse(self.blast.is_protein('>abc\nATCATCXXXX'))

    def test_multi_get_alphabet_nucl(self):
        self.assertTrue(self.blast.is_dna(self.multi_nucl_query))
        self.assertTrue(self.blast.is_protein(self.multi_nucl_query))
        self.assertFalse(self.blast.is_protein_and_not_dna(self.multi_nucl_query))
        self.assertTrue(self.blast.is_protein('>abc\nATCATCEDAGD\nATCATCEDAGD'))
        self.assertFalse(self.blast.is_protein('>abc\nATCATCXXXX\nATCATCXXXX'))

    def test_single_blastp(self):
        out = self.blast.blastp(fasta_string=self.single_prot_query, db=self.test_faa)
        first_match = out.split('\t')[1]
        first_match = first_match.split('|')[-1]
        self.assertEqual('FAM19036p_000425', first_match)

    def test_multi_blastp(self):
        out = self.blast.blastp(fasta_string=self.multi_prot_query, db=self.test_faa)
        first_match = out.split('\t')[1]
        first_match = first_match.split('|')[-1]
        self.assertEqual('FAM19036p_000425', first_match)

    def test_single_blastx(self):
        out = self.blast.blastx(fasta_string=self.single_nucl_query, db=self.test_faa)
        first_match = out.split('\t')[1]
        first_match = first_match.split('|')[-1]
        self.assertEqual('FAM19036p_000425', first_match)

    def test_multi_blastx(self):
        out = self.blast.blastx(fasta_string=self.multi_nucl_query, db=self.test_faa)
        first_match = out.split('\t')[1]
        first_match = first_match.split('|')[-1]
        self.assertEqual('FAM19036p_000425', first_match)

    def test_single_blastn_ffn(self):
        out = self.blast.blastn(fasta_string=self.single_nucl_query, db=self.test_ffn_1)
        first_match = out.split('\t')[1]
        first_match = first_match.split('|')[-1]
        self.assertEqual('FAM19036p_000425', first_match)

    def test_multi_blastn_ffn(self):
        out = self.blast.blastn(fasta_string=self.multi_nucl_query, db=self.test_ffn_1)
        first_match = out.split('\t')[1]
        first_match = first_match.split('|')[-1]
        self.assertEqual('FAM19036p_000425', first_match)

    def test_single_blastn_fna(self):
        out = self.blast.blastn(fasta_string=self.single_nucl_query, db=self.test_fna)
        match_quality = out.split('\t')[2]
        self.assertEqual(match_quality, '100.000')

    def test_multi_blastn_fna(self):
        out = self.blast.blastn(fasta_string=self.multi_nucl_query, db=self.test_fna)
        match_quality = out.split('\t')[2]
        self.assertEqual(match_quality, '100.000')

    def test_kwargs(self):
        out = self.blast.blastp(fasta_string=self.single_prot_query, db=self.test_faa, num_alignments=1)
        self.assertEqual(len(out.strip().split('\n')), 1)


class TestMultipleBlastDBs(TestCase):
    def setUp(self) -> None:
        self.blast = Blast()

        self.prot_file_extensions = ['.pdb', '.phr', '.pin', '.pot', '.psq', '.ptf', '.pto']
        self.nucl_file_extensions = ['.ndb', '.nhr', '.nin', '.not', '.nsq', '.ntf', '.nto']

        self.single_prot_query = open('../data/query_single_prot.fasta').read()
        self.single_nucl_query = open('../data/query_single_nucl.fasta').read()
        self.multi_prot_query = open('../data/query_multiple_prot.fasta').read()
        self.multi_nucl_query = open('../data/query_multiple_nucl.fasta').read()

        self.files = ['prot_nucl_1.ffn', 'prot_nucl_2.ffn']
        self.files = ['../data/ffn/' + file for file in self.files]

        for file in self.files:
            assert os.path.isfile(file) and ' ' not in file, file
            self.blast.mkblastdb(file, 'nucl', taxid=None, title=None, overwrite=True)

    def tearDown(self) -> None:
        for file in self.files:
            for ext in self.nucl_file_extensions:
                _file = file + ext
                if os.path.isfile(_file): os.remove(_file)

    def test_multiple_ffn(self):
        out = self.blast.blastn(self.single_nucl_query, self.files).strip()
        out = out.split('\n')
        self.assertEqual(2, len(out))
        line0 = out[0].split('\t')
        line1 = out[1].split('\t')
        matches = {line0[1], line1[1]}
        self.assertEqual({'FAM19036p_000425', 'FAKEENTRY'}, matches)

    def test_parse_kwarg_string(self):
        for kwarg_string in ['', '-evalue 0.01', '-evalue 1e-1', '-matrix BLOSUM80']:
            for delim in [' ', '=']:
                print(self.blast.parse_kwarg_string(kwarg_string.replace(' ', delim)))

        for kwarg_string in ['-arg /etc', '-arg "ab cdef"', '-arg $XXX', 'arg 0.4', 'ARG 0.4']:
            with self.assertRaises(AssertionError):
                print(self.blast.parse_kwarg_string(kwarg_string))

    def test_kwargs_as_list(self):
        self.assertEqual(self.blast.kwargs_as_list({'-a': '0', '-b': '1', '-c': '2'}),
                         ['-a', '0', '-b', '1', '-c', '2'])
        self.assertEqual(self.blast.kwargs_as_list({}), [])
