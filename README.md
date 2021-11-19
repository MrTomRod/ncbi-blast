# NcbiBlast

## Idea:

Simple Python wrapper for NCBI Blast.

## Setup:

You need to have ncbi-blast. You can download it using `get_latest_blast.py`

## Usage:

```python3
# load class
from ncbi_blast.Blast import Blast

blast = Blast(blast_path='/local/bin', outfmt=6)

queryfile = open('tests/data/query_single_nucl.fasta').read()
nucl_fasta1 = 'tests/data/ffn/prot_nucl_1.ffn'
nucl_fasta2 = 'tests/data/ffn/prot_nucl_2.ffn'

# create blast database for nucl_fastas
blast.mkblastdb(file=nucl_fasta1, dbtype='nucl')
blast.mkblastdb(file=nucl_fasta2, dbtype='nucl')

# blast database
output = blast.blastn(fasta_string=queryfile, db=nucl_fasta2)
# example output (outfmt=6): 'QUERY_ID\MATCH_ID\t100.000\t1677\t0\t0\t1\t1677\t1\t1677\t0.0\t3097\n'

# blast multiple databases
output = blast.blastn(fasta_string=queryfile, db=[nucl_fasta1, nucl_fasta2])
```