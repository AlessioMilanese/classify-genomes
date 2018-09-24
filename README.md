classify-genomes
========
This tool classify genomes sequences (as well as metagenomic assembled genome) according to the [mOTUs taxonomy](https://github.com/motu-tool/mOTUs_v2).

Pre-requisites
--------------
* Perl 5.24.0
* cdbtools
* Prodigal 2.6.3
* Python 2.7.12
* SAMtools 1.3.1
* HMMER 3.1b2
* vsearch

If you have [conda](https://conda.io/docs/), you can install the dependencies and create an enviroment (after cloning the directory, see Installation):
```bash
cd classify-genomes/env
conda env create -f classify-genomes.yaml
source activate classify-genomes-ENV
```
Note: type `source deactivate` to deactivate an active environment.

Installation
--------------
```bash
git clone https://github.com/AlessioMilanese/classify-genomes.git
cd classify-genomes
python setup.py
```

Note: in the following examples we assume that the python script ```classify-genomes``` is in the system path.


Simple examples
--------------

Download the genome with id [AWWC00000000](https://www.ncbi.nlm.nih.gov/nuccore/AWWC00000000.1) that was assembled in the HMP project and annotated as *Jonquetella sp. BV3C21*
```bash
wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/AW/WC/AWWC01/AWWC01.1.fsa_nt.gz
gunzip AWWC01.1.fsa_nt.gz
```
(if it doesn't work, download JPFJ01.1.fsa_nt.gz from https://www.ncbi.nlm.nih.gov/Traces/wgs/AWWC01?display=contigs)

Find the taxonomy annotation with classify-genomes:
```bash
classify-genomes AWWC01.1.fsa_nt
```

Which results in:
```
Extract genes
Map genes to mOTUs database
Find taxonomy

RESULT:
Fasta sequence: AWWC01.1.fsa_nt
Consensus taxonomy: 428712 Jonquetella anthropi [C]
Number of detected genes: 9
Percentage of agreeing genes: 100.0%

Gene annotation: AWWC01.1.fsa_nt
COG0016_1	428712 Jonquetella anthropi [C]
COG0012_1	428712 Jonquetella anthropi [C]
COG0541_1	428712 Jonquetella anthropi [C]
COG0172_1	428712 Jonquetella anthropi [C]
COG0552_1	428712 Jonquetella anthropi [C]
COG0525_1	428712 Jonquetella anthropi [C]
COG0533_1	428712 Jonquetella anthropi [C]
COG0495_1	428712 Jonquetella anthropi [C]
COG0018_1	428712 Jonquetella anthropi [C]
```

Command options
--------------

The tool expect as input a fasta file:
```
classify-genomes <fasta_file> [options]
```

The options are:
* **`-t`** number of threads (default 1)
* **`-v`** verbose level: 1=error, 2=warning, 3=message, 4+=debugging. Default is 3.
* **`-o`** save result to a file. Default is standard output.
* **`-m`** save the fasta sequence of the extracted marker genes.
* **`-s`** print the result in a single line. Useful if you want to analyse many genomes and cat all the taxonomies.
