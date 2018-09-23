classify-genomes
========
This tool classify genomes sequences (as well as metagenomic assembled genome) according to the [mOTUs taxonomy](https://github.com/motu-tool/mOTUs_v2). 

Pre-requisites
--------------
* Perl 5.24.0
* Prodigal 2.6.3
* Python 2.7.12
* SAMtools 1.3.1
* HMMER 3.1b2
* vsearch

If you have [conda](https://conda.io/docs/), you can install the dependencies and create an enviroment (after cloning the directory, see Installation):
```bash
cd classify-genomes/env
conda env create -f classify-genomes.yaml
source activate classify-genomes-env
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
