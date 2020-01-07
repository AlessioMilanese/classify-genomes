classify-genomes
========
This tool classify genomes sequences (as well as metagenomic assembled genome) according to the [mOTUs taxonomy](https://github.com/motu-tool/mOTUs_v2) (mOTUs version 2.0.0).

Pre-requisites
--------------
* Perl 5.24.0
* cdbtools
* Prodigal 2.6.3
* Python 2.7.12
* SAMtools 1.3.1
* HMMER 3.1b2
* vsearch

If you have [conda](https://conda.io/docs/), you can install the dependencies and create an environment (after cloning the directory, see Installation):
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
Consensus NCBI ID: 428712
Consensus taxonomy: Jonquetella anthropi [Jonquetella sp. BV3C21/Jonquetella anthropi]
Consensus mOTUs: ref_mOTU_v25_01951

Annotation: Consistent
Confidence: 99.9 %

Number of detected genes: 10
Number of mapped genes: 9
Number of agreeing genes: 9
Percentage of agreeing genes: 100.0%

Gene annotation:
COG0016_1	100.0	428712 Jonquetella anthropi [Jonquetella sp. BV3C21/Jonquetella anthropi]@ref_mOTU_v25_01951
COG0012_1	100.0	428712 Jonquetella anthropi [Jonquetella sp. BV3C21/Jonquetella anthropi]@ref_mOTU_v25_01951
COG0541_1	99.9	428712 Jonquetella anthropi [Jonquetella sp. BV3C21/Jonquetella anthropi]@ref_mOTU_v25_01951
COG0525_1	100.0	428712 Jonquetella anthropi [Jonquetella sp. BV3C21/Jonquetella anthropi]@ref_mOTU_v25_01951
COG0552_1	100.0	428712 Jonquetella anthropi [Jonquetella sp. BV3C21/Jonquetella anthropi]@ref_mOTU_v25_01951
COG0172_1	100.0	428712 Jonquetella anthropi [Jonquetella sp. BV3C21/Jonquetella anthropi]@ref_mOTU_v25_01951
COG0533_1	100.0	428712 Jonquetella anthropi [Jonquetella sp. BV3C21/Jonquetella anthropi]@ref_mOTU_v25_01951
COG0495_1	99.9	428712 Jonquetella anthropi [Jonquetella sp. BV3C21/Jonquetella anthropi]@ref_mOTU_v25_01951
COG0018_1	100.0	428712 Jonquetella anthropi [Jonquetella sp. BV3C21/Jonquetella anthropi]@ref_mOTU_v25_01951
```

Let's analyse the result: the genome in the file `AWWC01.1.fsa_nt` is annotated as `Jonquetella anthropi` with a NCBI taxonomy id equal to `428712`, which belongs to the mOTU `ref_mOTU_v25_01951`.

The genome is `Consistent` in its annotation (based on the marker genes), and we provide a confidence on the annotation (see more info on [Estimate correctness of the assignment](https://github.com/AlessioMilanese/classify-genomes#estimate-correctness-of-the-assignment)).

The genome contains `10` marker genes (MGs) out of 10. The marker genes have the property to be present in single-copy in bacterial genomes, hence if the tool extract more than 10 MGs there might be problems with the genome that you are analysing. After that there is the information of the number of genes that were mapped to the mOTUs database (`Number of mapped genes`) and the number of genes that support the consensus taxonomy `Number of agreeing genes` (as well as how much percentage of the found ones they represents, in this case `100%`). In the case above, `10` genes were detected in the genome and `9` mapped to the mOTUs database. All the `9` mapped genes agrees on the same annotation (`Number of agreeing genes`). Check the following example:

![alt text](https://github.com/AlessioMilanese/classify-genomes/blob/master/pics/explain_gene_numbers.png)

Finally, there is a list with all the identified genes, the target gene in the mOTUs database and the percentage identity.


Implementation details
--------------

The classification of a fasta genome is divided into four steps:

1. **Predict genes**, we use Prodigal (Hyatt et al. BMC Bioinf 2010) to predict all genes and proteins from the input fasta file.
2. **Extract Marker genes**, we extract the 10 marker genes used to build the mOTUs (Milanese et al. Nature Comm 2019) using fetchMG.
3. **Map marker genes to the mOTUs database**, we map the extracted marker genes to the marker genes in the 14,212 mOTUs. We use vsearch (Rognes et al. PeerJ 2016) with the command `--usearch_global` and option `--maxrejects 10000 --minqt 0.7`. We use marker gene specific thresholds that were calibrated in Milanese et al. Nature Comm. 2019 to distinguish between different species. For example, `COG0049` is a conserved gene with threshold `98.9`, while `COG0124` is evolutionary under less pressure having a threshold of `94.7`.
4. **Find taxonomy annotation for the genome**. From the annotation of the single genes we derive the annotation for the submitted genome, which can be:
* No genes, no marker genes have been detected in the fasta file;
* No match, marker genes have been detected in the input genome, but there is no match to the mOTUs;
* Inconsistent, marker genes have been detected and map to the mOTUs database, but half of the genes map to one mOTU and the other half to another mOTU;
* Consistent, the majority of the genes agree on one mOTU.


Estimate correctness of the assignment
--------------

To test the classification accuracy we did a simulation (details below) where we removed some genomes from the database and we used them to test the `classify-genomes`. You can see the results for the specificity of the classification (a proxy for the probability of correct assignment) in the following plot:

![alt text](https://github.com/AlessioMilanese/classify-genomes/blob/master/pics/prec_recall.png)

This means that, if you have:
```
Number of mapped genes: 5
Number of agreeing genes: 4
```
for your classification, then the estimated probability to be correct is `0.98` (or `98%`).

#### Details of the simulation

We remove some or all the genomes from 1/10th of the ref-mOTUs (~1,200 mOTUs).

For 75% of the ~1,200 mOTUs we remove all genomes from the database, and we used one genome to test if it was correctly not assign to any mOTU. These genomes can contribute only to the false positives in the evaluation.

For 25% of the ~1,200 mOTUs we removed half of the genomes. We chose one of the removed genomes and test if it is classified to the correct mOTU. Since for this set some genomes from the correct mOTU are available, it is possible to evaluate the true positives.

The simulation is quite conservative, having 3 times more possible false positives, than true positives.


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
* **`-S`** map to the specI database (default: map to the mOTUs database)
* **`-db`** choose the mOTUs database: nr=non-reduntant, cen=centroids. Default is nr.
