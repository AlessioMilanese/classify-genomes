
import os
import sys
import subprocess
import glob
import time
import argparse
from collections import defaultdict
import shlex
bindir = os.path.abspath(os.path.dirname(__file__))
#input: contigs fastas
#gene prediction (wuerde das mit prodigal machen, output ist ein gbk oder gff und die sequences)
#fetchMG
#dann fuer die fetchMG genes die +-100 bp extracten (shini kannst du jaime fragen ob wir sein script da benutzten koennen, haben ja das gbk)
#db file generieren (die daten sind ja alle da, ausser der motu zugehoerigkeit, aber die muss spaeter gemacht werden)

def readSingleColumnFileToList(infileName):

	infile = open(infileName, "r")
	listInfileData = []
	for strLine in infile:
		strLine = strLine.rstrip('\n')
		if (" " in strLine):
			arrLine = strLine.split(' ')
			strLine = arrLine[0]
		if ("\t" in strLine):
			arrLine = strLine.split('\t')
			strLine = arrLine[0]

		listInfileData.append(strLine)

	return(listInfileData)


def parse2columnFile(infileName):
	dictIn = defaultdict(list)
	dictInvOut = defaultdict(list)

	infile = open(infileName, "r")

	for strLine in infile:

		strLine = strLine.rstrip('\n')
		arrLine = strLine.split('\t')

		strColumn_1_entry = arrLine[0]
		strColumn_2_entry = arrLine[1]

		dictIn[strColumn_1_entry] = strColumn_2_entry
		dictInvOut[strColumn_2_entry].append(strColumn_1_entry)

	return(dictIn, dictInvOut)


def predictGenes(contigsFasta, metagenomicMode, transTableNumber, protOut, nuclOut, gffOut):

	if (metagenomicMode):
		mode = " -p meta"
	else:
		mode = " -p single"


	transTableParameter = ""
	if (transTableNumber != "11"):
		transTableParameter = " -g " + str(transTableNumber)
		#prodigal_cmd = subprocess.Popen(['prodigal', '-f gff', '-q', '-c', '-a ' + str(protOut), '-d ' + str(nuclOut), '-o ' + str(gffOut), transTableParameter, '-m', mode, '-i ' + contigsFasta])
	#else:
		#prodigal_cmd = subprocess.Popen(['prodigal', '-f gff', '-q', '-c', '-a ' + str(protOut), '-d ' + str(nuclOut), '-o ' + str(gffOut), '-m' , mode, '-i ' + contigsFasta])

	cmd = "prodigal " + "-a " + str(protOut) + " -d " + str(nuclOut) + " -f gff -o " + str(gffOut) + " -c" + str(transTableParameter) + " -q -m " + mode + " -i " + contigsFasta
	print (cmd)
	popenCMD = shlex.split(cmd)
	prodigal_cmd = subprocess.Popen(popenCMD)
	prodigal_cmd.wait()

	#cmd = "prodigal " + "-a " + protOut + " -d " + nuclOut + " -f -o " + gffOut + " -c " + transTableParameter + " -q -m " + mode + " " + contigsFasta
	#os.system(cmd)


def runFetchMGs(protFasta, genFasta, numThreads, outMGfolder):

	#if i want to use this I need to reformat the genomes
	#if (metagenomicMode):
	#	besthit = ""
	#else:
	#	mode = " -v"

	#cmd = bin + "MOCATFetchMGs05.pl " + "-t " + str(numThreads) " -d " + nuclFasta + " -o " + outMGfolder + " " + protFasta
	#os.system(cmd)
	cmd = "perl " + str(bindir) + "/fetchMG.pl -m extraction -t " + str(numThreads) + " -o " + str(outMGfolder) + " -d " + str(genFasta)  + " " + str(protFasta)
	print (cmd)
	popenCMD = shlex.split(cmd)
	fetchMGs_cmd = subprocess.Popen(popenCMD)
	fetchMGs_cmd.wait()


def parseFetchMGs(outMGfolder, OGType):
	setMGs = set()
	glob.glob('*.gif')

	OGs2use = ["all"]
	if (OGType == "mOTU"):
		OGs2use = ["COG0012", "COG0016", "COG0018", "COG0172", "COG0215", "COG0495", "COG0525", "COG0533", "COG0541", "COG0552"]

	currentPath= os.path.sep.join([outMGfolder, "*.all.marker_genes_scores.table"])
	for currentFileName in glob.glob(currentPath):
		with open(currentFileName,'r') as currentFile:
			for strLine in currentFile:
				strLine = strLine.rstrip('\n')
				arrLine = strLine.split('\t')

				mgID = arrLine[0]
				currOG = arrLine[2]
				if ( currOG in OGs2use or "all" in OGs2use):
					setMGs.add(mgID)


	return(setMGs)

def MakeContigLengthFile(contigsFasta, contigLengthsFileName):
#	cmd = 'seqtk comp ' + contigsFasta + ' | cut -f1,2 > ' + contigLengthsFileName
#	os.system(cmd)

# recoded above in pure python
	contigLengthsFile = open(contigLengthsFileName,"w")
	seqtk_comp = subprocess.Popen(['seqtk', 'comp', contigsFasta], stdout=subprocess.PIPE)

	cut = subprocess.Popen(['cut', '-f1,2'], stdin=seqtk_comp.stdout, stdout=contigLengthsFile)
	cut.wait()

#genomes/SAGs/assemblies3/SAR324_J029/SAR324_J029.contigs.fasta

#MakeContigLengthFile("genomes/SAGs/assemblies3/SAR324_J029/SAR324_J029.contigs.fasta", "temp.txt")
#cmd = seqtk comp genomes/SAGs/assemblies3/SAR324_J029/SAR324_J029.contigs.fasta | cut -f1,3 > temp2.txt


def parseGeneLocations_contigCoords(setMGs, contigCoordsFileName, dictContig2Length, padLength, paddedLocationsFileName, extractLocationsFileName, geneLengthFileName, paddedLengthFileName):

	contigCoordsFile = open(contigCoordsFileName, "r")

	paddedLocationsFile = open(paddedLocationsFileName, "w")
	geneLengthFile = open(geneLengthFileName, "w")
	extractLocationsFile = open(extractLocationsFileName, "w")
	paddedLengthFile = open(paddedLengthFileName, "w")

	paddedLocationsFileSpace = open(paddedLocationsFileName + ".ssv" , "w")

	dictExtractName2writeName = defaultdict(list)
	padLength = int(padLength)

	for strLine in contigCoordsFile:
		strLine = strLine.rstrip('\n')
		arrLine = strLine.split('\t')
		geneName = arrLine[0]
		contigID = arrLine[1]
		start = int(arrLine[2])
		start = max(start, 1)
		end = int(arrLine[3])
		length = end - start + 1

		if (geneName in setMGs or "all" in setMGs):

			if (contigID in dictContig2Length):

				contigLength = int(dictContig2Length[contigID])

				paddedStart = max(start - padLength, 1)
				paddedEnd = min(end + padLength, contigLength)

				newGeneStart = start - paddedStart + 1
				newGeneEnd = newGeneStart + length - 1

				#geneName contigID ID in padded file start of gene on padded seq
				writeline = geneName + "\t" + contigID + "\t" + str(newGeneStart) + "\t" + str(newGeneEnd) + "\n"
				paddedLocationsFile.write(writeline)
				writeline = geneName + " " + contigID + " " + str(newGeneStart) + " " + str(newGeneEnd) + "\n"
				paddedLocationsFileSpace.write(writeline)

				#lyrata:1-108
				extractLine = contigID + ":" + str(paddedStart) + "-" + str(paddedEnd)
				extractLocationsFile.write(extractLine + "\n")

				paddedGeneName = geneName + " " + str(newGeneStart) + " " + str(newGeneEnd)
				dictExtractName2writeName[extractLine].append(paddedGeneName)

				lengthLine = geneName + "\t" + str(length)
				geneLengthFile.write(lengthLine + "\n")

				lengthLine2 = geneName + "\t" + str(paddedEnd - paddedStart + 1)
				paddedLengthFile.write(lengthLine2 + "\n")
			else:
				print("Could not find length for " + str(contigID) + " . Please check fasta or length file.")
	return(dictExtractName2writeName)



def parseGeneLocations_gff(setMGs, gffFileName, padLength, paddedLocationsFileName, extractLocationsFileName, geneLengthFileName, paddedLengthFileName):

	gffFile = open(gffFileName, "r")

	paddedLocationsFile = open(paddedLocationsFileName, "w")
	extractLocationsFile = open(extractLocationsFileName, "w")
	geneLengthFile = open(geneLengthFileName, "w")
	paddedLengthFile = open(paddedLengthFileName, "w")
	paddedLocationsFileSpace = open(paddedLocationsFileName + ".ssv" , "w")
	dictExtractName2writeName = defaultdict(list)
	padLength = int(padLength)

	#SAR324_J029_c172        Prodigal_v2.6.1 CDS     33      836     91.7    +       0       ID=170_1;partial=00;start_type=ATG;rbs_motif=GGAGG;rbs_spacer=5-10bp;gc_cont=0.417;conf=100.00;score=91.69;cscore=79.91;sscore=11.78;rscore=8.57;uscore=-0.99;tscore=4.21;
	for strLine in gffFile:
		strLine = strLine.strip()
		dictContigInfo = {}
		# Sequence Data: seqnum=1;seqlen=4639675;seqhdr="NC_000913 # Escherichia coli str. K-12 substr. MG1655, complete genome."
		if(strLine.startswith("# Sequence Data:")):
			strLine = strLine.replace("# Sequence Data:", "")
			strLine = strLine.split("#")[0]
			strLine = strLine.strip()

			infoList = strLine.split(';')
			for infoPart in infoList:
				infoSplit = infoPart.split('=')
				dictContigInfo[infoSplit[0]] = infoSplit[1]

			contigLength = int(dictContigInfo["seqlen"])

		elif(not strLine.startswith("#")):
			arrLine = strLine.split('\t')
			contigID = arrLine[0]
			strColumn_2_entry = arrLine[1]

			start = int(arrLine[3])
			start = max(start, 1)
			end = int(arrLine[4])
			#strand = arrLine[6]
			length = end - start + 1

			dictInfo = {}
			info = arrLine[8]
			infoList = info.split(';')
			for infoPart in infoList:
				if ('=' in infoPart):
					infoSplit = infoPart.split('=')
					dictInfo[infoSplit[0]] = infoSplit[1]

			paddedStart = max(start - padLength, 1)
			paddedEnd = min(end + padLength, contigLength)

			newGeneStart = start - paddedStart + 1
			newGeneEnd = newGeneStart + length - 1

			if "name" in dictInfo:
				geneName = dictInfo[name]
			elif "ID" in dictInfo:
				geneNumber = dictInfo["ID"].split("_")[1]
				geneName = contigID + "_" + geneNumber
				#print geneName
			else:
				print("Cannot find ID in GFF, please provide gff in prodigal format")

			if (geneName in setMGs or "all" in setMGs):
				#geneName contigID ID in padded file start of gene on padded seq
				#print geneName
				writeline = geneName + "\t" + contigID + "\t" + str(newGeneStart) + "\t" + str(newGeneEnd) + "\n"
				paddedLocationsFile.write(writeline)
				writeline = geneName + " " + contigID + " " + str(newGeneStart) + " " + str(newGeneEnd) + "\n"
				paddedLocationsFileSpace.write(writeline)

				#lyrata:1-108
				extractLine = contigID + ":" + str(paddedStart) + "-" + str(paddedEnd)
				extractLocationsFile.write(extractLine + "\n")

				lengthLine = geneName + "\t" + str(length)
				geneLengthFile.write(lengthLine + "\n")

				lengthLine2 = geneName + "\t" + str(paddedEnd - paddedStart + 1)
				paddedLengthFile.write(lengthLine2 + "\n")

				paddedGeneName = geneName + " " + str(newGeneStart) + " " + str(newGeneEnd)
				dictExtractName2writeName[extractLine].append(paddedGeneName)

	return(dictExtractName2writeName)

#==> data/mOTU.v1.padded <==
#>1048834.TC41_3292 101 1561

#==> data/mOTU.v1.padded.coord <==
#1048834.TC41_3292 101 1561

#==> data/mOTU.v1.padded.len <==
#1048834.TC41_3292       1661
# ==> change to length of gene?

def getPaddedSequences(extractLocationsFileName, contigsFasta, paddedFastaFileName, dictExtractName2writeName):

	paddedFastaFile = open(paddedFastaFileName, "w")

	cat = subprocess.Popen(['cat', extractLocationsFileName],
							stdout=subprocess.PIPE,
							)

	samtools = subprocess.Popen(['xargs', 'samtools', 'faidx', contigsFasta],
							stdin=cat.stdout,
							stdout=subprocess.PIPE,
							)

	output = samtools.stdout

	for line in output:
		if (line.startswith('>')):
			line = line.strip()
			recordName = line.replace('>', '')
			newRecordName = dictExtractName2writeName[recordName].pop()
			newRecordName = newRecordName.strip()
			paddedFastaFile.write(">" + newRecordName + "\n")
		else:
			paddedFastaFile.write(line)



def runMakePaddedMG_DB(contigsFasta, outfolder, runPrefix, gffFileName, contigCoordsFileName, contigLengthsFileName, padLength, geneNamesFileName, proteinFasta, transTableNumber, metagenomicMode, numThreads, OGType):

	#actual input

	#contigsFasta
	#outfolder
	#gffFileName (optional)
	#contigCoordsFileName (optional)
	#runPrefix
	#padLength (optional)
	#geneNamesFileName (optional)
	#genesFasta (optional)
	#proteinFasta (optional)
	#transTableNumber (optional; only genepred)
	#metagenomicMode (optional; only genepred)
	#numThreads (optional; only fetchmgs)


	if not os.path.exists(outfolder):
		os.makedirs(outfolder)
	else:
		print("Outfolder exists, hopefully nothing important will be overwritten...waiting 5 seconds so you can emergency cancel")
		time.sleep(5)

	#internal (files are produced by script)
	outMGfolder = os.path.sep.join([outfolder, runPrefix + "_markerGenes"])
	extractLocationsFileName = os.path.sep.join([outfolder, runPrefix + ".extract.coords"])
	paddedFastaFileName = os.path.sep.join([outfolder, runPrefix + ".padded.fasta"])
	paddedLengthFileName = os.path.sep.join([outfolder, runPrefix + ".padded.len"])
	paddedLocationsFileName = os.path.sep.join([outfolder, runPrefix + ".padded.coords"])

	geneLengthFileName = os.path.sep.join([outfolder, runPrefix + ".genes.len"])

	contigLengthGiven = False
	if (contigLengthsFileName):
		contigLengthGiven = True
	else:
		contigLengthsFileName  = os.path.sep.join([outfolder, runPrefix + ".contigs.lengths"])



	#this section is about the genes for which to extract the padded sequence. either they are found using fetchMGs or they are preselected by the user using geneNamesFileName ('all' can be specified to get padded sequence for all genes)
	setMGs = set()

	#first step: predict genes or not?
	if ((not gffFileName) and (not contigCoordsFileName)):
		#print("Running prodigal...")
		gffFileName = os.path.sep.join([outfolder, runPrefix + ".all.genes.gff"])
		genesFasta = os.path.sep.join([outfolder, runPrefix + ".all.genes.fasta"])
		proteinFasta = os.path.sep.join([outfolder, runPrefix + ".all.proteins.fasta"])
		predictGenes(contigsFasta, metagenomicMode, transTableNumber, proteinFasta, genesFasta, gffFileName)
		#print("...done")

	#second step: search for MGs or not?
	if (proteinFasta and not geneNamesFileName):
		#print("Running fetchMGs...")
		#runFetchMGs(proteinFasta, genesFasta, numThreads, outMGfolder)
		runFetchMGs(proteinFasta, genesFasta, numThreads, outMGfolder)
		setMGs = parseFetchMGs(outMGfolder, OGType)
		#print("...done")

	if(geneNamesFileName == "all"):
		setMGs.add("all")
	elif(geneNamesFileName):
		listMGs = readSingleColumnFileToList(geneNamesFileName)
		setMGs = set(listMGs)
		#print(setMGs)

	if (len(setMGs) == 0):
		print("No genes provided for extraction, provide parameter 'XXX=all' to extract all padded sequences for all genes.")
		sys.exit(1)


	#this part is about getting the locations of the genes from either provided files or from gff files from geneprediction above
	if (gffFileName):
		dictExtractName2writeName = parseGeneLocations_gff(setMGs, gffFileName, padLength, paddedLocationsFileName, extractLocationsFileName, geneLengthFileName, paddedLengthFileName)

		#parseGeneLocations_gff(setMGs, gffFileName, padLength, paddedLocationsFileName, extractLocationsFileName, geneLengthFileName, paddedLengthFileName)

	elif(contigCoordsFileName and contigLengthsFileName):
		if (not contigLengthGiven):
			MakeContigLengthFile(contigsFasta, contigLengthsFileName)
		dictContig2Length, dictLength2Contig = parse2columnFile(contigLengthsFileName)
		dictExtractName2writeName = parseGeneLocations_contigCoords(setMGs, contigCoordsFileName, dictContig2Length, padLength, paddedLocationsFileName, extractLocationsFileName, geneLengthFileName, paddedLengthFileName)
	else:
		print("No gff or gene coordinate file provided; also gene calling disabled. Cannot parse gene locations on contigs. Exiting.")
		sys.exit(1)

	#And now the padded sequences are actually extracted
	if(len(dictExtractName2writeName.keys()) > 0):
		getPaddedSequences(extractLocationsFileName, contigsFasta, paddedFastaFileName, dictExtractName2writeName)

#needed files
# contigs fasta
# dependencies


def main(argv=None):

	parser = argparse.ArgumentParser(description='This script is used to extract padded sequences from contigs. Can predict genes and extract marker genes in the process. Part of the mOTUs v2 pipeline', add_help = True)
	parser.add_argument('contigsFasta', action="store", help='Fasta file with contig sequences')
	parser.add_argument('outfolder', action="store", help='Folder for output')
	parser.add_argument('runPrefix', action="store", help='Name for this run')
	parser.add_argument('--gffFile', '-g', action="store", dest='gffFileName', help='gff formatted file that provides the location of genes. Gene names should be encoded in the attribute section as either ID (ID=*_XX where contigname_XX is the gene name; this is prodigal style) or name (name=geneName). If neither --gffFile nor -contigCoordsFile provided, de novo gene prediction will be performed.')
	parser.add_argument('--contigCoordsFile', '-c',action="store", dest='contigCoordsFileName', help='Alternative to --gffFile. Tab delimited file with geneName\\tcontigName\\tgeneStart\\geneEnd. If neither --gffFile nor -contigCoordsFile provided, de novo gene prediction will be performed.')
	parser.add_argument('--contigLengthsFile', '-s',action="store", dest='contigLengthsFileName', help='Can be provided with --contigCoordsFile. Tab delimited file with contigName\\tcontigLength. If not given contig length will be calculated from fasta file.')
	parser.add_argument('--padLength', '-l', action='store', dest='padLength', type=int, default=100, help="Lengths of the padding at each side of the gene, default is 100bp.")
	parser.add_argument('--geneNamesFile', '-n', action='store', dest='geneNamesFileName', help="File with names of genes to be padded. If not set and --proteinFasta provided (or when gene prediction is run), fetchMGs will be run to extract the 40 single copy marker genes. If set to 'all', all gene will be padded.")
	parser.add_argument('--proteinFasta', '-p', action='store', dest='proteinFasta', help="If gene prediction is not run, but marker genes should be extracted using fetchMGs and padded, please provide the protein sequences using this parameter and do not set --geneNamesFile")
	parser.add_argument('--transTableNumber', '-t',action="store", dest='transTableNumber', type=int, default=11,help='Select a translation table for prodigal gene prediction, default is 11.')
	parser.add_argument('--metagenomicMode', '-m',action="store", dest='metagenomicMode', default="meta",help="Select a procedure for prodigal gene prediction, default is meta; alternative is 'single'.")
	parser.add_argument('--OGType', '-o',action="store", dest='OGType', default="all",help="If running fetchMGs, set to 'mOTU' to constrain to 10 mOTU MGs.")
	parser.add_argument('--numThreads', '-a',action="store", dest='numThreads', type=int, default=1, help="Number of threads for fetchMGs.")
	parser.add_argument('--version', action='version', version='%(prog)s 1.0')
	args = parser.parse_args()

	numThreads = str(args.numThreads)
	transTableNumber = str(args.transTableNumber)

	runMakePaddedMG_DB(args.contigsFasta, args.outfolder, args.runPrefix, args.gffFileName, args.contigCoordsFileName, args.contigLengthsFileName,args.padLength, args.geneNamesFileName, args.proteinFasta, transTableNumber, args.metagenomicMode, numThreads, args.OGType)

	return 0        # success

if __name__ == '__main__':
	status = main()
	sys.exit(status)
