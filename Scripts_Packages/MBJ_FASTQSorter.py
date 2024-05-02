import sys,os
from Bio import SeqIO
import datetime
import time
import argparse
import csv
import gzip
import pandas as pd
import re

# Adapted from:
#	Song X, Gao HY, Herrup K, Hart RP. Optimized splitting of mixed-species RNA sequencing data. 
#	J Bioinform Comput Biol. 2022 Apr;20(2):2250001. doi: 10.1142/S0219720022500019.
#
# Steps:
# 1. Read SAM input and parse into:
#    a. dictionary of read numbers pointing to tuple (immutable) required aligment info
#    b. list of unique read IDs pointing to list (updatable) of rows in dict
# 2. Report counts
# 3. Decide best genome:
#    a. If ID matches only one genome, that's it
#    b. If ID has multiple genome matches, choose primary alignment
#    c. (future) if any reads remaining, choose best alignment score
# 4. Output tables for splitting fastq files in a separate program

# Global variables with standard naming convension
VERSION = '1.0'
DESCRIPTION = '''\
    The goal is to separate BULK RNAseq fastq files into two species specific subsets
    We iterate through the mixed-genome SAM file to and choose the best score by genome.
    Read IDs matching both genomes are counted.
	Open SAM input, parse into a dictionary of alignment info and
	a dictionary of read IDs (pointing to a tuple of alignment records).
	'''
def checkOptions():
    parser = argparse.ArgumentParser(usage='samtools view <BAM> | python byPrim.py',description=DESCRIPTION)
    parser.add_argument('-v','--version', action='version', version=VERSION)
    parser.add_argument('-s','--sample',dest="sample_name",required=True,default=None,help='Sample ID for insertion into file name')
    parser.add_argument('-i','--FastIn',dest="fastq_in",required=True,default=None,help='Preprocessed FASTQ file input directory')
    parser.add_argument('-o','--FastOut',dest="fastq_out",required=True,default=None,help='Sorted FASTQ file output directory')

    args = parser.parse_args()
    #if args.sample_name:
        #args.sample_name = args.sample_name.upper()
    return args


def main():

	opt = checkOptions()
	# all options are accessible as opt.var
	# print(opt)
	sampleName = opt.sample_name
	fastqInput = opt.fastq_in
	fastqOutput = opt.fastq_out
	
	# Create the input variable for forward FASTQ sorting
	forInput = fastqInput + "/QualTrim_" + sampleName + "_R1.fastq.gz"		
	
	# Initializing temporary lists to hold the FASTQ reads
	fastTEMP1 = []
	fastTEMP2 = []
	MsLine = []
	HumLine = []
	RevFile = []
	MsRevReads = []
	HumRevReads = []
	lineno = 0
	
	# Load the relevant read ID tables
	print("\nLoading the read ID tables for " + sampleName + "...")
	ID = pd.read_csv("./" + sampleName + "_table.tsv.gz", sep='\t', header=None)
	MIX = pd.read_csv("./" + sampleName + "_mixtable.tsv.gz", sep='\t', header=None)
	
	# Sort the tables into the necessary read IDs
	SAPTab = ID[(ID.iloc[:,3] == 'Sapiens')].iloc[:,0]
	MIXTab = MIX.iloc[:,0]
	
	# Convert dataframes to list
	SAPTabFor = SAPTab.tolist()
	MIXTabFor = MIXTab.tolist()
	SAPTabRev = SAPTab.tolist()
	MIXTabRev = MIXTab.tolist()
	
	# Cleaning up dataframes that are no longer needed
	lst = [ID, MIX, SAPTab, MIXTab]
	del ID
	del MIX
	del SAPTab
	del MIXTab
	del lst
	
	# Load forward FASTQ file line by line
	print("Sorting the forward FASTQ file " + forInput + "...")
	with gzip.open(forInput, 'rt') as f:
		for line in f:
			lineno+=1
			if lineno%4000000 == 0:
				print("Forward FASTQ for " + sampleName + " is on read " + str(lineno/4))
			if lineno%4 == 1:
				readID = re.split(r'@|\s',line)[1]
				if (readID in MIXTabFor):
					flag = "Neither"
					MIXTabFor.remove(readID)
				if (readID in SAPTabFor):
					flag = "Human"
					SAPTabFor.remove(readID)
				else:
					flag = "Mouse"
			if flag == 'Mouse':
				fastTEMP1.append(line)
				MsLine.append(lineno-1)
			if flag == 'Human':
				fastTEMP2.append(line)
				HumLine.append(lineno-1)
			else:
				continue
		
	# Writing the mouse and human R1 FASTQ files
	print("Writing the HUMAN forward FASTQ file for " + sampleName + "...")
	with gzip.open(fastqOutput + "/Human/Sorted_Human_" + sampleName + "_R1.fastq.gz", "wt") as HuForOutput:
		for row in fastTEMP2:
			HuForOutput.write(str(row))
	print("Finished writing the HUMAN forward FASTQ file for " + sampleName + "...")
	print("Writing the MOUSE forward FASTQ file for " + sampleName + "...")
	with gzip.open(fastqOutput + "/Mouse/Sorted_Mouse_" + sampleName + "_R1.fastq.gz", "wt") as MsForOutput:
		for row in fastTEMP1:
			MsForOutput.write(str(row))
	print("Finished writing the MOUSE forward FASTQ file for " + sampleName + "...\n")
	
	# Clearing the temporary fastq lists for use with the reverse file
	fastTEMP1.clear()
	fastTEMP2.clear()

	# Create the input variable for reverse FASTQ sorting
	revInput = fastqInput + "/QualTrim_" + sampleName + "_R2.fastq.gz"
	
	# Loading the entire reverse file into memory
	print("Loading the reverse FASTQ file " + revInput + "...")
	with gzip.open(revInput, 'rt') as f:
		for line in f:
			RevFile.append(line)
	print("Finished loading the reverse FASTQ file for " + revInput)
	
	# Subsetting the Reverse file based on the forward file line indices		
	print("Subsetting the MOUSE reads for " + revInput + "...")
	MsRevReads = [RevFile[i] for i in MsLine]
	print ("Finished subsetting the HUMAN reads for " + revInput)
	print("Subsetting the MOUSE reads for " + revInput + "...")
	HumRevReads = [RevFile[j] for j in HumLine]
	print ("Finished subsetting the HUMAN reads for " + revInput)
	
	# Deleting the lists to free RAM
	print("Deleting the lists for " + revInput + "...")
	lst = [RevFile, MsLine, HumLine]
	del RevFile
	del MsLine
	del HumLine
	del lst
		
	# Writing the mouse and human R2 FASTQ files
	print("Writing the HUMAN reverse FASTQ file for " + sampleName + "...")
	with gzip.open(fastqOutput + "/Human/Sorted_Human_" + sampleName + "_R2.fastq.gz", "wt") as HuRevOutput:
		for row in HumRevReads:
			HuRevOutput.write(str(row))
	print("Finished writing the HUMAN reverse FASTQ file for " + sampleName + "...")
	
	print("Writing the MOUSE reverse FASTQ file for " + sampleName + "...")
	with gzip.open(fastqOutput + "/Mouse/Sorted_Mouse_" + sampleName + "_R2.fastq.gz", "wt") as MsRevOutput:
		for row in MsRevReads:
			MsRevOutput.write(str(row))
	print("Finished writing the MOUSE reverse FASTQ file for " + sampleName + "...\n")
	
	# Clearing the temporary fastq lists for use with the reverse file
	MsRevReads.clear()
	HumRevReads.clear()
	
	sys.exit(0)

if __name__ == '__main__':
	main()
