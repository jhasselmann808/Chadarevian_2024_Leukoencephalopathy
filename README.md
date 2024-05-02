# Chadarevian_2024_Leukoencephalopathy
The pipeline included in this repository is designed to run on Linux systems and in conjuction with the reference files located at https://drive.google.com/drive/folders/1gzf8VM6-btcsj0bY6BxgDQj6_EOcFGLD?usp=sharing. Combined, this will provide all of the necessary packages, scripts, metadata, and reference files needed to generate the analysis that was used in the Chadarevian et al. (Neuron, 2024) publication. The FASTQ files will be automatically downloaded from GEO (GSE253623) and renamed to work with the downstream analysis. If you have previously downloaded the FASTQ files, you will need to ensure that the file names align with the filenames listed in the "./Metadata/Chadarevian_2024_Leukoencephalopathy_Bulk_Metadata.tsv" file

This pipeline uses a number of tools made available by other people/groups. If you use parts of this pipeline for subsequent analysis, please be sure to properly cite the original creators of the tools.

	FASTQC
		Andrews, S. (2014). FastQC: A quality control tool for high throughput sequencing data.

	BBDUK
		Bushnell B. - sourceforge.net/projects/bbmap/

	byPrim.py
		Song X, Gao HY, Herrup K, Hart RP. Optimized splitting of mixed-species RNA sequencing data. 
		J Bioinform Comput Biol. 2022 Apr;20(2):2250001. doi: 10.1142/S0219720022500019.

	HISAT2
		Kim, D., Paggi, J.M., Park, C. et al. Graph-based genome alignment and genotyping with HISAT2 and 
		HISAT-genotype. Nat Biotechnol 37, 907–915 (2019). doi: 10.1038/s41587-019-0201-4

	Kallisto
		Nicolas L Bray, Harold Pimentel, Pall Melsted and Lior Pachter. Near-optimal
		probabilistic RNA-seq quantification. Nature Biotechnology. 34, 525-27 (2016)
		doi:10.1038/nbt.3519

	DESeq2/tximport
		Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and
		dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)
				    
		Soneson C, Love MI, Robinson MD (2015). “Differential analyses for RNA-seq: 
		transcript-level estimates improve gene-level inferences.” F1000Research, 4. 
		doi: 10.12688/f1000research.7563.1.
 	
	Tidyverse
		Wickham H, Averick M, Bryan J, et al. (2019). “Welcome to the tidyverse.” 
 		Journal of Open Source Software, 4(43), 1686. doi:10.21105/joss.01686. 
	
	ggplot2
		Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag 
		New York. ISBN 978-3-319-24277-4, https://ggplot2.tidyverse.org. 

To sucessfully run the pipeline on your system, perform the following:

	1. Download this repository and save it into your desired analysis directory (~1TB space needed for full analysis)
 	2. Download the "References.tar.xz" from https://drive.google.com/drive/folders/1gzf8VM6-btcsj0bY6BxgDQj6_EOcFGLD?usp=sharing and extract the directory to ./Chadarevian_2024_Leukoencephalopathy/
	2. Ensure that the scripts in the ./Chadarevian_2024_Leukoencephalopathy/Scripts_Packages/ directory are executable
 		- chmod -R 755 ./Chadarevian_2024_Leukoencephalopathy/Scripts_Packages
	3. Ensure the pigz package is installed on the system
		- Terminal > sudo apt install pigz
	4. Ensure R is installed on the system (https://cran.r-project.org/bin/linux/ubuntu/fullREADME.html)
		- The necessary R packages will be installed automatically
	5. Create a Python environment that has access to the following packages:
		- biopython
		- pandas
		- datetime
		- time
		- argparse
		- csv
		- gzip
		- re
	6. Review lines 91-97 of the ./Chadarevian_2024_Leukoencephalopathy/Scripts_Packages/Chadarevian_2024_Leukoencephalopathy_RNASeqProcessing.sh script to ensure settings work for your system.
	7. Open a Terminal and activate the Python environment
	8. Run the Chadarevian_2024_Leukoencephalopathy_RNASeqProcessing.sh script and wait, and wait, and wait..........
