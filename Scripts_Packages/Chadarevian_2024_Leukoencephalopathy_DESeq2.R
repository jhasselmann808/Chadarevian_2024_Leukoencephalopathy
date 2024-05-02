# ==========================================================================================================  
# ======================================= Loading Required Packages ========================================
# ========================================================================================================== 

# Loading required CRAN packages
# If packages are not available, they will be installed and loaded
CRANpkg <- c("BiocManager", "remotes", "ggplot2", "ggrepel", 
             "pheatmap", "reshape2", "data.table", "tidyverse", 
             "gtable", "grid")
CRANreq <- unlist(lapply(CRANpkg, require, character.only=T))
CRANneed <- CRANpkg[CRANreq==F]

if (length(CRANneed)>0){
  install.packages(CRANneed)
  lapply(CRANneed, require, character.only=T)
}

# Loading required Bioconductor packages
# If packages are not available, they will be installed and loaded
BIOCpkg <- c("DESeq2", "tximport")
BIOCreq <- unlist(lapply(BIOCpkg, require, character.only=T))
BIOCneed <- BIOCpkg[BIOCreq==F]

if (length(BIOCneed)>0){
  BiocManager::install(BIOCneed)
  lapply(BIOCneed, require, character.only=T)
  }

# Loading required GitHub packages
# If packages are not available, they will be installed and loaded
GITpkg <- c("threadr")
GITreq <- unlist(lapply(GITpkg, require, character.only=T))
GITneed <- GITpkg[GITreq==F]

if (length(GITneed)>0){
  remotes::install_github(GITneed)
  lapply(GITneed, require, character.only=T)
  }

# ==========================================================================================================  
# ======================================= Defining Custom Functions ========================================
# ========================================================================================================== 

# A modified PCA function that assigns shapes and colors by group
New_PCA <- function (vstObj, var1 = NULL, var2 = NULL, ntop = NULL, labs = NULL){
  
  pcaData <- plotPCA(vstObj, intgroup=c(var1, var2), ntop=ntop, returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  if (leg_pos == "none"){
    marg <- c(0.5, 0.5, 0, 0)
  } else {
    marg <- c(0,0,0,0)
  }
  
  if (!(labs)){
    if (is.null(var2)){
      pca <- ggplot(pcaData, aes(PC1, PC2, color=get(var1), label=NULL)) + 
        geom_point(size=6) +
        theme_bw() +
        theme(axis.line=element_line(linewidth=2), 
              axis.ticks = element_line(color='black', linewidth=2), 
              axis.ticks.length = unit(0.2, 'cm'), 
              axis.text = element_text(face='bold', color ='black', size=20),
              text = element_text(face='bold', color ='black', size=20),
              panel.border = element_blank(),
              plot.margin = unit(marg, "cm"),
              legend.position='right') +
        labs(color = var1) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance"))
      
    } else {
      
      pca <- ggplot(pcaData, aes(PC1, PC2, color=get(var1), shape=get(var2), label=NULL)) + 
        #scale_color_manual(values = c("#0072B2","red", "#D55E00")) +
        geom_point(size=6) +
        theme_bw() +
        theme(axis.line=element_line(linewidth=2), 
              axis.ticks = element_line(color='black', linewidth=2), 
              axis.ticks.length = unit(0.2, 'cm'), 
              axis.text = element_text(face='bold', color ='black', size=20),
              text = element_text(face='bold', color ='black', size=20),
              panel.border = element_blank(),
              plot.margin = unit(marg, "cm"),
              legend.position='right') +
        labs(color = var1, shape = var2) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance"))
      
    }
  } else if (labs) {
    
    if (is.null(var2)){
      
      pca <- ggplot(pcaData, aes(PC1, PC2, color=get(var1), label=name)) + 
        geom_point(size=6) +
        theme_bw() +
        geom_text_repel(size = 3) +
        theme(axis.line=element_line(linewidth=2), 
              axis.ticks = element_line(color='black', linewidth=2), 
              axis.ticks.length = unit(0.2, 'cm'), 
              axis.text = element_text(face='bold', color ='black', size=20),
              text = element_text(face='bold', color ='black', size=20),
              panel.border = element_blank(),
              plot.margin = unit(marg, "cm"),
              legend.position='right') +
        labs(color = var1) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance"))
      
    } else {
      
      pca <- ggplot(pcaData, aes(PC1, PC2, color=get(var1), shape=get(var2), label=name)) + 
        geom_point(size=6) +
        theme_bw() +
        geom_text_repel(size = 3) +
        theme(axis.line=element_line(linewidth=2), 
              axis.ticks = element_line(color='black', linewidth=2), 
              axis.ticks.length = unit(0.2, 'cm'), 
              axis.text = element_text(face='bold', color ='black', size=20),
              text = element_text(face='bold', color ='black', size=20),
              panel.border = element_blank(),
              plot.margin = unit(c(0,0,0,0), "cm"),
              legend.position='right') +
        labs(color = var1, shape = var2) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance"))
      
    }
  }
  return(pca)
  
}

# A function that creates a list of FDR and LFC cutoff genes from the DESeq results
Make_Folded_Genes <- function (object, j, dat, FDR, LFC){
  
  # Extract the results from the DESeq object depending on the type of comparison
  if (length(dat[[j]]) == 1) {
    res <- results(object, name = dat[[j]], alpha = FDR)
  } else if (length(dat[[j]]) == 2){
    res <- results(object, list( dat[[j]] ), alpha = FDR)
  }
  # Remove all NAs from results
  nonNA_res <- res[!is.na(res$padj),]
  pcutoff_res <- nonNA_res[(nonNA_res$padj<FDR),]
  # Apply a log-fold change cutoff to the results
  folded_res <- pcutoff_res[(pcutoff_res$log2FoldChange >= LFC | pcutoff_res$log2FoldChange <= -LFC),]
  folded_genes <- data.frame(folded_res)
  folded_genes$Gene <- rownames(folded_genes)
  folded_genes <- folded_genes[, c(7,1:6)]
  folded_genes <- folded_genes[rev(order(folded_genes$log2FoldChange)),]
  
  return(folded_genes)
  
}

# A function that generates and saves PCA plots in conjuction with New_PCA
savePCA <- function (object, wid=8, hei=6, reso=300) {
  # Setting the plot size based on whether there is a legend or not
  if (leg_pos == 'none'){
    wide <- 6
    tall <- 6
  }else {
    wide <- 8
    tall <- 6
  }
  
  # Perform VST on the DESeq2 object
  vstObj <- varianceStabilizingTransformation(object, blind = TRUE)
  # Generating the PCA plot
  pca <- New_PCA(vstObj, var1, var2, labs = labs, ntop = ntop) + theme(legend.position = leg_pos)
  #Displaying the PCA plot
  print(pca)
  
  if (!(labs)){ # Labels == FALSE
    if (is.null(var2)){ # var2 == NULL
      # Saving PCA plot without batch correction, without variable 2, and without labels
      
      tiff(filename = file.path(outDir, paste0("PCA_Color_", var1, "_Top", ntop, 
                                               "Genes_Legend_", leg_pos, "_NoLabels.tiff")),
           res=300, units='in', width = wide, height=tall)
      
      print(pca)
      
      dev.off()
    } 
    else { # Batch == FALSE, Labels == FALSE, but var2 =/= NULL
      # Saving PCA plot without batch correction and without labels
      tiff(filename = file.path(outDir, paste0("PCA_Color_", var1, "_Shape_", var2, 
                                               "_Top", ntop, "Genes_Legend_", leg_pos, "_NoLabels.tiff")),
           res=300, units='in', width=wide, height=tall)
      
      print(pca)
      
      dev.off()
    }
  } 
  else { # Batch == FALSE, Labels == TRUE
    if (is.null(var2)){ # var2 == NULL
      # Saving PCA plot without batch correction and without variable 2
      tiff(filename = file.path(outDir, paste0("PCA_Color_", var1, "_Top", ntop, 
                                               "Genes_Legend_", leg_pos, ".tiff")),
           res=300, units='in', width=wide, height=tall)
      
      print(pca)
      
      dev.off()
    }
    else { # Batch == FALSE, Labels == TRUE, var2 =/= NULL
      # Saving PCA plot without batch correction
      tiff(filename = file.path(outDir, paste0("PCA_Color_", var1, "_Shape_", var2, 
                                               "_Top", ntop, "Genes_Legend_", leg_pos, ".tiff")),
           res=300, units='in', width=wide, height=tall)
      
      print(pca)
      
      dev.off()
    }
  }
}

# A function that creates and saves a JSON file that contains the details of the analysis
bulkJSON <- function(samples = NULL){
  
  etad <- format(Sys.time(), "%a %d %b %Y")
  emit <- format(Sys.time(), "%H:%M:%S %p %Z")
  verR <- paste0(R.Version()[6], ".", R.Version()[7])
  verDES <- packageDescription("DESeq2")$Version
  fastqs <- list(lapply(samples$Kallisto_Filename, function(x) {paste(as.character(x), sep="\n")}))
  filt <- paste("genes <=", rem, "total reads")
  comps <- list(lapply(names(comparisons), function(x) {paste(x, sep="\n")}))
  
  metrics <- c(etad, emit, verR, verDES, 
               fastqs, form, filt, FDR, 
               LFC, Sample_Group, comps, outDir)
  
  names(metrics) <- c("Date", "Time", "R Version", "DESeq2 Version",
                      "Samples", "Design", "Genes Filtered", "FDR", 
                      "LFC", "Sample_Group", "Comparisons", "Output Directory")
  
  threadr::write_json(metrics, 
                      file = file.path(outDir, "AnalysisInfo.json"), 
                      auto_unbox = T)
}

# ==========================================================================================================  
# =========================================== Main Code Runs Below =========================================
# ========================================================================================================== 
# Import variables from the command line
argList <- commandArgs(TRUE)
base.path <- as.character(argList[1])

# Import deduplicated list of mouse/human homologue genes
homologues <- read.table(file.path(base.path, "References/R/Human_Mouse_GeneHomologues_Deduplicated.tsv"), header=T)

# Setting metadata file directory
metaDir <- file.path(base.path, "Metadata/Chadarevian_2024_Leukoencephalopathy_Bulk_Metadata.tsv")

# The metadata file is read into a dataframe
samples_full <- read.table(metaDir, header = TRUE, stringsAsFactors = T)

# ==========================================================================================================
# ======================= Import Mouse data and filter for 1:1 homologue genes =============================
# ==========================================================================================================

# Filtering metadata to only contain the mouse samples
samples_mouse <- samples_full[(samples_full$Species %in% "Mouse"),]

# Extracts file path information for each abundance file from the metadata
files_Ms <- file.path(base.path, paste0("Kallisto_", samples_mouse$Genome), samples_mouse$Kallisto_Filename, "abundance.h5")
names(files_Ms) <- paste0(samples_mouse$Sample_Name)

# Verifying that all specified files exist
if (all(file.exists(files_Ms))==FALSE) {
  for (i in 1:length(files_Ms)){
    if (file.exists(files_Ms[i])==FALSE){
      print(paste0("File on line ", i, " is ", file.exists(files_Ms[i])))
    }
  }
  stop("Please ensure the metadata file matches the samples in analysis directory")
} else {
  print("All files exits. Moving on.")
}

# Import Kallisto abundance files and summarize to the gene level using TxImport
gene_list_Ms <- read.csv(file.path(base.path, "References/Kallisto", samples_mouse$Genome[1], paste0("tx2gene_", samples_mouse$Genome[1], ".tsv")), 
                         sep='\t', 
                         header=TRUE)
txi.gene_Ms <- tximport(files_Ms, 
                        type = "kallisto", 
                        txIn = TRUE, 
                        txOut = FALSE, 
                        tx2gene = gene_list_Ms, 
                        countsFromAbundance = 'no' )

# Remove all genes that are not included in the homologues table
Ms_Counts <- txi.gene_Ms$counts
Ms_Counts <- subset(Ms_Counts, rownames(Ms_Counts) %in% homologues$Mouse)

# Round the counts and convert to integers
Ms_Counts <- as.matrix(round(Ms_Counts))
mode(Ms_Counts) <- "integer"

# Create the DESeq2 object
form <- "~Mouse_Cell"
dds_gene_Ms <- DESeqDataSetFromMatrix(countData=Ms_Counts, 
                                      colData=samples_mouse, 
                                      design= eval(parse(text=form)))

# Filtering out genes with less than an average of 10 reads per sample
rem <- 10*length(colnames(counts(dds_gene_Ms)))
keep <- rowSums(counts(dds_gene_Ms)) >= rem
print(paste0("Keeping ", nrow(dds_gene_Ms[keep,]), " out of ", nrow(dds_gene_Ms), " genes"))
dds_gene_Ms <- dds_gene_Ms[keep,]
remove(keep)

# Creating output directory variables
analysisDate <- paste0(format(Sys.time(), "%Y%m%d"), "_Analysis")
outDir_merged <- file.path(base.path, "DESeq_Results", analysisDate)

# Create the output directories if they do not exist
if (!(dir.exists(outDir_merged))){
  if (!(dir.exists(file.path(base.path, "DESeq_Results")))){
    dir.create(file.path(base.path, "DESeq_Results"))
  }
  dir.create(file.path(base.path, "DESeq_Results", analysisDate))
}

# Creating a list of genes that were retained in the mouse object
keptGenes <- data.frame(Gene = rownames(dds_gene_Ms))

# ==========================================================================================================
# ============================ Import human data and change to mouse genes =================================
# ==========================================================================================================

# Filtering metadata to only contain the human samples
samples_human <- samples_full[(samples_full$Species %in% "Human"),]

# Extracts file path information for each abundance file from the metadata
files_hu <- file.path(base.path, paste0("Kallisto_", samples_human$Genome), samples_human$Kallisto_Filename, "abundance.h5")
names(files_hu) <- paste0(samples_human$Sample_Name)

# Verifying that all specified files exist
if (all(file.exists(files_hu))==FALSE) {
  for (i in 1:length(files_hu)){
    if (file.exists(files_hu[i])==FALSE){
      print(paste0("File on line ", i, " is ", file.exists(files_hu[i])))
    }
  }
  stop("Please ensure the metadata file matches the samples in analysis directory")
} else {
  print("All files exits. Moving on.")
}

# Import Kallisto abundance files and summarize to the gene level using TxImport
gene_list_hu <- read.csv(file.path(base.path, "References/Kallisto", samples_human$Genome[1], paste0("tx2gene_", samples_human$Genome[1], ".tsv")), 
                         sep='\t', 
                         header=TRUE)
txi.gene_hu <- tximport(files_hu, 
                        type = "kallisto", 
                        txIn = TRUE, 
                        txOut = FALSE, 
                        tx2gene=gene_list_hu, 
                        countsFromAbundance = 'no' )

# Remove all genes that are not included in the deduplicated human/mouse homologue file
Hu_Counts <- txi.gene_hu$counts
Hu_Counts_filt <- as.data.frame(subset(Hu_Counts, rownames(Hu_Counts) %in% homologues$Human))

# Replace human gene names with mouse gene names
Hu_Counts_filt$Human <- rownames(Hu_Counts_filt)
Hu_Counts_Ms <- merge(Hu_Counts_filt, homologues, by = "Human")
rownames(Hu_Counts_Ms) <- Hu_Counts_Ms$Mouse
Hu_Counts_Ms$Human <- NULL
Hu_Counts_Ms$Mouse <- NULL

# Round the counts and convert to integers
Hu_Counts_Ms <- as.matrix(round(Hu_Counts_Ms))
mode(Hu_Counts_Ms) <- "integer"

# Create the DESeq2 object
form <- "~1"
dds_gene_Hu <- DESeqDataSetFromMatrix(countData=Hu_Counts_Ms, 
                                      colData=samples_human, 
                                      design= eval(parse(text=form)))

# Adjusting kept genes to remove genes that are not in the human dataset
Hu_genes <- data.frame(Gene = rownames(dds_gene_Hu))
keptGenes_filt <- subset(keptGenes, keptGenes$Gene %in% Hu_genes$Gene)

# Writing the table of genes that will be used in this analysis
write.table(keptGenes_filt, 
            file.path(outDir_merged, paste0("RetainedGenes_", analysisDate, "_", rem, "CumulativeCounts", ".tsv")),
            quote=F, 
            sep='\t', 
            row.names=F)

# Filtering Human DESeq2 object by keptGenes_filt list
dds_gene_Hu <- dds_gene_Hu[keptGenes_filt$Gene,]

# Writing the raw human counts into a table
raw_hu <- as.data.frame(txi.gene_hu[['counts']])
raw_hu$Gene <- rownames(raw_hu)
raw_hu <- raw_hu[,c(ncol(raw_hu), 1:(ncol(raw_hu)-1))]
write.table(raw_hu, 
            file.path(outDir_merged, paste0("RawCounts_", analysisDate, "_HumanOnly.tsv")),
            quote = F, 
            sep = '\t', 
            row.names = F)

# ==========================================================================================================
# ================================== Performing VST on Mouse and Human Objects =============================
# ==========================================================================================================

# Filtering the mouse DESeq2 object by the keptGenes_filt list
dds_gene_Ms <- dds_gene_Ms[keptGenes_filt$Gene,]

# Normalizing the mouse data using DESeq2's Variance Stabilizing Transformation and saving the VST file
vsdB_gene_Ms <- varianceStabilizingTransformation(dds_gene_Ms, blind = TRUE)
VST_Ms <- as.data.frame(assay(vsdB_gene_Ms))
VST_Ms$Gene <- rownames(VST_Ms)
VST_Ms <- VST_Ms[,c(ncol(VST_Ms), (1:ncol(VST_Ms)-1))]
write.table(VST_Ms, 
            file.path(outDir_merged, paste0("VSTCounts_", analysisDate, "_MouseOnly.tsv")),
            quote = F, 
            sep = '\t', 
            row.names = F)

# Normalizing the human data using DESeq2's Variance Stabilizing Transformation and saving the VST file
vsdB_gene_Hu <- varianceStabilizingTransformation(dds_gene_Hu, blind = TRUE)
VST_Hu <- as.data.frame(assay(vsdB_gene_Hu))
VST_Hu$Gene <- rownames(VST_Hu)
VST_Hu <- VST_Hu[,c(ncol(VST_Hu), (1:ncol(VST_Hu)-1))]
write.table(VST_Hu, file.path(outDir_merged, paste0("VSTCounts_", analysisDate, "_HumanOnly.tsv")),
            quote = F, sep = '\t', row.names = F)

# Generating and saving a mouse/human merged VST file
VST_Merge <- merge(VST_Ms, VST_Hu, by = 'Gene')
VST_Merge <- VST_Merge[, c(1:10, 15:18, 11:14)]
write.table(VST_Merge, file.path(outDir_merged, paste0("VSTCounts_", analysisDate, "_MergedMouse+Human.tsv")),
            quote = F, sep = '\t', row.names = F)

# ==========================================================================================================
# ================================ DESeq2 Analysis on the Mouse Samples ====================================
# ==========================================================================================================

# Setting control group
dds_gene_Ms$Mouse_Cell <- relevel(dds_gene_Ms$Mouse_Cell, ref="hCSF1_PBS")

# Perform DESeq analysis
dds_gene_Ms <- DESeq(dds_gene_Ms)

# Set FDR and LFC cutoff levels
FDR <- 0.01
LFC <- 2

# Print the comparison names contained within dds_gene_Ms
resultsNames(dds_gene_Ms)

# Name the comparison and create the output directory variable
Sample_Group <- "Mouse_Samples_Only"
outDir <- file.path(base.path, "DESeq_Results", analysisDate, Sample_Group)

# List of the comparisons to include in the outputs
comparisons <- list("Mouse_Cell_hFIRE_HPC_Brain_vs_hCSF1_PBS",
                    "Mouse_Cell_hFIRE_PBS_vs_hCSF1_PBS")

# List of names for the comparison output directories
names(comparisons) <- list("hFIRE_HPC_BrainvshCSF1_PBS",
                           "hFIRE_PBSvshCSF1_PBS")

# Create the output directories if they do not exist
if (!(dir.exists(outDir))){
  if (!(dir.exists(file.path(base.path, "DESeq_Results", analysisDate)))){
    if (!(dir.exists(file.path(base.path, "DESeq_Results")))){
      dir.create(file.path(base.path, "DESeq_Results"))
    }
    dir.create(file.path(base.path, "DESeq_Results", analysisDate))
  }
  dir.create(file.path(base.path, "DESeq_Results", analysisDate, Sample_Group))
}

# Output and save the differential expression analysis values for each comparison
for (i in 1:length(comparisons)){
  
  # Create the output directory for the comparison and gene files
  if (!(dir.exists(file.path(outDir, names(comparisons)[i], "Genes")))){
    if (!(dir.exists(file.path(outDir, names(comparisons)[i])))){
      dir.create(file.path(outDir, names(comparisons)[i]))
    }
    dir.create(file.path(outDir, names(comparisons)[i], "Genes"))
  }
  if (!(dir.exists(file.path(outDir, names(comparisons)[i], "GSEA")))){
    dir.create(file.path(outDir, names(comparisons)[i], "GSEA"))
  }
  
  # Setting the output directory for the results files
  geneOut <- file.path(outDir, names(comparisons)[i], "Genes")
  
  print(paste0("Viewing results for: ", comparisons[[i]]))
  
  # Extract the results from the DESeq object depending on the type of comparison
  if (length(comparisons[[i]]) == 1) {
    res <- results(dds_gene_Ms, name = comparisons[[i]], alpha = FDR)
  } else if (length(comparisons[[i]]) == 2){
    res <- results(dds_gene_Ms, list( comparisons[[i]] ), alpha = FDR)
  }
  
  # Remove all NAs from results
  nonNA_res <- res[!is.na(res$padj),]
  # Apply an FDR cutoff to the results
  pcutoff_res <- nonNA_res[(nonNA_res$padj<FDR),]
  # Apply a log-fold change cutoff to the results
  folded_res <- pcutoff_res[(pcutoff_res$log2FoldChange >= LFC | pcutoff_res$log2FoldChange <= -LFC),]
  summary(pcutoff_res) 
  summary(folded_res)
  
  # Saving the resulsts file with all genes in the dataset 
  All_genes <- data.frame(nonNA_res)
  All_genes$Gene <- rownames(All_genes)
  All_genes <- All_genes[, c(7,1:6)]
  All_genes <- All_genes[rev(order(All_genes$log2FoldChange)),]
  # Replacing any padj=0 with a value one order of magnitude lower than the lowest padj
  # This is so the -log10(padj) will not return Inf
  lowP <- min(subset(All_genes$padj, -log10(All_genes$padj) < Inf))
  All_genes$padj[All_genes$padj == 0] <- lowP*0.1
  All_genes$log10 <- -log10(All_genes$padj)
  All_genes$signed <- (All_genes$log2FoldChange * All_genes$log10)/abs(All_genes$log2FoldChange)
  write.table(All_genes, 
              file=file.path(geneOut, paste0("AllGenes_", names(comparisons)[i], "_", Sample_Group, ".tsv")), 
              row.names = F, 
              quote=FALSE, 
              sep='\t')
  
  # Saving the gene list for use with GSEA
  GSEA_genes <- All_genes[,c(1,ncol(All_genes))]
  GSEA_genes <- GSEA_genes[order(GSEA_genes$signed, decreasing = T),]
  write.table(GSEA_genes, 
              file=file.path(outDir, names(comparisons)[i], "GSEA", paste0("GSEA_", names(comparisons)[i], "_", Sample_Group, ".rnk")), 
              row.names = F, col.names = F, 
              quote=FALSE, 
              sep='\t')
  
  # Saving the results file with genes that passed the FDR cutoff
  diff_genes <- data.frame(pcutoff_res)
  diff_genes$Gene <- rownames(diff_genes)
  diff_genes <- diff_genes[, c(7,1:6)]
  diff_genes <- diff_genes[rev(order(diff_genes$log2FoldChange)),]
  write.table(diff_genes, 
              file=file.path(geneOut, paste0("SigGenes_", names(comparisons)[i], "_FDR", FDR, "_", Sample_Group, ".tsv")),
              row.names = F, 
              quote=FALSE, 
              sep='\t')
  
  # Saving the genes that passed the FDR and LFC cutoffs
  folded_genes <- data.frame(folded_res)
  folded_genes$Gene <- rownames(folded_genes)
  folded_genes <- folded_genes[, c(7,1:6)]
  folded_genes <- folded_genes[rev(order(folded_genes$log2FoldChange)),]
  write.table(folded_genes, 
              file=file.path(geneOut, paste0("SigFoldGenes_", names(comparisons)[i], "_FDR", FDR, "_LFC", LFC, "_", Sample_Group, ".tsv")), 
              row.names = F, 
              quote=FALSE, 
              sep='\t')
  
}

# Writing the raw mouse counts into a table
raw <- as.data.frame(txi.gene_Ms[['counts']])
raw$Gene <- rownames(raw)
raw <- raw[,c(ncol(raw), 1:(ncol(raw)-1))]
write.table(raw, 
            file.path(outDir, paste0("RawCounts_", analysisDate, "_", Sample_Group, ".tsv")),
            quote = F, 
            sep = '\t', 
            row.names = F)

# ==========================================================================================================  
# =============================================== PCA Plots ================================================
# ==========================================================================================================

# Set the variables for the PCA plot
var1 <- "Mouse_Type" # Must include for PCA generation (sets color)
var2 <- "Cell_Type" # Make NULL if not using an additional variable (sets shape)
ntop <- 2000 # Number of genes to include in PCA
labs <- T # Make FALSE or F if not including labels
leg_pos <- "right" # Legend position. Use "none" to remove legend

# Create and save the PCA plot
savePCA(dds_gene_Ms)

# Saving a JSON file with the analysis details
bulkJSON(samples = samples_mouse)

# ==========================================================================================================  
# ========================================== Mouse Only Heatmaps ===========================================
# ========================================================================================================== 

# Creating the Heatmap output directory
if (!(dir.exists(file.path(outDir, "Heatmaps")))){
  dir.create(file.path(outDir, "Heatmaps"))
}

# Setting FDR and LFC cutoffs for the heatmap               
FDR <- FDR
LFC <- LFC

# Setting the heatmap cell height
height <- 5
width <- 10

# Creating and saving a median-centered heatmap for the "SigFold" genes from each comparison
for (j in 1:length(comparisons)){
  
  # Creating color palette for heatmap
  colfunc <- colorRampPalette(c("midnightblue", "dodgerblue4", "lightskyblue1", "#1FA187FF", "darkgreen", "darkslategray"))
  
  # Make the list of genes to filter the heatmap
  folded_genes <- Make_Folded_Genes(dds_gene_Ms, j, comparisons, FDR, LFC)
  
  # Filtering the data by the chosen gene list and assaying
  vsdB_filt_Ms <- vsdB_gene_Ms[(rownames(vsdB_gene_Ms) %in% folded_genes$Gene),]
  vsdB_filt_assay_Ms <- assay(vsdB_filt_Ms)
  
  # Median centering the filtered data
  vsdB_filt_assay_Ms <- t(as.data.frame(scale(t(vsdB_filt_assay_Ms), scale=F, center=T)))
  
  # Color key for heatmap columns
  metaCols <- c("Mouse_Type", "Cell_Type")
  annot <- as.data.frame(colData(dds_gene_Ms)[,metaCols])
  rownames(annot) <- colnames(vsdB_filt_assay_Ms)
  annotation_colors <- list(Mouse_Type=c(hCSF1="#0072B2", hFIRE="red"),
                            Cell_Type=c(PBS="darkblue", HPC="darkred"))

  # Making a pretty heatmap
  map <- pheatmap(vsdB_filt_assay_Ms,
                  annotation_col = annot,
                  annotation_colors = annotation_colors,
                  border_color = NA,
                  cellwidth = width,
                  cellheight = height,
                  cluster_cols = T,
                  clustering_distance_cols = "euclidean",
                  cluster_rows = T,
                  clustering_distance_rows = "euclidean",
                  show_rownames = T,
                  show_colnames = T,
                  fontsize = height,
                  legend = T,
                  treeheight_row = 0,
                  treeheight_col = 10,
                  color = colfunc(100))
  
  map$gtable$grobs[[1]]$gp <- gpar(lwd = 4)
  
  tiff(filename = file.path(outDir, "Heatmaps", 
                            paste0("Heatmap_SigFoldGenes_", names(comparisons)[j], "_FDR", FDR, "_LFC", LFC, "_", Sample_Group, "_Centered.tiff")), 
       height = (height/72 * nrow(vsdB_filt_assay_Ms))+2, 
       width = (width/72 * ncol(vsdB_filt_assay_Ms)) + 2, 
       units ='in', 
       res = 300, 
       compression = 'lzw')
  
  grid.newpage()
  grid.draw(map$gtable)
  dev.off()
}


# ==========================================================================================================  
# ======================================= Mouse Only Volcano Plots =========================================
# ========================================================================================================== 

# Creating the Volcano output directory
if (!(dir.exists(file.path(outDir, "Volcano")))){
  dir.create(file.path(outDir, "Volcano"))
}

# Setting FDR cutoff for volcano plot                
FDR <- FDR
LFC <- LFC
numLabs <- 10 # Number of genes to label

for (k in 1:length(comparisons)){
  # Make the list of genes to filter the heatmap
  # Extract the results from the DESeq object depending on the type of comparison
  if (length(comparisons[[k]]) == 1) {
    res <- results(dds_gene_Ms, name = comparisons[[k]], alpha = FDR)
  } else if (length(comparisons[[k]]) == 2){
    res <- results(dds_gene_Ms, list( comparisons[[k]] ), alpha = FDR)
  }
  
  # Remove all NAs from results
  nonNA_res <- res[!is.na(res$padj),]
  # Replacing any padj=0 with a value one order of magnitude lower than the lowest padj
  lowP <- min(subset(nonNA_res$padj, -log10(nonNA_res$padj) < Inf))
  nonNA_res$padj[nonNA_res$padj == 0] <- as.numeric(lowP*0.1)
  # Apply an FDR cutoff to the results
  pcutoff_res <- nonNA_res[(nonNA_res$padj<FDR),]
  # Apply a log-fold change cutoff to the results
  folded_res <- pcutoff_res[(pcutoff_res$log2FoldChange >= LFC | pcutoff_res$log2FoldChange <= -LFC),]
  folded_genes <- data.frame(folded_res)
  folded_genes$Gene <- rownames(folded_genes)
  folded_genes <- folded_genes[, c(7,1:6)]
  
  # Calculating the ideal limit values for the x and y axes
  logMin <- min(nonNA_res$log2FoldChange)
  logMax <- max(nonNA_res$log2FoldChange)
  pMax <- max(-log10(nonNA_res$padj))
  
  print(paste0("logMin: ", logMin))
  print(paste0("logMax: ", logMax))
  print(paste0("pMax: ", pMax))
  
  logRound <- 2
  pRound <- 10
  
  xMin <- (logRound*round((logMin-(logRound/2))/logRound))
  xMax <- (logRound*round((logMax+(logRound/2))/logRound))
  xDivide <- 2
  yMax <- (pRound*round((pMax+(pRound/2))/pRound))
  yDivide <- 10
  
  if (TRUE) {
    # Generating and saving the volcano plot
    mutateddf <- mutate(as.data.frame(nonNA_res), 
                        Cutoff=rownames(nonNA_res) %in% rownames(folded_res),
                        Color=ifelse(nonNA_res$log2FoldChange>0, '#AD0505', '#0C0C80'),
                        Combo=nonNA_res$log2FoldChange*-log10(nonNA_res$padj)) #Will have different colors depending on significance
    mutateddf <- within(mutateddf, Color[Cutoff == "FALSE"] <- 'gray80')
    rownames(mutateddf) <- rownames(nonNA_res)
    
    input <- cbind(gene=rownames(mutateddf), mutateddf) #convert the rownames to a column
    input <- input[order(input$Combo, decreasing = T),]
    
    geneLabels <- subset(input, input$Cutoff == "TRUE")
    geneUp <- subset(geneLabels, geneLabels$log2FoldChange > 0) # Keeping only positive LFC genes
    geneDown <- subset(geneLabels, geneLabels$log2FoldChange < 0) # Keeping only negative LFC genes
    
    # Setting the color scale based on which subsets of genes are present
    if (nrow(geneUp)== 0 && nrow(geneDown) == 0){
      cols <- c("gray80")
    } else if (nrow(geneUp) == 0 && nrow(geneDown) > 0){
      cols <- c("#0C0C80", "gray80")
    } else if (nrow(geneUp) >0 && nrow(geneDown) == 0) {
      cols <- c("#AD0505", "gray80") 
    } else {
      cols <- c("#0C0C80", "#AD0505", "gray80")
    }
    
    # Checking if there are enough positive LFC genes to label the number set
    # by numLabs. If not, using all genes in the list
    if (nrow(geneUp) >= numLabs) {
      geneUp <- head(geneUp, numLabs)
    }
    
    # Checking if there are enough negative LFC genes to label the number set
    # by numLabs. If not, using all genes in the list
    if (nrow(geneDown) >= numLabs){
      geneDown <- tail(geneDown, numLabs)
    }
    
    volc <- ggplot(input, aes(log2FoldChange, -log10(padj))) + #volcano plot with log2Foldchange versus pvalue
      theme_classic() +
      theme(axis.line = element_line(linewidth=1), 
            axis.text = element_text(color = "black", size = 10, face = "bold"),
            axis.ticks = element_line(color = "black", linewidth = 1),
            axis.title = element_text(face = 'bold', size = 10),
            plot.title = element_text(face = 'bold'),
            legend.position = "NA") +
      geom_point(data=input, aes(col=Color)) + #add points colored by significance
      scale_color_manual(values = cols) + 
      xlim(xMin, xMax) +
      ylim(0, yMax) +
      scale_x_continuous(limits = c(xMin, xMax), breaks = seq(xMin, xMax, by=xDivide), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, yMax), breaks = seq(0, yMax, by=yDivide), expand = c(0, 0)) +
      geom_hline(yintercept = -log10(FDR), linetype = "dashed", linewidth = 1, color = "gray30") +
      geom_vline(xintercept = c(LFC, -LFC), linetype = "dashed", linewidth = 1, color = "gray30") +
      geom_label_repel(data=rbind(geneUp, geneDown), color="white", size=2, segment.color = 'black', 
                       aes(label=gene, fill = alpha(c(Color), 0.7))) +
      scale_fill_identity()
    
    # Saving the volcano plot  
    tiff(filename = file.path(outDir, "Volcano", 
                              paste0("Volcano_", names(comparisons)[k], "_FDR", FDR, "_LFC", LFC, "_", Sample_Group, ".tiff")),
         res=300, 
         units='in', 
         width=7.2, 
         height=6, 
         compression = 'lzw')
    
    plot(volc)
    
    dev.off()
    
  }
}

# ==========================================================================================================  
# ======================================== Merged Species Heatmaps =========================================
# ========================================================================================================== 

# Creating the merged heatmap output directory
if (!(dir.exists(file.path(outDir_merged, "Merged_Heatmaps")))){
  dir.create(file.path(outDir_merged, "Merged_Heatmaps"))
}

# Import the previously saved VST file generated from the mixed species DESeq2 object
vsdB_gene_merged <- read.table(file.path(outDir_merged, paste0("VSTCounts_", analysisDate, "_MergedMouse+Human.tsv")), 
                        header=T)
rownames(vsdB_gene_merged) <- vsdB_gene_merged$Gene
vsdB_gene_merged$Gene <- NULL

# Setting FDR and LFC cutoffs for the heatmap               
FDR <- 0.01
LFC <- 2

# Setting the heatmap cell height
height <- 5
width <- 10

# Creating and saving a median-centered heatmap for the "SigFold" genes from each comparison
for (j in 1:length(comparisons)){
  
  # Creating color palette for heatmap
  colfunc <- colorRampPalette(c("midnightblue", "dodgerblue4", "lightskyblue1", "#1FA187FF", "darkgreen", "darkslategray"))
  
  # Make the list of genes to filter the heatmap
  folded_genes <- Make_Folded_Genes(dds_gene_Ms, j, comparisons, FDR, LFC)
  
  # Filtering the data by the chosen gene list and assaying
  vsdB_filt_assay <- vsdB_gene_merged[(rownames(vsdB_gene_merged) %in% folded_genes$Gene),]
  
  # Median centering the filtered data
  vsdB_filt_assay <- t(as.data.frame(scale(t(vsdB_filt_assay), scale=F, center=T)))
  
  # Color key for heatmap columns
  metaCols <- c("Mouse_Type", "Cell_Type")
  annot <- as.data.frame(samples_full[, metaCols])
  rownames(annot) <- colnames(vsdB_filt_assay)
  annotation_colors <- list(Mouse_Type=c(hCSF1="#0072B2", hFIRE="red"),
                            Cell_Type=c(PBS="darkblue", HPC="darkred"))
  
  # Make a pretty heatmap
  map <- pheatmap(vsdB_filt_assay,
                  annotation_col = annot,
                  annotation_colors = annotation_colors,
                  border_color = NA,
                  cellwidth = width,
                  cellheight = height,
                  cluster_cols = F,
                  cluster_rows = T,
                  clustering_distance_rows = "euclidean",
                  show_rownames = T,
                  show_colnames = T,
                  fontsize = height,
                  legend = T,
                  treeheight_row = 0,
                  treeheight_col = 10,
                  color = colfunc(100))
  
  map$gtable$grobs[[1]]$gp <- gpar(lwd = 4)
  
  tiff(filename = file.path(outDir_merged, "Merged_Heatmaps", 
                            paste0("Heatmap_MergedSpecies_SigFoldGenes_", names(comparisons)[j], "_FDR", FDR, "_LFC", LFC, "_", Sample_Group, "_Centered.tiff")), 
       height = (height/72 * nrow(vsdB_filt_assay))+2, 
       width = (width/72 * ncol(vsdB_filt_assay)) + 2, 
       units ='in', 
       res = 300, 
       compression = 'lzw')
  
  grid.newpage()
  grid.draw(map$gtable)
  dev.off()
}

# ===========================================================================================================  
# ============================================== Merged Box Plots ===========================================
# ===========================================================================================================

# Creating the merged boxplot output directory
if (!(dir.exists(file.path(outDir_merged, "Merged_BoxPlots")))){
  dir.create(file.path(outDir_merged, "Merged_BoxPlots"))
}

# Creating a list of genes to generate boxplots with
genes <- list(c("Csf1r", "Cx3cr1", "P2ry12", "C1qa", "Stab1", "Aif1"),
              c("Tyrobp", "Trem2", "Gpr34", "Olfml3", "Irf8", "Spi1"),
              c("Fcrls", "Il1a", "Ccr6", "Treml2", "Il7r", "Ly9"))
names(genes) <- c("xMG_Fully_Recovered_Top", 
                  "xMG_Fully_Recovered_Middle", 
                  "xMG_Unrecovered_Bottom")
              
# Generate boxplots for each set of genes in the list
for (i in names(genes)){
  
  # Subset the vsdB_gene file to only contain genes in the current list and melt
  dat <- vsdB_gene_merged[(rownames(vsdB_gene_merged) %in% genes[[i]]),]
  dat$Gene <- rownames(dat)
  dat_melt <- melt(as.data.table(dat),
                   id.vars = "Gene",
                   variable.name = "ID", 
                   value.name = c("Expression"))
  
  # Generate a metadata file for the boxplots
  samples_full$ID <- samples_full$Sample_Name
  samples_filtered <- samples_full[, c("ID", "Cell_Type", "Mouse_Type", "Mouse_Cell")]
  
  # Merge gene counts and sample metadata
  full_dat <- merge(dat_melt, samples_filtered, by="ID")
  
  # Set the factor levels to order the samples on the plot
  full_dat$Mouse_Cell <- factor(full_dat$Mouse_Cell, 
                                  levels = c("hCSF1_PBS", "hFIRE_PBS",
                                             "hFIRE_HPC_Microglia", "hFIRE_HPC_Brain"))
  
  full_dat$Gene <- factor(full_dat$Gene, levels = genes[[i]])
  
  # Generate the box plot
  bp <- ggplot(full_dat, aes(x=Gene, y=Expression, fill=Mouse_Cell)) +
    scale_fill_manual(values=c(hCSF1_PBS="#55a0fb", hFIRE_PBS="#ccf9e8", 
                               hFIRE_HPC_Microglia="darkgreen", hFIRE_HPC_Brain="burlywood2")) +
    geom_boxplot(outlier.shape = NA, lwd = 0.1) +
    geom_point(position=position_jitterdodge(), size = 0.5, aes(shape=Mouse_Cell)) +
    scale_shape_manual(values=c(hCSF1_PBS=15, hFIRE_PBS=16, 
                                hFIRE_HPC_Microglia=17, hFIRE_HPC_Brain=18)) +
    facet_wrap(~Gene, scales = 'free_x', ncol = 7) +
    theme_linedraw() +
    ylab("VST Normalized Expression") +
    xlab(NULL) +
    scale_x_discrete(expand = c(0,0)) +
    theme(axis.title = element_text(face = 'bold', size = 10), 
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "right",
          legend.text = element_text(size= 10),
          legend.title = element_blank(),
          strip.text = element_text(face = 'bold', size = 10))
  
  
  # Saving the box plots
  tiff(filename = file.path(outDir_merged, "Merged_BoxPlots",
                            paste0("BoxPlot_", i, "_MergedMouse+Human_Legend.tiff")),
       res=300, 
       units='in', 
       width=8, 
       height=1.33, 
       compression = 'lzw')
  
    plot(bp)
  
  dev.off()    

}




