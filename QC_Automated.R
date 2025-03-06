#!/usr/bin/env Rscript

library(optparse)
library(dplyr)
library(SoupX)
library(Seurat)
library(Matrix)
library(ggplot2)
library(stringr)
library(clustree)
library(DropletQC)
library(sctransform)
#library(singleCellTK)

# Run: Rscript QC_Automated.R --DropletQC 1 --SoupX 1 --Scrublet 1 --Gene 200 --RNA 2000 --Mt 10 --Mitopattern MT --Sample "10X-Lautz-H1-20240320-3v3-1DI-hg/outs 10X-Lautz-H3-20240320-3v3-1DI-hg/outs 10X-Lautz-L1-20240320-3v3-1DI-hg/outs 10X-Lautz-L3-20240320-3v3-1DI-hg/outs"
## Run: QC_Automated.R --DropletQC 1 --SoupX 1 --Scrublet 1 --Gene 0 --RNA 0 --Mt 0 --Mitopattern mt --Sample "/data/BCS/Services/June_Nakamura/Hydrocephalus-SnRNA-seq/alignment/10X-Goto-MUT-20231020-3v3-1DI-mm/outs /data/BCS/Services/June_Nakamura/Hydrocephalus-SnRNA-seq/alignment/10X-Goto-WT-20231020-3v3-1DI-mm/outs /data/BCS/Services/June_Nakamura/Hydrocephalus-SnRNA-seq/alignment/10X-Nakamura-June-228-WT-20240326-3v3-1DI-mm/outs /data/BCS/Services/June_Nakamura/Hydrocephalus-SnRNA-seq/alignment/10X-Nakamura-June-259-MUT-20240326-3v3-1DI-mm/outs"

opt_list = list(make_option(c("--DropletQC"), type = "integer", default = 1), make_option(c("--SoupX"), type = "integer", default = 1), make_option(c("--Scrublet"), type = "integer", default = 1), make_option(c("--Gene"), type = "integer", default = 0), make_option(c("--RNA"), type = "integer", default = 0), make_option(c("--Mt"), type = "integer", default = 0), make_option(c("--Mitopattern"), type = "character", default = "MT"), make_option(c("--Sample"), type = "character"))
opt_parser = OptionParser(option_list = opt_list)
opt = parse_args(opt_parser); #print(opt);

#args = commandArgs(trailingOnly=TRUE); #print(args);

if (opt$DropletQC == 0) {emptydrops = 0} else {emptydrops = 1} 
if (opt$SoupX == 0) {ambient_rna = 0} else {ambient_rna = 1}
if (opt$Scrublet == 0) {doublets = 0} else {doublets = 1}
if (opt$Gene == 0) {gene = 0} else {gene = opt$Gene} 
if (opt$RNA == 0) {rna = 0} else {rna = opt$RNA}
if (opt$Mt == 0) {mt = 100} else {mt = opt$Mt}
if (opt$Mitopattern == "") {mpattern = "MT"} else {mpattern = opt$Mitopattern}


	
samples = str_split(opt$Sample, " "); #print(samples[[1]]);
sample_no = length(samples[[1]]); #print(sample_no);
sample_files = c(); folders = c();
for (s in 1:sample_no) {
	sf = samples[[1]][s]  # 10X-Lautz-H1-20240320-3v3-1DI-hg/outs
	sample_files[s] = sf
	bgz = paste0("gzip -c -d ",gsub('/outs','',sample_files[s]),"/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > ",gsub('/outs','',sample_files[s]),"/barcodes.tsv"); #print(bgz);
	fgz = paste0("gzip -c -d ",gsub('/outs','',sample_files[s]),"/outs/filtered_feature_bc_matrix/features.tsv.gz > ",gsub('/outs','',sample_files[s]),"/features.tsv"); #print(fgz);
	mtx = paste0("gzip -c -d ",gsub('/outs','',sample_files[s]),"/outs/filtered_feature_bc_matrix/matrix.mtx.gz > ",gsub('/outs','',sample_files[s]),"/matrix.mtx"); #print(mtx);
	system(bgz, wait=TRUE)
	system(fgz, wait=TRUE)
	system(mtx, wait=TRUE)
}

pnames = function(dirname) { 
	pname = gsub('/outs','',dirname)
}

read_dirs = function(dirname) {
	rawfeat = Read10X(paste0(dirname, "/raw_feature_bc_matrix"))
	filteredfeat = Read10X(data.dir = paste0(dirname, "/filtered_feature_bc_matrix"))
	return(list(rawfeat, filteredfeat))
}

dropletqc = function(dirname) {
	#print("DropletQC started detection of empty dropltes....")
	write("DropletQC started detection of empty dropltes....", fileConn2, append=TRUE)
	mut_data_comp = Read10X(data.dir = paste0(dirname, "/filtered_feature_bc_matrix")); #print(ncol(mut_data_comp)); print(nrow(mut_data_comp)); print(head(colSums(mut_data_comp))); paste(cat("\n\n"));
	# Path or Location of Bam file
	bamfile = "possorted_genome_bam.bam"
	bam_location = paste0(dirname,"/",bamfile)
	# Detecting empty droplets using DropletQC
	nf = nuclear_fraction_tags(bam=bam_location, bam_index=paste0(bam_location,".bai"), barcodes=paste0(dirname,"/filtered_feature_bc_matrix/barcodes.tsv.gz"), tiles=1, cores=1, verbose=FALSE); #print(head(nf));
	nf_filtered = nf[colnames(mut_data_comp), ]
	nf_umis = cbind(nf_filtered, as.data.frame(colSums(mut_data_comp))); #print(head(nf_umis));
	empty_drops = identify_empty_drops(nf_umi=nf_umis); #print(paste("Cell count before removing Empty Droplets:", nrow(empty_drops))); #print(head(empty_drops)); paste(cat("\n"));
	write(paste("Cell count before removing Empty Droplets:", nrow(empty_drops)), fileConn1, append=TRUE)
	
	filtered_drops = subset(empty_drops, cell_status == "cell"); #print(paste("Empty Droplet Count:", nrow(subset(empty_drops, cell_status == "empty_droplet")))); print(paste("Cell count after removing Empty Droplets:", nrow(filtered_drops)));
	write(paste("Empty Droplet Count:", nrow(subset(empty_drops, cell_status == "empty_droplet"))), fileConn1, append=TRUE)
	write(paste("Cell count after removing Empty Droplets:", nrow(filtered_drops)), fileConn1, append=TRUE)
	
	if (doublets == 1) {
		# Cell indices of empty drops and cells
		cell_indices_keep = c(); cell_indices_remove = c();
		for (index in c(rownames(filtered_drops))) {
			cell_indices_keep = c(cell_indices_keep, which(index == c(rownames(empty_drops))))
		}
		cell_remove = setdiff(c(rownames(empty_drops)), c(rownames(filtered_drops))); #print(cell_remove);
		for (indice in cell_remove) {
			cell_indices_remove = c(cell_indices_remove, which(indice == c(rownames(empty_drops))))		
		}
		#print(paste("Cell indices to keep:")); print(cell_indices_keep); paste(cat("\n\n")); print(paste("Cell indices to remove:")); print(cell_indices_remove);

		# write modified filtered data in a file for Scrublet
		#print(data.frame(cell_indices_remove))
		write.table(data.frame(cell_indices_remove), file=paste0(pname,"EmptyDroplets.txt"), sep="\n", row.names = FALSE, col.names=FALSE, quote = FALSE)
	}
	
	#print("DropletQC successfully completed detection of empty dropltes.")
	write("DropletQC successfully completed detection of empty dropltes.", fileConn2, append=TRUE)
	# Save RDS File
	saveRDS(filtered_drops, file = paste0(gsub('/outs','',dirname),"_DropletFiltered.rds"))
	
	return(filtered_drops)
}

# SCTransform normalization and cluster identification of Seurat object for providing as Input to SoupX
norm_clusters = function(seurat_obj) {
	#print("Performing SCTransform normalization and cluster identification....")
	write("Performing SCTransform normalization and cluster identification....", fileConn2, append=TRUE)
	srat = SCTransform(seurat_obj, verbose = F)
	srat = RunPCA(srat, verbose = F)
	srat = RunUMAP(srat, dims = 1:40, verbose = F)
	srat = FindNeighbors(srat, dims = 1:40, verbose = F)
	srat = FindClusters(srat, verbose = T); #print(head(srat));
	
	#print("Successfully completed SCTransform normalization and cluster identification.")
	write("Successfully completed SCTransform normalization and cluster identification.", fileConn2, append=TRUE)
	
	return(srat)
}

soupx_inputs = function(dirname) {
	#print("Determining SoupX inputs...."); print(dirname);
	input_dirs = read_dirs(dirname); #print(dirname); 
	tod =  input_dirs[[1]]; #print("tod:"); print(nrow(tod)); print(ncol(tod));
	toc = input_dirs[[2]]; #print("toc:"); print(nrow(toc)); print(ncol(toc));
	if (emptydrops == 1) {
		filtered_drops = dropletqc(dirname)
		#tod = tod[, rownames(filtered_drops)]; print("tod:"); print(nrow(tod)); print(ncol(tod));
		tod = tod; #print("tod:"); print(nrow(tod)); print(ncol(tod));
		toc = toc[, rownames(filtered_drops)]; #print("toc:"); print(nrow(toc)); print(ncol(toc));
	} else {
		tod = tod
		toc = toc
	}
	#print("Successfully determined inputs for SoupX.")
	write("Successfully determined inputs for SoupX.", fileConn2, append=TRUE)
	
	return(list(tod, toc))	
}

soupx = function(dirname) {
	##soup_input = soupx_inputs(dirname)
	##tod = soup_input[[1]]
	##toc = soup_input[[2]]
	
	input_dirs = read_dirs(dirname); #print(dirname); 
	tod =  input_dirs[[1]]; #print("tod:"); print(nrow(tod)); print(ncol(tod));
	toc = input_dirs[[2]]; #print("toc:"); print(nrow(toc)); print(ncol(toc));
	
	if (emptydrops == 1) {
		filtered_drops = dropletqc(dirname)
		#tod = tod[, rownames(filtered_drops)]; print("tod:"); print(nrow(tod)); print(ncol(tod));
		tod = tod; #print("tod:"); print(nrow(tod)); print(ncol(tod));
		toc = toc[, rownames(filtered_drops)]; #print("toc:"); print(nrow(toc)); print(ncol(toc));
	} else {
		tod = tod
		toc = toc
	}
	
	# Soup Channel
	sc = SoupChannel(tod, toc, calcSoupProfile = FALSE)
	# Calculate soup profile
	soupProf = data.frame(row.names = rownames(toc), est = rowSums(toc)/sum(toc), counts = rowSums(toc))
	sc = setSoupProfile(sc, soupProf)
	
	# Seurat objects before any filtering step
	input_dirs = read_dirs(dirname)
	fil_data = input_dirs[[2]]; pname = pnames(dirname);
	sc_mut_comp_obj = CreateSeuratObject(counts = fil_data, project=pname); #print(ncol(sc_mut_comp_obj)); # , min.cells = 3, min.features = 200
	# Mitochondrial contamination
	sc_mut_comp_obj[["percent.mt"]] = PercentageFeatureSet(sc_mut_comp_obj, pattern = paste0("^",mpattern,"-")); #print(sc_mut_cntrl[["percent.mt"]]); # For mouse: mt, Human: MT
	
	if (emptydrops == 1) {
		##filtered_drops = dropletqc(dirname)
		sc_mut_comp_obj = sc_mut_comp_obj[, rownames(filtered_drops)]
	}
	
	nclusters = norm_clusters(sc_mut_comp_obj)
	meta = nclusters@meta.data; #print(head(meta));
	umap = nclusters@reductions$umap@cell.embeddings; #print("UMAP:"); print(head(umap)); 
	
	# Set seurat clusters as SoupX cell clusters
	sc = setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
	#print(head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20))
	# Set dimensionality reduction info
	sc = setDR(sc, umap)
	#print("SoupX started detection of ambient RNA....")
	write("SoupX started detection of ambient RNA....", fileConn2, append=TRUE)
	
	# Estimate rho
	sc = autoEstCont(sc); #print(sc$fit$rhoEst[1]);
	write(paste0("global rho: ",sc$fit$rhoEst[1]), fileConn1, append=TRUE)
	# Clean the data
	scout = adjustCounts(sc); #print(head(scout));	
	
	#print("SoupX successfully completed detection of ambient RNA.")
	write("SoupX successfully completed detection of ambient RNA.", fileConn2, append=TRUE)
	
	# Save RDS File
	saveRDS(scout, file = paste0(gsub('/outs','',dirname),"_AmbientRNAFiltered.rds"))
	
	return(scout)
}

scrublets = function(dirname, scrublet_file) {
	pname = pnames(dirname)
	if (ambient_rna == 1) { 
		soup_out = soupx(dirname) 
		# Seurat objects
		sc_mut_comp = CreateSeuratObject(counts = soup_out, project=pname); #print(head(sc_mut_comp)); print(ncol(sc_mut_comp)); 
		# Mitochondrial contamination
		sc_mut_comp[["percent.mt"]] = PercentageFeatureSet(sc_mut_comp, pattern = paste0("^",mpattern,"-")); # For mouse: mt, Human: MT
	} else {
		input_dirs = read_dirs(dirname)
		fil_data = input_dirs[[2]]
		if (emptydrops == 1) {
			filtered_drops = dropletqc(dirname)
			fil_data = fil_data[, rownames(filtered_drops)]	
		}
		# Seurat objects
		sc_mut_comp = CreateSeuratObject(counts = fil_data, project=pname); #print(head(sc_mut_comp)); print(ncol(sc_mut_comp)); 
		# Mitochondrial contamination
		sc_mut_comp[["percent.mt"]] = PercentageFeatureSet(sc_mut_comp, pattern = "^MT-"); # For mouse: mt, Human: MT	
	}
	
	#print("Reading output from Srublet...")
	write("Reading output from Srublet...", fileConn2, append=TRUE)
	
	# Run python script from R
	system(paste('scrublet_doublet.py', gsub('/outs','',dirname)), wait = TRUE)
	
	# Read Scrublet output files that contains doublets
	scrublet_data = read.csv(scrublet_file, sep="\t", header=TRUE); #print(nrow(scrublet_data));
	doublet_data = subset(scrublet_data, Doublet=="True"); #print(paste("Doublet Count:", nrow(doublet_data))); 
	filtered_data = subset(scrublet_data, Doublet=="False"); #print(paste("Cell count after removing Doublets:", nrow(filtered_data))); #print(head(filtered_data)); print(head(filtered_data$Celltag));
	write(paste("Doublet Count:", nrow(doublet_data)), fileConn1, append=TRUE)
	write(paste("Cell count after removing Doublets:", nrow(filtered_data)), fileConn1, append=TRUE)
	
	# Removing Doublets
	sc_mut_comp_doublet_filtered =  sc_mut_comp[, filtered_data$Celltag]; #print(head(sc_mut_comp_doublet_filtered)); #print(ncol(sc_mut_comp_doublet_filtered)); 
	
	#print("Completed analyzing output from Srublet.")
	write("Completed analyzing output from Srublet.", fileConn2, append=TRUE)
	
	# Save RDS File
	saveRDS(sc_mut_comp_doublet_filtered, file = paste0(gsub('/outs','',dirname),"_DoubletFiltered.rds"))
	
	return(sc_mut_comp_doublet_filtered)
}


#sample_files = c("10X-Lautz-H1-20240320-3v3-1DI-hg/outs", "10X-Lautz-H3-20240320-3v3-1DI-hg/outs", "10X-Lautz-L1-20240320-3v3-1DI-hg/outs", "10X-Lautz-L3-20240320-3v3-1DI-hg/outs", "10X-Lautz-cardiac-spheroid-high-20240216-3v3-1DI-hg/outs")
print(sample_files)
for (x in 1:length(sample_files)) {
	fileConn1 = paste0(gsub('/outs','',sample_files[x]),"_QC_Log.out")
	fileConn2 = paste0(gsub('/outs','',sample_files[x]),"_QC_Log.err")
	#Check its existence
	if (file.exists(fileConn1)) {
	  #Delete file if it exists
	  file.remove(fileConn1)
	}
	if (file.exists(fileConn2)) {
	  #Delete file if it exists
	  file.remove(fileConn2)
	}
	#print(sample_files[x])
	#if (doublets == 1) {
	#	#print("Running Scrublet for doublet detection....")
	#	# Run python script from R
	#	system(paste('python3 scrublet_doublet.py', gsub('/outs','',sample_files[x])), wait = FALSE)
	#	#print("Scrublet completed doublet detection.")
	#}

	pname = pnames(sample_files[x]);
	if (doublets == 1) {
		#scrublet_files = c("10X-Lautz-H1-20240320-3v3-1DI-hgScrublet_Doublets_Scrublet.csv", "10X-Lautz-H3-20240320-3v3-1DI-hgScrublet_Doublets_Scrublet.csv", "10X-Lautz-L1-20240320-3v3-1DI-hgScrublet_Doublets_Scrublet.csv", "10X-Lautz-L3-20240320-3v3-1DI-hgScrublet_Doublets_Scrublet.csv", "10X-Lautz-cardiac-spheroid-high-20240216-3v3-1DI-hgScrublet_Doublets_Scrublet.csv")
		scrublet_files = paste0(gsub("/outs","",sample_files[x]),"Scrublet_Doublets_Scrublet.csv")
		sc_mut_comp_data = scrublets(sample_files[x], scrublet_files)
	} else {
		if (ambient_rna == 1) { 
			soup_out = soupx(sample_files[x]) 
			# Seurat objects
			sc_mut_comp_data = CreateSeuratObject(counts = soup_out, project=pname); #print(head(sc_mut_comp_data)); print(ncol(sc_mut_comp_data)); 
			# Mitochondrial contamination
			sc_mut_comp_data[["percent.mt"]] = PercentageFeatureSet(sc_mut_comp_data, pattern = "^MT-"); # For mouse: mt, Human: MT
		} else {
			input_dirs = read_dirs(sample_files[x])
			fil_data = input_dirs[[2]]
			if (emptydrops == 1) {
				filtered_drops = dropletqc(sample_files[x])
				fil_data = fil_data[, rownames(filtered_drops)]	
			}
			# Seurat objects
			sc_mut_comp_data = CreateSeuratObject(counts = fil_data, project=pname); #print(head(sc_mut_comp_data)); print(ncol(sc_mut_comp_data)); 
			# Mitochondrial contamination
			sc_mut_comp_data[["percent.mt"]] = PercentageFeatureSet(sc_mut_comp_data, pattern = "^",mpattern,"-"); # For mouse: mt, Human: MT	
		}
	}

	# Removing low quality cells
	sc_mut_comp_data = subset(sc_mut_comp_data, subset = nFeature_RNA > gene & nCount_RNA > rna & percent.mt < mt); #print(ncol(sc_mut_comp_data)); # ncol - gives the no of cells

	# Save RDS File
	saveRDS(sc_mut_comp_data, file = paste0(gsub('/outs','',sample_files[x]),".rds"))
}



