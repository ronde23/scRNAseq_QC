#!/usr/bin/env python

import os
import re
import sys
import glob
import scipy.io
import pathlib
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scrublet as scr

cmdlist = sys.argv; #print(sys.argv); # /data/BCS/Services/Marie_Filippi/ScRNAseq-021924/alignment/10X-Filppi-Keima-High-20230207-3v3-1DI-mm/outs/raw_feature_bc_matrix
# /data/BCS/Services/Marie_Filippi/ScRNAseq-021924/alignment/AS_SLAM_HHT/outs/raw_feature_bc_matrix

def DoubletsResults(fname):
    outfile = open(((fname.split('/'))[-1]).replace('.csv','_Scrublet.csv'), 'w')
    outfile.write('Celltag'+'\t'+'Doublet'+'\n') #+fname.replace('.csv', '')
    with open(fname, 'r') as infile: # fname = "Score.csv"
        lines = infile.readlines()
        for line in lines[:]:
            split_line = (line.rstrip()).split('\t')
            celltag = re.sub('[^a-zA-Z0-9-]','',split_line[0]); #print(celltag);
            outfile.write(celltag+'\t'+split_line[1]+'\n')
    outfile.close()
    # Deleting file
    path = pathlib.Path(fname)
    path.unlink()

empty_droplet_files = glob.glob('*EmptyDroplets.txt'); #print(empty_droplet_files);

for i in range(len(cmdlist)-1):
    input_dir = cmdlist[i+1] # '/data/BCS/Services/Marie_Filippi/ScRNAseq-021924/alignment/10X-Filppi-Keima-Low-20230207-3v3-1DI-mm/outs/raw_feature_bc_matrix'
    
    if (len(empty_droplet_files) == 0): # DropletQC has not been run
        barcodelist = []

    elif (len(empty_droplet_files) > 0):
        pname = input_dir; barcode_fname = pname+"EmptyDroplets.txt";
        
        barcodefl = open(barcode_fname, 'r')
        lines = barcodefl.readlines()
        barcodelist = []
        for line in lines:
            barcodelist.append(int(line.rstrip())-1)
        barcodefl.close()
    
    counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx').T.tocsc()  # tocsc()
    #print(counts_matrix); print(counts_matrix.shape); 
    
    data = scipy.sparse.csr_matrix(np.delete(counts_matrix.toarray(),barcodelist, axis=0)); #print(data.shape); print(data);
    
    genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1)) # numpy array
    features = ((np.array(pd.read_csv(input_dir + '/features.tsv', delimiter='\t', header=None)))[:,1]); # Extracting column 1
    #features = features[:, None] # Converting array of form (n,) to (n,1)
    #features_trans = features.T # Taking Transpose of a matrix
    cells = np.array(pd.read_csv(input_dir + '/barcodes.tsv', delimiter='\t', header=None)); 
   # print(cells); print(cells.shape); # print(features_trans); print(features_trans.shape); 
    #print(cells[1])
    #for x in range(len(barcodelist)):
    cells_data = np.delete(cells, barcodelist)
    #print(cells_data.shape); print(cells_data);
    
    if ("raw" in input_dir):
        outname = input_dir.replace("/outs/raw_feature_bc_matrix", "")
        #dirname_df = input_dir.replace("/outs/raw_feature_bc_matrix", "") + 'Scrublet_Data.csv' 
        dirname_scores = input_dir.replace("/outs/raw_feature_bc_matrix", "") + 'Scrublet_Scores.csv' 
        dirnames_doublets = input_dir.replace("/outs/raw_feature_bc_matrix", "") + 'Scrublet_Doublets.csv'
    elif("filtered" in input_dir):
        outname = input_dir.replace("/outs/filtered_feature_bc_matrix", "")
        #dirname_df = input_dir.replace("/outs/filtered_feature_bc_matrix", "") + 'Scrublet_Data.csv' 
        dirname_scores = input_dir.replace("/outs/filtered_feature_bc_matrix", "") + 'Scrublet_Scores.csv' 
        dirnames_doublets = input_dir.replace("/outs/filtered_feature_bc_matrix", "") + 'Scrublet_Doublets.csv'
    
    # You need to delete these following 3 lines if your matrix and tsv file are in raw or filtered folders.    
    dirname_scores = input_dir.replace("/outs/filtered_feature_bc_matrix", "") + 'Scrublet_Scores.csv' 
    dirnames_doublets = input_dir.replace("/outs/filtered_feature_bc_matrix", "") + 'Scrublet_Doublets.csv'
    #print(dirname_scores, dirnames_doublets)
    
    #print('Counts matrix shape: {} rows, {} columns'.format(data.shape[0], data.shape[1]))
    #print('Number of genes in gene list: {}'.format(len(genes)))
    scrub = scr.Scrublet(data, expected_doublet_rate=0.05) #### counts_matrix
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
    #print(doublet_scores, predicted_doublets); print(dirname_scores, dirnames_doublets);
    
    df_scores = pd.DataFrame(doublet_scores, index=cells_data)
    df_doublets = pd.DataFrame(predicted_doublets, index=cells_data); #print(df_doublets.head()); print(df_doublets.columns);
    df_doublets.columns = ['Doublet']; #print(df_doublets.head());
    df_scores.to_csv(dirname_scores, sep = '\t', header = False)
    df_doublets.to_csv(dirnames_doublets, sep = '\t', header = False)
    # Writing formatted results to files    
    DoubletsResults(dirname_scores)
    DoubletsResults(dirnames_doublets)
    
    #scrub_doublets = scrub.call_doublets(threshold=0.25)
    outname = input_dir.replace("/outs/filtered_feature_bc_matrix", "")
    scrub.plot_histogram()
    plt.savefig(outname+'doublet_score_histogram.png')
    #print('Running UMAP...')
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    #print('Done.')
    
    #print((pd.DataFrame(doublet_scores)).value_counts()); print((pd.DataFrame(predicted_doublets)).value_counts());
    
    scrub.plot_embedding('UMAP', order_points=True);
    plt.savefig(outname+'UMAP.png')
    




