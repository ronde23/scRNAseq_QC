################################## Readme File ##############################

Requirements:
1. Linux OS
2. Python version 3.x
3. R version >= 4.0.0

This single cell RNAseq quality control pipeline uses 3 tools - DropletQC for identifying empty droplets, SoupX to estimate ambient RNA in cells, Scrublet for identifying doublets.


To run the python scripts and generate output you need to run the following command:
If you want to use DropletQC, SoupX and Scrublet to filter empty droplets, ambient RNA, doublets followed by filtering of low quality cells (damaged cells) use the command below:
Run: Rscript QC_Automated.R --DropletQC 1 --SoupX 1 --Scrublet 1 --Gene 200 --RNA 2000 --Mt 10 --Mitopattern MT --Sample "Sample1/outs Sample2/outs Sample3/outs"

If you want to use DropletQC, SoupX and Scrublet to filter empty droplets, ambient RNA, doublets use the command below:
Run: Rscipt QC_Automated.R --DropletQC 1 --SoupX 1 --Scrublet 1 --Gene 0 --RNA 0 --Mt 0 --Mitopattern mt --Sample "Sample1/outs Sample2/outs Sample3/outs"

If you want to use DropletQCand Scrublet to filter empty droplets and doublets followed by filtering of low quality cells (damaged cells) use the command below:
Run: Rscript QC_Automated.R --DropletQC 1 --SoupX 0 --Scrublet 1 --Gene 200 --RNA 2000 --Mt 10 --Mitopattern MT --Sample "Sample1/outs Sample2/outs Sample3/outs"

If you want to use DropletQC and SoupX to filter empty droplets, ambient RNA, doublets followed by filtering of low quality cells (damaged cells) use the command below:
Run: Rscript QC_Automated.R --DropletQC 1 --SoupX 1 --Scrublet 0 --Gene 200 --RNA 2000 --Mt 10 --Mitopattern MT --Sample "Sample1/outs Sample2/outs Sample3/outs"

If you want to use SoupX and Scrublet to filter empty droplets, ambient RNA, doublets followed by filtering of low quality cells (damaged cells) use the command below:
Run: Rscript QC_Automated.R --DropletQC 0 --SoupX 1 --Scrublet 1 --Gene 200 --RNA 2000 --Mt 10 --Mitopattern MT --Sample "Sample1/outs Sample2/outs Sample3/outs"
