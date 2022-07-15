### ================ ###
### 12 Metabarcoding ###
### ================ ###


### Preparations

# Change into project directory
cd ~/meg_ss22   # change to path (directory location) on your system

# Update Git repository
git pull

# On error, remove old project directory and download again
# git clone https://github.com/mhelmkampf/meg_ss22



### QIIME 2

# All data files are found in meg_ss22/data/meta
# Drag and drop .qzv files at view.qiime2.org to view

## Demultiplexing

# 01_crassa-demux-primersremoved.qzv
# 02_crassa-demux-joined.qzv

## Denoising

# 03_crassa-deblur-stats.qzv
# 04_crassa-deblur-table.qzv
# 05_crassa-deblur-seqs.qzv

## Taxonomic classification

# 06_crassa-tax-all-bar.qzv
# 07_crasse-tax-cm-bar.qzv

## Diversity analyses

# 08_crassa-alpha-rarefaction.qzv
# 09_crassa-evenness-sign.qzv
# 10_bray-curtis-pcoa.qzv
# 11_jaccard-pcoa.qzv
