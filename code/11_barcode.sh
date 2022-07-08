### ================ ###
### 11 DNA barcoding ###
### ================ ###

### Preparations

# Remove old project directory

# Navigate to your location of choice (e.g. cd ~) and re-download repository
git clone https://github.com/mhelmkampf/meg_ss22

# On Windows, install and launch BioEdit
# zipped installation file provided in meg_ss22/apps
# or download from https://bioedit.software.informer.com/

# On Mac, launch 4Peaks
# executable provided in meg_ss22/apps



### Obtain barcode from Sanger trace files

# PCR / sequencing primers used for COI:
# M13-tailed primers (Ivanova et al. 2007, Molecluar Ecology Notes, lower case = tail)
# VF2_t1    tgtaaaacgacggccagtCAACCAACCACAAAGACATTGGCAC
# FishF2_t1 tgtaaaacgacggccagtCGACTAATCATAAAGATATCGGCAC
# FishR2_t1 caggaaacagctatgacACTTCAGGGTGACCGAAGAATCAGAA
# FR1d_t1   caggaaacagctatgacACCTCAGGGTGTCCGAARAAYCARAA
#> ~ 650 bp fragment in 5' region of COI

# View trace files with 4Peaks (Mac) or BioEdit (Windows)
# files are provided in meg_ss22/data/barcode
# each sample was sequenced in forward and reverse direction: *_F.ab1, *_R.ab1

# Trim low quality 5' and 3' positions from F and R file
# 4Peaks: select range of positions to keep, crop (Cmd + K) and export as plain text file
# BioEdit: select range of positions to keep, copy selection (Ctrl + C) to plain text file
# save as *_F_tr.fas, *_R_tr.fas

# Reverse complement R sequence
# http://www.reverse-complement.com/
# save as *_R_rc.fas

# Align sequences using Clustal Omega
# https://www.ebi.ac.uk/Tools/msa/clustalo/
# select DNA, output format Pearson/FASTA
# save as *.aln

# Create consensus sequence with Cons
# https://www.ebi.ac.uk/Tools/msa/emboss_cons/
# select DNA, change name to *_con
# save as *_con.fas

# Compile consensus sequences in one file
cd meg_ss22/data/barcode/fasta
cat *_con.fas > ../co1.fas



### Identify species of origin based on barcode

# Barcode of Life Data System (BOLD): http://www.boldsystems.org/index.php/
# select Identification, Animal Identification (COI), Species Level Barcode Records
# explore species, BIN and tree pages
# evaluate id quality based on similarity score, within-BIN and NN distances



### =============== Bonus material =============== ###

### Alternative to BOLD: identify barcode with BLAST and the NCBI nt database

# https://blast.ncbi.nlm.nih.gov/Blast.cgi
# select Nucleotide BLAST, Nucleotide collection (nr/nt), megablast



### Build phylogenetic tree of barcodes with MAFFT

# https://mafft.cbrc.jp/alignment/server/
# align consensus sequences with default settings
# select phylogenetic tree with NJ, Jukes-Cantor model, Bootstrap on settings
# view tree on Phylo.io



### ================== Solutions ================= ###

### Species IDs

# 10: Oncorhynchus keta                      100%
# 11: Gadus chalcogrammus                    100%
# 14: Platichthys flesus / Pleuronectes sp.  100%
# 15: Gadus chalcogrammus                    100%
# 19: Gadus chalcogrammus                    100%
# 20: Platichthys flesus / Pleuronectes sp.  100%
