### ================ ###
### 11 DNA barcoding ###
### ================ ###

### Preparations

# Remove old project directory

# Navigate to your location of choice (e.g. cd ~) and re-download repository
git clone https://github.com/mhelmkampf/meg_ss22

# On Windows, install BioEdit
# zipped installation file provided in meg_ss22/apps
# or download from https://bioedit.software.informer.com



### Obtain barcode from Sanger trace files

# PCR / sequencing primers used for COI:
# M13-tailed primers (Ivanova et al. 2007, Molecluar Ecology Notes, lower case = tail)
# VF2_t1    tgtaaaacgacggccagtCAACCAACCACAAAGACATTGGCAC
# FishF2_t1 tgtaaaacgacggccagtCGACTAATCATAAAGATATCGGCAC
# FishR2_t1 caggaaacagctatgacACTTCAGGGTGACCGAAGAATCAGAA
# FR1d_t1   caggaaacagctatgacACCTCAGGGTGTCCGAARAAYCARAA
#> ~ 650 bp fragment in 5' region of COI

# View trace files with 4Peaks (Mac) or BioEdit (Windows)
# found in meg_ss22/apps
# each sample has been sequenced in forward and reverse direction: *_F.ab1, *_R.ab1

# Trim low quality 5' and 3' positions from F and R file
# 4Peaks: select range of positions to keep, crop
# BioEdit:

# Export trimmed sequence to fasta file
#> *_F_tr.fas, *_R_tr.fas

# Reverse complement R sequence at http://www.reverse-complement.com
#> *_R_rc.fas

# Align sequences using Clustal Omega: https://www.ebi.ac.uk/Tools/msa/clustalo/
# settings: DNA, output format ClustalW with character counts
# repeat with settings output format Pearson/FASTA
#> *.aln

# Create consensus sequence with Cons: https://www.ebi.ac.uk/Tools/msa/emboss_cons/
# settings: DNA, change name to *_con
#> *_con.fas

# Compile consensus sequences in one file
cd meg_ss22/data/barcode/fasta
cat *_con.fas > ../co1.fas



### Identify species of origin based on barcode

# BOLD Identification System: http://www.boldsystems.org/index.php/IDS_OpenIdEngine
# set to Animal Identification (COI), Species Level Barcode Records
# explore species, BIN and tree pages
#? evaluate id quality based on similarity score, within-BIN and NN distances



### =============== Bonus material =============== ###

### BLAST

# https://blast.ncbi.nlm.nih.gov/Blast.cgi
