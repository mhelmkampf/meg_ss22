### =============================== ###
### 07 SNPs and population genomics ###
### =============================== ###


### Genome assembly QC

# Quality assessment with QUAST
# Web tool: http://cab.cc.spbu.ru/quast/
# Report for H. puella reference genome: http://cab.cc.spbu.ru/quast/reports/02_Jun_2022_10:19:29_143981/report.html
# Report for example data: http://cab.cc.spbu.ru/quast/reports/02_Jun_2022_10:18:17_569501/report.html


### Variant Call Format (in bash)

## Print (compressed) VCF file to screen
cd ~/meg_ss22/data/genome   # if project directory is in home (~), adjust if necessary
zcat < snps_lg12_phased_36.vcf.gz | head
zcat < snps_rad.vcf.gz | head


## Demo code to be run on the UOL HPC cluster: subset VCF
vcftools --gzvcf ../chapter2_phased_mac2.vcf.gz \
  --keep meg_36.ids \
  --chr LG12 \
  --mac 4 --max-missing 0.33 --remove-indels \
  --recode --stdout | gzip > snps_lg12_phased_36.vcf.gz


### Handling VCF files in R

## Preparations
install.packages("vcfR")

library(tidyverse)
library(vcfR)
library(adegenet)

# Set working directory to [...]/meg_ss22


## Read in VCF
lg12_vcf <- read.vcfR("data/snps/snps_lg12_phased_36.vcf.gz")
lg12_vcf

rad_vcf <- read.vcfR("data/snps/snps_rad.vcf.gz")
rad_vcf


## View VCF sections
lg12_vcf@meta
rad_vcf@meta

head(lg12_vcf@fix)
head(rad_vcf@fix)

head(lg12_vcf@gt)
head(rad_vcf@gt)

rad_vcf@gt[1:20, 1:6]   # [row, column] indices can be used to subset data


## Convert to genlight object
lg12_gl <- vcfR2genlight(lg12_vcf)
class(lg12_gl)


# Plot genlight object
glPlot(lg12_gl)   # warning: computationally intense


## Distribution of allele frequencies
freq_lg12 <- glMean(lg12_gl)   # compute mean of alternate alleles

hist(freq_lg12, breaks = 10, proba = TRUE, col = "grey",   # Shortcut for a quick histogram
     xlab = "Allele frequencies",
     main = "Distribution of alternate allele frequencies")


## Principal Component Analysis
pca_lg12 <- glPca(lg12_gl, nf = 2)
pca_lg12
pca_lg12$scores

scores_lg12 <- as.data.frame(pca_lg12$scores) %>%
  rownames_to_column("Sample") %>%
  as_tibble()

# If this step stalls, load with: scores_lg12 <- read_csv("data/tables/pca_lg12_scores.csv")


## Exercise: Plot PC scores


## Add species and location information, improve plot
scores_lg12SL <- scores_lg12 %>%
  mutate(Species = str_sub(Sample, -6, -4),
         Location = str_sub(Sample, -3, -1))

l <- ggplot(data = scores_lg12SL, aes(x = PC1, y = PC2, color = Species, shape = Location)) +
  geom_point(size = 4, alpha = 0.75) +
  labs(title = "LG12") +
  scale_color_manual(values = c("goldenrod2", "royalblue3", "grey20", "coral2", "grey80")) +
  theme_light(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.x = element_text(vjust = -1.5))
l


## Zooming into the plot
l + coord_fixed(xlim = c(-20, 20), ylim = c(-40, 0))


## Plot eigenvalues
pca_lg12$eig

var_lg12 <- pca_lg12$eig / sum(pca_lg12$eig)

barplot(var_lg12, main = "Proportion of variance explained", las = 2)   # Another shortcut for a quick barplot



### ================== Solutions ================= ###

## Plot PC scores (minimal solution)
ggplot(data = scores_lg12, aes(x = PC1, y = PC2)) +
  geom_point()



### =============== Bonus material =============== ###
  
### Analyze RADseq example data

## Convert to genlight object (adegenet from here on)
ids <- getID(rad_vcf)
length(unique(myID, incomparables = NA)) == length(ids)
rad_tmp <- rad_vcf[!duplicated(ids, incomparables = NA), ]   # hack to fix duplicate ids in file

rad_gl <- vcfR2genlight(rad_tmp)


## Plot genlight object
glPlot(rad_gl)


## Distribution of allele frequencies
freq_rad <- glMean(rad_gl)   # compute mean of alternate alleles

hist(freq_rad, breaks = 10, proba = TRUE, col = "grey", 
     xlab = "Allele frequencies",
     main = "Distribution of alternate allele frequencies")


## Principal Component Analysis
pca_rad <- glPca(rad_gl, nf = 2)
pca_rad$scores

scores_rad <- as.data.frame(pca_rad$scores) %>%
  rownames_to_column("Sample") %>%
  as_tibble()


## Plot PC scores (minimal solution)
ggplot(data = scores_rad, aes(x = PC1, y = PC2)) +
  geom_point()


## Add species and location information, improve plot
scores_radSL <- scores_rad %>%
  mutate(Species = str_sub(Sample, -8, -6),
         Location = str_sub(Sample, -5, -3))

r <- ggplot(data = scores_radSL, aes(x = PC1, y = PC2, color = Species, shape = Location)) +
  geom_point(size = 4, alpha = 0.75) +
  labs(title = "RADseq") +
  scale_color_manual(values = c("grey20", "coral2", "grey80")) +
  theme_light(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.x = element_text(vjust = -1.5))
r


## Plot eigenvalues
pca_rad$eig

var_rad <- pca_rad$eig / sum(pca_rad$eig)

barplot(var_rad, main = "Proportion of variance explained", las = 2)   # Another shortcut for a quick barplot
