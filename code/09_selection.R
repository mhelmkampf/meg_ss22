### ============ ###
### 09 Selection ###
### ============ ###


### Preparations
instlg12.packages("vcfR")

library(tidyverse)
library(vcfR)
library(adegenet)

# Set working directory to [...]/meg_ss22



### Exercise 1: Estimate Ne from LD with NeEstimator

## Software available in /software/NeEstimator
## To run, execute in terminal Ne2-1.exe on Windows, Ne2M on Mac, or Ne2L on Linux
## Provide path to input file when prompted, e.g. ../../data/puella_caribbean.gen)
## For method choose linkage disequilibrium (1)
## Results will be saved in /other/software/NeEstimator



### Exercise 2: Simulating the effect of selection on allele frequencies

## Go to: https://www.biologysimulations.com/population-genetics
## Set population size to 500, number of generations to 100, and frequency of red allele to 0.05
## Choose a combination of parameters regarding Natural Selection for the following scenarios:
# - Selection for a dominant phenotype (Red)
# - Selection against a recessive phenotype (Red)
# - Heterozygote advantage
## Run the simulation a few times and observe



### Exercise 3: Identifying genomic regions potentially under selection

## Fst calculations were done with vcftools, output and log files are in data/selection)

# Calculate joint Fst per SNP
# vcftools --gzvcf snps_lg12_phased_36.vcf.gz \
# --weir-fst-pop pop.gumhon.txt \
# --weir-fst-pop pop.indbel.txt \
# --weir-fst-pop pop.nigbel.txt \
# --weir-fst-pop pop.nighon.txt \
# --weir-fst-pop pop.puebel.txt \
# --weir-fst-pop pop.unibel.txt \
# --stdout 1> fst_snp_lg12.tsv 2> fst_snp_lg12.log

# Calculate joint Fst per 50kb window
# vcftools --gzvcf snps_lg12_phased_36.vcf.gz \
# --weir-fst-pop pop.gumhon.txt \
# --weir-fst-pop pop.indbel.txt \
# --weir-fst-pop pop.nigbel.txt \
# --weir-fst-pop pop.nighon.txt \
# --weir-fst-pop pop.puebel.txt \
# --weir-fst-pop pop.unibel.txt \
# --fst-window-step 5000 \
# --fst-window-size 50000 \
# --stdout 1> fst_50kb_lg12.tsv 2> fst_50kb_lg12.log


## Plot joint Fst, per SNP
snp <- read_tsv("data/selection/fst_snp_lg12.tsv")

(p <- ggplot(data = snp, aes(x = POS, y = WEIR_AND_COCKERHAM_FST)) +
  geom_point(size = 0.25, alpha = 0.5) +
  labs(x = "Position", y = "Fst") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
)


## Add chromosome-wide Fst (value found in log files)
p + geom_hline()


## Plot joint Fst, along 50kb windows
win <- read_tsv("data/selection/fst_50kb_lg12.tsv")

(w <- ggplot(data = win, aes(x = BIN_START, y = WEIGHTED_FST)) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_hline(yintercept = genome_wide, color = "blue") +
    labs(x = "Window", y = "Mean Fst") +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
)


## Identify 99.5% quantile
quantile(win$WEIGHTED_FST, probs = 0.995)

vin <- win %>%
  mutate(OUTLIER = case_when(
    WEIGHTED_FST > quantile(WEIGHTED_FST, probs = 0.995) ~ "yes",
    TRUE ~ "no")
  )


## Plot 99.5% quantile
(v <- ggplot(data = vin, aes(x = BIN_START, y = WEIGHTED_FST, color = OUTLIER)) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_hline(yintercept = genome_wide, color = "blue") +
    labs(x = "Window", y = "Mean Fst") +
    scale_color_manual(values = c("black", "red")) +
    guides(color = "none") +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
)


## Retrieve window positions of peaks
max(win$WEIGHTED_FST)
print(vin %>% filter(OUTLIER == "yes"), n = nrow(vin))


## Genomic regions of interest
# Region1: 20160000, 20285000 (125 kb)
# Region2: 22150000, 22245000  (95 kb)


# Extract genomic regions of interest from vcf
# vcftools --gzvcf snps_lg12_phased_36.vcf.gz \
# --chr LG12 \
# --from-bp 20160000 \
# --to-bp 20285000 \
# --recode \
# --stdout | gzip > region1.vcf.gz
#  
# vcftools --gzvcf snps_lg12_phased_36.vcf.gz \
# --chr LG12 \
# --from-bp 22150000 \
# --to-bp 22245000 \
# --recode \
# --stdout | gzip > region2.vcf.gz


## PCA based on genomic regions of interest (region 1)
region1_vcf <- read.vcfR("data/selection/region1.vcf.gz")

region1_gl <- vcfR2genlight(region1_vcf)

pca_region1 <- glPca(region1_gl, nf = 2)

scores_region1 <- as.data.frame(pca_region1$scores) %>%
  rownames_to_column("Sample") %>%
  mutate(Species = str_sub(Sample, -6, -4)) %>%
  as_tibble()


## If the code above is not working, load data from file
# scores_region1 <- read_csv("data/selection/pca_region1_scores.csv")


## Plot region-specific PCA
(r1 <- ggplot(data = scores_region1, aes(x = PC1, y = PC2, color = Species)) +
  geom_point(size = 4, alpha = 0.75) +
  labs(title = "Region 1") +
  scale_color_manual(values = c("goldenrod2", "royalblue3", "grey20", "coral2", "grey80")) +
  theme_light(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.x = element_text(vjust = -1.5))
)



### ================== Solutions ================= ###

## Simulating the effect of selection on allele frequencies
# Selection for a dominant phenotype: e.g. red 1, purple 1, blue 0.7
# -> Red quickly goes to fixation (directional selection)
# Selection against a recessive phenotype: e.g. red 0.7, purple 1, blue 1
# -> Red remains a low frequency (stabilizing selection)
# Heterozygote advantage: e.g. red 0.7, purple 1, blue 0.7
# -> red and blue allele converge at 0.5, heterozygote excess (balancing selection)
# General observation: changing the strength of selection affects time scale


## Add chromosome-wide Fst
genome_wide <- 0.052
p + geom_hline(yintercept = genome_wide, color = "blue")



### =============== Bonus material =============== ###
  
### Identify region2 by BLAST

## Go to: https://blast.ncbi.nlm.nih.gov/Blast.cgi
## Choose blastx
## Enter region2_sub.fas as query and run BLAST with default settings


### Smooth (spline interpolation)
smooth <- as.data.frame(spline(vin$BIN_START, vin$WEIGHTED_FST))

(s <- ggplot() +
    geom_line(data = smooth, aes(x = x, y = y), color = "gray20") +
    labs(x = "Window position", y = "Mean Fst") +
    geom_hline(yintercept = genome_wide, color = "blue") +
    scale_y_continuous(limits=c(0, 0.6)) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
)


## Highlight window with max Fst
max(win$WEIGHTED_FST)
(peak <- win[which(win$WEIGHTED_FST == max(win$WEIGHTED_FST)), "BIN_START"] %>% pull)

s + geom_vline(xintercept = peak, color = "red", alpha = 0.5)

s + xlim(1.8e+07, 2.2e+07)
s + xlim(2.00e+07, 2.05e+07)


## PCA based on genomic regions of interest (region 2)
region2_vcf <- read.vcfR("data/selection/region2.vcf.gz")

region2_gl <- vcfR2genlight(region2_vcf)

pca_region2 <- glPca(region2_gl, nf = 2)

scores_region2 <- as.data.frame(pca_region2$scores) %>%
  rownames_to_column("Sample") %>%
  mutate(Species = str_sub(Sample, -6, -4)) %>%
  as_tibble()


## If the code above is not working, load data from file
# scores_region2 <- read_csv("data/selection/pca_region2_scores.csv")


## Plot region-specific PCA
(r2 <- ggplot(data = scores_region2, aes(x = PC1, y = PC2, color = Species)) +
    geom_point(size = 4, alpha = 0.75) +
    labs(title = "Region 2") +
    scale_color_manual(values = c("goldenrod2", "royalblue3", "grey20", "coral2", "grey80")) +
    theme_light(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.title.x = element_text(vjust = -1.5))
)
