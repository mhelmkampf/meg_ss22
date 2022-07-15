### ======================= ###
### 04 Population structure ###
### ======================= ###


## Working directory
# getwd()   use to see whether your current directory is meg_ss22
# setwd()   set current directory to meg_ss22 if necessary

## Install packages
# install.packages("adegenet")
# install.packages("genepop")
# install.packages("pegas")
# install.packages("NEff")
# install.packages("tidyverse")
install.packages("hierfstat")


### Exercise 1: Population model

library(NEff)

## Explore the influence of demographic parameters on N, Ne, and heterozygosity (change parameters carefully)
population(max.repeat = 1, max.time = 50, quiet = FALSE, N.allele = 20, N.loci = 10,
           Na = 100,                               # starting population
           n.recruit.fem = 302,                    # no. offspring per female
           delay = 2,                              # years to sexual maturity
           sex.ratio = 0.5,                        # proportion of females
           surv.ad.fem = 0.3, surv.ad.mal = 0.3,   # adult female and male survival rate
           surv.juv = 0.015,                       # juvenile survival rate
           skew.recruit.fem = 0,                   # female variance in reproduction (variance around mean)
           skew.recruit.mal = 0,                   # male variance in reproduction (variance around mean)
           const = FALSE)


### Exercise 2: Estimate Ne from genetic data with NeEstimator

## Software available in apps/NeEstimator
## To run, execute in terminal Ne2-1.exe on Windows, Ne2M on Mac, or Ne2L on Linux
## Provide path to input file when prompted, e.g. ../../data/msats/puella_caribbean.gen)
## For method choose heterozygote excess (2)
## Results will be saved in apps/NeEstimator


### Exercise 3: Calculate global Fst and test for differentiation

## Load packages
library(adegenet)
library(genepop)
library(pegas)


## Read in data (file extensions: Windows = .txt, macOS = .gen)
caribbean <- read.genepop("data/msats/puella_caribbean.gen", ncode = 3)
temporal <- read.genepop("data/msats/puella_temporal.gen", ncode = 3)
hamlets <- read.genepop("data/msats/hamlets_caribbean.gen", ncode = 3)


## Calculate global Fst (per-locus and overall)
pegas::Fst(as.loci(caribbean))
pegas::Fst(as.loci(temporal))
pegas::Fst(as.loci(hamlets))

genepop::Fst("data/msats/puella_caribbean.txt", outputFile = "local/caribbean_Fst.txt")
genepop::Fst("data/msats/puella_temporal.txt", outputFile = "local/temporal_Fst.txt")
genepop::Fst("data/msats/hamlets_caribbean.txt", outputFile = "local/hamlets_Fst.txt")


## Test for differentiation (G-test)
test_diff("data/msats/puella_caribbean.txt", outputFile = "local/caribbean_diff.txt")   # how to? help(genepop)
test_diff("data/msats/puella_temporal.txt", outputFile = "local/temporal_diff.txt")
test_diff("data/msats/hamlets_caribbean.txt", outputFile = "local/hamlets_diff.txt")


### Exercise 4: Calculate population-specific Fst

## Load packages
library(tidyverse)
library(hierfstat)


## Population-specific Fst
betas(hamlets)


## Convert output to tidyverse data frame (tibble)
b <- betas(hamlets)$betaiovl %>%                   # extracts Fsts from betas object
  bind_rows() %>%                                  # converts to tibble
  pivot_longer(cols = everything()) %>%            # transforms data from wide to long format
  rename("Population" = name, "Fst" = value) %>%   # renames columns
  arrange(desc(Fst))                               # sorts data by descending Fst


## Plot with ggplot
ggplot(data = b, aes(x = reorder(Population, -Fst), y = Fst)) +        # determines how data is mapped
  geom_bar(stat = "identity", color = "grey20", fill = "skyblue3") +   # uses barplot geom, sets outline and fill color
  labs(x = "Population") +                                             # changes x-axis label
  theme_minimal(base_size = 16) +                                      # uses basic theme "minimal" (e.g. no frame, white background)
  theme(
    panel.grid.minor = element_blank(),                                # removes minor grid lines
    panel.grid.major.x = element_blank(),                              # removes major grid lines intersecting x-axis
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.25),   # adjusts angle and position of x-axis labels
    )


## Test differentiation between populations
test_diff("data/msats/hamlets_caribbean.txt", pairs = TRUE, outputFile = "local/hamlets_diffpairs.txt", iterations = 1000)


## Clean up temporary files (may not work on Windows)
system("rm fichier.in cmdline.txt")



### =============== Bonus material =============== ###

### PCoA of genetic distances

d <- genet.dist(hamlets, method = "Ds")   # Nei's genetic distance
pcoa(as.matrix(d))
