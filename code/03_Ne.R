### ============================ ###
### 03 Effective population size ###
### ============================ ###

## Working directory
# getwd()   use to see whether your current directory is meg_ss22
# setwd()   set current directory to meg_ss22 if necessary

## Install packages
install.packages("tidyverse")
install.packages("NEff")

## Load packages
library(adegenet)
library(genepop)
library(pegas)


### Exercise 1: Test whether self-fertilization occurs in hamlets

## Read in data
barbados <- read.genepop("data/puella_barbados.gen", ncode = 3)   # file extensions: Windows = .txt, macOS = .gen
caribbean <- read.genepop("data/puella_caribbean.gen", ncode = 3)   # file extensions: Windows = .txt, macOS = .gen

## Calculate Fis
pegas::Fst(as.loci(barbados))
pegas::Fst(as.loci(caribbean))

## Test for heterozygote deficiency
help(test_HW)

test_HW("data/puella_barbados.txt", which = "deficit", outputFile = "local/barbados_Hdef.txt")
test_HW("data/puella_caribbean.txt", which = "global deficit", outputFile = "local/caribbean_Hdef.txt")

## Clean up (may not work on Windows)
system("rm fichier.in cmdline.txt")


### Exercise 2: Fluctuating population size

N <- c(100, 120, 80, 110, 10)
N
Na <- mean(N)                 # arithmetic mean: Na <- (100 + 120 + 80 + 110 + 10) / 5
Na
Ne <- 1 / mean(1 / N)         # harmonic mean: Ne <- 5 / (1/100 + 1/120 + 1/80 + 1/110 + 1/10)
Ne
Ne / Na                       # ratio

## Plot Ne over time after bottleneck
N <- c(100, 120, 80, 110, 10, rep(100, 45))

for (i in 1:length(N)) {
  Ne[i] <- 1 / mean(1 / N[1:i])
}

df <- data.frame(Generation = 1:length(Ne), Ne)

plot(x = df$Generation, y = df$Ne)

## Optional: fancier plotting with ggplot
library(tidyverse)

ggplot(data = df, aes(x = Generation, y = Ne)) + 
  geom_point() + geom_line() +
  theme_bw(base_size = 16)


### Exercise 3: Estimate Ne from genetic data with NeEstimator

## Software available in other/software or at http://www.molecularfisherieslaboratory.com.au/neestimator-software/
## Run command line version (e.g. Ne2-1.exe on Windows) or GUI version (NeEstimator2x1.jar)
## Requires multi-locus diploid genotypes in Genepop format
## Relative path to data files is ../../../data/
## For method choose heterozygote excess
## Results will be saved in /other/software/NeEstimator


### Exercise 4: Explore population model (homework)

library(NEff)

## Explore the influence of demographic parameters on N, Ne, and heterozygosity (change parameters carefully)
population(max.repeat = 1, max.time = 50, quiet = FALSE, N.allele = 20, N.loci = 10,
           Na = 100,                                # starting population
           n.recruit.fem = 300,                     # no. offspring per female
           delay = 2,                               # years to sexual maturity
           sex.ratio = 0.5,                         # proportion of females
           surv.ad.fem = 0.3, surv.ad.mal = 0.3,    # adult female and male survival rate
           surv.juv = 0.015, const = FALSE,         # juvenile survival rate
           skew.recruit.fem = 0,                    # female variance in reproduction (variance around mean)
           skew.recruit.mal = 0                     # male variance in reproduction (variance around mean)
)



### =============== Bonus material =============== ###

### Uneven breeding sex ratio (elephant seal example)

N <- 625
Nm <- 1:625
Nf <- N - Nm

Ne <- (4 * Nm * Nf) / (Nm + Nf)

df <- data.frame(Nm, Ne)

g <- ggplot(data = df, aes(x = Nm, y = Ne)) +
  geom_line(size = 1) +
  theme_bw(base_size = 16)
g

g + geom_vline(xintercept = 21, color = "red", size = 1) 


### Missing data (poppr)

all_gc <- as.genclone(caribbean)
info_table(all_gc, type = "missing", plot = TRUE)   # save manually as 10 x 10
