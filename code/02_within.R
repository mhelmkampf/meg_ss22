### ============================= ###
### 02 Within population genetics ###
### ============================= ###


## Check working directory if unclear
getwd()

## Install packages (if not already installed)
# install.packages("adegenet")
# install.packages("pegas")
# install.packages("genepop")
# install.packages("poppr")


### Hardy-Weinberg example calculation (1 locus, 2 alleles)

## Allele frequencies
n <- 10  # no. individuals
y <- ((6 * 2) + 1) / (n * 2)  # yellow
b <- ((3 * 2) + 1) / (n * 2)  # blue

## Heterozygosity
Ho <- 1 / n
He <- 2 * y * b

## Fixation index
Fis = (He - Ho) / He

## Expected genotype frequencies
yye <- y^2
ybe <- 2 * y * b
bbe <- b^2

## Prepare data frame
dat <- data.frame(row.names = c("YY", "YB", "BB"),
                  "observed" = c(6, 1, 3),
                  "expected" = c(yye, ybe, bbe)
                 )

## Perform chi-squared test of goodness of fit. What does this mean?
test <- chisq.test(dat$observed, p = dat$expected)
test
pchisq(test$statistic, df = 1, lower.tail = FALSE) 

# Note: The Chi-square test is normally not recommended for such low sample sizes
# chisq.test() does not allow changing the degrees of freedom through a parameter, 
# so we need to recalculate the p-value with pchisq()


### Encoding example data in Genepop format and using R packages

library(adegenet)
library(pegas)
library(genepop)
library(poppr)

help(read.genepop)
yellowblue <- read.genepop("data/msats/yellowblue.gen", ncode = 1)

yellowblue
yellowblue@tab
yellowblue@loc.n.all
yellowblue@all.names

hw.test(yellowblue)


### Real-world data: microsatellite data of a reef fish population from Barbados

barbados <- read.genepop("data/msats/puella_barbados.txt", ncode = 3)   # file extensions: Windows = .txt, macOS = .gen
barbados

## What is the most / least diverse locus in terms of number of alleles?
barbados@loc.n.all   # Number of loci and alleles
barbados@all.names   # Names of alleles for each locus

## Locus summary using poppr (no. alleles, heterozygosity, evenness)
locus_table(barbados)

## Is this population in HWE? What does this tell us?
test_HW("data/msats/puella_barbados.txt", outputFile = "local/barbados_HW.txt")
hw.test(barbados)  # alternative using pegas package

## Compare to Cayo de Media Luna population
medialuna <- read.genepop("data/msats/puella_medialuna.txt", ncode = 3)   # file extensions: Windows = .txt, macOS = .gen
test_HW("data/msats/puella_medialuna.txt", outputFile = "local/medialuna_HW.txt")
# hw.test(medialuna)


### All populations: what patterns can you identify? (locus / pop outliers)
caribbean <- read.genepop("data/msats/puella_caribbean.txt", ncode = 3)   # file extensions: Windows = .txt, macOS = .gen
test_HW("data/msats/puella_caribbean.txt", outputFile = "local/caribbean_HW.txt")
# hw.test(caribbean)


### Clean up (may not work on Windows)
system("rm fichier.in cmdline.txt")
