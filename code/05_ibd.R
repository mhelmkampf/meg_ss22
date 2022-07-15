### ======================== ###
### 05 Isolation by distance ###
### ======================== ###


## Working directory
# getwd()   use to see whether your current directory is meg_ss22
# setwd()   set current directory to meg_ss22 if necessary

## Install packages
# install.packages("tidyverse")
# install.packages("adegenet")
# install.packages("hierfstat")
# install.packages("genepop")
install.packages("ape")


### Exercise 1: Population-specific Fst

## Load packages
library(tidyverse)
library(adegenet)
library(hierfstat)
library(genepop)


## Read in data (file extensions: Windows = .txt, macOS = .gen)
hamlets <- read.genepop("data/msats/hamlets_caribbean.gen", ncode = 3)


## Population-specific Fst (average for each population)
(bh <- betas(hamlets))


## Convert output to tidyverse data frame (tibble)
fh <- bh$betaiovl %>%                              # extracts Fsts from betas object
  bind_rows() %>%                                  # converts to tibble
  pivot_longer(cols = everything()) %>%            # transforms data from wide to long format
  rename("Population" = name, "Fst" = value) %>%   # renames columns
  arrange(desc(Fst))                               # sorts data by descending Fst


## Plot with ggplot
ggplot(data = fh, aes(x = reorder(Population, -Fst), y = Fst)) +       # determines how data is mapped
  geom_bar(stat = "identity", color = "grey20", fill = "skyblue3") +   # uses barplot geom, sets outline and fill color
  labs(x = "Population") +                                             # changes x-axis label
  theme_minimal(base_size = 14) +                                      # uses basic theme "minimal" (e.g. no frame, white background)
  theme(                                                               # more specific theme elements follow:
    panel.grid.minor = element_blank(),                                # removes minor grid lines
    panel.grid.major.x = element_blank(),                              # removes major grid lines intersecting x-axis
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.25),   # adjusts angle and position of x-axis labels
    )


### Exercise 2: Pairwise Fst

## Calculate pairwise Fst between all populations
Fst("data/msats/puella_caribbean.txt", pairs = TRUE, outputFile = "local/puella_FstPairs.txt")


## Test for differentiation between populations (exact G test)
test_diff("data/msats/puella_caribbean.txt", pairs = TRUE, outputFile = "local/puella_diffPairs.txt", iterations = 1000)


## Clean up temporary files (may not work on Windows)
system("rm fichier.in cmdline.txt")


### Exercise 3: Isolation by distance

## Read in data
geo <- read_csv("data/msats/puella_geo.csv", col_names = TRUE, col_types = "cd")   # pairwise geographical distances in km
fst <- read_csv("data/msats/puella_fst.csv", col_names = TRUE, col_types = "cd")   # pairwise Fst between populations


## Convert data from wide to long format
geo_long <- pivot_longer(geo,
                         cols = 2:15,
                         names_to = "popA",
                         values_to = "distance",
                         values_drop_na = TRUE)

fst_long <- pivot_longer(fst,
                         cols = 2:15,
                         names_to = "popA",
                         values_to = "distance",
                         values_drop_na = TRUE)


## Create data frame with both geographical and genetic distances 
all_dist <- full_join(geo_long, fst_long, by = c("popB", "popA")) %>%
  rename(distance = distance.x, fst = distance.y) %>%
  select(popA, popB, distance, fst)


## Plot geographical over genetic distance
a <- ggplot(data = all_dist, aes(x = distance, y = fst)) +
  geom_point() +
  labs(title = "Caribbean", x = "Distance [km]", y = "Fst") +
  theme_minimal(base_size = 14) +
  theme(axis.title.x = element_text(vjust = -1.75),
        axis.title.y = element_text(vjust = 2)
        )
a


## Linear regression analysis
mod_a <- lm(formula = fst ~ distance, data = all_dist)
summary(mod_a)

a + geom_smooth(method = "lm", se = FALSE, color = "red")


## Guna Yala only
guna_dist <- all_dist %>%
  filter(popA %in% c(12:16) & popB %in% c(12:16))

(g <- ggplot(data = guna_dist, aes(x = distance, y = fst)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Guna Yala", x = "Distance [km]", y = "Fst") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = -0.17),
        axis.title.x = element_text(vjust = -1.75),
        axis.title.y = element_text(vjust = 2)
  )
)

mod_g <- lm(formula = fst ~ distance, data = guna_dist)
summary(mod_g)


## Honduras and Belize only (Exercise)
honbel_dist <- all_dist %>%
  filter(popA %in% c(1:5) & popB %in% c(1:5))

(h <- ggplot(data = honbel_dist, aes(x = distance, y = fst)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    labs(title = "Honduras & Belize", x = "Distance [km]", y = "Fst") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = -0.17),
          axis.title.x = element_text(vjust = -1.75),
          axis.title.y = element_text(vjust = 2)
    )
)

# mod_h <- lm(formula = fst ~ distance, data = honbel_dist)
# summary(mod_h)


## Honduras, Belize & Panama
hbp_dist <- all_dist %>%
  filter(popA %in% c(1:5, 7, 9:16) & popB %in% c(1:5, 7, 9:16))

(p <- ggplot(data = hbp_dist, aes(x = distance, y = fst)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    labs(title = "Honduras, Belize & Panama", x = "Distance [km]", y = "Fst") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = -0.17),
          axis.title.x = element_text(vjust = -1.75),
          axis.title.y = element_text(vjust = 2)
    )
)

mod_p <- lm(formula = fst ~ distance, data = hbp_dist)
summary(mod_p)


## Mantel test for IBD
library(ape)

mgeo <- as.matrix(geo[, -1])   # convert to matrix (discard first column)
mfst <- as.matrix(fst[, -1])   # convert to matrix (discard first column)

mgeo[upper.tri(mgeo)] <- t(mgeo)[upper.tri(mgeo)]   # copy bottom matrix triangle to top
mfst[upper.tri(mfst)] <- t(mfst)[upper.tri(mfst)]   # copy bottom matrix triangle to top

mantel.test(mgeo, mfst, permutations = 10000)   # perform Mantel test



### =============== Bonus material =============== ###

## Plot as panels
# install.packages("gridExtra")
# library(gridExtra)
# grid.arrange(g, h, p, a)
