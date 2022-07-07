### ================================ ###
### 01 Introduction to R and RStudio ###
### Marine Ecological Genetics SS22  ###
### ================================ ###


## Simple math
3 + 5 / 2
(3 + 5) / 2


## Creating and assigning objects: object_name <- ...
a <- 1
a
b

a <- 2
a
a + 5

b <- a + 5
b

c <- "hello"
a + c


## Example data frame (table)
iris


## Basic R syntax: function(object, parameters)
head(iris)

help(head)   # Getting help with a function

head(iris, n = 5)

View(iris)


## Extracting data and basic stats
iris$Sepal.Width

# Exercise: calculate maximum, mean  and standard deviation of sepal width
# Hint: search online (e.g. stackoverflow.com/questions)
max(iris$Sepal.Width)
mean(iris$Sepal.Width)
sd(iris$Sepal.Width)

# Exercise: assign petal length of first 20 rows to new object
# Possible solutions:
pl20 <- head(iris$Petal.Length, n = 20)
pl20 <- iris$Petal.Length[1:20]


## Simple plotting: histogram
hist(iris$Sepal.Width)

# Exercise: drop plot title, fix x axis label, and color columns blue:
# Hint: use "help(hist)"
# Possible solution:
hist(iris$Sepal.Width,
     main = "",
     xlab = "Sepal Width",
     ylab = "Frequency",
     col = "blue")


## Plotting and stats example 2: Scatter plot with regression line
plot(x = iris$Sepal.Length,
     y = iris$Petal.Length)

reg <- lm(formula = Petal.Length ~ Sepal.Length,
          data = iris)

abline(reg, col = "red")


## Packages
install.packages("tidyverse")

library("tidyverse")

sessionInfo()


## Preparation for next time: set up Rproject with git
# File | New Project ... | Version Control | Git
# Repository URL: https://github.com/mhelmkampf/meg_ss22
# Project directory name: meg_ss22
# Subdirectory: ~


## Links and resources
# R cheatsheet: https://cran.r-project.org/doc/contrib/Short-refcard.pdf
# R for Beginners: YaRrr! The Pirate's Guide to R, https://bookdown.org/ndphillips/YaRrr/
# Advanced R: R for Data Science, https://r4ds.had.co.nz/index.html
