## Analysis script for the probabilioty of agreement methodology

## Import required source code.
source("PAgreement_Source_Code.R")

## Read in the data:
rr <- read.csv(file = "RRdataGSvsCB.csv", header = T)
## Note that the dataframe must have four columns with headers 'Subject', 'MS', 'Repeat' and 'Measurements', 
## and N = sum{m_ij} rows. The 'Subject' column contains subject codes and the entries of 'MS' must be 1's and 
## 2's indicating which measurement system a given row corresponds to. The entries of 'Repeat' indicate which 
## repeated measurement a given row corrsponds to, and the entries of the 'Measurements' column correspond to 
## the repeated measurements by the given system on each of the n subjects. Note that the observations must be 
## ordered by subject with repeated measurements appearing consecutively with MS1 measurements appearing first 
## and MS2 second.

## Construct relevant plots:

## Scatter plot (ignoring repeated measurements)
scatter.plot(data = rr, n = 21)
## Scatter plot (accounting for repeated measurements)
rep.scatter.plot(data = rr, n = 21)
## Repeatability plot
repeat.plot(data = rr, n = 21)
## Modified QQ-plot
mod.qqplot(data = rr, n = 21)

## Probability of agreement analysis:
p <- prob.agree(n = 21, delta = 5, data = rr)

## Inspect parameter estimates and standard errors. (Order of return: mu, alpha, beta, sigma_s, sigma_1, sigma_2).
p$Estimates
p$St.Errors

## Inspect probability of agreement estimates.
p$PAgree
