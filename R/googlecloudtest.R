
#### Install and load packages ####
# install.packages("remotes")
# install.packages("devtools")
# install.packages("rstudioapi") # Get location of the script
# install.packages("readxl")
# install.packages("PBSadmb")
# devtools::install_github("cmlegault/ASAPplots", build_vignettes = TRUE)
# remotes::install_github("r4ss/r4ss", ref="development")
# remotes::install_github("Bai-Li-NOAA/Age_Structured_Stock_Assessment_Model_Comparison")
# install.packages("doParallel")
# install.packages("foreach")

# setwd("C:/Users/bai.li/Documents/Github/Age_Structured_Stock_Assessment_Model_Comparison/")
# devtools::load_all()

library(readxl)
library(PBSadmb)
library(ASAPplots)
library(r4ss)
library(ASSAMC)
library(rstudioapi)
library(parallel)
library(doParallel)
library(foreach)

#### Set working directory ####
dir <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(dir)
getwd()

## Setup working directory, number of simulations, and seed number
maindir <- file.path(dir, "results")
em_input_filenames <- read.csv(file.path(maindir, "em_input", "em_input_filenames.csv"))

## Downloading stock assessment packages 
mydir <- file.path(maindir, "em_input")

url <- "https://vlab.ncep.noaa.gov/documents/259399/11683431/ss/2bd2329d-4623-c420-a360-abd6c3dd52c1?version=1.0&t=1599077106790&download=true"
download.file(file.path(url, "ss"), file.path(mydir, "ss"), mode="wb")
# chmod +x ss (in terminal)

#### Operating model setup ####
om_sim_num <- 5 # total number of iterations per case
keep_sim_num <- 3 # number of kept iterations per case
figure_number <- 3 # number of individual iteration to plot

seed_num <- 9924

#### Life-history parameters ####
year <- 1:30
ages <- 1:12   #Age structure of the popn

initial_equilibrium_F <- TRUE
median_R0 <- 1000000 #Average annual unfished recruitment (scales the popn)
median_h <- 0.75 #Steepness of the Beverton-Holt spawner-recruit relationship.
mean_R0 <- NULL
mean_h <- NULL
SRmodel <- 1 # 1=Beverton-Holt; 2=Ricker
M <- 0.2       #Age-invariant natural mortality

Linf <- 800	  #Asymptotic average length
K <- 0.18     	#Growth coefficient
a0 <- -1.36    #Theoretical age at size 0
a.lw <- 0.000000025  #Length-weight coefficient
b.lw <- 3.0     	    #Length-weight exponent
A50.mat <- 2.25      #Age at 50% maturity
slope.mat <- 3       #Slope of maturity ogive
pattern.mat <- 1     #Simple logistic maturity
female.proportion <- 0.5   #Sex ratio

#### Fleet settings ####
fleet_num <- 1

#CV of landings for OM
cv.L <- list()
cv.L$fleet1 <- 0.005

#Input CV of landings for EMs
input.cv.L <- list()
input.cv.L$fleet1 <- 0.01

#Annual sample size (nfish) of age comp samples
n.L <- list()
n.L$fleet1 <- 200

#Define fleet selectivity
sel_fleet <- list()

sel_fleet$fleet1$pattern <- 1
sel_fleet$fleet1$A50.sel1 <- 2
sel_fleet$fleet1$slope.sel1 <- 1

#### Survey settings ####
survey_num <- 1

#CV of surveys for OM
cv.survey <- list()
cv.survey$survey1 <- 0.1

#Input CV of surveys for EMs
input.cv.survey <- list()
input.cv.survey$survey1 <- 0.2

#Annual sample size (nfish) of age comp samples
n.survey <- list()
n.survey$survey1 <- 200

#Define survey selectivity
sel_survey <- list()

sel_survey$survey1$pattern <- 1
sel_survey$survey1$A50.sel1 <- 1.5
sel_survey$survey1$slope.sel1 <- 2

#### Other settings ####
logf_sd <- 0.2
f_dev_change <- FALSE
f_pattern <- 1
start_val <- 0.01
middle_val <- NULL
end_val <- 0.39
f_val <- NULL
start_year <- 1
middle_year <- NULL

logR_sd <- 0.6
r_dev_change <- TRUE

om_bias_cor <- TRUE
bias_cor_method <- "median_unbiased" #Options: "none", "median_unbiased", and "mean_unbiased"
em_bias_cor <- TRUE

#### Case 1: B-H median-unbiased model and mismatch ####
## No adhoc bias adjustment in AMAK and ASAP
## BAM input steepness = median-unbiased steepness = 0.75
## SS input steepness = median-unbiased steepness = 0.75

case1 <- save_initial_input(base_case=TRUE, case_name = "C1")
run_om(input_list=case1, show_iter_num=F)

source(file.path(dir, "googlecloudtest", "run_ss.R"))


