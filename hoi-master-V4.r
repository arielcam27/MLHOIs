# required packages
library(plotly)       # plotting
require(ggpubr)       # plotting
library("FME")        # fitting models to data
require(growthrates)  # fitting models to data 
library(caret)        # ML models
library(pROC)         # ML models
library(e1071)        # ML models
library(nlme)         # Statistics
require(sets)         # Data structure
require(magrittr)     # pipe operator

# IMPORTANTE: JuliaLang and JuliaCall required for SDEs solver.
library(JuliaCall)
julia_setup(JULIA_HOME = "/home/user/julia-0.6.4/bin") # change dir to local Julia
diffeqr::diffeq_setup()

# Scripts.
source("get_single_V4.r")   # estimate r, K from single-species
source("get_pairwise_V4.r") # estimate a_ij from pairwise-species
source("get_three_V4.r")    # estimate b_ijk from threewise-species
source("gen_samples_V4.r")  # generate HOI samples
source("machine_V4.r")      # ML models
source("bender_V4.r")       # Bender-Case test
source("wootton_V4.r")      # Wootton test

#-----------

# STEP 1: Read data and estimate r, K from 1-species data.

# Input: Three csv files for each species.
# Current version requires each csv file to have this structure:
# ____________________________
# |week|replicate|individuals|
# ----------------------------

# Current version considerates hard-coded 7 weeks for each replicate.

message("Estimating r, K from single species data...")

# INPUT: csv file.
# OUTPUT: df with time-series ("data") and estimated y0, r, K ("pars").

neb1 <- get_single_week1("./data/1-neb.csv") # Nebulosa 1-species data.
sim1 <- get_single_week1("./data/2-sim.csv") # Simulans 1-species data.
ebo1 <- get_single_week1("./data/3-ebo.csv") # Ebony 1-species data.
# TODO: Adapt for n 1-species data.

#------------------

simulateODE1 <- function(df, replicateNum) { 
# Function to visualize 1-species data vs deterministic logistic growth.
# df: dataframe with data (parameters and splitted data).
# replicateNum: replicate number.
  
  f <- function(u,p,t) {
  # Logistic growth ODE.
    r = p[1]
    K = p[2]
    du = r*u*(1.0 - u/K)
    return(c(du))
  }
  
  # parameter vector p=(r, K, y0).
  p <- c(as.numeric(df$pars$r[replicateNum]), 
         as.numeric(df$pars$K[replicateNum]), 
         as.numeric(df$pars$y0[replicateNum]))
  
  # Initial condition.
  u0 = c(p[3])
  
  # Time interval (hard coded, 7 weeks).
  tspan <- list(1, 7)
  
  sol = diffeqr::ode.solve(f, u0, tspan, p=p, saveat=1.0)
  
  udf = as.data.frame(sol$u)
  colnames(udf) <- c("sol")
  tt <- df$data[[1]]$time
  yy <- lapply(df$data[[replicateNum]], as.double)
  udf <- cbind(tt, udf)
  udf <- cbind(udf, yy)
  print(udf)
  fig <- plotly::plot_ly(udf, x = ~time, y = ~sol, type = 'scatter', mode = 'lines')
  fig <- plotly::add_trace(fig, x = ~time, y = ~x1, mode = 'markers')
  fig
}

# Example: Plot Ebony: replicate #1, data and estimated ODE solution.
simulateODE1(ebo1, 1)

#-----------------------------------------

get_mean <- function(df) { 
# Function to retrive 1 value for parameters r, K, y0.
# Used function: Median.
  
  pars_mean = c(median(unlist(lapply(df$pars$r, as.numeric))),
                median(unlist(lapply(df$pars$K, as.numeric))),
                median(unlist(lapply(df$pars$y0, as.numeric))))
  return(pars_mean)
}

# Obtain 1 value for parameters r, K, y0 for each species.
neb_mean = get_mean(neb1) # Nebulosa
sim_mean = get_mean(sim1) # Simulans
ebo_mean = get_mean(ebo1) # Ebony

#-----------------------------------------

# STEP 2: Read data and estimate a_ij from 2-species data.

# Input: Three csv files for each pair of species.
# Current version requires each csv file to have this structure:
# ____________________________________
# |week|replicate|species 1|species 2|
# ------------------------------------

# Current version considerates hard-coded 7 weeks for each replicate.

message("Estimating aij from two-species data...")

amin = -20.0
amax =  20.0

# INPUT: csv file, 2 dataframes with 1-species parameters, and lower-upper limits for aij.
# OUTPUT: df with time-series ("data") and estimated a12, a21, x10, x20.

# Estimate Nebulosa vs Simulans aij.
df_neb_sim0 <- get_pairwise_week1("./data/4-neb_sim.csv", neb_mean, sim_mean, amin, amax)
# Estimate Nebulosa vs Ebony aij.
df_neb_ebo0 <- get_pairwise_week1("./data/5-neb_ebo.csv", neb_mean, ebo_mean, amin, amax)
# Estimate Simulans vs Ebony aij.
df_sim_ebo0 <- get_pairwise_week1("./data/6-sim_ebo.csv", sim_mean, ebo_mean, amin, amax)


simulateODE2 <- function(df1, df2, df3, speciesN, replicateNum) { 
# Function to visualize 2-species data vs deterministic Lotka-Volterra model.
# df: dataframe with data (parameters and splitted data).
# replicateNum: replicate number.  
  
  f <- function(u,p,t) {
  # SDE model: deterministic part.
    
    r1  = p[1]
    r2  = p[2]
    K1  = p[3]
    K2  = p[4]
    a12 = p[5]
    a21 = p[6]
    
    du1 = r1*u[1]*(1 - u[1]/K1 + a12*u[2]/K2)
    du2 = r2*u[2]*(1 - u[2]/K2 + a21*u[1]/K1)
    return(c(du1, du2))
  }
  
  # Parameter vector p=(r1,r2,K1,K2,a12,a21,x10,x20)
  p <- c(as.numeric(df1$pars$r[replicateNum]),
         as.numeric(df2$pars$r[replicateNum]),
         as.numeric(df1$pars$K[replicateNum]),
         as.numeric(df2$pars$K[replicateNum]), 
         as.numeric(df3$pars$a12[replicateNum]),
         as.numeric(df3$pars$a21[replicateNum]),
         as.numeric(df3$pars$x10[replicateNum]),
         as.numeric(df3$pars$x20[replicateNum]))
  
  # Initial condition.
  u0 = c(p[7], p[8])
  # Time interval (hard coded, 7 weeks).
  tspan <- list(1, 7)
  
  sol = diffeqr::ode.solve(f, u0, tspan, p=p, saveat=1.0)
  
  udf = as.data.frame(sol$u)
  tt <- df3$data[[replicateNum]]$time
  yy1 <- df3$data[[replicateNum]]$x1
  
  yy2 <- df3$data[[replicateNum]]$x2
  yy <- data.frame(yy1, yy2)
  udf <- cbind(tt, udf)
  udf <- cbind(udf, yy)
  
  fig <- plotly::plot_ly(x = 1:7, y = ~sol$u[,speciesN], type = 'scatter', mode = 'lines')
  fig <- plotly::add_trace(fig, x = 1:7, y = yy[, speciesN], mode = 'markers')
  fig
}

# Example: Plot Nebulosa vs Simulans: replicate #3, plot species #2 (Simulans) data 
# and estimated ODE solution.
simulateODE2(neb1, ebo1, df_neb_ebo0, 2, 3)

get_mean2 <- function(df) { 
  # Function to retrive 1 value for parameters a12,a21,x10,x20.
  # Used function: Median.
  pars_mean = c(median(unlist(lapply(df$pars$a12, as.numeric))),
                median(unlist(lapply(df$pars$a21, as.numeric))),
                median(unlist(lapply(df$pars$x10, as.numeric))),
                median(unlist(lapply(df$pars$x20, as.numeric))))
  return(pars_mean)
}

# Obtain 1 value for parameters a12, a21, x10, x20 for each pair of species.
neb_sim_mean = get_mean2(df_neb_sim0) # Nebulosa vs Simulans.
neb_ebo_mean = get_mean2(df_neb_ebo0) # Nebulosa vs Ebony.
sim_ebo_mean = get_mean2(df_sim_ebo0) # Simulans vs Ebony.

#------------------------------------------


# STEP 2: Read data and estimate b_ijk from 3-species data.

# Input: One csv file.
# Current version requires each csv file to have this structure:
# ______________________________________________
# |week|replicate|species 1|species 2|species 3|
# ----------------------------------------------

# Current version considerates hard-coded 7 weeks for each replicate.

message("Estimating bijk from three-species data...")

# Hard coded min and max values for bijk parameters.
bmin = -5.0
bmax =  5.0

# 1-species parameters.
parSingle <- cbind(neb_mean, sim_mean, ebo_mean)
# 2-species parameters.
parPair   <- cbind(neb_sim_mean, neb_ebo_mean, sim_ebo_mean)

print(parSingle)
print(parPair)

# INPUT: csv file, 1-species parameters, 2-species parameters, and lower-upper limits for bijk.
# OUTPUT: df with time-series ("data") and estimated b123, b231, b312.

# Estimate bijk.
df_neb_sim_ebo <- get_three_week("./data/7-neb_sim_ebo.csv", parSingle, parPair, bmin, bmax)

# Get 1 value for b123, b231, b312.
get_mean3 <- function(df) { 
  pars_mean = c(median(unlist(lapply(df$pars$b123, as.numeric))),
                median(unlist(lapply(df$pars$b231, as.numeric))),
                median(unlist(lapply(df$pars$b312, as.numeric))))
  return(pars_mean)
}

parThree = get_mean3(df_neb_sim_ebo)
print(parThree)

simulateODE3 <- function(parSingle, parPair, df3, speciesN, replicateNum) { 
  # Function to visualize 3-species data vs deterministic HOI model.
  # df: dataframe with data (parameters and splitted data).
  # replicateNum: replicate number.  
  
  f <- function(u,p,t) {
  # SDE model: deterministic part.
    r1 = parSingle[1,1]
    r2 = parSingle[1,2]
    r3 = parSingle[1,3]
    
    K1 = parSingle[2,1]
    K2 = parSingle[2,2]
    K3 = parSingle[2,3]
    
    a12 = parPair[1,1]
    a21 = parPair[2,1]
    
    a23 = parPair[1,3]
    a32 = parPair[2,3]
    
    a13 = parPair[1,2]
    a31 = parPair[2,2]
    
    b123 = df3$pars[["b123"]][1]
    b231 = df3$pars[["b231"]][1]
    b312 = df3$pars[["b312"]][1]
    
    # b123 = 0
    # b231 = 0
    # b312 = 0
    
    x1 = u[1]
    x2 = u[2]
    x3 = u[3]
    
    du1 = r1*x1*(1 - x1/K1 + a12*x2/K2 + a13*x3/K3 + b123*x2*x3/(K2*K3))
    du2 = r2*x2*(1 - x2/K2 + a23*x3/K3 + a21*x1/K1 + b231*x3*x1/(K3*K1))
    du3 = r3*x3*(1 - x3/K3 + a31*x1/K1 + a32*x2/K2 + b312*x1*x2/(K1*K2))
    return(c(du1, du2, du3))
  }
  
  
  p <- 0
  
  # Initial condition (hard coded, 5 individuals per species in this problem).
  u0 = c(5, 5, 5)
  
  # Time interval (hard coded, 7 weeks).
  tspan <- list(1, 7)
  
  sol = diffeqr::ode.solve(f, u0, tspan, saveat=1.0)
  
  udf = as.data.frame(sol$u)
  tt <- df3$data[[replicateNum]]$time
  yy1 <- df3$data[[replicateNum]]$x1
  #print(yy1)
  yy2 <- df3$data[[replicateNum]]$x2
  yy3 <- df3$data[[replicateNum]]$x3
  yy <- data.frame(yy1, yy2, yy3)
  udf <- cbind(tt, udf)
  udf <- cbind(udf, yy)
  #print(yy)
  fig <- plotly::plot_ly(x = 1:7, y = ~sol$u[,speciesN], type = 'scatter', mode = 'lines')
  fig <- plotly::add_trace(fig, x = 1:7, y = yy[, speciesN], mode = 'markers')
  fig
}

# Example: Plot Nebulosa vs Simulans vs Ebony: 
# replicate #3, plot species #2 (Simulans) data 
simulateODE3(parSingle, parPair, df_neb_sim_ebo, 2, 3)

#--------------------------------


# STEP 3: Generate testing samples.

message("Generating testing samples...")
# gen_samples_V2: Lotka-Volteraa parameters as is. 
# Random initial conditions, bijk random from -10 to 10 (hard coded).

number = 1000

noise = 0.01
df_res1 = gen_samples_V2(number, noise, parSingle, parPair, parThree)

noise = 0.1
df_res2 = gen_samples_V2(number, noise, parSingle, parPair, parThree)

noise = 0.25
df_res3 = gen_samples_V2(number, noise, parSingle, parPair, parThree)

noise = 0.5
df_res4 = gen_samples_V2(number, noise, parSingle, parPair, parThree)

plotHOIs <- function(df) {
plot(1:21, df$df_HOI[1,1:21])
for (ii in 1:number){
  lines(1:21, df$df_HOI[ii,1:21])
}
}

plotNONHOIs <- function(df) {
  plot(1:21, df$df_NONHOI[1,1:21])
  for (ii in 1:number){
    lines(1:21, df$df_NONHOI[ii,1:21])
  }
}

plotHOIs(df_res4)
plotNONHOIs(df_res4)

plot(1:21, df_res4$df_HOI123SDE1[1,1:21])
for (ii in 1:number){
  lines(1:21, df_res4$df_HOI123SDE1[ii,1:21], col='blue')
  lines(1:21, df_res4$df_NONHOI123SDE1[ii,1:21], col='red')
}

plot(1:21, df_res4$df_HOI1SDE1[1,1:21])
for (ii in 1:number){
  lines(1:21, df_res4$df_HOI1SDE1[ii,1:21])
}

  #---------------

# STEP 4: Perform classical tests.

change_ABC <- function(df1, df2, df3, row, column) {
  Out <- data.frame(0, 0, 0)
  colnames(Out) <- c("week", "replicate", "x1")
  Out <- Out[-1,]
  for (ii in 1:7){
    aux <- data.frame(ii, 1, df1[row, 7 * (column - 1) + ii])
    colnames(aux) <- c("week", "replicate", "x1")
    Out <- rbind(Out, aux)
  }

  for (ii in 1:7){
    aux <- data.frame(ii, 2, df2[row, 7 * (column - 1) + ii])
    colnames(aux) <- c("week", "replicate", "x1")
    Out <- rbind(Out, aux)
  }

  for (ii in 1:7){
    aux <- data.frame(ii, 3, df3[row, 7 * (column - 1) + ii])
    colnames(aux) <- c("week", "replicate", "x1")
    Out <- rbind(Out, aux)
  }
  return(Out)
}

print(change_ABC(df_res1$df_HOI13SDE1, df_res1$df_HOI13SDE2, df_res1$df_HOI13SDE3, 1, 1))

check_HOI_bender <- function(df_res, row) {
  aux <- 0
  
  data_A <- change_ABC(df_res$df_HOI1SDE1, df_res$df_HOI1SDE2, df_res$df_HOI1SDE3, row, 1)
  data_AB <- change_ABC(df_res$df_HOI12SDE1, df_res$df_HOI12SDE2, df_res$df_HOI12SDE3, row, 1)
  data_AC <- change_ABC(df_res$df_HOI13SDE1, df_res$df_HOI13SDE2, df_res$df_HOI13SDE3, row, 1)
  data_ABC <- change_ABC(df_res$df_HOI123SDE1, df_res$df_HOI123SDE2, df_res$df_HOI123SDE3, row, 1)
  aux <- aux + bendertest(data_A, data_AB, data_AC, data_ABC)
  
  data_A <- change_ABC(df_res$df_HOI2SDE1, df_res$df_HOI2SDE2, df_res$df_HOI2SDE3, row, 2)
  data_AB <- change_ABC(df_res$df_HOI12SDE1, df_res$df_HOI12SDE2, df_res$df_HOI12SDE3, row, 2)
  data_AC <- change_ABC(df_res$df_HOI23SDE1, df_res$df_HOI23SDE2, df_res$df_HOI23SDE3, row, 2)
  data_ABC <- change_ABC(df_res$df_HOI123SDE1, df_res$df_HOI123SDE2, df_res$df_HOI123SDE3, row, 2)
  aux <- aux + bendertest(data_A, data_AB, data_AC, data_ABC)
  
  data_A <- change_ABC(df_res$df_HOI3SDE1, df_res$df_HOI3SDE2, df_res$df_HOI3SDE3, row, 3)
  data_AB <- change_ABC(df_res$df_HOI13SDE1, df_res$df_HOI13SDE2, df_res$df_HOI13SDE3, row, 3)
  data_AC <- change_ABC(df_res$df_HOI23SDE1, df_res$df_HOI23SDE2, df_res$df_HOI23SDE3, row, 3)
  data_ABC <- change_ABC(df_res$df_HOI123SDE1, df_res$df_HOI123SDE2, df_res$df_HOI123SDE3, row, 3)
  aux <- aux + bendertest(data_A, data_AB, data_AC, data_ABC)
  return(aux)
}

check_NONHOI_bender <- function(df_res, row) {
  aux <- 0

  data_A <- change_ABC(df_res$df_NONHOI1SDE1, df_res$df_NONHOI1SDE2, df_res$df_NONHOI1SDE3, row, 1)
  data_AB <- change_ABC(df_res$df_NONHOI12SDE1, df_res$df_NONHOI12SDE2, df_res$df_NONHOI12SDE3, row, 1)
  data_AC <- change_ABC(df_res$df_NONHOI13SDE1, df_res$df_NONHOI13SDE2, df_res$df_NONHOI13SDE3, row, 1)
  data_ABC <- change_ABC(df_res$df_NONHOI123SDE1, df_res$df_NONHOI123SDE2, df_res$df_NONHOI123SDE3, row, 1)
  aux <- aux + bendertest(data_A, data_AB, data_AC, data_ABC)

  data_A <- change_ABC(df_res$df_NONHOI2SDE1, df_res$df_NONHOI2SDE2, df_res$df_NONHOI2SDE3, row, 2)
  data_AB <- change_ABC(df_res$df_NONHOI12SDE1, df_res$df_NONHOI12SDE2, df_res$df_NONHOI12SDE3, row, 2)
  data_AC <- change_ABC(df_res$df_NONHOI23SDE1, df_res$df_NONHOI23SDE2, df_res$df_NONHOI23SDE3, row, 2)
  data_ABC <- change_ABC(df_res$df_NONHOI123SDE1, df_res$df_NONHOI123SDE2, df_res$df_NONHOI123SDE3, row, 2)
  aux <- aux + bendertest(data_A, data_AB, data_AC, data_ABC)

  data_A <- change_ABC(df_res$df_NONHOI3SDE1, df_res$df_NONHOI3SDE2, df_res$df_NONHOI3SDE3, row, 3)
  data_AB <- change_ABC(df_res$df_NONHOI13SDE1, df_res$df_NONHOI13SDE2, df_res$df_NONHOI13SDE3, row, 3)
  data_AC <- change_ABC(df_res$df_NONHOI23SDE1, df_res$df_NONHOI23SDE2, df_res$df_NONHOI23SDE3, row, 3)
  data_ABC <- change_ABC(df_res$df_NONHOI123SDE1, df_res$df_NONHOI123SDE2, df_res$df_NONHOI123SDE3, row, 3)
  aux <- aux + bendertest(data_A, data_AB, data_AC, data_ABC)
  return(aux)
}

check_HOI_bender(df_res3, 4)
check_NONHOI_bender(df_res3, 10)

message("Check: Bender...")
totalHOIbender1 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "\n") }
  if (check_HOI_bender(df_res1, ii) > 1) {
    totalHOIbender1 <- totalHOIbender1 + 1
  }
}

totalNONHOIbender1 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "\n") }
  aux = check_NONHOI_bender(df_res1, ii)
  if (aux < 2) {
    totalNONHOIbender1 <- totalNONHOIbender1 + 1
  }
}

totalHOIbender2 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "\n") }
  if (check_HOI_bender(df_res2, ii) > 1) {
    totalHOIbender2 <- totalHOIbender2 + 1
  }
}

totalHOIbender3 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "\n") }
  if (check_HOI_bender(df_res3, ii) > 1) {
    totalHOIbender3 <- totalHOIbender3 + 1
  }
}

totalHOIbender4 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "\n") }
  if (check_HOI_bender(df_res4, ii) > 1) {
    totalHOIbender4 <- totalHOIbender4 + 1
  }
}

totalNONHOIbender2 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "\n") }
  aux = check_NONHOI_bender(df_res2, ii)
  if (aux < 2) {
    totalNONHOIbender2 <- totalNONHOIbender2 + 1
  }
}

totalNONHOIbender3 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "\n") }
  aux = check_NONHOI_bender(df_res3, ii)
  if (aux < 2) {
    totalNONHOIbender3 <- totalNONHOIbender3 + 1
  }
}

totalNONHOIbender4 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "\n") }
  aux = check_NONHOI_bender(df_res4, ii)
  if (aux < 2) {
    totalNONHOIbender4 <- totalNONHOIbender4 + 1
  }
}

#----


check_HOI_wootton <- function(df_res, row) {
  aux <- 0

  data_A <- change_ABC(df_res$df_HOI1SDE1, df_res$df_HOI1SDE2, df_res$df_HOI1SDE3, row, 1)
  data_AB <- change_ABC(df_res$df_HOI12SDE1, df_res$df_HOI12SDE2, df_res$df_HOI12SDE3, row, 1)
  data_AC <- change_ABC(df_res$df_HOI13SDE1, df_res$df_HOI13SDE2, df_res$df_HOI13SDE3, row, 1)
  data_ABC <- change_ABC(df_res$df_HOI123SDE1, df_res$df_HOI123SDE2, df_res$df_HOI123SDE3, row, 1)
  aux <- aux + woottontest(data_A, data_AB, data_AC, data_ABC)

  data_A <- change_ABC(df_res$df_HOI2SDE1, df_res$df_HOI2SDE2, df_res$df_HOI2SDE3, row, 2)
  data_AB <- change_ABC(df_res$df_HOI12SDE1, df_res$df_HOI12SDE2, df_res$df_HOI12SDE3, row, 2)
  data_AC <- change_ABC(df_res$df_HOI23SDE1, df_res$df_HOI23SDE2, df_res$df_HOI23SDE3, row, 2)
  data_ABC <- change_ABC(df_res$df_HOI123SDE1, df_res$df_HOI123SDE2, df_res$df_HOI123SDE3, row, 2)
  aux <- aux + woottontest(data_A, data_AB, data_AC, data_ABC)

  data_A <- change_ABC(df_res$df_HOI3SDE1, df_res$df_HOI3SDE2, df_res$df_HOI3SDE3, row, 3)
  data_AB <- change_ABC(df_res$df_HOI13SDE1, df_res$df_HOI13SDE2, df_res$df_HOI13SDE3, row, 3)
  data_AC <- change_ABC(df_res$df_HOI23SDE1, df_res$df_HOI23SDE2, df_res$df_HOI23SDE3, row, 3)
  data_ABC <- change_ABC(df_res$df_HOI123SDE1, df_res$df_HOI123SDE2, df_res$df_HOI123SDE3, row, 3)
  aux <- aux + woottontest(data_A, data_AB, data_AC, data_ABC)
  return(aux)
}

check_NONHOI_wootton <- function(df_res, row) {
  aux <- 0

  data_A <- change_ABC(df_res$df_NONHOI1SDE1, df_res$df_NONHOI1SDE2, df_res$df_NONHOI1SDE3, row, 1)
  data_AB <- change_ABC(df_res$df_NONHOI12SDE1, df_res$df_NONHOI12SDE2, df_res$df_NONHOI12SDE3, row, 1)
  data_AC <- change_ABC(df_res$df_NONHOI13SDE1, df_res$df_NONHOI13SDE2, df_res$df_NONHOI13SDE3, row, 1)
  data_ABC <- change_ABC(df_res$df_NONHOI123SDE1, df_res$df_NONHOI123SDE2, df_res$df_NONHOI123SDE3, row, 1)
  aux <- aux + woottontest(data_A, data_AB, data_AC, data_ABC)

  data_A <- change_ABC(df_res$df_NONHOI2SDE1, df_res$df_NONHOI2SDE2, df_res$df_NONHOI2SDE3, row, 2)
  data_AB <- change_ABC(df_res$df_NONHOI12SDE1, df_res$df_NONHOI12SDE2, df_res$df_NONHOI12SDE3, row, 2)
  data_AC <- change_ABC(df_res$df_NONHOI23SDE1, df_res$df_NONHOI23SDE2, df_res$df_NONHOI23SDE3, row, 2)
  data_ABC <- change_ABC(df_res$df_NONHOI123SDE1, df_res$df_NONHOI123SDE2, df_res$df_NONHOI123SDE3, row, 2)
  aux <- aux + woottontest(data_A, data_AB, data_AC, data_ABC)

  data_A <- change_ABC(df_res$df_NONHOI3SDE1, df_res$df_NONHOI3SDE2, df_res$df_NONHOI3SDE3, row, 3)
  data_AB <- change_ABC(df_res$df_NONHOI13SDE1, df_res$df_NONHOI13SDE2, df_res$df_NONHOI13SDE3, row, 3)
  data_AC <- change_ABC(df_res$df_NONHOI23SDE1, df_res$df_NONHOI23SDE2, df_res$df_NONHOI23SDE3, row, 3)
  data_ABC <- change_ABC(df_res$df_NONHOI123SDE1, df_res$df_NONHOI123SDE2, df_res$df_NONHOI123SDE3, row, 3)
  aux <- aux + woottontest(data_A, data_AB, data_AC, data_ABC)
  return(aux)
}

check_HOI_wootton(df_res2, 2)

message("Check: Wootton...")

totalHOIwootton1 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "\n") }
  if (check_HOI_wootton(df_res1, ii) > 0) {
    totalHOIwootton1 <- totalHOIwootton1 + 1
  }
}

totalHOIwootton2 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "\n") }
  if (check_HOI_wootton(df_res2, ii) > 0) {
    totalHOIwootton2 <- totalHOIwootton2 + 1
  }
}

totalHOIwootton3 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "\n") }
  if (check_HOI_wootton(df_res3, ii) > 0) {
    totalHOIwootton3 <- totalHOIwootton3 + 1
  }
}

totalHOIwootton4 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "\n") }
  if (check_HOI_wootton(df_res4, ii) > 0) {
    totalHOIwootton4 <- totalHOIwootton4 + 1
  }
}

totalNONHOIwootton1 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "\n") }
  if (check_NONHOI_wootton(df_res1, ii) == 0) {
    totalNONHOIwootton1 <- totalNONHOIwootton1 + 1
  }
}

totalNONHOIwootton2 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "\n") }
  if (check_NONHOI_wootton(df_res2, ii) == 0) {
    totalNONHOIwootton2 <- totalNONHOIwootton2 + 1
  }
}

totalNONHOIwootton3 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "\n") }
  if (check_NONHOI_wootton(df_res3, ii) == 0) {
    totalNONHOIwootton3 <- totalNONHOIwootton3 + 1
  }
}

totalNONHOIwootton4 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "\n") }
  if (check_NONHOI_wootton(df_res4, ii) == 0) {
    totalNONHOIwootton4 <- totalNONHOIwootton4 + 1
  }
}


#-----

# STEP 5: Generate training samples for ML models.

message("Creating training sets")

# Number of desired sambples for training ML models.
number_train = 10000

noise = 0.0
df_train1 = gen_samples_training_V2(number_train, noise, parSingle, parPair, parThree)

noise = 0.05
df_train2 = gen_samples_training_V2(number_train, noise, parSingle, parPair, parThree)

noise = 0.1
df_train3 = gen_samples_training_V2(number_train, noise, parSingle, parPair, parThree)

noise = 0.15
df_train4 = gen_samples_training_V2(number_train, noise, parSingle, parPair, parThree)

noise = 0.2
df_train5 = gen_samples_training_V2(number_train, noise, parSingle, parPair, parThree)

HOI1 = df_train1$df_HOI
HOI2 = df_train2$df_HOI
HOI3 = df_train3$df_HOI
HOI4 = df_train4$df_HOI
HOI5 = df_train5$df_HOI
HOI = HOI1
HOI = rbind(HOI, HOI2)
HOI = rbind(HOI, HOI3)
HOI = rbind(HOI, HOI4)
HOI = rbind(HOI, HOI5)

NONHOI1 = df_train1$df_NONHOI
NONHOI2 = df_train2$df_NONHOI
NONHOI3 = df_train3$df_NONHOI
NONHOI4 = df_train4$df_NONHOI
NONHOI5 = df_train5$df_NONHOI
NONHOI = NONHOI1
NONHOI = rbind(NONHOI, NONHOI2)
NONHOI = rbind(NONHOI, NONHOI3)
NONHOI = rbind(NONHOI, NONHOI4)
NONHOI = rbind(NONHOI, NONHOI5)

df_train = list(df_HOI = HOI, df_NONHOI = NONHOI)

plot(1:21, df_train$df_HOI[1,1:21])
for (ii in 1:number_train){
  lines(1:21, df_train$df_HOI[ii,1:21], col='blue')
  lines(1:21, df_train$df_NONHOI[ii,1:21], col='red')
}

dataset <- rbind(df_train$df_HOI, df_train$df_NONHOI)

message("Training: GLM...")
list_glm   <- machine_lin_HOI(dataset)
message("Training: KNN...")
list_knn   <- machine_knn_HOI(dataset)
message("Training: SVM...")
list_svm   <- machine_svm_HOI(dataset)

print(list_glm$acc)
print(list_knn$acc)
print(list_svm$acc)

#plot(list_glm$ROC, print.thres="best", print.thres.best.method="closest.topleft")


#------------

# STEP 6: Challenge ML models.


message("Check: GLM...")

totalHOIglm1 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "
") }
  aux <- 0
  if (predict(list_glm$fit, df_res1$df_HOI123SDE1[ii,]) == "hoi") { aux <- aux + 1 }
  if (predict(list_glm$fit, df_res1$df_HOI123SDE2[ii,]) == "hoi") { aux <- aux + 1 }
  if (predict(list_glm$fit, df_res1$df_HOI123SDE3[ii,]) == "hoi") { aux <- aux + 1 }
  if (aux > 1) { totalHOIglm1 <- totalHOIglm1 + 1 }
}

totalHOIglm2 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "
") }
  aux <- 0
  if (predict(list_glm$fit, df_res2$df_HOI123SDE1[ii,]) == "hoi") { aux <- aux + 1 }
  if (predict(list_glm$fit, df_res2$df_HOI123SDE2[ii,]) == "hoi") { aux <- aux + 1 }
  if (predict(list_glm$fit, df_res2$df_HOI123SDE3[ii,]) == "hoi") { aux <- aux + 1 }
  if (aux > 1) { totalHOIglm2 <- totalHOIglm2 + 1 }
}

totalHOIglm3 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "
") }
  aux <- 0
  if (predict(list_glm$fit, df_res3$df_HOI123SDE1[ii,]) == "hoi") { aux <- aux + 1 }
  if (predict(list_glm$fit, df_res3$df_HOI123SDE2[ii,]) == "hoi") { aux <- aux + 1 }
  if (predict(list_glm$fit, df_res3$df_HOI123SDE3[ii,]) == "hoi") { aux <- aux + 1 }
  if (aux > 1) { totalHOIglm3 <- totalHOIglm3 + 1 }
}

totalHOIglm4 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "
") }
  aux <- 0
  if (predict(list_glm$fit, df_res4$df_HOI123SDE1[ii,]) == "hoi") { aux <- aux + 1 }
  if (predict(list_glm$fit, df_res4$df_HOI123SDE2[ii,]) == "hoi") { aux <- aux + 1 }
  if (predict(list_glm$fit, df_res4$df_HOI123SDE3[ii,]) == "hoi") { aux <- aux + 1 }
  if (aux > 1) { totalHOIglm4 <- totalHOIglm4 + 1 }
}

totalNONHOIglm1 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "
") }
  aux <- 0
  if (predict(list_glm$fit, df_res1$df_NONHOI123SDE1[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (predict(list_glm$fit, df_res1$df_NONHOI123SDE2[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (predict(list_glm$fit, df_res1$df_NONHOI123SDE3[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (aux > 1) { totalNONHOIglm1 <- totalNONHOIglm1 + 1 }
}

totalNONHOIglm2 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "
") }
  aux <- 0
  if (predict(list_glm$fit, df_res2$df_NONHOI123SDE1[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (predict(list_glm$fit, df_res2$df_NONHOI123SDE2[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (predict(list_glm$fit, df_res2$df_NONHOI123SDE3[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (aux > 1) { totalNONHOIglm2 <- totalNONHOIglm2 + 1 }
}

totalNONHOIglm3 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "
") }
  aux <- 0
  if (predict(list_glm$fit, df_res3$df_NONHOI123SDE1[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (predict(list_glm$fit, df_res3$df_NONHOI123SDE2[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (predict(list_glm$fit, df_res3$df_NONHOI123SDE3[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (aux > 1) { totalNONHOIglm3 <- totalNONHOIglm3 + 1 }
}

totalNONHOIglm4 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "
") }
  aux <- 0
  if (predict(list_glm$fit, df_res4$df_NONHOI123SDE1[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (predict(list_glm$fit, df_res4$df_NONHOI123SDE2[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (predict(list_glm$fit, df_res4$df_NONHOI123SDE3[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (aux > 1) { totalNONHOIglm4 <- totalNONHOIglm4 + 1 }
}


#----

message("Check: KNN...")

totalHOIknn1 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "
") }
  aux <- 0
  if (predict(list_knn$fit, df_res1$df_HOI123SDE1[ii,]) == "hoi") { aux <- aux + 1 }
  if (predict(list_knn$fit, df_res1$df_HOI123SDE2[ii,]) == "hoi") { aux <- aux + 1 }
  if (predict(list_knn$fit, df_res1$df_HOI123SDE3[ii,]) == "hoi") { aux <- aux + 1 }
  if (aux > 1) { totalHOIknn1 <- totalHOIknn1 + 1 }
}

totalHOIknn2 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "
") }
  aux <- 0
  if (predict(list_knn$fit, df_res2$df_HOI123SDE1[ii,]) == "hoi") { aux <- aux + 1 }
  if (predict(list_knn$fit, df_res2$df_HOI123SDE2[ii,]) == "hoi") { aux <- aux + 1 }
  if (predict(list_knn$fit, df_res2$df_HOI123SDE3[ii,]) == "hoi") { aux <- aux + 1 }
  if (aux > 1) { totalHOIknn2 <- totalHOIknn2 + 1 }
}


totalHOIknn3 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "
") }
  aux <- 0
  if (predict(list_knn$fit, df_res3$df_HOI123SDE1[ii,]) == "hoi") { aux <- aux + 1 }
  if (predict(list_knn$fit, df_res3$df_HOI123SDE2[ii,]) == "hoi") { aux <- aux + 1 }
  if (predict(list_knn$fit, df_res3$df_HOI123SDE3[ii,]) == "hoi") { aux <- aux + 1 }
  if (aux > 1) { totalHOIknn3 <- totalHOIknn3 + 1 }
}


totalHOIknn4 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "
") }
  aux <- 0
  if (predict(list_knn$fit, df_res4$df_HOI123SDE1[ii,]) == "hoi") { aux <- aux + 1 }
  if (predict(list_knn$fit, df_res4$df_HOI123SDE2[ii,]) == "hoi") { aux <- aux + 1 }
  if (predict(list_knn$fit, df_res4$df_HOI123SDE3[ii,]) == "hoi") { aux <- aux + 1 }
  if (aux > 1) { totalHOIknn4 <- totalHOIknn4 + 1 }
}


totalNONHOIknn1 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "
") }
  aux <- 0
  if (predict(list_knn$fit, df_res1$df_NONHOI123SDE1[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (predict(list_knn$fit, df_res1$df_NONHOI123SDE2[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (predict(list_knn$fit, df_res1$df_NONHOI123SDE3[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (aux > 1) { totalNONHOIknn1 <- totalNONHOIknn1 + 1 }
}

totalNONHOIknn2 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "
") }
  aux <- 0
  if (predict(list_knn$fit, df_res2$df_NONHOI123SDE1[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (predict(list_knn$fit, df_res2$df_NONHOI123SDE2[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (predict(list_knn$fit, df_res2$df_NONHOI123SDE3[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (aux > 1) { totalNONHOIknn2 <- totalNONHOIknn2 + 1 }
}

totalNONHOIknn3 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "
") }
  aux <- 0
  if (predict(list_knn$fit, df_res3$df_NONHOI123SDE1[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (predict(list_knn$fit, df_res3$df_NONHOI123SDE2[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (predict(list_knn$fit, df_res3$df_NONHOI123SDE3[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (aux > 1) { totalNONHOIknn3 <- totalNONHOIknn3 + 1 }
}

totalNONHOIknn4 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "
") }
  aux <- 0
  if (predict(list_knn$fit, df_res4$df_NONHOI123SDE1[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (predict(list_knn$fit, df_res4$df_NONHOI123SDE2[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (predict(list_knn$fit, df_res4$df_NONHOI123SDE3[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (aux > 1) { totalNONHOIknn4 <- totalNONHOIknn4 + 1 }
}

#---

message("Check: SVM...")
  
totalHOIsvm1 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "
") }
  aux <- 0
  if (predict(list_svm$fit, df_res1$df_HOI123SDE1[ii,]) == "hoi") { aux <- aux + 1 }
  if (predict(list_svm$fit, df_res1$df_HOI123SDE2[ii,]) == "hoi") { aux <- aux + 1 }
  if (predict(list_svm$fit, df_res1$df_HOI123SDE3[ii,]) == "hoi") { aux <- aux + 1 }
  if (aux > 1) { totalHOIsvm1 <- totalHOIsvm1 + 1 }
}

totalHOIsvm2 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "
") }
  aux <- 0
  if (predict(list_svm$fit, df_res2$df_HOI123SDE1[ii,]) == "hoi") { aux <- aux + 1 }
  if (predict(list_svm$fit, df_res2$df_HOI123SDE2[ii,]) == "hoi") { aux <- aux + 1 }
  if (predict(list_svm$fit, df_res2$df_HOI123SDE3[ii,]) == "hoi") { aux <- aux + 1 }
  if (aux > 1) { totalHOIsvm2 <- totalHOIsvm2 + 1 }
}

totalHOIsvm3 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "
") }
  aux <- 0
  if (predict(list_svm$fit, df_res3$df_HOI123SDE1[ii,]) == "hoi") { aux <- aux + 1 }
  if (predict(list_svm$fit, df_res3$df_HOI123SDE2[ii,]) == "hoi") { aux <- aux + 1 }
  if (predict(list_svm$fit, df_res3$df_HOI123SDE3[ii,]) == "hoi") { aux <- aux + 1 }
  if (aux > 1) { totalHOIsvm3 <- totalHOIsvm3 + 1 }
}

totalHOIsvm4 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "
") }
  aux <- 0
  if (predict(list_svm$fit, df_res4$df_HOI123SDE1[ii,]) == "hoi") { aux <- aux + 1 }
  if (predict(list_svm$fit, df_res4$df_HOI123SDE2[ii,]) == "hoi") { aux <- aux + 1 }
  if (predict(list_svm$fit, df_res4$df_HOI123SDE3[ii,]) == "hoi") { aux <- aux + 1 }
  if (aux > 1) { totalHOIsvm4 <- totalHOIsvm4 + 1 }
}

totalNONHOIsvm1 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "
") }
  aux <- 0
  if (predict(list_svm$fit, df_res1$df_NONHOI123SDE1[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (predict(list_svm$fit, df_res1$df_NONHOI123SDE2[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (predict(list_svm$fit, df_res1$df_NONHOI123SDE3[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (aux > 1) { totalNONHOIsvm1 <- totalNONHOIsvm1 + 1 }
}

totalNONHOIsvm2 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "
") }
  aux <- 0
  if (predict(list_svm$fit, df_res2$df_NONHOI123SDE1[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (predict(list_svm$fit, df_res2$df_NONHOI123SDE2[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (predict(list_svm$fit, df_res2$df_NONHOI123SDE3[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (aux > 1) { totalNONHOIsvm2 <- totalNONHOIsvm2 + 1 }
}


totalNONHOIsvm3 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "
") }
  aux <- 0
  if (predict(list_svm$fit, df_res3$df_NONHOI123SDE1[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (predict(list_svm$fit, df_res3$df_NONHOI123SDE2[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (predict(list_svm$fit, df_res3$df_NONHOI123SDE3[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (aux > 1) { totalNONHOIsvm3 <- totalNONHOIsvm3 + 1 }
}

totalNONHOIsvm4 <- 0
for (ii in 1:number) {
  if (ii %% floor(number/10) == 0) { cat("Progress:", ii, "
") }
  aux <- 0
  if (predict(list_svm$fit, df_res4$df_NONHOI123SDE1[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (predict(list_svm$fit, df_res4$df_NONHOI123SDE2[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (predict(list_svm$fit, df_res4$df_NONHOI123SDE3[ii,]) == "nonhoi") { aux <- aux + 1 }
  if (aux > 1) { totalNONHOIsvm4 <- totalNONHOIsvm4 + 1 }
}

#----------------------------------------

# STEP 7: Predict HOI/NON-HOI label for experimental data.

dataOut_A <-read.table("./data/1-neb.csv",sep=",")
colnames(dataOut_A) <- c("week", "replicate", "x1")
dataOut_AB <-read.table("./data/4-neb_sim.csv",sep=",")
colnames(dataOut_AB) <- c("week", "replicate", "x1", "x2")
dataOut_AC <-read.table("./data/5-neb_ebo.csv",sep=",")
colnames(dataOut_AC) <- c("week", "replicate", "x1", "x2")
dataOut_ABC <-read.table("./data/7-neb_sim_ebo.csv",sep=",")
colnames(dataOut_ABC) <- c("week", "replicate", "x1", "x2", "x3")
bendertest(dataOut_A, dataOut_AB, dataOut_AC, dataOut_ABC)

dataOut_A <-read.table("./data/2-sim.csv",sep=",")
colnames(dataOut_A) <- c("week", "replicate", "x1")
dataOut_AB <-read.table("./data/4-neb_sim.csv",sep=",")
colnames(dataOut_AB) <- c("week", "replicate", "x2", "x1")
dataOut_AC <-read.table("./data/6-sim_ebo.csv",sep=",")
colnames(dataOut_AC) <- c("week", "replicate", "x1", "x2")
dataOut_ABC <-read.table("./data/7-neb_sim_ebo.csv",sep=",")
colnames(dataOut_ABC) <- c("week", "replicate", "x2", "x1", "x3")
bendertest(dataOut_A, dataOut_AB, dataOut_AC, dataOut_ABC)

dataOut_A <-read.table("./data/3-ebo.csv",sep=",")
colnames(dataOut_A) <- c("week", "replicate", "x1")
dataOut_AB <-read.table("./data/5-neb_ebo.csv",sep=",")
colnames(dataOut_AB) <- c("week", "replicate", "x2", "x1")
dataOut_AC <-read.table("./data/6-sim_ebo.csv",sep=",")
colnames(dataOut_AC) <- c("week", "replicate", "x2", "x1")
dataOut_ABC <-read.table("./data/7-neb_sim_ebo.csv",sep=",")
colnames(dataOut_ABC) <- c("week", "replicate", "x2", "x3", "x1")
bendertest(dataOut_A, dataOut_AB, dataOut_AC, dataOut_ABC)

#---

dataOut_A <-read.table("./data/1-neb.csv",sep=",")
colnames(dataOut_A) <- c("week", "replicate", "x1")
dataOut_AB <-read.table("./data/4-neb_sim.csv",sep=",")
colnames(dataOut_AB) <- c("week", "replicate", "x1", "x2")
dataOut_AC <-read.table("./data/5-neb_ebo.csv",sep=",")
colnames(dataOut_AC) <- c("week", "replicate", "x1", "x2")
dataOut_ABC <-read.table("./data/7-neb_sim_ebo.csv",sep=",")
colnames(dataOut_ABC) <- c("week", "replicate", "x1", "x2", "x3")
woottontest(dataOut_A, dataOut_AB, dataOut_AC, dataOut_ABC)

dataOut_A <-read.table("./data/2-sim.csv",sep=",")
colnames(dataOut_A) <- c("week", "replicate", "x1")
dataOut_AB <-read.table("./data/4-neb_sim.csv",sep=",")
colnames(dataOut_AB) <- c("week", "replicate", "x2", "x1")
dataOut_AC <-read.table("./data/6-sim_ebo.csv",sep=",")
colnames(dataOut_AC) <- c("week", "replicate", "x1", "x2")
dataOut_ABC <-read.table("./data/7-neb_sim_ebo.csv",sep=",")
colnames(dataOut_ABC) <- c("week", "replicate", "x2", "x1", "x3")
woottontest(dataOut_A, dataOut_AB, dataOut_AC, dataOut_ABC)

dataOut_A <-read.table("./data/3-ebo.csv",sep=",")
colnames(dataOut_A) <- c("week", "replicate", "x1")
dataOut_AB <-read.table("./data/5-neb_ebo.csv",sep=",")
colnames(dataOut_AB) <- c("week", "replicate", "x2", "x1")
dataOut_AC <-read.table("./data/6-sim_ebo.csv",sep=",")
colnames(dataOut_AC) <- c("week", "replicate", "x2", "x1")
dataOut_ABC <-read.table("./data/7-neb_sim_ebo.csv",sep=",")
colnames(dataOut_ABC) <- c("week", "replicate", "x2", "x3", "x1")
woottontest(dataOut_A, dataOut_AB, dataOut_AC, dataOut_ABC)

getSampleHOI <- function(ii) {
  dataSample <-read.table("./data/7-neb_sim_ebo.csv",sep=",")
  colnames(dataSample) <- c("week", "replicate", "x1", "x2", "x3")
  splitted.data <- multisplit(dataSample, c("replicate"))
  dataSample1 <- splitted.data[[ii]]
  dataSample2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0)
  colnames(dataSample2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                             "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                             "x31", "x32", "x33", "x34", "x35", "x36", "x37")
  dataSample2[1, 1:7] <- dataSample1$x1[-8]
  dataSample2[1, 8:14] <- dataSample1$x2[-8]
  dataSample2[1, 15:21] <- dataSample1$x3[-8]
  return(dataSample2)
}

sample1 <- getSampleHOI(1)
sample2 <- getSampleHOI(2)
sample3 <- getSampleHOI(3)

predict(list_glm$fit, sample1)
predict(list_glm$fit, sample2)
predict(list_glm$fit, sample3)

predict(list_knn$fit, sample1)
predict(list_knn$fit, sample2)
predict(list_knn$fit, sample3)

predict(list_svm$fit, sample1)
predict(list_svm$fit, sample2)
predict(list_svm$fit, sample3)