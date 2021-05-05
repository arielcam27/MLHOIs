get_pairwise_week1 <- function(file, df1, df2, amin, amax)  {
  # Function to estimate parameters from 2-species data.
  # INPUT: csv file, species #1 parameters, species #2 parameters, min value 
  # for aij, and max value for aij.
  # OUTPUT: Dataframe with parameters "p" (x0, y0, a12, a21) and splitted data "data".
  
  dataOut <-read.table(file,sep=",")
  colnames(dataOut) <- c("time", "replicate", "x1", "x2")
  
  splitted.data <- multisplit(dataOut, c("replicate"))
  
  model <- function(u,par,t) {
    with(as.list(par), {
        r1 = df1[1]
        r2 = df2[1]
        K1 = df1[2]
        K2 = df2[2]
        
        a12 = par[3]
        a21 = par[4]
        
        du1 = r1*u[1]*(1 - u[1]/K1 + a12*u[2]/K2)
        du2 = r2*u[2]*(1 - u[2]/K2 + a21*u[1]/K1)
        return(c(du1, du2))
    })
  }
  
  # Hard coded mid value for aij guess.
  amid = 0.01 * amin
  
  # Intervals hard coded. Customize for another problem.
  low <- c(x10 = 5.00, x20 = 5.00, a12 = amin, a21 = amin)
  par <- c(x10 = 15.0, x20 = 15.0, a12 = amid, a21 = amid)
  upp <- c(x10 = 20.0, x20 = 20.0, a12 = amax, a21 = amax)
  
  ModelCost1 <- function(P) {
    # Cost function for model fitting.
    
    tspan <- list(0, 7)
    x10 = P["x10"]
    x20 = P["x20"]
    u0 = c(x10, x20)
    par[3] = P["a12"]
    par[4] = P["a21"]
    out = diffeqr::ode.solve(f = model, u0 = u0, tspan = tspan, p = par, saveat=1.0)
    out2 = data.frame(time = out$t, x1 = out$u[,1], x2 = out$u[,2])
    out2 = out2[-1,]
    colnames(out2) <- c("time", "x1", "x2")
    return(modCost(out2, Data))
  }
  
  df_pars <- data.frame("replicate" = 1,
                        "x10" = 0,
                        "x20" = 0, 
                        "a12"  = 0, 
                        "a21"  = 0, 
                        stringsAsFactors = FALSE)
  
  df_pars <- df_pars[-1,]
  
  for (ii in 1:3) {
    Data <- splitted.data[[ii]]
    Data <- Data[,-2]
    Fit <- modFit(f = ModelCost1, lower = low, p = par, upper = upp, method = "Marq")
    df_pars[nrow(df_pars) +1,] = c(0, Fit[[1]][[1]], Fit[[1]][[2]], Fit[[1]][[3]], Fit[[1]][[4]]) 
  }
  
  return(list("pars" = df_pars, "data" = splitted.data))
}
