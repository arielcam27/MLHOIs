
get_single_week1 <- function(file)  {

  dataOut <-read.table(file,sep=",")
  colnames(dataOut) <- c("time", "replicate", "x1")

  splitted.data <- multisplit(dataOut, c("replicate"))
  
  model <- function(u,par,t) {
    with(as.list(par), {
      r = par[2]
      K = par[3]
      du = r*u*(1.0 - u/K)
      return(du)
    })
  }
    
  low <- c(y0 = 10.00, r = 0.10, K = 100)
  par <- c(y0 = 25.00, r = 1.00, K = 500)
  upp <- c(y0 = 30.00, r = 10.0, K = 1500)
  
  ModelCost1 <- function(P) {
    tspan <- list(0, 7)
    u0 = P["y0"]
    par[2] = P["r"]
    par[3] = P["K"]
    out = diffeqr::ode.solve(f = model, u0 = u0, tspan = tspan, p = par, saveat=1.0)
    out2 = data.frame(time = out$t, x1 = out$u)
    out2 = out2[-1,]
    colnames(out2) <- c("time", "x1")
    return(modCost(out2, Data))
  }
  
  df_pars <- data.frame("replicate" = 1,
                        "y0" = 0, 
                        "r"  = 0, 
                        "K"  = 0, 
                        stringsAsFactors = FALSE)
  
  df_pars <- df_pars[-1,]

  for (ii in 1:3) {
    Data <- splitted.data[[ii]]
    Data <- Data[,-2]
    Fit <- modFit(f = ModelCost1, lower = low, p = par, upper = upp, method = "Marq")
    df_pars[nrow(df_pars) +1,] = c(0, Fit[[1]][[1]], Fit[[1]][[2]], Fit[[1]][[3]]) 
  }
  
  return(list("pars" = df_pars, "data" = splitted.data))
}
