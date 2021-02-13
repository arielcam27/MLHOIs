

get_three_week <- function(file, parSingle, parPair, bmin, bmax)  {
  dataOut <-read.table(file,sep=",")
  colnames(dataOut) <- c("time", "replicate", "x1", "x2", "x3")
  
  splitted.data <- multisplit(dataOut, c("replicate"))

  model <- function(u,par,t) {
    with(as.list(par), {
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
      
      b123 = par[4]
      b231 = par[5]
      b312 = par[6]
      
      x1 = u[1]
      x2 = u[2]
      x3 = u[3]
      
      du1 = r1*x1*(1.0 - x1/K1 + a12*x2/K2 + a13*x3/K3 + b123*x2*x3/(K2*K3))
      du2 = r2*x2*(1.0 - x2/K2 + a23*x3/K3 + a21*x1/K1 + b231*x3*x1/(K3*K1))
      du3 = r3*x3*(1.0 - x3/K3 + a31*x1/K1 + a32*x2/K2 + b312*x1*x2/(K1*K2))
      
      return(c(du1,du2,du3))}
    )
  }
  
  # initial "guess"
  bmid = 0.01 * amin
  
  low <-   c(x10 = 1.0,  x20 = 1.0,  x30 = 1.0,  b123=bmin, b231=bmin, b312=bmin)
  par <-   c(x10 = 5.0,  x20 = 5.0,  x30 = 5.0,  b123=bmid, b231=bmid, b312=bmid)
  up <-    c(x10 = 10.0, x20 = 10.0, x30 = 10.0, b123=bmax, b231=bmax, b312=bmax)
  
  # model cost,
  ModelCost3 <- function(P) {
    tspan <- list(0, 7)
    x10 = P["x10"]
    x20 = P["x20"]
    x30 = P["x30"]
    u0 = c(x10, x20, x30)
    par[4] = P["b123"]
    par[5] = P["b231"]
    par[6] = P["b312"]
    out = diffeqr::ode.solve(f = model, u0 = u0, tspan = tspan, p = par, saveat=1.0)
    out2 = data.frame(time = out$t, x1 = out$u[,1], x2 = out$u[,2], x3 = out$u[,3])
    out2 = out2[-1,]
    colnames(out2) <- c("time", "x1", "x2", "x3")
    return(modCost(out2, Data))  # object of class modCost
  }
  
  df_pars <- data.frame("replicate" = 1,
                        "x10" = 0, 
                        "x20" = 0, 
                        "x30" = 0, 
                        "b123" = 0,
                        "b231" = 0,
                        "b312" = 0,
                        stringsAsFactors = FALSE)
  
  df_pars <- df_pars[-1,]
  
  for (ii in 1:3) {
    Data <- splitted.data[[ii]]
    Data <- Data[,-2]
    Fit <- modFit(f = ModelCost3, 
                  p = par, 
                  lower = low,
                  upper = up)
    
    df_pars[nrow(df_pars) +1,] = c(ii, 
                                   Fit[[1]][[1]], 
                                   Fit[[1]][[2]], 
                                   Fit[[1]][[3]], 
                                   Fit[[1]][[4]], 
                                   Fit[[1]][[5]],
                                   Fit[[1]][[6]])
  }
  
  return(list("pars" = df_pars, "data" = splitted.data))
}