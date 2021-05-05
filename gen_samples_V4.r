gen_samples <- function(number, parSingle, parPair, parThree) {
# Function to generate synthetic samples to test models.
# INPUT: number of samples, 1-species parameters, 2-species parameters,
# 3-species paramenters.
# OUTPUT: Dataframes with HOI and NON-HOI samples.
# TODO: Remove hard coded dataframes with noisy samples.
  
  modelHOI= function(u, p, t) {
    # Deterministic HOI model.
    
    r1 = p[[1]]
    r2 = p[[2]]
    r3 = p[[3]]
    
    K1 = p[[4]]
    K2 = p[[5]]
    K3 = p[[6]]
    
    a12 = p[[7]]
    a13 = p[[8]]
    a21 = p[[9]]
    a23 = p[[10]]
    a31 = p[[11]]
    a32 = p[[12]]
    
    b123 = p[[13]]
    b231 = p[[14]]
    b312 = p[[15]]
    
    x1 = u[1]
    x2 = u[2]
    x3 = u[3]
    
    dx1 = r1*x1*(1 - x1/K1 + a12*x2/K2 + a13*x3/K3 + b123*x2*x3/(K2*K3))
    dx2 = r2*x2*(1 - x2/K2 + a23*x3/K3 + a21*x1/K1 + b231*x3*x1/(K3*K1))
    dx3 = r3*x3*(1 - x3/K3 + a31*x1/K1 + a32*x2/K2 + b312*x1*x2/(K1*K2))
    return(c(dx1, dx2, dx3))
  }
  
  modelNONHOI= function(u, p, t) {
    # Deterministic NON-HOI model.
    
    r1 = p[[1]]
    r2 = p[[2]]
    r3 = p[[3]]
    
    K1 = p[[4]]
    K2 = p[[5]]
    K3 = p[[6]]
    
    a12 = p[[7]]
    a13 = p[[8]]
    a21 = p[[9]]
    a23 = p[[10]]
    a31 = p[[11]]
    a32 = p[[12]]
    
    b123 = 0
    b231 = 0
    b312 = 0
    
    x1 = u[1]
    x2 = u[2]
    x3 = u[3]
    
    dx1 = r1*x1*(1 - x1/K1 + a12*x2/K2 + a13*x3/K3 + b123*x2*x3/(K2*K3))
    dx2 = r2*x2*(1 - x2/K2 + a23*x3/K3 + a21*x1/K1 + b231*x3*x1/(K3*K1))
    dx3 = r3*x3*(1 - x3/K3 + a31*x1/K1 + a32*x2/K2 + b312*x1*x2/(K1*K2))
    
    return(c(dx1, dx2, dx3))
  }
  
  stochastic <- function(u, p, t) {
    #return(0.2^2*c(sqrt(u[1]), sqrt(u[2]), sqrt(u[3])))
    return(-0.5^2 * c(u[1], u[2], u[3]))
  }
  
  # Empty dataframes.
  
  dataOutHOI <- data.frame(0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0,
                           "hoi")
  colnames(dataOutHOI) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                            "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                            "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                            "isHOI")
  
  dataOutHOIsingle1SDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0,
                           "hoi")
  colnames(dataOutHOIsingle1SDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                            "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                            "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                            "isHOI")
  dataOutHOIsingle1SDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      "hoi")
  colnames(dataOutHOIsingle1SDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                       "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                       "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                       "isHOI")
  dataOutHOIsingle1SDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      "hoi")
  colnames(dataOutHOIsingle1SDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                       "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                       "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                       "isHOI")
  
  dataOutHOIsingle2SDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      "hoi")
  colnames(dataOutHOIsingle2SDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                       "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                       "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                       "isHOI")
  dataOutHOIsingle2SDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      "hoi")
  colnames(dataOutHOIsingle2SDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                       "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                       "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                       "isHOI")
  dataOutHOIsingle2SDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      "hoi")
  colnames(dataOutHOIsingle2SDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                       "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                       "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                       "isHOI")
  
  dataOutHOIsingle3SDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      "hoi")
  colnames(dataOutHOIsingle3SDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                       "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                       "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                       "isHOI")
  dataOutHOIsingle3SDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      "hoi")
  colnames(dataOutHOIsingle3SDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                       "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                       "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                       "isHOI")
  dataOutHOIsingle3SDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      "hoi")
  colnames(dataOutHOIsingle3SDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                       "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                       "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                       "isHOI")
  
  dataOutHOIpairwise1SDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      "hoi")
  colnames(dataOutHOIpairwise1SDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                       "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                       "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                       "isHOI")
  dataOutHOIpairwise1SDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        "hoi")
  colnames(dataOutHOIpairwise1SDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                         "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                         "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                         "isHOI")
  dataOutHOIpairwise1SDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        "hoi")
  colnames(dataOutHOIpairwise1SDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                         "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                         "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                         "isHOI")
  
  dataOutHOIpairwise2SDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        "hoi")
  colnames(dataOutHOIpairwise2SDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                         "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                         "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                         "isHOI")
  dataOutHOIpairwise2SDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        "hoi")
  colnames(dataOutHOIpairwise2SDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                         "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                         "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                         "isHOI")
  dataOutHOIpairwise2SDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        "hoi")
  colnames(dataOutHOIpairwise2SDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                         "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                         "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                         "isHOI")
  
  dataOutHOIpairwise3SDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        "hoi")
  colnames(dataOutHOIpairwise3SDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                         "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                         "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                         "isHOI")
  dataOutHOIpairwise3SDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        "hoi")
  colnames(dataOutHOIpairwise3SDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                         "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                         "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                         "isHOI")
  dataOutHOIpairwise3SDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        "hoi")
  colnames(dataOutHOIpairwise3SDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                         "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                         "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                         "isHOI")
  
  
  dataOutHOISDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0,
                           "hoi")
  colnames(dataOutHOISDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                "isHOI")
  dataOutHOISDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0,
                               "hoi")
  colnames(dataOutHOISDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                "isHOI")
  dataOutHOISDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0,
                               "hoi")
  colnames(dataOutHOISDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                "isHOI")
  
  dataOutNONHOI <- data.frame(0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0,
                              "nonhoi")
  colnames(dataOutNONHOI) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                               "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                               "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                               "isHOI")
  
  dataOutNONHOIsingle1SDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         "nonhoi")
  colnames(dataOutNONHOIsingle1SDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                          "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                          "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                          "isHOI")
  dataOutNONHOIsingle1SDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         "nonhoi")
  colnames(dataOutNONHOIsingle1SDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                          "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                          "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                          "isHOI")
  dataOutNONHOIsingle1SDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         "nonhoi")
  colnames(dataOutNONHOIsingle1SDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                          "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                          "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                          "isHOI")
  
  dataOutNONHOIsingle2SDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         "nonhoi")
  colnames(dataOutNONHOIsingle2SDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                          "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                          "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                          "isHOI")
  dataOutNONHOIsingle2SDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         "nonhoi")
  colnames(dataOutNONHOIsingle2SDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                          "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                          "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                          "isHOI")
  dataOutNONHOIsingle2SDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         "nonhoi")
  colnames(dataOutNONHOIsingle2SDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                          "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                          "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                          "isHOI")
  
  dataOutNONHOIsingle3SDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         "nonhoi")
  colnames(dataOutNONHOIsingle3SDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                          "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                          "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                          "isHOI")
  dataOutNONHOIsingle3SDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         "nonhoi")
  colnames(dataOutNONHOIsingle3SDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                          "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                          "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                          "isHOI")
  dataOutNONHOIsingle3SDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         "nonhoi")
  colnames(dataOutNONHOIsingle3SDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                          "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                          "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                          "isHOI")
  
  dataOutNONHOIpairwise1SDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           "nonhoi")
  colnames(dataOutNONHOIpairwise1SDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                            "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                            "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                            "isHOI")
  dataOutNONHOIpairwise1SDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           "nonhoi")
  colnames(dataOutNONHOIpairwise1SDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                            "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                            "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                            "isHOI")
  dataOutNONHOIpairwise1SDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           "nonhoi")
  colnames(dataOutNONHOIpairwise1SDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                            "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                            "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                            "isHOI")
  
  dataOutNONHOIpairwise2SDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           "nonhoi")
  colnames(dataOutNONHOIpairwise2SDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                            "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                            "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                            "isHOI")
  dataOutNONHOIpairwise2SDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           "nonhoi")
  colnames(dataOutNONHOIpairwise2SDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                            "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                            "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                            "isHOI")
  dataOutNONHOIpairwise2SDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           "nonhoi")
  colnames(dataOutNONHOIpairwise2SDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                            "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                            "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                            "isHOI")
  
  dataOutNONHOIpairwise3SDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           "nonhoi")
  colnames(dataOutNONHOIpairwise3SDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                            "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                            "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                            "isHOI")
  dataOutNONHOIpairwise3SDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           "nonhoi")
  colnames(dataOutNONHOIpairwise3SDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                            "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                            "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                            "isHOI")
  dataOutNONHOIpairwise3SDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           "nonhoi")
  colnames(dataOutNONHOIpairwise3SDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                            "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                            "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                            "isHOI")
  
  dataOutNONHOISDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0,
                              "nonhoi")
  colnames(dataOutNONHOISDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                   "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                   "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                   "isHOI")
  dataOutNONHOISDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0,
                                  "nonhoi")
  colnames(dataOutNONHOISDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                   "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                   "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                   "isHOI")
  dataOutNONHOISDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0,
                                  "nonhoi")
  colnames(dataOutNONHOISDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                   "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                   "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                   "isHOI")

  # Index for keeping track of number of samples generated.
  index <- 1
  # Tolerance level to accept "sufficiently different" HOI and NON-HOI samples.
  TOL <- 50
  
  while (index < number + 1) {
    # Print progress.
    if (index == 1) { cat("Progress:", index, "\n") }
    if (index %% 100 == 0) { cat("Progress:", index, "\n") }
    
    # Generate random initial conditions (hard coded).
    start = c(x1=runif(1, min=2, max=10), 
              x2=runif(1, min=2, max=10), 
              x3=runif(1, min=2, max=10))
    
    # Prepare initial conditions for 1- and 2-species.
    start1 <- c(start[1], 0, 0)
    start2 <- c(0, start[2], 0)
    start3 <- c(0, 0, start[3])
    start4 <- c(start[1], start[2], 0)
    start5 <- c(start[1], 0, start[3])
    start6 <- c(0, start[2], start[3])
    
    # Generate random model parameters (hard coded).
    r1 = parSingle[1,1] * (1 + runif(1, -0.2, 0.2))
    r2 = parSingle[1,2] * (1 + runif(1, -0.2, 0.2))
    r3 = parSingle[1,3] * (1 + runif(1, -0.2, 0.2))
    
    K1 = parSingle[2,1] * (1 + runif(1, -0.2, 0.2))
    K2 = parSingle[2,2] * (1 + runif(1, -0.2, 0.2))
    K3 = parSingle[2,3] * (1 + runif(1, -0.2, 0.2))
    
    a12 = parPair[1,1] * (1 + runif(1, -0.2, 0.2))
    a21 = parPair[2,1] * (1 + runif(1, -0.2, 0.2))
    
    a23 = parPair[1,3] * (1 + runif(1, -0.2, 0.2))
    a32 = parPair[2,3] * (1 + runif(1, -0.2, 0.2))
    
    a13 = parPair[1,2] * (1 + runif(1, -0.2, 0.2))
    a31 = parPair[2,2] * (1 + runif(1, -0.2, 0.2))  
    
    j1 <- parThree[1] * (1 + runif(1, -1.0, 4.0))
    j2 <- parThree[2] * (1 + runif(1, -1.0, 4.0))
    j3 <- parThree[3] * (1 + runif(1, -1.0, 4.0))
    
    # Check if HOI parameters are "big enough" (hard coded).
    flag1 = abs(j1) > 0.00001
    flag2 = abs(j2) > 0.00001
    flag3 = abs(j3) > 0.00001
    
    if (flag1 | flag2 | flag3) {
      
      b123 <- j1
      b231 <- j2
      b312 <- j3
      
      # Parameter vector p.
      par = c(r1=r1, r2=r2, r3=r3, 
              K1=K1, K2=K2, K3=K3, 
              a12=a12, a13=a13, a21=a21, a23=a23, a31=a31, a32=a32,
              b123=b123, b231=b231, b312=b312)
      
      tryCatch(
        {
          # HOI sample.
          outHOI = diffeqr::ode.solve(f = modelHOI, u0 = start, tspan = list(0, 7), p = par, saveat=1.0)
          out2HOI = data.frame(time = outHOI$t, x1 = outHOI$u[,1], x2 = outHOI$u[,2], x3 = outHOI$u[,3])
          out2HOI = out2HOI[-1,]
          colnames(out2HOI) <- c("time", "x1", "x2", "x3")
          
          # NON-HOI sample.
          outNONHOI = diffeqr::ode.solve(f = modelNONHOI, u0 = start, tspan = list(0, 7), p = par, saveat=1.0)
          out2NONHOI = data.frame(time = outNONHOI$t, x1 = outNONHOI$u[,1], x2 = outNONHOI$u[,2], x3 = outNONHOI$u[,3])
          out2NONHOI = out2NONHOI[-1,]
          colnames(out2NONHOI) <- c("time", "x1", "x2", "x3")
          
          # Check tolerance level condition.
          flag4 = dist(rbind(unlist(out2HOI), unlist(out2NONHOI))) > TOL
          
          if (flag4) {
            # Save HOI sample.
            dataOutHOI[index, 1:7]   <- out2HOI$x1
            dataOutHOI[index, 8:14]  <- out2HOI$x2
            dataOutHOI[index, 15:21] <- out2HOI$x3
            dataOutHOI[index, 22]    <- "hoi"
            # Save NON-HOI sample.
            dataOutNONHOI[index, 1:7]   <- out2NONHOI$x1
            dataOutNONHOI[index, 8:14]  <- out2NONHOI$x2
            dataOutNONHOI[index, 15:21] <- out2NONHOI$x3
            dataOutNONHOI[index, 22]    <- "nonhoi"
            
            # Generate and save HOI noisy sample (3-species, three replicates).
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOISDE1[index, 1:7]   <- df$x1
            dataOutHOISDE1[index, 8:14]  <- df$x2
            dataOutHOISDE1[index, 15:21] <- df$x3
            dataOutHOISDE1[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOISDE2[index, 1:7]   <- df$x1
            dataOutHOISDE2[index, 8:14]  <- df$x2
            dataOutHOISDE2[index, 15:21] <- df$x3
            dataOutHOISDE2[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOISDE3[index, 1:7]   <- df$x1
            dataOutHOISDE3[index, 8:14]  <- df$x2
            dataOutHOISDE3[index, 15:21] <- df$x3
            dataOutHOISDE3[index, 22]    <- "hoi"
            
            # Generate and save HOI noisy sample (1-species, species #1, three replicates).
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start1, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIsingle1SDE1[index, 1:7]   <- df$x1
            dataOutHOIsingle1SDE1[index, 8:14]  <- df$x2
            dataOutHOIsingle1SDE1[index, 15:21] <- df$x3
            dataOutHOIsingle1SDE1[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start1, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIsingle1SDE2[index, 1:7]   <- df$x1
            dataOutHOIsingle1SDE2[index, 8:14]  <- df$x2
            dataOutHOIsingle1SDE2[index, 15:21] <- df$x3
            dataOutHOIsingle1SDE2[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start1, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIsingle1SDE3[index, 1:7]   <- df$x1
            dataOutHOIsingle1SDE3[index, 8:14]  <- df$x2
            dataOutHOIsingle1SDE3[index, 15:21] <- df$x3
            dataOutHOIsingle1SDE3[index, 22]    <- "hoi"
            
            # Generate and save HOI noisy sample (1-species, species #2, three replicates).
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start2, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIsingle2SDE1[index, 1:7]   <- df$x1
            dataOutHOIsingle2SDE1[index, 8:14]  <- df$x2
            dataOutHOIsingle2SDE1[index, 15:21] <- df$x3
            dataOutHOIsingle2SDE1[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start2, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIsingle2SDE2[index, 1:7]   <- df$x1
            dataOutHOIsingle2SDE2[index, 8:14]  <- df$x2
            dataOutHOIsingle2SDE2[index, 15:21] <- df$x3
            dataOutHOIsingle2SDE2[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start2, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIsingle2SDE3[index, 1:7]   <- df$x1
            dataOutHOIsingle2SDE3[index, 8:14]  <- df$x2
            dataOutHOIsingle2SDE3[index, 15:21] <- df$x3
            dataOutHOIsingle2SDE3[index, 22]    <- "hoi"
            
            # Generate and save HOI noisy sample (1-species, species #3, three replicates).
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start3, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIsingle3SDE1[index, 1:7]   <- df$x1
            dataOutHOIsingle3SDE1[index, 8:14]  <- df$x2
            dataOutHOIsingle3SDE1[index, 15:21] <- df$x3
            dataOutHOIsingle3SDE1[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start3, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIsingle3SDE2[index, 1:7]   <- df$x1
            dataOutHOIsingle3SDE2[index, 8:14]  <- df$x2
            dataOutHOIsingle3SDE2[index, 15:21] <- df$x3
            dataOutHOIsingle3SDE2[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start3, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIsingle3SDE3[index, 1:7]   <- df$x1
            dataOutHOIsingle3SDE3[index, 8:14]  <- df$x2
            dataOutHOIsingle3SDE3[index, 15:21] <- df$x3
            dataOutHOIsingle3SDE3[index, 22]    <- "hoi"
            
            # Generate and save HOI noisy sample (2-species, species #1 vs #2, three replicates).
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start4, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIpairwise1SDE1[index, 1:7]   <- df$x1
            dataOutHOIpairwise1SDE1[index, 8:14]  <- df$x2
            dataOutHOIpairwise1SDE1[index, 15:21] <- df$x3
            dataOutHOIpairwise1SDE1[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start4, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIpairwise1SDE2[index, 1:7]   <- df$x1
            dataOutHOIpairwise1SDE2[index, 8:14]  <- df$x2
            dataOutHOIpairwise1SDE2[index, 15:21] <- df$x3
            dataOutHOIpairwise1SDE2[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start4, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIpairwise1SDE3[index, 1:7]   <- df$x1
            dataOutHOIpairwise1SDE3[index, 8:14]  <- df$x2
            dataOutHOIpairwise1SDE3[index, 15:21] <- df$x3
            dataOutHOIpairwise1SDE3[index, 22]    <- "hoi"
            
            # Generate and save HOI noisy sample (2-species, species #1 vs #3, three replicates).
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start5, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIpairwise2SDE1[index, 1:7]   <- df$x1
            dataOutHOIpairwise2SDE1[index, 8:14]  <- df$x2
            dataOutHOIpairwise2SDE1[index, 15:21] <- df$x3
            dataOutHOIpairwise2SDE1[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start5, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIpairwise2SDE2[index, 1:7]   <- df$x1
            dataOutHOIpairwise2SDE2[index, 8:14]  <- df$x2
            dataOutHOIpairwise2SDE2[index, 15:21] <- df$x3
            dataOutHOIpairwise2SDE2[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start5, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIpairwise2SDE3[index, 1:7]   <- df$x1
            dataOutHOIpairwise2SDE3[index, 8:14]  <- df$x2
            dataOutHOIpairwise2SDE3[index, 15:21] <- df$x3
            dataOutHOIpairwise2SDE3[index, 22]    <- "hoi"
            
            # Generate and save HOI noisy sample (2-species, species #2 vs #3, three replicates).
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start6, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIpairwise3SDE1[index, 1:7]   <- df$x1
            dataOutHOIpairwise3SDE1[index, 8:14]  <- df$x2
            dataOutHOIpairwise3SDE1[index, 15:21] <- df$x3
            dataOutHOIpairwise3SDE1[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start6, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIpairwise3SDE2[index, 1:7]   <- df$x1
            dataOutHOIpairwise3SDE2[index, 8:14]  <- df$x2
            dataOutHOIpairwise3SDE2[index, 15:21] <- df$x3
            dataOutHOIpairwise3SDE2[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start6, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIpairwise3SDE3[index, 1:7]   <- df$x1
            dataOutHOIpairwise3SDE3[index, 8:14]  <- df$x2
            dataOutHOIpairwise3SDE3[index, 15:21] <- df$x3
            dataOutHOIpairwise3SDE3[index, 22]    <- "hoi"
            
            # Generate NON-HOI samples.
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start, tspan = list(0, 7), p = par, saveat=1.0)
            out2NONHOISDE1 = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            out2NONHOISDE1 = out2NONHOISDE1[-1,]
            colnames(out2NONHOISDE1) <- c("time", "x1", "x2", "x3")
            dataOutNONHOISDE1[index, 1:7]   <- out2NONHOISDE1$x1
            dataOutNONHOISDE1[index, 8:14]  <- out2NONHOISDE1$x2
            dataOutNONHOISDE1[index, 15:21] <- out2NONHOISDE1$x3
            dataOutNONHOISDE1[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start, tspan = list(0, 7), p = par, saveat=1.0)
            out2NONHOISDE2 = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            out2NONHOISDE2 = out2NONHOISDE2[-1,]
            colnames(out2NONHOISDE2) <- c("time", "x1", "x2", "x3")
            dataOutNONHOISDE2[index, 1:7]   <- out2NONHOISDE2$x1
            dataOutNONHOISDE2[index, 8:14]  <- out2NONHOISDE2$x2
            dataOutNONHOISDE2[index, 15:21] <- out2NONHOISDE2$x3
            dataOutNONHOISDE2[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start, tspan = list(0, 7), p = par, saveat=1.0)
            out2NONHOISDE3 = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            out2NONHOISDE3 = out2NONHOISDE3[-1,]
            colnames(out2NONHOISDE3) <- c("time", "x1", "x2", "x3")
            dataOutNONHOISDE3[index, 1:7]   <- out2NONHOISDE3$x1
            dataOutNONHOISDE3[index, 8:14]  <- out2NONHOISDE3$x2
            dataOutNONHOISDE3[index, 15:21] <- out2NONHOISDE3$x3
            dataOutNONHOISDE3[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start1, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIsingle1SDE1[index, 1:7]   <- df$x1
            dataOutNONHOIsingle1SDE1[index, 8:14]  <- df$x2
            dataOutNONHOIsingle1SDE1[index, 15:21] <- df$x3
            dataOutNONHOIsingle1SDE1[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start1, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIsingle1SDE2[index, 1:7]   <- df$x1
            dataOutNONHOIsingle1SDE2[index, 8:14]  <- df$x2
            dataOutNONHOIsingle1SDE2[index, 15:21] <- df$x3
            dataOutNONHOIsingle1SDE2[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start1, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIsingle1SDE3[index, 1:7]   <- df$x1
            dataOutNONHOIsingle1SDE3[index, 8:14]  <- df$x2
            dataOutNONHOIsingle1SDE3[index, 15:21] <- df$x3
            dataOutNONHOIsingle1SDE3[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start2, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIsingle2SDE1[index, 1:7]   <- df$x1
            dataOutNONHOIsingle2SDE1[index, 8:14]  <- df$x2
            dataOutNONHOIsingle2SDE1[index, 15:21] <- df$x3
            dataOutNONHOIsingle2SDE1[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start2, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIsingle2SDE2[index, 1:7]   <- df$x1
            dataOutNONHOIsingle2SDE2[index, 8:14]  <- df$x2
            dataOutNONHOIsingle2SDE2[index, 15:21] <- df$x3
            dataOutNONHOIsingle2SDE2[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start2, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIsingle2SDE3[index, 1:7]   <- df$x1
            dataOutNONHOIsingle2SDE3[index, 8:14]  <- df$x2
            dataOutNONHOIsingle2SDE3[index, 15:21] <- df$x3
            dataOutNONHOIsingle2SDE3[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start3, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIsingle3SDE1[index, 1:7]   <- df$x1
            dataOutNONHOIsingle3SDE1[index, 8:14]  <- df$x2
            dataOutNONHOIsingle3SDE1[index, 15:21] <- df$x3
            dataOutNONHOIsingle3SDE1[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start3, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIsingle3SDE2[index, 1:7]   <- df$x1
            dataOutNONHOIsingle3SDE2[index, 8:14]  <- df$x2
            dataOutNONHOIsingle3SDE2[index, 15:21] <- df$x3
            dataOutNONHOIsingle3SDE2[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start3, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIsingle3SDE3[index, 1:7]   <- df$x1
            dataOutNONHOIsingle3SDE3[index, 8:14]  <- df$x2
            dataOutNONHOIsingle3SDE3[index, 15:21] <- df$x3
            dataOutNONHOIsingle3SDE3[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start4, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIpairwise1SDE1[index, 1:7]   <- df$x1
            dataOutNONHOIpairwise1SDE1[index, 8:14]  <- df$x2
            dataOutNONHOIpairwise1SDE1[index, 15:21] <- df$x3
            dataOutNONHOIpairwise1SDE1[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start4, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIpairwise1SDE2[index, 1:7]   <- df$x1
            dataOutNONHOIpairwise1SDE2[index, 8:14]  <- df$x2
            dataOutNONHOIpairwise1SDE2[index, 15:21] <- df$x3
            dataOutNONHOIpairwise1SDE2[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start4, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIpairwise1SDE3[index, 1:7]   <- df$x1
            dataOutNONHOIpairwise1SDE3[index, 8:14]  <- df$x2
            dataOutNONHOIpairwise1SDE3[index, 15:21] <- df$x3
            dataOutNONHOIpairwise1SDE3[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start5, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIpairwise2SDE1[index, 1:7]   <- df$x1
            dataOutNONHOIpairwise2SDE1[index, 8:14]  <- df$x2
            dataOutNONHOIpairwise2SDE1[index, 15:21] <- df$x3
            dataOutNONHOIpairwise2SDE1[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start5, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIpairwise2SDE2[index, 1:7]   <- df$x1
            dataOutNONHOIpairwise2SDE2[index, 8:14]  <- df$x2
            dataOutNONHOIpairwise2SDE2[index, 15:21] <- df$x3
            dataOutNONHOIpairwise2SDE2[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start5, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIpairwise2SDE3[index, 1:7]   <- df$x1
            dataOutNONHOIpairwise2SDE3[index, 8:14]  <- df$x2
            dataOutNONHOIpairwise2SDE3[index, 15:21] <- df$x3
            dataOutNONHOIpairwise2SDE3[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start6, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIpairwise3SDE1[index, 1:7]   <- df$x1
            dataOutNONHOIpairwise3SDE1[index, 8:14]  <- df$x2
            dataOutNONHOIpairwise3SDE1[index, 15:21] <- df$x3
            dataOutNONHOIpairwise3SDE1[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start6, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIpairwise3SDE2[index, 1:7]   <- df$x1
            dataOutNONHOIpairwise3SDE2[index, 8:14]  <- df$x2
            dataOutNONHOIpairwise3SDE2[index, 15:21] <- df$x3
            dataOutNONHOIpairwise3SDE2[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start6, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIpairwise3SDE3[index, 1:7]   <- df$x1
            dataOutNONHOIpairwise3SDE3[index, 8:14]  <- df$x2
            dataOutNONHOIpairwise3SDE3[index, 15:21] <- df$x3
            dataOutNONHOIpairwise3SDE3[index, 22]    <- "nonhoi"

            index <- index + 1
          }
        },
        error=function(cond) {
          #index <- index + 1
          #print("error?")
          #message("ERROR: Possibly NaN")
        }
      )
      
      
    }
  }
  return(list(df_HOI = dataOutHOI, 
              df_HOI1SDE1 = dataOutHOIsingle1SDE1,
              df_HOI1SDE2 = dataOutHOIsingle1SDE2,
              df_HOI1SDE3 = dataOutHOIsingle1SDE3,
              df_HOI2SDE1 = dataOutHOIsingle2SDE1,
              df_HOI2SDE2 = dataOutHOIsingle2SDE2,
              df_HOI2SDE3 = dataOutHOIsingle2SDE3,
              df_HOI3SDE1 = dataOutHOIsingle3SDE1,
              df_HOI3SDE2 = dataOutHOIsingle3SDE2,
              df_HOI3SDE3 = dataOutHOIsingle3SDE3,
              df_HOI12SDE1 = dataOutHOIpairwise1SDE1,
              df_HOI12SDE2 = dataOutHOIpairwise1SDE2,
              df_HOI12SDE3 = dataOutHOIpairwise1SDE3,
              df_HOI13SDE1 = dataOutHOIpairwise2SDE1,
              df_HOI13SDE2 = dataOutHOIpairwise2SDE2,
              df_HOI13SDE3 = dataOutHOIpairwise2SDE3,
              df_HOI23SDE1 = dataOutHOIpairwise3SDE1,
              df_HOI23SDE2 = dataOutHOIpairwise3SDE2,
              df_HOI23SDE3 = dataOutHOIpairwise3SDE3,
              df_HOI123SDE1 = dataOutHOISDE1,
              df_HOI123SDE2 = dataOutHOISDE2,
              df_HOI123SDE3 = dataOutHOISDE3,
              df_NONHOI = dataOutNONHOI, 
              df_NONHOI1SDE1 = dataOutNONHOIsingle1SDE1,
              df_NONHOI1SDE2 = dataOutNONHOIsingle1SDE2,
              df_NONHOI1SDE3 = dataOutNONHOIsingle1SDE3,
              df_NONHOI2SDE1 = dataOutNONHOIsingle2SDE1,
              df_NONHOI2SDE2 = dataOutNONHOIsingle2SDE2,
              df_NONHOI2SDE3 = dataOutNONHOIsingle2SDE3,
              df_NONHOI3SDE1 = dataOutNONHOIsingle3SDE1,
              df_NONHOI3SDE2 = dataOutNONHOIsingle3SDE2,
              df_NONHOI3SDE3 = dataOutNONHOIsingle3SDE3,
              df_NONHOI12SDE1 = dataOutNONHOIpairwise1SDE1,
              df_NONHOI12SDE2 = dataOutNONHOIpairwise1SDE2,
              df_NONHOI12SDE3 = dataOutNONHOIpairwise1SDE3,
              df_NONHOI13SDE1 = dataOutNONHOIpairwise2SDE1,
              df_NONHOI13SDE2 = dataOutNONHOIpairwise2SDE2,
              df_NONHOI13SDE3 = dataOutNONHOIpairwise2SDE3,
              df_NONHOI23SDE1 = dataOutNONHOIpairwise3SDE1,
              df_NONHOI23SDE2 = dataOutNONHOIpairwise3SDE2,
              df_NONHOI23SDE3 = dataOutNONHOIpairwise3SDE3,
              df_NONHOI123SDE1 = dataOutNONHOISDE1,
              df_NONHOI123SDE2 = dataOutNONHOISDE2,
              df_NONHOI123SDE3 = dataOutNONHOISDE3
              ))
}



















gen_samples_V2 <- function(number, noise, parSingle, parPair, parThree) {
  # Function to generate synthetic samples to test models.
  # INPUT: number of samples, 1-species parameters, 2-species parameters,
  # 3-species paramenters.
  # OUTPUT: Dataframes with HOI and NON-HOI samples.
  # TODO: Remove hard coded dataframes with noisy samples.
  
  modelHOI= function(u, p, t) {
    # Deterministic HOI model.
    
    r1 = p[[1]]
    r2 = p[[2]]
    r3 = p[[3]]
    K1 = p[[4]]
    K2 = p[[5]]
    K3 = p[[6]]
    a12 = p[[7]]
    a13 = p[[8]]
    a21 = p[[9]]
    a23 = p[[10]]
    a31 = p[[11]]
    a32 = p[[12]]
    b123 = p[[13]]
    b231 = p[[14]]
    b312 = p[[15]]
    
    x1 = u[1]
    x2 = u[2]
    x3 = u[3]
    
    dx1 = r1*x1*(1 - x1/K1 + a12*x2/K2 + a13*x3/K3 + b123*x2*x3/(K2*K3))
    dx2 = r2*x2*(1 - x2/K2 + a23*x3/K3 + a21*x1/K1 + b231*x3*x1/(K3*K1))
    dx3 = r3*x3*(1 - x3/K3 + a31*x1/K1 + a32*x2/K2 + b312*x1*x2/(K1*K2))
    return(c(dx1, dx2, dx3))
  }
  
  modelNONHOI= function(u, p, t) {
    # Deterministic NON-HOI model.
    
    r1 = p[[1]]
    r2 = p[[2]]
    r3 = p[[3]]
    K1 = p[[4]]
    K2 = p[[5]]
    K3 = p[[6]]
    a12 = p[[7]]
    a13 = p[[8]]
    a21 = p[[9]]
    a23 = p[[10]]
    a31 = p[[11]]
    a32 = p[[12]]
    
    x1 = u[1]
    x2 = u[2]
    x3 = u[3]
    
    dx1 = r1*x1*(1 - x1/K1 + a12*x2/K2 + a13*x3/K3)
    dx2 = r2*x2*(1 - x2/K2 + a23*x3/K3 + a21*x1/K1)
    dx3 = r3*x3*(1 - x3/K3 + a31*x1/K1 + a32*x2/K2)
    
    return(c(dx1, dx2, dx3))
  }
  
  stochastic <- function(u, p, t) {
    #return(0.2^2*c(sqrt(u[1]), sqrt(u[2]), sqrt(u[3])))
    return(- noise^2 * c(u[1], u[2], u[3]))
  }
  
  # Empty dataframes.
  
  dataOutHOI <- data.frame(0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0,
                           "hoi")
  colnames(dataOutHOI) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                            "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                            "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                            "isHOI")
  
  dataOutHOIsingle1SDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      "hoi")
  colnames(dataOutHOIsingle1SDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                       "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                       "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                       "isHOI")
  dataOutHOIsingle1SDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      "hoi")
  colnames(dataOutHOIsingle1SDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                       "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                       "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                       "isHOI")
  dataOutHOIsingle1SDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      "hoi")
  colnames(dataOutHOIsingle1SDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                       "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                       "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                       "isHOI")
  
  dataOutHOIsingle2SDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      "hoi")
  colnames(dataOutHOIsingle2SDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                       "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                       "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                       "isHOI")
  dataOutHOIsingle2SDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      "hoi")
  colnames(dataOutHOIsingle2SDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                       "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                       "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                       "isHOI")
  dataOutHOIsingle2SDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      "hoi")
  colnames(dataOutHOIsingle2SDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                       "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                       "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                       "isHOI")
  
  dataOutHOIsingle3SDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      "hoi")
  colnames(dataOutHOIsingle3SDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                       "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                       "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                       "isHOI")
  dataOutHOIsingle3SDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      "hoi")
  colnames(dataOutHOIsingle3SDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                       "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                       "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                       "isHOI")
  dataOutHOIsingle3SDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0,
                                      "hoi")
  colnames(dataOutHOIsingle3SDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                       "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                       "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                       "isHOI")
  
  dataOutHOIpairwise1SDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        "hoi")
  colnames(dataOutHOIpairwise1SDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                         "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                         "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                         "isHOI")
  dataOutHOIpairwise1SDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        "hoi")
  colnames(dataOutHOIpairwise1SDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                         "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                         "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                         "isHOI")
  dataOutHOIpairwise1SDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        "hoi")
  colnames(dataOutHOIpairwise1SDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                         "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                         "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                         "isHOI")
  
  dataOutHOIpairwise2SDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        "hoi")
  colnames(dataOutHOIpairwise2SDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                         "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                         "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                         "isHOI")
  dataOutHOIpairwise2SDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        "hoi")
  colnames(dataOutHOIpairwise2SDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                         "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                         "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                         "isHOI")
  dataOutHOIpairwise2SDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        "hoi")
  colnames(dataOutHOIpairwise2SDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                         "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                         "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                         "isHOI")
  
  dataOutHOIpairwise3SDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        "hoi")
  colnames(dataOutHOIpairwise3SDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                         "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                         "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                         "isHOI")
  dataOutHOIpairwise3SDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        "hoi")
  colnames(dataOutHOIpairwise3SDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                         "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                         "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                         "isHOI")
  dataOutHOIpairwise3SDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0,
                                        "hoi")
  colnames(dataOutHOIpairwise3SDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                         "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                         "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                         "isHOI")
  dataOutHOISDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0,
                               "hoi")
  colnames(dataOutHOISDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                "isHOI")
  dataOutHOISDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0,
                               "hoi")
  colnames(dataOutHOISDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                "isHOI")
  dataOutHOISDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0,
                               "hoi")
  colnames(dataOutHOISDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                "isHOI")
  
  dataOutNONHOI <- data.frame(0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0,
                              "nonhoi")
  colnames(dataOutNONHOI) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                               "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                               "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                               "isHOI")
  
  dataOutNONHOIsingle1SDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         "nonhoi")
  colnames(dataOutNONHOIsingle1SDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                          "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                          "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                          "isHOI")
  dataOutNONHOIsingle1SDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         "nonhoi")
  colnames(dataOutNONHOIsingle1SDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                          "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                          "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                          "isHOI")
  dataOutNONHOIsingle1SDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         "nonhoi")
  colnames(dataOutNONHOIsingle1SDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                          "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                          "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                          "isHOI")
  
  dataOutNONHOIsingle2SDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         "nonhoi")
  colnames(dataOutNONHOIsingle2SDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                          "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                          "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                          "isHOI")
  dataOutNONHOIsingle2SDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         "nonhoi")
  colnames(dataOutNONHOIsingle2SDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                          "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                          "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                          "isHOI")
  dataOutNONHOIsingle2SDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         "nonhoi")
  colnames(dataOutNONHOIsingle2SDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                          "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                          "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                          "isHOI")
  
  dataOutNONHOIsingle3SDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         "nonhoi")
  colnames(dataOutNONHOIsingle3SDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                          "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                          "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                          "isHOI")
  dataOutNONHOIsingle3SDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         "nonhoi")
  colnames(dataOutNONHOIsingle3SDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                          "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                          "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                          "isHOI")
  dataOutNONHOIsingle3SDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0,
                                         "nonhoi")
  colnames(dataOutNONHOIsingle3SDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                          "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                          "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                          "isHOI")
  
  dataOutNONHOIpairwise1SDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           "nonhoi")
  colnames(dataOutNONHOIpairwise1SDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                            "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                            "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                            "isHOI")
  dataOutNONHOIpairwise1SDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           "nonhoi")
  colnames(dataOutNONHOIpairwise1SDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                            "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                            "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                            "isHOI")
  dataOutNONHOIpairwise1SDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           "nonhoi")
  colnames(dataOutNONHOIpairwise1SDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                            "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                            "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                            "isHOI")
  
  dataOutNONHOIpairwise2SDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           "nonhoi")
  colnames(dataOutNONHOIpairwise2SDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                            "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                            "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                            "isHOI")
  dataOutNONHOIpairwise2SDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           "nonhoi")
  colnames(dataOutNONHOIpairwise2SDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                            "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                            "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                            "isHOI")
  dataOutNONHOIpairwise2SDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           "nonhoi")
  colnames(dataOutNONHOIpairwise2SDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                            "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                            "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                            "isHOI")
  
  dataOutNONHOIpairwise3SDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           "nonhoi")
  colnames(dataOutNONHOIpairwise3SDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                            "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                            "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                            "isHOI")
  dataOutNONHOIpairwise3SDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           "nonhoi")
  colnames(dataOutNONHOIpairwise3SDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                            "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                            "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                            "isHOI")
  dataOutNONHOIpairwise3SDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0,
                                           "nonhoi")
  colnames(dataOutNONHOIpairwise3SDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                            "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                            "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                            "isHOI")
  
  dataOutNONHOISDE1 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0,
                                  "nonhoi")
  colnames(dataOutNONHOISDE1) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                   "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                   "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                   "isHOI")
  dataOutNONHOISDE2 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0,
                                  "nonhoi")
  colnames(dataOutNONHOISDE2) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                   "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                   "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                   "isHOI")
  dataOutNONHOISDE3 <- data.frame(0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0,
                                  "nonhoi")
  colnames(dataOutNONHOISDE3) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                   "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                                   "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                                   "isHOI")
  
  # Index for keeping track of number of samples generated.
  index <- 1
  # Tolerance level to accept "sufficiently different" HOI and NON-HOI samples.
  TOL <- 1
  
  while (index < number + 1) {
    # Print progress.
    if (index == 1) { cat("Progress:", index, "\n") }
    if (index %% 10 == 0) { cat("Progress:", index, "\n") }
    
    # Generate random initial conditions (hard coded).
    start = c(x1=runif(1, min=2, max=10), 
              x2=runif(1, min=2, max=10), 
              x3=runif(1, min=2, max=10))
    
    # Prepare initial conditions for 1- and 2-species.
    start1 <- c(start[1], 0, 0)
    start2 <- c(0, start[2], 0)
    start3 <- c(0, 0, start[3])
    start4 <- c(start[1], start[2], 0)
    start5 <- c(start[1], 0, start[3])
    start6 <- c(0, start[2], start[3])
    
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
    
    # Generate random bijk parameters (hard coded).
    j1 <- parThree[1] * runif(1, -5.0, 5.0)
    j2 <- parThree[2] * runif(1, -5.0, 5.0)
    j3 <- parThree[3] * runif(1, -5.0, 5.0)
    
    # Check if HOI parameters are "big enough" (hard coded).
    flag1 = abs(j1) > 0.001
    flag2 = abs(j2) > 0.001
    flag3 = abs(j3) > 0.001
    
    if (flag1 | flag2 | flag3) {
      
      b123 <- j1
      b231 <- j2
      b312 <- j3
      
      # Parameter vector par.
      par = c(r1=r1, r2=r2, r3=r3, 
              K1=K1, K2=K2, K3=K3, 
              a12=a12, a13=a13, a21=a21, a23=a23, a31=a31, a32=a32,
              b123=b123, b231=b231, b312=b312)
      
      tryCatch(
        {
          # HOI sample.
          outHOI = diffeqr::ode.solve(f = modelHOI, u0 = start, tspan = list(0, 7), p = par, saveat=1.0)
          out2HOI = data.frame(time = outHOI$t, x1 = outHOI$u[,1], x2 = outHOI$u[,2], x3 = outHOI$u[,3])
          out2HOI = out2HOI[-1,]
          colnames(out2HOI) <- c("time", "x1", "x2", "x3")
          
          # NON-HOI sample.
          outNONHOI = diffeqr::ode.solve(f = modelNONHOI, u0 = start, tspan = list(0, 7), p = par, saveat=1.0)
          out2NONHOI = data.frame(time = outNONHOI$t, x1 = outNONHOI$u[,1], x2 = outNONHOI$u[,2], x3 = outNONHOI$u[,3])
          out2NONHOI = out2NONHOI[-1,]
          colnames(out2NONHOI) <- c("time", "x1", "x2", "x3")
          
          # Check tolerance level condition.
          flag4 = dist(rbind(unlist(out2HOI), unlist(out2NONHOI))) > TOL
          
          if (flag4) {
            # Save HOI sample.
            dataOutHOI[index, 1:7]   <- out2HOI$x1
            dataOutHOI[index, 8:14]  <- out2HOI$x2
            dataOutHOI[index, 15:21] <- out2HOI$x3
            dataOutHOI[index, 22]    <- "hoi"
            # Save NON-HOI sample.
            dataOutNONHOI[index, 1:7]   <- out2NONHOI$x1
            dataOutNONHOI[index, 8:14]  <- out2NONHOI$x2
            dataOutNONHOI[index, 15:21] <- out2NONHOI$x3
            dataOutNONHOI[index, 22]    <- "nonhoi"
            # Generate and save HOI noisy sample (3-species, three replicates).
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOISDE1[index, 1:7]   <- df$x1
            dataOutHOISDE1[index, 8:14]  <- df$x2
            dataOutHOISDE1[index, 15:21] <- df$x3
            dataOutHOISDE1[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOISDE2[index, 1:7]   <- df$x1
            dataOutHOISDE2[index, 8:14]  <- df$x2
            dataOutHOISDE2[index, 15:21] <- df$x3
            dataOutHOISDE2[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOISDE3[index, 1:7]   <- df$x1
            dataOutHOISDE3[index, 8:14]  <- df$x2
            dataOutHOISDE3[index, 15:21] <- df$x3
            dataOutHOISDE3[index, 22]    <- "hoi"
            # Generate and save HOI noisy sample (1-species, species #1, three replicates).
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start1, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIsingle1SDE1[index, 1:7]   <- df$x1
            dataOutHOIsingle1SDE1[index, 8:14]  <- df$x2
            dataOutHOIsingle1SDE1[index, 15:21] <- df$x3
            dataOutHOIsingle1SDE1[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start1, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIsingle1SDE2[index, 1:7]   <- df$x1
            dataOutHOIsingle1SDE2[index, 8:14]  <- df$x2
            dataOutHOIsingle1SDE2[index, 15:21] <- df$x3
            dataOutHOIsingle1SDE2[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start1, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIsingle1SDE3[index, 1:7]   <- df$x1
            dataOutHOIsingle1SDE3[index, 8:14]  <- df$x2
            dataOutHOIsingle1SDE3[index, 15:21] <- df$x3
            dataOutHOIsingle1SDE3[index, 22]    <- "hoi"
            # Generate and save HOI noisy sample (1-species, species #2, three replicates).
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start2, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIsingle2SDE1[index, 1:7]   <- df$x1
            dataOutHOIsingle2SDE1[index, 8:14]  <- df$x2
            dataOutHOIsingle2SDE1[index, 15:21] <- df$x3
            dataOutHOIsingle2SDE1[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start2, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIsingle2SDE2[index, 1:7]   <- df$x1
            dataOutHOIsingle2SDE2[index, 8:14]  <- df$x2
            dataOutHOIsingle2SDE2[index, 15:21] <- df$x3
            dataOutHOIsingle2SDE2[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start2, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIsingle2SDE3[index, 1:7]   <- df$x1
            dataOutHOIsingle2SDE3[index, 8:14]  <- df$x2
            dataOutHOIsingle2SDE3[index, 15:21] <- df$x3
            dataOutHOIsingle2SDE3[index, 22]    <- "hoi"
            # Generate and save HOI noisy sample (1-species, species #3, three replicates).
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start3, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIsingle3SDE1[index, 1:7]   <- df$x1
            dataOutHOIsingle3SDE1[index, 8:14]  <- df$x2
            dataOutHOIsingle3SDE1[index, 15:21] <- df$x3
            dataOutHOIsingle3SDE1[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start3, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIsingle3SDE2[index, 1:7]   <- df$x1
            dataOutHOIsingle3SDE2[index, 8:14]  <- df$x2
            dataOutHOIsingle3SDE2[index, 15:21] <- df$x3
            dataOutHOIsingle3SDE2[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start3, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIsingle3SDE3[index, 1:7]   <- df$x1
            dataOutHOIsingle3SDE3[index, 8:14]  <- df$x2
            dataOutHOIsingle3SDE3[index, 15:21] <- df$x3
            dataOutHOIsingle3SDE3[index, 22]    <- "hoi"
            # Generate and save HOI noisy sample (2-species, species #1 vs #2, three replicates).
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start4, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIpairwise1SDE1[index, 1:7]   <- df$x1
            dataOutHOIpairwise1SDE1[index, 8:14]  <- df$x2
            dataOutHOIpairwise1SDE1[index, 15:21] <- df$x3
            dataOutHOIpairwise1SDE1[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start4, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIpairwise1SDE2[index, 1:7]   <- df$x1
            dataOutHOIpairwise1SDE2[index, 8:14]  <- df$x2
            dataOutHOIpairwise1SDE2[index, 15:21] <- df$x3
            dataOutHOIpairwise1SDE2[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start4, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIpairwise1SDE3[index, 1:7]   <- df$x1
            dataOutHOIpairwise1SDE3[index, 8:14]  <- df$x2
            dataOutHOIpairwise1SDE3[index, 15:21] <- df$x3
            dataOutHOIpairwise1SDE3[index, 22]    <- "hoi"
            # Generate and save HOI noisy sample (2-species, species #1 vs #3, three replicates).
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start5, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIpairwise2SDE1[index, 1:7]   <- df$x1
            dataOutHOIpairwise2SDE1[index, 8:14]  <- df$x2
            dataOutHOIpairwise2SDE1[index, 15:21] <- df$x3
            dataOutHOIpairwise2SDE1[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start5, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIpairwise2SDE2[index, 1:7]   <- df$x1
            dataOutHOIpairwise2SDE2[index, 8:14]  <- df$x2
            dataOutHOIpairwise2SDE2[index, 15:21] <- df$x3
            dataOutHOIpairwise2SDE2[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start5, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIpairwise2SDE3[index, 1:7]   <- df$x1
            dataOutHOIpairwise2SDE3[index, 8:14]  <- df$x2
            dataOutHOIpairwise2SDE3[index, 15:21] <- df$x3
            dataOutHOIpairwise2SDE3[index, 22]    <- "hoi"
            # Generate and save HOI noisy sample (2-species, species #2 vs #3, three replicates).
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start6, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIpairwise3SDE1[index, 1:7]   <- df$x1
            dataOutHOIpairwise3SDE1[index, 8:14]  <- df$x2
            dataOutHOIpairwise3SDE1[index, 15:21] <- df$x3
            dataOutHOIpairwise3SDE1[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start6, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIpairwise3SDE2[index, 1:7]   <- df$x1
            dataOutHOIpairwise3SDE2[index, 8:14]  <- df$x2
            dataOutHOIpairwise3SDE2[index, 15:21] <- df$x3
            dataOutHOIpairwise3SDE2[index, 22]    <- "hoi"
            
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start6, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutHOIpairwise3SDE3[index, 1:7]   <- df$x1
            dataOutHOIpairwise3SDE3[index, 8:14]  <- df$x2
            dataOutHOIpairwise3SDE3[index, 15:21] <- df$x3
            dataOutHOIpairwise3SDE3[index, 22]    <- "hoi"
            # Generate NON-HOI samples.
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start, tspan = list(0, 7), p = par, saveat=1.0)
            out2NONHOISDE1 = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            out2NONHOISDE1 = out2NONHOISDE1[-1,]
            colnames(out2NONHOISDE1) <- c("time", "x1", "x2", "x3")
            dataOutNONHOISDE1[index, 1:7]   <- out2NONHOISDE1$x1
            dataOutNONHOISDE1[index, 8:14]  <- out2NONHOISDE1$x2
            dataOutNONHOISDE1[index, 15:21] <- out2NONHOISDE1$x3
            dataOutNONHOISDE1[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start, tspan = list(0, 7), p = par, saveat=1.0)
            out2NONHOISDE2 = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            out2NONHOISDE2 = out2NONHOISDE2[-1,]
            colnames(out2NONHOISDE2) <- c("time", "x1", "x2", "x3")
            dataOutNONHOISDE2[index, 1:7]   <- out2NONHOISDE2$x1
            dataOutNONHOISDE2[index, 8:14]  <- out2NONHOISDE2$x2
            dataOutNONHOISDE2[index, 15:21] <- out2NONHOISDE2$x3
            dataOutNONHOISDE2[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start, tspan = list(0, 7), p = par, saveat=1.0)
            out2NONHOISDE3 = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            out2NONHOISDE3 = out2NONHOISDE3[-1,]
            colnames(out2NONHOISDE3) <- c("time", "x1", "x2", "x3")
            dataOutNONHOISDE3[index, 1:7]   <- out2NONHOISDE3$x1
            dataOutNONHOISDE3[index, 8:14]  <- out2NONHOISDE3$x2
            dataOutNONHOISDE3[index, 15:21] <- out2NONHOISDE3$x3
            dataOutNONHOISDE3[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start1, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIsingle1SDE1[index, 1:7]   <- df$x1
            dataOutNONHOIsingle1SDE1[index, 8:14]  <- df$x2
            dataOutNONHOIsingle1SDE1[index, 15:21] <- df$x3
            dataOutNONHOIsingle1SDE1[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start1, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIsingle1SDE2[index, 1:7]   <- df$x1
            dataOutNONHOIsingle1SDE2[index, 8:14]  <- df$x2
            dataOutNONHOIsingle1SDE2[index, 15:21] <- df$x3
            dataOutNONHOIsingle1SDE2[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start1, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIsingle1SDE3[index, 1:7]   <- df$x1
            dataOutNONHOIsingle1SDE3[index, 8:14]  <- df$x2
            dataOutNONHOIsingle1SDE3[index, 15:21] <- df$x3
            dataOutNONHOIsingle1SDE3[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start2, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIsingle2SDE1[index, 1:7]   <- df$x1
            dataOutNONHOIsingle2SDE1[index, 8:14]  <- df$x2
            dataOutNONHOIsingle2SDE1[index, 15:21] <- df$x3
            dataOutNONHOIsingle2SDE1[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start2, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIsingle2SDE2[index, 1:7]   <- df$x1
            dataOutNONHOIsingle2SDE2[index, 8:14]  <- df$x2
            dataOutNONHOIsingle2SDE2[index, 15:21] <- df$x3
            dataOutNONHOIsingle2SDE2[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start2, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIsingle2SDE3[index, 1:7]   <- df$x1
            dataOutNONHOIsingle2SDE3[index, 8:14]  <- df$x2
            dataOutNONHOIsingle2SDE3[index, 15:21] <- df$x3
            dataOutNONHOIsingle2SDE3[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start3, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIsingle3SDE1[index, 1:7]   <- df$x1
            dataOutNONHOIsingle3SDE1[index, 8:14]  <- df$x2
            dataOutNONHOIsingle3SDE1[index, 15:21] <- df$x3
            dataOutNONHOIsingle3SDE1[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start3, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIsingle3SDE2[index, 1:7]   <- df$x1
            dataOutNONHOIsingle3SDE2[index, 8:14]  <- df$x2
            dataOutNONHOIsingle3SDE2[index, 15:21] <- df$x3
            dataOutNONHOIsingle3SDE2[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start3, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIsingle3SDE3[index, 1:7]   <- df$x1
            dataOutNONHOIsingle3SDE3[index, 8:14]  <- df$x2
            dataOutNONHOIsingle3SDE3[index, 15:21] <- df$x3
            dataOutNONHOIsingle3SDE3[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start4, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIpairwise1SDE1[index, 1:7]   <- df$x1
            dataOutNONHOIpairwise1SDE1[index, 8:14]  <- df$x2
            dataOutNONHOIpairwise1SDE1[index, 15:21] <- df$x3
            dataOutNONHOIpairwise1SDE1[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start4, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIpairwise1SDE2[index, 1:7]   <- df$x1
            dataOutNONHOIpairwise1SDE2[index, 8:14]  <- df$x2
            dataOutNONHOIpairwise1SDE2[index, 15:21] <- df$x3
            dataOutNONHOIpairwise1SDE2[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start4, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIpairwise1SDE3[index, 1:7]   <- df$x1
            dataOutNONHOIpairwise1SDE3[index, 8:14]  <- df$x2
            dataOutNONHOIpairwise1SDE3[index, 15:21] <- df$x3
            dataOutNONHOIpairwise1SDE3[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start5, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIpairwise2SDE1[index, 1:7]   <- df$x1
            dataOutNONHOIpairwise2SDE1[index, 8:14]  <- df$x2
            dataOutNONHOIpairwise2SDE1[index, 15:21] <- df$x3
            dataOutNONHOIpairwise2SDE1[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start5, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIpairwise2SDE2[index, 1:7]   <- df$x1
            dataOutNONHOIpairwise2SDE2[index, 8:14]  <- df$x2
            dataOutNONHOIpairwise2SDE2[index, 15:21] <- df$x3
            dataOutNONHOIpairwise2SDE2[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start5, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIpairwise2SDE3[index, 1:7]   <- df$x1
            dataOutNONHOIpairwise2SDE3[index, 8:14]  <- df$x2
            dataOutNONHOIpairwise2SDE3[index, 15:21] <- df$x3
            dataOutNONHOIpairwise2SDE3[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start6, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIpairwise3SDE1[index, 1:7]   <- df$x1
            dataOutNONHOIpairwise3SDE1[index, 8:14]  <- df$x2
            dataOutNONHOIpairwise3SDE1[index, 15:21] <- df$x3
            dataOutNONHOIpairwise3SDE1[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start6, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIpairwise3SDE2[index, 1:7]   <- df$x1
            dataOutNONHOIpairwise3SDE2[index, 8:14]  <- df$x2
            dataOutNONHOIpairwise3SDE2[index, 15:21] <- df$x3
            dataOutNONHOIpairwise3SDE2[index, 22]    <- "nonhoi"
            
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start6, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            colnames(df) <- c("time", "x1", "x2", "x3")
            dataOutNONHOIpairwise3SDE3[index, 1:7]   <- df$x1
            dataOutNONHOIpairwise3SDE3[index, 8:14]  <- df$x2
            dataOutNONHOIpairwise3SDE3[index, 15:21] <- df$x3
            dataOutNONHOIpairwise3SDE3[index, 22]    <- "nonhoi"
            
            index <- index + 1
          }
        },
        error=function(cond) {
          #index <- index + 1
          #print("error?")
          #message("ERROR: Possibly NaN")
        }
      )
      
      
    }
  }
  return(list(df_HOI = dataOutHOI, 
              df_HOI1SDE1 = dataOutHOIsingle1SDE1,
              df_HOI1SDE2 = dataOutHOIsingle1SDE2,
              df_HOI1SDE3 = dataOutHOIsingle1SDE3,
              df_HOI2SDE1 = dataOutHOIsingle2SDE1,
              df_HOI2SDE2 = dataOutHOIsingle2SDE2,
              df_HOI2SDE3 = dataOutHOIsingle2SDE3,
              df_HOI3SDE1 = dataOutHOIsingle3SDE1,
              df_HOI3SDE2 = dataOutHOIsingle3SDE2,
              df_HOI3SDE3 = dataOutHOIsingle3SDE3,
              df_HOI12SDE1 = dataOutHOIpairwise1SDE1,
              df_HOI12SDE2 = dataOutHOIpairwise1SDE2,
              df_HOI12SDE3 = dataOutHOIpairwise1SDE3,
              df_HOI13SDE1 = dataOutHOIpairwise2SDE1,
              df_HOI13SDE2 = dataOutHOIpairwise2SDE2,
              df_HOI13SDE3 = dataOutHOIpairwise2SDE3,
              df_HOI23SDE1 = dataOutHOIpairwise3SDE1,
              df_HOI23SDE2 = dataOutHOIpairwise3SDE2,
              df_HOI23SDE3 = dataOutHOIpairwise3SDE3,
              df_HOI123SDE1 = dataOutHOISDE1,
              df_HOI123SDE2 = dataOutHOISDE2,
              df_HOI123SDE3 = dataOutHOISDE3,
              df_NONHOI = dataOutNONHOI, 
              df_NONHOI1SDE1 = dataOutNONHOIsingle1SDE1,
              df_NONHOI1SDE2 = dataOutNONHOIsingle1SDE2,
              df_NONHOI1SDE3 = dataOutNONHOIsingle1SDE3,
              df_NONHOI2SDE1 = dataOutNONHOIsingle2SDE1,
              df_NONHOI2SDE2 = dataOutNONHOIsingle2SDE2,
              df_NONHOI2SDE3 = dataOutNONHOIsingle2SDE3,
              df_NONHOI3SDE1 = dataOutNONHOIsingle3SDE1,
              df_NONHOI3SDE2 = dataOutNONHOIsingle3SDE2,
              df_NONHOI3SDE3 = dataOutNONHOIsingle3SDE3,
              df_NONHOI12SDE1 = dataOutNONHOIpairwise1SDE1,
              df_NONHOI12SDE2 = dataOutNONHOIpairwise1SDE2,
              df_NONHOI12SDE3 = dataOutNONHOIpairwise1SDE3,
              df_NONHOI13SDE1 = dataOutNONHOIpairwise2SDE1,
              df_NONHOI13SDE2 = dataOutNONHOIpairwise2SDE2,
              df_NONHOI13SDE3 = dataOutNONHOIpairwise2SDE3,
              df_NONHOI23SDE1 = dataOutNONHOIpairwise3SDE1,
              df_NONHOI23SDE2 = dataOutNONHOIpairwise3SDE2,
              df_NONHOI23SDE3 = dataOutNONHOIpairwise3SDE3,
              df_NONHOI123SDE1 = dataOutNONHOISDE1,
              df_NONHOI123SDE2 = dataOutNONHOISDE2,
              df_NONHOI123SDE3 = dataOutNONHOISDE3
  ))
}




















gen_samples_training_V2 <- function(number, noise, parSingle, parPair, parThree) {
  # Function to generate synthetic samples to train ML models.
  # INPUT: number of samples, noise level, 1-species parameters, 
  # 2-species parameters, 3-species paramenters.
  # OUTPUT: Dataframes with HOI and NON-HOI samples.
  # TODO: Remove hard coded dataframes with noisy samples.
  
  modelHOI= function(u, p, t) {
    r1 = p[[1]]
    r2 = p[[2]]
    r3 = p[[3]]
    
    K1 = p[[4]]
    K2 = p[[5]]
    K3 = p[[6]]
    
    a12 = p[[7]]
    a13 = p[[8]]
    a21 = p[[9]]
    a23 = p[[10]]
    a31 = p[[11]]
    a32 = p[[12]]
    
    b123 = p[[13]]
    b231 = p[[14]]
    b312 = p[[15]]
    
    x1 = u[1]
    x2 = u[2]
    x3 = u[3]
    
    dx1 = r1*x1*(1 - x1/K1 + a12*x2/K2 + a13*x3/K3 + b123*x2*x3/(K2*K3))
    dx2 = r2*x2*(1 - x2/K2 + a23*x3/K3 + a21*x1/K1 + b231*x3*x1/(K3*K1))
    dx3 = r3*x3*(1 - x3/K3 + a31*x1/K1 + a32*x2/K2 + b312*x1*x2/(K1*K2))
    return(c(dx1, dx2, dx3))
  }
  
  modelNONHOI= function(u, p, t) {
    r1 = p[[1]]
    r2 = p[[2]]
    r3 = p[[3]]
    
    K1 = p[[4]]
    K2 = p[[5]]
    K3 = p[[6]]
    
    a12 = p[[7]]
    a13 = p[[8]]
    a21 = p[[9]]
    a23 = p[[10]]
    a31 = p[[11]]
    a32 = p[[12]]
    
    x1 = u[1]
    x2 = u[2]
    x3 = u[3]
    
    dx1 = r1*x1*(1 - x1/K1 + a12*x2/K2 + a13*x3/K3)
    dx2 = r2*x2*(1 - x2/K2 + a23*x3/K3 + a21*x1/K1)
    dx3 = r3*x3*(1 - x3/K3 + a31*x1/K1 + a32*x2/K2)
    
    return(c(dx1, dx2, dx3))
  }

  stochastic <- function(u, p, t) {
    #return(0.2^2*c(sqrt(u[1]), sqrt(u[2]), sqrt(u[3])))
    return(-noise^2 * c(u[1], u[2], u[3]))
  }
  
  # Empty dataframes.
  dataOutHOI <- data.frame(0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0,
                           "hoi")
  colnames(dataOutHOI) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                            "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                            "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                            "isHOI")
  
  dataOutNONHOI <- data.frame(0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0,
                           "nonhoi")
  colnames(dataOutNONHOI) <- c("x11", "x12", "x13", "x14", "x15", "x16", "x17",
                            "x21", "x22", "x23", "x24", "x25", "x26", "x27",
                            "x31", "x32", "x33", "x34", "x35", "x36", "x37",
                            "isHOI")
  
  # Index for keeping track of number of samples generated.
  index <- 1
  # Tolerance level to accept "sufficiently different" HOI and NON-HOI samples.
  TOL <- 1
  # Muliplier factor to modify bijk limits.
  jMul <- 1
  
  while (index < number + 1) {
    # Print progress.
    if (index == 1) { cat("Progress:", index, "\n") }
    if (index %% 10 == 0) { cat("Progress:", index, "\n") }
    
    # Generate random initial conditions (hard coded).
    start = c(x1=runif(1, min=2, max=10), 
              x2=runif(1, min=2, max=10), 
              x3=runif(1, min=2, max=10))
    
    # Prepare initial conditions for 1- and 2-species.
    start1 <- c(start[1], 0, 0)
    start2 <- c(0, start[2], 0)
    start3 <- c(0, 0, start[3])
    start4 <- c(start[1], start[2], 0)
    start5 <- c(start[1], 0, start[3])
    start6 <- c(0, start[2], start[3])
    
    # Generate random model parameters (hard coded).
    r1 = parSingle[1,1] * (1 + runif(1, -0.2, 0.2))
    r2 = parSingle[1,2] * (1 + runif(1, -0.2, 0.2))
    r3 = parSingle[1,3] * (1 + runif(1, -0.2, 0.2))
    
    K1 = parSingle[2,1] * (1 + runif(1, -0.2, 0.2))
    K2 = parSingle[2,2] * (1 + runif(1, -0.2, 0.2))
    K3 = parSingle[2,3] * (1 + runif(1, -0.2, 0.2))
    
    a12 = parPair[1,1] * (1 + runif(1, -0.2, 0.2))
    a21 = parPair[2,1] * (1 + runif(1, -0.2, 0.2))
    
    a23 = parPair[1,3] * (1 + runif(1, -0.2, 0.2))
    a32 = parPair[2,3] * (1 + runif(1, -0.2, 0.2))
    
    a13 = parPair[1,2] * (1 + runif(1, -0.2, 0.2))
    a31 = parPair[2,2] * (1 + runif(1, -0.2, 0.2))  
    
    j1 <- parThree[1] * runif(1, -2.0 * jMul, 2.0 * jMul)
    j2 <- parThree[2] * runif(1, -2.0 * jMul, 2.0 * jMul)
    j3 <- parThree[3] * runif(1, -2.0 * jMul, 2.0 * jMul)
    
    # Check if HOI parameters are "big enough" (hard coded).
    flag1 = abs(j1) > 0.00001
    flag2 = abs(j2) > 0.00001
    flag3 = abs(j3) > 0.00001
    
    if (flag1 | flag2 | flag3) {
      
      b123 <- j1
      b231 <- j2
      b312 <- j3
      
      # Parameter vector par.
      par = c(r1=r1, r2=r2, r3=r3, 
              K1=K1, K2=K2, K3=K3, 
              a12=a12, a13=a13, a21=a21, a23=a23, a31=a31, a32=a32,
              b123=b123, b231=b231, b312=b312)
      
      tryCatch(
        {
          # HOI sample.
          outHOI = diffeqr::ode.solve(f = modelHOI, u0 = start, tspan = list(0, 7), p = par, saveat=1.0)
          out2HOI = data.frame(time = outHOI$t, x1 = outHOI$u[,1], x2 = outHOI$u[,2], x3 = outHOI$u[,3])
          out2HOI = out2HOI[-1,]
          colnames(out2HOI) <- c("time", "x1", "x2", "x3")
          # NON-HOI sample.
          outNONHOI = diffeqr::ode.solve(f = modelNONHOI, u0 = start, tspan = list(0, 7), p = par, saveat=1.0)
          out2NONHOI = data.frame(time = outNONHOI$t, x1 = outNONHOI$u[,1], x2 = outNONHOI$u[,2], x3 = outNONHOI$u[,3])
          out2NONHOI = out2NONHOI[-1,]
          colnames(out2NONHOI) <- c("time", "x1", "x2", "x3")
          
          # Check tolerance level condition.
          norm_HOINONHOI <- dist(rbind(unlist(out2HOI), unlist(out2NONHOI)))
          #print(norm_HOINONHOI)
          flag4 <- (norm_HOINONHOI > TOL)
          #print(flag4)
          
          if (flag4) {
            # Save HOI sample.
            sol = diffeqr::sde.solve(f = modelHOI, g = stochastic, u0 = start, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            dataOutHOI[index, 1:7]   <- df$x1
            dataOutHOI[index, 8:14]  <- df$x2
            dataOutHOI[index, 15:21] <- df$x3
            dataOutHOI[index, 22]    <- "hoi"
            
            # Save NON-HOI sample.
            sol = diffeqr::sde.solve(f = modelNONHOI, g = stochastic, u0 = start, tspan = list(0, 7), p = par, saveat=1.0)
            df = data.frame(time = sol$t, x1 = sol$u[,1], x2 = sol$u[,2], x3 = sol$u[,3])
            df = df[-1,]
            dataOutNONHOI[index, 1:7]   <- df$x1
            dataOutNONHOI[index, 8:14]  <- df$x2
            dataOutNONHOI[index, 15:21] <- df$x3
            dataOutNONHOI[index, 22]    <- "nonhoi"
            
            index = index + 1

            jMul <- 1
          }
     else { jMul <- jMul + 0.1 # increase jMul to increase noise level por bijk
          }
        },
        error=function(cond) {
          #index <- index + 1
          #print("error?")
          #message("ERROR: Possibly NaN")
        }
      )
      
      
    }
  }
  return(list(df_HOI = dataOutHOI, 
              df_NONHOI = dataOutNONHOI
              )
         )
}
