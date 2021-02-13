estimate <- function(dataOut) {
  
  splitted.data <- multisplit(dataOut, c("replicate"))
  
  lower1 <- c(y0 = 5,  mumax = -5.0)
  p1     <- c(y0 = 10,  mumax = 0.5)
  upper1 <- c(y0 = 15, mumax = 5.0)
   
  r_A = c()
  pgr_A = c()
   
  for (ii in as.set(dataOut$replicate)) {
   dat <- splitted.data[[ii]]
   #print(dat)
     
    # try(
    #   fit1 <- fit_growthmodel(FUN = grow_exponential,
    #                           p = p1,
    #                           dat$week,
    #                           dat$x1,
    #                           lower = lower1,
    #                           upper = upper1)
    # )
    # 
    # r_A[ii] = coef(fit1)[2]
    
    # 3.5.3 from "A Primer in Ecology with R", M Henry H Stevens
    n <- nrow(dat)
    N.change <- dat$x1[-1]/dat$x1[-n]
    #cat("Nchange:", N.change, "\n")
    interval <- diff(dat$week)
    #cat("interval:", interval, "\n")
    pgr <- log(N.change)/interval
    #cat("pgr:", pgr, "\n")
    #pgr <- na.omit(pgr)
    Nt <- dat$x1[-n]
    #plot(pgr ~ Nt)
    tryCatch(
      {
        mod1 <- lm(pgr ~ Nt)
        r <- coef(mod1)[[1]]
        #pgr_A[ii] <- r
        pgr_A[ii] <- mean(pgr)
        K <- -1/coef(mod1)[[2]]
      }
    )
    #summary(mod1)
  }
  return(list("r"=r_A, 
              "pgr"=pgr_A)
              #"pgrTOTAL"=pgr)
         )
}

bendertest <- function(single, pairwise1, pairwise2, three) {

# single
list <- estimate(single)
r_A <- c()
r_A = rbind(r_A, list$r)
pgr_A = c()
pgr_A = rbind(pgr_A, list$pgr)

# pairwise
list <- estimate(pairwise1)
r_AB <- c()
r_AB = rbind(r_AB, list$r)
pgr_AB = c()
pgr_AB = rbind(pgr_AB, list$pgr)
list <- estimate(pairwise2)
r_AC <- c()
r_AC = rbind(r_AC, list$r)
pgr_AC = c()
pgr_AC = rbind(pgr_AC, list$pgr)

# three
list <- estimate(three)
r_ABC <- c()
r_ABC = rbind(r_ABC, list$r)
pgr_ABC = c()
pgr_ABC = rbind(pgr_ABC, list$pgr)

HOIcount <- 0

# Test <- t.test(r_A + r_ABC, r_AB + r_AC, var.equal = FALSE)
# if (Test$p.value < 0.05) {
#   #print("HOI")
#   HOIcount <- HOIcount + 1
# } else {
#     #print("NON-HOI")
#   }

# #t.test(pgr_LHS, pgr_RHS)
#print(mean(pgr_A + pgr_ABC))
#print(mean(pgr_AB + pgr_AC))
Test <- t.test(pgr_A + pgr_ABC, pgr_AB + pgr_AC, var.equal = FALSE, alternative = "two.sided", paired = TRUE)
if (Test$p.value < 0.01) {
  #print("HOI")
  HOIcount <- HOIcount + 1
} else {
  #print("NON-HOI")
}

return(HOIcount)
}
