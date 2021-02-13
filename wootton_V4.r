woottontest <- function(single, pairwise1, pairwise2, three) {
  
  get_LHS <- function(ii, jj, tt) {
    return(three$x1[7*(ii-1)+tt] * single$x1[7*(jj-1)+tt])
  }
  
  get_RHS <- function(ii, jj, tt) {
    return(pairwise1$x1[7*(ii-1)+tt] * pairwise2$x1[7*(jj-1)+tt])
  }

  LHS = c()
  RHS = c()
  
  t_fix <- 2

  for (ii in 1:3){
      LHS[ii] <- get_LHS(ii, ii, t_fix)
      RHS[ii] <- get_RHS(ii, ii, t_fix)
  }
  
  #print(LHS)
  #print(RHS)
  
  HOIcount <- 0
  
  Test <- t.test(LHS, RHS, var.equal = FALSE)

  if (!is.nan(Test$p.value) & Test$p.value < 0.05) { 
    #print("HOI") 
    HOIcount <- HOIcount + 1
  } else { 
    #print("NON-HOI") 
  }

  return(HOIcount)
}