#########################################################
## Stat 202A - Homework 6
## Author: 
## Date : 
## Description: This script implements QR decomposition
## and linear regression based on QR
#########################################################

#############################################################
## INSTRUCTIONS: Please fill in the missing lines of code
## only where specified. Do not change function names, 
## function inputs or outputs. You can add examples at the
## end of the script (in the "Optional examples" section) to 
## double-check your work, but MAKE SURE TO COMMENT OUT ALL 
## OF YOUR EXAMPLES BEFORE SUBMITTING.
##
## Very important: Do not use the function "setwd" anywhere
## in your code. If you do, I will be unable to grade your 
## work since R will attempt to change my working directory
## to one that does not exist.
#############################################################

##################################
## Function 1: QR decomposition ##
##################################

myQR <- function(A){
  
  ## Perform QR factorization on the matrix A
  ## FILL IN CODE HERE ##
  
  n = dim(A)[1];
  m = dim(A)[2];
  
  R = A
  Q = diag(n);
  
  for (k in 1:m)
  {
    x = array(0,c(n,1));
    x[k:n,1] = R[k:n,k];
    v = x;
    v[k] = x[k] + sign(x[k,1]) * norm(x);
    
    s = norm(v);
    if(s != 0){
      u = v/s;
      R = R - 2*(u%*%(t(u)%*%R));
      Q = Q - 2*(u%*%(t(u)%*%Q));
    }
  }
  
  ## Function should output a list with Q.transpose and R
  return(list("Q" = t(Q), "R" = R))
  
}

###############################################
## Function 2: Linear regression based on QR ##
###############################################

myLM <- function(X, Y){
  
  ## X is an n x p matrix of explanatory variables
  ## Y is an n dimensional vector of responses
  ## Do NOT simulate data in this function. n and p
  ## should be determined by X.
  ## Use myQR inside of this function
  
  n = dim(X)[1];
  p = dim(X)[2];
    
  ## FILL CODE HERE ##
  ## Assuming Y is a n X 1 vector
  Z = cbind(matrix(rep(1, n), nrow = n), X, Y);
  list_qr = myQR(Z);
  
  R = list_qr[[2]]
  
  R1 = R[1:(p+1), 1:(p+1)];
  Y1 = R[1:(p+1), p+2];
  
  beta_ls = solve(R1, Y1);
  ## Function returns beta_ls, the least squares
  ## solution vector
  return(beta_ls)
}

#n= 10
#p = 5
#X = matrix(rnorm(n*p), nrow=n)
#beta_true = matrix(rep(0, p), nrow = p)
#Y = X %*% beta_true + rnorm(n)
#print(myLM(X,Y))

