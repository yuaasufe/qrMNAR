qrMNAR

R package "qrMNAR" for estimation of quantile regression coefficients with nonignorable missing data. Provide an algorithm  to implement the proposed method to estimate quantile regression coefficients with nonignorable missing data. It also provides a bootstrap method to estimate the standard error of the  proposed estimation.



Installation

Requirements

- Rtools
  

Package Install

    # install.packages("quantreg")
    # install.packages("Rcpp")
    # install.packages("foreach")
    # install.packages("snow")
    # install.packages("doSNOW")
    
    ## download "qrMNAR_0.1.0.tar.gz", then use the following code to install this R package
    install.packages("path/to/qrMNAR_0.1.0.tar.gz", repos = NULL, type = "source")



Example

    library(qrMNAR)
    
    #### This package has two main functions "msIpwQr" and "estimate_sd.boot". "msIpwQr" is used to estimate quantile regression coefficients; "estimate_sd.boot" is used to estimate standard errors of the proposed estimates via weighted bootstrapping and function "estimate_sd.boot" can using parallel computing for bootstrapping. More detailed information about these two functions can be obtained through
    
    help("msIpwQr")
    help("estimate_sd.boot")
    
    #### Example 1: Call function "msIpwQr" to estimate quantile regression coefficients
    ## generate sim_data for example
    nsize <- 1000
    beta <- c(1,-2,2,0.5) ; gam <- c(0.5, 0.5, 0.2, 0); theta <- c(-2, 0.5,0.5,0.5)
    x1 <- rbinom(nsize, 1, 0.5); x2 <- rnorm(nsize, mean = 2, sd = 0.5);
    x3 <- runif(nsize, min = 0, max = 2); xx <- cbind(1, x1, x2, x3);
    err <- xx%*%gam*rnorm( nsize, mean = 0, sd = 0.5); yy <- xx %*% beta + err
    xy <- cbind(1, x1, x2, yy); prob.m <- 1/(1+ exp(- xy %*% theta))
    delta.ind <- sapply(prob.m, function(p){rbinom(n=1, size=1, prob=p)})
    yy[delta.ind==0] <- NA
    set.seed(100)
    sim_data <- data.frame(x1, x2, x3, yy)
    p1<- 1; p2<- 1; q1<- 0; q2 <- 1
    
    ## set the initial value related to the algorithm
    thresh = 1e-04; max_iter.ms <- 50;  m <- 10
    tau_seq.ms  <- seq(1/50, 49/50, by = 1/50); tau_seq.out <- seq(0.25, 0.75, by = 0.25)
    
    ## using sim_data to estimate quantile regression coefficients
    out <- msIpwQr(sim_data, p1, p2, q1, q2, thresh, max_iter.ms, tau_seq.out, tau_seq.ms, m )
    # true value of quantile regression coefficients
    beta_tau_true <- beta + gam %*% t(qnorm( tau_seq.out, mean = 0, sd = 0.5)) 
    # quantile regression coefficients estimates
    out[[3]] 
    # summary of quantile regression model
    out[[7]] 
    
    
    ####  Example 2: Call function "estimate_sd.boot" to get a bootstrap estimate of the standard errors of resulting estimators; this function implement bootstrap through parallel computing
    ## set the initial value related to the algorithm
    library(parallel)
    num_cores <- detectCores()
    B = 100
     
    ## to obtain standard errors of quantile regression coefficients estimates via bootstrapping by parallel computing
    out_sd <- estimate_sd.boot(sim_data, p1, p2, q1, q2, thresh, max_iter.ms,                                                   tau_seq.out,  tau_seq.ms, m, B, num_cores)
    # standard errors of quantile regression coefficients estimates via bootstrapping
    out_sd[[3]] 


