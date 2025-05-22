# CausalSemiComp::SimDataWeibFrail
# shape = k (alpha), scale = 1/lambda
SimDataWeibFrail = function(n.sample, params, no.protected = T, no.large = T,
                            cens.exp.rate = 0.1, cens, cens.admin = 100, 
                            round.times = T, X = NULL, dim_x, beta_A, tau){
  if (!is.null(X) & (no.large | no.protected)) 
    stop("Can't specify X and ask no large or no protected")
  list2env(params, envir = environment())
  
  gamma.scale <- params$theta
  gamma.common.shape <- params$rho / params$theta
  gamma.each.shape <- (1 - params$rho) / params$theta
  
  n.sample.temp <- n.sample
  cond.sample <- F # this argument determines the while loop
  while (cond.sample == F) {
    if (no.protected | no.large) {
      n.sample.temp <- n.sample.temp * 4
    }
    if (is.null(X)) {
      x1 <- rbinom(n = n.sample.temp, size = 1, prob = 0.5)
      x2 <- rnorm(n.sample.temp)
      if(dim_x==2){X <- cbind(x1, x2)}
      if(dim_x==1){X <- x2}
    }else {
      n.times <- ceiling(n.sample.temp/nrow(X))
      X <- do.call(rbind, replicate(n.times, X, simplify = F)) %>% 
        as.matrix
      n.sample.temp <- nrow(X)
    }
    # frailty generation
    # gamma.scale = 1/beta in original gamma dist
    gamma.common <- rgamma(n.sample.temp, shape = gamma.common.shape, 
                                    scale = gamma.scale)
    gamma0 <- gamma.common + rgamma(n.sample.temp, shape = gamma.each.shape, 
                                    scale = gamma.scale)
    gamma1 <- gamma.common + rgamma(n.sample.temp, shape = gamma.each.shape, 
                                    scale = gamma.scale)
    gamma.out <- cbind(gamma0, gamma1)

    U1.0 <- runif(n.sample.temp)
    U1.1 <- runif(n.sample.temp)
    U12.0 <- runif(n.sample.temp)
    U2.0 <- runif(n.sample.temp)
    U2.1 <- runif(n.sample.temp)
    U12.1 <- runif(n.sample.temp)
    
    # X %*% params$beta when there 2 and more covariaes. * when there is one covariate
    # DN: expb001.gamma <- gamma0 * exp(X %*% params$beta.a0.01) 
    #TODO  make this part more general so it could fit any dim_x 
    # matrix(X) %*% params$beta.a0.01[2] (with dim_x=1) 
    # and matrix(X) %*% params$beta.a0.01 (dim_x=2) should work 
    
    expb001.gamma <- gamma0 * exp(as.matrix(X) %*% params$beta.a0.01[2]) # [2]
    expb002.gamma <- gamma0 * exp(as.matrix(X) %*% params$beta.a0.02[2])
    expb012.gamma <- gamma0 * exp(as.matrix(X) %*% params$beta.a0.12[2])
    expb101.gamma <- gamma1 * exp(as.matrix(X) %*% params$beta.a1.01[2])
    expb102.gamma <- gamma1 * exp(as.matrix(X) %*% params$beta.a1.02[2])
    expb112.gamma <- gamma1 * exp(as.matrix(X) %*% params$beta.a1.12[2])
    
    # generating 4 PO times: [ -log(U) / (base_haz * lambda^alpha)] ^ (1/alpha)
    # note that U ~ 1-U
    T1.0 <- (-log(U1.0)/(expb001.gamma * 
       params$base.weib.scale.a0.01^(-params$base.weib.shape.a0.01)))^
      (1/params$base.weib.shape.a0.01)
    T1.1 <- (-log(U1.1)/(expb101.gamma * params$base.weib.scale.a1.01^
      (-params$base.weib.shape.a1.01)))^
      (1/params$base.weib.shape.a1.01)
    T2.0 <- (-log(U2.0)/(expb002.gamma * params$base.weib.scale.a0.02^
      (-params$base.weib.shape.a0.02)))^
      (1/params$base.weib.shape.a0.02)
    T2.1 <- (-log(U2.1)/(expb102.gamma * params$base.weib.scale.a1.02^
       (-params$base.weib.shape.a1.02)))^
      (1/params$base.weib.shape.a1.02)
    if (round.times == T) {
      T1.0 <- round(T1.0, 1)
      T1.1 <- round(T1.1, 1)
      T2.0 <- round(T2.0, 1)
      T2.1 <- round(T2.1, 1)
    }
    # if infected first, change T2(0) according expb012.gamma, and T2(1) according expb112.gamma,
    # due to the transformation 1 -> 2
    for (i in 1:n.sample.temp) {
      if (T2.0[i] >= T1.0[i]) {
        U12.0i <- U12.0[i]
        T2.0[i] <- (-log(U12.0i)/(expb012.gamma[i] * 
            params$base.weib.scale.a0.12^(-params$base.weib.shape.a0.12)) + 
            T1.0[i]^params$base.weib.shape.a0.12)^(1/params$base.weib.shape.a0.12)
        if (round.times == T) {
          T2.0[i] <- round(T2.0[i], 1)
          if (T2.0[i] == T1.0[i]) {
            T2.0[i] <- round(T1.0[i] + runif(1, 0.1, 1), 1)
          }
        }
      }
      if (T2.1[i] >= T1.1[i]) {
        U12.1i <- U12.1[i]
        T2.1[i] <- (-log(U12.1i)/(expb112.gamma[i] * 
                                    params$base.weib.scale.a1.12^(-params$base.weib.shape.a1.12)) + 
                      T1.1[i]^params$base.weib.shape.a1.12)^(1/params$base.weib.shape.a1.12)
        if (round.times == T) {
          T2.1[i] <- round(T2.1[i], 1)
          if (T2.1[i] == T1.1[i]) {
            T2.1[i] <- round(T1.1[i] + runif(1, 0.1, 
                                             1), 1)
          }
        }
      }
    }
    # something with protected (those not infected under treatment, but does under ctr) and large
    # out sums no.protected and no.large, 
    # by removing protected and those with high T2.0 and high T2.1, if needed
    out.protected <- (T1.0 < T2.0) & (T1.1 > T2.1)
    out.large <- T2.0 > 50 | T2.1 > 50
    if (no.protected == T & no.large == T) {
      out <- out.protected | out.large
    }
    if (no.protected == T & no.large == F) {
      out <- out.protected
    }
    if (no.protected == F & no.large == T) {
      out <- out.large
    }
    if (any(no.protected, no.large)) {
      T1.0 <- T1.0[!out]
      T1.1 <- T1.1[!out]
      T2.0 <- T2.0[!out]
      T2.1 <- T2.1[!out]
      X <- X[!out, ]
      gamma.out <- gamma.out[!out, ]
    }
    n.sample.real <- length(T1.0)
    if (n.sample.real >= n.sample) {
      cond.sample <- T
    }
  } # end of while loop 
  # determined by cond.sample, ends where cond.sample <- T, i.e. n.sample.real >= n.sample
  
  # truncate time vectors, X, and gamma in n.sample
  T1.0 <- T1.0[1:n.sample]
  T1.1 <- T1.1[1:n.sample]
  T2.0 <- T2.0[1:n.sample]
  T2.1 <- T2.1[1:n.sample]
  
  if (round.times == T) {
    T2.0[T1.0 == 0] <- pmax(T2.0[T1.0 == 0], 0.1)
    T2.1[T1.1 == 0] <- pmax(T2.1[T1.1 == 0], 0.1)
    T1.0[T1.0 == 0] <- 0.1
    T1.1[T1.1 == 0] <- 0.1
    T2.0[T2.0 == 0] <- 0.1
    T2.1[T2.1 == 0] <- 0.1
  }
  X <- as.matrix(X)[1:n.sample, ] # X[1:n.sample, ]
  gamma.out <- gamma.out[1:n.sample, ]
  
  if(is.null(params$ceffs_A) == TRUE){
    # randomized treatment assignment
    A <- rbinom(n = n.sample, size = 1, prob = 0.5)
  }else{
    #20.01.25 non-randomized treatment assignment
    A_probs = expit(as.matrix(cbind(rep(1, n.sample), X)) %*% params$ceffs_A)
    A <- rbinom(n = n.sample, size = 1, prob = A_probs)
  }
  
  # censure 
  #TODO TZ:22.01.25 censure is used for observed times, not POs.
  C <- rexp(n.sample, rate = cens.exp.rate)
  if (round.times == T) {
    C <- round(C, 1)
    C[C == T1.0 | C == T2.0 | C == T1.1 | C == T2.1] <- C[C == 
              T1.0 | C == T2.0 | C == T1.1 | C == T2.1] + 0.05
  }
  
  # update times (tilde in the paper) and delta1, delta2 
  if(cens == T){
    T1 <- T2 <- delta1 <- delta2 <- vector(length = n.sample)
    T1[A == 0] <- pmin(T1.0[A == 0], T2.0[A == 0], C[A == 0], cens.admin)
    T1[A == 1] <- pmin(T1.1[A == 1], T2.1[A == 1], C[A == 1], cens.admin)
    T2[A == 0] <- pmin(T2.0[A == 0], C[A == 0], cens.admin)
    T2[A == 1] <- pmin(T2.1[A == 1], C[A == 1], cens.admin)
  }
  if(cens == F){
    T1 <- T2 <- delta1 <- delta2 <- vector(length = n.sample)
    T1[A == 0] <- pmin(T1.0[A == 0], T2.0[A == 0], cens.admin) #, cens.admin
    T1[A == 1] <- pmin(T1.1[A == 1], T2.1[A == 1], cens.admin) #, cens.admin
    T2[A == 0] <- pmin(T2.0[A == 0], cens.admin) #, cens.admin
    T2[A == 1] <- pmin(T2.1[A == 1], cens.admin) #, cens.admin
  }
  
  # delta1[A == 0] <- T1[A == 0] == T1.0[A == 0]
  # delta1[A == 1] <- T1[A == 1] == T1.1[A == 1]
  # delta2[A == 0] <- T2[A == 0] == T2.0[A == 0]
  # delta2[A == 1] <- T2[A == 1] == T2.1[A == 1]
  delta1.0 <- ifelse(T1.0 <= tau & T1.0 <= T2.0, 1, 0)
  delta1.1 <- ifelse(T1.1 <= tau & T1.1 <= T2.1, 1, 0)
  delta2.0 <- ifelse(T2.0 <= tau, 1, 0)
  delta2.1 <- ifelse(T2.1 <= tau, 1, 0)
  
  # change T2 for patients with T1 occurence (infection), but only after cens.admin
  T2[delta1 == 1 & T1 == cens.admin] <- cens.admin + 0.05
  list.to.return <- list(T1.0 = T1.0, T1.1 = T1.1, T2.0 = T2.0, T2.1 = T2.1, X = X 
     #,T1 = T1, T2 = T2, delta1 = delta1, delta2 = delta2
     , A = A, C = C
     , delta1.0=delta1.0, delta1.1=delta1.1, delta2.0=delta2.0, delta2.1=delta2.1
     , gamma.out = gamma.out)
}
