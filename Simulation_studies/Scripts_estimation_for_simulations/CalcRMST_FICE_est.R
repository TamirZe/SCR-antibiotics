# CausalSemiComp::CalcRMST

CalcRMST_one_formula = function (rho, tau, n.gamma.vals, n.sample.pers
                     #, population, 
                      ,Xnames, Xnames_formula, data, res, list.out = T, detailed = F){
  #n.X.vals <- nrow(population)
  n.X.vals = nrow(data)
  m.A1 <- mean(data$A == 1)
  
  #Xmat <- as.matrix(subset(data, select = Xnames))
  formula = as.formula(paste0("Surv(T1, delta1) ~", paste0(Xnames_formula, collapse = "+")))
  fit <- survival::coxph(formula, data = data)
  Xmat = model.matrix(fit)
  
  exp.b001 <- exp(Xmat %*% coef(res$fit.list$fit.a0.01))
  exp.b002 <- exp(Xmat %*% coef(res$fit.list$fit.a0.02))
  exp.b012 <- exp(Xmat %*% coef(res$fit.list$fit.a0.12))
  exp.b101 <- exp(Xmat %*% coef(res$fit.list$fit.a1.01))
  exp.b102 <- exp(Xmat %*% coef(res$fit.list$fit.a1.02))
  exp.b112 <- exp(Xmat %*% coef(res$fit.list$fit.a1.12))
  
  
  #$beta.a0.01
  #[1]  0.0000000 -0.6931472
  
  #$beta.a0.02
  #[1] 0.0000000 0.6931472
  
  #$beta.a0.12
  #[1] -0.6931472  0.6931472
  
  #$beta.a1.01
  #[1] -1.386294  1.098612
  
  #$beta.a1.02
  #[1] -0.2876821  0.4054651
  
  #$beta.a1.12
  #[1] 0 0

  #exp.b001 <- exp(Xmat %*% c( 0, -0.6931472))
  #exp.b002 <- exp(Xmat %*% c(0, 0.6931472))
  #exp.b012 <- exp(Xmat %*% c(-0.6931472,  0.6931472))
  #exp.b101 <- exp(Xmat %*% c(-1.386294,  1.098612))
  #exp.b102 <- exp(Xmat %*% c(-0.2876821,  0.4054651))
  #exp.b112 <- exp(Xmat %*% c(0, 0))

  step.A0T1 <- res$H.step.funcs$step.A0T1
  step.A0T2 <- res$H.step.funcs$step.A0T2
  step.A0T12 <- res$H.step.funcs$step.A0T12
  step.A1T1 <- res$H.step.funcs$step.A1T1
  step.A1T2 <- res$H.step.funcs$step.A1T2
  step.A1T12 <- res$H.step.funcs$step.A1T12
  
  A0T1.times <- pmin(knots(step.A0T1), tau) %>% unique
  A0T2.times <- pmin(knots(step.A0T2), tau) %>% unique
  A0T12.times <- pmin(knots(step.A0T12), tau) %>% unique
  A1T1.times <- pmin(knots(step.A1T1), tau) %>% unique
  A1T2.times <- pmin(knots(step.A1T2), tau) %>% unique
  A1T12.times <- pmin(knots(step.A1T12), tau) %>% unique
  
  #theta.est = 3
  theta.est <- (1 - m.A1) * res$thetas[1] + m.A1 * res$thetas[2]
  gamma.scale <- theta.est
  gamma.common.shape <- rho/theta.est 
  gamma.each.shape <- (1 - rho)/theta.est
  T1.0.sim <- T1.1.sim <- T2.0.sim <- T2.1.sim <- vector(length = n.X.vals * n.gamma.vals * n.sample.pers)
  
  for (i in 1:n.X.vals) {
    gamma.common <- rgamma(n.gamma.vals, shape = gamma.common.shape, # gamma.common.shape must be positive
                           scale = gamma.scale)
    gamma0 <- gamma.common + rgamma(n.gamma.vals, shape = gamma.each.shape, 
                                    scale = gamma.scale)
    gamma1 <- gamma.common + rgamma(n.gamma.vals, shape = gamma.each.shape, 
                                    scale = gamma.scale)
    for (j in 1:n.gamma.vals) {
      #print(paste0(i, "  ,  " , j))
      exp.b001.gamma0.sim <- exp.b001[i] * gamma0[j]
      exp.b002.gamma0.sim <- exp.b002[i] * gamma0[j]
      exp.b012.gamma0.sim <- exp.b012[i] * gamma0[j]
      exp.b101.gamma1.sim <- exp.b101[i] * gamma1[j]
      exp.b102.gamma1.sim <- exp.b102[i] * gamma1[j]
      exp.b112.gamma1.sim <- exp.b112[i] * gamma1[j]
      
      step.S.a0.01j <- function(t) {
        exp(-step.A0T1(t) * exp.b001.gamma0.sim)
      }
      step.S.a0.02j <- function(t) {
        exp(-step.A0T2(t) * exp.b002.gamma0.sim)
      }
      step.S.a1.01j <- function(t) {
        exp(-step.A1T1(t) * exp.b101.gamma1.sim)
      }
      step.S.a1.02j <- function(t) {
        exp(-step.A1T2(t) * exp.b102.gamma1.sim)
      }
      
      pr.A0T1 <- -diff(c(1, step.S.a0.01j(A0T1.times)))
      # pr.A0T1[length(pr.A0T1)] <- pr.A0T1[length(pr.A0T1)] + 1 - sum(pr.A0T1)
      #(cbind(A0T1.times, pr.A0T1))
      pr.A0T1 <- c( pr.A0T1, 1 - sum(pr.A0T1) )
      #View(cbind(c(A0T1.times,365), pr.A0T1))
      
      pr.A0T2 <- -diff(c(1, step.S.a0.02j(A0T2.times)))
      # pr.A0T2[length(pr.A0T2)] <- pr.A0T2[length(pr.A0T2)] + 1 - sum(pr.A0T2)
      pr.A0T2 <- c( pr.A0T2, 1 - sum(pr.A0T2) ) 
      
      pr.A1T1 <- -diff(c(1, step.S.a1.01j(A1T1.times)))
      # pr.A1T1[length(pr.A1T1)] <- pr.A1T1[length(pr.A1T1)] +   1 - sum(pr.A1T1)
      pr.A1T1 <- c( pr.A1T1, 1 - sum(pr.A1T1) ) 
      
      pr.A1T2 <- -diff(c(1, step.S.a1.02j(A1T2.times)))
      # pr.A1T2[length(pr.A1T2)] <- pr.A1T2[length(pr.A1T2)] + 1 - sum(pr.A1T2)
      pr.A1T2 <- c( pr.A1T2, 1 - sum(pr.A1T2) )
      
      #a = sample(1:5, size = 100000, replace = T, prob = c(0.78,0.1,0.1,0.01,0.01))
      #table(a) / length(a)
      
      T1.0.sim.temp <- sample(c(A0T1.times, tau), n.sample.pers, replace = T, prob = pr.A0T1)
      T2.0.sim.temp <- sample(c(A0T2.times, tau), n.sample.pers, replace = T, prob = pr.A0T2)
      T1.1.sim.temp <- sample(c(A1T1.times, tau), n.sample.pers, replace = T, prob = pr.A1T1)
      T2.1.sim.temp <- sample(c(A1T2.times, tau), n.sample.pers, replace = T, prob = pr.A1T2)
      disease.0 <- (T1.0.sim.temp < T2.0.sim.temp) # < pmin(T2.0.sim.temp, tau)
      disease.1 <- (T1.1.sim.temp < T2.1.sim.temp) # < pmin(T2.1.sim.temp, tau)
      T1.0.sim.temp[!disease.0] <- Inf
      T1.1.sim.temp[!disease.1] <- Inf
      # T1.0.tau <- T1.0.sim.temp == tau
      # T1.1.tau <- T1.1.sim.temp == tau
      # T2.0.sim.temp[disease.0 & T1.0.tau] <- tau
      # T2.1.sim.temp[disease.1 & T1.1.tau] <- tau
      # disease.0.not.tau <- (disease.0 & !T1.0.tau)
      # disease.1.not.tau <- (disease.1 & !T1.1.tau)
      s0 <- sum(disease.0)
      s1 <- sum(disease.1)
      if (s0 > 0) {
        for (k in 1:s0) {
          step.S.a0.12k <- function(t) {
            exp(-(step.A0T12(t) - step.A0T12(T1.0.sim.temp[disease.0][k])) * 
                  exp.b012.gamma0.sim)
          }
          A0T12.times.k <- A0T12.times[A0T12.times >= T1.0.sim.temp[disease.0][k]]
          pr.A0T12 <- -diff(c(1, step.S.a0.12k(A0T12.times.k)))
          if (length(pr.A0T12) <= 1) {
             T2.0.sim.temp[disease.0][k] <- tau
          }
           else {
            # pr.A0T12[length(pr.A0T12)] <- pr.A0T12[length(pr.A0T12)] + 
            #   1 - sum(pr.A0T12)
            pr.A0T12 <- c( pr.A0T12, 1 - sum(pr.A0T12) )
            T2.0.sim.temp[disease.0][k] <- sample(c(A0T12.times.k, tau), 1, replace = T, prob = pr.A0T12)
          }
        }
      }
      if (s1 > 0) {
        for (k in 1:s1) {
          step.S.a1.12k <- function(t) {
            exp(-(step.A1T12(t) - step.A1T12(T1.1.sim.temp[disease.1][k])) * 
                  exp.b112.gamma1.sim)
          }
          A1T12.times.k <- A1T12.times[A1T12.times >= T1.1.sim.temp[disease.1][k]]
          pr.A1T12 <- -diff(c(1, step.S.a1.12k(A1T12.times.k)))
          if (length(pr.A1T12) <= 1) {
            T2.1.sim.temp[disease.1][k] <- tau
           }
          else {
            # pr.A1T12[length(pr.A1T12)] <- pr.A1T12[length(pr.A1T12)] + 
            #   1 - sum(pr.A1T12)
            pr.A1T12 <- c( pr.A1T12, 1 - sum(pr.A1T12) )
            T2.1.sim.temp[disease.1][k] <- sample(c(A1T12.times.k, tau), 1, replace = T, prob = pr.A1T12)
          }
        }
      }
      st <- (i - 1) * n.gamma.vals * n.sample.pers + (j - 1) * n.sample.pers + 1
      en <- (i - 1) * n.gamma.vals * n.sample.pers + j * n.sample.pers
      T1.0.sim[st:en] <- T1.0.sim.temp
      T2.0.sim[st:en] <- T2.0.sim.temp
      T1.1.sim[st:en] <- T1.1.sim.temp
      T2.1.sim[st:en] <- T2.1.sim.temp
    }
  }
  
  #T_lst_cefa_new  = list(T1.0.sim=T1.0.sim, T2.0.sim=T2.0.sim, T1.1.sim=T1.1.sim, T2.1.sim=T2.1.sim)
  #saveRDS(T_lst_cefa_new  , file = "T_lst_cefa_new.rds")
  
  
  ret <- list(T1.0.sim = T1.0.sim, T1.1.sim = T1.1.sim, T2.0.sim = T2.0.sim, T2.1.sim = T2.1.sim)
    
  return(ret)
}

####################################################################################################33
CalcRMST = function (rho, tau, n.gamma.vals, n.sample.pers
                     #, population, 
                     ,Xnames, Xnames_formula, data, res, list.out = T, detailed = F){
  X_sample_location = c("not_taken", "urine", "wound", "blood", "sputum", "several_sources")
  #n.X.vals <- nrow(population)
  n.X.vals = nrow(data)
  m.A1 <- mean(data$A == 1)
  
  #Xmat <- as.matrix(subset(data, select = Xnames))
  formula0 = as.formula(paste0("Surv(T1, delta1) ~",
       paste0(Xnames_formula, collapse = "+")))
  fit0 <- survival::coxph(formula0, data = data)
  Xmat0 = model.matrix(fit0)
  
  formula1 = formula0
  if(any(X_sample_location %in% Xnames_formula)){
    formula1 = as.formula(paste0("Surv(T1, delta1) ~",
           paste0(Xnames_formula[-match(X_sample_location, Xnames_formula)], collapse = "+")))
  }
 # if("covariate_antib_resist_ceftazidime_90" %in% Xnames_formula){
 #   formula1 = as.formula(paste0("Surv(T1, delta1) ~",
 #      paste0(Xnames_formula[-grep("covariate_antib_resist_ceftazidime_90", Xnames_formula)], collapse = "+")))
 # }
  fit1 <- survival::coxph(formula1, data = data)
  Xmat1 = model.matrix(fit1)
  
  # @IMPO note that fit0/fit1, Xmat0/Xmat1 and formula0/formula1 DO NOT refer to TREATMENT groups 0 and 1!!!. 
  # They refer to the different states - from STATE 0 and STATE 1  
  
  exp.b001 <- exp(Xmat0 %*% coef(res$fit.list$fit.a0.01))
  exp.b002 <- exp(Xmat0 %*% coef(res$fit.list$fit.a0.02))
  exp.b012 <- exp(Xmat1 %*% coef(res$fit.list$fit.a0.12))
  exp.b101 <- exp(Xmat0 %*% coef(res$fit.list$fit.a1.01))
  exp.b102 <- exp(Xmat0 %*% coef(res$fit.list$fit.a1.02))
  exp.b112 <- exp(Xmat1 %*% coef(res$fit.list$fit.a1.12))
  
  step.A0T1 <- res$H.step.funcs$step.A0T1
  step.A0T2 <- res$H.step.funcs$step.A0T2
  step.A0T12 <- res$H.step.funcs$step.A0T12
  step.A1T1 <- res$H.step.funcs$step.A1T1
  step.A1T2 <- res$H.step.funcs$step.A1T2
  step.A1T12 <- res$H.step.funcs$step.A1T12
  
  A0T1.times <- pmin(knots(step.A0T1), tau) %>% unique
  A0T2.times <- pmin(knots(step.A0T2), tau) %>% unique
  A0T12.times <- pmin(knots(step.A0T12), tau) %>% unique
  A1T1.times <- pmin(knots(step.A1T1), tau) %>% unique
  A1T2.times <- pmin(knots(step.A1T2), tau) %>% unique
  A1T12.times <- pmin(knots(step.A1T12), tau) %>% unique
  
  theta.est <- (1 - m.A1) * res$thetas[1] + m.A1 * res$thetas[2]
  gamma.scale <- theta.est
  gamma.common.shape <- rho/theta.est 
  gamma.each.shape <- (1 - rho)/theta.est
  T1.0.sim <- T1.1.sim <- T2.0.sim <- T2.1.sim <- vector(length = n.X.vals * n.gamma.vals * n.sample.pers)
  
  for (i in 1:n.X.vals) {
    gamma.common <- rgamma(n.gamma.vals, shape = gamma.common.shape, # gamma.common.shape must be positive
                           scale = gamma.scale)
    gamma0 <- gamma.common + rgamma(n.gamma.vals, shape = gamma.each.shape, 
                                    scale = gamma.scale)
    gamma1 <- gamma.common + rgamma(n.gamma.vals, shape = gamma.each.shape, 
                                    scale = gamma.scale)
    for (j in 1:n.gamma.vals) {
      #print(paste0(i, "  ,  " , j))
      exp.b001.gamma0.sim <- exp.b001[i] * gamma0[j]
      exp.b002.gamma0.sim <- exp.b002[i] * gamma0[j]
      exp.b012.gamma0.sim <- exp.b012[i] * gamma0[j]
      exp.b101.gamma1.sim <- exp.b101[i] * gamma1[j]
      exp.b102.gamma1.sim <- exp.b102[i] * gamma1[j]
      exp.b112.gamma1.sim <- exp.b112[i] * gamma1[j]
      
      step.S.a0.01j <- function(t) {
        exp(-step.A0T1(t) * exp.b001.gamma0.sim)
      }
      step.S.a0.02j <- function(t) {
        exp(-step.A0T2(t) * exp.b002.gamma0.sim)
      }
      step.S.a1.01j <- function(t) {
        exp(-step.A1T1(t) * exp.b101.gamma1.sim)
      }
      step.S.a1.02j <- function(t) {
        exp(-step.A1T2(t) * exp.b102.gamma1.sim)
      }
      
      pr.A0T1 <- -diff(c(1, step.S.a0.01j(A0T1.times)))
      # pr.A0T1[length(pr.A0T1)] <- pr.A0T1[length(pr.A0T1)] + 1 - sum(pr.A0T1)
      #(cbind(A0T1.times, pr.A0T1))
      pr.A0T1 <- c( pr.A0T1, 1 - sum(pr.A0T1) )
      #View(cbind(c(A0T1.times,365), pr.A0T1))
      
      pr.A0T2 <- -diff(c(1, step.S.a0.02j(A0T2.times)))
      # pr.A0T2[length(pr.A0T2)] <- pr.A0T2[length(pr.A0T2)] + 1 - sum(pr.A0T2)
      pr.A0T2 <- c( pr.A0T2, 1 - sum(pr.A0T2) ) 
      
      pr.A1T1 <- -diff(c(1, step.S.a1.01j(A1T1.times)))
      # pr.A1T1[length(pr.A1T1)] <- pr.A1T1[length(pr.A1T1)] +   1 - sum(pr.A1T1)
      pr.A1T1 <- c( pr.A1T1, 1 - sum(pr.A1T1) ) 
      
      pr.A1T2 <- -diff(c(1, step.S.a1.02j(A1T2.times)))
      # pr.A1T2[length(pr.A1T2)] <- pr.A1T2[length(pr.A1T2)] + 1 - sum(pr.A1T2)
      pr.A1T2 <- c( pr.A1T2, 1 - sum(pr.A1T2) )
      
      #a = sample(1:5, size = 100000, replace = T, prob = c(0.78,0.1,0.1,0.01,0.01))
      #table(a) / length(a)
      
      T1.0.sim.temp <- sample(c(A0T1.times, tau), n.sample.pers, replace = T, prob = pr.A0T1)
      T2.0.sim.temp <- sample(c(A0T2.times, tau), n.sample.pers, replace = T, prob = pr.A0T2)
      T1.1.sim.temp <- sample(c(A1T1.times, tau), n.sample.pers, replace = T, prob = pr.A1T1)
      T2.1.sim.temp <- sample(c(A1T2.times, tau), n.sample.pers, replace = T, prob = pr.A1T2)
      disease.0 <- (T1.0.sim.temp < T2.0.sim.temp) # < pmin(T2.0.sim.temp, tau)
      disease.1 <- (T1.1.sim.temp < T2.1.sim.temp) # < pmin(T2.1.sim.temp, tau)
      T1.0.sim.temp[!disease.0] <- Inf
      T1.1.sim.temp[!disease.1] <- Inf
      # T1.0.tau <- T1.0.sim.temp == tau
      # T1.1.tau <- T1.1.sim.temp == tau
      # T2.0.sim.temp[disease.0 & T1.0.tau] <- tau
      # T2.1.sim.temp[disease.1 & T1.1.tau] <- tau
      # disease.0.not.tau <- (disease.0 & !T1.0.tau)
      # disease.1.not.tau <- (disease.1 & !T1.1.tau)
      s0 <- sum(disease.0)
      s1 <- sum(disease.1)
      if (s0 > 0) {
        for (k in 1:s0) {
          step.S.a0.12k <- function(t) {
            exp(-(step.A0T12(t) - step.A0T12(T1.0.sim.temp[disease.0][k])) * 
                  exp.b012.gamma0.sim)
          }
          A0T12.times.k <- A0T12.times[A0T12.times >= T1.0.sim.temp[disease.0][k]]
          pr.A0T12 <- -diff(c(1, step.S.a0.12k(A0T12.times.k)))
          if (length(pr.A0T12) <= 1) {
            T2.0.sim.temp[disease.0][k] <- tau
          }
          else {
            # pr.A0T12[length(pr.A0T12)] <- pr.A0T12[length(pr.A0T12)] + 
            #   1 - sum(pr.A0T12)
            pr.A0T12 <- c( pr.A0T12, 1 - sum(pr.A0T12) )
            T2.0.sim.temp[disease.0][k] <- sample(c(A0T12.times.k, tau), 1, replace = T, prob = pr.A0T12)
          }
        }
      }
      if (s1 > 0) {
        for (k in 1:s1) {
          step.S.a1.12k <- function(t) {
            exp(-(step.A1T12(t) - step.A1T12(T1.1.sim.temp[disease.1][k])) * 
                  exp.b112.gamma1.sim)
          }
          A1T12.times.k <- A1T12.times[A1T12.times >= T1.1.sim.temp[disease.1][k]]
          pr.A1T12 <- -diff(c(1, step.S.a1.12k(A1T12.times.k)))
          if (length(pr.A1T12) <= 1) {
            T2.1.sim.temp[disease.1][k] <- tau
          }
          else {
            # pr.A1T12[length(pr.A1T12)] <- pr.A1T12[length(pr.A1T12)] + 
            #   1 - sum(pr.A1T12)
            pr.A1T12 <- c( pr.A1T12, 1 - sum(pr.A1T12) )
            T2.1.sim.temp[disease.1][k] <- sample(c(A1T12.times.k, tau), 1, replace = T, prob = pr.A1T12)
          }
        }
      }
      st <- (i - 1) * n.gamma.vals * n.sample.pers + (j - 1) * n.sample.pers + 1
      en <- (i - 1) * n.gamma.vals * n.sample.pers + j * n.sample.pers
      T1.0.sim[st:en] <- T1.0.sim.temp
      T2.0.sim[st:en] <- T2.0.sim.temp
      T1.1.sim[st:en] <- T1.1.sim.temp
      T2.1.sim[st:en] <- T2.1.sim.temp
    }
  }
  
  #T_lst_cefa_new  = list(T1.0.sim=T1.0.sim, T2.0.sim=T2.0.sim, T1.1.sim=T1.1.sim, T2.1.sim=T2.1.sim)
  #saveRDS(T_lst_cefa_new  , file = "T_lst_cefa_new.rds")
  
  
  ret <- list(T1.0.sim = T1.0.sim, T1.1.sim = T1.1.sim, T2.0.sim = T2.0.sim, T2.1.sim = T2.1.sim)
  
  return(ret)
}
####################################################################################################

####################################################################################################
CalcRMST_tweak = function (rho, tau, n.gamma.vals, n.sample.pers
                     #, population, 
                     ,Xnames, Xnames_formula, data, res, list.out = T, detailed = F){

  m.A1 <- mean(data$A == 1)
  X_tweak <- as.matrix(qnorm(c(0.2, 0.4, 0.6, 0.8)))
  #n.X.vals = nrow(data)
  n.X.vals = nrow(X_tweak)
  
  #Xmat <- as.matrix(subset(data, select = Xnames))
  formula0 = as.formula(paste0("Surv(T1, delta1) ~",
                               paste0(Xnames_formula, collapse = "+")))
  fit0 <- survival::coxph(formula0, data = data)
  #Xmat0 = model.matrix(fit0)
  Xmat0 = X_tweak
  
  formula1 = formula0
  fit1 <- survival::coxph(formula1, data = data)
  #Xmat1 = model.matrix(fit1)
  Xmat1 = X_tweak
  
  # @IMPO note that fit0/fit1, Xmat0/Xmat1 and formula0/formula1 DO NOT refer to TREATMENT groups 0 and 1!!!. 
  # They refer to the different states - from STATE 0 and STATE 1  
  
  exp.b001 <- exp(Xmat0 %*% coef(res$fit.list$fit.a0.01))
  exp.b002 <- exp(Xmat0 %*% coef(res$fit.list$fit.a0.02))
  exp.b012 <- exp(Xmat1 %*% coef(res$fit.list$fit.a0.12))
  exp.b101 <- exp(Xmat0 %*% coef(res$fit.list$fit.a1.01))
  exp.b102 <- exp(Xmat0 %*% coef(res$fit.list$fit.a1.02))
  exp.b112 <- exp(Xmat1 %*% coef(res$fit.list$fit.a1.12))
  
  #exp.b001 <- exp(Xmat0 * c(-0.6931472))
  #exp.b002 <- exp(Xmat0 * c(0.6931472))
  #exp.b012 <- exp(Xmat1 * c(0.6931472))
  #exp.b101 <- exp(Xmat0 * c(1.098612))
  #exp.b102 <- exp(Xmat0 * c(0.4054651))
  #exp.b112 <- exp(Xmat1 * c(0))
  
  
  step.A0T1 <- res$H.step.funcs$step.A0T1
  step.A0T2 <- res$H.step.funcs$step.A0T2
  step.A0T12 <- res$H.step.funcs$step.A0T12
  step.A1T1 <- res$H.step.funcs$step.A1T1
  step.A1T2 <- res$H.step.funcs$step.A1T2
  step.A1T12 <- res$H.step.funcs$step.A1T12
  
  A0T1.times <- pmin(knots(step.A0T1), tau) %>% unique
  A0T2.times <- pmin(knots(step.A0T2), tau) %>% unique
  A0T12.times <- pmin(knots(step.A0T12), tau) %>% unique
  A1T1.times <- pmin(knots(step.A1T1), tau) %>% unique
  A1T2.times <- pmin(knots(step.A1T2), tau) %>% unique
  A1T12.times <- pmin(knots(step.A1T12), tau) %>% unique
  
  #theta.est = 2
  theta.est <- (1 - m.A1) * res$thetas[1] + m.A1 * res$thetas[2]
  gamma.scale <- theta.est
  gamma.common.shape <- rho/theta.est 
  gamma.each.shape <- (1 - rho)/theta.est
  T1.0.sim <- T1.1.sim <- T2.0.sim <- T2.1.sim <- vector(length = n.X.vals * n.gamma.vals * n.sample.pers)
  
  for (i in 1:n.X.vals) {
    gamma.common <- rgamma(n.gamma.vals, shape = gamma.common.shape, # gamma.common.shape must be positive
                           scale = gamma.scale)
    gamma0 <- gamma.common + rgamma(n.gamma.vals, shape = gamma.each.shape, 
                                    scale = gamma.scale)
    gamma1 <- gamma.common + rgamma(n.gamma.vals, shape = gamma.each.shape, 
                                    scale = gamma.scale)
    for (j in 1:n.gamma.vals) {
      #print(paste0(i, "  ,  " , j))
      exp.b001.gamma0.sim <- exp.b001[i] * gamma0[j]
      exp.b002.gamma0.sim <- exp.b002[i] * gamma0[j]
      exp.b012.gamma0.sim <- exp.b012[i] * gamma0[j]
      exp.b101.gamma1.sim <- exp.b101[i] * gamma1[j]
      exp.b102.gamma1.sim <- exp.b102[i] * gamma1[j]
      exp.b112.gamma1.sim <- exp.b112[i] * gamma1[j]
      
      step.S.a0.01j <- function(t) {
        exp(-step.A0T1(t) * exp.b001.gamma0.sim)
      }
      step.S.a0.02j <- function(t) {
        exp(-step.A0T2(t) * exp.b002.gamma0.sim)
      }
      step.S.a1.01j <- function(t) {
        exp(-step.A1T1(t) * exp.b101.gamma1.sim)
      }
      step.S.a1.02j <- function(t) {
        exp(-step.A1T2(t) * exp.b102.gamma1.sim)
      }
      
      pr.A0T1 <- -diff(c(1, step.S.a0.01j(A0T1.times)))
      # pr.A0T1[length(pr.A0T1)] <- pr.A0T1[length(pr.A0T1)] + 1 - sum(pr.A0T1)
      #(cbind(A0T1.times, pr.A0T1))
      pr.A0T1 <- c( pr.A0T1, 1 - sum(pr.A0T1) )
      #View(cbind(c(A0T1.times,365), pr.A0T1))
      
      pr.A0T2 <- -diff(c(1, step.S.a0.02j(A0T2.times)))
      # pr.A0T2[length(pr.A0T2)] <- pr.A0T2[length(pr.A0T2)] + 1 - sum(pr.A0T2)
      pr.A0T2 <- c( pr.A0T2, 1 - sum(pr.A0T2) ) 
      
      pr.A1T1 <- -diff(c(1, step.S.a1.01j(A1T1.times)))
      # pr.A1T1[length(pr.A1T1)] <- pr.A1T1[length(pr.A1T1)] +   1 - sum(pr.A1T1)
      pr.A1T1 <- c( pr.A1T1, 1 - sum(pr.A1T1) ) 
      
      pr.A1T2 <- -diff(c(1, step.S.a1.02j(A1T2.times)))
      # pr.A1T2[length(pr.A1T2)] <- pr.A1T2[length(pr.A1T2)] + 1 - sum(pr.A1T2)
      pr.A1T2 <- c( pr.A1T2, 1 - sum(pr.A1T2) )
      
      #a = sample(1:5, size = 100000, replace = T, prob = c(0.78,0.1,0.1,0.01,0.01))
      #table(a) / length(a)
      
      T1.0.sim.temp <- sample(c(A0T1.times, tau), n.sample.pers, replace = T, prob = pr.A0T1)
      T2.0.sim.temp <- sample(c(A0T2.times, tau), n.sample.pers, replace = T, prob = pr.A0T2)
      T1.1.sim.temp <- sample(c(A1T1.times, tau), n.sample.pers, replace = T, prob = pr.A1T1)
      T2.1.sim.temp <- sample(c(A1T2.times, tau), n.sample.pers, replace = T, prob = pr.A1T2)
      disease.0 <- (T1.0.sim.temp < T2.0.sim.temp) # < pmin(T2.0.sim.temp, tau)
      disease.1 <- (T1.1.sim.temp < T2.1.sim.temp) # < pmin(T2.1.sim.temp, tau)
      T1.0.sim.temp[!disease.0] <- Inf
      T1.1.sim.temp[!disease.1] <- Inf
      # T1.0.tau <- T1.0.sim.temp == tau
      # T1.1.tau <- T1.1.sim.temp == tau
      # T2.0.sim.temp[disease.0 & T1.0.tau] <- tau
      # T2.1.sim.temp[disease.1 & T1.1.tau] <- tau
      # disease.0.not.tau <- (disease.0 & !T1.0.tau)
      # disease.1.not.tau <- (disease.1 & !T1.1.tau)
      s0 <- sum(disease.0)
      s1 <- sum(disease.1)
      if (s0 > 0) {
        for (k in 1:s0) {
          step.S.a0.12k <- function(t) {
            exp(-(step.A0T12(t) - step.A0T12(T1.0.sim.temp[disease.0][k])) * 
                  exp.b012.gamma0.sim)
          }
          A0T12.times.k <- A0T12.times[A0T12.times >= T1.0.sim.temp[disease.0][k]]
          pr.A0T12 <- -diff(c(1, step.S.a0.12k(A0T12.times.k)))
          if (length(pr.A0T12) <= 1) {
            T2.0.sim.temp[disease.0][k] <- tau
          }
          else {
            # pr.A0T12[length(pr.A0T12)] <- pr.A0T12[length(pr.A0T12)] + 
            #   1 - sum(pr.A0T12)
            pr.A0T12 <- c( pr.A0T12, 1 - sum(pr.A0T12) )
            T2.0.sim.temp[disease.0][k] <- sample(c(A0T12.times.k, tau), 1, replace = T, prob = pr.A0T12)
          }
        }
      }
      if (s1 > 0) {
        for (k in 1:s1) {
          step.S.a1.12k <- function(t) {
            exp(-(step.A1T12(t) - step.A1T12(T1.1.sim.temp[disease.1][k])) * 
                  exp.b112.gamma1.sim)
          }
          A1T12.times.k <- A1T12.times[A1T12.times >= T1.1.sim.temp[disease.1][k]]
          pr.A1T12 <- -diff(c(1, step.S.a1.12k(A1T12.times.k)))
          if (length(pr.A1T12) <= 1) {
            T2.1.sim.temp[disease.1][k] <- tau
          }
          else {
            # pr.A1T12[length(pr.A1T12)] <- pr.A1T12[length(pr.A1T12)] + 
            #   1 - sum(pr.A1T12)
            pr.A1T12 <- c( pr.A1T12, 1 - sum(pr.A1T12) )
            T2.1.sim.temp[disease.1][k] <- sample(c(A1T12.times.k, tau), 1, replace = T, prob = pr.A1T12)
          }
        }
      }
      st <- (i - 1) * n.gamma.vals * n.sample.pers + (j - 1) * n.sample.pers + 1
      en <- (i - 1) * n.gamma.vals * n.sample.pers + j * n.sample.pers
      T1.0.sim[st:en] <- T1.0.sim.temp
      T2.0.sim[st:en] <- T2.0.sim.temp
      T1.1.sim[st:en] <- T1.1.sim.temp
      T2.1.sim[st:en] <- T2.1.sim.temp
    }
  }
  
  ret <- list(T1.0.sim = T1.0.sim, T1.1.sim = T1.1.sim, T2.0.sim = T2.0.sim, T2.1.sim = T2.1.sim, theta.est = theta.est)
  
  return(ret)
}
####################################################################################################




effect_calc_outer <- function(t, stratum_indices, T1.0.sim, T1.1.sim, T2.0.sim, T2.1.sim) {
  F1_0 = mean(T1.0.sim[stratum_indices] <= t)
  F1_1 = mean(T1.1.sim[stratum_indices] <= t)
  #F2_0 = mean(T2.0.sim[stratum_indices] <= t)
  #F2_1 = mean(T2.1.sim[stratum_indices] <= t)
  effect_diff = mean(F1_1) - mean(F1_0)
  effect_ratio = mean(F1_1) / mean(F1_0)
  return(list(F1_1=F1_1, F1_0=F1_0
              #, F2_1=F2_1, F2_0=F2_0
              , effect_diff=effect_diff
              , effect_ratio=effect_ratio))
}



CalcRMST_calculations = function(T1.0.sim, T1.1.sim, T2.0.sim, T2.1.sim,
                                 time_points=c(1:365)){
    last_point = max(time_points)
    ad <- (T1.0.sim < T2.0.sim) & (T1.1.sim < T2.1.sim)
    #ad2 <- (T1.0.sim <= T2.0.sim & T1.0.sim !=T2.0.sim) & (T1.1.sim <= T2.1.sim & T1.1.sim !=T2.1.sim)
    nd <- (T1.0.sim >= T2.0.sim) & (T1.1.sim >= T2.1.sim)
    prop.ad <- mean(ad)
    prop.nd <- mean(nd)
    
    as <- (T2.0.sim == last_point) & (T2.1.sim == last_point)
    prop.as <- mean(as)
    
    ios0 <- ((T2.0.sim == last_point) | (T1.0.sim < T2.0.sim)) 
    ios1 <-  ((T2.1.sim == last_point) | (T1.1.sim < T2.1.sim))
    ios <- ios0 & ios1
    prop.ios0 <- mean(ios0)
    prop.ios1 <- mean(ios1)
    prop.ios <- mean(ios)
    
    # FICE_1 = mean(T1.1.sim[ios] < 365) - mean(T1.0.sim[ios] < 365)
    # SACE_1 = mean(T1.1.sim[as] < 365) - mean(T1.0.sim[as] < 365)
    # AICE_1 = mean(T1.1.sim[ad] < 100) - mean(T1.0.sim[ad] < 100)
    
    effect_calc <- function(t, stratum_indices) {
      F1_1 = mean(T1.1.sim[stratum_indices] <= t)
      F1_0 = mean(T1.0.sim[stratum_indices] <= t)
      F2_1 = mean(T2.1.sim[stratum_indices] <= t)
      F2_0 = mean(T2.0.sim[stratum_indices] <= t)
      effect_diff = mean(F1_1) - mean(F1_0)
      effect_ratio = mean(F1_1) / mean(F1_0)
      return(list(F1_1=F1_1, F1_0=F1_0, F2_1=F2_1, F2_0=F2_0, 
                  effect_diff=effect_diff, effect_ratio=effect_ratio))
    }
    
    #time_points = c(1:365) # seq(0,365,5) # c(1:365)
    print("start calculating FICE_1_vec in MC")
    FICE_1_vec = sapply(time_points, effect_calc, ios)
    
    #FICE_1_vec_cefa = FICE_1_vec
    #saveRDS(FICE_1_vec_cefa  , file = "FICE_1_vec_cefa.rds")
    
    d = data.frame(time = time_points, t(as.data.frame(FICE_1_vec)))
    d = sapply(d,as.numeric) %>% data.frame()
    dd <- melt(data.table(d), id.vars="time")
    
    
    return(list(nd=nd, ad=ad, as=as, ios0=ios0, ios1=ios1, ios=ios, 
                prop.nd=prop.nd, prop.ad=prop.ad, prop.as=prop.as, 
                prop.ios0=prop.ios0, prop.ios1=prop.ios1, prop.ios=prop.ios, 
                FICE_1_vec=FICE_1_vec, d=d, dd=dd))
}

