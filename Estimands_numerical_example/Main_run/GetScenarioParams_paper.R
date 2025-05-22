GetScenarioParams_new = function(scenario.num, rho){
  
  #TODO scenario A, reflects proportions, similar FICE and AICE
  if (scenario.num == 5) {
    base.weib.scale.a1.01 = 2
    base.weib.scale.a1.02 = 2.75
    base.weib.scale.a1.12 = 2.25
    base.weib.scale.a0.01 = 2.5
    base.weib.scale.a0.02 = 2.25
    base.weib.scale.a0.12 = 2.75
    
    base.weib.shape.a1.01 = 2.5
    base.weib.shape.a1.02 = 2.1
    base.weib.shape.a1.12 = 2.1
    base.weib.shape.a0.01 = 2.5
    base.weib.shape.a0.02 = 2.1
    base.weib.shape.a0.12 = 2.1
    
    beta.a1.01 <- log(c(0.25, 3))
    beta.a1.02 <- log(c(0.75, 1.5))
    beta.a1.12 <- log(c(1, 1))
    beta.a0.01 <- -log(c(1, 2))
    beta.a0.02 <- log(c(1, 2))
    beta.a0.12 <- log(c(0.5, 2))
    k.tau <- 0.25
    theta <- 3 * k.tau/(1 - k.tau)
    rho = 0
  }
  
  if (scenario.num == 9) {
    base.weib.scale.a1.01 = 2
    base.weib.scale.a1.02 = 2.75
    base.weib.scale.a1.12 = 2.25
    base.weib.scale.a0.01 = 2.5
    base.weib.scale.a0.02 = 2.25
    base.weib.scale.a0.12 = 2.75
    
    base.weib.shape.a1.01 = 2.5
    base.weib.shape.a1.02 = 2.1
    base.weib.shape.a1.12 = 2.1
    base.weib.shape.a0.01 = 2.5
    base.weib.shape.a0.02 = 2.1
    base.weib.shape.a0.12 = 2.1
    
    beta.a1.01 <- log(c(0.25, 3))
    beta.a1.02 <- log(c(0.75, 1.5))
    beta.a1.12 <- log(c(1, 1))
    beta.a0.01 <- -log(c(1, 2))
    beta.a0.02 <- log(c(1, 2))
    beta.a0.12 <- log(c(0.5, 2))
    k.tau <- 0.25
    theta <- 9 * k.tau/(1 - k.tau)
    rho = 0
  }
  
  #TODO scenario B, DO NOT reflect strata proportions, different FICE and AICE
  if (scenario.num == 11) {
    base.weib.scale.a1.01 = 0.1
    base.weib.scale.a1.02 = 1
    base.weib.scale.a1.12 = 1
    base.weib.scale.a0.01 = 0.1
    base.weib.scale.a0.02 = 1
    base.weib.scale.a0.12 = 1
    
    base.weib.shape.a1.01 = 1
    base.weib.shape.a1.02 = 3
    base.weib.shape.a1.12 = 3
    base.weib.shape.a0.01 = 1
    base.weib.shape.a0.02 = 0.5
    base.weib.shape.a0.12 = 0.5
    
    beta.a1.01 <- log(c(1, 2))
    beta.a1.02 <- log(c(1, 2))
    beta.a1.12 <- log(c(1, 2))
    beta.a0.01 <- -log(c(1, 2))
    beta.a0.02 <- -log(c(1, 2))
    beta.a0.12 <- -log(c(1, 2))
    
    k.tau <- 0.25
    theta <- 3 * k.tau/(1 - k.tau)
    rho = 0
  }
  
  if (scenario.num == 13) {
    base.weib.scale.a1.01 = 0.1
    base.weib.scale.a1.02 = 1
    base.weib.scale.a1.12 = 1
    base.weib.scale.a0.01 = 0.1
    base.weib.scale.a0.02 = 1
    base.weib.scale.a0.12 = 1
    
    base.weib.shape.a1.01 = 1
    base.weib.shape.a1.02 = 3
    base.weib.shape.a1.12 = 3
    base.weib.shape.a0.01 = 1
    base.weib.shape.a0.02 = 0.5
    base.weib.shape.a0.12 = 0.5
    
    beta.a1.01 <- log(c(1, 2))
    beta.a1.02 <- log(c(1, 2))
    beta.a1.12 <- log(c(1, 2))
    beta.a0.01 <- -log(c(1, 2))
    beta.a0.02 <- -log(c(1, 2))
    beta.a0.12 <- -log(c(1, 2))
    
    k.tau <- 0.25
    theta <- 9 * k.tau/(1 - k.tau)
    rho = 0
  }
  
  params <- mget(ls(environment()))
  return(params)
}
