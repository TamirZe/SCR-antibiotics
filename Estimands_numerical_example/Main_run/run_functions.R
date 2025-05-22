Data_generation_func = function(scen, rho, tau, n.sample, dim_x){
  ##############################################################
  # set parameters ####
  # GetScenarioParams_new from GetScenarioParams3.R
  params <- GetScenarioParams_new(scenario.num = scen)
  params$rho = rho
  R <- 100
  #lapply(params,round,2)
  # high scale (lambda = "1/lambda" in exponential dist) - high mean T
  # high alpha (shape) - low mean T
  # increase T1.0 (base.weib.scale.a0.01) - increasing the effect (hi I(0))
  # increase T1.1 (base.weib.scale.a1.01) - reduce the effect (low I(1))
  ##############################################################
  
  ##############################################################
  # generate the dataset ####
  Daniel::CatIndex(1)
  sim.df <- SimDataWeibFrail(n.sample = n.sample, params = params,
                             no.protected = F, no.large = F,
                             cens.exp.rate = 0.0295, cens = F, cens.admin = 100, 
                             round.times = F, dim_x = dim_x, tau=tau)
  
  my.data = data.frame(
    T1.0.orig = sim.df$T1.0, T2.0.orig = sim.df$T2.0
    ,T1.1.orig = sim.df$T1.1, T2.1.orig = sim.df$T2.1
    ,T1.0 = sim.df$T1.0, T2.0 = sim.df$T2.0
    ,T1.1 = sim.df$T1.1, T2.1 = sim.df$T2.1
    ,delta1.0 = sim.df$delta1.0, delta2.0 = sim.df$delta2.0
    ,delta1.1 = sim.df$delta1.1, delta2.1 = sim.df$delta2.1
    #,T1 = sim.df$T1, T2 = sim.df$T2,
    #,delta1 = sim.df$delta1, delta2 = sim.df$delta2
    ,A = sim.df$A, X = sim.df$X
    ,I_0 = ifelse(sim.df$T1.0 <= tau & sim.df$T1.0 <= sim.df$T2.0, 1, 0), S_0 = ifelse(sim.df$T2.0 >= tau, 1, 0)
    ,I_1 = ifelse(sim.df$T1.1 <= tau & sim.df$T1.1 <= sim.df$T2.1, 1, 0), S_1 = ifelse(sim.df$T2.1 >= tau, 1, 0))
  
  #adaptations of T1 and T2 ####
  # inf
  my.data$T1.0[my.data$T1.0 > my.data$T2.0 | my.data$T1.0 > 1] <- Inf
  my.data$T1.1[my.data$T1.1 > my.data$T2.1 | my.data$T1.1 > 1] <- Inf
  # tau "RMST"
  if (exists("tau")) {
    my.data$T2.0 <- pmin(my.data$T2.0, tau)
    my.data$T2.1 <- pmin(my.data$T2.1, tau)
  }
  
  
  # stratum ####
  my.data$ai = as.numeric(  (my.data$T1.0 <= my.data$T2.0) & (my.data$T1.1 <= my.data$T2.1) )
  # ai = as.numeric( (my.data$T1.0 < 1) & (my.data$T1.1 < 1) )
  my.data$as = as.numeric( (my.data$T2.0 >= 1) & (my.data$T2.1 >= 1) )
  my.data$ios = as.numeric( ( my.data$T1.0 <= 1 | my.data$T2.0 >= 1 ) & 
                              ( my.data$T1.1 <= 1 | my.data$T2.1 >= 1 ) )
  identical(my.data$S_0 * my.data$S_1, my.data$as)
  identical(my.data$I_0 * my.data$I_1, my.data$ai)
  ##############################################################
  
  my.data$two_pt = my.data$as == 0 &  my.data$ai == 0 & my.data$ios == 1
  return(list(my.data=my.data, params=params))
}
##############################################################

##############################################################
# estimands ####
create_estimands_mat = function(scen, params, my.data, r_vec = seq(1,1,0.5)){
  rho = params$rho
  theta=params$theta
  
  # SACE_estimands
  SACE_mat = NULL
  # iterate over several r values 
  for (j in 1:length(r_vec)) {
    SACE_mat = rbind(SACE_mat, 
                     TV_SACE_FICE_estimand_func(dataset=my.data, r=r_vec[j], len=0.01, estimand="SACE"))
  }
  SACE_mat = data.frame(SACE_mat)
  fixed_time_SACE = filter(SACE_mat, t==r)
  point_SACE_t_r1 = filter(SACE_mat, r==1)
  
  # FICE_estimands
  FICE_mat = NULL
  for (j in 1:length(r_vec)) {
    FICE_mat = rbind(FICE_mat, 
                     TV_SACE_FICE_estimand_func(dataset=my.data, r=r_vec[j], len=0.01, estimand="FICE"))
  }
  FICE_mat = data.frame(FICE_mat)
  fixed_time_FICE = filter(FICE_mat, t==r)
  point_FICE_t_r1 = filter(FICE_mat, r==1)
  
  # two_pt_estimands
  two_pt_mat = NULL
  for (j in 1:length(r_vec)) {
    two_pt_mat = rbind(two_pt_mat, 
                       TV_SACE_FICE_estimand_func(dataset=my.data, r=r_vec[j], len=0.01, estimand="two_pt"))
  }
  two_pt_mat = data.frame(two_pt_mat)
  fixed_time_two_pt = filter(two_pt_mat, t==r)
  point_two_pt_t_r1 = filter(two_pt_mat, r==1)
  
  # AICE_estimands
  AICE_mat = AICE_estimand_func(dataset = my.data, t_max=1, len=0.01) %>% # , r=1
    data.frame()
  AICE_t1 = filter(AICE_mat, t==1)
  
  
  # all_estimands_mat
  all_estimands_mat_scen = 
    data.frame( Scen=as.factor(scen), theta=theta, theta_label=as.factor(paste0("theta == ", round(theta, 2)))
                ,rho=rho, rho_label=as.factor(paste0("rho == ", rho))
                ,rbind(data.frame(SACE_mat, estimand="SACE") 
                       ,data.frame(FICE_mat, estimand="FICE")
                       ,data.frame(AICE_mat, estimand="AICE")
                       ,data.frame(two_pt_mat, estimand="two_pt")) )
  #all_estimands_r1_mat = all_estimands_mat %>% filter(r %in% c(0,1))
  return(all_estimands_mat_scen = all_estimands_mat_scen)
}
##############################################################

##############################################################
sum_data_func = function(my.data){
  return(list(a = summary(my.data),
              b = my.data %>% summarise_all("mean"), #  %>% group_by(A)
              c = my.data %>% group_by(A) %>% summarise_all("mean") %>% print(width=Inf)))
}
##############################################################