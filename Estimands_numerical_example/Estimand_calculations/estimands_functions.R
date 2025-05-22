# estimands functions
# SACE_FICE ####  
point_SACE_FICE_estimand_func = function(dataset, t, r, estimand){ # r<=t
  if(estimand=="SACE"){
    dataset$stratum_r = (dataset$T2.0>=r) & (dataset$T2.1>=r)
  }
  if(estimand=="FICE"){
    dataset$stratum_r = ( dataset$T1.0<=r | dataset$T2.0>=r ) & 
                        ( dataset$T1.1<=r | dataset$T2.1>=r )
    # dataset$stratum_r = ( dataset$T1.0<dataset$T2.0 | dataset$T2.0>r ) & 
    #                     ( dataset$T1.1<dataset$T2.1 | dataset$T2.1>r )
    # dataset$stratum_r = ( (dataset$T1.0<r & dataset$T1.0<dataset$T2.0) | dataset$T2.0>r ) & 
    #                     ( (dataset$T1.1<r & dataset$T1.1<dataset$T2.1) | dataset$T2.1>r )
  }
  if(estimand=="two_utype"){
    dataset$stratum_r = dataset$two_utype
  }
    
    
  prop_r = mean(dataset$stratum_r)
  inf1_r = mean(dataset$T1.1[dataset$stratum_r] <= t)
  inf0_r = mean(dataset$T1.0[dataset$stratum_r] <= t)
  estimand_r = inf1_r - inf0_r
  nonsurv1_r = mean(dataset$T2.1[dataset$stratum_r] <= t)
  nonsurv0_r = mean(dataset$T2.0[dataset$stratum_r] <= t)
  effect_on_death_r = nonsurv1_r - nonsurv0_r
  return(c(t=t, r=r, prop_r=prop_r, 
           inf1_r=inf1_r, inf0_r=inf0_r, estimand_r=estimand_r,
           nonsurv1_r=nonsurv1_r, nonsurv0_r=nonsurv0_r, effect_on_death_r=effect_on_death_r))
}

TV_SACE_FICE_estimand_func = function(dataset, r, len, estimand){ # fixed r value
  names_col = c("t", "r", "prop", "inf1", "inf0", "effect",
                "T2_1", "T2_0", "effect_on_death")
  t_points = seq(0, r, len)
  df_estimand_fixed_r = matrix(nrow = length(t_points), ncol = length(names_col))
  for (i in 1:length(t_points)) { # run on t values
    df_estimand_fixed_r[i, ] = 
      point_SACE_FICE_estimand_func(dataset=dataset, t=t_points[i], r=r, estimand=estimand)
  }
  colnames(df_estimand_fixed_r) = names_col
  return(df_estimand_fixed_r)
}

# AICE #### 
AICE_estimand_func = function(dataset, t_max, len, r=NULL){ 
  names_col = c("t", "r", "prop", "inf1", "inf0", "effect",
                "T2_1", "T2_0", "effect_on_death")
  if(is.null(r)){
    dataset$ai = (dataset$T1.0<dataset$T2.0) & (dataset$T1.1<dataset$T2.1)
  }else{
    dataset$ai = (dataset$T1.0<r & dataset$T1.0<dataset$T2.0) & 
      (dataset$T1.1<r & dataset$T1.1<dataset$T2.1) 
  }
  prop_ai = mean(dataset$ai)
  t_points = seq(0, t_max, len)
  df_AICE = matrix(nrow = length(t_points), ncol = length(names_col))
  for (i in 1:length(t_points)) { # run on t values
    t = t_points[i]
    inf1_ai = mean(dataset$T1.1[dataset$ai] <= t)
    inf0_ai = mean(dataset$T1.0[dataset$ai] <= t)
    AICE = inf1_ai - inf0_ai
    nonsurv1_ai = mean(dataset$T2.1[dataset$ai] <= t)
    nonsurv0_ai = mean(dataset$T2.0[dataset$ai] <= t)
    effect_on_death_ai = nonsurv1_ai - nonsurv0_ai
    df_AICE[i, ] = c(t=t, r=0, prop_ai=prop_ai, 
     inf1_ai=inf1_ai, inf0_ai=inf0_ai, AICE=AICE,
     nonsurv1=nonsurv1_ai, nonsurv0=nonsurv0_ai, effect_on_death=effect_on_death_ai)
  }
  colnames(df_AICE) = names_col
  return(df_AICE)
}

