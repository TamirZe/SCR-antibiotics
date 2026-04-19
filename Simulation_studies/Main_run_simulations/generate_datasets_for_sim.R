library(dplyr)
library(survival)
library(devtools)
library(ggplot2)
library(CausalSemiComp)

# set the path according to the configuration of your local system
setwd("~/R projects/AAA PhD/Causal-effects-on-non-terminal-event-time-with-application")
source("Estimands_numerical_example/Main_run/GetScenarioParams_paper.R") 
source("Estimands_numerical_example/Main_run/run_functions.R") 
source("Estimands_numerical_example/Data_generation/SimDataWeibFrail.R")
source("Estimands_numerical_example/Estimand_calculations/estimands_functions.R")
source("Estimands_numerical_example/Estimand_calculations/CalcTrueCausalParams.R")

scen_seed = 50001; lett = "a"
scen_ind = 7 # 1-3: scen 5, 4-6: scen 9, 7-9: scen 7
tau = 1; n.sample = 4000; dim_x = 1 # n.sample = 10000 # n.sample = 30000000 
set.seed(scen_seed + match(lett,letters))

#scen_values = c(5,9) # main text
scen_values = c(5,9,7) # go to GetScenarioParams_new to see the details of these scenarios- 5:theta=1, 9:theta=3, 7:theta=2
rho_values = c(0,0.5,1)
scen_vec = rep(scen_values, each=length(rho_values))
rho_vec = rep(rho_values, length(unique(scen_vec)))
estimands_mat = NULL
coeff_A = c(0.2, 0.2) # coeff_A = NULL

##########################################################
data_lst <- param_lst <- list()
for (i in 1:500) {
  print(i)
  #params <- GetScenarioParams_new(scenario.num = scen) inside my.data_generation_tweak
  set.seed(101 + i)
  my.data_generation = Data_generation_func(scen=scen_vec[scen_ind], rho=rho_vec[scen_ind], tau=tau,  # scen=7
                                            n.sample=n.sample, dim_x=dim_x, coeff_A=coeff_A)
  my.data = my.data_generation$my.data
  
  params = my.data_generation$params
  data_lst[[i]] <- my.data
  param_lst[[i]] <- params
}

saveRDS(param_lst, file = paste0("Estimands_numerical_example/Estimation_clstr/data_lsts/Output_clstr/S", scen_ind,
                                 "/simulation_data_params_lst_n", n.sample, "_S", scen_ind, ".rds"))
saveRDS(data_lst, file = paste0("Estimands_numerical_example/Estimation_clstr/data_lsts/Output_clstr/S", scen_ind,
                                "/simulation_data_lst_n", n.sample, "_S", scen_ind, ".rds"))

