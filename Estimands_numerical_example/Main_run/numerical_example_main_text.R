library(dplyr)
library(survival)
library(devtools)
library(ggplot2)
library(CausalSemiComp)

#library(Daniel)
#library(ggpubr)
#library(reshape2)
#library(boot)

ddd= data.frame(a=1:3, b=4:6)
ccc = ddd %>% filter(a==1)

# set the path according to the configuration of your local system
setwd("~/R projects/AAA PhD/Causal-effects-on-non-terminal-event-time-with-application")
source("Estimands_numerical_example/Main_run/GetScenarioParams_paper.R") 
source("Estimands_numerical_example/Main_run/run_functions.R") 
source("Estimands_numerical_example/Data_generation/SimDataWeibFrail.R")
source("Estimands_numerical_example/Estimand_calculations/estimands_functions.R")
source("Estimands_numerical_example/Estimand_calculations/CalcTrueCausalParams.R")

#Sys.unsetenv("GITHUB_PAT")
#install_github("daniel258/CausalSemiComp")


scen_seed = 50001; lett = "a"
scen = 1; tau = 1; n.sample = 10000; dim_x = 1 # n.sample = 30000000 
set.seed(scen_seed + match(lett,letters))

scen_values = c(5,9) # main text
rho_values = c(0,0.5,1)
scen_vec = rep(scen_values, each=length(rho_values))
rho_vec = rep(rho_values, length(unique(scen_vec)))
estimands_mat = NULL
# sum_my.data_lst: parameters set -> rho value -> a,b,c
sum_my.data_lst <- replicate(length(scen_values), setNames(replicate(length(rho_values), vector("list", 0), 
           simplify = FALSE), paste0("rho=", rho_values)), simplify = FALSE)

for (i in 1:length(scen_vec)) {
  print(i)
  my.data_generation = Data_generation_func(scen=scen_vec[i], rho=rho_vec[i], tau=tau, n.sample=n.sample, dim_x=dim_x)
  my.data = my.data_generation$my.data
  params = my.data_generation$params
  
  all_estimands_mat_scen = create_estimands_mat(scen=as.factor(scen_vec[i]), params=params,
                                                my.data=my.data, r_vec=seq(1,1,0.5))
  estimands_mat = rbind(estimands_mat, all_estimands_mat_scen)
  
  # first argument - scen (theta in our case)
  # second argument - rho
  sum_my.data_lst[[ match(scen_vec, scen_values)[i] ]][[ first(which(rho_vec==rho_vec[i])) ]] = 
    sum_data_func(my.data=my.data)
}

main_text_sum_my.data_lst = sum_my.data_lst
#save(main_text_sum_my.data_lst, file = "main_text_sum_my.data_lst.RData")
summary_my.data_lst <- lapply(main_text_sum_my.data_lst, function(sublist) {
  lapply(sublist, function(innerlist) {
    innerlist[[2]]  # extract the 2nd element
  })
})

# all estimands
all_estimands_mat_fig = estimands_mat %>% 
  subset(select = c("Scen", "theta", "theta_label", "rho", "rho_label", "t", "r", "estimand", "effect")) %>% 
  filter(estimand %in% c("SACE", "AICE", "FICE") & r %in% c(0,1))

##############################################################
# plot estimands 
all_estimands_selected_r_first =  
  ggplot(all_estimands_mat_fig, aes(x = t, y = effect)) + 
  geom_line(aes(color = estimand), size=1.5) + 
  facet_grid(theta_label ~ rho_label, labeller = label_parsed) +
  labs(y = "Causal effect on infection") + 
  xlab("Time from baseline (years)") +
  ylim(c(min(all_estimands_mat_fig$effect), max(all_estimands_mat_fig$effect))) + 
  scale_colour_manual(name = "Estimand"
                      , values=c("tomato", "springgreen3", "deepskyblue2") # slateblue1
                      , labels=c("AICE(t)","FICE(t)", "SACE(t)"))  +
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  theme(plot.title = element_text(size=18),
          axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 10),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          legend.text = element_text(size = 11), 
          legend.title = element_text(size = 14), 
          legend.key.size = unit(2, 'cm'),
          strip.text.x = element_text(size = 16, face = "bold"),
          strip.text.y = element_text(size = 16, face = "bold"))  
all_estimands_selected_r_first
##############################################################

