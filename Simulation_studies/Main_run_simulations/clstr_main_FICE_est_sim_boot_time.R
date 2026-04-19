args = commandArgs(trailingOnly = TRUE)
print(args)
print(length(args))

job_id = NULL
if(length(args)>0){
  job_id = as.numeric(args[1])
}else{
  stop(' ERROR no job id given!!!')
}

print("one")
library(dplyr)
library(data.table)
library(survival)
library(reshape)
library(future.apply)
library(progressr)

# Set up parallel processing
n_plan <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
#plan(multicore, workers = n_plan)
future::plan(future::multisession, workers = n_plan) # workers = n_plan # 6
#future::plan(future::multicore, workers = n_plan) 
print("multisession")

# parameter rho
###############################################################
n_sample=4000
scen = "S7"
###############################################################

#print("two")
###############################################################
# source for Simulations_studies
main_path = paste0("/scratch200/tamirzehavi/AntibioticsSCR/")
path_estimation_code = paste0(main_path, "Estimation/Code/Scripts_estimation_for_simulations/")
#/scratch200/tamirzehavi/AntibioticsSCR/Estimation/Data/Output/S1
#path_estimation_data_out = paste0(main_path, "Estimation/Data/Output/", scen, "/tweak_n", n_sample, "/")
source(paste0(path_estimation_code, "EMcausalSC_GH_current_FICE_est_time.R"))
source(paste0(path_estimation_code, "CalcRMST_FICE_est.R"))
# --- FIX: ?????? ???? ?????? ??? ---
path_estimation_data_out <- paste0(
  "/scratch200/tamirzehavi/AntibioticsSCR/Estimation/Data/Output/",
  scen, "/tweak_n", n_sample, "/"
)
#dir.create(path_estimation_data_out, recursive = TRUE, showWarnings = FALSE)
###############################################################


# generaed by numerical_example_main_text_estimation
simulation_data_lst <- readRDS(paste0(main_path, "Estimation/Data/Output/", scen, "/simulation_data_lst_n", n_sample, "_", scen, ".rds"))
param_lst           <- readRDS(paste0(main_path, "Estimation/Data/Output/", scen, "/simulation_data_params_lst_n", n_sample, "_", scen, ".rds"))
orig_data <- simulation_data_lst[[job_id]] 
rho = param_lst[[1]]$rho

nboot = 200 # 200 # < length(simulation_data_lst)
message("DEBUG NBOOT=", nboot)
Xnames <- Xnames_formula <- "X"
match_bool = FALSE
#EM
max.iter = 10000; eps.conv = 0.00001  #  0.0001 #  0.00001 # 0.000001
# eps.conv = 0.0001
#MC
n_gamma = 100; n_quad_PO = 100 # 100
reduced_return = TRUE
###############################################################

# Bootstrap/main sim estimation function

run_estimation <- function(boot_iter=NULL, data, rho, Xnames, Xnames_formula, 
                           match_bool = FALSE, 
                           max.iter = 10000, eps.conv = 0.0001,
                           n_gamma, n_quad_PO, mc_tweak, reduced_return,
                           em_time = TRUE, log_fun = NULL, max_time_sec = Inf, check_timeout = NULL) {
  
  
  # Create bootstrapped dataset or use original data
  if(is.null(boot_iter)){ 
    data_b = data
  } else{
    # bootstrap - set seed for reproducibility
    set.seed(101 + boot_iter)
    if(match_bool == TRUE) {
      pairs_b = sort(sample(unique(data$pair), length(unique(data$pair)), replace = TRUE))
      data_b = data[unlist(sapply(pairs_b, function(p) which(data$pair == p))), ] %>% arrange(pair)
    } else {
      rows_b = sample(1:nrow(data), nrow(data), replace = TRUE)
      data_b = data[rows_b, ]
    }
  }
  
  # max.iter = 10000, w = NULL, eps.conv = 0.0001
  # EM estimation
  stage <- "init"
  t_boot0 <- proc.time()[3]
  check_timeout <- function(stage = "") {
  elapsed <- proc.time()[3] - t_boot0
  if (elapsed > max_time_sec) {
      stop(sprintf("TIMEOUT: %.1fs exceeded at %s", elapsed, stage))
    }
  }
  print("EM")
  stage <- "EM"
  set.seed(101)
  res = EMcausalSC_FICE_est(data = data_b, Xnames = Xnames, Xnames_formula = Xnames_formula, Lname = NULL,
                            max.iter = max.iter, w = NULL, eps.conv = eps.conv, 
                            init.thetas = c(1, 1),
                            track_time = em_time, log_fun = log_fun, max_time_sec = max_time_sec, check_timeout = check_timeout)
  
  #num_EM_iterations = res$iter
  beta_EM = res$betas
  thetas_EM = res$thetas
  
  check_timeout("after EM / before MC")
  
  # MCMC simulation
  print("MC")
  stage <- "MC"
  t_mc0 <- proc.time()[3]
  Calc <- if (mc_tweak==F) CalcRMST else CalcRMST_tweak
  set.seed(101)
  sim_data = Calc(rho = rho, tau = 1, n.gamma.vals = n_gamma, n.sample.pers = n_quad_PO, 
                      Xnames = Xnames, Xnames_formula = Xnames_formula, data = data, 
                      res = res, list.out = TRUE, detailed = FALSE)
                      if (!is.null(log_fun)) {
  log_fun("MC finished in %.1f sec", proc.time()[3] - t_mc0)}
  stage <- "post-MC"
  check_timeout("after MC")

  
  # Extract simulation results 
  T1.0.sim = sim_data$T1.0.sim
  T1.1.sim = sim_data$T1.1.sim
  T2.0.sim = sim_data$T2.0.sim
  T2.1.sim = sim_data$T2.1.sim

  if(reduced_return == TRUE){
    return(list(
    boot_iter = boot_iter,
    num_EM_iterations = res$iter,
    em_timing = res$em_timing,
    beta_EM = beta_EM,
    thetas_EM = thetas_EM,
    theta.est = sim_data$theta.est,
    T1.0.sim = T1.0.sim,
    T1.1.sim = T1.1.sim,
    T2.0.sim = T2.0.sim,
    T2.1.sim = T2.1.sim
    ))
  }
  # Calculate estimators
  print("CalcRMST_calculations")
  estimators = CalcRMST_calculations(T1.0.sim = T1.0.sim, T1.1.sim = T1.1.sim,
                                     T2.0.sim = T2.0.sim, T2.1.sim = T2.1.sim,
                                     time_points = c(1:365)/365) # time_points = c(1:365) # seq(0,1, by=(1/365))
  
  # Calculate proportions
  as = estimators$as; prop.as = estimators$prop.as; ad = estimators$ad
  prop.ad = estimators$prop.ad; 
  ad0 <- (T1.0.sim < T2.0.sim); ad1 <- (T1.1.sim < T2.1.sim); 
  prop.ad0 = mean(ad0); prop.ad1 = mean(ad1)
  ad_try = ad0 & ad1; prop.ad_try = mean(ad_try)
  prop.nd = estimators$prop.nd
  
  ios = estimators$ios
  ios0 = T1.0.sim <= T2.0.sim | T2.0.sim >= 1
  ios1 = T1.1.sim <= T2.1.sim | T2.1.sim >= 1
  prop.ios0 = estimators$prop.ios0
  prop.ios1 = estimators$prop.ios1
  prop.ios = estimators$prop.ios
  
  ut0 = ((T1.0.sim <= T2.0.sim & T2.0.sim < 1) & (T1.1.sim > T2.1.sim & T2.1.sim >= 1))
  ut1 = ((T1.1.sim <= T2.1.sim & T2.1.sim < 1) & (T1.0.sim > T2.0.sim & T2.0.sim >= 1))
  ut = ut0 | ut1
  prop.ut0 = mean(ut0)
  prop.ut1 = mean(ut1)
  prop.ut = prop.ut0 + prop.ut1
  
  frailty_proportions_lst = list(
    prop.ios0 = prop.ios0, prop.ios1 = prop.ios1, prop.ios = prop.ios,
    prop.as = prop.as, prop.ad0 = prop.ad0, prop.ad1 = prop.ad1, prop.ad = prop.ad,
    prop.ut0 = prop.ut0, prop.ut1 = prop.ut1, prop.ut = prop.ut
  )
  
  # Calculate one-year estimators
  FICE_1_vec = estimators$FICE_1_vec
  T_2_ios_effect = estimators$d$F2_1 - estimators$d$F2_0
  # < 365
  FICE_1 = mean(T1.1.sim[ios] < 1) - mean(T1.0.sim[ios] < 1)
  SACE_1 = mean(T1.1.sim[as] < 1) - mean(T1.0.sim[as] < 1)
  AICE_1 = mean(T1.1.sim[ad] < 1) - mean(T1.0.sim[ad] < 1)
  
  ut0_CE = mean(T1.1.sim[ut0] < 1) - mean(T1.0.sim[ut0] < 1)
  ut1_CE = mean(T1.1.sim[ut1] < 1) - mean(T1.0.sim[ut1] < 1)
  ut_CE = mean(T1.1.sim[ut] < 1) - mean(T1.0.sim[ut] < 1)
  ut_CE_direct = (prop.ut1 - prop.ut0) / (prop.ut)
  
  frailty_estimators_one_year_lst = list(
    FICE_1_vec = FICE_1_vec, T_2_ios_effect = T_2_ios_effect, FICE_1 = FICE_1,
    SACE_1 = SACE_1, AICE_1 = AICE_1,
    ut0_CE = ut0_CE, ut1_CE = ut1_CE, ut_CE = ut_CE, ut_CE_direct = ut_CE_direct
  )
  
  # Calculate time-varying estimators
  times_calc = c(1:365) / 365
  est_ad_lst = list()
  est_as_lst = list()
  est_ut_lst = list()
  est_ut0_lst = list()
  est_ut1_lst = list()
  
  for (i in 1:length(times_calc)) {
    est_ad_lst[[i]] = effect_calc_outer(times_calc[i], ad, T1.0.sim, T1.1.sim, T2.0.sim, T2.1.sim)
    est_as_lst[[i]] = effect_calc_outer(times_calc[i], as, T1.0.sim, T1.1.sim, T2.0.sim, T2.1.sim)
    est_ut_lst[[i]] = effect_calc_outer(times_calc[i], ut, T1.0.sim, T1.1.sim, T2.0.sim, T2.1.sim)
    est_ut0_lst[[i]] = effect_calc_outer(times_calc[i], ut0, T1.0.sim, T1.1.sim, T2.0.sim, T2.1.sim)
    est_ut1_lst[[i]] = effect_calc_outer(times_calc[i], ut1, T1.0.sim, T1.1.sim, T2.0.sim, T2.1.sim)
  }
  
  from_lst_to_mat_for_figure = function(times, est_lst) {
    df = data.frame(time = times, t(do.call("cbind", est_lst)))
    df = sapply(df, as.numeric) %>% data.frame()
    ddf <- melt(data.table(df), id.vars = "time")
    return(ddf)
  }
  
  ddf_ios = estimators$dd
  ddf_ad = from_lst_to_mat_for_figure(times = times_calc, est_lst = est_ad_lst)
  ddf_as = from_lst_to_mat_for_figure(times = times_calc, est_lst = est_as_lst)
  ddf_ut = from_lst_to_mat_for_figure(times = times_calc, est_lst = est_ut_lst)
  ddf_ut0 = from_lst_to_mat_for_figure(times = times_calc, est_lst = est_ut0_lst)
  ddf_ut1 = from_lst_to_mat_for_figure(times = times_calc, est_lst = est_ut1_lst)
  
  ddf_ios$variable = paste0("ios_", ddf_ios$variable)
  ddf_as$variable = paste0("as_", ddf_as$variable)
  ddf_ad$variable = paste0("ad_", ddf_ad$variable)
  ddf_ut$variable = paste0("ut_", ddf_ut$variable)
  
  all_estimands = data.frame(rbind(ddf_ios, ddf_as, ddf_ad, ddf_ut))
  all_estimands_diff_plt = all_estimands %>%
    filter(variable %in% c("ios_effect_diff", "as_effect_diff", "ad_effect_diff", "ut_effect_diff"))
  all_estimands_ratio_plt = all_estimands %>%
    filter(variable %in% c("ios_effect_ratio", "as_effect_ratio", "ad_effect_ratio", "ut_effect_ratio"))
  
  frailty_estimators_lst = list(
    ddf_ios = ddf_ios, ddf_as = ddf_as, ddf_ad = ddf_ad, ddf_ut = ddf_ut,
    all_estimands = all_estimands,
    all_estimands_diff_plt = all_estimands_diff_plt,
    all_estimands_ratio_plt = all_estimands_ratio_plt
  )
  
  # Calculate total effects
  T1_F_mean_vec_0 = sapply(times_calc, function(x) mean(T1.0.sim <= x & T1.0.sim < T2.0.sim))
  T1_F_mean_vec_1 = sapply(times_calc, function(x) mean(T1.1.sim <= x & T1.1.sim < T2.1.sim))
  
  total_effect_ratio = T1_F_mean_vec_1 / T1_F_mean_vec_0
  total_effect_diff = T1_F_mean_vec_1 - T1_F_mean_vec_0
  
  T2_F_mean_vec_0 = sapply(times_calc, function(x) mean(T2.0.sim <= x))
  T2_F_mean_vec_1 = sapply(times_calc, function(x) mean(T2.1.sim <= x))
  
  T2_effect_ratio = T2_F_mean_vec_1 / T2_F_mean_vec_0
  T2_effect_diff = T2_F_mean_vec_1 - T2_F_mean_vec_0
  
  # Calculate bounds
  T2_F_among_infected = mean(T2.0.sim[ad0] < 1)
  T2_F_among_non_infected = mean(T2.0.sim[!ad0] < 1)
  
  # Bounds without ORP
  U_pi_ios_wout_ORP = min(prop.ios0, prop.ios1)
  L_pi_ios_wout_ORP = max(0, prop.ios0 + prop.ios1 - 1)
  
  U_P1_t_wout_ORP = pmin(T1_F_mean_vec_1, prop.ios0) / L_pi_ios_wout_ORP
  L_P1_t_wout_ORP = pmax(0, (T1_F_mean_vec_1 + prop.ios0 - 1) / U_pi_ios_wout_ORP)
  
  U_P0_t_wout_ORP = pmin(T1_F_mean_vec_0, prop.ios1) / L_pi_ios_wout_ORP
  L_P0_t_wout_ORP = pmax(0, (T1_F_mean_vec_0 + prop.ios1 - 1) / U_pi_ios_wout_ORP)
  
  U_t_wout_ORP_diff = U_P1_t_wout_ORP - L_P0_t_wout_ORP
  L_t_wout_ORP_diff = L_P1_t_wout_ORP - U_P0_t_wout_ORP
  
  U_t_wout_ORP_ratio = (pmin(T1_F_mean_vec_1, prop.ios0)) / (pmax(0, T1_F_mean_vec_0 + prop.ios1 - 1))
  L_t_wout_ORP_ratio = (pmax(0, T1_F_mean_vec_1 + prop.ios0 - 1)) / (pmin(T1_F_mean_vec_0, prop.ios1))
  
  u_t_wout_ORP_diff = (pmin(T1_F_mean_vec_1, prop.ios0)) - (pmax(0, T1_F_mean_vec_0 + prop.ios1 - 1))
  l_t_wout_ORP_diff = (pmax(0, prop.ios0 - 1 + T1_F_mean_vec_1)) - (pmin(T1_F_mean_vec_0, prop.ios1))
  
  U2_t_wout_ORP_diff = u_t_wout_ORP_diff / ifelse(u_t_wout_ORP_diff >= 0, L_pi_ios_wout_ORP, U_pi_ios_wout_ORP)
  L2_t_wout_ORP_diff = l_t_wout_ORP_diff / ifelse(l_t_wout_ORP_diff >= 0, U_pi_ios_wout_ORP, L_pi_ios_wout_ORP)
  
  # Bounds under weak-ORP
  U_pi_ios_weak_ORP = pmin(prop.ios0, prop.ios1 + (T1_F_mean_vec_0[365] * T2_F_among_infected))
  L_pi_ios_weak_ORP = max(T1_F_mean_vec_0[365], prop.ios0 + prop.ios1 - 1)
  
  U_P0_t_weak_ORP = pmin(1, T1_F_mean_vec_0 / L_pi_ios_weak_ORP)
  L_P0_t_weak_ORP = T1_F_mean_vec_0 / U_pi_ios_weak_ORP
  U_P1_t_weak_ORP = pmin(T1_F_mean_vec_1, prop.ios0) / L_pi_ios_weak_ORP
  L_P1_t_weak_ORP = pmax(0, T1_F_mean_vec_1 + prop.ios0 - 1) / U_pi_ios_weak_ORP
  
  U_t_weak_ORP_diff = U_P1_t_weak_ORP - L_P0_t_weak_ORP
  L_t_weak_ORP_diff = L_P1_t_weak_ORP - U_P0_t_weak_ORP
  
  U_t_weak_ORP_ratio = (pmin(T1_F_mean_vec_1, prop.ios0) / T1_F_mean_vec_0)
  L_t_weak_ORP_ratio = pmax(0, T1_F_mean_vec_1 + prop.ios0 - 1) / T1_F_mean_vec_0
  
  u_t_weak_ORP_diff = (pmin(T1_F_mean_vec_1, prop.ios0)) - (T1_F_mean_vec_0)
  l_t_weak_ORP_diff = (pmax(0, prop.ios0 - 1 + T1_F_mean_vec_1)) - (T1_F_mean_vec_0)
  
  U2_t_weak_ORP_diff = u_t_weak_ORP_diff / ifelse(u_t_weak_ORP_diff >= 0, L_pi_ios_weak_ORP, U_pi_ios_weak_ORP)
  L2_t_weak_ORP_diff = l_t_weak_ORP_diff / ifelse(l_t_weak_ORP_diff >= 0, U_pi_ios_weak_ORP, L_pi_ios_weak_ORP)
  
  # Bounds under ios-ORP
  pi_ios_ORP_IOS = prop.ios0
  P_0_t_ios_ORP = T1_F_mean_vec_0 / prop.ios0
  U_P1_t_ios_ORP = pmin(1, T1_F_mean_vec_1 / prop.ios0)
  L_P1_t_ios_ORP = pmax(0, 1 - ((1 - T1_F_mean_vec_1) / prop.ios0))
  U_t_ios_ORP_diff = U_P1_t_ios_ORP - P_0_t_ios_ORP
  L_t_ios_ORP_diff = L_P1_t_ios_ORP - P_0_t_ios_ORP
  
  bounds_proportions_lst = list(
    L_pi_ios_wout_ORP = L_pi_ios_wout_ORP, U_pi_ios_wout_ORP = U_pi_ios_wout_ORP,
    L_pi_ios_weak_ORP = L_pi_ios_weak_ORP, U_pi_ios_weak_ORP = U_pi_ios_weak_ORP,
    pi_ios_ORP_IOS = pi_ios_ORP_IOS,
    prop.ios0 = prop.ios0, prop.ios1 = prop.ios1
  )
  
  bounds_estimators_diff_lst = list(
    L_t_wout_ORP_diff = L_t_wout_ORP_diff, U_t_wout_ORP_diff = U_t_wout_ORP_diff,
    L2_t_wout_ORP_diff = L2_t_wout_ORP_diff, U2_t_wout_ORP_diff = U2_t_wout_ORP_diff,
    L_t_weak_ORP_diff = L_t_weak_ORP_diff, U_t_weak_ORP_diff = U_t_weak_ORP_diff,
    L2_t_weak_ORP_diff = L2_t_weak_ORP_diff, U2_t_weak_ORP_diff = U2_t_weak_ORP_diff,
    L_t_ios_ORP_diff = L_t_ios_ORP_diff, U_t_ios_ORP_diff = U_t_ios_ORP_diff
  )
  
  bounds_estimators_ratio_lst = list(
    L_t_wout_ORP_ratio = L_t_wout_ORP_ratio, U_t_wout_ORP_ratio = U_t_wout_ORP_ratio,
    L_t_weak_ORP_ratio = L_t_weak_ORP_ratio, U_t_weak_ORP_ratio = U_t_weak_ORP_ratio
  )
  
  # Return all results
  list(
    boot_iter = boot_iter,
    num_EM_iterations = res$iter,
    em_timing = res$em_timing,
    beta_EM = beta_EM,
    thetas_EM = thetas_EM,
    frailty_proportions_lst = frailty_proportions_lst,
    frailty_estimators_one_year_lst = frailty_estimators_one_year_lst,
    frailty_estimators_lst = frailty_estimators_lst,
    total_effect_ratio = total_effect_ratio,
    total_effect_diff = total_effect_diff,
    T2_effect_ratio = T2_effect_ratio,
    T2_effect_diff = T2_effect_diff,
    bounds_proportions_lst = bounds_proportions_lst,
    bounds_estimators_diff_lst = bounds_estimators_diff_lst,
    bounds_estimators_ratio_lst = bounds_estimators_ratio_lst
  )
}

# Point estmation in one ieration of simulated dataset
sim_results = run_estimation(boot_iter = NULL, 
               data = orig_data,  
               rho = rho,
               Xnames = Xnames, 
               Xnames_formula = Xnames_formula,
               match_bool = FALSE,
               max.iter = max.iter, eps.conv = eps.conv,
               n_gamma = n_gamma, n_quad_PO = n_quad_PO, 
               mc_tweak = TRUE, reduced_return = reduced_return)


f_main <- paste0(path_estimation_data_out, "sim_results_sim", job_id, "_scen", scen, ".rds")
message("OUTDIR=", path_estimation_data_out)
message("f_main=", f_main)
message("getwd=", getwd())
flush(stdout())
#saveRDS(sim_results, f_main)
#cat("WROTE", f_main, "size", file.info(f_main)$size, "\n")
saveRDS(sim_results, f_main)
message("SAVED: ", f_main)
stopifnot(file.exists(f_main), file.info(f_main)$size > 0)
#saveRDS(sim_results, 
#        file = paste0(path_estimation_data_out, 'sim_results_', 'sim', job_id, '_scen', scen, '.rds')) 
        
#f_main <- paste0(path_estimation_data_out, 'sim_results_', 'sim', job_id, '_scen', scen, '.rds')
#stopifnot(file.exists(f_main))
      
 
#stop("Only main sim WOUT BOOTSTRAP")      
       
# Run bootstrap with progress bar
print(paste0("Starting ", nboot, " bootstrap iterations with parallel processing..."))

# Use progressr for progress tracking
handlers(global = TRUE)
handlers("void")
#handlers("txtprogressbar")


#with_progress({
#  p <- progressor(steps = nboot)
#  bootstrap_results <- future_lapply(1:nboot, function(j) {
#    p(sprintf("Bootstrap iteration %d", j))
#    run_estimation(boot_iter = j, 
#                  data = orig_data,  
#                  rho = rho,
#                  Xnames = Xnames, 
#                  Xnames_formula = Xnames_formula,
#                  match_bool = FALSE,
#                  n_gamma=n_gamma, n_quad_PO=n_quad_PO)
#  }, future.seed = TRUE)
#})





log_out <- stdout()  # this is exactly the file set by: #SBATCH -o .../%x_%A_%a.out

log_line <- function(fmt, ...) {
  cat(sprintf(fmt, ...), "\n", file = log_out)
  flush.console()  # usually harmless; helps when stdout is buffered
}

log_line("[%s] SLURM job start | job_id=%s array_task=%s host=%s",
         Sys.time(),
         Sys.getenv("SLURM_JOB_ID"),
         Sys.getenv("SLURM_ARRAY_TASK_ID"),
         Sys.info()[["nodename"]])


with_progress({
  p <- progressor(along = 1:nboot)  # one tick per bootstrap

  bootstrap_results <- future_lapply(
    1:nboot,
    function(j) {
    library(survival)
      t0 <- proc.time()[3]

      # per-bootstrap logger
      log_fun <- function(fmt, ...) {
        log_line("[%s] BOOT %d/%d %s",
                 Sys.time(), j, nboot, sprintf(fmt, ...))
      }

      log_fun("START")
      p()  # ok (handlers("void") = no spam)

      res <- tryCatch(
        run_estimation(
          boot_iter = j,
          data = orig_data,
          rho = rho,
          Xnames = Xnames,
          Xnames_formula = Xnames_formula,
          match_bool = FALSE,
          n_gamma = n_gamma,
          n_quad_PO = n_quad_PO,
          max.iter = max.iter,
          eps.conv = eps.conv,
          mc_tweak = TRUE,
          reduced_return = reduced_return,
          # NEW:
          em_time = TRUE,
          log_fun = log_fun,
          max_time_sec = 60 * 60,      # e.g. 30 minutes per bootstrap
                         # collect + log EM iter timing
          check_timeout = NULL
        ),
        #error = function(e) {
        #  log_fun("FAILED: %s", conditionMessage(e))
        #  return(NULL)
        #}
     error = function(e) {

        now <- proc.time()[3]
      
        if (exists("t0")) {
          log_fun("FAILED after %.1f sec | stage=%s",
                  now - t0,
                  if (exists("stage")) stage else "unknown")
        } else {
          log_fun("FAILED | stage=%s",
                  if (exists("stage")) stage else "unknown")
        }
      
        if (exists("t_mc0") && exists("stage") && stage == "MC") {
          log_fun("MC elapsed before failure: %.1f sec",
                  now - t_mc0)
          log_fun("EM elapsed before MC: %.1f sec",
                  t_mc0 - t0)
        }
      
        log_fun("ERROR: %s", conditionMessage(e))
        log_fun("ERR_CLASS: %s", paste(class(e), collapse = ","))
      
        cl <- conditionCall(e)
        if (!is.null(cl))
          log_fun("ERR_CALL: %s", paste(deparse(cl), collapse=" "))
      
        if (requireNamespace("rlang", quietly = TRUE)) {
          tb <- rlang::trace_back()
          tb_txt <- paste(utils::capture.output(print(tb)), collapse = " | ")
          log_fun("TRACE: %s", tb_txt)
        } else {
          calls_txt <- paste(
            vapply(sys.calls(),
                   function(x) paste(deparse(x), collapse=" "),
                   ""),
            collapse=" | ")
          log_fun("SYS.CALLS: %s", calls_txt)
        }
      
        return(NULL)
  }

      )

      t1 <- proc.time()[3]
      #dt_min <- as.numeric(difftime(t1, t0, units = "mins"))
      dt_min <- (t1 - t0) / 60
      log_fun("END (%.2f min)", dt_min)

      res
    },
    future.seed = TRUE
  )
})



log_line("[%s] SLURM job end | job_id=%s array_task=%s",
         Sys.time(),
         Sys.getenv("SLURM_JOB_ID"),
         Sys.getenv("SLURM_ARRAY_TASK_ID"))



print("Bootstrap completed!")

# Save bootstrap results
#saveRDS(bootstrap_results, 
#        file = paste0(path_estimation_data_out, 'bootstrap_results_', 'sim', job_id, '_scen', scen, '.rds'))


#f_boot <- paste0(path_estimation_data_out, "bootstrap_results_sim", job_id, "_scen", scen, ".rds")
#message("OUTDIR=", path_estimation_data_out)
#message("f_main=", f_boot)
#message("getwd=", getwd())
#flush(stdout())

#saveRDS(bootstrap_results, f_boot)
#cat("WROTE", f_boot, "size", file.info(f_boot)$size, "\n")

#saveRDS(bootstrap_results, f_boot)
#message("SAVED: ", f_boot)
#stopifnot(file.exists(f_boot), file.info(f_boot)$size > 0)

# ===== ????? bootstrap?? ?????? ??? ??? =====


print("Bootstrap completed!")

f_boot <- paste0(
  path_estimation_data_out,
  "bootstrap_results_sim", job_id, "_scen", scen, ".rds"
)

bootstrap_results_clean <-
  bootstrap_results[!vapply(bootstrap_results, is.null, logical(1))]

message("BOOTSTRAP: kept ",
        length(bootstrap_results_clean),
        " out of ", length(bootstrap_results))

if (length(bootstrap_results_clean) > 0) {
  saveRDS(bootstrap_results_clean, f_boot)
  message("SAVED: ", f_boot)
  stopifnot(file.exists(f_boot), file.info(f_boot)$size > 0)
} else {
  message("NOT SAVED: all bootstrap iterations failed")
}


#saveRDS(bootstrap_results, f_boot)
#message("SAVED: ", f_boot)
#stopifnot(file.exists(f_boot), file.info(f_boot)$size > 0)
#print("All results saved!")