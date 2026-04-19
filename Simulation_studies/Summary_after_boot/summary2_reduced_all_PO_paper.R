library(dplyr)
library(data.table)
library(survival)
library(reshape)
library(future.apply)
library(progressr)
library(ggplot2)

# source from numerical_example_main_text_estimation
setwd("~/R projects/AAA PhD/Causal-effects-on-non-terminal-event-time-with-application")
base_path <- "~/R projects/AAA PhD/Causal-effects-on-non-terminal-event-time-with-application/Estimands_numerical_example/Estimation_clstr/data_lsts/Output_clstr"

source("Estimands_numerical_example/Main_run/GetScenarioParams_paper.R") 
source("Estimands_numerical_example/Main_run/run_functions.R") 
source("Estimands_numerical_example/Data_generation/SimDataWeibFrail.R")
source("Estimands_numerical_example/Estimand_calculations/estimands_functions.R")
source("Estimands_numerical_example/Estimand_calculations/CalcTrueCausalParams.R")
source("Estimands_numerical_example/Simulation_studies/Summary_after_boot/functions_calculations_sim.R")
# estimation simulations results

########################################################
n.sample = 4000
scen_ind = 7

n_sim = 30000000
X = as.matrix(qnorm(c(0.2, 0.4, 0.6, 0.8)))

# parameters from numerical_example_main_text_estimation
scen_seed = 50001; lett = "a"
tau = 1; dim_x = dim(X)[2] 
set.seed(scen_seed + match(lett,letters))

scen_values = c(5,9,7) # main text # 5: theta=1, 7: theta=2, 9: theta=3
#scen_values = c(5,7,9)
rho_values = c(0,0.5,1)
scen_vec = rep(scen_values, each=length(rho_values))
rho_vec = rep(rho_values, length(unique(scen_vec)))
coeff_A = c(0.2, 0.2) # coeff_A = NULL
########################################################

# true values for each time point (from large population)
times_calc=seq(0,365,14)/365
times_calc = unique(c(times_calc, 1))


#######################################################################
setwd("~/R projects/AAA PhD/Causal-effects-on-non-terminal-event-time-with-application")
base_path <- "~/R projects/AAA PhD/Causal-effects-on-non-terminal-event-time-with-application/Estimands_numerical_example/Estimation_clstr/data_lsts/Output_clstr"

# loop over scenarios and create datasets
scen_inds <- c(1:3, 7:9)

for (scen_ind in scen_inds) {
  message("=== S", scen_ind, " ===")
  
  setwd(file.path(base_path, paste0("S", scen_ind, "/tweak_n", n.sample)))
  files <- list.files(pattern = "\\.rds$")
  main_files <- files[grepl("^sim_results", files)]
  
  #### Population true values for each time point (from large population)
  set.seed(101)
  my.data_generation_tweak = Data_generation_func(scen=scen_vec[scen_ind], rho=rho_vec[scen_ind], tau=tau,
                                                  n.sample=n_sim, X=X, dim_x=dim_x, coeff_A=coeff_A)
  my.data_tweak = my.data_generation_tweak$my.data
  dat_T.sim_tweak = data.frame(T1.0.sim = my.data_tweak$T1.0, T1.1.sim = my.data_tweak$T1.1,
                               T2.0.sim = my.data_tweak$T2.0, T2.1.sim = my.data_tweak$T2.1)
  #dat_T.sim_tweak = readRDS(file.path(base_path, paste0("S", scen_ind, "/tweak_n", n.sample), 
  #                                     paste0("dat_T.sim_tweak_S", scen_ind, ".rds")))
  target_param_by_t <- sapply(times_calc, function(tc)
    calc_3(d = dat_T.sim_tweak, last_point = max(times_calc), t_cut = tc)
  )
  
  saveRDS(dat_T.sim_tweak, file = paste0("dat_T.sim_tweak_S", scen_ind, ".rds"))
  rm(my.data_generation_tweak); rm(my.data_tweak); rm(dat_T.sim_tweak); gc()
  saveRDS(target_param_by_t, file = paste0("target_param_by_t_S", scen_ind, ".rds"))
  rm(target_param_by_t)
  
  # Estimation - point_est_list
  
  point_est_list <- lapply(times_calc, function(tc)
    matrix(NA_real_, nrow = length(main_files), ncol = 7,
           dimnames = list(NULL, c("FICE","SACE","AICE","TE","pi_ios","pi_as","pi_ad")))
  )
  names(point_est_list) <- times_calc
  
  for (i in seq_along(main_files)) {
    x <- readRDS(main_files[i])
    for (j in seq_along(times_calc)) {
      point_est_list[[j]][i, ] <- calc_3(d = x, last_point = max(times_calc), t_cut = times_calc[j])
    }
    rm(x); gc()
  }
  saveRDS(point_est_list, file = paste0("point_est_list_S", scen_ind, ".rds"))
  rm(point_est_list)
  
  # build boot_list at selected time points
  # boot_files <- files[grepl("^bootstrap_results", files)]
  #t_boot <- sapply(c(0.5, 1), function(tt) times_calc[which.min(abs(times_calc - tt))])
  #boot_list <- build_boot_list(boot_files = boot_files, t_boot = t_boot, last_point = max(times_calc))
  # rm(boot_list); gc()

  all_results <- lapply(setNames(main_files, gsub("\\.rds$", "", main_files)), function(f) readRDS(f))
  beta_mat  <- do.call(rbind, lapply(all_results, `[[`, "beta_EM"))
  theta_mat <- do.call(rbind, lapply(all_results, `[[`, "thetas_EM"))
  beta_est = colMeans(beta_mat,  na.rm = TRUE)
  beta_se = apply(beta_mat, 2, sd, na.rm = TRUE)
  theta_est = colMeans(theta_mat, na.rm = TRUE)
  theta_se = apply(theta_mat, 2, sd, na.rm = TRUE)
  theta_MC_unq <- unique(do.call(rbind, lapply(all_results, `[[`, "theta.est")))
  theta_MC_mean = mean(do.call(rbind, lapply(all_results, `[[`, "theta.est")))
  theta_MC_se = sd(do.call(rbind, lapply(all_results, `[[`, "theta.est")))
  hist_num_EM_iterations = hist(do.call(rbind, lapply(all_results, `[[`, "num_EM_iterations")))
  # EM iterations and time, per sim iteration (then we can average)
  num_EM_iterations  <- do.call(rbind, lapply(all_results, `[[`, "num_EM_iterations"))
  em_timing  <- do.call(rbind, lapply(all_results, `[[`, "em_timing"))
  last_iter_elapsed_sec <- sapply(all_results, function(x) {
    tail(x$em_timing$iter_elapsed_sec, 1)
  })
  saveRDS(list(beta_est = beta_est, beta_se = beta_se, theta_est = theta_est, theta_se = theta_se
               ,theta_MC_unq = theta_MC_unq, theta_MC_mean = theta_MC_mean, theta_MC_se = theta_MC_se
               ,hist_num_EM_iterations = hist_num_EM_iterations
               ,num_EM_iterations=num_EM_iterations, last_iter_elapsed_sec=last_iter_elapsed_sec)
          ,file = paste0("summary_results_EM_S", scen_ind, ".rds"))
  
}

# ============================================================
# Bootstrap list + summary with empirical/boot SD + coverage
# ============================================================

build_boot_list <- function(boot_files, t_boot, last_point) {
  cols <- c("FICE_mean","FICE_sd","FICE_nok","SACE_mean","SACE_sd","SACE_nok",
            "AICE_mean","AICE_sd","AICE_nok","TE_mean","TE_sd","TE_nok",
            "pi_ios_mean","pi_ios_sd","pi_as_mean","pi_as_sd",
            "pi_ad_mean","pi_ad_sd","n_null","n_tot")
  boot_list <- lapply(t_boot, function(tc)
    matrix(NA_real_, nrow = length(boot_files), ncol = length(cols),
           dimnames = list(NULL, cols)))
  names(boot_list) <- t_boot
  
  for (i in seq_along(boot_files)) {
    br <- readRDS(boot_files[i])
    for (j in seq_along(t_boot)) {
      boot_list[[j]][i, ] <- calculations_per_file(boots      = br,
                                                   last_point = last_point,
                                                   t_cut      = t_boot[j])
    }
    rm(br)
  }
  boot_list
}

summarize_scen <- function(scen_name, point_est_list, boot_list, target_param_by_t, tt) {
  tc  <- as.numeric(names(point_est_list))
  tcb <- as.numeric(names(boot_list))
  j   <- which.min(abs(tc  - tt))
  jb  <- which.min(abs(tcb - tt))
  
  estimands <- colnames(point_est_list[[j]])
  z <- qnorm(0.975)
  
  est    <- colMeans(point_est_list[[j]], na.rm = TRUE)
  emp_sd <- apply(point_est_list[[j]], 2, sd, na.rm = TRUE)
  
  boot_SE <- sapply(estimands, function(e) {
    col <- paste0(e, "_sd")
    if (col %in% colnames(boot_list[[jb]])) mean(boot_list[[jb]][, col], na.rm = TRUE) else NA_real_
  })
  
  coverage <- sapply(estimands, function(e) {
    col <- paste0(e, "_sd")
    if (!(col %in% colnames(boot_list[[jb]]))) return(NA_real_)
    pt  <- point_est_list[[j]][, e]
    sdi <- boot_list[[jb]][, col]
    tru <- target_param_by_t[e, j]
    mean((pt - z * sdi <= tru) & (tru <= pt + z * sdi), na.rm = TRUE)
  })
  
  truth <- target_param_by_t[estimands, j]
  
  data.frame(
    scen     = scen_name,
    t        = round(tc[j], 3),
    estimand = estimands,
    truth    = round(truth,    4),
    est      = round(est,      4),
    bias     = round(est - truth, 4),
    bias_pct = round(100 * (est - truth) / truth, 2),
    emp_sd   = round(emp_sd,   4),
    boot_SE  = round(boot_SE,  4),
    se_bias  = round(boot_SE - emp_sd, 4),
    se_bias_pct = round(100 * (boot_SE - emp_sd) / emp_sd, 2),
    coverage = round(coverage, 4),
    row.names = NULL)
}

# ---- build and save boot_list for each scenario (run once) ----

###############################
'''rm(my.data_tweak, my.data_generation_tweak, dat_T.sim_tweak); gc()
for (scen_ind in scen_inds) {
  message("=== building boot_list S", scen_ind, " ===")
  
  setwd(file.path(base_path, paste0("S", scen_ind, "/tweak_n", n.sample)))
  files <- list.files(pattern = "\\.rds$")
  boot_files <- files[grepl("^bootstrap_results", files)]
  t_boot <- sapply(c(0.5, 1), function(tt) times_calc[which.min(abs(times_calc - tt))])
  boot_list <- build_boot_list(boot_files = boot_files,
                               t_boot     = t_boot,
                               last_point = max(times_calc))
  saveRDS(boot_list, file = paste0("boot_list_S", scen_ind, ".rds"))
  rm(boot_list)
}'''
###############################

###############################
library(future.apply)
build_boot_list_par <- function(boot_files, t_boot, last_point) {
  res <- future_lapply(boot_files, function(f) {
    br <- readRDS(f)
    out <- t(sapply(t_boot, function(tc)
      calculations_per_file(boots = br, last_point = last_point, t_cut = tc)))
    rm(br)
    out  # matrix: length(t_boot) x 20
  }, future.seed = NULL)
  
  cols <- colnames(res[[1]])
  boot_list <- lapply(seq_along(t_boot), function(j) {
    m <- do.call(rbind, lapply(res, function(x) x[j, ]))
    colnames(m) <- cols
    m
  })
  names(boot_list) <- t_boot
  boot_list
}

plan(multisession, workers = 4)  # adjust
for (scen_ind in scen_inds) {
  message("=== building boot_list S", scen_ind, " ===")
  dir_s <- file.path(base_path, paste0("S", scen_ind, "/tweak_n", n.sample))
  boot_files <- list.files(path = dir_s, pattern = "^bootstrap_results.*\\.rds$", full.names = TRUE)
  
  t_boot <- sapply(c(0.5, 1), function(tt) times_calc[which.min(abs(times_calc - tt))])
  
  boot_list <- build_boot_list_par(boot_files = boot_files,
                                   t_boot     = t_boot,
                                   last_point = max(times_calc))
  saveRDS(boot_list, file = file.path(dir_s, paste0("boot_list_S", scen_ind, ".rds")))
  rm(boot_list); gc()
}
plan(sequential)
###############################

# ---- summary loop (reads only) ---- start running fast
###############################
summary_tbl <- list()
summary_results_EM_lst <- list() 
for (scen_ind in scen_inds) {
  dir_s <- file.path(base_path, paste0("S", scen_ind, "/tweak_n", n.sample))
  point_est_list    <- readRDS(file.path(dir_s, paste0("point_est_list_S",    scen_ind, ".rds")))
  boot_list         <- readRDS(file.path(dir_s, paste0("boot_list_S",         scen_ind, ".rds")))
  target_param_by_t <- readRDS(file.path(dir_s, paste0("target_param_by_t_S", scen_ind, ".rds")))
  summary_results_EM <- readRDS(file.path(dir_s, paste0("summary_results_EM_S", scen_ind, ".rds")))

  summary_tbl[[paste0("S", scen_ind)]] <- do.call(rbind, lapply(c(0.5, 1), function(tt)
    summarize_scen(scen_name         = paste0("S", scen_ind),
                   point_est_list    = point_est_list,
                   boot_list         = boot_list,
                   target_param_by_t = target_param_by_t,
                   tt                = tt)))
  # EM info
  summary_results_EM_lst[[paste0("S", scen_ind)]] <- list(
    beta_est = summary_results_EM$beta_est
    ,beta_se  = summary_results_EM$beta_se
    ,theta_est = summary_results_EM$theta_est
    ,theta_se  = summary_results_EM$theta_se
    ,theta_MC_unq = summary_results_EM$theta_MC_unq
    ,theta_MC_mean = summary_results_EM$theta_MC_mean
    ,theta_MC_se   = summary_results_EM$theta_MC_se
    ,hist_num_EM_iterations = summary_results_EM$hist_num_EM_iterations
    ,num_EM_iterations = summary_results_EM$num_EM_iterations
    ,mean_num_EM_iterations = mean(summary_results_EM$num_EM_iterations)
    ,last_iter_elapsed_sec = summary_results_EM$last_iter_elapsed_sec
    ,mean_last_iter_elapsed_sec = mean(summary_results_EM$last_iter_elapsed_sec)
  )
}

summary_all <- do.call(rbind, summary_tbl)
rownames(summary_all) <- NULL
summary_all$group <- ifelse(summary_all$scen %in% paste0("S", 1:3), "S1-S3",
                            ifelse(summary_all$scen %in% paste0("S", 7:9), "S7-S9", NA))
target_estimands <- c("FICE","SACE","TE")

max_bias_tbl <- do.call(rbind, lapply(c("S1-S3","S7-S9"), function(g) {
  do.call(rbind, lapply(unique(summary_all$t), function(tt) {
    do.call(rbind, lapply(target_estimands, function(e) {
      sub <- summary_all[summary_all$group == g &
                           summary_all$estimand == e &
                           summary_all$t == tt, ]
      k_bias <- which.max(abs(sub$bias))
      k_se   <- which.max(abs(sub$se_bias))
      data.frame(
        group              = g,
        t                  = tt,
        estimand           = e,
        max_abs_bias       = sprintf("%.4f (%.2f%%)",
                                     abs(sub$bias[k_bias]),
                                     abs(sub$bias_pct[k_bias])),
        max_abs_bias_pct   = max(abs(sub$bias_pct),    na.rm = TRUE),
        max_abs_se_bias    = sprintf("%.4f (%.2f%%)",
                                     abs(sub$se_bias[k_se]),
                                     abs(sub$se_bias_pct[k_se])),
        max_abs_se_bias_pct= max(abs(sub$se_bias_pct), na.rm = TRUE))
    }))
  }))
}))

min_bias_tbl <- do.call(rbind, lapply(c("S1-S3","S7-S9"), function(g) {
  do.call(rbind, lapply(unique(summary_all$t), function(tt) {
    do.call(rbind, lapply(target_estimands, function(e) {
      sub <- summary_all[summary_all$group == g &
                           summary_all$estimand == e &
                           summary_all$t == tt, ]
      k_bias <- which.min(abs(sub$bias))
      k_se   <- which.min(abs(sub$se_bias))
      data.frame(
        group              = g,
        t                  = tt,
        estimand           = e,
        min_abs_bias       = sprintf("%.4f (%.2f%%)",
                                     abs(sub$bias[k_bias]),
                                     abs(sub$bias_pct[k_bias])),
        min_abs_bias_pct   = min(abs(sub$bias_pct),    na.rm = TRUE),
        min_abs_se_bias    = sprintf("%.4f (%.2f%%)",
                                     abs(sub$se_bias[k_se]),
                                     abs(sub$se_bias_pct[k_se])),
        min_abs_se_bias_pct= min(abs(sub$se_bias_pct), na.rm = TRUE))
    }))
  }))
}))

sapply(summary_results_EM_lst, `[[`, "mean_num_EM_iterations")
sapply(summary_results_EM_lst, `[[`, "mean_last_iter_elapsed_sec")

##############################################################
# ============================================================
# Tables for SM paper - FICE only, with bias/emp SD/boot SE/coverage, per scenario and time point
# ============================================================

# Filter to FICE only, select relevant columns
fice_tbl <- do.call(rbind, lapply(summary_tbl, function(df)
  df[df$estimand == "FICE", c("scen","t","estimand","truth","est","bias","emp_sd","boot_SE","coverage")]
))
rownames(fice_tbl) <- NULL

library(xtable)

fice_tbl <- do.call(rbind, lapply(summary_tbl, function(df)
  df[df$estimand == "FICE", c("scen","t","truth","est","bias","emp_sd","boot_SE","coverage")]
))
fice_tbl <- merge(x = fice_tbl, y = scen_meta, by.x = "scen", by.y = "name", sort = FALSE)
fice_tbl <- fice_tbl[order(fice_tbl$theta, fice_tbl$rho, fice_tbl$t), ]

body_rows  <- character(0)
prev_theta <- NULL
theta_n    <- tapply(fice_tbl$theta, fice_tbl$theta, length)

for (k in seq_len(nrow(fice_tbl))) {
  r  <- fice_tbl[k, ]
  th <- sprintf("\\multirow{%d}{*}{$%g$}", theta_n[as.character(r$theta)], r$theta)
  ro <- sprintf("\\multirow{2}{*}{$%g$}", r$rho)
  
  is_new_theta <- is.null(prev_theta) || r$theta != prev_theta
  is_new_rho   <- k == 1 || r$rho != fice_tbl$rho[k-1] || is_new_theta
  
  if (is_new_theta && !is.null(prev_theta)) body_rows <- c(body_rows, "\\midrule")
  prev_theta <- r$theta
  
  body_rows <- c(body_rows, sprintf(
    "%s & %s & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f \\\\",
    if (is_new_theta) th else "",
    if (is_new_rho)   ro else "",
    r$t, r$truth, r$est, r$bias, r$emp_sd, r$boot_SE, r$coverage
  ))
}

cat(paste(c(
  "\\begin{table}[ht]", "\\centering",
  "\\caption{FICE results. Est = mean estimate; Bias = Est $-$ Truth; Emp SD = empirical SD; Boot SE = mean bootstrap SE; Coverage = empirical 95\\% CI coverage.}",
  "\\label{tab:fice}",
  "\\begin{tabular}{ccrrrrrrr}", "\\toprule",
  "$\\theta$ & $\\rho$ & $t$ & Truth & Est & Bias & Emp SD & Boot SE & Coverage \\\\",
  "\\midrule", body_rows, "\\bottomrule", "\\end{tabular}", "\\end{table}"
), collapse = "\n"), "\n")
##############################################################

# ============================================================
# plots of curves with/wout sd ribbon, per scenario
# ============================================================

plot_curves_sd <- function(point_list, target_by_t, times_calc,
                           estimands = c("FICE","SACE","TE"),
                           title = NULL, error_bar = FALSE) {
  colors          <- c("FICE" = "green3", "SACE" = "steelblue", "TE" = "red")
  labels_estimand <- c("FICE" = "FICE",   "SACE" = "SACE",      "TE" = "Total effect")
  
  df <- do.call(rbind, lapply(seq_along(times_calc), function(j)
    do.call(rbind, lapply(estimands, function(est) {
      pt <- point_list[[j]][, est]
      data.frame(t = times_calc[j], estimand = est, truth = target_by_t[est, j],
                 mean_est = mean(pt, na.rm = TRUE), sd_est = sd(pt, na.rm = TRUE))
    }))
  ))
  df$estimand     <- factor(df$estimand,     levels = c("TE","FICE","SACE"))
  df$estimand_lab <- factor(labels_estimand[as.character(df$estimand)],
                            levels = labels_estimand[c("TE","FICE","SACE")])
  
  base_theme <- list(
    theme_bw(),
    scale_color_manual(values = colors, labels = labels_estimand),
    scale_fill_manual( values = colors, labels = labels_estimand, guide = "none"),
    scale_linetype_manual(values = c("Estimate" = "22", "Truth" = "solid"),
                          labels = c("Estimate" = "Estimate", "Truth" = "Estimand")),
    labs(x = "Time after baseline in years (t)", y = "Effect - Infection",
         title = title, color = NULL, linetype = NULL),
    theme(legend.position = "right", axis.title = element_text(size = 11),
          plot.title = element_text(hjust = 0.5), legend.key.width = unit(1, "cm"))
  )
  
  uncertainty_layer <- if (error_bar)
    geom_errorbar(data = df_long[df_long$type == "Estimate", ],
                  aes(ymin = value - sd_est, ymax = value + sd_est),
                  width = 0.02, linewidth = 0.4, alpha = 0.6, color = "black")
  else
    geom_ribbon(aes(ymin = mean_est - sd_est, ymax = mean_est + sd_est, fill = estimand),
                alpha = 0.25, color = NA)
  
  df_long <- rbind(
    data.frame(df[, c("t","estimand","sd_est")], value = df$mean_est, type = "Estimate"),
    data.frame(df[, c("t","estimand","sd_est")], value = df$truth,    type = "Truth")
  )
  df_long$estimand <- factor(df_long$estimand, levels = c("TE","FICE","SACE"))
  df_long$type     <- factor(df_long$type, levels = c("Estimate","Truth"))
  
  p_main <- ggplot(df_long, aes(x = t, y = value, color = estimand, linetype = type)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey40", linewidth = 0.5) +
    uncertainty_layer +
    geom_line(linewidth = 0.9) +
    base_theme
  
  p_with_sd <- ggplot(df, aes(x = t, color = estimand, fill = estimand)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey40", linewidth = 0.5) +
    uncertainty_layer +
    geom_line(aes(y = mean_est, linetype = "Estimate"), linewidth = 0.9) +
    geom_line(aes(y = truth,    linetype = "Truth"),    linewidth = 0.9) +
    facet_wrap(~ estimand_lab, nrow = 1, scales = "free_y") +
    base_theme
  
  list(main = p_main, with_sd = p_with_sd)
}

plot_grid_2x3 <- function(point_est_lists, target_by_t_lists, times_calc, scen_meta,
                          estimands = c("FICE","SACE","TE"), error_bar = FALSE, show_color_legend = TRUE){
  colors <- c("FICE" = ifelse(length(estimands) == 1, "black", "green3"),
              "SACE" = "steelblue", "TE" = "red")
  
  df <- do.call(rbind, lapply(seq_len(nrow(scen_meta)), function(k) {
    nm <- scen_meta$name[k]; pl <- point_est_lists[[nm]]; tb <- target_by_t_lists[[nm]]
    do.call(rbind, lapply(seq_along(times_calc), function(j)
      do.call(rbind, lapply(estimands, function(est) {
        pt <- pl[[j]][, est]
        data.frame(t = times_calc[j], estimand = est, truth = tb[est, j],
                   mean_est = mean(pt, na.rm = TRUE), sd_est = sd(pt, na.rm = TRUE),
                   rho = scen_meta$rho[k], theta = scen_meta$theta[k])
      }))
    ))
  }))
  df$estimand <- factor(df$estimand, levels = c("TE","FICE","SACE"))
  
  df_long <- rbind(
    data.frame(df[, c("t","estimand","sd_est","rho","theta")], value = df$mean_est, type = "Estimate"),
    data.frame(df[, c("t","estimand","sd_est","rho","theta")], value = df$truth,    type = "Truth")
  )
  df_long$type <- factor(df_long$type, levels = c("Estimate","Truth"))
  
  uncertainty_layer <- if (error_bar)
    geom_errorbar(data = df_long[df_long$type == "Estimate", ],
                  aes(ymin = value - sd_est, ymax = value + sd_est, group = estimand),
                  color = "black", width = 0.04, linewidth = 0.5, alpha = 0.5)
  else
    geom_ribbon(data = df_long[df_long$type == "Estimate", ],
                aes(ymin = value - sd_est, ymax = value + sd_est, fill = estimand),
                alpha = 0.25, color = NA)
  
  ggplot(df_long, aes(x = t, y = value, linetype = type)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey40", linewidth = 0.5) +
    uncertainty_layer +
    geom_line(aes(color = estimand), linewidth = 0.9) +
    facet_grid(rows = vars(theta), cols = vars(rho),
               labeller = label_bquote(rows = theta == .(theta), cols = rho == .(rho))) +
    scale_color_manual(values = colors,
                       labels = c(TE = "Total effect", FICE = "FICE", SACE = "SACE"),
                       guide = if (show_color_legend) "legend" else "none") +
    scale_fill_manual(values = colors, guide = "none") +
    scale_linetype_manual(values = c("Estimate" = "22", "Truth" = "solid"),
                          labels = c("Estimate" = "Mean estimate", "Truth" = "True")) +
    labs(x = "Time after baseline in years (t)",
         y = if (length(estimands) == 1) paste(estimands, "(t)") else "Causal effect on infection",
         color = NULL, linetype = NULL) +
    theme_bw() +
    theme(legend.key.width = unit(1, "cm"),
          legend.text  = element_text(size = 14),
          strip.text   = element_text(size = 16),
          axis.title   = element_text(size = 14),
          axis.text    = element_text(size = 12))
}

plots_by_scenario <- list()
for (scen_ind in scen_inds) {
  point_est_list = readRDS(file.path(base_path, paste0("S", scen_ind, "/tweak_n", n.sample), 
                                     paste0("point_est_list_S", scen_ind, ".rds")))
  target_param_by_t = readRDS(file.path(base_path, paste0("S", scen_ind, "/tweak_n", n.sample), 
                                        paste0("target_param_by_t_S", scen_ind, ".rds")))
  
  plots_by_scenario[[paste0("S", scen_ind)]] = 
    plot_curves_sd(point_list = point_est_list, target_by_t = target_param_by_t, 
                   times_calc = times_calc, title = NULL) # paste0("S", scen_ind)
}

point_est_lists   <- setNames(lapply(scen_inds, function(s)
  readRDS(file.path(base_path, paste0("S", s, "/tweak_n", n.sample),
                    paste0("point_est_list_S", s, ".rds")))), paste0("S", scen_inds))
target_by_t_lists <- setNames(lapply(scen_inds, function(s)
  readRDS(file.path(base_path, paste0("S", s, "/tweak_n", n.sample),
                    paste0("target_param_by_t_S", s, ".rds")))), paste0("S", scen_inds))

scen_meta <- data.frame(name = paste0("S", scen_inds), rho = c(0, 0.5, 1, 0, 0.5, 1), theta = c(1, 1, 1, 2, 2, 2))

# main text
plot_grid_2x3(point_est_lists = point_est_lists, target_by_t_lists = target_by_t_lists,
  times_calc = times_calc, scen_meta = scen_meta, estimands = "FICE", error_bar = TRUE, show_color_legend = FALSE)

# SM
plot_grid_2x3(point_est_lists, target_by_t_lists, times_calc, scen_meta)
