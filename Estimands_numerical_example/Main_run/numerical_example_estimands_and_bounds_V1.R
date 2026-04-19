# numerical example - bounds examination

#####################################################################################
# ------------------------------------------------------------------
# numerical example - bounds examination
library(dplyr)
library(ggplot2)
library(future.apply)
library(purrr)

setwd("~/R projects/AAA PhD/Causal-effects-on-non-terminal-event-time-with-application")
source("Estimands_numerical_example/Main_run/GetScenarioParams_paper.R")
source("Estimands_numerical_example/Main_run/run_functions.R")
source("Estimands_numerical_example/Data_generation/SimDataWeibFrail.R")
source("Estimands_numerical_example/Estimand_calculations/estimands_functions.R")
source("Estimands_numerical_example/Estimand_calculations/CalcTrueCausalParams.R")


# ------------------------------------------------------------------
# params
# ------------------------------------------------------------------
tau              <- 1
scen_seed        <- 50001; lett <- "a"

main_bool = TRUE # "main" or "SM"
scen_values_main      <- c(5, 7) # Main text (5: theta=1, 7: theta=2, 9 [original version]: theta=3) 
scen_values_SM = c(11, 12) #SM (11: theta=1, 12: theta=2, 13 [original version]: theta=3)
scen_values       <- if (main_bool) scen_values_main else scen_values_SM

X_numerical_ex   <- as.matrix(qnorm(seq(0.1, 0.9, 0.1)))
dim_x            <- ncol(X_numerical_ex)
coeff_A          <- c(0.2, 0.2)
times_calc       <- seq(0, 365, 5) / 365
criterion_vec    <- c("ios_ORP", "weak_ORP", "none")
rho_values       <- c(0, 0.5, 1) # c(0) # c(0, 0.5, 1)
ylims_diff       <- c(-0.15, 0.25)
ylims_ratio      <- c(0.5, 2.5) 
legend_text_size <- 11
n_sim            <- 30000000

set.seed(scen_seed + match(lett, letters))

# helper functions from functions_calculations_sim
#we can instead do source("Estimands_numerical_example/Simulation_studies/Summary_after_boot/functions_calculations_sim.R")
########################################################
calc_3 <- function(d, last_point=1, t_cut=1) {
  if (is.null(d)) return(c(FICE=NA_real_, SACE=NA_real_, AICE=NA_real_))
  
  T1.0 <- d$T1.0; T1.1 <- d$T1.1
  T2.0 <- d$T2.0; T2.1 <- d$T2.1
  
  ios <- ((T2.0==last_point) | (T1.0<T2.0)) & ((T2.1==last_point) | (T1.1<T2.1))
  ad  <- (T1.0<T2.0) & (T1.1<T2.1)
  as  <- (T2.0==last_point) & (T2.1==last_point)
  
  diff_mean <- function(mask){
    if (!any(mask, na.rm=TRUE)) return(NA_real_)
    mean(T1.1[mask] < t_cut, na.rm=TRUE) - mean(T1.0[mask] < t_cut, na.rm=TRUE)
  }
  
  c(FICE=diff_mean(ios), SACE=diff_mean(as), AICE=diff_mean(ad)
    , TE = mean(T1.1 < t_cut, na.rm=TRUE) - mean(T1.0 < t_cut, na.rm=TRUE)
    ,pi_ios=mean(ios, na.rm=T), pi_as=mean(as, na.rm=T), pi_ad=mean(ad, na.rm=T)
  )
}
calculations_per_file <- function(boots, last_point=1, t_cut=1) {
  #m <- t(vapply(boots, function(b) calc_3(b, last_point, t_cut), numeric(3)))
  m <- t(sapply(boots, function(b) calc_3(b, last_point, t_cut)))
  c(
    FICE_mean=mean(m[,"FICE"], na.rm=TRUE), FICE_sd=sd(m[,"FICE"], na.rm=TRUE), FICE_nok=sum(!is.na(m[,"FICE"]))
    ,SACE_mean=mean(m[,"SACE"], na.rm=TRUE), SACE_sd=sd(m[,"SACE"], na.rm=TRUE), SACE_nok=sum(!is.na(m[,"SACE"]))
    ,AICE_mean=mean(m[,"AICE"], na.rm=TRUE), AICE_sd=sd(m[,"AICE"], na.rm=TRUE), AICE_nok=sum(!is.na(m[,"AICE"]))
    ,TE_mean=mean(m[,"TE"], na.rm=TRUE), TE_sd=sd(m[,"TE"], na.rm=TRUE), TE_nok=sum(!is.na(m[,"TE"]))
    ,pi_ios_mean=mean(m[,"pi_ios"], na.rm=TRUE), pi_ios_sd=sd(m[,"pi_ios"], na.rm=TRUE)
    ,pi_as_mean =mean(m[,"pi_as"],  na.rm=TRUE), pi_as_sd=sd(m[,"pi_as"], na.rm=TRUE)
    ,pi_ad_mean =mean(m[,"pi_ad"],  na.rm=TRUE), pi_ad_sd=sd(m[,"pi_ad"], na.rm=TRUE)
    ,n_null=sum(vapply(boots, is.null, logical(1)))
    ,n_tot=length(boots)
  )
}
########################################################

# ------------------------------------------------------------------
# run one scenario x one criterion
# returns ###list(bounds, estimands)### with all 3 strata (the estimands) and both scales
# ------------------------------------------------------------------
run_one <- function(scen_ind, criterion, n_sim, times_calc, rho, save_data_files=FALSE) {
  dat <- Data_generation_func(
    scen = scen_values[scen_ind], rho = rho, tau = tau,
    n.sample = n_sim, X = X_numerical_ex, dim_x = dim_x,
    coeff_A = coeff_A, criterion = criterion
  )
  
  theta      <- dat$params$theta
  print(table(dat$my.data$X))
  d          <- dat$my.data
  
  if(save_data_files == TRUE){
    saveRDS(d, paste0("d_scen_", scen_ind, "_", scen_values[scen_ind], "_theta", theta, "_rho", rho, "_", criterion, ".RDS"))
  }
  
  T1.0 <- d$T1.0; T1.1 <- d$T1.1
  T2.0 <- d$T2.0; T2.1 <- d$T2.1
  rm(d, dat); gc()
  
  last_point <- max(times_calc)
  last_ind   <- length(times_calc)
  
  F0 <- sapply(times_calc, function(t) mean(T1.0 <= t & T1.0 < T2.0))
  F1 <- sapply(times_calc, function(t) mean(T1.1 <= t & T1.1 < T2.1))
  
  ad0  <- T1.0 < T2.0
  ios0 <- (T2.0 == last_point) | (T1.0 < T2.0)
  ios1 <- (T2.1 == last_point) | (T1.1 < T2.1)
  ios  <- ios0 & ios1
  ad   <- (T1.0 < T2.0) & (T1.1 < T2.1)
  as  <- (T2.0 == last_point) & (T2.1 == last_point)
  
  p0 <- mean(ios0); p1 <- mean(ios1)
  T2_inf_dead_0 <- mean(T2.0[ad0] < 1)
  
  # -- bounds diff scale --
  L_pi_w <- max(0, p0 + p1 - 1);  U_pi_w <- min(p0, p1)
  if (L_pi_w > 0) {
    # L_wout_diff <- (pmax(0, F1 + p0 - 1) / U_pi_w) - (pmin(F0, p1) / L_pi_w)
    # U_wout_diff <- (pmin(F1, p0)          / L_pi_w) - (pmax(0, F0 + p1 - 1) / U_pi_w)
    #TODO
    l_wout_diff = pmax(0, F1 + p0 - 1) - pmin(F0, p1)
    u_wout_diff = pmin(F1, p0) - pmax(0, F0 + p1 - 1)
    L_wout_diff <- l_wout_diff / ifelse(l_wout_diff >= 0, U_pi_w, L_pi_w)
    U_wout_diff <- u_wout_diff / ifelse(u_wout_diff >= 0, L_pi_w, U_pi_w)
    
    L_wout_ratio <- pmax(0, F1 + p0 - 1) / pmin(F0, p1)
    U_wout_ratio <- pmin(F1, p0) / pmax(0, F0 + p1 - 1)
  } else {
    L_wout_diff  <- rep(NA_real_, length(times_calc))
    U_wout_diff  <- rep(NA_real_, length(times_calc))
    L_wout_ratio <- rep(NA_real_, length(times_calc))
    U_wout_ratio <- rep(NA_real_, length(times_calc))
  }
  
  L_pi_wk <- max(F0[last_ind], p0 + p1 - 1)
  U_pi_wk <- pmin(p0, p1 + F0[last_ind] * T2_inf_dead_0)
  if (L_pi_wk > 0) {
    # L_weak_diff  <- (pmax(0, F1 + p0 - 1) / U_pi_wk) - pmin(1, F0 / L_pi_wk)
    # U_weak_diff  <- (pmin(F1, p0)          / L_pi_wk) - (F0 / U_pi_wk)
    #TODO
    l_weak_diff = pmax(0, F1 + p0 - 1) - F0
    u_weak_diff = pmin(F1, p0) - F0
    L_weak_diff <- l_weak_diff / ifelse(l_weak_diff >= 0, U_pi_wk, L_pi_wk)
    U_weak_diff <- u_weak_diff / ifelse(u_weak_diff >= 0, L_pi_wk, U_pi_wk)
    
    L_weak_ratio <- pmax(0, F1 + p0 - 1) / F0
    U_weak_ratio <- pmin(F1, p0) / F0
  } else {
    L_weak_diff  <- rep(NA_real_, length(times_calc))
    U_weak_diff  <- rep(NA_real_, length(times_calc))
    L_weak_ratio <- rep(NA_real_, length(times_calc))
    U_weak_ratio <- rep(NA_real_, length(times_calc))
  }
  
  L_ios_diff  <- pmax(0, 1 - (1 - F1) / p0) - F0 / p0
  U_ios_diff  <- pmin(1, F1 / p0)            - F0 / p0
  L_ios_ratio <- pmax(0, 1 - (1 - F1) / p0) / (F0 / p0)
  U_ios_ratio <- pmin(1, F1 / p0)            / (F0 / p0)
  
  # -- true estimands: FICE (ios), AICE (ad), SACE (as), both scales --
  calc_effects <- function(stratum) {
    F1_1 <- sapply(times_calc, function(t) mean(T1.1[stratum] <= t))
    F1_0 <- sapply(times_calc, function(t) mean(T1.0[stratum] <= t))
    data.frame(time = times_calc, diff = F1_1 - F1_0, ratio = F1_1 / F1_0)
  }
  
  eff <- rbind(
    cbind(estimand = "FICE", calc_effects(ios)),
    cbind(estimand = "AICE", calc_effects(ad)),
    cbind(estimand = "SACE", calc_effects(as))
  )
  eff$theta <- theta; eff$rho <- rho; eff$dgm <- criterion
  
  true_diff  <- eff$diff[ eff$estimand == "FICE"]
  true_ratio <- eff$ratio[eff$estimand == "FICE"]
  
  bounds_df = NULL
  if(rho == 0){
    make_bounds_rows <- function(assumption, L_diff, U_diff, L_ratio, U_ratio) {
      rbind(
        data.frame(assumption=assumption, bound="Lower", scale="diff",  time=times_calc, value=L_diff),
        data.frame(assumption=assumption, bound="Upper", scale="diff",  time=times_calc, value=U_diff),
        data.frame(assumption=assumption, bound="Lower", scale="ratio", time=times_calc, value=L_ratio),
        data.frame(assumption=assumption, bound="Upper", scale="ratio", time=times_calc, value=U_ratio)
      )
    }
    
    bounds_df <- rbind(
      make_bounds_rows("Without ORP", L_wout_diff, U_wout_diff, L_wout_ratio, U_wout_ratio),
      make_bounds_rows("Weak-ORP",    L_weak_diff, U_weak_diff, L_weak_ratio, U_weak_ratio),
      make_bounds_rows("ios-ORP",     L_ios_diff,  U_ios_diff,  L_ios_ratio,  U_ios_ratio)
    )
    bounds_df$theta <- theta
    bounds_df$rho   <- rho
    bounds_df$dgm   <- criterion
    bounds_df$true_diff  <- true_diff[ match(bounds_df$time, times_calc)]
    bounds_df$true_ratio <- true_ratio[match(bounds_df$time, times_calc)]
  }
  
  list(bounds = bounds_df
       ,estimands = eff
       #,my.data = d
       )
}


combos <- expand.grid(rho = rho_values, scen_ind = seq_along(scen_values), stringsAsFactors = FALSE)
combos$scen <- scen_values[combos$scen_ind]

# ------------------------------------------------------------------
# save datasets for all scenarios under no ORP assumptions (Section 3.4)
# ------------------------------------------------------------------
########################################################################
setwd(paste0("~/R projects/AAA PhD/Causal-effects-on-non-terminal-event-time-with-application/Estimands_numerical_example/Main_run/Bounds/",
      ifelse(main_bool, "Main_text", "SM"), "/DFs"))
plan(multisession, workers = 2)
raw_save_datasets <- future_lapply(
  seq_len(nrow(combos)),
  function(i)
    run_one(
      scen_ind   = combos$scen_ind[i],
      criterion  = "none",
      n_sim      = n_sim,
      times_calc = times_calc,
      rho        = combos$rho[i],
      save_data_files=TRUE
    ),
  future.seed = 101
)
plan(sequential)

# summarizing results for all scenarios
#path = paste0("~/R projects/AAA PhD/Causal-effects-on-non-terminal-event-time-with-application/Estimands_numerical_example/Main_run/Bounds/",
#             ifelse(main_bool, "Main_text", "SM"), "/DFs")
files <- list.files(getwd(), full.names = TRUE, pattern = "\\.RDS$|\\.csv$|\\.RData$")
results <- map_dfr(files, function(f) {
  df <- readRDS(f)  
  name <- tools::file_path_sans_ext(basename(f))
  n <- nrow(df)
  tibble(
    dataset  = name,
    ai       = mean(df$ai),
    as       = mean(df$as),
    ios      = mean(df$ios),
    two_pt   = mean(df$two_pt),
    n        = n
  )
})
print(results, n = Inf)
########################################################################

########################################################################
# ------------------------------------------------------------------
# build master_df (bounds) and estimands_df
# ------------------------------------------------------------------

# we can remove  "none", and add raw_save_datasets to raw
plan(multisession, workers = 1)
raw <- future_lapply(
  seq_len(nrow(combos)),
  function(i) lapply(criterion_vec, function(cr)
    run_one(scen_ind   = combos$scen_ind[i],
            criterion  = cr,
            n_sim      = n_sim,
            times_calc = times_calc,
            rho        = combos$rho[i])),
  future.seed = 101
)
plan(sequential)

master_df    <- do.call(rbind, lapply(unlist(raw, recursive = FALSE), `[[`, "bounds"))
estimands_df <- do.call(rbind, lapply(unlist(raw, recursive = FALSE), `[[`, "estimands"))

estimands_df$theta_label <- factor(paste0("theta==", estimands_df$theta),
                                   levels = paste0("theta==", sort(unique(estimands_df$theta))))
estimands_df$rho_label   <- factor(paste0("rho==", estimands_df$rho),
                                   levels = paste0("rho==", sort(unique(estimands_df$rho))))
saveRDS(estimands_df, "estimands_df.RDS")

# factors
master_df$assumption  <- factor(master_df$assumption,
                                levels = c("ios-ORP", "Weak-ORP", "Without ORP"))
dgm_labels            <- c(ios_ORP = "ios-ORP", weak_ORP = "Weak-ORP", none = "Without ORP")
master_df$dgm_label   <- factor(dgm_labels[master_df$dgm],
                                levels = c("ios-ORP", "Weak-ORP", "Without ORP"))
master_df$theta_label <- factor(paste0("theta==", master_df$theta),
                                levels = paste0("theta==", sort(unique(master_df$theta))))
master_df$line_id     <- paste(master_df$dgm, master_df$assumption,
                               master_df$bound, master_df$theta, master_df$rho, sep = "_")
saveRDS(master_df, "master_df.RDS")
########################################################################

# ------------------------------------------------------------------
# plot: bounds grid (rows = theta, cols = DGM)
# ------------------------------------------------------------------
plot_grid_bounds <- function(df, scale = "diff", ylims = NULL, legend_text_size = 11) {
  df <- df[df$scale == scale, ]
  # Before the ggplot call, shift red lower bounds slightly
  df$value[df$assumption == "Without ORP" & grepl("Lower", df$bound, ignore.case = TRUE)] <-
    df$value[df$assumption == "Without ORP" & grepl("Lower", df$bound, ignore.case = TRUE)] - 0.005
  
  true_col <- if (scale == "diff") "true_diff" else "true_ratio"
  df_true          <- unique(df[, c("theta_label", "dgm_label", "time", true_col)])
  names(df_true)[names(df_true) == true_col] <- "true_val"
  df_true$group_id <- paste(df_true$theta_label, df_true$dgm_label)
  clrs <- c("ios-ORP" = "gray", "Weak-ORP" = "deepskyblue1", "Without ORP" = "red")
  ref  <- if (scale == "diff") 0 else 1
  
  p <- ggplot(df, aes(x = time, y = value, color = assumption, group = line_id)) +
    geom_line(linewidth = 1.5) +
    #geom_point(size = 2) +
    geom_point(size = 0, alpha = 0) + # for dots in the legend
    geom_line(data = df_true,
              aes(x = time, y = true_val, group = group_id),
              color = "black", linewidth = 1.5, linetype = "dashed", inherit.aes = FALSE) +
    geom_hline(yintercept = ref, linetype = "dashed", color = "black", linewidth = 0.7) +
    scale_color_manual(name = "Assumption", values = clrs) +
    facet_grid(theta_label ~ dgm_label,
               labeller = labeller(theta_label = label_parsed)
               #,scales = "free_y"
               ) +
    #scale_y_continuous(breaks = seq(-1, 1, by = 0.05)) + 
    xlab("Time after baseline in years (t)") +
    ylab(if (scale == "diff") "FICE(t) bounds (risk difference)" else "Bounds (risk ratio)") +
    #guides(color = guide_legend(override.aes = list(size = 6))) +
    guides(color = guide_legend(override.aes = list(size = 4, shape = 16, alpha = 1, linetype = 0))) +
    theme_bw() +
    theme(
      strip.text      = element_text(size = 13, face = "bold"),
      axis.title      = element_text(size = 16),
      axis.text       = element_text(size = 11),
      panel.spacing   = unit(0.8, "lines"),
      legend.position = "right",
      legend.title    = element_text(size = 16, face = "bold"),
      legend.text     = element_text(size=rel(1.75)) # size = legend_text_size + 1
    )
  
  if (!is.null(ylims)) p <- p + coord_cartesian(ylim = ylims)
  p
}

plot_estimands <- function(dgm_name, scale = "diff", ylims = NULL, legend_text_size = 11) {
  sub             <- estimands_df[estimands_df$dgm == dgm_name, ]
  sub$theta_label <- factor(paste0("theta==", sub$theta),
                            levels = paste0("theta==", sort(unique(sub$theta))))
  sub$rho_label   <- factor(paste0("rho==", sub$rho),
                            levels = paste0("rho==", sort(unique(sub$rho))))
  sub$estimand    <- factor(sub$estimand, levels = c("FICE", "SACE", "AICE"))
  sub$group_id    <- paste(sub$estimand, sub$theta, sub$rho)
  yvar <- if (scale == "diff") "diff" else "ratio"
  ref  <- if (scale == "diff") 0 else 1
  
  p <- ggplot(sub, aes_string("time", yvar, color = "estimand", group = "group_id")) +
    geom_line(linewidth = 1.5) +
    geom_hline(yintercept = ref, linetype = "dashed", color = "black", linewidth = 0.7) +
    scale_color_manual(name = "Estimand",
                       values = c("FICE" = "green3", "SACE" = "steelblue", "AICE" = "red")) +
    facet_grid(theta_label ~ rho_label,
               labeller = labeller(theta_label = label_parsed, rho_label = label_parsed)) +
    xlab("Time after baseline in years (t)") +
    ylab(if (scale == "diff") "Causal effect on infection" else "Risk ratio") +
    theme_bw() +
    theme(
      strip.text      = element_text(size = 13, face = "bold"),
      axis.title      = element_text(size = 16),
      axis.text       = element_text(size = 16),
      axis.text.x     = element_text(size = 10),
      axis.text.y     = element_text(size = 16),
      panel.spacing   = unit(0.8, "lines"),
      legend.position = "right",
      legend.title = element_text(size = legend_text_size + 4, face = "bold"),  # add face="bold"
      legend.text     = element_text(size = legend_text_size + 4)
    )
  
  if (!is.null(ylims)) p <- p + coord_cartesian(ylim = ylims)
  p
}


master_df_bounds = master_df %>% filter(rho==0)
fig_ratio_bounds <- 
  plot_grid_bounds(master_df_bounds, scale = "ratio", ylims = ylims_ratio, legend_text_size = legend_text_size)

fig_diff_bounds_SM  <- 
  plot_grid_bounds(master_df_bounds, scale = "diff",  ylims = ylims_diff,  legend_text_size = legend_text_size)
fig_diff_bounds_main_text <- plot_grid_bounds(df = master_df_bounds %>% filter(dgm == "ios_ORP" & theta == 1),
                                  scale = "diff", ylims = c(-0.15, 0.25), legend_text_size = legend_text_size)
fig_diff_bounds_main_text + facet_null() + theme(strip.text = element_blank())


fig_est_ios_ORP  <- plot_estimands("ios_ORP")
fig_est_weak_ORP <- plot_estimands("weak_ORP")
fig_est_wout_ORP <- plot_estimands("none")


# check DGMs
out_b <- do.call(rbind, lapply(times_calc, function(tc) {
  sub_b <- master_df[abs(master_df$time - tc) < 1e-6 & master_df$scale == "diff" & master_df$rho == 0,
                     c("time","dgm","theta","assumption","bound","value")]
  sub_b$time <- round(tc, 4)
  sub_b
}))
