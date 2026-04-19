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