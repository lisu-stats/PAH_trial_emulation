tmle_msm_fit_once <- function(data_sub, regimes_sub) {
  ltmleMSM(
    data                     = data_sub,
    Anodes                   = Anodes,
    Cnodes                   = Cnodes,
    Lnodes                   = Lnodes_filtered,
    Ynodes                   = Ynodes,
    final.Ynodes             = final.Ynodes,
    survivalOutcome          = TRUE,
    regimes                  = regimes_sub,
    summary.measures         = summary.measures,
    working.msm              = working.msm,
    Qform                    = Qform,
    gform                    = gform,
    deterministic.g.function = MaintainTreatment,
    SL.library               = mylibrary2,
    gbounds                  = c(0.01,1),
    variance.method          = "ic",
    iptw.only                = FALSE
  )
}

gj_worker <- function(g, data, regimes, groups, coef_names) {
  library(ltmle)
  
  keep <- which(groups != g)
  data_g    <- data[keep, , drop = FALSE]
  regimes_g <- regimes[keep, , , drop = FALSE]
  
  set.seed(g)
  fit_g <- try(tmle_msm_fit_once(data_g, regimes_g), silent = TRUE)
  if (inherits(fit_g, "try-error") || is.null(fit_g$beta))
    return(rep(NA_real_, length(coef_names)))
  
  bi <- fit_g$beta
  names(bi)<-coef_names
  bi
  print(bi)
}


# G=20
# n=743
# groups <- sample(rep(1:G, length.out = n))
# 
# 
# gj_worker(1, wide_df_ltmle, regimes, groups, coef.names )

grouped_jackknife_ltmleMSM <- function(data, regimes, G = NULL, ncores = 4,
                                       coef_names = NULL) {
  
  n <- nrow(data)
  stopifnot(dim(regimes)[1] == n)
  
  # Full-sample MSM fit (deterministic if SL library is deterministic + seed set)
  set.seed(12345)
  fit_full <- tmle_msm_fit_once(data, regimes)
  beta_full <- fit_full$beta
  names(beta_full)<-coef_names
  print(beta_full)
  
  if (is.null(coef_names)) coef_names <- names(beta_full)
  p <- length(coef_names)
  
  # Choose number of groups (default: ≈ sqrt(n), min 10 max 25)
  if (is.null(G)) G <- max(10, min(25, round(sqrt(n))))
  
  # random partition into G groups (set seed for reproducibility)
  set.seed(12345)
  groups <- sample(rep(1:G, length.out = n))
  
 

  
  # cluster for Windows
  cl <- parallel::makeCluster(ncores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  # export needed objects
  parallel::clusterExport(
    cl,
    varlist = c(
      "gj_worker", "tmle_msm_fit_once",
      "Anodes","Cnodes","Lnodes_filtered","Ynodes","final.Ynodes",
      "summary.measures","working.msm","Qform","gform",
      "mylibrary2", "coef.names"
    ),
    envir = globalenv()
  )
  
  
  
  # load package once on each worker
  parallel::clusterEvalQ(cl, { library(ltmle); NULL })

  # run grouped jackknife
  gj_list <- parallel::parLapply(
    cl,
    X = 1:G,
    fun = gj_worker,
    data = data,
    regimes = regimes,
    groups = groups,
    coef_names = coef_names
  )
  
  # bind into matrix (G × p)
  B <- do.call(rbind, gj_list)
  B <- as.data.frame(B)
  names(B)<-coef_names
  
  # delete-a-group jackknife variance:
  # Var = (G - 1)/G * sum( (θ_g - mean(θ_g))^2 )
  B_mean <- colMeans(B, na.rm = TRUE)
  gj_var <- (G - 1) / G * colSums((t(t(B) - B_mean))^2, na.rm = TRUE)
  gj_se  <- sqrt(gj_var)
  
  list(
    beta_full = beta_full[coef_names],
    beta_groups = B,
    gj_mean  = B_mean,
    gj_se    = gj_se,
    gj_var   = gj_var,
    groups   = groups
  )
}



#------------------(!!!!!run group jackknife!!!!!!!)-----------------------------------------

#--------------------------------- (1) main model -------------------------------------------------
horizons = c(5,10,15,20)


###create coef.names 
coef.names<-NULL
for (k in horizons) {
  coef.names       <- c(coef.names, paste0("year", k))
}
  coef.names       <- c(coef.names, "assigned_treatment")


  working.msm<-"Y ~ -1 + year5  + year10 + year15  + year20 + assigned_treatment"
  
  

GJK.main <- grouped_jackknife_ltmleMSM(
  data    = wide_df_ltmle,
  regimes = regimes,
  G       = 20,        # recommended
  ncores  = 20,
  coef_names = coef.names
)

print(GJK.main$beta_full)
print(GJK.main$gj_se)

save(GJK.main, file='./results/Group20/group_jackknife_msm_main.rdata')






get.summary.main<-function(GJK.main)
{  
# logistic expit function
expit <- function(x) 1 / (1 + exp(-x))


#############################################################
## 1. Compute baseline risks (never-treated) and risk differences
##    for the MSM with no intercept and treatment × year interaction
#############################################################

# INPUT:
#   beta      = named vector of MSM coefficients (full or group-delete)
#   horizons  = c(5,10,15,20)
#
# OUTPUT:
#   list with:
#      $risk_never  = baseline risks at each horizon
#      $RD          = risk differences (treated - never)

compute_risk_and_RD <- function(beta, horizons = c(5,10,15,20)) {
  
  risk_never <- numeric(length(horizons))
  RD         <- numeric(length(horizons))
  
  for (j in seq_along(horizons)) {
    
    k <- horizons[j]
    
    # names of coefficients in MSM
    b0_name <- paste0("year", k)
    b1_name <- "assigned_treatment"
    
    # logit for never-treated
    b0 <- beta[b0_name]
    
    # logit for treated
    b1 <- beta[b1_name]
    
    # risks
    p_never <- expit(b0)
    p_init  <- expit(b0 + b1)
    
    # store
    risk_never[j] <- p_never
    RD[j]         <- p_init - p_never
  }
  
  names(risk_never) <- paste0("risk_never_year", horizons)
  names(RD)         <- paste0("RD_year",          horizons)
  
  list(risk_never = risk_never,
       RD         = RD)
}


#############################################################
## 2. Compute full-sample baseline risks and RDs
#############################################################

horizons <- c(5,10,15,20)

full_vals <- compute_risk_and_RD(GJK.main$beta_full, horizons)

risk_never_full <- full_vals$risk_never
RD_full         <- full_vals$RD


#############################################################
## 3. Compute grouped jackknife estimates for each group
#############################################################

G <- nrow(GJK.main$beta_groups)

risk_never_groups <- matrix(NA, nrow = length(horizons), ncol = G,
                            dimnames = list(paste0("year", horizons),
                                            paste0("g", 1:G)))

RD_groups <- matrix(NA, nrow = length(horizons), ncol = G,
                    dimnames = list(paste0("year", horizons),
                                    paste0("g", 1:G)))


for (g in 1:G) {
  
  tmp <- compute_risk_and_RD(GJK.main$beta_groups[g,], horizons)
  
  risk_never_groups[, g] <- unlist(tmp$risk_never)
  RD_groups[, g]         <- unlist(tmp$RD)
}


#############################################################
## 4. Compute grouped jackknife variance and SE
#############################################################

# means across delete-g groups
risk_never_mean <- rowMeans(risk_never_groups, na.rm = TRUE)
RD_mean         <- rowMeans(RD_groups, na.rm = TRUE)

# grouped jackknife variance formula:
#   Var = (G - 1)/G * sum_g (theta_g - theta_mean)^2
risk_never_var <- (G - 1)/G * rowSums((risk_never_groups - risk_never_mean)^2,
                                      na.rm = TRUE)
RD_var         <- (G - 1)/G * rowSums((RD_groups - RD_mean)^2,
                                      na.rm = TRUE)

risk_never_se <- sqrt(risk_never_var)
RD_se         <- sqrt(RD_var)


#############################################################
## 5. 95% confidence intervals
#############################################################

z <- qt(0.975, df=G-1)

risk_never_lo <- risk_never_full - z * risk_never_se
risk_never_hi <- risk_never_full + z * risk_never_se

RD_lo <- RD_full - z * RD_se
RD_hi <- RD_full + z * RD_se


#############################################################
## 6. Final summary table
#############################################################

RD_summary <- data.frame(
  horizon       = horizons,
  
  # baseline risks
  risk_never      = risk_never_full,
  risk_never_SE   = risk_never_se,
  risk_never_LO   = risk_never_lo,
  risk_never_HI   = risk_never_hi,
  
  # risk differences
  RD              = RD_full,
  RD_SE           = RD_se,
  RD_LO           = RD_lo,
  RD_HI           = RD_hi
)

print(round(RD_summary, 6))

return(RD_summary)

}

get.surv.summary.main <- function(GJK.main) {
  
  # logistic expit
  expit <- function(x) 1 / (1 + exp(-x))
  
  #############################################################
  ## 1. Compute risks and survival probabilities from MSM
  #############################################################
  
  compute_survival <- function(beta, horizons = c(5,10,15,20)) {
    
    risk_never <- numeric(length(horizons))
    risk_trt   <- numeric(length(horizons))
    
    surv_never <- numeric(length(horizons))
    surv_trt   <- numeric(length(horizons))
    
    for (j in seq_along(horizons)) {
      
      k <- horizons[j]
      
      b0 <- beta[paste0("year", k)]
      b1 <- beta["assigned_treatment"]
      
      # cumulative risks
      p_never <- expit(b0)
      p_trt   <- expit(b0 + b1)
      
      # survival probabilities
      risk_never[j] <- p_never
      risk_trt[j]   <- p_trt
      surv_never[j] <- 1 - p_never
      surv_trt[j]   <- 1 - p_trt
    }
    
    list(
      risk_never = risk_never,
      risk_trt   = risk_trt,
      surv_never = surv_never,
      surv_trt   = surv_trt
    )
  }
  
  
  #############################################################
  ## 2. Full-sample estimates
  #############################################################
  
  horizons <- c(5,10,15,20)
  
  full_vals <- compute_survival(GJK.main$beta_full, horizons)
  
  surv_never_full <- full_vals$surv_never
  surv_trt_full   <- full_vals$surv_trt
  
  
  #############################################################
  ## 3. Grouped jackknife estimates
  #############################################################
  
  G <- nrow(GJK.main$beta_groups)
  
  surv_never_groups <- matrix(NA, nrow = length(horizons), ncol = G)
  surv_trt_groups   <- matrix(NA, nrow = length(horizons), ncol = G)
  
  for (g in 1:G) {
    
    tmp <- compute_survival(GJK.main$beta_groups[g,], horizons)
    
    surv_never_groups[, g] <- unlist(tmp$surv_never)
    surv_trt_groups[, g]   <- unlist(tmp$surv_trt)
  }
  
  
  #############################################################
  ## 4. Grouped jackknife variance and SE
  #############################################################
  
  surv_never_mean <- rowMeans(surv_never_groups, na.rm = TRUE)
  surv_trt_mean   <- rowMeans(surv_trt_groups, na.rm = TRUE)
  
  surv_never_var <- (G - 1) / G *
    rowSums((surv_never_groups - surv_never_mean)^2, na.rm = TRUE)
  
  surv_trt_var <- (G - 1) / G *
    rowSums((surv_trt_groups - surv_trt_mean)^2, na.rm = TRUE)
  
  surv_never_se <- sqrt(surv_never_var)
  surv_trt_se   <- sqrt(surv_trt_var)
  
  
  #############################################################
  ## 5. 95% confidence intervals
  #############################################################
  
  z <- qt(0.975, df=G-1)
  
  surv_never_lo <- surv_never_full - z * surv_never_se
  surv_never_hi <- surv_never_full + z * surv_never_se
  
  surv_trt_lo <- surv_trt_full - z * surv_trt_se
  surv_trt_hi <- surv_trt_full + z * surv_trt_se
  
  
  #############################################################
  ## 6. Final summary table
  #############################################################
  
  summary_tab <- data.frame(
    horizon = horizons,
    
    # never treated
    surv_never    = surv_never_full,
    surv_never_SE = surv_never_se,
    surv_never_LO = surv_never_lo,
    surv_never_HI = surv_never_hi,
    
    # treated
    surv_trt    = surv_trt_full,
    surv_trt_SE = surv_trt_se,
    surv_trt_LO = surv_trt_lo,
    surv_trt_HI = surv_trt_hi
  )
  
  print(round(summary_tab, 6))
  
  return(summary_tab)
}


print('## (1) main model: working.msm<-"Y ~ -1 + year5  + year10 + year15  + year20 + assigned_treatment"')


print('## 20 jackknife groups; risk for never initiate regime and risk differences by years')

load(file='./results/Group20/group_jackknife_msm_main.rdata')

RD_summary_main_g20<-get.summary.main(GJK.main)

print('## 20 jackknife groups; survivial probabilities for never initiate regime and initiate regime by years')

Surv_summary_main_g20<-get.surv.summary.main(GJK.main)





#---------------------------- (2) treatment by year -----------------------------------------------------------

###create coef.names 
coef.names<-NULL
for (k in horizons) {
  coef.names       <- c(coef.names, paste0("year", k))
}
for (k in horizons) {
  coef.names       <- c(coef.names, paste0("assigned_treatment:year", k))
}

working.msm<-"Y ~ -1 + year5 + assigned_treatment:year5 + year10 + assigned_treatment:year10 + year15 + assigned_treatment:year15 + year20 + assigned_treatment:year20"


 GJK <- grouped_jackknife_ltmleMSM(
   data    = wide_df_ltmle,
   regimes = regimes,
   G       = 20,        # recommended
   ncores  = 20,
   coef_names = coef.names
 )

 print(GJK$beta_full)
 print(GJK$gj_se)


 save(GJK, file='./results/Group20/group_jackknife_msm_year.rdata')







get.summary.trt.year<-function(GJK)
{  
  
#---------------------summarize results, year*treatment effects------------
#############################################################
## 0. Helpers
#############################################################

# logistic expit function
expit <- function(x) 1 / (1 + exp(-x))


#############################################################
## 1. Compute baseline risks (never-treated) and risk differences
##    for the MSM with no intercept and treatment × year interaction
#############################################################

# INPUT:
#   beta      = named vector of MSM coefficients (full or group-delete)
#   horizons  = c(5,10,15,20)
#
# OUTPUT:
#   list with:
#      $risk_never  = baseline risks at each horizon
#      $RD          = risk differences (treated - never)

compute_risk_and_RD <- function(beta, horizons = c(5,10,15,20)) {
  
  risk_never <- numeric(length(horizons))
  RD         <- numeric(length(horizons))
  
  for (j in seq_along(horizons)) {
    
    k <- horizons[j]
    
    # names of coefficients in MSM
    b0_name <- paste0("year", k)
    b1_name <- paste0("assigned_treatment:year", k)
    
    # logit for never-treated
    b0 <- beta[b0_name]
    
    # logit for treated
    b1 <- beta[b1_name]
    
    # risks
    p_never <- expit(b0)
    p_init  <- expit(b0 + b1)
    
    # store
    risk_never[j] <- p_never
    RD[j]         <- p_init - p_never
  }
  
  names(risk_never) <- paste0("risk_never_year", horizons)
  names(RD)         <- paste0("RD_year",          horizons)
  
  list(risk_never = risk_never,
       RD         = RD)
}


#############################################################
## 2. Compute full-sample baseline risks and RDs
#############################################################

horizons <- c(5,10,15,20)

full_vals <- compute_risk_and_RD(GJK$beta_full, horizons)

risk_never_full <- full_vals$risk_never
RD_full         <- full_vals$RD


#############################################################
## 3. Compute grouped jackknife estimates for each group
#############################################################

G <- nrow(GJK$beta_groups)

risk_never_groups <- matrix(NA, nrow = length(horizons), ncol = G,
                            dimnames = list(paste0("year", horizons),
                                            paste0("g", 1:G)))

RD_groups <- matrix(NA, nrow = length(horizons), ncol = G,
                    dimnames = list(paste0("year", horizons),
                                    paste0("g", 1:G)))


for (g in 1:G) {
  
  tmp <- compute_risk_and_RD(GJK$beta_groups[g,], horizons)
  
  risk_never_groups[, g] <- unlist(tmp$risk_never)
  RD_groups[, g]         <- unlist(tmp$RD)
}


#############################################################
## 4. Compute grouped jackknife variance and SE
#############################################################

# means across delete-g groups
risk_never_mean <- rowMeans(risk_never_groups, na.rm = TRUE)
RD_mean         <- rowMeans(RD_groups, na.rm = TRUE)

# grouped jackknife variance formula:
#   Var = (G - 1)/G * sum_g (theta_g - theta_mean)^2
risk_never_var <- (G - 1)/G * rowSums((risk_never_groups - risk_never_mean)^2,
                                      na.rm = TRUE)
RD_var         <- (G - 1)/G * rowSums((RD_groups - RD_mean)^2,
                                      na.rm = TRUE)

risk_never_se <- sqrt(risk_never_var)
RD_se         <- sqrt(RD_var)


#############################################################
## 5. 95% confidence intervals
#############################################################

z <- qt(0.975, df=G-1)

risk_never_lo <- risk_never_full - z * risk_never_se
risk_never_hi <- risk_never_full + z * risk_never_se

RD_lo <- RD_full - z * RD_se
RD_hi <- RD_full + z * RD_se


#############################################################
## 6. Final summary table
#############################################################

RD_summary <- data.frame(
  horizon       = horizons,
  
  # baseline risks
  risk_never      = risk_never_full,
  risk_never_SE   = risk_never_se,
  risk_never_LO   = risk_never_lo,
  risk_never_HI   = risk_never_hi,
  
  # risk differences
  RD              = RD_full,
  RD_SE           = RD_se,
  RD_LO           = RD_lo,
  RD_HI           = RD_hi
)

print(round(RD_summary, 6))

return(RD_summary)
}

get.surv.summary.trt.year <- function(GJK) {
  
  #############################################################
  ## 0. Helpers
  #############################################################
  
  expit <- function(x) 1 / (1 + exp(-x))
  
  
  #############################################################
  ## 1. Compute risks and survival probabilities
  ##    MSM: year-specific baselines + treatment × year effects
  #############################################################
  
  compute_survival <- function(beta, horizons = c(5,10,15,20)) {
    
    risk_never <- numeric(length(horizons))
    risk_trt   <- numeric(length(horizons))
    surv_never <- numeric(length(horizons))
    surv_trt   <- numeric(length(horizons))
    
    for (j in seq_along(horizons)) {
      
      k <- horizons[j]
      
      b0 <- beta[paste0("year", k)]
      b1 <- beta[paste0("assigned_treatment:year", k)]
      
      # cumulative risks
      p_never <- expit(b0)
      p_trt   <- expit(b0 + b1)
      
      # survival probabilities
      risk_never[j] <- p_never
      risk_trt[j]   <- p_trt
      surv_never[j] <- 1 - p_never
      surv_trt[j]   <- 1 - p_trt
    }
    
    list(
      surv_never = surv_never,
      surv_trt   = surv_trt
    )
  }
  
  
  #############################################################
  ## 2. Full-sample survival estimates
  #############################################################
  
  horizons <- c(5,10,15,20)
  
  full_vals <- compute_survival(GJK$beta_full, horizons)
  
  surv_never_full <- full_vals$surv_never
  surv_trt_full   <- full_vals$surv_trt
  
  
  #############################################################
  ## 3. Grouped jackknife estimates
  #############################################################
  
  G <- nrow(GJK$beta_groups)
  
  surv_never_groups <- matrix(
    NA_real_,
    nrow = length(horizons),
    ncol = G,
    dimnames = list(paste0("year", horizons), paste0("g", seq_len(G)))
  )
  
  surv_trt_groups <- matrix(
    NA_real_,
    nrow = length(horizons),
    ncol = G,
    dimnames = list(paste0("year", horizons), paste0("g", seq_len(G)))
  )
  
  for (g in seq_len(G)) {
    
    tmp <- compute_survival(GJK$beta_groups[g, ], horizons)
    
    surv_never_groups[, g] <- unlist(tmp$surv_never)
    surv_trt_groups[, g]   <- unlist(tmp$surv_trt)
  }
  
  
  #############################################################
  ## 4. Grouped jackknife variance and SE
  #############################################################
  
  surv_never_mean <- rowMeans(surv_never_groups, na.rm = TRUE)
  surv_trt_mean   <- rowMeans(surv_trt_groups, na.rm = TRUE)
  
  surv_never_var <- (G - 1) / G *
    rowSums((surv_never_groups - surv_never_mean)^2, na.rm = TRUE)
  
  surv_trt_var <- (G - 1) / G *
    rowSums((surv_trt_groups - surv_trt_mean)^2, na.rm = TRUE)
  
  surv_never_se <- sqrt(surv_never_var)
  surv_trt_se   <- sqrt(surv_trt_var)
  
  
  #############################################################
  ## 5. 95% confidence intervals
  #############################################################
  
  z <- qt(0.975, df=G-1)
  
  surv_never_lo <- surv_never_full - z * surv_never_se
  surv_never_hi <- surv_never_full + z * surv_never_se
  
  surv_trt_lo <- surv_trt_full - z * surv_trt_se
  surv_trt_hi <- surv_trt_full + z * surv_trt_se
  
  
  #############################################################
  ## 6. Final summary table
  #############################################################
  
  surv_summary <- data.frame(
    horizon = horizons,
    
    # never treated
    surv_never    = surv_never_full,
    surv_never_SE = surv_never_se,
    surv_never_LO = surv_never_lo,
    surv_never_HI = surv_never_hi,
    
    # treated
    surv_trt    = surv_trt_full,
    surv_trt_SE = surv_trt_se,
    surv_trt_LO = surv_trt_lo,
    surv_trt_HI = surv_trt_hi
  )
  
  print(round(surv_summary, 6))
  
  return(surv_summary)
}



print('
## (2)treatment by year model: working.msm<-"Y ~ -1 + year5 + assigned_treatment:year5 + year10 + assigned_treatment:year10 + year15 + assigned_treatment:year15 + year20 + assigned_treatment:year20"
')

print('## 20 jackknife groups; risk for never initiate regime and risk differences by years')

load(file='./results/Group20/group_jackknife_msm_year.rdata')
RD_summary_trtyear_g20<-get.summary.trt.year(GJK)

print('## 20 jackknife groups; survivial probabilities for never initiate regime and initiate regime by years')

Surv_summary_trtyear_g20<-get.surv.summary.trt.year(GJK)





#(3) -------------------------------stratify GRIPHON_eligibility------------------------ 



working.msm<-"Y ~ -1 + year5 + year10 +  year15 + year20 + assigned_treatment + year5:GRIPHON_eligibility__yes +year10:GRIPHON_eligibility__yes+ year15:GRIPHON_eligibility__yes + year20:GRIPHON_eligibility__yes + assigned_treatment:GRIPHON_eligibility__yes "


coef.names<-NULL
horizons = c(5,10,15,20)
for (k in horizons) {
  coef.names       <- c(coef.names, paste0("year", k))
}
coef.names       <- c(coef.names, "assigned_treatment")

for (k in horizons) {
  coef.names       <- c(coef.names, paste0("GRIPHON_eligibility__yes:year", k))
}

coef.names       <- c(coef.names, "GRIPHON_eligibility__yes:assigned_treatment")



GJK.GRIPHON <- grouped_jackknife_ltmleMSM(
  data    = wide_df_ltmle,
  regimes = regimes,
  G       = 20,        # recommended
  ncores  = 20,
  coef_names = coef.names
)

print(GJK.GRIPHON$beta_full)
print(GJK.GRIPHON$gj_se)

save(GJK.GRIPHON, file='./results/Group20/group_jackknife_msm_GRIPHON.rdata')





get.summary.GRIPHON<-function(GJK.GRIPHON)
{
  
## -------------------------------------------
## 0) Setup
## -------------------------------------------

# horizons we want to report
horizons <- c(5, 10, 15, 20)

# expit helper
expit <- function(x) 1 / (1 + exp(-x))


## -------------------------------------------
## 1) Compute risks & RDs for GRIPHON strata (no-intercept MSM)
##
## For each year k:
##   GRIPHON = 0:
##     logit_never = beta[year_k]
##     logit_init  = beta[year_k] + beta[trt:year_k]
##
##   GRIPHON = 1:
##     logit_never = beta[year_k] + beta[GRIPHON:year_k]
##     logit_init  = beta[year_k]  + beta[GRIPHON:year_k]
##                + beta[trt:year_k] + beta[GRIPHON:trt:year_k]
##
## We return:
##   risk_never_g0, risk_never_g1, RD_g0, RD_g1  (vectors over horizons)
## -------------------------------------------

compute_stratified_metrics <- function(beta, horizons = c(5,10,15,20)) {
  
  risk_never_g0 <- RD_g0 <- numeric(length(horizons))
  risk_never_g1 <- RD_g1 <- numeric(length(horizons))
  
  for (j in seq_along(horizons)) {
    k <- horizons[j]
    
    year_k         <- paste0("year", k)
    trt_year_k     <- "assigned_treatment"
    mod_year_k     <- paste0("GRIPHON_eligibility__yes:year", k)
    mod_trt_year_k <- "GRIPHON_eligibility__yes:assigned_treatment"
    
    # GRIPHON=0 (never / init)
    logit_never_g0 <- beta[year_k]
    logit_init_g0  <- beta[year_k] + beta[trt_year_k]
    
    # GRIPHON=1 (never / init) -- note: GRIPHON main-column "GRIPHON_eligibility__yes"
    logit_never_g1 <- beta[year_k] + beta[mod_year_k]
    logit_init_g1  <- beta[year_k] + + beta[mod_year_k] + beta[trt_year_k] + beta[mod_trt_year_k]
    
    # probabilities
    p_never_g0 <- expit(logit_never_g0)
    p_init_g0  <- expit(logit_init_g0)
    p_never_g1 <- expit(logit_never_g1)
    p_init_g1  <- expit(logit_init_g1)
    
    # store
    risk_never_g0[j] <- p_never_g0
    risk_never_g1[j] <- p_never_g1
    RD_g0[j]         <- p_init_g0 - p_never_g0
    RD_g1[j]         <- p_init_g1 - p_never_g1
  }
  
  names(risk_never_g0) <- paste0("risk_never_g0_year", horizons)
  names(risk_never_g1) <- paste0("risk_never_g1_year", horizons)
  names(RD_g0)         <- paste0("RD_g0_year",          horizons)
  names(RD_g1)         <- paste0("RD_g1_year",          horizons)
  
  list(
    risk_never_g0 = risk_never_g0,
    risk_never_g1 = risk_never_g1,
    RD_g0         = RD_g0,
    RD_g1         = RD_g1
  )
}


## -------------------------------------------
## 2) Full-sample estimates from GJK.GRIPHON$beta_full
## -------------------------------------------

full_vals <- compute_stratified_metrics(GJK.GRIPHON$beta_full, horizons)

risk_never_g0_full <- full_vals$risk_never_g0
risk_never_g1_full <- full_vals$risk_never_g1
RD_g0_full         <- full_vals$RD_g0
RD_g1_full         <- full_vals$RD_g1


## -------------------------------------------
## 3) Grouped jackknife (delete-a-group) estimates for each group
##     GJK.GRIPHON$beta_groups is assumed to be G × p
## -------------------------------------------

G <- nrow(GJK.GRIPHON$beta_groups)

risk_never_g0_grp <- matrix(NA_real_, nrow = length(horizons), ncol = G)
risk_never_g1_grp <- matrix(NA_real_, nrow = length(horizons), ncol = G)
RD_g0_grp         <- matrix(NA_real_, nrow = length(horizons), ncol = G)
RD_g1_grp         <- matrix(NA_real_, nrow = length(horizons), ncol = G)

for (g in seq_len(G)) {
  tmp <- compute_stratified_metrics(GJK.GRIPHON$beta_groups[g, ], horizons)
  risk_never_g0_grp[, g] <- unlist(tmp$risk_never_g0)
  risk_never_g1_grp[, g] <- unlist(tmp$risk_never_g1)
  RD_g0_grp[, g]         <- unlist(tmp$RD_g0)
  RD_g1_grp[, g]         <- unlist(tmp$RD_g1)
}

rownames(risk_never_g0_grp) <- rownames(risk_never_g1_grp) <- paste0("year", horizons)
rownames(RD_g0_grp)         <- rownames(RD_g1_grp)         <- paste0("year", horizons)
colnames(risk_never_g0_grp) <- colnames(risk_never_g1_grp) <- paste0("g", seq_len(G))
colnames(RD_g0_grp)         <- colnames(RD_g1_grp)         <- paste0("g", seq_len(G))


## -------------------------------------------
## 4) Grouped jackknife SEs and 95% CIs
##     Var(theta) = (G-1)/G * sum_g (theta_g - mean(theta_g))^2
## -------------------------------------------

# means across groups
risk_never_g0_mean <- rowMeans(risk_never_g0_grp, na.rm = TRUE)
risk_never_g1_mean <- rowMeans(risk_never_g1_grp, na.rm = TRUE)
RD_g0_mean         <- rowMeans(RD_g0_grp,         na.rm = TRUE)
RD_g1_mean         <- rowMeans(RD_g1_grp,         na.rm = TRUE)

# variances
risk_never_g0_var <- (G - 1)/G * rowSums((risk_never_g0_grp - risk_never_g0_mean)^2, na.rm = TRUE)
risk_never_g1_var <- (G - 1)/G * rowSums((risk_never_g1_grp - risk_never_g1_mean)^2, na.rm = TRUE)
RD_g0_var         <- (G - 1)/G * rowSums((RD_g0_grp         - RD_g0_mean        )^2, na.rm = TRUE)
RD_g1_var         <- (G - 1)/G * rowSums((RD_g1_grp         - RD_g1_mean        )^2, na.rm = TRUE)

# SEs
risk_never_g0_se <- sqrt(risk_never_g0_var)
risk_never_g1_se <- sqrt(risk_never_g1_var)
RD_g0_se         <- sqrt(RD_g0_var)
RD_g1_se         <- sqrt(RD_g1_var)

# 95% CIs
z <- qt(0.975, df=G-1)

risk_never_g0_lo <- risk_never_g0_full - z * risk_never_g0_se
risk_never_g0_hi <- risk_never_g0_full + z * risk_never_g0_se
risk_never_g1_lo <- risk_never_g1_full - z * risk_never_g1_se
risk_never_g1_hi <- risk_never_g1_full + z * risk_never_g1_se

RD_g0_lo <- RD_g0_full - z * RD_g0_se
RD_g0_hi <- RD_g0_full + z * RD_g0_se
RD_g1_lo <- RD_g1_full - z * RD_g1_se
RD_g1_hi <- RD_g1_full + z * RD_g1_se


## -------------------------------------------
## 5) Final tidy summary table
##     (one row per horizon, with GRIPHON=0 and GRIPHON=1 columns)
## -------------------------------------------

summary_GRIPHON <- data.frame(
  horizon = horizons,
  
  # Baseline risks (never) for GRIPHON=0
  risk_never_g0     = risk_never_g0_full,
  risk_never_g0_SE  = risk_never_g0_se,
  risk_never_g0_lo  = risk_never_g0_lo,
  risk_never_g0_hi  = risk_never_g0_hi,
  
  # Risk differences (treated - never) for GRIPHON=0
  RD_g0             = RD_g0_full,
  RD_g0_SE          = RD_g0_se,
  RD_g0_lo          = RD_g0_lo,
  RD_g0_hi          = RD_g0_hi,
  
  # Baseline risks (never) for GRIPHON=1
  risk_never_g1     = risk_never_g1_full,
  risk_never_g1_SE  = risk_never_g1_se,
  risk_never_g1_lo  = risk_never_g1_lo,
  risk_never_g1_hi  = risk_never_g1_hi,
  
  # Risk differences (treated - never) for GRIPHON=1
  RD_g1             = RD_g1_full,
  RD_g1_SE          = RD_g1_se,
  RD_g1_lo          = RD_g1_lo,
  RD_g1_hi          = RD_g1_hi
)

print(round(summary_GRIPHON, 6))

return(summary_GRIPHON)
}


get.surv.summary.GRIPHON <- function(GJK.GRIPHON) {
  
  ## -------------------------------------------
  ## 0) Setup
  ## -------------------------------------------
  
  horizons <- c(5, 10, 15, 20)
  expit <- function(x) 1 / (1 + exp(-x))
  
  
  ## -------------------------------------------
  ## 1) Compute survival probabilities by GRIPHON
  ##
  ## MSM with no intercept, year effects, treatment effects,
  ## GRIPHON × year and GRIPHON × treatment interactions
  ## -------------------------------------------
  
  compute_stratified_survival <- function(beta, horizons = c(5,10,15,20)) {
    
    surv_never_g0 <- surv_trt_g0 <- numeric(length(horizons))
    surv_never_g1 <- surv_trt_g1 <- numeric(length(horizons))
    
    for (j in seq_along(horizons)) {
      
      k <- horizons[j]
      
      year_k         <- paste0("year", k)
      trt_year_k     <- "assigned_treatment"
      mod_year_k     <- paste0("GRIPHON_eligibility__yes:year", k)
      mod_trt_year_k <- "GRIPHON_eligibility__yes:assigned_treatment"
      
      ## GRIPHON = 0
      logit_never_g0 <- beta[year_k]
      logit_trt_g0   <- beta[year_k] + beta[trt_year_k]
      
      ## GRIPHON = 1
      logit_never_g1 <- beta[year_k] + beta[mod_year_k]
      logit_trt_g1   <- beta[year_k] + beta[mod_year_k] +
        beta[trt_year_k] + beta[mod_trt_year_k]
      
      ## Convert to survival probabilities
      surv_never_g0[j] <- 1 - expit(logit_never_g0)
      surv_trt_g0[j]   <- 1 - expit(logit_trt_g0)
      surv_never_g1[j] <- 1 - expit(logit_never_g1)
      surv_trt_g1[j]   <- 1 - expit(logit_trt_g1)
    }
    
    list(
      surv_never_g0 = surv_never_g0,
      surv_trt_g0   = surv_trt_g0,
      surv_never_g1 = surv_never_g1,
      surv_trt_g1   = surv_trt_g1
    )
  }
  
  
  ## -------------------------------------------
  ## 2) Full-sample estimates
  ## -------------------------------------------
  
  full_vals <- compute_stratified_survival(GJK.GRIPHON$beta_full, horizons)
  
  surv_never_g0_full <- full_vals$surv_never_g0
  surv_trt_g0_full   <- full_vals$surv_trt_g0
  surv_never_g1_full <- full_vals$surv_never_g1
  surv_trt_g1_full   <- full_vals$surv_trt_g1
  
  
  ## -------------------------------------------
  ## 3) Grouped jackknife estimates
  ## -------------------------------------------
  
  G <- nrow(GJK.GRIPHON$beta_groups)
  
  surv_never_g0_grp <- matrix(NA_real_, length(horizons), G)
  surv_trt_g0_grp   <- matrix(NA_real_, length(horizons), G)
  surv_never_g1_grp <- matrix(NA_real_, length(horizons), G)
  surv_trt_g1_grp   <- matrix(NA_real_, length(horizons), G)
  
  for (g in seq_len(G)) {
    tmp <- compute_stratified_survival(GJK.GRIPHON$beta_groups[g, ], horizons)
    surv_never_g0_grp[, g] <- unlist(tmp$surv_never_g0)
    surv_trt_g0_grp[, g]   <- unlist(tmp$surv_trt_g0)
    surv_never_g1_grp[, g] <- unlist(tmp$surv_never_g1)
    surv_trt_g1_grp[, g]   <- unlist(tmp$surv_trt_g1)
  }
  
  
  ## -------------------------------------------
  ## 4) Grouped jackknife SEs
  ## -------------------------------------------
  
  jk_var <- function(mat) {
    m <- rowMeans(mat, na.rm = TRUE)
    (G - 1) / G * rowSums((mat - m)^2, na.rm = TRUE)
  }
  
  surv_never_g0_se <- sqrt(jk_var(surv_never_g0_grp))
  surv_trt_g0_se   <- sqrt(jk_var(surv_trt_g0_grp))
  surv_never_g1_se <- sqrt(jk_var(surv_never_g1_grp))
  surv_trt_g1_se   <- sqrt(jk_var(surv_trt_g1_grp))
  
  
  ## -------------------------------------------
  ## 5) 95% confidence intervals
  ## -------------------------------------------
  
  z <- qt(0.975, df=G-1)
  
  ci <- function(est, se) {
    cbind(lo = est - z * se, hi = est + z * se)
  }
  
  ci_never_g0 <- ci(surv_never_g0_full, surv_never_g0_se)
  ci_trt_g0   <- ci(surv_trt_g0_full,   surv_trt_g0_se)
  ci_never_g1 <- ci(surv_never_g1_full, surv_never_g1_se)
  ci_trt_g1   <- ci(surv_trt_g1_full,   surv_trt_g1_se)
  
  
  ## -------------------------------------------
  ## 6) Final summary table
  ## -------------------------------------------
  
  summary_GRIPHON <- data.frame(
    horizon = horizons,
    
    surv_never_g0    = surv_never_g0_full,
    surv_never_g0_SE = surv_never_g0_se,
    surv_never_g0_lo = ci_never_g0[, "lo"],
    surv_never_g0_hi = ci_never_g0[, "hi"],
    
    surv_trt_g0    = surv_trt_g0_full,
    surv_trt_g0_SE = surv_trt_g0_se,
    surv_trt_g0_lo = ci_trt_g0[, "lo"],
    surv_trt_g0_hi = ci_trt_g0[, "hi"],
    
    surv_never_g1    = surv_never_g1_full,
    surv_never_g1_SE = surv_never_g1_se,
    surv_never_g1_lo = ci_never_g1[, "lo"],
    surv_never_g1_hi = ci_never_g1[, "hi"],
    
    surv_trt_g1    = surv_trt_g1_full,
    surv_trt_g1_SE = surv_trt_g1_se,
    surv_trt_g1_lo = ci_trt_g1[, "lo"],
    surv_trt_g1_hi = ci_trt_g1[, "hi"]
  )
  
  print(round(summary_GRIPHON, 6))
  return(summary_GRIPHON)
}


print('
## (3)main model modified by GRIPHON_eligibility: working.msm<-"Y ~ -1 + year5 + year10 +  year15 + year20 + assigned_treatment + year5:GRIPHON_eligibility__yes +year10:GRIPHON_eligibility__yes+ year15:GRIPHON_eligibility__yes + year20:GRIPHON_eligibility__yes + assigned_treatment:GRIPHON_eligibility__yes

')

print('## 20 jackknife groups; risk for never initiate regime and risk differences by years')

load(file='./results/Group20/group_jackknife_msm_GRIPHON.rdata')
RD_summary_GRIPHON_g20<-get.summary.GRIPHON(GJK.GRIPHON)


print('## 20 jackknife groups; survivial probabilities for never initiate regime and initiate regime by years')
Surv_summary_GRIPHON_g20<-get.surv.summary.GRIPHON(GJK.GRIPHON)



# (4)-------------------------------stratify by prevalent or incident cases------------------------------------------------


working.msm<-"Y ~ -1 + year5 + year10 +  year15 + year20 + assigned_treatment + year5:incident_or_prevalent__prevalent +year10:incident_or_prevalent__prevalent+ year15:incident_or_prevalent__prevalent + year20:incident_or_prevalent__prevalent + assigned_treatment:incident_or_prevalent__prevalent "



coef.names<-NULL
horizons = c(5,10,15,20)
for (k in horizons) {
  coef.names       <- c(coef.names, paste0("year", k))
}
coef.names       <- c(coef.names, "assigned_treatment")

for (k in horizons) {
  coef.names       <- c(coef.names, paste0("incident_or_prevalent__prevalent:year", k))
}

coef.names       <- c(coef.names, "incident_or_prevalent__prevalent:assigned_treatment")




GJK.prevalent <- grouped_jackknife_ltmleMSM(
  data    = wide_df_ltmle,
  regimes = regimes,
  G       = 20,        # recommended
  ncores  = 20,
  coef_names = coef.names
)

print(GJK.prevalent$beta_full)
print(GJK.prevalent$gj_se)

save(GJK.prevalent, file='./results/Group20/group_jackknife_msm_prevalent.rdata')







get.summary.prevalent<-function(GJK.prevalent)
  
{

## -------------------------------------------
## 0) Setup
## -------------------------------------------

# horizons we want to report
horizons <- c(5, 10, 15, 20)

# expit helper
expit <- function(x) 1 / (1 + exp(-x))


## -------------------------------------------
## 1) Compute risks & RDs for prevalent strata (no-intercept MSM)
##
## For each year k:
##   prevalent = 0:
##     logit_never = beta[year_k]
##     logit_init  = beta[year_k] + beta[trt:year_k]
##
##   prevalent = 1:
##     logit_never = beta[year_k] + beta[prevalent:year_k]
##     logit_init  = beta[year_k]  + beta[prevalent:year_k]
##                + beta[trt:year_k] + beta[prevalent:trt:year_k]
##
## We return:
##   risk_never_g0, risk_never_g1, RD_g0, RD_g1  (vectors over horizons)
## -------------------------------------------

compute_stratified_metrics <- function(beta, horizons = c(5,10,15,20)) {
  
  risk_never_g0 <- RD_g0 <- numeric(length(horizons))
  risk_never_g1 <- RD_g1 <- numeric(length(horizons))
  
  for (j in seq_along(horizons)) {
    k <- horizons[j]
    
    year_k         <- paste0("year", k)
    trt_year_k     <-"assigned_treatment"
    mod_year_k     <- paste0("incident_or_prevalent__prevalent:year", k)
    mod_trt_year_k <- "incident_or_prevalent__prevalent:assigned_treatment"
    
    # prevalent=0 (never / init)
    logit_never_g0 <- beta[year_k]
    logit_init_g0  <- beta[year_k] + beta[trt_year_k]
    
    # prevalent=1 (never / init) -- note: prevalent main-column "prevalent_eligibility__yes"
    logit_never_g1 <- beta[year_k] + beta[mod_year_k]
    logit_init_g1  <- beta[year_k] + + beta[mod_year_k] + beta[trt_year_k] + beta[mod_trt_year_k]
    
    # probabilities
    p_never_g0 <- expit(logit_never_g0)
    p_init_g0  <- expit(logit_init_g0)
    p_never_g1 <- expit(logit_never_g1)
    p_init_g1  <- expit(logit_init_g1)
    
    # store
    risk_never_g0[j] <- p_never_g0
    risk_never_g1[j] <- p_never_g1
    RD_g0[j]         <- p_init_g0 - p_never_g0
    RD_g1[j]         <- p_init_g1 - p_never_g1
  }
  
  names(risk_never_g0) <- paste0("risk_never_g0_year", horizons)
  names(risk_never_g1) <- paste0("risk_never_g1_year", horizons)
  names(RD_g0)         <- paste0("RD_g0_year",          horizons)
  names(RD_g1)         <- paste0("RD_g1_year",          horizons)
  
  list(
    risk_never_g0 = risk_never_g0,
    risk_never_g1 = risk_never_g1,
    RD_g0         = RD_g0,
    RD_g1         = RD_g1
  )
}


## -------------------------------------------
## 2) Full-sample estimates from GJK.prevalent$beta_full
## -------------------------------------------

full_vals <- compute_stratified_metrics(GJK.prevalent$beta_full, horizons)

risk_never_g0_full <- full_vals$risk_never_g0
risk_never_g1_full <- full_vals$risk_never_g1
RD_g0_full         <- full_vals$RD_g0
RD_g1_full         <- full_vals$RD_g1


## -------------------------------------------
## 3) Grouped jackknife (delete-a-group) estimates for each group
##     GJK.prevalent$beta_groups is assumed to be G × p
## -------------------------------------------

G <- nrow(GJK.prevalent$beta_groups)

risk_never_g0_grp <- matrix(NA_real_, nrow = length(horizons), ncol = G)
risk_never_g1_grp <- matrix(NA_real_, nrow = length(horizons), ncol = G)
RD_g0_grp         <- matrix(NA_real_, nrow = length(horizons), ncol = G)
RD_g1_grp         <- matrix(NA_real_, nrow = length(horizons), ncol = G)

for (g in seq_len(G)) {
  tmp <- compute_stratified_metrics(GJK.prevalent$beta_groups[g, ], horizons)
  risk_never_g0_grp[, g] <- unlist(tmp$risk_never_g0)
  risk_never_g1_grp[, g] <- unlist(tmp$risk_never_g1)
  RD_g0_grp[, g]         <- unlist(tmp$RD_g0)
  RD_g1_grp[, g]         <- unlist(tmp$RD_g1)
}

rownames(risk_never_g0_grp) <- rownames(risk_never_g1_grp) <- paste0("year", horizons)
rownames(RD_g0_grp)         <- rownames(RD_g1_grp)         <- paste0("year", horizons)
colnames(risk_never_g0_grp) <- colnames(risk_never_g1_grp) <- paste0("g", seq_len(G))
colnames(RD_g0_grp)         <- colnames(RD_g1_grp)         <- paste0("g", seq_len(G))


## -------------------------------------------
## 4) Grouped jackknife SEs and 95% CIs
##     Var(theta) = (G-1)/G * sum_g (theta_g - mean(theta_g))^2
## -------------------------------------------

# means across groups
risk_never_g0_mean <- rowMeans(risk_never_g0_grp, na.rm = TRUE)
risk_never_g1_mean <- rowMeans(risk_never_g1_grp, na.rm = TRUE)
RD_g0_mean         <- rowMeans(RD_g0_grp,         na.rm = TRUE)
RD_g1_mean         <- rowMeans(RD_g1_grp,         na.rm = TRUE)

# variances
risk_never_g0_var <- (G - 1)/G * rowSums((risk_never_g0_grp - risk_never_g0_mean)^2, na.rm = TRUE)
risk_never_g1_var <- (G - 1)/G * rowSums((risk_never_g1_grp - risk_never_g1_mean)^2, na.rm = TRUE)
RD_g0_var         <- (G - 1)/G * rowSums((RD_g0_grp         - RD_g0_mean        )^2, na.rm = TRUE)
RD_g1_var         <- (G - 1)/G * rowSums((RD_g1_grp         - RD_g1_mean        )^2, na.rm = TRUE)

# SEs
risk_never_g0_se <- sqrt(risk_never_g0_var)
risk_never_g1_se <- sqrt(risk_never_g1_var)
RD_g0_se         <- sqrt(RD_g0_var)
RD_g1_se         <- sqrt(RD_g1_var)

# 95% CIs
z <- qt(0.975, df=G-1)

risk_never_g0_lo <- risk_never_g0_full - z * risk_never_g0_se
risk_never_g0_hi <- risk_never_g0_full + z * risk_never_g0_se
risk_never_g1_lo <- risk_never_g1_full - z * risk_never_g1_se
risk_never_g1_hi <- risk_never_g1_full + z * risk_never_g1_se

RD_g0_lo <- RD_g0_full - z * RD_g0_se
RD_g0_hi <- RD_g0_full + z * RD_g0_se
RD_g1_lo <- RD_g1_full - z * RD_g1_se
RD_g1_hi <- RD_g1_full + z * RD_g1_se


## -------------------------------------------
## 5) Final tidy summary table
##     (one row per horizon, with prevalent=0 and prevalent=1 columns)
## -------------------------------------------

summary_prevalent <- data.frame(
  horizon = horizons,
  
  # Baseline risks (never) for prevalent=0
  risk_never_g0     = risk_never_g0_full,
  risk_never_g0_SE  = risk_never_g0_se,
  risk_never_g0_lo  = risk_never_g0_lo,
  risk_never_g0_hi  = risk_never_g0_hi,
  
  # Risk differences (treated - never) for prevalent=0
  RD_g0             = RD_g0_full,
  RD_g0_SE          = RD_g0_se,
  RD_g0_lo          = RD_g0_lo,
  RD_g0_hi          = RD_g0_hi,
  
  # Baseline risks (never) for prevalent=1
  risk_never_g1     = risk_never_g1_full,
  risk_never_g1_SE  = risk_never_g1_se,
  risk_never_g1_lo  = risk_never_g1_lo,
  risk_never_g1_hi  = risk_never_g1_hi,
  
  # Risk differences (treated - never) for prevalent=1
  RD_g1             = RD_g1_full,
  RD_g1_SE          = RD_g1_se,
  RD_g1_lo          = RD_g1_lo,
  RD_g1_hi          = RD_g1_hi
)

print(round(summary_prevalent, 6))

return(summary_prevalent)
}


get.surv.summary.prevalent <- function(GJK.prevalent) {
  
  ## -------------------------------------------
  ## 0) Setup
  ## -------------------------------------------
  
  horizons <- c(5, 10, 15, 20)
  expit <- function(x) 1 / (1 + exp(-x))
  
  
  ## -------------------------------------------
  ## 1) Compute survival probabilities by prevalent strata
  ## -------------------------------------------
  
  compute_stratified_survival <- function(beta, horizons = c(5,10,15,20)) {
    
    surv_never_g0 <- surv_trt_g0 <- numeric(length(horizons))
    surv_never_g1 <- surv_trt_g1 <- numeric(length(horizons))
    
    for (j in seq_along(horizons)) {
      
      k <- horizons[j]
      
      year_k         <- paste0("year", k)
      trt_year_k     <- "assigned_treatment"
      mod_year_k     <- paste0("incident_or_prevalent__prevalent:year", k)
      mod_trt_year_k <- "incident_or_prevalent__prevalent:assigned_treatment"
      
      ## incident (prevalent = 0)
      logit_never_g0 <- beta[year_k]
      logit_trt_g0   <- beta[year_k] + beta[trt_year_k]
      
      ## prevalent = 1
      logit_never_g1 <- beta[year_k] + beta[mod_year_k]
      logit_trt_g1   <- beta[year_k] + beta[mod_year_k] +
        beta[trt_year_k] + beta[mod_trt_year_k]
      
      ## convert to survival probabilities
      surv_never_g0[j] <- 1 - expit(logit_never_g0)
      surv_trt_g0[j]   <- 1 - expit(logit_trt_g0)
      surv_never_g1[j] <- 1 - expit(logit_never_g1)
      surv_trt_g1[j]   <- 1 - expit(logit_trt_g1)
    }
    
    list(
      surv_never_g0 = surv_never_g0,
      surv_trt_g0   = surv_trt_g0,
      surv_never_g1 = surv_never_g1,
      surv_trt_g1   = surv_trt_g1
    )
  }
  
  
  ## -------------------------------------------
  ## 2) Full-sample estimates
  ## -------------------------------------------
  
  full_vals <- compute_stratified_survival(GJK.prevalent$beta_full, horizons)
  
  surv_never_g0_full <- full_vals$surv_never_g0
  surv_trt_g0_full   <- full_vals$surv_trt_g0
  surv_never_g1_full <- full_vals$surv_never_g1
  surv_trt_g1_full   <- full_vals$surv_trt_g1
  
  
  ## -------------------------------------------
  ## 3) Grouped jackknife estimates
  ## -------------------------------------------
  
  G <- nrow(GJK.prevalent$beta_groups)
  
  surv_never_g0_grp <- matrix(NA_real_, length(horizons), G)
  surv_trt_g0_grp   <- matrix(NA_real_, length(horizons), G)
  surv_never_g1_grp <- matrix(NA_real_, length(horizons), G)
  surv_trt_g1_grp   <- matrix(NA_real_, length(horizons), G)
  
  for (g in seq_len(G)) {
    
    tmp <- compute_stratified_survival(GJK.prevalent$beta_groups[g, ], horizons)
    
    surv_never_g0_grp[, g] <- unlist(tmp$surv_never_g0)
    surv_trt_g0_grp[, g]   <- unlist(tmp$surv_trt_g0)
    surv_never_g1_grp[, g] <- unlist(tmp$surv_never_g1)
    surv_trt_g1_grp[, g]   <- unlist(tmp$surv_trt_g1)
  }
  
  
  ## -------------------------------------------
  ## 4) Grouped jackknife SEs
  ## -------------------------------------------
  
  jk_var <- function(mat) {
    m <- rowMeans(mat, na.rm = TRUE)
    (G - 1) / G * rowSums((mat - m)^2, na.rm = TRUE)
  }
  
  surv_never_g0_se <- sqrt(jk_var(surv_never_g0_grp))
  surv_trt_g0_se   <- sqrt(jk_var(surv_trt_g0_grp))
  surv_never_g1_se <- sqrt(jk_var(surv_never_g1_grp))
  surv_trt_g1_se   <- sqrt(jk_var(surv_trt_g1_grp))
  
  
  ## -------------------------------------------
  ## 5) 95% confidence intervals
  ## -------------------------------------------
  
  z <- qt(0.975, df=G-1)
  
  ci <- function(est, se) {
    cbind(lo = est - z * se, hi = est + z * se)
  }
  
  ci_never_g0 <- ci(surv_never_g0_full, surv_never_g0_se)
  ci_trt_g0   <- ci(surv_trt_g0_full,   surv_trt_g0_se)
  ci_never_g1 <- ci(surv_never_g1_full, surv_never_g1_se)
  ci_trt_g1   <- ci(surv_trt_g1_full,   surv_trt_g1_se)
  
  
  ## -------------------------------------------
  ## 6) Final summary table
  ## -------------------------------------------
  
  summary_prevalent <- data.frame(
    horizon = horizons,
    
    surv_never_g0    = surv_never_g0_full,
    surv_never_g0_SE = surv_never_g0_se,
    surv_never_g0_lo = ci_never_g0[, "lo"],
    surv_never_g0_hi = ci_never_g0[, "hi"],
    
    surv_trt_g0    = surv_trt_g0_full,
    surv_trt_g0_SE = surv_trt_g0_se,
    surv_trt_g0_lo = ci_trt_g0[, "lo"],
    surv_trt_g0_hi = ci_trt_g0[, "hi"],
    
    surv_never_g1    = surv_never_g1_full,
    surv_never_g1_SE = surv_never_g1_se,
    surv_never_g1_lo = ci_never_g1[, "lo"],
    surv_never_g1_hi = ci_never_g1[, "hi"],
    
    surv_trt_g1    = surv_trt_g1_full,
    surv_trt_g1_SE = surv_trt_g1_se,
    surv_trt_g1_lo = ci_trt_g1[, "lo"],
    surv_trt_g1_hi = ci_trt_g1[, "hi"]
  )
  
  print(round(summary_prevalent, 6))
  return(summary_prevalent)
}

print('
## (4)main model modified by prevalent or incident case: working.msm<-"Y ~ -1 + year5 + year10 +  year15 + year20 + assigned_treatment + year5:incident_or_prevalent__prevalent +year10:incident_or_prevalent__prevalent+ year15:incident_or_prevalent__prevalent + year20:incident_or_prevalent__prevalent + assigned_treatment:incident_or_prevalent__prevalent "

')
print('## 20 jackknife groups; risk for never initiate regime and risk differences by years')

load(file='./results/Group20/group_jackknife_msm_prevalent.rdata')

RD_summary_prevalent_g20<-get.summary.prevalent(GJK.prevalent)

print('## 20 jackknife groups; survivial probabilities for never initiate regime and initiate regime by years')
Surv_summary_prevalent_g20<-get.surv.summary.prevalent(GJK.prevalent)



# -------------------------------stratify by RPH or non-RPH------------------------------------------------

working.msm<-"Y ~ -1 + year5 + year10 +  year15 + year20 + assigned_treatment + year5:RPH +year10:RPH+ year15:RPH + year20:RPH + assigned_treatment:RPH "



coef.names<-NULL
horizons = c(5,10,15,20)
for (k in horizons) {
  coef.names       <- c(coef.names, paste0("year", k))
}
coef.names       <- c(coef.names, "assigned_treatment")

for (k in horizons) {
  coef.names       <- c(coef.names, paste0("RPH:year", k))
}

coef.names       <- c(coef.names, "RPH:assigned_treatment")




GJK.RPH <- grouped_jackknife_ltmleMSM(
  data    = wide_df_ltmle,
  regimes = regimes,
  G       = 20,        # recommended
  ncores  = 20,
  coef_names = coef.names
)

print(GJK.RPH$beta_full)
print(GJK.RPH$gj_se)

save(GJK.RPH, file='./results/Group20/group_jackknife_msm_RPH.rdata')





get.summary.RPH<-function(GJK.RPH)
  
{
  
  ## -------------------------------------------
  ## 0) Setup
  ## -------------------------------------------
  
  # horizons we want to report
  horizons <- c(5, 10, 15, 20)
  
  # expit helper
  expit <- function(x) 1 / (1 + exp(-x))
  
  
  ## -------------------------------------------
  ## 1) Compute risks & RDs for RPH strata (no-intercept MSM)
  ##
  ## For each year k:
  ##   RPH = 0:
  ##     logit_never = beta[year_k]
  ##     logit_init  = beta[year_k] + beta[trt:year_k]
  ##
  ##   RPH = 1:
  ##     logit_never = beta[year_k] + beta[RPH:year_k]
  ##     logit_init  = beta[year_k]  + beta[RPH:year_k]
  ##                + beta[trt:year_k] + beta[RPH:trt:year_k]
  ##
  ## We return:
  ##   risk_never_g0, risk_never_g1, RD_g0, RD_g1  (vectors over horizons)
  ## -------------------------------------------
  
  compute_stratified_metrics <- function(beta, horizons = c(5,10,15,20)) {
    
    risk_never_g0 <- RD_g0 <- numeric(length(horizons))
    risk_never_g1 <- RD_g1 <- numeric(length(horizons))
    
    for (j in seq_along(horizons)) {
      k <- horizons[j]
      
      year_k         <- paste0("year", k)
      trt_year_k     <-"assigned_treatment"
      mod_year_k     <- paste0("RPH:year", k)
      mod_trt_year_k <- "RPH:assigned_treatment"
      
      # RPH=0 (never / init)
      logit_never_g0 <- beta[year_k]
      logit_init_g0  <- beta[year_k] + beta[trt_year_k]
      
      # RPH=1 (never / init) -- note: RPH main-column "RPH_eligibility__yes"
      logit_never_g1 <- beta[year_k] + beta[mod_year_k]
      logit_init_g1  <- beta[year_k] + + beta[mod_year_k] + beta[trt_year_k] + beta[mod_trt_year_k]
      
      # probabilities
      p_never_g0 <- expit(logit_never_g0)
      p_init_g0  <- expit(logit_init_g0)
      p_never_g1 <- expit(logit_never_g1)
      p_init_g1  <- expit(logit_init_g1)
      
      # store
      risk_never_g0[j] <- p_never_g0
      risk_never_g1[j] <- p_never_g1
      RD_g0[j]         <- p_init_g0 - p_never_g0
      RD_g1[j]         <- p_init_g1 - p_never_g1
    }
    
    names(risk_never_g0) <- paste0("risk_never_g0_year", horizons)
    names(risk_never_g1) <- paste0("risk_never_g1_year", horizons)
    names(RD_g0)         <- paste0("RD_g0_year",          horizons)
    names(RD_g1)         <- paste0("RD_g1_year",          horizons)
    
    list(
      risk_never_g0 = risk_never_g0,
      risk_never_g1 = risk_never_g1,
      RD_g0         = RD_g0,
      RD_g1         = RD_g1
    )
  }
  
  
  ## -------------------------------------------
  ## 2) Full-sample estimates from GJK.RPH$beta_full
  ## -------------------------------------------
  
  full_vals <- compute_stratified_metrics(GJK.RPH$beta_full, horizons)
  
  risk_never_g0_full <- full_vals$risk_never_g0
  risk_never_g1_full <- full_vals$risk_never_g1
  RD_g0_full         <- full_vals$RD_g0
  RD_g1_full         <- full_vals$RD_g1
  
  
  ## -------------------------------------------
  ## 3) Grouped jackknife (delete-a-group) estimates for each group
  ##     GJK.RPH$beta_groups is assumed to be G × p
  ## -------------------------------------------
  
  G <- nrow(GJK.RPH$beta_groups)
  
  risk_never_g0_grp <- matrix(NA_real_, nrow = length(horizons), ncol = G)
  risk_never_g1_grp <- matrix(NA_real_, nrow = length(horizons), ncol = G)
  RD_g0_grp         <- matrix(NA_real_, nrow = length(horizons), ncol = G)
  RD_g1_grp         <- matrix(NA_real_, nrow = length(horizons), ncol = G)
  
  for (g in seq_len(G)) {
    tmp <- compute_stratified_metrics(GJK.RPH$beta_groups[g, ], horizons)
    risk_never_g0_grp[, g] <- unlist(tmp$risk_never_g0)
    risk_never_g1_grp[, g] <- unlist(tmp$risk_never_g1)
    RD_g0_grp[, g]         <- unlist(tmp$RD_g0)
    RD_g1_grp[, g]         <- unlist(tmp$RD_g1)
  }
  
  rownames(risk_never_g0_grp) <- rownames(risk_never_g1_grp) <- paste0("year", horizons)
  rownames(RD_g0_grp)         <- rownames(RD_g1_grp)         <- paste0("year", horizons)
  colnames(risk_never_g0_grp) <- colnames(risk_never_g1_grp) <- paste0("g", seq_len(G))
  colnames(RD_g0_grp)         <- colnames(RD_g1_grp)         <- paste0("g", seq_len(G))
  
  
  ## -------------------------------------------
  ## 4) Grouped jackknife SEs and 95% CIs
  ##     Var(theta) = (G-1)/G * sum_g (theta_g - mean(theta_g))^2
  ## -------------------------------------------
  
  # means across groups
  risk_never_g0_mean <- rowMeans(risk_never_g0_grp, na.rm = TRUE)
  risk_never_g1_mean <- rowMeans(risk_never_g1_grp, na.rm = TRUE)
  RD_g0_mean         <- rowMeans(RD_g0_grp,         na.rm = TRUE)
  RD_g1_mean         <- rowMeans(RD_g1_grp,         na.rm = TRUE)
  
  # variances
  risk_never_g0_var <- (G - 1)/G * rowSums((risk_never_g0_grp - risk_never_g0_mean)^2, na.rm = TRUE)
  risk_never_g1_var <- (G - 1)/G * rowSums((risk_never_g1_grp - risk_never_g1_mean)^2, na.rm = TRUE)
  RD_g0_var         <- (G - 1)/G * rowSums((RD_g0_grp         - RD_g0_mean        )^2, na.rm = TRUE)
  RD_g1_var         <- (G - 1)/G * rowSums((RD_g1_grp         - RD_g1_mean        )^2, na.rm = TRUE)
  
  # SEs
  risk_never_g0_se <- sqrt(risk_never_g0_var)
  risk_never_g1_se <- sqrt(risk_never_g1_var)
  RD_g0_se         <- sqrt(RD_g0_var)
  RD_g1_se         <- sqrt(RD_g1_var)
  
  # 95% CIs
  z <- qt(0.975, df=G-1)
  
  risk_never_g0_lo <- risk_never_g0_full - z * risk_never_g0_se
  risk_never_g0_hi <- risk_never_g0_full + z * risk_never_g0_se
  risk_never_g1_lo <- risk_never_g1_full - z * risk_never_g1_se
  risk_never_g1_hi <- risk_never_g1_full + z * risk_never_g1_se
  
  RD_g0_lo <- RD_g0_full - z * RD_g0_se
  RD_g0_hi <- RD_g0_full + z * RD_g0_se
  RD_g1_lo <- RD_g1_full - z * RD_g1_se
  RD_g1_hi <- RD_g1_full + z * RD_g1_se
  
  
  ## -------------------------------------------
  ## 5) Final tidy summary table
  ##     (one row per horizon, with RPH=0 and RPH=1 columns)
  ## -------------------------------------------
  
  summary_RPH <- data.frame(
    horizon = horizons,
    
    # Baseline risks (never) for RPH=0
    risk_never_g0     = risk_never_g0_full,
    risk_never_g0_SE  = risk_never_g0_se,
    risk_never_g0_lo  = risk_never_g0_lo,
    risk_never_g0_hi  = risk_never_g0_hi,
    
    # Risk differences (treated - never) for RPH=0
    RD_g0             = RD_g0_full,
    RD_g0_SE          = RD_g0_se,
    RD_g0_lo          = RD_g0_lo,
    RD_g0_hi          = RD_g0_hi,
    
    # Baseline risks (never) for RPH=1
    risk_never_g1     = risk_never_g1_full,
    risk_never_g1_SE  = risk_never_g1_se,
    risk_never_g1_lo  = risk_never_g1_lo,
    risk_never_g1_hi  = risk_never_g1_hi,
    
    # Risk differences (treated - never) for RPH=1
    RD_g1             = RD_g1_full,
    RD_g1_SE          = RD_g1_se,
    RD_g1_lo          = RD_g1_lo,
    RD_g1_hi          = RD_g1_hi
  )
  
  print(round(summary_RPH, 6))
  
  return(summary_RPH)
}


get.surv.summary.RPH <- function(GJK.RPH) {
  
  ## -------------------------------------------
  ## 0) Setup
  ## -------------------------------------------
  
  horizons <- c(5, 10, 15, 20)
  expit <- function(x) 1 / (1 + exp(-x))
  
  
  ## -------------------------------------------
  ## 1) Compute survival probabilities by RPH strata
  ## -------------------------------------------
  
  compute_stratified_survival <- function(beta, horizons = c(5,10,15,20)) {
    
    surv_never_g0 <- surv_trt_g0 <- numeric(length(horizons))
    surv_never_g1 <- surv_trt_g1 <- numeric(length(horizons))
    
    for (j in seq_along(horizons)) {
      
      k <- horizons[j]
      
      year_k         <- paste0("year", k)
      trt_year_k     <- "assigned_treatment"
      mod_year_k     <- paste0("RPH:year", k)
      mod_trt_year_k <- "RPH:assigned_treatment"
      
      ## incident (RPH = 0)
      logit_never_g0 <- beta[year_k]
      logit_trt_g0   <- beta[year_k] + beta[trt_year_k]
      
      ## RPH = 1
      logit_never_g1 <- beta[year_k] + beta[mod_year_k]
      logit_trt_g1   <- beta[year_k] + beta[mod_year_k] +
        beta[trt_year_k] + beta[mod_trt_year_k]
      
      ## convert to survival probabilities
      surv_never_g0[j] <- 1 - expit(logit_never_g0)
      surv_trt_g0[j]   <- 1 - expit(logit_trt_g0)
      surv_never_g1[j] <- 1 - expit(logit_never_g1)
      surv_trt_g1[j]   <- 1 - expit(logit_trt_g1)
    }
    
    list(
      surv_never_g0 = surv_never_g0,
      surv_trt_g0   = surv_trt_g0,
      surv_never_g1 = surv_never_g1,
      surv_trt_g1   = surv_trt_g1
    )
  }
  
  
  ## -------------------------------------------
  ## 2) Full-sample estimates
  ## -------------------------------------------
  
  full_vals <- compute_stratified_survival(GJK.RPH$beta_full, horizons)
  
  surv_never_g0_full <- full_vals$surv_never_g0
  surv_trt_g0_full   <- full_vals$surv_trt_g0
  surv_never_g1_full <- full_vals$surv_never_g1
  surv_trt_g1_full   <- full_vals$surv_trt_g1
  
  
  ## -------------------------------------------
  ## 3) Grouped jackknife estimates
  ## -------------------------------------------
  
  G <- nrow(GJK.RPH$beta_groups)
  
  surv_never_g0_grp <- matrix(NA_real_, length(horizons), G)
  surv_trt_g0_grp   <- matrix(NA_real_, length(horizons), G)
  surv_never_g1_grp <- matrix(NA_real_, length(horizons), G)
  surv_trt_g1_grp   <- matrix(NA_real_, length(horizons), G)
  
  for (g in seq_len(G)) {
    
    tmp <- compute_stratified_survival(GJK.RPH$beta_groups[g, ], horizons)
    
    surv_never_g0_grp[, g] <- unlist(tmp$surv_never_g0)
    surv_trt_g0_grp[, g]   <- unlist(tmp$surv_trt_g0)
    surv_never_g1_grp[, g] <- unlist(tmp$surv_never_g1)
    surv_trt_g1_grp[, g]   <- unlist(tmp$surv_trt_g1)
  }
  
  
  ## -------------------------------------------
  ## 4) Grouped jackknife SEs
  ## -------------------------------------------
  
  jk_var <- function(mat) {
    m <- rowMeans(mat, na.rm = TRUE)
    (G - 1) / G * rowSums((mat - m)^2, na.rm = TRUE)
  }
  
  surv_never_g0_se <- sqrt(jk_var(surv_never_g0_grp))
  surv_trt_g0_se   <- sqrt(jk_var(surv_trt_g0_grp))
  surv_never_g1_se <- sqrt(jk_var(surv_never_g1_grp))
  surv_trt_g1_se   <- sqrt(jk_var(surv_trt_g1_grp))
  
  
  ## -------------------------------------------
  ## 5) 95% confidence intervals
  ## -------------------------------------------
  
  z <- qt(0.975, df=G-1)
  
  ci <- function(est, se) {
    cbind(lo = est - z * se, hi = est + z * se)
  }
  
  ci_never_g0 <- ci(surv_never_g0_full, surv_never_g0_se)
  ci_trt_g0   <- ci(surv_trt_g0_full,   surv_trt_g0_se)
  ci_never_g1 <- ci(surv_never_g1_full, surv_never_g1_se)
  ci_trt_g1   <- ci(surv_trt_g1_full,   surv_trt_g1_se)
  
  
  ## -------------------------------------------
  ## 6) Final summary table
  ## -------------------------------------------
  
  summary_RPH <- data.frame(
    horizon = horizons,
    
    surv_never_g0    = surv_never_g0_full,
    surv_never_g0_SE = surv_never_g0_se,
    surv_never_g0_lo = ci_never_g0[, "lo"],
    surv_never_g0_hi = ci_never_g0[, "hi"],
    
    surv_trt_g0    = surv_trt_g0_full,
    surv_trt_g0_SE = surv_trt_g0_se,
    surv_trt_g0_lo = ci_trt_g0[, "lo"],
    surv_trt_g0_hi = ci_trt_g0[, "hi"],
    
    surv_never_g1    = surv_never_g1_full,
    surv_never_g1_SE = surv_never_g1_se,
    surv_never_g1_lo = ci_never_g1[, "lo"],
    surv_never_g1_hi = ci_never_g1[, "hi"],
    
    surv_trt_g1    = surv_trt_g1_full,
    surv_trt_g1_SE = surv_trt_g1_se,
    surv_trt_g1_lo = ci_trt_g1[, "lo"],
    surv_trt_g1_hi = ci_trt_g1[, "hi"]
  )
  
  print(round(summary_RPH, 6))
  return(summary_RPH)
}

print('
## (4)main model modified by RPH : working.msm<-"Y ~ -1 + year5 + year10 +  year15 + year20 + assigned_treatment + year5:RPH +year10:RPH+ year15:RPH + year20:RPH + assigned_treatment:RPH "

')
print('## 20 jackknife groups; risk for never initiate regime and risk differences by years')

load(file='./results/Group20/group_jackknife_msm_RPH.rdata')

RD_summary_RPH_g20<-get.summary.RPH(GJK.RPH)

print('## 20 jackknife groups; survivial probabilities for never initiate regime and initiate regime by years')
Surv_summary_RPH_g20<-get.surv.summary.RPH(GJK.RPH)





# -----------------------------
#  (6) ------ Jackknife CIs for MSM betas--------------------------------
# -----------------------------
# GROUPED JK (delete-a-group)


jk_ci_msm_group <- function(beta_full, beta_groups, level = 0.95) {
  stopifnot(is.numeric(beta_full), is.matrix(beta_groups) || is.data.frame(beta_groups))
  beta_groups <- as.matrix(beta_groups)
  
  # align coefficient names across full and grouped fits
  common <- intersect(names(beta_full), colnames(beta_groups))
  if (length(common) == 0L) stop("No common coefficient names between beta_full and beta_groups.")
  if (length(common) < length(beta_full)) {
    warning("Some coefficients in beta_full are not present in beta_groups; returning common subset only.")
  }
  
  est <- beta_full[common]
  B   <- beta_groups[, common, drop = FALSE]
  
  G <- nrow(B)
  if (G < 2L) stop("beta_groups must have at least 2 rows (groups).")
  
  # grouped-jackknife variance: (G-1)/G * sum_g (theta_g - mean)^2
  Bmean <- colMeans(B, na.rm = TRUE)
  var_jk <- (G - 1) / G * colSums((t(t(B) - Bmean))^2, na.rm = TRUE)
  se_jk  <- sqrt(var_jk)
  
  # CIs and p-values (Wald, normal approx)
  alpha <- 1 - level
  z     <- qt(1 - alpha/2, df=G-1)
  lo    <- est - z * se_jk
  hi    <- est + z * se_jk
  pval  <- 2 * pt(-abs(est / se_jk), df=G-1)
  
  out <- data.frame(
    coefficient = common,
    estimate    = as.numeric(est),
    SE_jk       = as.numeric(se_jk),
    CI_lo       = as.numeric(lo),
    CI_hi       = as.numeric(hi),
    p_value     = as.numeric(pval),
    row.names   = NULL
  )
  out[order(out$coefficient), ]
}

# main model results
load(file='./results/Group20/group_jackknife_msm_main.rdata')

msm_ci_group <- jk_ci_msm_group(
  beta_full   = GJK.main$beta_full,
  beta_groups = GJK.main$beta_groups,
  level       = 0.95
)


print('## 20 jackknife groups; MSM coefficients for main model')

print(msm_ci_group, digits = 2)



load(file='./results/Group20/group_jackknife_msm_year.rdata')
## treatment by year results
msm_ci_group <- jk_ci_msm_group(
  beta_full   = GJK$beta_full,
  beta_groups = GJK$beta_groups,
  level       = 0.95
)


print('## 20 jackknife groups; MSM coefficients for treatment*year model')
print(msm_ci_group, digits = 2)






load(file='./results/Group20/group_jackknife_msm_GRIPHON.rdata')

msm_ci_group <- jk_ci_msm_group(
  beta_full   = GJK.GRIPHON$beta_full,
  beta_groups = GJK.GRIPHON$beta_groups,
  level       = 0.95
)


print('## 20 jackknife groups; MSM coefficients for  model stratified by GRIPHON')
print(msm_ci_group, digits = 2)




load(file='./results/Group20/group_jackknife_msm_prevalent.rdata')

msm_ci_group.prevalent <- jk_ci_msm_group(
  beta_full   = GJK.prevalent$beta_full,
  beta_groups = GJK.prevalent$beta_groups,
  level       = 0.95
)

print('## 20 jackknife groups; MSM coefficients for  model stratified by prevalent case')

print(msm_ci_group.prevalent, digits = 2)


load(file='./results/Group20/group_jackknife_msm_RPH.rdata')

msm_ci_group.RPH <- jk_ci_msm_group(
  beta_full   = GJK.RPH$beta_full,
  beta_groups = GJK.RPH$beta_groups,
  level       = 0.95
)

print('## 20 jackknife groups; MSM coefficients for  model stratified by RPH case')

print(msm_ci_group.RPH, digits = 2)


#------------------test whether treatment effects vary by time -------------------

#' Test whether treatment effects vary over time (F-distribution adjusted)
#'
#' @param GJK Output object from grouped_jackknife_ltmleMSM
#' @param trt_year_prefix String prefix for treatment-by-time interaction coefficients
#' @param years Vector of time horizons (e.g., c(5, 10, 15, 20))
#' @param ref_year Reference year to compare all other treatment effects against
#'
#' @return A list containing the F-statistic, degrees of freedom, p-value, and contrast table
test.const.trt <- function(GJK, 
                           trt_year_prefix = "assigned_treatment:year",
                           years = c(5, 10, 15, 20), 
                           ref_year = 5) {
  
  ## 1. Helper: Compute jackknife covariance & G_eff handling NAs
  jk_cov_info <- function(beta_groups, required_cols) {
    B <- as.matrix(beta_groups[, required_cols, drop = FALSE])
    
    # Filter complete rows for the required parameters
    valid_rows <- complete.cases(B)
    B_valid    <- B[valid_rows, , drop = FALSE]
    G_eff      <- nrow(B_valid)
    
    if (G_eff < 2) stop("Insufficient valid jackknife replicates to compute covariance.")
    
    B_centered <- sweep(B_valid, 2, colMeans(B_valid), FUN = "-")
    V_jk       <- ((G_eff - 1) / G_eff) * (t(B_centered) %*% B_centered)
    
    list(V_jk = V_jk, G_eff = G_eff)
  }
  
  ## 2. Helper: Build contrast matrix R and compute F-test
  test_constant_trt_effect <- function(beta_full, beta_groups, trt_year_prefix, years, ref_year) {
    
    nm <- names(beta_full)
    trt_names <- paste0(trt_year_prefix, years)
    
    # Check that specified interaction coefficients exist
    missing_coefs <- setdiff(trt_names, nm)
    if (length(missing_coefs) > 0) {
      stop("The following treatment terms were not found in beta_full: ", 
           paste(missing_coefs, collapse = ", "))
    }
    
    # Reference contrast setup
    ref    <- paste0(trt_year_prefix, ref_year)
    others <- setdiff(trt_names, ref)
    q      <- length(others) # Number of contrasts being tested
    
    # Construct contrast matrix R (q x p)
    p <- length(beta_full)
    R <- matrix(0, nrow = q, ncol = p, dimnames = list(others, nm))
    for (i in seq_len(q)) {
      R[i, others[i]] <- 1
      R[i, ref]       <- -1
    }
    
    # Compute covariance for required terms
    cov_info <- jk_cov_info(beta_groups, required_cols = nm)
    V        <- cov_info$V_jk
    G_eff    <- cov_info$G_eff
    
    if (G_eff <= q) {
      stop(sprintf("Insufficient groups (G_eff = %d) to test q = %d contrasts. Require G_eff > q.", G_eff, q))
    }
    
    # Linear contrast variance: R * V * R^T
    beta <- beta_full[nm]
    RVRT <- R %*% V %*% t(R)
    
    # Invert contrast covariance matrix
    RVRT_inv <- tryCatch(
      solve(RVRT),
      error = function(e) {
        warning("Used generalized inverse for Wald test due to singular contrast matrix.")
        MASS::ginv(RVRT)
      }
    )
    
    diff_vec <- as.numeric(R %*% beta)
    
    # Hotelling's T^2 statistic
    T2 <- drop(t(diff_vec) %*% RVRT_inv %*% diff_vec)
    
    # Convert T^2 to F-statistic: F ~ F(df1 = q, df2 = G_eff - q)
    df1    <- q
    df2    <- G_eff - q
    F_stat <- ((df2) / (q * (G_eff - 1))) * T2
    pval   <- pf(F_stat, df1 = df1, df2 = df2, lower.tail = FALSE)
    
    contrasts_df <- data.frame(
      contrast = paste0(others, " - ", ref),
      diff     = diff_vec,
      row.names = NULL
    )
    
    list(
      F_statistic = F_stat,
      T2_statistic = T2,
      df1         = df1,
      df2         = df2,
      G_eff       = G_eff,
      p.value     = pval,
      contrasts   = contrasts_df
    )
  }
  
  ## 3. Run test and print results
  wald_out <- test_constant_trt_effect(
    beta_full       = GJK$beta_full,
    beta_groups     = GJK$beta_groups,
    trt_year_prefix = trt_year_prefix,
    years           = years,
    ref_year        = ref_year
  )
  
  cat(sprintf(
    "\n=== Grouped Jackknife F-Test for Constant Treatment Effect ===\nF(%d, %d) = %.3f (T^2 = %.3f)\np-value = %.4g (G_eff = %d)\n\n",
    wald_out$df1, wald_out$df2, wald_out$F_statistic, wald_out$T2_statistic, wald_out$p.value, wald_out$G_eff
  ))
  
  print(wald_out$contrasts)
  
  return(invisible(wald_out))
}




print('##### test if treatement effects are constant over time')
load(file='./results/Group20/group_jackknife_msm_year.rdata')
print('20 jackknife groups')
test.const.trt(GJK)




