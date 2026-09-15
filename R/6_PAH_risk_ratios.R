get.summary.main.RR <- function(GJK.main) {
  
  # logistic expit function
  expit <- function(x) 1 / (1 + exp(-x))
  
  #############################################################
  ## 1. Compute baseline risks, RD, and RR
  #############################################################
  
  compute_risk_RD_RR <- function(beta, horizons = c(5,10,15,20)) {
    
    risk_never <- numeric(length(horizons))
    RD         <- numeric(length(horizons))
    RR         <- numeric(length(horizons))
    
    for (j in seq_along(horizons)) {
      
      k <- horizons[j]
      
      b0 <- beta[paste0("year", k)]
      b1 <- beta["assigned_treatment"]
      
      p_never <- expit(b0)
      p_init  <- expit(b0 + b1)
      
      risk_never[j] <- p_never
      RD[j]         <- p_init - p_never
      RR[j]         <- ifelse(p_never > 0, p_init / p_never, NA)
    }
    
    names(risk_never) <- paste0("risk_never_year", horizons)
    names(RD)         <- paste0("RD_year",         horizons)
    names(RR)         <- paste0("RR_year",         horizons)
    
    list(
      risk_never = risk_never,
      RD         = RD,
      RR         = RR
    )
  }
  
  #############################################################
  ## 2. Full-sample estimates
  #############################################################
  
  horizons <- c(5,10,15,20)
  
  full_vals <- compute_risk_RD_RR(GJK.main$beta_full, horizons)
  
  risk_never_full <- full_vals$risk_never
  RD_full         <- full_vals$RD
  RR_full         <- full_vals$RR
  
  #############################################################
  ## 3. Grouped jackknife estimates
  #############################################################
  
  G <- nrow(GJK.main$beta_groups)
  
  risk_never_groups <- matrix(NA, length(horizons), G)
  RD_groups         <- matrix(NA, length(horizons), G)
  RR_groups         <- matrix(NA, length(horizons), G)
  
  for (g in 1:G) {
    
    tmp <- compute_risk_RD_RR(GJK.main$beta_groups[g,], horizons)
    
    risk_never_groups[, g] <- unlist(tmp$risk_never)
    RD_groups[, g]         <- unlist(tmp$RD)
    RR_groups[, g]         <- unlist(tmp$RR)
  }
  
  #############################################################
  ## 4. Variance + SE (jackknife)
  #############################################################
  
  risk_never_mean <- rowMeans(risk_never_groups, na.rm = TRUE)
  RD_mean         <- rowMeans(RD_groups, na.rm = TRUE)
  RR_mean         <- rowMeans(RR_groups, na.rm = TRUE)
  
  risk_never_var <- (G - 1)/G * rowSums((risk_never_groups - risk_never_mean)^2)
  RD_var         <- (G - 1)/G * rowSums((RD_groups - RD_mean)^2)
  RR_var         <- (G - 1)/G * rowSums((RR_groups - RR_mean)^2)
  
  risk_never_se <- sqrt(risk_never_var)
  RD_se         <- sqrt(RD_var)
  RR_se         <- sqrt(RR_var)
  
  #############################################################
  ## 5. Confidence intervals
  #############################################################
  
  z <- qt(0.975, df=G-1)
  
  risk_never_lo <- risk_never_full - z * risk_never_se
  risk_never_hi <- risk_never_full + z * risk_never_se
  
  RD_lo <- RD_full - z * RD_se
  RD_hi <- RD_full + z * RD_se
  
  RR_lo <- RR_full - z * RR_se
  RR_hi <- RR_full + z * RR_se
  
  #############################################################
  ## 6. Final summary
  #############################################################
  
  summary <- data.frame(
    horizon = horizons,
    
    # risks
    risk_never    = risk_never_full,
    risk_never_SE = risk_never_se,
    risk_never_LO = risk_never_lo,
    risk_never_HI = risk_never_hi,
    
    # RD
    RD    = RD_full,
    RD_SE = RD_se,
    RD_LO = RD_lo,
    RD_HI = RD_hi,
    
    # RR
    RR    = RR_full,
    RR_SE = RR_se,
    RR_LO = RR_lo,
    RR_HI = RR_hi
  )
  
  print(round(summary, 6))
  return(summary)
}


print('## (1) main model: working.msm<-"Y ~ -1 + year5  + year10 + year15  + year20 + assigned_treatment"')




print('## 20 jackknife groups; risk for never initiate regime and risk ratios by years')

load(file='./results/Group20/group_jackknife_msm_main.rdata')

RR_summary_main_g20<-get.summary.main.RR(GJK.main)



