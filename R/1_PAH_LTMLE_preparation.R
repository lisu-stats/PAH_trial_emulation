#----------------------------------------------------------------------------------------------------
#----------------------------prepare LTMLE arguments for model fitting ---------------------------------------------
#----------------------------------------------------------------------------------------------------

## =========================================================
## 0) Packages
## =========================================================
library(ltmle)  

## =========================================================
## 1) Create dummy variables for baseline categorical variables---------------------
## =========================================================



## Baseline W (factors + continuous we listed)

# Baseline factor variables we currently include in Wvars:
W_factor <- c("sex","diagnosis_verified","vasoresponder","GRIPHON_eligibility", "incident_or_prevalent")

# Turn character → factor and add an explicit NA level "UNK" (unknown)
to_factor_with_unknown <- function(x, unk = "UNK") {
  xf <- addNA(factor(x), ifany = TRUE)   # adds NA as a factor level
  lv <- levels(xf)
  lv[is.na(lv)] <- unk                   # rename the NA level
  levels(xf) <- lv
  xf
}

# Example for one variable:
wide_df_scaled$incident_or_prevalent <- to_factor_with_unknown(wide_df_scaled$incident_or_prevalent)

# Build dummies for each, drop the first level 
make_dummies <- function(x, prefix) {
  x <- factor(x)                         # ensure factor
  mm <- model.matrix(~ x - 1)            # one column per level
  colnames(mm) <- paste0(prefix, "__", sub("^x", "", colnames(mm)))
  if (ncol(mm) >= 2) mm <- mm[, -1, drop = FALSE]  # drop 1st level as reference
  mm
}

# Apply to our final data used for ltmle (wide_df_scaled) or earlier if we prefer:
for (w in W_factor) {
  if (w %in% names(wide_df_scaled)) {
    mm <- make_dummies(wide_df_scaled[[w]], w)
    for (j in seq_len(ncol(mm))) {
      wide_df_scaled[[colnames(mm)[j]]] <- mm[, j]
    }
  }
}

# Remove the original factor columns from the data and from Wvars
wide_df_scaled <- wide_df_scaled[, !names(wide_df_scaled) %in% W_factor]
W_dummies <- grep(paste0("^(", paste(W_factor, collapse="|"), ")__"), 
                  names(wide_df_scaled), value = TRUE)

# New Wvars to use in Qform/gform:
Wvars_cont   <- c("CCI","age_diagnosis","time_to_mono_or_dual_therapy")
Wvars <- c("id_is_OC", "RPH", W_dummies, Wvars_cont)





## =========================================================
## 2) Regimes (n x K* x 2): observed A-path vs. never initiate-------
##    ltmleMSM requires regimes as array [n x numACnodes x numRegimes]
## =========================================================

Anodes <- paste0("A", 1:K)
Cnodes <- paste0("C", 1:(K-1))
Ynodes <- paste0("Y", 1:K)


n <- nrow(wide_df_scaled)
regimes <- array(0L, dim = c(n, K, 2))

# Regime 1: "initiate as eligible" (assume observed data follow this rule) => use observed A path
Aobs <- as.matrix(wide_df_scaled[, Anodes, drop = FALSE])
Aobs[!is.na(Aobs) & Aobs != 0 & Aobs != 1] <- as.numeric(Aobs[!is.na(Aobs) & Aobs != 0 & Aobs != 1] >= 0.5)
regimes[, , 1] <- Aobs

# Regime 2: "never initiate" (regimes[,,2] already zeros)








## =========================================================
##----------- 3) Build  Lnodes = numeric time-varying + functional_class dummies------
## =========================================================


## ---- Collect numeric L candidates: any numeric column ending in _<t> that is not A/C/Y ----
ends_with_time <- grepl("_[0-9]+$", names(wide_df_scaled))
is_num <- sapply(wide_df_scaled, is.numeric)
is_A   <- names(wide_df_scaled) %in% Anodes
is_Y   <- names(wide_df_scaled) %in% Ynodes
is_C   <- names(wide_df_scaled) %in% Cnodes

L_candidates <- names(wide_df_scaled)[ends_with_time & is_num & !(is_A | is_Y | is_C)]

## ---- Build Lnodes by exact time order (1..K) ----
Lnodes_ordered <- unlist(lapply(1:K, function(t) {
  # any numeric time-varying var that ends with _t (e.g., bs_bmi_7, reveal_lite_2_7)
  L_num_t <- L_candidates[grepl(paste0("_", t, "$"), L_candidates)]
  
  L_num_t
}))


# Final Lnodes (numeric + FC dummies)
Lnodes <- Lnodes_ordered

## =========================================================
## Fix Cnodes: force to factor(censored, uncensored)
## =========================================================


Cnodes <- paste0("C", 1:(K-1))
C_present <- Cnodes[Cnodes %in% names(wide_df_scaled)]

coerce_C <- function(x) {
  # If factor, convert to character to clean
  if (is.factor(x)) x <- as.character(x)
  x <- trimws(as.character(x))
  
  # Numeric-like? Convert quietly
  suppressWarnings({
    xnum <- as.numeric(x)
  })
  # Map numeric 0/1
  out <- x
  out[!is.na(xnum) & xnum == 1] <- "uncensored"
  out[!is.na(xnum) & xnum == 0] <- "censored"
  
  # Map common strings
  lower <- tolower(out)
  out[lower %in% c("uncensored","true","t","yes","y")] <- "uncensored"
  out[lower %in% c("censored","false","f","no","n")]   <- "censored"
  
  # Anything else becomes NA (ltmle allows NA in censoring after event)
  out[!(out %in% c("censored","uncensored"))] <- NA_character_
  
  factor(out, levels = c("censored","uncensored"))
}

# Apply to all Cnodes
for (cn in C_present) {
  wide_df_scaled[[cn]] <- coerce_C(wide_df_scaled[[cn]])
}

# Quick assert
stopifnot(all(sapply(C_present, function(cn) {
  x <- wide_df_scaled[[cn]]
  is.factor(x) && identical(levels(x), c("censored","uncensored"))
})))



## =========================================================
#---------(4) FIX COLUMN ORDERING FOR LTMLE / LTMLEMSM----------------
## =========================================================

reorder_for_ltmle <- function(df, Kstar, Wvars, Anodes, Cnodes, Ynodes, Lnodes) {
  
  # 1. Baseline variables (W)
  baseline_block <- Wvars[Wvars %in% names(df)]
  
  # 2. Build time-block ordering: C_t, A_t, L_t..., Y_t
  time_blocks <- list()
  
  for (t in 1:Kstar) {
    Ct <- Cnodes[t]
    At <- Anodes[t]
    Yt <- Ynodes[t]
    
    # Lnodes for that time t = any that end with _t or __t
    Lt <- Lnodes[grepl(paste0("_", t, "$"), Lnodes) |
                   grepl(paste0("_", t, "__"), Lnodes)]
    
    time_blocks[[t]] <- c(Lt,  At, Yt, Ct)
  }
  
  ordered_cols <- c(baseline_block, unlist(time_blocks))
  
  # keep only columns that exist
  ordered_cols <- ordered_cols[ordered_cols %in% names(df)]
  
  # return reordered df
  df[, ordered_cols, drop = FALSE]
}

## === CALL THE REORDERING RIGHT BEFORE ltmleMSM() ===

wide_df_ltmle <- reorder_for_ltmle(
  df = wide_df_scaled,
  Kstar = K,
  Wvars = Wvars,
  Anodes = Anodes,
  Cnodes = Cnodes,
  Ynodes = Ynodes,
  Lnodes = Lnodes
)


remove_L1 <- c(
  "bs_bmi_1",
  "hb_pvr_calc_1",
  "reveal_lite_2_1",
  "functional_class__1_1",
  "functional_class__2_1",
  "functional_class__4_1"
)

Lnodes_filtered <- setdiff(Lnodes_ordered, remove_L1)



#------------(5) bulid Qform and gform ------------------------------------------------
remove_Wvars<-c( "incident_or_prevalent__UNK","incident_or_prevalent__prevalent")
Wvars_filtered <- setdiff(Wvars, remove_Wvars)
   

# 1. Extract the exact Q-node order from our data frame
Qnodes_in_data <- names(wide_df_ltmle)[names(wide_df_ltmle) %in% c(Lnodes_filtered, Ynodes)]

# 2. Initialize Qform with correct names and length
Qform <- setNames(rep(NA_character_, length(Qnodes_in_data)), Qnodes_in_data)

# 3. Lnodes get "~ 1" (our request)
for (L in Lnodes_filtered) {
  Qform[L] <- "Q.kplus1 ~ 1"
}

# 4. Build the Ynode formulas
L_for_t <- function(t) Lnodes_filtered[grepl(paste0("_", t, "$"), Lnodes_filtered)]

for (t in seq_along(Ynodes)) {
  yname <- Ynodes[t]

  # All L variables we want included for Y_t
  Lt_all <- c(
    if (t == 1) remove_L1 else character(0),
    L_for_t(t)
  )

  rhs <- c(paste0("A", t), Lt_all, Wvars_filtered)

  Qform[yname] <- paste0("Q.kplus1 ~ ", paste(rhs, collapse = " + "))
}

# # 5. Validation
stopifnot(identical(names(Qform), Qnodes_in_data))  # <- KEY 
stopifnot(!any(is.na(Qform)))


gform <- c()

for (t in 1:K) {

  rhs_base <- c(
    Wvars_filtered,
    if (t == 1) remove_L1 else NULL,     # include the removed L1 variables in A1, C1
    Lnodes_filtered[grepl(paste0("_", t, "$"), Lnodes_filtered)]
  )

  # A-model
  rhsA <- if (t == 1) rhs_base else c(rhs_base, paste0("A", t-1))
  
  gform[paste0("A", t)] <- paste0("A", t, " ~ ", paste(rhsA, collapse = " + "))

  # C-model (if needed)
  if (paste0("C", t) %in% Cnodes) {
    rhsC <- c(rhs_base, paste0("A", t))
  if (t<K)  gform[paste0("C", t)] <- paste0("C", t, " ~ ", paste(rhsC, collapse = " + "))
  }
}

#-----(6) create time horizons for counterfactual outcomes and summary measures to be used in MSM-----------------------------
# 2 regimes x 1 summary measure (indicator for "never") x |k_list| horizons
k_list <- seq(5, 20, 5)  ### 5, 10, 15, 20 years since mono/dual therapies
k_list <- k_list[k_list <= K]
final.Ynodes <- paste0("Y", k_list)


summary.measures <- array(NA_real_, dim = c(2, 5, length(k_list)))
summary.measures[1, 1, ] <- 1  # regime 1 (initiate)
summary.measures[2, 1, ] <- 0  # regime 2 (never)


num.regimes<-2
for (i in 1:num.regimes) {
  
  for (k in 1:length(final.Ynodes))
  {
    summary.measures[i, 2, k] <- as.numeric(k==1) #indicator for the years
    summary.measures[i, 3, k] <- as.numeric(k==2)  
    summary.measures[i, 4, k] <- as.numeric(k==3)
    summary.measures[i, 5, k] <- as.numeric(k==4)
  }
  
}

colnames(summary.measures) <- c("assigned_treatment", "year5","year10", "year15", "year20")






## Then call ltmleMSM using wide_df_ltmle:

## =========================================================
## 6) Fit ltmleMSM
##    - survivalOutcome=TRUE for one-time event Y  (monotone 0→1)
##    - gbounds to stabilize small/large propensity tails
##    - SL.library: start simple; can expand later
## =========================================================


# choose simple libraries due to treatment and outcome data sparsity
mylibrary2 <- list(Q=c("SL.mean","SL.glm", "SL.gam", "SL.bayesglm"),
                  g=c("SL.mean","SL.glm","SL.bayesglm"))




