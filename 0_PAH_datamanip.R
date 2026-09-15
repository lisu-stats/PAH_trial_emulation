library(dplyr)

load(file='baseline_data_v2.rdata')
load(file='followup_data.rdata')
load(file='RPH_followup_data_v3.rdata')
RPH_followup_data<-RPH_follow_up_data_V3

#-----------------------------(1). pick up the columns, format the visit dates ---------------------------------------------


# Remove row-number column
baseline <- subset(baseline_data_v2, select = setdiff(names(baseline_data_v2), "X"))
followup <- subset(followup_data, select = setdiff(names(followup_data), c("...1", "diagnosis_date", "time_from_diagnosis_to_visit")))

names(baseline)
names(followup)
names(RPH_followup_data)<-c("id","study_visit" , "visit_date" , "bs_bmi" ,"hb_pvr_calc","functional_class","reveal_lite_2")
names(RPH_followup_data)

## 1) Coerce to character first 
baseline$diagnosis_date            <- as.character(baseline$diagnosis_date)
followup$visit_date                <- as.character(followup$visit_date)
RPH_followup_data$visit_date       <- as.character(RPH_followup_data$visit_date)

## 2) Parse with explicit formats

# baseline: "01/12/2013" ->  DD/MM/YYYY (1 December 2013)
baseline$diagnosis_date <- as.Date(baseline$diagnosis_date, format = "%d/%m/%Y")

# followup: already ISO "YYYY-MM-DD"; as.Date() is fine directly
followup$visit_date     <- as.Date(followup$visit_date)  # or format = "%Y-%m-%d"

# RPH: "14/02/2019" -> DD/MM/YYYY
RPH_followup_data$visit_date <- as.Date(RPH_followup_data$visit_date, format = "%d/%m/%Y")


#-----

# A) Show first 10 parsed values
print(baseline$diagnosis_date[1:10])
print(followup$visit_date[1:10])
print(RPH_followup_data$visit_date[1:10])

# B) Count NAs introduced by parsing
cat("NAs in baseline$diagnosis_date:", sum(is.na(baseline$diagnosis_date)), "\n")
cat("NAs in followup$visit_date    :", sum(is.na(followup$visit_date)), "\n")
cat("NAs in RPH$visit_date         :", sum(is.na(RPH_followup_data$visit_date)), "\n")

# C) Sanity windows (years should look sensible, e.g., 1980–2035)
range_baseline <- range(baseline$diagnosis_date, na.rm = TRUE)
range_followup <- range(followup$visit_date,     na.rm = TRUE)
range_rph      <- range(RPH_followup_data$visit_date, na.rm = TRUE)
print(range_baseline); print(range_followup); print(range_rph)


#--------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------(2). impute missing baseline covariates and stack baseline and follow-up data ---------------------------------------------

# create a NA category for functional_class, add a missing data indicator for bs_bmi, hb_pvr_calc, reveal_lite_2
# and then impute the missing values in the original column with sex*age group medians






## =========================
## Baseline median imputation (simple, combined)
## =========================

## -- 1) SETTINGS: add any other baseline numeric vars we want imputed
vars_numeric <- c("bs_bmi", "hb_pvr_calc", "reveal_lite_2")
# e.g., vars_numeric <- c("bs_bmi","hb_pvr_calc","reveal_lite_2","sbp","dbp")

## -- 2) Build Sex × Age‑quartile groups once
baseline$sex <- toupper(trimws(as.character(baseline$sex)))

# robust quartiles for age_diagnosis
q <- unique(quantile(baseline$age_diagnosis, probs = c(0,.25,.5,.75,1), na.rm = TRUE))
if (length(q) < 5L) {
  # widen bounds slightly if ties collapse cutpoints
  amin <- min(baseline$age_diagnosis, na.rm = TRUE) - 1e-9
  amax <- max(baseline$age_diagnosis, na.rm = TRUE) + 1e-9
  q <- sort(unique(c(amin, q, amax)))
}
baseline$age_quartile <- cut(baseline$age_diagnosis, breaks = q,
                             include.lowest = TRUE,
                             labels = paste0("Q", seq_len(length(q) - 1)))

grp_key <- interaction(baseline$sex, baseline$age_quartile, drop = TRUE, lex.order = TRUE)

## -- 3) Helper function: median impute within groups (with overall fallback)
impute_by_group <- function(x, grp) {
  miss <- is.na(x)
  # median per group (NA if group has no observed values)
  med_by_grp <- tapply(x, grp, function(z) if (all(is.na(z))) NA_real_ else median(z, na.rm = TRUE))
  overall <- median(x, na.rm = TRUE)
  
  # fill by group
  if (is.factor(grp)) {
    for (g in levels(grp)) {
      idx <- miss & (grp == g)
      med <- med_by_grp[[g]]
      if (is.null(med) || is.na(med)) med <- overall
      x[idx] <- med
    }
  }
  # fallback when grp is NA (e.g., missing sex/age)
  idx_na_grp <- miss & is.na(grp)
  if (any(idx_na_grp)) x[idx_na_grp] <- overall
  
  x
}

## -- 4) Create missingness indicators + impute for all numeric vars
for (v in vars_numeric) {
  # (a) missingness indicator (pre‑imputation)
  baseline[[paste0(v, "_missing")]] <- as.integer(is.na(baseline[[v]]))
  # (b) impute
  x <- suppressWarnings(as.numeric(baseline[[v]]))
  baseline[[v]] <- impute_by_group(x, grp_key)
}

## -- 5) functional_class: indicator + median impute as ordinal 1..4 + set ref="3"
# (a) original missingness indicator first
baseline$functional_class_missing <- as.integer(is.na(baseline$functional_class))

# (b) coerce to numeric codes 1..4 
map_fc_to_num <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("I","1")]   <- "1"
  x[x %in% c("II","2")]  <- "2"
  x[x %in% c("III","3")] <- "3"
  x[x %in% c("IV","4")]  <- "4"
  suppressWarnings(as.numeric(x))
}
fc_num <- map_fc_to_num(baseline$functional_class)

# (c) median impute within Sex × Age‑quartile; then round/cap to {1,2,3,4}
fc_num <- impute_by_group(fc_num, grp_key)
fc_num <- round(fc_num)
fc_num[fc_num < 1] <- 1
fc_num[fc_num > 4] <- 4



# (d) put back as factor with reference "3"
baseline$functional_class <- factor(as.character(fc_num), levels = c("1","2","3","4"))
baseline$functional_class <- stats::relevel(baseline$functional_class, ref = "3")


stopifnot(all(na.omit(as.numeric(as.character(baseline$functional_class))) %in% 1:4))
table(baseline$functional_class, useNA = "ifany")  ## all impuated as 3



#---------------(3). stack baseline and followup data-------------------------------------------------------------------------------------


## 0) Clean inputs
if ("...1" %in% names(baseline)) baseline[["...1"]] <- NULL
if ("...1" %in% names(followup)) followup[["...1"]] <- NULL

baseline$id <- as.character(baseline$id)
followup$id <- as.character(followup$id)
RPH_followup_data$id <- as.character(RPH_followup_data$id)


## 1) Create baseline visit rows (visit 0)
baseline_visit <- data.frame(
  id = baseline$id,
  study_visit = 0L,
  visit_date = baseline$diagnosis_date,
  bs_bmi = baseline$bs_bmi,
  hb_pvr_calc = baseline$hb_pvr_calc,
  functional_class = baseline$functional_class,
  reveal_lite_2 = baseline$reveal_lite_2,
  stringsAsFactors = FALSE
)

# prevent factor mismatch
if ("functional_class" %in% names(baseline_visit)) {
  baseline_visit$functional_class <- as.character(baseline_visit$functional_class)
}
if ("functional_class" %in% names(followup)) {
  if (is.factor(followup$functional_class)==F) {
    followup$functional_class <- as.character(followup$functional_class)
  }
}

if ("functional_class" %in% names(RPH_followup_data)) {
  if (is.factor(RPH_followup_data$functional_class)==F) {
    RPH_followup_data$functional_class <- as.character(RPH_followup_data$functional_class)
  }
}


## 2) Align columns for rbind
cols <- union(names(baseline_visit), names(followup))

baseline_visit[setdiff(cols, names(baseline_visit))] <- NA
followup[setdiff(cols, names(followup))] <- NA
RPH_followup_data[setdiff(cols, names(RPH_followup_data))] <- NA

baseline_visit <- baseline_visit[, cols]
followup       <- followup[, cols]
RPH_followup_data <- RPH_followup_data[, cols]


## 3) Stack baseline + follow-up
visits_long <- rbind(baseline_visit, followup,RPH_followup_data)

## 4) Prepare static baseline variables
visit_vars <- intersect(names(baseline), names(followup))

# IMPORTANT UPDATE:
# static_vars now includes missing indicators such as:
# bs_bmi_missing, hb_pvr_calc_missing, reveal_lite_2_missing
static_vars <- setdiff(names(baseline), c("id", visit_vars))

baseline_static <- baseline[, c("id", static_vars), drop = FALSE]

if (any(duplicated(baseline_static$id))) {
  baseline_static <- baseline_static[!duplicated(baseline_static$id), , drop = FALSE]
}

## 5) Merge baseline static variables (including missing indicators)
long_analytic <- merge(visits_long, baseline_static, by = "id", all.x = TRUE)



## 7) Derived times
long_analytic$visit_date<-as.Date(long_analytic$visit_date)
long_analytic$diagnosis_date<-as.Date(long_analytic$diagnosis_date)

long_analytic$time_from_diagnosis_to_visit<-as.numeric((long_analytic$visit_date-long_analytic$diagnosis_date))/365.2425

numeric_cols <- c("time_from_diagnosis_to_visit",
                  "time_to_mono_or_dual_therapy",
                  "time_to_selexipag_from_diagnosis")

for (cc in numeric_cols) {
  if (cc %in% names(long_analytic)) {
    suppressWarnings(long_analytic[[cc]] <- as.numeric(long_analytic[[cc]]))
  }
}

long_analytic$time_from_backgroundT_to_visit <-
  long_analytic$time_from_diagnosis_to_visit - long_analytic$time_to_mono_or_dual_therapy

long_analytic$time_from_backgroundT_to_selexipag <-
  long_analytic$time_to_selexipag_from_diagnosis - long_analytic$time_to_mono_or_dual_therapy

long_analytic$time_from_backgroundT_to_death<-long_analytic$surv_time - long_analytic$time_to_mono_or_dual_therapy

long_analytic$status<-ifelse(long_analytic$sub_cause =='death', 1, 0)

names(long_analytic)

#----------------(4). create gap times between visits------------------------------------------------------------------------------------
# note the large gaps for prevalent cases!!!!!!!!!


# Ensure data are ordered correctly
long_analytic <- long_analytic[order(long_analytic$id,
                                     long_analytic$study_visit,
                                     long_analytic$visit_date), ]

# Create a new variable for gap times
long_analytic$gap_time <- NA

# Compute within each id
ids <- unique(long_analytic$id)

for (id_value in ids) {
  
  idx <- which(long_analytic$id == id_value)
  
  # Extract the vector of visit times for this id
  tvec <- long_analytic$time_from_backgroundT_to_visit[idx]
  
  # Compute diffs: gap between consecutive visits
  gaps <- c(NA, diff(tvec))
  
  # Assign back to the dataset
  long_analytic$gap_time[idx] <- gaps
}  

#----------------------------------------------------------------------------------------------------
#-----------------------------(5). create discrete-time (yearly) survival indicators---------------------------------------------
#----------------------------------------------------------------------------------------------------


############################################################
# SETTINGS
############################################################

# length of interval (in original time units; change to 1 if already in years)
year_len <- 1

# time-varying covariates to include
tv_vars <- c("bs_bmi", "hb_pvr_calc", "functional_class", "reveal_lite_2")

# baseline covariates (add/remove depending on  data)
baseline_vars <- c("sex","CCI","age_diagnosis","centre", "time_to_mono_or_dual_therapy", "diagnosis_verified" , "incident_or_prevalent", "vasoresponder","GRIPHON_eligibility"  )

############################################################
# PREPARE
############################################################

# ensure proper types
long_analytic$id <- as.character(long_analytic$id)
long_analytic$time_from_backgroundT_to_visit  <- as.numeric(long_analytic$time_from_backgroundT_to_visit)
long_analytic$time_from_backgroundT_to_death  <- as.numeric(long_analytic$time_from_backgroundT_to_death)
long_analytic$status <- as.numeric(long_analytic$status)

# maximum number of intervals in cohort
K <- ceiling(max(long_analytic$time_from_backgroundT_to_death, na.rm = TRUE) / year_len)

############################################################
# FUNCTION: MAKE WIDE ROW FOR ONE ID
############################################################

make_row <- function(pid) {
  
  di <- long_analytic[long_analytic$id == pid, , drop = FALSE]
  di <- di[order(di$time_from_backgroundT_to_visit), ]
  
  # survival time and event indicator
  Ti <- di$time_from_backgroundT_to_death[1]
  ev <- di$status[1]
  
  # censoring time is same as Ti but meaning depends on ev
  censor_time <- Ti
  
  # baseline values from visit 0
  base_row <- di[di$study_visit == 0, ]
  if (nrow(base_row) == 0) base_row <- di[1, ]  # fallback
  
  out <- list(id = pid)
  
  # add baseline
  for (bv in baseline_vars) out[[bv]] <- base_row[[bv]]
  
  # construct interval nodes
  for (t in 1:K) {
    
    t_start <- (t - 1) * year_len
    t_end   <- min(t * year_len, Ti)
    
    #### CENSORING NODE
    cens_status <- if (t_end < censor_time) "uncensored" else "censored"
    out[[paste0("C", t)]] <- factor(cens_status,
                                    levels = c("censored","uncensored"))
    
    #### TREATMENT NODE (A_t)
    t_selex <- di$time_from_backgroundT_to_selexipag[1]
    
    if (!is.na(t_selex) && t_end >= t_selex) {
      A_t <- 1L
    } else {
      A_t <- 0L
    }
    
    out[[paste0("A", t)]] <- A_t
    
    #### TIME-VARYING COVARIATES (L_t) — per-variable closest non-missing ≤ t_end, treatment initiation and event
    
    # For convenience
    t_visit <- di$time_from_backgroundT_to_visit
    
    for (v in tv_vars) {
      # Candidate rows: observed (non-missing) for this variable
      obs_mask <- !is.na(di[[v]]) & !is.na(t_visit)
      
      
      # cutoff for covariates: do not use post-treatment values
      cov_cutoff <- if (!is.na(t_selex)) min(t_end, t_selex) else t_end
      
      prior_mask <- obs_mask & (t_visit <= cov_cutoff)
      
      
      if (any(prior_mask)) {
        # Pick the observation with the largest time ≤ t_end (i.e., closest to interval end)
        pick_v <- which.max(t_visit[prior_mask])
        idx <- which(prior_mask)[pick_v]
        out[[paste0(v, "_", t)]] <- di[[v]][idx]
      } else if (any(obs_mask)) {
        # No prior non-missing: fallback to the earliest non-missing overall (e.g., baseline)
        pick_v <- which.min(t_visit[obs_mask])
        idx <- which(obs_mask)[pick_v]
        out[[paste0(v, "_", t)]] <- di[[v]][idx]
      } else {
        # No non-missing anywhere for this variable for this person
        out[[paste0(v, "_", t)]] <- NA
      }
    }
    #### OUTCOME NODE (Y_t)
    y_t <- if (ev == 1 && Ti <= t * year_len) 1L else 0L
    if (t > 1) y_t <- max(out[[paste0("Y", t-1)]], y_t)
    out[[paste0("Y", t)]] <- y_t
  }
  
  # ensure last Y is 1 for event
  if (ev == 1) out[[paste0("Y", K)]] <- 1L
  
  return(out)
}

############################################################
# APPLY TO ALL IDs
############################################################

ids <- unique(long_analytic$id)
wide_list <- lapply(ids, make_row)

# convert to data.frame
wide_df <- do.call(rbind.data.frame, wide_list)
row.names(wide_df) <- NULL






#----------------------------------------------------------------------------------------------------
#-----------------------------(6). inspect wide format data against raw data ---------------------------------------------
#----------------------------------------------------------------------------------------------------

############################################################
# 0. SETUP
############################################################

# IDs to inspect
ids_to_check <- c("OC15TS", "OC22DR", "OC097R", "OC04GT", "OC02E4", "RPH4", "RPH32", "RPH23", "OC0648")

# Time-varying covariates
tv_vars <- c("bs_bmi","hb_pvr_calc","functional_class","reveal_lite_2")

# Infer number of intervals K from the WIDE dataset
infer_K <- function(df) {
  ycols <- grep("^Y[0-9]+$", names(df), value = TRUE)
  max(as.integer(sub("^Y","", ycols)))
}

K <- infer_K(wide_df)

# Generate node names
Ynodes <- paste0("Y", 1:K)
Anodes <- paste0("A", 1:K)
Cnodes <- paste0("C", 1:K)

############################################################
# 1. RAW DATA EXTRACTOR
############################################################

get_raw_data <- function(id_value) {
  raw <- subset(long_analytic, id == id_value)
  raw <- raw[order(raw$time_from_backgroundT_to_visit, raw$study_visit), ]
  raw
}

############################################################
# 2. WIDE DATA EXTRACTOR
############################################################

get_wide_row <- function(id_value) {
  wide_df[wide_df$id == id_value, , drop = FALSE]
}

############################################################
# 3. UPDATED INTERVAL COMPARISON (variable-specific LOCF)
############################################################

interval_compare <- function(raw, wide, K, tv_vars, year_len = 1) {
  
  df <- data.frame(
    interval = 1:K,
    t_start  = (0:(K-1)) * year_len,
    t_end    = (1:K) * year_len
  )
  
  # Add wide C,A,Y
  for (t in 1:K) {
    df[[paste0("C",t)]] <- as.character(wide[[paste0("C",t)]])
    df[[paste0("A",t)]] <- wide[[paste0("A",t)]]
    df[[paste0("Y",t)]] <- wide[[paste0("Y",t)]]
  }
  
  # Prepare raw and wide Lnodes columns
  for (v in tv_vars) {
    df[[paste0("RAW_", v)]]  <- NA
    df[[paste0("WIDE_", v)]] <- NA
  }
  
  # Loop through intervals
  for (i in 1:K) {
    
    tend <- df$t_end[i]
    t_visit <- raw$time_from_backgroundT_to_visit
    
    for (v in tv_vars) {
      
      obs_mask   <- !is.na(raw[[v]]) & !is.na(t_visit)
      prior_mask <- obs_mask & (t_visit <= tend)
      
      if (any(prior_mask)) {
        # Closest non-missing BEFORE interval end
        pick_v <- which.max(t_visit[prior_mask])
        idx <- which(prior_mask)[pick_v]
        raw_val <- raw[[v]][idx]
      } else if (any(obs_mask)) {
        # Fallback: earliest non-missing
        pick_v <- which.min(t_visit[obs_mask])
        idx <- which(obs_mask)[pick_v]
        raw_val <- raw[[v]][idx]
      } else {
        raw_val <- NA
      }
      
      # Wide value for interval i
      wide_val <- wide[[paste0(v, "_", i)]]
      
      df[[paste0("RAW_", v)]][i]  <- raw_val
      df[[paste0("WIDE_", v)]][i] <- wide_val
    }
  }
  
  df
}

############################################################
# 4. STRUCTURAL CHECKS
############################################################

check_wide_row <- function(row, K) {
  
  y <- as.numeric(row[Ynodes])
  a <- as.numeric(row[Anodes])
  c <- as.character(row[Cnodes])
  
  list(
    Y_monotone = all(diff(y) >= 0, na.rm = TRUE),
    
    A_switch_once = {
      d <- diff(a)
      all(d >= 0, na.rm = TRUE) && sum(d > 0, na.rm = TRUE) <= 1
    },
    
    C_stays_censored = {
      first_c <- which(c == "censored")[1]
      is.na(first_c) || all(c[first_c:length(c)] == "censored")
    }
  )
}

############################################################
# 5. MASTER PATIENT CHECKER FUNCTION
############################################################

check_patient <- function(id_value) {
  
  cat("\n", strrep("=", 70), "\nID:", id_value,
      "\n", strrep("-", 70), "\n\n", sep = "")
  
  raw  <- get_raw_data(id_value)
  wide <- get_wide_row(id_value)
  
  if (nrow(raw)==0) { cat("No raw data.\n"); return(NULL) }
  if (nrow(wide)==0){ cat("No wide data.\n"); return(NULL) }
  
  cat("RAW DATA (baseline + follow-up):\n")
  print(raw, row.names = FALSE, max = 500)
  
  cat("\nWIDE DATA:\n")
  print(wide, row.names = FALSE, max = 500)
  
  cat("\nINTERVAL-BY-INTERVAL COMPARISON (variable-specific LOCF):\n")
  comp <- interval_compare(raw, wide, K, tv_vars)
  print(comp, row.names = FALSE, max = 2000)
  
  cat("\nSTRUCTURAL CHECKS:\n")
  print(check_wide_row(wide, K))
  
  invisible(list(raw = raw, wide = wide, comparison = comp))
}

############################################################
# 6. RUN FOR ALL FIVE PATIENTS
############################################################
sink("data_inspection_September_15_2026.txt")

results <- lapply(ids_to_check, check_patient)
names(results) <- ids_to_check

sink()
#----------------------------------------------------------------------------------------------------
#-----------------------------(7). standardise non-factor variables ---------------------------------------------
#----------------------------------------------------------------------------------------------------

## ---- 1) Infer K and define node sets ----
infer_K <- function(df) {
  ycols <- grep("^Y[0-9]+$", names(df), value = TRUE)
  if (length(ycols) == 0) stop("No Y nodes found in wide_df.")
  max(as.integer(sub("^Y", "", ycols)))
}

K <- infer_K(wide_df)
Anodes <- paste0("A", 1:K)
Ynodes <- paste0("Y", 1:K)
Cnodes <- paste0("C", 1:K)   # usually factors ("censored"/"uncensored")

## ---- 2) Identify columns to scale ----
# Time-varying L-nodes typically have pattern "<var>_<t>"
# We’ll pick any column ending with "_<integer>" that is numeric and not A/Y/C nodes.
ends_with_time <- grepl("_[0-9]+$", names(wide_df))
is_num  <- sapply(wide_df, is.numeric)
is_fact <- sapply(wide_df, is.factor)

# Baseline numeric variables we want scaled
baseline_to_scale <- intersect(c("CCI", "age_diagnosis", "time_to_mono_or_dual_therapy"),
                               names(wide_df))

# Exclusions
is_id     <- names(wide_df) %in% "id"
is_A      <- names(wide_df) %in% Anodes
is_Y      <- names(wide_df) %in% Ynodes
is_C      <- names(wide_df) %in% Cnodes         # Cnodes are factors; will be excluded by !is_num too
is_meta   <- grepl("^\\.", names(wide_df))      # if we kept .Ti, .status etc.
is_flag   <- grepl("_missing$", names(wide_df)) # keep 0/1 flags unscaled

# Time-varying numeric L-node candidates:
L_numeric <- ends_with_time & is_num & !(is_A | is_Y | is_C)

# Final set to scale = (L_numeric) OR (baseline_to_scale), excluding id/A/Y/meta/flags/factors
to_scale <- (L_numeric | (names(wide_df) %in% baseline_to_scale)) &
  !is_id & !is_A & !is_Y & !is_meta & !is_flag & !is_fact

## ---- 3) Safe z-scoring helper ----
zscale <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(x - m)   # avoid division by zero; centers only
  (x - m) / s
}

## ---- 4) Make a scaled copy ----
wide_df_scaled <- wide_df
wide_df_scaled[to_scale] <- lapply(wide_df_scaled[to_scale], zscale)

## ---- 5) Quick sanity checks (optional) ----
# Treatment and outcome untouched

# A nodes: must be {0,1,NA}
avals <- suppressWarnings(as.numeric(unlist(wide_df_scaled[Anodes[Anodes %in% names(wide_df_scaled)]], use.names = FALSE)))
stopifnot(all(is.na(avals) | avals %in% c(0,1)))

# Y nodes: must be {0,1,NA}
yvals <- suppressWarnings(as.numeric(unlist(wide_df_scaled[Ynodes[Ynodes %in% names(wide_df_scaled)]], use.names = FALSE)))
stopifnot(all(is.na(yvals) | yvals %in% c(0,1)))




# Peek at means/sds for a few scaled columns
check_cols <- c(intersect(baseline_to_scale, names(wide_df_scaled)),
                names(wide_df_scaled)[which(L_numeric)[1:min(126, sum(L_numeric))]])
summary_stats <- t(vapply(wide_df_scaled[check_cols], function(x) c(mean=mean(x, na.rm=TRUE),
                                                                    sd=sd(x, na.rm=TRUE)),
                          numeric(2)))
print(round(summary_stats, 3))


# 3. Create clean FC dummies

# =========================================================
# REMOVE ORIGINAL FUNCTIONAL CLASS COLUMNS 
# =========================================================
orig_fc <- grep("^functional_class_[0-9]+$", names(wide_df_scaled), value = TRUE)

if (length(orig_fc) > 0) {
  wide_df_scaled <- wide_df_scaled[, !names(wide_df_scaled) %in% orig_fc]
}


fc_dummy <- c()

for (t in 1:K) {
  
  # extract from unscaled wide_df to avoid interference from scaling
  fc_raw_name <- paste0("functional_class_", t)
  if (!(fc_raw_name %in% names(wide_df))) next
  
  fc <- factor(as.character(wide_df[[fc_raw_name]]))
  
  # set "3" as reference
  if ("3" %in% levels(fc)) {
    fc <- stats::relevel(fc, ref = "3")
  }
  
  # non-reference levels
  other_lvls <- setdiff(levels(fc), "3")
  
  for (lvl in other_lvls) {
    # SAFE naming: functional_class__<level>_<time>
    newcol <- paste0("functional_class__", lvl, "_", t)
    
    wide_df_scaled[[newcol]] <- as.integer(fc == lvl)
    fc_dummy <- c(fc_dummy, newcol)
  }
}


#####create indicator for Original cohort or RPH

wide_df_scaled$id_is_OC <- as.numeric(grepl("^(OC)", wide_df_scaled$id))

wide_df_scaled$RPH<-ifelse(wide_df_scaled$centre=='Papworth',1,0)

names(wide_df_scaled)


############################################################
# DONE — wide_df is READY for LTMLE
############################################################

# We will later define:
#   Anodes <- paste0("A", 1:K)
#   Cnodes <- paste0("C", 1:K)
#   Ynodes <- paste0("Y", 1:K)
#   Lnodes <- as.vector(t(outer(tv_vars, 1:K, paste, sep="_")))
#
# and then call ltmle().




