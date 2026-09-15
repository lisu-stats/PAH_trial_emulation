
###############################################################################
# E-value calculations directly from the RD_summary object
# RD_summary must contain:
#   horizon, risk_never, RD, RD_LO, RD_HI
###############################################################################

# ---- 1. Take the RD_summary object exactly as is ----
load(file='../results/Group20/group_jackknife_msm_main.rdata')


df <- get.summary.main(GJK.main)

# ---- 2. Compute treated risks and risk ratios ----

df$p1     <- df$risk_never + df$RD
df$p1_lo  <- df$risk_never + df$RD_LO       # lower RD = stronger effect
df$p1_hi  <- df$risk_never + df$RD_HI       # upper RD = closer to null

df$RR     <- df$p1    / df$risk_never
df$RR_lo  <- df$p1_lo / df$risk_never
df$RR_hi  <- df$p1_hi / df$risk_never       # this RR_hi is used for E-value CI

# ---- 3. E-value function for protective effects (RR < 1) ----

Evalue <- function(RR) {
  inv <- 1 / RR
  inv + sqrt(inv * (inv - 1))
}

# ---- 4. Apply E-value formulas ----

df$E_point <- Evalue(df$RR)
df$E_CI    <- Evalue(df$RR_hi)   # CI bound closest to null

# ---- 5. Print a clean summary table ----

out <- df[, c("horizon", "RR", "E_point", "RR_hi", "E_CI")]
round(out, 3)






