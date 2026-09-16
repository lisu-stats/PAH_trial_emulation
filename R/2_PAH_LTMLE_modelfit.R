# main model

working.msm<-"Y ~ -1 + year5  + year10 + year15  + year20 + assigned_treatment"


set.seed(12345)
fit_msm <- ltmleMSM(
  data             = wide_df_ltmle,
  Anodes           = Anodes,
  Cnodes           = Cnodes,
  Lnodes           = Lnodes_filtered,
  Ynodes           = Ynodes,
  final.Ynodes     = final.Ynodes,
  survivalOutcome  = TRUE,
  regimes          = regimes,
  summary.measures = summary.measures,
  working.msm      = working.msm,
  Qform            = Qform,
  gform            = gform,
  deterministic.g.function = MaintainTreatment,
  SL.library       = mylibrary2,
  gbounds          = c(0.01, 1),
  variance.method = "ic"
)

## =========================================================
## Results
## =========================================================
sm <- summary(fit_msm)
print(sm)




#treatment effects by year

working.msm<-"Y ~ -1 + year5 + assigned_treatment:year5 + year10 + assigned_treatment:year10 + year15 + assigned_treatment:year15 + year20 + assigned_treatment:year20"


set.seed(12345)
fit_msm2 <- ltmleMSM(
  data             = wide_df_ltmle,
  Anodes           = Anodes,
  Cnodes           = Cnodes,
  Lnodes           = Lnodes_filtered,
  Ynodes           = Ynodes,
  final.Ynodes     = final.Ynodes,
  survivalOutcome  = TRUE,
  regimes          = regimes,
  summary.measures = summary.measures,
  working.msm      = working.msm,
  Qform            = Qform,
  gform            = gform,
  deterministic.g.function = MaintainTreatment,
  SL.library       = mylibrary2,
  gbounds          = c(0.01, 1),
  variance.method = "ic"
)

sm2 <- summary(fit_msm2)
print(sm2)





#treatment effects by years, stratified by GRIPHON_eligibility

#working.msm<-"Y ~ -1 + year5 + assigned_treatment:year5 + year10 + assigned_treatment:year10 + year15 + assigned_treatment:year15 + year20 + assigned_treatment:year20 + year5:GRIPHON_eligibility__yes +year10:GRIPHON_eligibility__yes+ year15:GRIPHON_eligibility__yes + year20:GRIPHON_eligibility__yes + year5:assigned_treatment:GRIPHON_eligibility__yes +year10:assigned_treatment:GRIPHON_eligibility__yes+ year15:assigned_treatment:GRIPHON_eligibility__yes + year20:assigned_treatment:GRIPHON_eligibility__yes"

working.msm<-"Y ~ -1 + year5 + year10 +  year15 + year20 + assigned_treatment + year5:GRIPHON_eligibility__yes +year10:GRIPHON_eligibility__yes+ year15:GRIPHON_eligibility__yes + year20:GRIPHON_eligibility__yes + assigned_treatment:GRIPHON_eligibility__yes "



set.seed(12345)
fit_msm3 <- ltmleMSM(
  data             = wide_df_ltmle,
  Anodes           = Anodes,
  Cnodes           = Cnodes,
  Lnodes           = Lnodes_filtered,
  Ynodes           = Ynodes,
  final.Ynodes     = final.Ynodes,
  survivalOutcome  = TRUE,
  regimes          = regimes,
  summary.measures = summary.measures,
  working.msm      = working.msm,
  Qform            = Qform,
  gform            = gform,
  deterministic.g.function = MaintainTreatment,
  SL.library       = mylibrary2,
  gbounds          = c(0.01, 1),
  variance.method = "ic"
)


sm3 <- summary(fit_msm3)
print(sm3)



#treatment effects by years, stratified by prevalent or incident case

working.msm<-"Y ~ -1 + year5 + year10  + year15  + year20 + year5:incident_or_prevalent__prevalent +year10:incident_or_prevalent__prevalent+ year15:incident_or_prevalent__prevalent + year20:incident_or_prevalent__prevalent + assigned_treatment +  assigned_treatment:incident_or_prevalent__prevalent"

set.seed(12345)
fit_msm4 <- ltmleMSM(
  data             = wide_df_ltmle,
  Anodes           = Anodes,
  Cnodes           = Cnodes,
  Lnodes           = Lnodes_filtered,
  Ynodes           = Ynodes,
  final.Ynodes     = final.Ynodes,
  survivalOutcome  = TRUE,
  regimes          = regimes,
  summary.measures = summary.measures,
  working.msm      = working.msm,
  Qform            = Qform,
  gform            = gform,
  deterministic.g.function = MaintainTreatment,
  SL.library       = mylibrary2,
  gbounds          = c(0.01, 1),
  variance.method = "ic"
)


sm4 <- summary(fit_msm4)
print(sm4)


#treatment effects by   RPH or non-RPH------------------------------------------------

working.msm<-"Y ~ -1 + year5 + year10 +  year15 + year20 + assigned_treatment + year5:RPH +year10:RPH+ year15:RPH + year20:RPH + assigned_treatment:RPH "

set.seed(12345)
fit_msm5 <- ltmleMSM(
  data             = wide_df_ltmle,
  Anodes           = Anodes,
  Cnodes           = Cnodes,
  Lnodes           = Lnodes_filtered,
  Ynodes           = Ynodes,
  final.Ynodes     = final.Ynodes,
  survivalOutcome  = TRUE,
  regimes          = regimes,
  summary.measures = summary.measures,
  working.msm      = working.msm,
  Qform            = Qform,
  gform            = gform,
  deterministic.g.function = MaintainTreatment,
  SL.library       = mylibrary2,
  gbounds          = c(0.01, 1),
  variance.method = "ic"
)


sm5 <- summary(fit_msm5)
print(sm5)
