####clean dataset----
incident_CVD <- data_baseline
dput(names(incident_CVD))
incident_CVD_1 <- subset(incident_CVD, select = c("CALC_all", "IPH_all", "LRNC_all", "Max_IMT_complete", "Stenosis_complete", "ergoid"))

incident_CVD_sup <- data_secondscan_2
incident_CVD_sup <- subset(incident_CVD_sup, select = c("ergoid", "ergocar_date", "ergocarfup_date", "First_age", "Sex"))
incident_CVD_sup <- incident_CVD_sup %>%
  distinct(ergoid, .keep_all = TRUE)

incident_CVD_2 <- left_join(incident_CVD_1, incident_CVD_sup, by = "ergoid")
incident_CVD_2 <- left_join(incident_CVD_2, data_PA, by = "ergoid")
incident_CVD_2 <- subset(incident_CVD_2, !is.na(total_MET))

summary(incident_CVD_2)

#######cvd outcome
Stroke <- read_sav('Stroke_data_Luo_Oct2024.sav')
names(Stroke)
Stroke_0 <- subset(Stroke, select = c("ergoid","IC", "prevalent_stroke","Baseline_date","incident_stroke","stroke_date","stroke_subtype" ))
Stroke_1 <- rename(Stroke_0, ic_Stroke = IC, prev_stroke_date = Baseline_date, prev_Stroke = prevalent_stroke, inc_stroke = incident_stroke, inc_stroke_date = stroke_date)
Stroke_1 <- subset(Stroke_1, ic_Stroke == 1)

incident_CVD_3 <- left_join(incident_CVD_2, Stroke_1[c("ergoid", "prev_Stroke", "inc_stroke", "inc_stroke_date","stroke_subtype")], by = "ergoid")
table(incident_CVD_3$prev_Stroke)

incident_CVD_5 <- subset(incident_CVD_3, prev_Stroke == 0)
#end date end day of 2020
incident_CVD_5$inc_stroke_date[is.na(incident_CVD_5$inc_stroke_date)] <- as.Date("2020-12-31")

#mortality
death <- read_sav('fp_VitalStatus_2022-42.sav')
death_0 <- subset(death, select = c("ergoid", "fp_vitalstatus", "fp_censordate"))
death_1 <- rename(death_0, death = fp_vitalstatus, death_date = fp_censordate)
death_1$death[is.na(death_1$death)] <- 0
incident_CVD_5 <- left_join(incident_CVD_5, death_1, by = "ergoid")
incident_CVD_5$death_check <- incident_CVD_5$death
table(incident_CVD_5$death_check)
incident_CVD_5$death_check[incident_CVD_5$death_date > as.Date("2020-12-31")] <- 0
incident_CVD_5$death_time <- time_length(difftime(incident_CVD_5$death_date, incident_CVD_5$ergocar_date), "years")
incident_CVD_5$death_check[incident_CVD_5$death_date > incident_CVD_5$inc_stroke_date] <- 0
incident_CVD_5$inc_stroke_date[
  (incident_CVD_5$inc_stroke_date > incident_CVD_5$death_date)&
    (incident_CVD_5$death_check == 1)
] <- incident_CVD_5$death_date[
  (incident_CVD_5$inc_stroke_date > incident_CVD_5$death_date)&
    (incident_CVD_5$death_check == 1)
]
incident_CVD_5$stroke_time <- time_length(difftime(incident_CVD_5$inc_stroke_date, incident_CVD_5$ergocar_date), "years")
table(incident_CVD_5$inc_stroke)
summary(incident_CVD_5)

incident_CVD_5 <- mutate(incident_CVD_5, 
                         stat = case_when(inc_stroke == 1 ~ 1,
                                          inc_stroke == 0 & death_check == 1 ~ 2,
                                          inc_stroke == 0 & death_check == 0 ~ 0))

table(incident_CVD_5$stat)
summary(incident_CVD_5$stroke_time)


dput(names(incident_CVD_5))
incident_CVD_check <- incident_CVD_5[,c("prev_Stroke", "inc_stroke", "inc_stroke_date","death", 
                                        "death_date", "death_check", "death_time", 
                                        "stroke_time", "stat")]
incident_CVD_5 <- left_join(incident_CVD_5, income_edu, by="ergoid")
incident_CVD_5 <- left_join(incident_CVD_5, energy, by="ergoid")
incident_CVD_5 <- left_join(incident_CVD_5, fhis, by = "ergoid")#fhis, familiy history of CVD

#协变量填补
Scan_covar_use <- subset(Scan_covar, select = -c(Sex,First_age))
incident_CVD_5 <- left_join(incident_CVD_5, Scan_covar_use, by = "ergoid")
summary(incident_CVD_5)

dput(names(incident_CVD_5))
imp0 <- mice(incident_CVD_5, maxit=0, defaultMethod = c("pmm","logreg","polr","polyreg"))# pmm for any variable, logreg for binary, polr for ordered, polyreg for unordered
meth <- imp0$method
meth[c("ergoid", "ergocar_date", "ergocarfup_date", "METh_modvig", "income", 
       "rs_cohort", "alcohol", "CRP", "ethnicity", "stroke_time", "CVD_his", "stroke_subtype")] <- "" #variables not being imputed will also not be used as predictors
# predictor matrix
pred <- imp0$predictorMatrix
pred[ ,c("ergoid", "ergocar_date", "ergocarfup_date", "METh_modvig", "income",
         "rs_cohort", "alcohol", "CRP", "ethnicity", "stroke_time", "CVD_his", "stroke_subtype")] <- 0 #will not be used as predictor in any model
incident_CVD_5_0 <- mice(incident_CVD_5, m=30, maxit=30, meth=meth, seed=619, predictorMatrix = pred)
incident_CVD_5_long <- complete(incident_CVD_5_0, action = 'long', include = T)
incident_CVD_5_long_imputed30 <- complete(incident_CVD_5_0, action = 'long')
summary(incident_CVD_5_long_imputed30)

#cleaning the covariates
incident_CVD_5_long$rs_cohort <- factor(incident_CVD_5_long$rs_cohort)
#binary smoking
incident_CVD_5_long <- mutate(incident_CVD_5_long,
                              smokingstatus_2 = case_when(
                                smokingstatus == "Never_smoke" ~ 0,
                                smokingstatus == "Past_smoke" ~ 1, 
                                smokingstatus == "Current_smoke" ~ 1))
incident_CVD_5_long$smokingstatus_2 <- factor(incident_CVD_5_long$smokingstatus_2, levels = c(0,1), labels = c("Never_Past", "Current_smoke"))
#binary education
incident_CVD_5_long <- mutate(incident_CVD_5_long,
                              edu_2 = case_when(
                                edu == "primary education"  | edu == "lower/intermediate general education OR lower vocational education"   | edu == "intermediate vocational education OR higher general education" ~ 0,
                                edu == "higher vocational education OR university" ~ 1))
incident_CVD_5_long$edu_2 <- factor(incident_CVD_5_long$edu_2, levels = c(0,1), labels = c("low", "high"))
#hypertension
incident_CVD_5_long <- mutate(incident_CVD_5_long, hypertension_1 = case_when(SBP >= 140 | DBP >= 90 | antihypertensive == "User" ~ 1,TRUE ~ 0))
incident_CVD_5_long$hypertension_1 <- factor(incident_CVD_5_long$hypertension_1, labels = c("normal","hypertension"))
#lipid lowering
incident_CVD_5_long <- mutate(incident_CVD_5_long, hyperchol = case_when(TC >= 6.2 | lipidreducing == "User" ~ 1,TRUE ~ 0))
incident_CVD_5_long$hyperchol <- factor(incident_CVD_5_long$hyperchol, labels = c("normal","Hypercholesterolaemia"))

incident_CVD_5_long$IPH_all <- factor(incident_CVD_5_long$IPH_all, labels = c("no_IPH", "IPH"))
incident_CVD_5_long$LRNC_all <- factor(incident_CVD_5_long$LRNC_all, labels = c("no_LRNC", "LRNC"))

incident_CVD_5_long <- mutate(incident_CVD_5_long, one_component_all = case_when(
  IPH_all == "IPH" | LRNC_all == "LRNC" ~ 1,
  TRUE ~ 0
))
incident_CVD_5_long$one_component_all <- factor(incident_CVD_5_long$one_component_all, labels = c("no_one_component", "one_component"))
data_CVD_nonlinear <- subset(incident_CVD_5_long, .imp == 1)
####stroke case number----
table(data_CVD_nonlinear$stat)
table(data_CVD_nonlinear$stat, data_CVD_nonlinear$one_component_all)


incident_CVD_5_long_IPH <- subset(incident_CVD_5_long, IPH_all == "IPH")
incident_CVD_5_long_IPH <- as.mids(incident_CVD_5_long_IPH)
incident_CVD_5_long_no_IPH <- subset(incident_CVD_5_long, IPH_all == "no_IPH")
incident_CVD_5_long_no_IPH <- as.mids(incident_CVD_5_long_no_IPH)
incident_CVD_5_long_LRNC <- subset(incident_CVD_5_long, LRNC_all == "LRNC")
incident_CVD_5_long_LRNC <- as.mids(incident_CVD_5_long_LRNC)
incident_CVD_5_long_no_LRNC <- subset(incident_CVD_5_long, LRNC_all == "no_LRNC")
incident_CVD_5_long_no_LRNC <- as.mids(incident_CVD_5_long_no_LRNC)
incident_CVD_5_long_one_component <- subset(incident_CVD_5_long, one_component_all == "one_component")
incident_CVD_5_long_one_component <- as.mids(incident_CVD_5_long_one_component)
incident_CVD_5_long_no_one_component <- subset(incident_CVD_5_long, one_component_all == "no_one_component")
incident_CVD_5_long_no_one_component <- as.mids(incident_CVD_5_long_no_one_component)
incident_CVD_5_long <- as.mids(incident_CVD_5_long)
incident_CVD_5_check <- complete(incident_CVD_5_long, action = 'long')

####function to get the output
format_pooled_results <- function(pooled_model){
  fit_summary <- pooled_model
  fit_summary$low <- fit_summary$estimate - 1.96 * fit_summary$std.error
  fit_summary$high <- fit_summary$estimate + 1.96 * fit_summary$std.error
  fit_subset <- fit_summary[, c("term", "estimate", "low", "high")]
  fit_subset$estimate <- round(exp(fit_subset$estimate), 2)
  fit_subset$low <- round(exp(fit_subset$low), 2)
  fit_subset$high <- round(exp(fit_subset$high), 2)
  formatted_result <- paste0(fit_subset$estimate, " (", fit_subset$low, ", ", fit_subset$high, ")")
  return(formatted_result)
}

##############################one_component-------------------------------
#####total PA----
#####stroke
#nonliear
rg1 <- coxph(Surv(stroke_time, stat == 1) ~
               total_MET + one_component_all +
               Sex + First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
               hyperchol + hypertension_1 + dm + Max_IMT_complete, data = data_CVD_nonlinear)
rg2 <- coxph(Surv(stroke_time, stat == 1) ~
               ns(total_MET, 3) + one_component_all +
               Sex + First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
               hyperchol + hypertension_1 + dm + Max_IMT_complete, data = data_CVD_nonlinear)
print(anova(rg1, rg2), digits = 4)
#full cohort
rg1 <- with(incident_CVD_5_long,
            coxph(Surv(stroke_time, stat == 1) ~
                    I(total_MET/20) +
                    Sex + First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
                    hyperchol + hypertension_1 + dm + Max_IMT_complete))
summary(pool(rg1))
full_cohort <- format_pooled_results(summary(pool(rg1))[1,])
full_cohort
#no_one_component
rg1 <- with(incident_CVD_5_long_no_one_component,
            coxph(Surv(stroke_time, stat == 1) ~
                    I(total_MET/20) +
                    Sex + First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
                    hyperchol + hypertension_1 + dm + Max_IMT_complete))
summary(pool(rg1))
no_one_component <- format_pooled_results(summary(pool(rg1))[1,])
no_one_component
#one_component
rg1 <- with(incident_CVD_5_long_one_component,
            coxph(Surv(stroke_time, stat == 1) ~
                    I(total_MET/20) +
                    Sex + First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
                    hyperchol + hypertension_1 + dm + Max_IMT_complete))
summary(pool(rg1))
one_component <- format_pooled_results(summary(pool(rg1))[1,])
one_component
#interaction
rg1 <- with(incident_CVD_5_long,
            coxph(Surv(stroke_time, stat == 1) ~
                    I(total_MET/20)*one_component_all +
                    Sex + First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
                    hyperchol + hypertension_1 + dm + Max_IMT_complete))
summary(pool(rg1))
fit_rg0 <- summary(pool(rg1))[14,]
p_interaction <- c(round(fit_rg0$p.value, digits = 4))
p_interaction
total_stroke <- data.frame(full_cohort, no_one_component, one_component, p_interaction)
total_stroke

#####modvig PA----
#####stroke
#nonliear
rg1 <- coxph(Surv(stroke_time, stat == 1) ~
               METh_modvig + one_component_all +
               Sex + First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
               hyperchol + hypertension_1 + dm + Max_IMT_complete, data = data_CVD_nonlinear)
rg2 <- coxph(Surv(stroke_time, stat == 1) ~
               ns(METh_modvig, 3) + one_component_all +
               Sex + First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
               hyperchol + hypertension_1 + dm + Max_IMT_complete, data = data_CVD_nonlinear)
print(anova(rg1, rg2), digits = 4)
#full cohort
rg1 <- with(incident_CVD_5_long,
            coxph(Surv(stroke_time, stat == 1) ~
                    I(METh_modvig/20) +
                    Sex + First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
                    hyperchol + hypertension_1 + dm + Max_IMT_complete))
summary(pool(rg1))
full_cohort <- format_pooled_results(summary(pool(rg1))[1,])
full_cohort
#no_one_component
rg1 <- with(incident_CVD_5_long_no_one_component,
            coxph(Surv(stroke_time, stat == 1) ~
                    I(METh_modvig/20) +
                    Sex + First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
                    hyperchol + hypertension_1 + dm + Max_IMT_complete))
summary(pool(rg1))
no_one_component <- format_pooled_results(summary(pool(rg1))[1,])
no_one_component
#one_component
rg1 <- with(incident_CVD_5_long_one_component,
            coxph(Surv(stroke_time, stat == 1) ~
                    I(METh_modvig/20) +
                    Sex + First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
                    hyperchol + hypertension_1 + dm + Max_IMT_complete))
summary(pool(rg1))
one_component <- format_pooled_results(summary(pool(rg1))[1,])
one_component
#interaction
rg1 <- with(incident_CVD_5_long,
            coxph(Surv(stroke_time, stat == 1) ~
                    I(METh_modvig/20)*one_component_all +
                    Sex + First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
                    hyperchol + hypertension_1 + dm + Max_IMT_complete))
summary(pool(rg1))
fit_rg0 <- summary(pool(rg1))[14,]
p_interaction <- c(round(fit_rg0$p.value, digits = 4))
p_interaction
modvig_stroke <- data.frame(full_cohort, no_one_component, one_component, p_interaction)
modvig_stroke


#####moderate PA----
#####stroke
#nonliear
rg1 <- coxph(Surv(stroke_time, stat == 1) ~
               METh_moderate + one_component_all +
               Sex + First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
               hyperchol + hypertension_1 + dm + Max_IMT_complete, data = data_CVD_nonlinear)
rg2 <- coxph(Surv(stroke_time, stat == 1) ~
               ns(METh_moderate, 3) + one_component_all +
               Sex + First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
               hyperchol + hypertension_1 + dm + Max_IMT_complete, data = data_CVD_nonlinear)
print(anova(rg1, rg2), digits = 4)
#full cohort
rg1 <- with(incident_CVD_5_long,
            coxph(Surv(stroke_time, stat == 1) ~
                    I(METh_moderate/20) +
                    Sex + First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
                    hyperchol + hypertension_1 + dm + Max_IMT_complete))
summary(pool(rg1))
full_cohort <- format_pooled_results(summary(pool(rg1))[1,])
full_cohort
#no_one_component
rg1 <- with(incident_CVD_5_long_no_one_component,
            coxph(Surv(stroke_time, stat == 1) ~
                    I(METh_moderate/20) +
                    Sex + First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
                    hyperchol + hypertension_1 + dm + Max_IMT_complete))
summary(pool(rg1))
no_one_component <- format_pooled_results(summary(pool(rg1))[1,])
no_one_component
#one_component
rg1 <- with(incident_CVD_5_long_one_component,
            coxph(Surv(stroke_time, stat == 1) ~
                    I(METh_moderate/20) +
                    Sex + First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
                    hyperchol + hypertension_1 + dm + Max_IMT_complete))
summary(pool(rg1))
one_component <- format_pooled_results(summary(pool(rg1))[1,])
one_component
#interaction
rg1 <- with(incident_CVD_5_long,
            coxph(Surv(stroke_time, stat == 1) ~
                    I(METh_moderate/20)*one_component_all +
                    Sex + First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
                    hyperchol + hypertension_1 + dm + Max_IMT_complete))
summary(pool(rg1))
fit_rg0 <- summary(pool(rg1))[14,]
p_interaction <- c(round(fit_rg0$p.value, digits = 4))
p_interaction
moderate_stroke <- data.frame(full_cohort, no_one_component, one_component, p_interaction)
moderate_stroke


#####vigorous PA----
#####stroke
#nonliear
rg1 <- coxph(Surv(stroke_time, stat == 1) ~
               METh_vigorous + one_component_all +
               Sex + First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
               hyperchol + hypertension_1 + dm + Max_IMT_complete, data = data_CVD_nonlinear)
rg2 <- coxph(Surv(stroke_time, stat == 1) ~
               ns(METh_vigorous, knots = c(0, 5, 10)) + one_component_all +
               Sex + First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
               hyperchol + hypertension_1 + dm + Max_IMT_complete, data = data_CVD_nonlinear)
print(anova(rg1, rg2), digits = 4)
#full cohort
rg1 <- with(incident_CVD_5_long,
            coxph(Surv(stroke_time, stat == 1) ~
                    I(METh_vigorous/20) +
                    Sex + First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
                    hyperchol + hypertension_1 + dm + Max_IMT_complete))
summary(pool(rg1))
full_cohort <- format_pooled_results(summary(pool(rg1))[1,])
full_cohort
#no_one_component
rg1 <- with(incident_CVD_5_long_no_one_component,
            coxph(Surv(stroke_time, stat == 1) ~
                    I(METh_vigorous/20) +
                    Sex + First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
                    hyperchol + hypertension_1 + dm + Max_IMT_complete))
summary(pool(rg1))
no_one_component <- format_pooled_results(summary(pool(rg1))[1,])
no_one_component
#one_component
rg1 <- with(incident_CVD_5_long_one_component,
            coxph(Surv(stroke_time, stat == 1) ~
                    I(METh_vigorous/20) +
                    Sex + First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
                    hyperchol + hypertension_1 + dm + Max_IMT_complete))
summary(pool(rg1))
one_component <- format_pooled_results(summary(pool(rg1))[1,])
one_component
#interaction
rg1 <- with(incident_CVD_5_long,
            coxph(Surv(stroke_time, stat == 1) ~
                    I(METh_vigorous/20)*one_component_all +
                    Sex + First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
                    hyperchol + hypertension_1 + dm + Max_IMT_complete))
summary(pool(rg1))
fit_rg0 <- summary(pool(rg1))[14,]
p_interaction <- c(round(fit_rg0$p.value, digits = 4))
p_interaction
vigorous_stroke <- data.frame(full_cohort, no_one_component, one_component, p_interaction)
vigorous_stroke


final_one_component <- rbind(total_stroke,modvig_stroke,
                             moderate_stroke,vigorous_stroke)
final_one_component

####sex difference----
#total
rg0 <- with(incident_CVD_5_long,
            coxph(Surv(stroke_time, stat == 1) ~
                    I(total_MET/20)*Sex +
                    First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
                    hyperchol + hypertension_1 + dm + Max_IMT_complete))
summary(pool(rg0))
#modvig
rg0 <- with(incident_CVD_5_long,
            coxph(Surv(stroke_time, stat == 1) ~
                    I(METh_modvig/20)*Sex +
                    First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
                    hyperchol + hypertension_1 + dm + Max_IMT_complete))
summary(pool(rg0))
#moderate
rg0 <- with(incident_CVD_5_long,
            coxph(Surv(stroke_time, stat == 1) ~
                    I(METh_vigorous/20)*Sex +
                    First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
                    hyperchol + hypertension_1 + dm + Max_IMT_complete))
summary(pool(rg0))
#vigorous
rg0 <- with(incident_CVD_5_long,
            coxph(Surv(stroke_time, stat == 1) ~
                    I(METh_moderate/20)*Sex +
                    First_age + BMI + smokingstatus_2 + rs_cohort + edu_2 + 
                    hyperchol + hypertension_1 + dm + Max_IMT_complete))
summary(pool(rg0))

####table1----
incident_CVD_5_check_1 <- subset(incident_CVD_5_check, .imp == 1)
dput(names(incident_CVD_5_check_1))
incident_CVD_5_check_1$Sex <- factor(incident_CVD_5_check_1$Sex, labels = c("Man", "Woman"))
incident_CVD_5_check_1 <- left_join(incident_CVD_5_check_1, ID_follow, by = "ergoid")
incident_CVD_5_check_1$follow[is.na(incident_CVD_5_check_1$follow)] <- 0
render.median.IQR <- function(x, ...) {
  c('', 
    `Mean (SD)` = sprintf("%s (%s)", round(mean(x), 1), round(sd(x), 1)),
    `Median [IQR]` = sprintf("%s [%s, %s]",  round(median(x), 1), 
                             round(quantile(x, 0.25), 1), round(quantile(x, 0.75), 1)))
}
table_1_1 <- table1(~ Sex + First_age + BMI + smokingstatus + edu_2 + 
                      hyperchol + hypertension_1 + dm
                    | follow,
                    data = incident_CVD_5_check_1,
                    render.continuous=c(.="Mean (SD)"))
table_1_1 <- as.data.frame(table_1_1)
table_1_1
table_1_1 <- table1(~ total_MET + METh_modvig + METh_moderate + METh_vigorous
                    | follow,
                    data = incident_CVD_5_check_1,
                    render= render.median.IQR)
table_1_1 <- as.data.frame(table_1_1)
table_1_1

#cognitive function
MMSE_ej <- read_sav('ergo_jong_MMSE_Luoshiyuan.sav')
dput(names(MMSE_ej))
MMSE_ej <- subset(MMSE_ej, select = c("ergoid", "dmmse"))
MMSE_ej <- rename(MMSE_ej, mmse = dmmse)
incident_CVD_5_check_1_ej <- subset(incident_CVD_5_check_1, rs_cohort == 3)
incident_CVD_5_check_1_ej <- left_join(incident_CVD_5_check_1_ej, MMSE_ej, by = "ergoid")

MMSE_e5 <- read_sav('ergo_5_MMSE_Luoshiyuan.sav')
dput(names(MMSE_e5))
MMSE_e5 <- subset(MMSE_e5, select = c("ergoid", "emmse"))
MMSE_e5 <- rename(MMSE_e5, mmse = emmse)
incident_CVD_5_check_1_e5 <- subset(incident_CVD_5_check_1, rs_cohort != 3)
incident_CVD_5_check_1_e5 <- left_join(incident_CVD_5_check_1_e5, MMSE_e5, by = "ergoid")
incident_CVD_5_check_1_mmse <- rbind(incident_CVD_5_check_1_ej, incident_CVD_5_check_1_e5)
incident_CVD_5_check_1_mmse$mmse[is.na(incident_CVD_5_check_1_mmse$mmse)] <- 0
incident_CVD_5_check_1_mmse <- mutate(incident_CVD_5_check_1_mmse, mmse_dys = case_when(
  mmse < 24 ~ 1,
  TRUE ~ 0
))
incident_CVD_5_check_1_mmse$mmse_dys <- factor(incident_CVD_5_check_1_mmse$mmse_dys)
summary(incident_CVD_5_check_1_mmse)

table_mmse <- table1(~  mmse + mmse_dys
                    | follow,
                    data = incident_CVD_5_check_1_mmse,
                    render.continuous=c(.="Mean (SD)"))
table_mmse <- as.data.frame(table_mmse)
table_mmse

table_mmse <- table1(~  mmse
                     | follow,
                     data = incident_CVD_5_check_1_mmse,
                     render= render.median.IQR)
table_mmse <- as.data.frame(table_mmse)
table_mmse

dput(names(incident_CVD_5_check_1_mmse))
cognitive_func <- subset(incident_CVD_5_check_1_mmse, select = c("ergoid","mmse","mmse_dys"))


ID_total_check <- subset(incident_CVD_5_check_1, select = "ergoid")
ID_total_check$total_check <- 1
incident_CVD_5_check_1_perartery <- data_secondscan_5
incident_CVD_5_check_1_perartery <- subset(incident_CVD_5_check_1_perartery, select = c("Cal_1", "IPH_1", "Lipid_1", "Maximum_1", "ergoid", "side"))
incident_CVD_5_check_1_perartery <- left_join(incident_CVD_5_check_1_perartery, ID_total_check, by = "ergoid")
incident_CVD_5_check_1_perartery <- left_join(incident_CVD_5_check_1_perartery,ID_use,by= c("ergoid", "side"))
incident_CVD_5_check_1_perartery$follow[is.na(incident_CVD_5_check_1_perartery$follow)] <- 0
incident_CVD_5_check_1_perartery <- subset(incident_CVD_5_check_1_perartery,!is.na(total_check)&!is.na(Cal_1))
incident_CVD_5_check_1_perartery$Lipid_1 <- factor(incident_CVD_5_check_1_perartery$Lipid_1)
incident_CVD_5_check_1_perartery$IPH_1 <- factor(incident_CVD_5_check_1_perartery$IPH_1)

table_1_1 <- table1(~ Maximum_1+IPH_1+Lipid_1
                    | follow,
                    data = incident_CVD_5_check_1_perartery)
table_1_1 <- as.data.frame(table_1_1)
table_1_1


####sex difference table1----
incident_CVD_5_check_1 <- subset(incident_CVD_5_check, .imp == 1)
dput(names(incident_CVD_5_check_1))
incident_CVD_5_check_1$Sex <- factor(incident_CVD_5_check_1$Sex, labels = c("Man", "Woman"))
incident_CVD_5_check_1 <- left_join(incident_CVD_5_check_1, ID_follow, by = "ergoid")
incident_CVD_5_check_1$follow[is.na(incident_CVD_5_check_1$follow)] <- 0
render.median.IQR <- function(x, ...) {
  c('', 
    `Mean (SD)` = sprintf("%s (%s)", round(mean(x), 1), round(sd(x), 1)),
    `Median [IQR]` = sprintf("%s [%s, %s]",  round(median(x), 1), 
                             round(quantile(x, 0.25), 1), round(quantile(x, 0.75), 1)))
}
incident_CVD_5_check_1$follow <- factor(incident_CVD_5_check_1$follow)
table_1_1 <- table1(~ follow + First_age + BMI + smokingstatus + edu_2 + 
                      hyperchol + hypertension_1 + dm
                    | Sex,
                    data = incident_CVD_5_check_1,
                    render.continuous=c(.="Mean (SD)"))
table_1_1 <- as.data.frame(table_1_1)
table_1_1
table_1_1 <- table1(~ total_MET + METh_modvig + METh_moderate + METh_vigorous
                    | Sex,
                    data = incident_CVD_5_check_1,
                    render= render.median.IQR)
table_1_1 <- as.data.frame(table_1_1)
table_1_1

#ognitive function
MMSE_ej <- read_sav('ergo_jong_MMSE_Luoshiyuan.sav')
dput(names(MMSE_ej))
MMSE_ej <- subset(MMSE_ej, select = c("ergoid", "dmmse"))
MMSE_ej <- rename(MMSE_ej, mmse = dmmse)
incident_CVD_5_check_1_ej <- subset(incident_CVD_5_check_1, rs_cohort == 3)
incident_CVD_5_check_1_ej <- left_join(incident_CVD_5_check_1_ej, MMSE_ej, by = "ergoid")

MMSE_e5 <- read_sav('ergo_5_MMSE_Luoshiyuan.sav')
dput(names(MMSE_e5))
MMSE_e5 <- subset(MMSE_e5, select = c("ergoid", "emmse"))
MMSE_e5 <- rename(MMSE_e5, mmse = emmse)
incident_CVD_5_check_1_e5 <- subset(incident_CVD_5_check_1, rs_cohort != 3)
incident_CVD_5_check_1_e5 <- left_join(incident_CVD_5_check_1_e5, MMSE_e5, by = "ergoid")
incident_CVD_5_check_1_mmse <- rbind(incident_CVD_5_check_1_ej, incident_CVD_5_check_1_e5)
incident_CVD_5_check_1_mmse$mmse[is.na(incident_CVD_5_check_1_mmse$mmse)] <- 0
incident_CVD_5_check_1_mmse <- mutate(incident_CVD_5_check_1_mmse, mmse_dys = case_when(
  mmse < 24 ~ 1,
  TRUE ~ 0
))
incident_CVD_5_check_1_mmse$mmse_dys <- factor(incident_CVD_5_check_1_mmse$mmse_dys)

table_mmse <- table1(~  mmse + mmse_dys
                     | Sex,
                     data = incident_CVD_5_check_1_mmse,
                     render.continuous=c(.="Mean (SD)"))
table_mmse <- as.data.frame(table_mmse)
table_mmse


ID_total_check <- subset(incident_CVD_5_check_1, select = "ergoid")
ID_total_check$total_check <- 1
incident_CVD_5_check_1_perartery <- data_secondscan_5
incident_CVD_5_check_1_perartery <- subset(incident_CVD_5_check_1_perartery, select = c("Cal_1", "IPH_1", "Lipid_1", "Maximum_1", "ergoid", "side", "Sex"))
incident_CVD_5_check_1_perartery <- left_join(incident_CVD_5_check_1_perartery, ID_total_check, by = "ergoid")
incident_CVD_5_check_1_perartery <- left_join(incident_CVD_5_check_1_perartery,ID_use,by= c("ergoid", "side"))
incident_CVD_5_check_1_perartery$follow[is.na(incident_CVD_5_check_1_perartery$follow)] <- 0
incident_CVD_5_check_1_perartery <- subset(incident_CVD_5_check_1_perartery,!is.na(total_check)&!is.na(Cal_1))
incident_CVD_5_check_1_perartery$Lipid_1 <- factor(incident_CVD_5_check_1_perartery$Lipid_1)
incident_CVD_5_check_1_perartery$IPH_1 <- factor(incident_CVD_5_check_1_perartery$IPH_1)

table_1_1 <- table1(~ Maximum_1+IPH_1+Lipid_1
                    | Sex,
                    data = incident_CVD_5_check_1_perartery)
table_1_1 <- as.data.frame(table_1_1)
table_1_1



