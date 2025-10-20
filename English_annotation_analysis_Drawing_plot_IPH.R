extract_summary_row <- function(tbl) {
  new_row <- paste0(tbl["Sum", ], "(", tbl["IPH", ], ")")
  return(new_row)
}

####IPH incident----
incident_followup_1_long_1 <- incident_followup_1_long[incident_followup_1_long$scan == 1,]
ergoid <- incident_followup_1_long_1$ergoid[incident_followup_1_long_1$IPH == 1]
ergoid <- unique(ergoid)
ergoid_0 <- unique(incident_followup_1_long_1$ergoid)
ergoid <- setdiff(ergoid_0,ergoid)
ergoid <- as.data.frame(ergoid)
ergoid$mark <- 1
IPH_ID <- data.frame(ergoid)
IPH_ID <- IPH_ID %>%
  distinct(ergoid, .keep_all = TRUE)

####inverse probability weighting----
incident_CVD_5_check_1 <- subset(incident_CVD_5_check, .imp == 1)#The response variable is 'follow'; 0 = lost to follow-up, 1 = followed.
dput(names(incident_CVD_5_check_1))
incident_CVD_5_check_1$Sex <- factor(incident_CVD_5_check_1$Sex, labels = c("Man", "Woman"))
IPH_ID_1 <- rename(IPH_ID, follow = mark)
incident_CVD_5_check_1 <- left_join(incident_CVD_5_check_1, IPH_ID_1, by = "ergoid")
incident_CVD_5_check_1$follow[is.na(incident_CVD_5_check_1$follow)] <- 0
incident_CVD_5_check_1 <- subset(incident_CVD_5_check_1, IPH_all == "no_IPH")####At baseline, participants without IPH form the base population for this analysis.
#Age and rs_cohort are too dominant; including them removes any associations.
incident_CVD_5_check_1 <- mutate(incident_CVD_5_check_1, age_group=case_when(
  First_age < median(incident_CVD_5_check_1$First_age) ~ 1,
  TRUE ~ 2))
incident_CVD_5_check_1$age_group <- factor(incident_CVD_5_check_1$age_group)
incident_CVD_5_check_1 <- left_join(incident_CVD_5_check_1, cognitive_func)
propensity_score_model <- glm(follow ~  Sex + age_group + BMI + smokingstatus_2 + edu_2 + hyperchol + hypertension_1 + dm + mmse_dys,
                              family = binomial, data = incident_CVD_5_check_1)

#Get raw weights#Get raw weights
incident_CVD_5_check_1$ps <- predict(propensity_score_model, type = "response")
#Get stabilized weights
p_follow <- mean(incident_CVD_5_check_1$follow)
incident_CVD_5_check_1 <- mutate(incident_CVD_5_check_1, IPW_stabilized = case_when(
  follow == 1 ~ p_follow / ps,
  follow == 0 ~ (1 - p_follow) / (1 - ps)),
  IPW = 1/ps)
IPW_final <- subset(incident_CVD_5_check_1, select = c("ergoid", "IPW_stabilized"))

####Merge IPW---
incident_followup_1_long_2 <- incident_followup_1_long[incident_followup_1_long$scan == 2,]
incident_followup_1_long_2 <- left_join(incident_followup_1_long_2, IPW_final, by = "ergoid")
incident_followup_1_long_2 <- left_join(incident_followup_1_long_2, IPH_ID, by = c("ergoid"))
incident_followup_1_long_3 <- incident_followup_1_long_2[!is.na(incident_followup_1_long_2$mark),]
incident_followup_1_long_IPH_case <- subset(incident_followup_1_long_3, .imp == 1 & scan == 2)
table(incident_followup_1_long_IPH_case$IPH)#147/937
incident_followup_1_long_3 <- as.mids(incident_followup_1_long_3, .imp = ".imp", .id = ".id")

#####Regression results for continuous variables----
extract_exp_ci_pval <- function(model_summary) {
  # Compute 95% CI on the log scale
  model_summary$low <- model_summary$estimate - 1.96 * model_summary$std.error
  model_summary$high <- model_summary$estimate + 1.96 * model_summary$std.error
  # Subset relevant columns
  result <- subset(model_summary, select = c("term", "estimate", "low", "high", "p.value"))
  # Exponentiate and round
  result$estimate <- round(exp(result$estimate), digits = 2)
  result$low <- round(exp(result$low), digits = 2)
  result$high <- round(exp(result$high), digits = 2)
  result$p.value <- round(result$p.value, digits = 3)
  return(result)
}

####Weighted
fit_continuous <- data.frame()
#total
rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + I(total_MET/20)
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
fit_rg1 <- extract_exp_ci_pval(summary(pool(rg0))[14,])
fit_continuous <- rbind(fit_continuous, fit_rg1)
#modvig
rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + I(METh_modvig/20)
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
fit_rg1 <- extract_exp_ci_pval(summary(pool(rg0))[14,])
fit_continuous <- rbind(fit_continuous, fit_rg1)
#moderate
rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + I(METh_moderate/20)
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
fit_rg1 <- extract_exp_ci_pval(summary(pool(rg0))[14,])
fit_continuous <- rbind(fit_continuous, fit_rg1)
#vigorous
rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + I(METh_vigorous/20)
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
fit_rg1 <- extract_exp_ci_pval(summary(pool(rg0))[14,])
fit_continuous_IPH <- rbind(fit_continuous, fit_rg1)
fit_continuous_IPH$extract <- paste0(fit_continuous_IPH$estimate, " (", fit_continuous_IPH$low, ", ", fit_continuous_IPH$high, ")")
fit_continuous_IPH

####sex difference----
#total
rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex*I(total_MET/20) + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
#modvig
rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex*I(METh_modvig/20) + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
#moderate
rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex*I(METh_moderate/20) + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
#vigorous
rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex*I(METh_vigorous/20) + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))


#####Nonlinear----
data_nonlinear <- subset(incident_followup_1_long_2, .imp == 1)
rg1 <- geeglm(IPH
              ~ Sex + First_age + followuptime + BMI + smokingstatus 
              + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
              + total_MET
              , family = binomial, id = ergoid, corstr = "independence", data = data_nonlinear,  weights = IPW_stabilized)
rg2 <- geeglm(IPH
              ~ Sex + First_age + followuptime + BMI + smokingstatus 
              + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
              + ns(total_MET, 3)
              , family = binomial, id = ergoid, corstr = "independence", data = data_nonlinear,  weights = IPW_stabilized)
print(anova(rg1, rg2, test = "Wald"), digits = 4)
QIC(rg1, rg2)


rg1 <- geeglm(IPH
              ~ Sex + First_age + followuptime + BMI + smokingstatus
              + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
              + METh_modvig
              , family = binomial, id = ergoid, corstr = "independence", data = data_nonlinear,  weights = IPW_stabilized)
rg2 <- geeglm(IPH
              ~ Sex + First_age + followuptime + BMI + smokingstatus
              + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
              + ns(METh_modvig, 3)
              , family = binomial, id = ergoid, corstr = "independence", data = data_nonlinear,  weights = IPW_stabilized)
print(anova(rg1, rg2, test = "Wald"), digits = 4)
QIC(rg1, rg2)


rg1 <- geeglm(IPH
              ~ Sex + First_age + followuptime + BMI + smokingstatus
              + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
              + METh_moderate
              , family = binomial, id = ergoid, corstr = "independence", data = data_nonlinear,  weights = IPW_stabilized)
rg2 <- geeglm(IPH
              ~ Sex + First_age + followuptime + BMI + smokingstatus
              + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
              + ns(METh_moderate, 3)
              , family = binomial, id = ergoid, corstr = "independence", data = data_nonlinear,  weights = IPW_stabilized)
print(anova(rg1, rg2, test = "Wald"), digits = 4)
QIC(rg1, rg2)


rg1 <- geeglm(IPH
              ~ Sex + First_age + followuptime + BMI + smokingstatus
              + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
              + METh_vigorous
              , family = binomial, id = ergoid, corstr = "independence", data = data_nonlinear,  weights = IPW_stabilized)
rg2 <- geeglm(IPH
              ~ Sex + First_age + followuptime + BMI + smokingstatus 
              + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
              + ns(METh_vigorous, knots = c(1, 5, 10))
              , family = binomial, id = ergoid, corstr = "independence", data = data_nonlinear,  weights = IPW_stabilized)
print(anova(rg1, rg2, test = "Wald"), digits = 4)
QIC(rg1, rg2)



####Check population proportions----
incident_followup_1_long_4 <- complete(incident_followup_1_long_3, action = 'long')
incident_followup_1_long_4 <- incident_followup_1_long_4[incident_followup_1_long_4$.imp == 1,]
incident_followup_1_long_4$IPH <- factor(incident_followup_1_long_4$IPH, labels = c("no_IPH", "IPH"))
incident_followup_1_long_4$Sex <- factor(incident_followup_1_long_4$Sex, labels = c("Man", "Woman"))

tb1 <- addmargins(table(incident_followup_1_long_4$IPH, incident_followup_1_long_4$total_guideline), 1)
tb2 <- addmargins(table(incident_followup_1_long_4$IPH, incident_followup_1_long_4$total_MET_3), 1)
tb3 <- addmargins(table(incident_followup_1_long_4$IPH, incident_followup_1_long_4$total_MET_4), 1)
tb4 <- addmargins(table(incident_followup_1_long_4$IPH, incident_followup_1_long_4$total_MET_5), 1)
tblist <- list(tb1, tb2, tb3, tb4)
summary_rows_IPH1 <- lapply(tblist, extract_summary_row)
summary_rows_IPH1
tb1 <- addmargins(table(incident_followup_1_long_4$IPH, incident_followup_1_long_4$METh_modvig_guideline), 1)
tb2 <- addmargins(table(incident_followup_1_long_4$IPH, incident_followup_1_long_4$METh_modvig_3), 1)
tb3 <- addmargins(table(incident_followup_1_long_4$IPH, incident_followup_1_long_4$METh_modvig_4), 1)
tb4 <- addmargins(table(incident_followup_1_long_4$IPH, incident_followup_1_long_4$METh_modvig_5), 1)
tblist <- list(tb1, tb2, tb3, tb4)
summary_rows_IPH2 <- lapply(tblist, extract_summary_row)
summary_rows_IPH2
tb1 <- addmargins(table(incident_followup_1_long_4$IPH, incident_followup_1_long_4$METh_moderate_3), 1)
tb2 <- addmargins(table(incident_followup_1_long_4$IPH, incident_followup_1_long_4$METh_moderate_4), 1)
tb3 <- addmargins(table(incident_followup_1_long_4$IPH, incident_followup_1_long_4$METh_moderate_5), 1)
tblist <- list(tb1, tb2, tb3)
summary_rows_IPH3 <- lapply(tblist, extract_summary_row)
summary_rows_IPH3
tb1 <- addmargins(table(incident_followup_1_long_4$IPH, incident_followup_1_long_4$Vigorous_3), 1)
tb2 <- addmargins(table(incident_followup_1_long_4$IPH, incident_followup_1_long_4$Vigorous_4), 1)
tb3 <- addmargins(table(incident_followup_1_long_4$IPH, incident_followup_1_long_4$Vigorous_5), 1)
tblist <- list(tb1, tb2, tb3)
summary_rows_IPH4 <- lapply(tblist, extract_summary_row)
summary_rows_IPH4


######Coefficients for quantile comparisons----
#function of extract coefficients
format_pooled_results_component <- function(pooled_model) {
  fit_summary <- pooled_model
  # Calculate 95% CI on the log scale
  fit_summary$low <- fit_summary$estimate - 1.96 * fit_summary$std.error
  fit_summary$high <- fit_summary$estimate + 1.96 * fit_summary$std.error
  # Subset and exponentiate
  fit_subset <- fit_summary[, c("term", "estimate", "low", "high", "p.value")]
  fit_subset$estimate <- round(exp(fit_subset$estimate), 2)
  fit_subset$low <- round(exp(fit_subset$low), 2)
  fit_subset$high <- round(exp(fit_subset$high), 2)
    # Round p-values to 3 decimals (or scientific notation for small values)
  fit_subset$p.value <- ifelse(fit_subset$p.value < 0.001,
                               "<0.001",
                               sprintf("%.3f", fit_subset$p.value))
    # Construct the formatted string
  formatted_result <- paste0(
    fit_subset$estimate, " (", fit_subset$low, ", ", fit_subset$high, "); P = ", fit_subset$p.value
  )
  # Add reference row at the beginning
  formatted_result <- c("Reference", formatted_result)
  return(formatted_result)
}

####total PA----
summary_rows_IPH1_final <- summary_rows_IPH1
#guideline
rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + total_guideline
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
summary_rows_IPH1_final <- append(summary_rows_IPH1_final, list(format_pooled_results_component(summary(pool(rg0))[14:15,])), after = 1)
summary_rows_IPH1_final
#Tertile
rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + total_MET_3
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
summary_rows_IPH1_final <- append(summary_rows_IPH1_final, list(format_pooled_results_component(summary(pool(rg0))[14:15,])), after = 3)
summary_rows_IPH1_final

#Quartile
rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + total_MET_4
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
summary_rows_IPH1_final <- append(summary_rows_IPH1_final, list(format_pooled_results_component(summary(pool(rg0))[14:16,])), after = 5)
summary_rows_IPH1_final

#Quintile
rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + total_MET_5
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
summary_rows_IPH1_final <- append(summary_rows_IPH1_final, list(format_pooled_results_component(summary(pool(rg0))[14:17,])), after = 7)
summary_rows_IPH1_final

#algin output table
final_list <- summary_rows_IPH1_final
total_IPH_output <- lapply(seq(1, length(final_list), by = 2), function(i){
  t(data.frame(Numeric = final_list[[i]], Characteristic = final_list[[i + 1]]))
})
total_IPH <- capture.output(print(total_IPH_output, quote = TRUE, row.names = FALSE))
total_IPH


####modvig PA----
summary_rows_IPH2_final <- summary_rows_IPH2
#guideline
rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + METh_modvig_guideline
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
summary_rows_IPH2_final <- append(summary_rows_IPH2_final, list(format_pooled_results_component(summary(pool(rg0))[14:15,])), after = 1)
summary_rows_IPH2_final
#Tertile
rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + METh_modvig_3
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
summary_rows_IPH2_final <- append(summary_rows_IPH2_final, list(format_pooled_results_component(summary(pool(rg0))[14:15,])), after = 3)
summary_rows_IPH2_final

#Quartile
rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + METh_modvig_4
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
summary_rows_IPH2_final <- append(summary_rows_IPH2_final, list(format_pooled_results_component(summary(pool(rg0))[14:16,])), after = 5)
summary_rows_IPH2_final

#Quintile
rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + METh_modvig_5
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
summary_rows_IPH2_final <- append(summary_rows_IPH2_final, list(format_pooled_results_component(summary(pool(rg0))[14:17,])), after = 7)
summary_rows_IPH2_final

#algin output table
final_list <- summary_rows_IPH2_final
modvig_IPH_output <- lapply(seq(1, length(final_list), by = 2), function(i){
  t(data.frame(Numeric = final_list[[i]], Characteristic = final_list[[i + 1]]))
})
modvig_IPH <- capture.output(print(modvig_IPH_output, quote = TRUE, row.names = FALSE))
modvig_IPH


####moderate----
summary_rows_IPH3_final <- summary_rows_IPH3
#Tertile
rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + METh_moderate_3
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
summary_rows_IPH3_final <- append(summary_rows_IPH3_final, list(format_pooled_results_component(summary(pool(rg0))[14:15,])), after = 1)
summary_rows_IPH3_final


#Quartile
rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + METh_moderate_4
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
summary_rows_IPH3_final <- append(summary_rows_IPH3_final, list(format_pooled_results_component(summary(pool(rg0))[14:16,])), after = 3)
summary_rows_IPH3_final

#Quintile
rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + METh_moderate_5
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
summary_rows_IPH3_final <- append(summary_rows_IPH3_final, list(format_pooled_results_component(summary(pool(rg0))[14:17,])), after = 5)
summary_rows_IPH3_final

#algin output table
final_list <- summary_rows_IPH3_final
moderate_IPH_output <- lapply(seq(1, length(final_list), by = 2), function(i){
  t(data.frame(Numeric = final_list[[i]], Characteristic = final_list[[i + 1]]))
})
moderate_IPH <- capture.output(print(moderate_IPH_output, quote = TRUE, row.names = FALSE))
moderate_IPH


####vigorous----
summary_rows_IPH4_final <- summary_rows_IPH4
#Tertile
rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + Vigorous_3
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
summary_rows_IPH4_final <- append(summary_rows_IPH4_final, list(format_pooled_results_component(summary(pool(rg0))[14:15,])), after = 1)
summary_rows_IPH4_final


#Quartile
rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + Vigorous_4
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
summary_rows_IPH4_final <- append(summary_rows_IPH4_final, list(format_pooled_results_component(summary(pool(rg0))[14:16,])), after = 3)
summary_rows_IPH4_final

#Quintile
rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + Vigorous_5
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
summary_rows_IPH4_final <- append(summary_rows_IPH4_final, list(format_pooled_results_component(summary(pool(rg0))[14:17,])), after = 5)
summary_rows_IPH4_final

#algin output table
final_list <- summary_rows_IPH4_final
vigorous_IPH_output <- lapply(seq(1, length(final_list), by = 2), function(i){
  t(data.frame(Numeric = final_list[[i]], Characteristic = final_list[[i + 1]]))
})
vigorous_IPH <- capture.output(print(vigorous_IPH_output, quote = TRUE, row.names = FALSE))
vigorous_IPH

####Final table output----
fit_continuous_IPH

total_IPH
modvig_IPH
moderate_IPH
vigorous_IPH
IPH_output <- list(total_IPH, modvig_IPH, moderate_IPH, vigorous_IPH)
IPH_output


#####################################################Plotting: Total PA------------------------------
######################guideline cutoff-------------------
#total_guideline tertile
dput(names(incident_followup))

cutvalue <- incident_followup %>%
  group_by(total_guideline) %>%
  summarize(quant0 = quantile(total_MET, probs = 0),
            quant25 = quantile(total_MET, probs = 0.25),
            quant50 = quantile(total_MET, probs = 0.50),
            quant75 = quantile(total_MET, probs = 0.75),
            quant100 = quantile(total_MET, probs = 1.00))
cutvalue
cutvalue_low1 <- cutvalue$quant0[1]
cutvalue_low2 <- cutvalue$quant0[2]
cutvalue_low3 <- cutvalue$quant0[3]
cutvalue_high1 <- cutvalue$quant100[1]
cutvalue_high2 <- cutvalue$quant100[2]
cutvalue_high3 <- 135
cutvalue_1 <- (cutvalue_low2 + cutvalue_high1)/2
cutvalue_2 <- (cutvalue_low3 + cutvalue_high2)/2
cutvalue_3 <- cutvalue_high3
cutvalue_1
cutvalue_2
cutvalue_3
X_position_1 <- (0 + cutvalue_1) / 2
X_position_2 <- (cutvalue_1 + cutvalue_2) / 2
X_position_3 <- (cutvalue_2 + cutvalue_3) / 2
X_position_1
X_position_2
X_position_3
X_position <- c(X_position_1, X_position_2, X_position_3)

#guideline
rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + total_guideline
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
fit_rg0 <- summary(pool(rg0))
fit_rg0$low <- fit_rg0$estimate - 1.96 * fit_rg0$std.error
fit_rg0$high <- fit_rg0$estimate + 1.96 * fit_rg0$std.error
fit_rg1 <- subset(fit_rg0, select = c("term", "estimate", "low", "high"))[14:15,]
fit_rg1$estimate <- round(exp(fit_rg1$estimate), digits = 2)
fit_rg1$low <- round(exp(fit_rg1$low), digits = 2)
fit_rg1$high <- round(exp(fit_rg1$high), digits = 2)
contrast_1 <- data.frame(term = "contrast", estimate = 1, low = 1, high = 1)
fit_rg1 <- rbind(contrast_1, fit_rg1)
fit_rg1
fit_rg2 <- cbind(fit_rg1, X_position)
dput(names(fit_rg2))

#Extract OR, CI, and p-value
fit_final <- subset(fit_rg0, select = c("term", "estimate", "low", "high", "p.value"))[14:15,]
fit_final$estimate <- round(exp(fit_final$estimate), digits = 2)
fit_final$low <- round(exp(fit_final$low), digits = 2)
fit_final$high <- round(exp(fit_final$high), digits = 2)
fit_final$p.value <- round(fit_final$p.value, digits = 3)
fit_final_IPH <- fit_final

####Define the X-axis position
Y_location <- 1

####Define colors
color_point <- "rgb(80, 29, 138)"
color <- "rgba(80, 29, 138, 1)"
color_trans <- "rgba(80, 29, 138, 0.8)"
color_vertex <- c(80, 29, 138)

####Plot points
fit_rg2
fit_rg2$neg_error <- fit_rg2$estimate - fit_rg2$low
fit_rg2$pos_error <- fit_rg2$high - fit_rg2$estimate
fig1_0 <- plot_ly() %>%
  add_markers(data = fit_rg2, 
              x = fit_rg2$X_position, 
              y = Y_location, 
              z = fit_rg2$estimate,
              type = "scatter3d", 
              mode = "markers",
              marker = list(size = 1.5, color = color_point, opacity = 1),
              showlegend = FALSE,
              error_z = list(type = "data", 
                             array = ~pos_error, 
                             arrayminus = ~neg_error, 
                             color = color_point,
                             thickness = 5,#Since this is a 3D plot, thickness has no effect here; consider adding lines manually for better control.
                             width = 5)) %>% 
  layout(
    scene = list(
      xaxis = list(title = "Total physical activity<br>(unit: METh-week)", 
                   range = c(0, cutvalue_high3), 
                   fixedrange = TRUE, 
                   showticklabels = FALSE, 
                   showgrid = FALSE, 
                   zeroline = FALSE, 
                   showline = FALSE),
      yaxis = list(title = "", 
                   range = c(-0.5, 16.5), 
                   fixedrange = TRUE, 
                   showticklabels = FALSE, 
                   showgrid = FALSE, 
                   zeroline = FALSE, 
                   showline = FALSE), 
      zaxis = list(title = "", 
                   range = c(0, 6), 
                   fixedrange = TRUE, 
                   showticklabels = FALSE, 
                   showgrid = FALSE, 
                   zeroline = FALSE, 
                   showline = FALSE))
  ) %>%
  add_trace(
    x = c(0), y = c(0), z = c(2.2),  # Adjust Z-axis title position
    type = "scatter3d",
    mode = "text",
    text = "Odds ratio of <br> incident IPH",
    textfont = list(size = 12),
    showlegend = FALSE
  )
fig1_0

####Draw the corresponding X-axis
X_axis <- data.frame(x = c(0, cutvalue_high3),
                     y = c(Y_location, Y_location),
                     z = c(0, 0))
#Add an arrow at the end of the X-axis
arrow_x <- cutvalue_high3
arrow_y <- Y_location
arrow_z <- 0
fig1_1 <- fig1_0 %>% 
  add_trace(
    x = X_axis$x,
    y = X_axis$y,
    z = X_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3, dash = "dash"),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(1),  #Arrow points along the X-axis
    v = list(0),
    w = list(0),
    sizemode = "absolute",
    sizeref = 5,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
fig1_1
# Add quantile labels on the X-axis
fig1_1 <- fig1_1 %>%
  add_trace(
    x = cutvalue_1,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_1, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = cutvalue_2,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_2, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig1_1

####Draw the corresponding Z-axis
Z_axis <- data.frame(x = c(0, 0),
                     y = c(Y_location, Y_location),
                     z = c(0, 2))
# Add an arrow at the end of the Z-axis
arrow_x <- 0
arrow_y <- Y_location
arrow_z <- 2
fig1_2 <- fig1_1 %>% 
  add_trace(
    x = Z_axis$x,
    y = Z_axis$y,
    z = Z_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(0),  
    v = list(0),
    w = list(1),#Arrow points along the Z-axis
    sizemode = "absolute",
    sizeref = 0.2,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add labels on the Z-axis
fig1_3 <- fig1_2 %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = 0,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 1,
    type = "scatter3d",
    mode = "text",
    text = 1,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 2,
    type = "scatter3d",
    mode = "text",
    text = 2,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig1_3


####Add a visual aid plane for the reference line
square_x <- c(0, 0, cutvalue_high3, cutvalue_high3)
square_z <- c(0, 1, 1, 0)
square_y1 <- c(Y_location, Y_location, Y_location, Y_location)
i_set <- c(0,2)
j_set <- c(1,3)
k_set <- c(2,0)
vertex_colors <- matrix(rep(color_vertex, each = 4), nrow = 4)/255
fig1_4 <- fig1_3  %>%
  add_trace(
    x = square_x, 
    z = square_z, 
    y = square_y1, 
    i = i_set, 
    j = j_set, 
    k = k_set,
    type = "mesh3d",
    vertexcolor = vertex_colors,
    opacity = 0.2,
    showscale = FALSE)
fig1_4

####Add auxiliary lines at the quantiles
line_x1 <- c(cutvalue_1, cutvalue_1)
line_z1 <- c(0, 1)
line_y1 <- c(Y_location, Y_location)
line_x2 <- c(cutvalue_2, cutvalue_2)
line_z2 <- c(0, 1)
line_y2 <- c(Y_location, Y_location)
fig1_5 <- fig1_4  %>% 
  add_trace(
    x = line_x1,
    y = line_y1,
    z = line_z1,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE) %>% 
  add_trace(
    x = line_x2,
    y = line_y2,
    z = line_z2,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE)
fig1_5

####Label the grouping method at the end of the X-axis
fig1_6 <- fig1_5 %>% 
  add_trace(
    x = cutvalue_high3,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = "Literature-based<br>cutoff",
    textposition = "right",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig1_6



######################Tertiles -------------------
#total_tertiles
cutvalue <- incident_followup %>%
  group_by(total_MET_3) %>%
  summarize(quant0 = quantile(total_MET, probs = 0),
            quant25 = quantile(total_MET, probs = 0.25),
            quant50 = quantile(total_MET, probs = 0.50),
            quant75 = quantile(total_MET, probs = 0.75),
            quant100 = quantile(total_MET, probs = 1.00))
cutvalue
cutvalue_low1 <- cutvalue$quant0[1]
cutvalue_low2 <- cutvalue$quant0[2]
cutvalue_low3 <- cutvalue$quant0[3]
cutvalue_high1 <- cutvalue$quant100[1]
cutvalue_high2 <- cutvalue$quant100[2]
cutvalue_high3 <- 135
cutvalue_1 <- (cutvalue_low2 + cutvalue_high1)/2
cutvalue_2 <- (cutvalue_low3 + cutvalue_high2)/2
cutvalue_3 <- cutvalue_high3
cutvalue_1
cutvalue_2
cutvalue_3
X_position_1 <- (0 + cutvalue_1) / 2
X_position_2 <- (cutvalue_1 + cutvalue_2) / 2
X_position_3 <- (cutvalue_2 + cutvalue_3) / 2
X_position_1
X_position_2
X_position_3
X_position <- c(X_position_1, X_position_2, X_position_3)

rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + total_MET_3
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
fit_rg0 <- summary(pool(rg0))
fit_rg0$low <- fit_rg0$estimate - 1.96 * fit_rg0$std.error
fit_rg0$high <- fit_rg0$estimate + 1.96 * fit_rg0$std.error
fit_rg1 <- subset(fit_rg0, select = c("term", "estimate", "low", "high"))[14:15,]
fit_rg1$estimate <- round(exp(fit_rg1$estimate), digits = 2)
fit_rg1$low <- round(exp(fit_rg1$low), digits = 2)
fit_rg1$high <- round(exp(fit_rg1$high), digits = 2)

contrast_1 <- data.frame(term = "contrast", estimate = 1, low = 1, high = 1)
fit_rg1 <- rbind(contrast_1, fit_rg1)
fit_rg1
fit_rg2 <- cbind(fit_rg1, X_position)
dput(names(fit_rg2))

#Extract OR, CI, and p-value
fit_final <- subset(fit_rg0, select = c("term", "estimate", "low", "high", "p.value"))[14:15,]
fit_final$estimate <- round(exp(fit_final$estimate), digits = 2)
fit_final$low <- round(exp(fit_final$low), digits = 2)
fit_final$high <- round(exp(fit_final$high), digits = 2)
fit_final$p.value <- round(fit_final$p.value, digits = 3)
fit_final_IPH <- rbind(fit_final_IPH, fit_final)


####Define the X-axis position
Y_location <- 5

####Define colors
color_point <- "rgb(28, 128, 65)"
color <- "rgba(28, 128, 65, 1)"
color_trans <- "rgba(28, 128, 65, 0.8)"
color_vertex <- c(28, 128, 65)

####Plot points
fit_rg2
fit_rg2$neg_error <- fit_rg2$estimate - fit_rg2$low
fit_rg2$pos_error <- fit_rg2$high - fit_rg2$estimate
fig2_0 <- fig1_6 %>%
  add_markers(data = fit_rg2, 
              x = fit_rg2$X_position, 
              y = Y_location, 
              z = fit_rg2$estimate,
              type = "scatter3d", 
              mode = "markers",
              marker = list(size = 1.5, color = color_point, opacity = 1),
              showlegend = FALSE,
              error_z = list(type = "data", 
                             array = ~pos_error, 
                             arrayminus = ~neg_error, 
                             color = color_point,
                             thickness = 5,#Since this is a 3D plot, thickness has no effect here; consider adding lines manually for better control.
                             width = 5))
fig2_0

####Draw the corresponding X-axis
X_axis <- data.frame(x = c(0, cutvalue_high3),
                     y = c(Y_location, Y_location),
                     z = c(0, 0))
#Add an arrow at the end of the X-axis
arrow_x <- cutvalue_high3
arrow_y <- Y_location
arrow_z <- 0
fig2_1 <- fig2_0 %>% 
  add_trace(
    x = X_axis$x,
    y = X_axis$y,
    z = X_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3, dash = "dash"),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(1),  #Arrow points along the X-axis
    v = list(0),
    w = list(0),
    sizemode = "absolute",
    sizeref = 5,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add quantile labels on the X-axis
fig2_1 <- fig2_1 %>%
  add_trace(
    x = cutvalue_1,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_1, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = cutvalue_2,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_2, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig2_1

####Draw the corresponding Z-axis
Z_axis <- data.frame(x = c(0, 0),
                     y = c(Y_location, Y_location),
                     z = c(0, 2))
# Add an arrow at the end of the Z-axis
arrow_x <- 0
arrow_y <- Y_location
arrow_z <- 2
fig2_2 <- fig2_1 %>% 
  add_trace(
    x = Z_axis$x,
    y = Z_axis$y,
    z = Z_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(0),  
    v = list(0),
    w = list(1),#Arrow points along the Z-axis
    sizemode = "absolute",
    sizeref = 0.2,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add labels on the Z-axis
fig2_3 <- fig2_2 %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = 0,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 1,
    type = "scatter3d",
    mode = "text",
    text = 1,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 2,
    type = "scatter3d",
    mode = "text",
    text = 2,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig2_3


####Add a visual aid plane for the reference line
square_x <- c(0, 0, cutvalue_high3, cutvalue_high3)
square_z <- c(0, 1, 1, 0)
square_y1 <- c(Y_location, Y_location, Y_location, Y_location)
i_set <- c(0,2)
j_set <- c(1,3)
k_set <- c(2,0)
vertex_colors <- matrix(rep(color_vertex, each = 4), nrow = 4)/255
fig2_4 <- fig2_3  %>%
  add_trace(
    x = square_x, 
    z = square_z, 
    y = square_y1, 
    i = i_set, 
    j = j_set, 
    k = k_set,
    type = "mesh3d",
    vertexcolor = vertex_colors,
    opacity = 0.2,
    showscale = FALSE)
fig2_4

####Add auxiliary lines at the quantiles
line_x1 <- c(cutvalue_1, cutvalue_1)
line_z1 <- c(0, 1)
line_y1 <- c(Y_location, Y_location)
line_x2 <- c(cutvalue_2, cutvalue_2)
line_z2 <- c(0, 1)
line_y2 <- c(Y_location, Y_location)
fig2_5 <- fig2_4  %>% 
  add_trace(
    x = line_x1,
    y = line_y1,
    z = line_z1,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE) %>% 
  add_trace(
    x = line_x2,
    y = line_y2,
    z = line_z2,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE)
fig2_5

####Label the grouping method at the end of the X-axis
fig2_6 <- fig2_5 %>% 
  add_trace(
    x = cutvalue_high3,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = "Tertile",
    textposition = "right",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig2_6




######################Quartiles -------------------
#total_quartiles
cutvalue <- incident_followup %>%
  group_by(total_MET_4) %>%
  summarize(quant0 = quantile(total_MET, probs = 0),
            quant25 = quantile(total_MET, probs = 0.25),
            quant50 = quantile(total_MET, probs = 0.50),
            quant75 = quantile(total_MET, probs = 0.75),
            quant100 = quantile(total_MET, probs = 1.00))
cutvalue
cutvalue_low1 <- cutvalue$quant0[1]
cutvalue_low2 <- cutvalue$quant0[2]
cutvalue_low3 <- cutvalue$quant0[3]
cutvalue_low4 <- cutvalue$quant0[4]
cutvalue_high1 <- cutvalue$quant100[1]
cutvalue_high2 <- cutvalue$quant100[2]
cutvalue_high3 <- cutvalue$quant100[3]
cutvalue_high4 <- 135
cutvalue_1 <- (cutvalue_low2 + cutvalue_high1)/2
cutvalue_2 <- (cutvalue_low3 + cutvalue_high2)/2
cutvalue_3 <- (cutvalue_low4 + cutvalue_high3)/2
cutvalue_4 <- cutvalue_high4
cutvalue_1
cutvalue_2
cutvalue_3
cutvalue_4
X_position_1 <- (0 + cutvalue_1) / 2
X_position_2 <- (cutvalue_1 + cutvalue_2) / 2
X_position_3 <- (cutvalue_2 + cutvalue_3) / 2
X_position_4 <- (cutvalue_3 + cutvalue_4) / 2
X_position_1
X_position_2
X_position_3
X_position_4
X_position <- c(X_position_1, X_position_2, X_position_3, X_position_4)

rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + total_MET_4
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
fit_rg0 <- summary(pool(rg0))
fit_rg0$low <- fit_rg0$estimate - 1.96 * fit_rg0$std.error
fit_rg0$high <- fit_rg0$estimate + 1.96 * fit_rg0$std.error
fit_rg1 <- subset(fit_rg0, select = c("term", "estimate", "low", "high"))[14:16,]
fit_rg1$estimate <- round(exp(fit_rg1$estimate), digits = 2)
fit_rg1$low <- round(exp(fit_rg1$low), digits = 2)
fit_rg1$high <- round(exp(fit_rg1$high), digits = 2)

contrast_1 <- data.frame(term = "contrast", estimate = 1, low = 1, high = 1)
fit_rg1 <- rbind(contrast_1, fit_rg1)
fit_rg1
fit_rg2 <- cbind(fit_rg1, X_position)
dput(names(fit_rg2))
#Extract OR, CI, and p-value
fit_final <- subset(fit_rg0, select = c("term", "estimate", "low", "high", "p.value"))[14:16,]
fit_final$estimate <- round(exp(fit_final$estimate), digits = 2)
fit_final$low <- round(exp(fit_final$low), digits = 2)
fit_final$high <- round(exp(fit_final$high), digits = 2)
fit_final$p.value <- round(fit_final$p.value, digits = 3)
fit_final_IPH <- rbind(fit_final_IPH, fit_final)

####Define the X-axis position
Y_location <- 10

####Define colors
color_point <- "rgb(229, 87, 9)"
color <- "rgba(229, 87, 9, 1)"
color_trans <- "rgba(229, 87, 9, 0.8)"
color_vertex <- c(229, 87, 9)

####Plot points
fit_rg2
fit_rg2$neg_error <- fit_rg2$estimate - fit_rg2$low
fit_rg2$pos_error <- fit_rg2$high - fit_rg2$estimate
fig3_0 <- fig2_6 %>%
  add_markers(data = fit_rg2, 
              x = fit_rg2$X_position, 
              y = Y_location, 
              z = fit_rg2$estimate,
              type = "scatter3d", 
              mode = "markers",
              marker = list(size = 1.5, color = color_point, opacity = 1),
              showlegend = FALSE,
              error_z = list(type = "data", 
                             array = ~pos_error, 
                             arrayminus = ~neg_error, 
                             color = color_point,
                             thickness = 5,#Since this is a 3D plot, thickness has no effect here; consider adding lines manually for better control.
                             width = 5))
fig3_0

####Draw the corresponding X-axis
X_axis <- data.frame(x = c(0, cutvalue_high4),
                     y = c(Y_location, Y_location),
                     z = c(0, 0))
#Add an arrow at the end of the X-axis
arrow_x <- cutvalue_high4
arrow_y <- Y_location
arrow_z <- 0
fig3_1 <- fig3_0 %>% 
  add_trace(
    x = X_axis$x,
    y = X_axis$y,
    z = X_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3, dash = "dash"),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(1),  #Arrow points along the X-axis
    v = list(0),
    w = list(0),
    sizemode = "absolute",
    sizeref = 5,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add quantile labels on the X-axis
fig3_1 <- fig3_1 %>%
  add_trace(
    x = cutvalue_1,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_1, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = cutvalue_2,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_2, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = cutvalue_3,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_3, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig3_1

####Draw the corresponding Z-axis
Z_axis <- data.frame(x = c(0, 0),
                     y = c(Y_location, Y_location),
                     z = c(0, 2))
# Add an arrow at the end of the Z-axis
arrow_x <- 0
arrow_y <- Y_location
arrow_z <- 2
fig3_2 <- fig3_1 %>% 
  add_trace(
    x = Z_axis$x,
    y = Z_axis$y,
    z = Z_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(0),  
    v = list(0),
    w = list(1),#Arrow points along the Z-axis
    sizemode = "absolute",
    sizeref = 0.2,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add labels on the Z-axis
fig3_3 <- fig3_2 %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = 0,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 1,
    type = "scatter3d",
    mode = "text",
    text = 1,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 2,
    type = "scatter3d",
    mode = "text",
    text = 2,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig3_3


####Add a visual aid plane for the reference line
square_x <- c(0, 0, cutvalue_high4, cutvalue_high4)
square_z <- c(0, 1, 1, 0)
square_y1 <- c(Y_location, Y_location, Y_location, Y_location)
i_set <- c(0,2)
j_set <- c(1,3)
k_set <- c(2,0)
vertex_colors <- matrix(rep(color_vertex, each = 4), nrow = 4)/255
fig3_4 <- fig3_3  %>%
  add_trace(
    x = square_x, 
    z = square_z, 
    y = square_y1, 
    i = i_set, 
    j = j_set, 
    k = k_set,
    type = "mesh3d",
    vertexcolor = vertex_colors,
    opacity = 0.2,
    showscale = FALSE)
fig3_4

####Add auxiliary lines at the quantiles
line_x1 <- c(cutvalue_1, cutvalue_1)
line_z1 <- c(0, 1)
line_y1 <- c(Y_location, Y_location)
line_x2 <- c(cutvalue_2, cutvalue_2)
line_z2 <- c(0, 1)
line_y2 <- c(Y_location, Y_location)
line_x3 <- c(cutvalue_3, cutvalue_3)
line_z3 <- c(0, 1)
line_y3 <- c(Y_location, Y_location)
fig3_5 <- fig3_4  %>% 
  add_trace(
    x = line_x1,
    y = line_y1,
    z = line_z1,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE) %>% 
  add_trace(
    x = line_x2,
    y = line_y2,
    z = line_z2,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE) %>% 
  add_trace(
    x = line_x3,
    y = line_y3,
    z = line_z3,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE)
fig3_5

####Label the grouping method at the end of the X-axis
fig3_6 <- fig3_5 %>% 
  add_trace(
    x = cutvalue_high4,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = "Quartile",
    textposition = "right",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig3_6



######################Quintiles -------------------
#total_quintiles
cutvalue <- incident_followup %>%
  group_by(total_MET_5) %>%
  summarize(quant0 = quantile(total_MET, probs = 0),
            quant25 = quantile(total_MET, probs = 0.25),
            quant50 = quantile(total_MET, probs = 0.50),
            quant75 = quantile(total_MET, probs = 0.75),
            quant100 = quantile(total_MET, probs = 1.00))
cutvalue
cutvalue_low1 <- cutvalue$quant0[1]
cutvalue_low2 <- cutvalue$quant0[2]
cutvalue_low3 <- cutvalue$quant0[3]
cutvalue_low4 <- cutvalue$quant0[4]
cutvalue_low5 <- cutvalue$quant0[5]
cutvalue_high1 <- cutvalue$quant100[1]
cutvalue_high2 <- cutvalue$quant100[2]
cutvalue_high3 <- cutvalue$quant100[3]
cutvalue_high4 <- cutvalue$quant100[4]
cutvalue_high5 <- 135
cutvalue_1 <- (cutvalue_low2 + cutvalue_high1)/2
cutvalue_2 <- (cutvalue_low3 + cutvalue_high2)/2
cutvalue_3 <- (cutvalue_low4 + cutvalue_high3)/2
cutvalue_4 <- (cutvalue_low5 + cutvalue_high4)/2
cutvalue_5 <- cutvalue_high5
cutvalue_1
cutvalue_2
cutvalue_3
cutvalue_4
cutvalue_5
X_position_1 <- (0 + cutvalue_1) / 2
X_position_2 <- (cutvalue_1 + cutvalue_2) / 2
X_position_3 <- (cutvalue_2 + cutvalue_3) / 2
X_position_4 <- (cutvalue_3 + cutvalue_4) / 2
X_position_5 <- (cutvalue_5 + cutvalue_4) / 2
X_position_1
X_position_2
X_position_3
X_position_4
X_position_5

X_position <- c(X_position_1, X_position_2, X_position_3, X_position_4, X_position_5)

rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + total_MET_5
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
fit_rg0 <- summary(pool(rg0))
fit_rg0$low <- fit_rg0$estimate - 1.96 * fit_rg0$std.error
fit_rg0$high <- fit_rg0$estimate + 1.96 * fit_rg0$std.error
fit_rg1 <- subset(fit_rg0, select = c("term", "estimate", "low", "high"))[14:17,]
fit_rg1$estimate <- round(exp(fit_rg1$estimate), digits = 2)
fit_rg1$low <- round(exp(fit_rg1$low), digits = 2)
fit_rg1$high <- round(exp(fit_rg1$high), digits = 2)
contrast_1 <- data.frame(term = "contrast", estimate = 1, low = 1, high = 1)
fit_rg1 <- rbind(contrast_1, fit_rg1)
fit_rg1
fit_rg2 <- cbind(fit_rg1, X_position)
dput(names(fit_rg2))
fit_rg2
#Extract OR, CI, and p-value
fit_final <- subset(fit_rg0, select = c("term", "estimate", "low", "high", "p.value"))[14:17,]
fit_final$estimate <- round(exp(fit_final$estimate), digits = 2)
fit_final$low <- round(exp(fit_final$low), digits = 2)
fit_final$high <- round(exp(fit_final$high), digits = 2)
fit_final$p.value <- round(fit_final$p.value, digits = 3)
fit_final_IPH <- rbind(fit_final_IPH, fit_final)
####Define the X-axis position
Y_location <- 16

####Define colors
color_point <- "rgb(43, 106, 153)"
color <- "rgba(43, 106, 153, 1)"
color_trans <- "rgba(43, 106, 153, 0.8)"
color_vertex <- c(43, 106, 153)

####Plot points
fit_rg2
fit_rg2$neg_error <- fit_rg2$estimate - fit_rg2$low
fit_rg2$pos_error <- fit_rg2$high - fit_rg2$estimate
fig4_0 <- fig3_6 %>%
  add_markers(data = fit_rg2, 
              x = fit_rg2$X_position, 
              y = Y_location, 
              z = fit_rg2$estimate,
              type = "scatter3d", 
              mode = "markers",
              marker = list(size = 1.5, color = color_point, opacity = 1),
              showlegend = FALSE,
              error_z = list(type = "data", 
                             array = ~pos_error, 
                             arrayminus = ~neg_error, 
                             color = color_point,
                             thickness = 5,#Since this is a 3D plot, thickness has no effect here; consider adding lines manually for better control.
                             width = 5))
fig4_0

####Draw the corresponding X-axis
X_axis <- data.frame(x = c(0, cutvalue_high5),
                     y = c(Y_location, Y_location),
                     z = c(0, 0))
#Add an arrow at the end of the X-axis
arrow_x <- cutvalue_high5
arrow_y <- Y_location
arrow_z <- 0
fig4_1 <- fig4_0 %>% 
  add_trace(
    x = X_axis$x,
    y = X_axis$y,
    z = X_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3, dash = "dash"),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(1),  #Arrow points along the X-axis
    v = list(0),
    w = list(0),
    sizemode = "absolute",
    sizeref = 5,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add quantile labels on the X-axis
fig4_1 <- fig4_1 %>%
  add_trace(
    x = cutvalue_1,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_1, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = cutvalue_2,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_2, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = cutvalue_3,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_3, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = cutvalue_4,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_4, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig4_1

####Draw the corresponding Z-axis
Z_axis <- data.frame(x = c(0, 0),
                     y = c(Y_location, Y_location),
                     z = c(0, 2))
# Add an arrow at the end of the Z-axis
arrow_x <- 0
arrow_y <- Y_location
arrow_z <- 2
fig4_2 <- fig4_1 %>% 
  add_trace(
    x = Z_axis$x,
    y = Z_axis$y,
    z = Z_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(0),  
    v = list(0),
    w = list(1),#Arrow points along the Z-axis
    sizemode = "absolute",
    sizeref = 0.2,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add labels on the Z-axis
fig4_3 <- fig4_2 %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = 0,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 1,
    type = "scatter3d",
    mode = "text",
    text = 1,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 2,
    type = "scatter3d",
    mode = "text",
    text = 2,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig4_3


####Add a visual aid plane for the reference line
square_x <- c(0, 0, cutvalue_high5, cutvalue_high5)
square_z <- c(0, 1, 1, 0)
square_y1 <- c(Y_location, Y_location, Y_location, Y_location)
i_set <- c(0,2)
j_set <- c(1,3)
k_set <- c(2,0)
vertex_colors <- matrix(rep(color_vertex, each = 4), nrow = 4)/255
fig4_4 <- fig4_3  %>%
  add_trace(
    x = square_x, 
    z = square_z, 
    y = square_y1, 
    i = i_set, 
    j = j_set, 
    k = k_set,
    type = "mesh3d",
    vertexcolor = vertex_colors,
    opacity = 0.2,
    showscale = FALSE)
fig4_4

####Add auxiliary lines at the quantiles
line_x1 <- c(cutvalue_1, cutvalue_1)
line_z1 <- c(0, 1)
line_y1 <- c(Y_location, Y_location)
line_x2 <- c(cutvalue_2, cutvalue_2)
line_z2 <- c(0, 1)
line_y2 <- c(Y_location, Y_location)
line_x3 <- c(cutvalue_3, cutvalue_3)
line_z3 <- c(0, 1)
line_y3 <- c(Y_location, Y_location)
line_x4 <- c(cutvalue_4, cutvalue_4)
line_z4 <- c(0, 1)
line_y4 <- c(Y_location, Y_location)
fig4_5 <- fig4_4  %>% 
  add_trace(
    x = line_x1,
    y = line_y1,
    z = line_z1,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE) %>% 
  add_trace(
    x = line_x2,
    y = line_y2,
    z = line_z2,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE) %>% 
  add_trace(
    x = line_x3,
    y = line_y3,
    z = line_z3,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE) %>% 
  add_trace(
    x = line_x4,
    y = line_y4,
    z = line_z4,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE)
fig4_5

####Label the grouping method at the end of the X-axis
fig4_6 <- fig4_5 %>% 
  add_trace(
    x = cutvalue_high5,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = "Quintile",
    textposition = "right",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig4_6

fig4_7 <- fig4_6 %>% layout(
  scene = list(camera = list(eye = list(x = 0.59097400646812713, y = -1.878239943166165, z = 1.309212875732691),
              aspectmode = "manual",
              aspectratio = list(x = 1, y = 2, z = 1)  # Y-axis is stretched
))
)
###Add the PA guideline-recommended threshold
fig4_7 <- fig4_7 %>% 
  add_text(x = 15, y = 1, z = 0,
           text = '*', textposition = 'middle center',
           textfont = list(size = 20, color = 'red'))%>% 
  add_text(x = 15, y = 5, z = 0,
           text = '*', textposition = 'middle center',
           textfont = list(size = 20, color = 'red'))%>% 
  add_text(x = 15, y = 10, z = 0,
           text = '*', textposition = 'middle center',
           textfont = list(size = 20, color = 'red'))%>% 
  add_text(x = 15, y = 16, z = 0,
           text = '*', textposition = 'middle center',
           textfont = list(size = 20, color = 'red'))%>% 
  add_text(x = 15, y = -0.5, z = 0,
           text = '     15  <br>Guideline<br>Recomm\'s', textposition = 'middle center',
           textfont = list(size = 12, color = 'red')) %>%
  layout(showlegend = FALSE)
  
fig4_7

saveNetwork(fig4_7, "plot.html")
webshot2::webshot("plot.html", 
                  "IPH_total.png",
                  vwidth = 1200, vheight = 800, zoom = 3)


#Crop the PNG so it can be stitched directly later
img_imager <- load.image("IPH_total.png")  # Load image
plot(img_imager)
dim(img_imager)
cropped_img <- imsub(img_imager, x %in% 1050:2700, y %in% 720:1950)
plot(cropped_img)
imager::save.image(cropped_img, "IPH_total.png")




#####################################################Plotting: Moderate-to-vigorous PA------------------------------
######################guideline cutoff-------------------
#METh_modvig_guideline tertile
cutvalue <- incident_followup %>%
  group_by(METh_modvig_guideline) %>%
  summarize(quant0 = quantile(METh_modvig, probs = 0),
            quant25 = quantile(METh_modvig, probs = 0.25),
            quant50 = quantile(METh_modvig, probs = 0.50),
            quant75 = quantile(METh_modvig, probs = 0.75),
            quant100 = quantile(METh_modvig, probs = 1.00))
cutvalue
cutvalue_low1 <- cutvalue$quant0[1]
cutvalue_low2 <- cutvalue$quant0[2]
cutvalue_low3 <- cutvalue$quant0[3]
cutvalue_high1 <- cutvalue$quant100[1]
cutvalue_high2 <- cutvalue$quant100[2]
cutvalue_high3 <- 110
cutvalue_1 <- (cutvalue_low2 + cutvalue_high1)/2
cutvalue_2 <- (cutvalue_low3 + cutvalue_high2)/2
cutvalue_3 <- cutvalue_high3
cutvalue_1
cutvalue_2
cutvalue_3
X_position_1 <- (0 + cutvalue_1) / 2
X_position_2 <- (cutvalue_1 + cutvalue_2) / 2
X_position_3 <- (cutvalue_2 + cutvalue_3) / 2
X_position_1
X_position_2
X_position_3
X_position <- c(X_position_1, X_position_2, X_position_3)

#guideline
rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + METh_modvig_guideline
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
fit_rg0 <- summary(pool(rg0))
fit_rg0$low <- fit_rg0$estimate - 1.96 * fit_rg0$std.error
fit_rg0$high <- fit_rg0$estimate + 1.96 * fit_rg0$std.error
fit_rg1 <- subset(fit_rg0, select = c("term", "estimate", "low", "high"))[14:15,]
fit_rg1$estimate <- round(exp(fit_rg1$estimate), digits = 2)
fit_rg1$low <- round(exp(fit_rg1$low), digits = 2)
fit_rg1$high <- round(exp(fit_rg1$high), digits = 2)

contrast_1 <- data.frame(term = "contrast", estimate = 1, low = 1, high = 1)
fit_rg1 <- rbind(contrast_1, fit_rg1)
fit_rg1
fit_rg2 <- cbind(fit_rg1, X_position)
dput(names(fit_rg2))
#Extract OR, CI, and p-value
fit_final <- subset(fit_rg0, select = c("term", "estimate", "low", "high", "p.value"))[14:15,]
fit_final$estimate <- round(exp(fit_final$estimate), digits = 2)
fit_final$low <- round(exp(fit_final$low), digits = 2)
fit_final$high <- round(exp(fit_final$high), digits = 2)
fit_final$p.value <- round(fit_final$p.value, digits = 3)
fit_final_IPH <- rbind(fit_final_IPH, fit_final)


####Define the X-axis position
Y_location <- 1

####Define colors
color_point <- "rgb(80, 29, 138)"
color <- "rgba(80, 29, 138, 1)"
color_trans <- "rgba(80, 29, 138, 0.8)"
color_vertex <- c(80, 29, 138)

####Plot points
fit_rg2
fit_rg2$neg_error <- fit_rg2$estimate - fit_rg2$low
fit_rg2$pos_error <- fit_rg2$high - fit_rg2$estimate
fig1_0 <- plot_ly() %>%
  add_markers(data = fit_rg2, 
              x = fit_rg2$X_position, 
              y = Y_location, 
              z = fit_rg2$estimate,
              type = "scatter3d", 
              mode = "markers",
              marker = list(size = 1.5, color = color_point, opacity = 1),
              showlegend = FALSE,
              error_z = list(type = "data", 
                             array = ~pos_error, 
                             arrayminus = ~neg_error, 
                             color = color_point,
                             thickness = 5,#Since this is a 3D plot, thickness has no effect here; consider adding lines manually for better control.
                             width = 5)) %>% 
  layout(
    scene = list(
      xaxis = list(title = "Moderate to vigorous<br>physical activity<br>(unit: METh-week)", 
                   range = c(0, cutvalue_high3), 
                   fixedrange = TRUE, 
                   showticklabels = FALSE, 
                   showgrid = FALSE, 
                   zeroline = FALSE, 
                   showline = FALSE),
      yaxis = list(title = "", 
                   range = c(-0.5, 16.5), 
                   fixedrange = TRUE, 
                   showticklabels = FALSE, 
                   showgrid = FALSE, 
                   zeroline = FALSE, 
                   showline = FALSE), 
      zaxis = list(title = "", 
                   range = c(0, 6), 
                   fixedrange = TRUE, 
                   showticklabels = FALSE, 
                   showgrid = FALSE, 
                   zeroline = FALSE, 
                   showline = FALSE))
  ) %>%
  add_trace(
    x = c(0), y = c(0), z = c(2.2),  # Adjust Z-axis title position
    type = "scatter3d",
    mode = "text",
    text = "Odds ratio of <br> incident IPH",
    textfont = list(size = 12),
    showlegend = FALSE
  )
fig1_0

####Draw the corresponding X-axis
X_axis <- data.frame(x = c(0, cutvalue_high3),
                     y = c(Y_location, Y_location),
                     z = c(0, 0))
#Add an arrow at the end of the X-axis
arrow_x <- cutvalue_high3
arrow_y <- Y_location
arrow_z <- 0
fig1_1 <- fig1_0 %>% 
  add_trace(
    x = X_axis$x,
    y = X_axis$y,
    z = X_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3, dash = "dash"),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(1),  #Arrow points along the X-axis
    v = list(0),
    w = list(0),
    sizemode = "absolute",
    sizeref = 5,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add quantile labels on the X-axis
fig1_1 <- fig1_1 %>%
  add_trace(
    x = cutvalue_1,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_1, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = cutvalue_2,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_2, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig1_1

####Draw the corresponding Z-axis
Z_axis <- data.frame(x = c(0, 0),
                     y = c(Y_location, Y_location),
                     z = c(0, 2))
# Add an arrow at the end of the Z-axis
arrow_x <- 0
arrow_y <- Y_location
arrow_z <- 2
fig1_2 <- fig1_1 %>% 
  add_trace(
    x = Z_axis$x,
    y = Z_axis$y,
    z = Z_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(0),  
    v = list(0),
    w = list(1),#Arrow points along the Z-axis
    sizemode = "absolute",
    sizeref = 0.2,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add labels on the Z-axis
fig1_3 <- fig1_2 %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = 0,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 1,
    type = "scatter3d",
    mode = "text",
    text = 1,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 2,
    type = "scatter3d",
    mode = "text",
    text = 2,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) 
fig1_3


####Add a visual aid plane for the reference line
square_x <- c(0, 0, cutvalue_high3, cutvalue_high3)
square_z <- c(0, 1, 1, 0)
square_y1 <- c(Y_location, Y_location, Y_location, Y_location)
i_set <- c(0,2)
j_set <- c(1,3)
k_set <- c(2,0)
vertex_colors <- matrix(rep(color_vertex, each = 4), nrow = 4)/255
fig1_4 <- fig1_3  %>%
  add_trace(
    x = square_x, 
    z = square_z, 
    y = square_y1, 
    i = i_set, 
    j = j_set, 
    k = k_set,
    type = "mesh3d",
    vertexcolor = vertex_colors,
    opacity = 0.2,
    showscale = FALSE)
fig1_4

####Add auxiliary lines at the quantiles
line_x1 <- c(cutvalue_1, cutvalue_1)
line_z1 <- c(0, 1)
line_y1 <- c(Y_location, Y_location)
line_x2 <- c(cutvalue_2, cutvalue_2)
line_z2 <- c(0, 1)
line_y2 <- c(Y_location, Y_location)
fig1_5 <- fig1_4  %>% 
  add_trace(
    x = line_x1,
    y = line_y1,
    z = line_z1,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE) %>% 
  add_trace(
    x = line_x2,
    y = line_y2,
    z = line_z2,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE)
fig1_5

####Label the grouping method at the end of the X-axis
fig1_6 <- fig1_5 %>% 
  add_trace(
    x = cutvalue_high3,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = "Literature-based<br>cutoff",
    textposition = "right",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig1_6



######################Tertiles -------------------
#total_tertiles
cutvalue <- incident_followup %>%
  group_by(METh_modvig_3) %>%
  summarize(quant0 = quantile(METh_modvig, probs = 0),
            quant25 = quantile(METh_modvig, probs = 0.25),
            quant50 = quantile(METh_modvig, probs = 0.50),
            quant75 = quantile(METh_modvig, probs = 0.75),
            quant100 = quantile(METh_modvig, probs = 1.00))
cutvalue
cutvalue_low1 <- cutvalue$quant0[1]
cutvalue_low2 <- cutvalue$quant0[2]
cutvalue_low3 <- cutvalue$quant0[3]
cutvalue_high1 <- cutvalue$quant100[1]
cutvalue_high2 <- cutvalue$quant100[2]
cutvalue_high3 <- 110
cutvalue_1 <- (cutvalue_low2 + cutvalue_high1)/2
cutvalue_2 <- (cutvalue_low3 + cutvalue_high2)/2
cutvalue_3 <- cutvalue_high3
cutvalue_1
cutvalue_2
cutvalue_3
X_position_1 <- (0 + cutvalue_1) / 2
X_position_2 <- (cutvalue_1 + cutvalue_2) / 2
X_position_3 <- (cutvalue_2 + cutvalue_3) / 2
X_position_1
X_position_2
X_position_3
X_position <- c(X_position_1, X_position_2, X_position_3)

rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + METh_modvig_3
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
fit_rg0 <- summary(pool(rg0))
fit_rg0$low <- fit_rg0$estimate - 1.96 * fit_rg0$std.error
fit_rg0$high <- fit_rg0$estimate + 1.96 * fit_rg0$std.error
fit_rg1 <- subset(fit_rg0, select = c("term", "estimate", "low", "high"))[14:15,]
fit_rg1$estimate <- round(exp(fit_rg1$estimate), digits = 2)
fit_rg1$low <- round(exp(fit_rg1$low), digits = 2)
fit_rg1$high <- round(exp(fit_rg1$high), digits = 2)

contrast_1 <- data.frame(term = "contrast", estimate = 1, low = 1, high = 1)
fit_rg1 <- rbind(contrast_1, fit_rg1)
fit_rg1
fit_rg2 <- cbind(fit_rg1, X_position)
dput(names(fit_rg2))
#Extract OR, CI, and p-value
fit_final <- subset(fit_rg0, select = c("term", "estimate", "low", "high", "p.value"))[14:15,]
fit_final$estimate <- round(exp(fit_final$estimate), digits = 2)
fit_final$low <- round(exp(fit_final$low), digits = 2)
fit_final$high <- round(exp(fit_final$high), digits = 2)
fit_final$p.value <- round(fit_final$p.value, digits = 3)
fit_final_IPH <- rbind(fit_final_IPH, fit_final)

####Define the X-axis position
Y_location <- 5

####Define colors
color_point <- "rgb(28, 128, 65)"
color <- "rgba(28, 128, 65, 1)"
color_trans <- "rgba(28, 128, 65, 0.8)"
color_vertex <- c(28, 128, 65)

####Plot points
fit_rg2
fit_rg2$neg_error <- fit_rg2$estimate - fit_rg2$low
fit_rg2$pos_error <- fit_rg2$high - fit_rg2$estimate
fig2_0 <- fig1_6 %>%
  add_markers(data = fit_rg2, 
              x = fit_rg2$X_position, 
              y = Y_location, 
              z = fit_rg2$estimate,
              type = "scatter3d", 
              mode = "markers",
              marker = list(size = 1.5, color = color_point, opacity = 1),
              showlegend = FALSE,
              error_z = list(type = "data", 
                             array = ~pos_error, 
                             arrayminus = ~neg_error, 
                             color = color_point,
                             thickness = 5,#Since this is a 3D plot, thickness has no effect here; consider adding lines manually for better control.
                             width = 5))
fig2_0

####Draw the corresponding X-axis
X_axis <- data.frame(x = c(0, cutvalue_high3),
                     y = c(Y_location, Y_location),
                     z = c(0, 0))
#Add an arrow at the end of the X-axis
arrow_x <- cutvalue_high3
arrow_y <- Y_location
arrow_z <- 0
fig2_1 <- fig2_0 %>% 
  add_trace(
    x = X_axis$x,
    y = X_axis$y,
    z = X_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3, dash = "dash"),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(1),  #Arrow points along the X-axis
    v = list(0),
    w = list(0),
    sizemode = "absolute",
    sizeref = 5,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add quantile labels on the X-axis
fig2_1 <- fig2_1 %>%
  add_trace(
    x = cutvalue_1,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_1, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = cutvalue_2,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_2, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig2_1

####Draw the corresponding Z-axis
Z_axis <- data.frame(x = c(0, 0),
                     y = c(Y_location, Y_location),
                     z = c(0, 2))
# Add an arrow at the end of the Z-axis
arrow_x <- 0
arrow_y <- Y_location
arrow_z <- 2
fig2_2 <- fig2_1 %>% 
  add_trace(
    x = Z_axis$x,
    y = Z_axis$y,
    z = Z_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(0),  
    v = list(0),
    w = list(1),#Arrow points along the Z-axis
    sizemode = "absolute",
    sizeref = 0.2,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add labels on the Z-axis
fig2_3 <- fig2_2 %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = 0,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 1,
    type = "scatter3d",
    mode = "text",
    text = 1,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 2,
    type = "scatter3d",
    mode = "text",
    text = 2,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)


####Add a visual aid plane for the reference line
square_x <- c(0, 0, cutvalue_high3, cutvalue_high3)
square_z <- c(0, 1, 1, 0)
square_y1 <- c(Y_location, Y_location, Y_location, Y_location)
i_set <- c(0,2)
j_set <- c(1,3)
k_set <- c(2,0)
vertex_colors <- matrix(rep(color_vertex, each = 4), nrow = 4)/255
fig2_4 <- fig2_3  %>%
  add_trace(
    x = square_x, 
    z = square_z, 
    y = square_y1, 
    i = i_set, 
    j = j_set, 
    k = k_set,
    type = "mesh3d",
    vertexcolor = vertex_colors,
    opacity = 0.2,
    showscale = FALSE)
fig2_4

####Add auxiliary lines at the quantiles
line_x1 <- c(cutvalue_1, cutvalue_1)
line_z1 <- c(0, 1)
line_y1 <- c(Y_location, Y_location)
line_x2 <- c(cutvalue_2, cutvalue_2)
line_z2 <- c(0, 1)
line_y2 <- c(Y_location, Y_location)
fig2_5 <- fig2_4  %>% 
  add_trace(
    x = line_x1,
    y = line_y1,
    z = line_z1,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE) %>% 
  add_trace(
    x = line_x2,
    y = line_y2,
    z = line_z2,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE)
fig2_5

####Label the grouping method at the end of the X-axis
fig2_6 <- fig2_5 %>% 
  add_trace(
    x = cutvalue_high3,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = "Tertile",
    textposition = "right",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig2_6




######################Quartiles -------------------
#total_quartiles
cutvalue <- incident_followup %>%
  group_by(METh_modvig_4) %>%
  summarize(quant0 = quantile(METh_modvig, probs = 0),
            quant25 = quantile(METh_modvig, probs = 0.25),
            quant50 = quantile(METh_modvig, probs = 0.50),
            quant75 = quantile(METh_modvig, probs = 0.75),
            quant100 = quantile(METh_modvig, probs = 1.00))
cutvalue
cutvalue_low1 <- cutvalue$quant0[1]
cutvalue_low2 <- cutvalue$quant0[2]
cutvalue_low3 <- cutvalue$quant0[3]
cutvalue_low4 <- cutvalue$quant0[4]
cutvalue_high1 <- cutvalue$quant100[1]
cutvalue_high2 <- cutvalue$quant100[2]
cutvalue_high3 <- cutvalue$quant100[3]
cutvalue_high4 <- 110
cutvalue_1 <- (cutvalue_low2 + cutvalue_high1)/2
cutvalue_2 <- (cutvalue_low3 + cutvalue_high2)/2
cutvalue_3 <- (cutvalue_low4 + cutvalue_high3)/2
cutvalue_4 <- cutvalue_high4
cutvalue_1
cutvalue_2
cutvalue_3
cutvalue_4
X_position_1 <- (0 + cutvalue_1) / 2
X_position_2 <- (cutvalue_1 + cutvalue_2) / 2
X_position_3 <- (cutvalue_2 + cutvalue_3) / 2
X_position_4 <- (cutvalue_3 + cutvalue_4) / 2
X_position_1
X_position_2
X_position_3
X_position_4
X_position <- c(X_position_1, X_position_2, X_position_3, X_position_4)

rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + METh_modvig_4
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
fit_rg0 <- summary(pool(rg0))
fit_rg0$low <- fit_rg0$estimate - 1.96 * fit_rg0$std.error
fit_rg0$high <- fit_rg0$estimate + 1.96 * fit_rg0$std.error
fit_rg1 <- subset(fit_rg0, select = c("term", "estimate", "low", "high"))[14:16,]
fit_rg1$estimate <- round(exp(fit_rg1$estimate), digits = 2)
fit_rg1$low <- round(exp(fit_rg1$low), digits = 2)
fit_rg1$high <- round(exp(fit_rg1$high), digits = 2)

contrast_1 <- data.frame(term = "contrast", estimate = 1, low = 1, high = 1)
fit_rg1 <- rbind(contrast_1, fit_rg1)
fit_rg1
fit_rg2 <- cbind(fit_rg1, X_position)
dput(names(fit_rg2))

#Extract OR, CI, and p-value
fit_final <- subset(fit_rg0, select = c("term", "estimate", "low", "high", "p.value"))[14:16,]
fit_final$estimate <- round(exp(fit_final$estimate), digits = 2)
fit_final$low <- round(exp(fit_final$low), digits = 2)
fit_final$high <- round(exp(fit_final$high), digits = 2)
fit_final$p.value <- round(fit_final$p.value, digits = 3)
fit_final_IPH <- rbind(fit_final_IPH, fit_final)
####Define the X-axis position
Y_location <- 10

####Define colors
color_point <- "rgb(229, 87, 9)"
color <- "rgba(229, 87, 9, 1)"
color_trans <- "rgba(229, 87, 9, 0.8)"
color_vertex <- c(229, 87, 9)

####Plot points
fit_rg2
fit_rg2$neg_error <- fit_rg2$estimate - fit_rg2$low
fit_rg2$pos_error <- fit_rg2$high - fit_rg2$estimate
fig3_0 <- fig2_6 %>%
  add_markers(data = fit_rg2, 
              x = fit_rg2$X_position, 
              y = Y_location, 
              z = fit_rg2$estimate,
              type = "scatter3d", 
              mode = "markers",
              marker = list(size = 1.5, color = color_point, opacity = 1),
              showlegend = FALSE,
              error_z = list(type = "data", 
                             array = ~pos_error, 
                             arrayminus = ~neg_error, 
                             color = color_point,
                             thickness = 5,#Since this is a 3D plot, thickness has no effect here; consider adding lines manually for better control.
                             width = 5))
fig3_0

####Draw the corresponding X-axis
X_axis <- data.frame(x = c(0, cutvalue_high4),
                     y = c(Y_location, Y_location),
                     z = c(0, 0))
#Add an arrow at the end of the X-axis
arrow_x <- cutvalue_high4
arrow_y <- Y_location
arrow_z <- 0
fig3_1 <- fig3_0 %>% 
  add_trace(
    x = X_axis$x,
    y = X_axis$y,
    z = X_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3, dash = "dash"),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(1),  #Arrow points along the X-axis
    v = list(0),
    w = list(0),
    sizemode = "absolute",
    sizeref = 5,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add quantile labels on the X-axis
fig3_1 <- fig3_1 %>%
  add_trace(
    x = cutvalue_1,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_1, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = cutvalue_2,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_2, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = cutvalue_3,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_3, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig3_1

####Draw the corresponding Z-axis
Z_axis <- data.frame(x = c(0, 0),
                     y = c(Y_location, Y_location),
                     z = c(0, 2))
# Add an arrow at the end of the Z-axis
arrow_x <- 0
arrow_y <- Y_location
arrow_z <- 2
fig3_2 <- fig3_1 %>% 
  add_trace(
    x = Z_axis$x,
    y = Z_axis$y,
    z = Z_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(0),  
    v = list(0),
    w = list(1),#Arrow points along the Z-axis
    sizemode = "absolute",
    sizeref = 0.2,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add labels on the Z-axis
fig3_3 <- fig3_2 %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = 0,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 1,
    type = "scatter3d",
    mode = "text",
    text = 1,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 2,
    type = "scatter3d",
    mode = "text",
    text = 2,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig3_3


####Add a visual aid plane for the reference line
square_x <- c(0, 0, cutvalue_high4, cutvalue_high4)
square_z <- c(0, 1, 1, 0)
square_y1 <- c(Y_location, Y_location, Y_location, Y_location)
i_set <- c(0,2)
j_set <- c(1,3)
k_set <- c(2,0)
vertex_colors <- matrix(rep(color_vertex, each = 4), nrow = 4)/255
fig3_4 <- fig3_3  %>%
  add_trace(
    x = square_x, 
    z = square_z, 
    y = square_y1, 
    i = i_set, 
    j = j_set, 
    k = k_set,
    type = "mesh3d",
    vertexcolor = vertex_colors,
    opacity = 0.2,
    showscale = FALSE)
fig3_4

####Add auxiliary lines at the quantiles
line_x1 <- c(cutvalue_1, cutvalue_1)
line_z1 <- c(0, 1)
line_y1 <- c(Y_location, Y_location)
line_x2 <- c(cutvalue_2, cutvalue_2)
line_z2 <- c(0, 1)
line_y2 <- c(Y_location, Y_location)
line_x3 <- c(cutvalue_3, cutvalue_3)
line_z3 <- c(0, 1)
line_y3 <- c(Y_location, Y_location)
fig3_5 <- fig3_4  %>% 
  add_trace(
    x = line_x1,
    y = line_y1,
    z = line_z1,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE) %>% 
  add_trace(
    x = line_x2,
    y = line_y2,
    z = line_z2,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE) %>% 
  add_trace(
    x = line_x3,
    y = line_y3,
    z = line_z3,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE)
fig3_5

####Label the grouping method at the end of the X-axis
fig3_6 <- fig3_5 %>% 
  add_trace(
    x = cutvalue_high4,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = "Quartile",
    textposition = "right",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig3_6





######################Quintiles -------------------
#total_quintiles
cutvalue <- incident_followup %>%
  group_by(METh_modvig_5) %>%
  summarize(quant0 = quantile(METh_modvig, probs = 0),
            quant25 = quantile(METh_modvig, probs = 0.25),
            quant50 = quantile(METh_modvig, probs = 0.50),
            quant75 = quantile(METh_modvig, probs = 0.75),
            quant100 = quantile(METh_modvig, probs = 1.00))
cutvalue
cutvalue_low1 <- cutvalue$quant0[1]
cutvalue_low2 <- cutvalue$quant0[2]
cutvalue_low3 <- cutvalue$quant0[3]
cutvalue_low4 <- cutvalue$quant0[4]
cutvalue_low5 <- cutvalue$quant0[5]
cutvalue_high1 <- cutvalue$quant100[1]
cutvalue_high2 <- cutvalue$quant100[2]
cutvalue_high3 <- cutvalue$quant100[3]
cutvalue_high4 <- cutvalue$quant100[4]
cutvalue_high5 <- 110
cutvalue_1 <- (cutvalue_low2 + cutvalue_high1)/2
cutvalue_2 <- (cutvalue_low3 + cutvalue_high2)/2
cutvalue_3 <- (cutvalue_low4 + cutvalue_high3)/2
cutvalue_4 <- (cutvalue_low5 + cutvalue_high4)/2
cutvalue_5 <- cutvalue_high5
cutvalue_1
cutvalue_2
cutvalue_3
cutvalue_4
cutvalue_5
X_position_1 <- (0 + cutvalue_1) / 2
X_position_2 <- (cutvalue_1 + cutvalue_2) / 2
X_position_3 <- (cutvalue_2 + cutvalue_3) / 2
X_position_4 <- (cutvalue_3 + cutvalue_4) / 2
X_position_5 <- (cutvalue_5 + cutvalue_4) / 2
X_position_1
X_position_2
X_position_3
X_position_4
X_position_5

X_position <- c(X_position_1, X_position_2, X_position_3, X_position_4, X_position_5)

rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + METh_modvig_5
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
fit_rg0 <- summary(pool(rg0))
fit_rg0$low <- fit_rg0$estimate - 1.96 * fit_rg0$std.error
fit_rg0$high <- fit_rg0$estimate + 1.96 * fit_rg0$std.error
fit_rg1 <- subset(fit_rg0, select = c("term", "estimate", "low", "high"))[14:17,]
fit_rg1$estimate <- round(exp(fit_rg1$estimate), digits = 2)
fit_rg1$low <- round(exp(fit_rg1$low), digits = 2)
fit_rg1$high <- round(exp(fit_rg1$high), digits = 2)

contrast_1 <- data.frame(term = "contrast", estimate = 1, low = 1, high = 1)
fit_rg1 <- rbind(contrast_1, fit_rg1)
fit_rg1
fit_rg2 <- cbind(fit_rg1, X_position)
dput(names(fit_rg2))
fit_rg2
#Extract OR, CI, and p-value
fit_final <- subset(fit_rg0, select = c("term", "estimate", "low", "high", "p.value"))[14:17,]
fit_final$estimate <- round(exp(fit_final$estimate), digits = 2)
fit_final$low <- round(exp(fit_final$low), digits = 2)
fit_final$high <- round(exp(fit_final$high), digits = 2)
fit_final$p.value <- round(fit_final$p.value, digits = 3)
fit_final_IPH <- rbind(fit_final_IPH, fit_final)
####Define the X-axis position
Y_location <- 16

####Define colors
color_point <- "rgb(43, 106, 153)"
color <- "rgba(43, 106, 153, 1)"
color_trans <- "rgba(43, 106, 153, 0.8)"
color_vertex <- c(43, 106, 153)

####Plot points
fit_rg2
fit_rg2$neg_error <- fit_rg2$estimate - fit_rg2$low
fit_rg2$pos_error <- fit_rg2$high - fit_rg2$estimate
fig4_0 <- fig3_6 %>%
  add_markers(data = fit_rg2, 
              x = fit_rg2$X_position, 
              y = Y_location, 
              z = fit_rg2$estimate,
              type = "scatter3d", 
              mode = "markers",
              marker = list(size = 1.5, color = color_point, opacity = 1),
              showlegend = FALSE,
              error_z = list(type = "data", 
                             array = ~pos_error, 
                             arrayminus = ~neg_error, 
                             color = color_point,
                             thickness = 5,#Since this is a 3D plot, thickness has no effect here; consider adding lines manually for better control.
                             width = 5))
fig4_0

####Draw the corresponding X-axis
X_axis <- data.frame(x = c(0, cutvalue_high5),
                     y = c(Y_location, Y_location),
                     z = c(0, 0))
#Add an arrow at the end of the X-axis
arrow_x <- cutvalue_high5
arrow_y <- Y_location
arrow_z <- 0
fig4_1 <- fig4_0 %>% 
  add_trace(
    x = X_axis$x,
    y = X_axis$y,
    z = X_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3, dash = "dash"),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(1),  #Arrow points along the X-axis
    v = list(0),
    w = list(0),
    sizemode = "absolute",
    sizeref = 5,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add quantile labels on the X-axis
fig4_1 <- fig4_1 %>%
  add_trace(
    x = cutvalue_1,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_1, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = cutvalue_2,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_2, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = cutvalue_3,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_3, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = cutvalue_4,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_4, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig4_1

####Draw the corresponding Z-axis
Z_axis <- data.frame(x = c(0, 0),
                     y = c(Y_location, Y_location),
                     z = c(0, 2))
# Add an arrow at the end of the Z-axis
arrow_x <- 0
arrow_y <- Y_location
arrow_z <- 2
fig4_2 <- fig4_1 %>% 
  add_trace(
    x = Z_axis$x,
    y = Z_axis$y,
    z = Z_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(0),  
    v = list(0),
    w = list(1),#Arrow points along the Z-axis
    sizemode = "absolute",
    sizeref = 0.2,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add labels on the Z-axis
fig4_3 <- fig4_2 %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = 0,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 1,
    type = "scatter3d",
    mode = "text",
    text = 1,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 2,
    type = "scatter3d",
    mode = "text",
    text = 2,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig4_3


####Add a visual aid plane for the reference line
square_x <- c(0, 0, cutvalue_high5, cutvalue_high5)
square_z <- c(0, 1, 1, 0)
square_y1 <- c(Y_location, Y_location, Y_location, Y_location)
i_set <- c(0,2)
j_set <- c(1,3)
k_set <- c(2,0)
vertex_colors <- matrix(rep(color_vertex, each = 4), nrow = 4)/255
fig4_4 <- fig4_3  %>%
  add_trace(
    x = square_x, 
    z = square_z, 
    y = square_y1, 
    i = i_set, 
    j = j_set, 
    k = k_set,
    type = "mesh3d",
    vertexcolor = vertex_colors,
    opacity = 0.2,
    showscale = FALSE)
fig4_4

####Add auxiliary lines at the quantiles
line_x1 <- c(cutvalue_1, cutvalue_1)
line_z1 <- c(0, 1)
line_y1 <- c(Y_location, Y_location)
line_x2 <- c(cutvalue_2, cutvalue_2)
line_z2 <- c(0, 1)
line_y2 <- c(Y_location, Y_location)
line_x3 <- c(cutvalue_3, cutvalue_3)
line_z3 <- c(0, 1)
line_y3 <- c(Y_location, Y_location)
line_x4 <- c(cutvalue_4, cutvalue_4)
line_z4 <- c(0, 1)
line_y4 <- c(Y_location, Y_location)
fig4_5 <- fig4_4  %>% 
  add_trace(
    x = line_x1,
    y = line_y1,
    z = line_z1,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE) %>% 
  add_trace(
    x = line_x2,
    y = line_y2,
    z = line_z2,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE) %>% 
  add_trace(
    x = line_x3,
    y = line_y3,
    z = line_z3,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE) %>% 
  add_trace(
    x = line_x4,
    y = line_y4,
    z = line_z4,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE)
fig4_5

####Label the grouping method at the end of the X-axis
fig4_6 <- fig4_5 %>% 
  add_trace(
    x = cutvalue_high5,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = "Quintile",
    textposition = "right",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig4_6
fig4_7 <- fig4_6 %>% layout(
  scene = list(camera = list(eye = list(x = 0.59097400646812713, y = -1.878239943166165, z = 1.309212875732691),
                             aspectmode = "manual",
                             aspectratio = list(x = 1, y = 2, z = 1)  # Y-axis is stretched
  ))
)
###Add the PA guideline-recommended threshold
fig4_7 <- fig4_7 %>% 
  add_text(x = 15, y = 1, z = 0,
           text = '*', textposition = 'middle center',
           textfont = list(size = 20, color = 'red'))%>% 
  add_text(x = 15, y = 5, z = 0,
           text = '*', textposition = 'middle center',
           textfont = list(size = 20, color = 'red'))%>% 
  add_text(x = 15, y = 10, z = 0,
           text = '*', textposition = 'middle center',
           textfont = list(size = 20, color = 'red'))%>% 
  add_text(x = 15, y = 16, z = 0,
           text = '*', textposition = 'middle center',
           textfont = list(size = 20, color = 'red'))%>% 
  add_text(x = 15, y = -0.5, z = 0,
           text = '     15  <br>Guideline<br>Recomm\'s', textposition = 'middle center',
           textfont = list(size = 12, color = 'red')) %>%
  layout(showlegend = FALSE)

fig4_7

saveNetwork(fig4_7, "plot.html")
webshot2::webshot("plot.html", 
                  "IPH_modvig.png",
                  vwidth = 1200, vheight = 800, zoom = 3)


#Crop the PNG so it can be stitched directly later
img_imager <- load.image("IPH_modvig.png")  # Load image
plot(img_imager)
dim(img_imager)
cropped_img <- imsub(img_imager, x %in% 1050:2700, y %in% 720:1950)
plot(cropped_img)
imager::save.image(cropped_img, "IPH_modvig.png")





#####################################################Plotting: Moderate PA------------------------------
######################Tertiles -------------------
cutvalue <- incident_followup %>%
  group_by(METh_moderate_3) %>%
  summarize(quant0 = quantile(METh_moderate, probs = 0),
            quant25 = quantile(METh_moderate, probs = 0.25),
            quant50 = quantile(METh_moderate, probs = 0.50),
            quant75 = quantile(METh_moderate, probs = 0.75),
            quant100 = quantile(METh_moderate, probs = 1.00))
cutvalue
cutvalue_low1 <- cutvalue$quant0[1]
cutvalue_low2 <- cutvalue$quant0[2]
cutvalue_low3 <- cutvalue$quant0[3]
cutvalue_high1 <- cutvalue$quant100[1]
cutvalue_high2 <- cutvalue$quant100[2]
cutvalue_high3 <- 85
cutvalue_1 <- (cutvalue_low2 + cutvalue_high1)/2
cutvalue_2 <- (cutvalue_low3 + cutvalue_high2)/2
cutvalue_3 <- cutvalue_high3
cutvalue_1
cutvalue_2
cutvalue_3
X_position_1 <- (0 + cutvalue_1) / 2
X_position_2 <- (cutvalue_1 + cutvalue_2) / 2
X_position_3 <- (cutvalue_2 + cutvalue_3) / 2
X_position_1
X_position_2
X_position_3
X_position <- c(X_position_1, X_position_2, X_position_3)

rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + METh_moderate_3
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
fit_rg0 <- summary(pool(rg0))
fit_rg0$low <- fit_rg0$estimate - 1.96 * fit_rg0$std.error
fit_rg0$high <- fit_rg0$estimate + 1.96 * fit_rg0$std.error
fit_rg1 <- subset(fit_rg0, select = c("term", "estimate", "low", "high"))[14:15,]
fit_rg1$estimate <- round(exp(fit_rg1$estimate), digits = 2)
fit_rg1$low <- round(exp(fit_rg1$low), digits = 2)
fit_rg1$high <- round(exp(fit_rg1$high), digits = 2)

contrast_1 <- data.frame(term = "contrast", estimate = 1, low = 1, high = 1)
fit_rg1 <- rbind(contrast_1, fit_rg1)
fit_rg1
fit_rg2 <- cbind(fit_rg1, X_position)
dput(names(fit_rg2))
#Extract OR, CI, and p-value
fit_final <- subset(fit_rg0, select = c("term", "estimate", "low", "high", "p.value"))[14:15,]
fit_final$estimate <- round(exp(fit_final$estimate), digits = 2)
fit_final$low <- round(exp(fit_final$low), digits = 2)
fit_final$high <- round(exp(fit_final$high), digits = 2)
fit_final$p.value <- round(fit_final$p.value, digits = 3)
fit_final_IPH <- rbind(fit_final_IPH, fit_final)

####Define the X-axis position
Y_location <- 1

####Define colors
color_point <- "rgb(28, 128, 65)"
color <- "rgba(28, 128, 65, 1)"
color_trans <- "rgba(28, 128, 65, 0.8)"
color_vertex <- c(28, 128, 65)

####Plot points
fit_rg2
fit_rg2$neg_error <- fit_rg2$estimate - fit_rg2$low
fit_rg2$pos_error <- fit_rg2$high - fit_rg2$estimate
fig1_0 <- plot_ly() %>%
  add_markers(data = fit_rg2, 
              x = fit_rg2$X_position, 
              y = Y_location, 
              z = fit_rg2$estimate,
              type = "scatter3d", 
              mode = "markers",
              marker = list(size = 1.5, color = color_point, opacity = 1),
              showlegend = FALSE,
              error_z = list(type = "data", 
                             array = ~pos_error, 
                             arrayminus = ~neg_error, 
                             color = color_point,
                             thickness = 5,#Since this is a 3D plot, thickness has no effect here; consider adding lines manually for better control.
                             width = 5)) %>% 
  layout(
    scene = list(
      xaxis = list(title = "Moderate physical activity<br>(unit: METh-week)", 
                   range = c(0, cutvalue_high3), 
                   fixedrange = TRUE, 
                   showticklabels = FALSE, 
                   showgrid = FALSE, 
                   zeroline = FALSE, 
                   showline = FALSE),
      yaxis = list(title = "", 
                   range = c(-0.5, 11.5), 
                   fixedrange = TRUE, 
                   showticklabels = FALSE, 
                   showgrid = FALSE, 
                   zeroline = FALSE, 
                   showline = FALSE), 
      zaxis = list(title = "", 
                   range = c(0, 6), 
                   fixedrange = TRUE, 
                   showticklabels = FALSE, 
                   showgrid = FALSE, 
                   zeroline = FALSE, 
                   showline = FALSE))
  ) %>%
  add_trace(
    x = c(0), y = c(0), z = c(2.2),  # Adjust Z-axis title position
    type = "scatter3d",
    mode = "text",
    text = "Odds ratio of <br> incident IPH",
    textfont = list(size = 12),
    showlegend = FALSE
  )
fig1_0

####Draw the corresponding X-axis
X_axis <- data.frame(x = c(0, cutvalue_high3),
                     y = c(Y_location, Y_location),
                     z = c(0, 0))
#Add an arrow at the end of the X-axis
arrow_x <- cutvalue_high3
arrow_y <- Y_location
arrow_z <- 0
fig2_1 <- fig1_0 %>% 
  add_trace(
    x = X_axis$x,
    y = X_axis$y,
    z = X_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3, dash = "dash"),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(1),  #Arrow points along the X-axis
    v = list(0),
    w = list(0),
    sizemode = "absolute",
    sizeref = 5,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add quantile labels on the X-axis
fig2_1 <- fig2_1 %>%
  add_trace(
    x = cutvalue_1,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_1, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = cutvalue_2,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_2, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig2_1

####Draw the corresponding Z-axis
Z_axis <- data.frame(x = c(0, 0),
                     y = c(Y_location, Y_location),
                     z = c(0, 2))
# Add an arrow at the end of the Z-axis
arrow_x <- 0
arrow_y <- Y_location
arrow_z <- 2
fig2_2 <- fig2_1 %>% 
  add_trace(
    x = Z_axis$x,
    y = Z_axis$y,
    z = Z_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(0),  
    v = list(0),
    w = list(1),#Arrow points along the Z-axis
    sizemode = "absolute",
    sizeref = 0.35,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add labels on the Z-axis
fig2_3 <- fig2_2 %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = 0,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 1,
    type = "scatter3d",
    mode = "text",
    text = 1,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 2,
    type = "scatter3d",
    mode = "text",
    text = 2,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) 
fig2_3


####Add a visual aid plane for the reference line
square_x <- c(0, 0, cutvalue_high3, cutvalue_high3)
square_z <- c(0, 1, 1, 0)
square_y1 <- c(Y_location, Y_location, Y_location, Y_location)
i_set <- c(0,2)
j_set <- c(1,3)
k_set <- c(2,0)
vertex_colors <- matrix(rep(color_vertex, each = 4), nrow = 4)/255
fig2_4 <- fig2_3  %>%
  add_trace(
    x = square_x, 
    z = square_z, 
    y = square_y1, 
    i = i_set, 
    j = j_set, 
    k = k_set,
    type = "mesh3d",
    vertexcolor = vertex_colors,
    opacity = 0.2,
    showscale = FALSE)
fig2_4

####Add auxiliary lines at the quantiles
line_x1 <- c(cutvalue_1, cutvalue_1)
line_z1 <- c(0, 1)
line_y1 <- c(Y_location, Y_location)
line_x2 <- c(cutvalue_2, cutvalue_2)
line_z2 <- c(0, 1)
line_y2 <- c(Y_location, Y_location)
fig2_5 <- fig2_4  %>% 
  add_trace(
    x = line_x1,
    y = line_y1,
    z = line_z1,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE) %>% 
  add_trace(
    x = line_x2,
    y = line_y2,
    z = line_z2,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE)
fig2_5

####Label the grouping method at the end of the X-axis
fig2_6 <- fig2_5 %>% 
  add_trace(
    x = cutvalue_high3,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = "Tertile",
    textposition = "right",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig2_6




######################Quartiles -------------------
#total_quartiles
cutvalue <- incident_followup %>%
  group_by(METh_moderate_4) %>%
  summarize(quant0 = quantile(METh_moderate, probs = 0),
            quant25 = quantile(METh_moderate, probs = 0.25),
            quant50 = quantile(METh_moderate, probs = 0.50),
            quant75 = quantile(METh_moderate, probs = 0.75),
            quant100 = quantile(METh_moderate, probs = 1.00))
cutvalue
cutvalue_low1 <- cutvalue$quant0[1]
cutvalue_low2 <- cutvalue$quant0[2]
cutvalue_low3 <- cutvalue$quant0[3]
cutvalue_low4 <- cutvalue$quant0[4]
cutvalue_high1 <- cutvalue$quant100[1]
cutvalue_high2 <- cutvalue$quant100[2]
cutvalue_high3 <- cutvalue$quant100[3]
cutvalue_high4 <- 85
cutvalue_1 <- (cutvalue_low2 + cutvalue_high1)/2
cutvalue_2 <- (cutvalue_low3 + cutvalue_high2)/2
cutvalue_3 <- (cutvalue_low4 + cutvalue_high3)/2
cutvalue_4 <- cutvalue_high4
cutvalue_1
cutvalue_2
cutvalue_3
cutvalue_4
X_position_1 <- (0 + cutvalue_1) / 2
X_position_2 <- (cutvalue_1 + cutvalue_2) / 2
X_position_3 <- (cutvalue_2 + cutvalue_3) / 2
X_position_4 <- (cutvalue_3 + cutvalue_4) / 2
X_position_1
X_position_2
X_position_3
X_position_4
X_position <- c(X_position_1, X_position_2, X_position_3, X_position_4)

rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + METh_moderate_4
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
fit_rg0 <- summary(pool(rg0))
fit_rg0$low <- fit_rg0$estimate - 1.96 * fit_rg0$std.error
fit_rg0$high <- fit_rg0$estimate + 1.96 * fit_rg0$std.error
fit_rg1 <- subset(fit_rg0, select = c("term", "estimate", "low", "high"))[14:16,]
fit_rg1$estimate <- round(exp(fit_rg1$estimate), digits = 2)
fit_rg1$low <- round(exp(fit_rg1$low), digits = 2)
fit_rg1$high <- round(exp(fit_rg1$high), digits = 2)

contrast_1 <- data.frame(term = "contrast", estimate = 1, low = 1, high = 1)
fit_rg1 <- rbind(contrast_1, fit_rg1)
fit_rg1
fit_rg2 <- cbind(fit_rg1, X_position)
dput(names(fit_rg2))
#Extract OR, CI, and p-value
fit_final <- subset(fit_rg0, select = c("term", "estimate", "low", "high", "p.value"))[14:16,]
fit_final$estimate <- round(exp(fit_final$estimate), digits = 2)
fit_final$low <- round(exp(fit_final$low), digits = 2)
fit_final$high <- round(exp(fit_final$high), digits = 2)
fit_final$p.value <- round(fit_final$p.value, digits = 3)
fit_final_IPH <- rbind(fit_final_IPH, fit_final)

####Define the X-axis position
Y_location <- 5

####Define colors
color_point <- "rgb(229, 87, 9)"
color <- "rgba(229, 87, 9, 1)"
color_trans <- "rgba(229, 87, 9, 0.8)"
color_vertex <- c(229, 87, 9)

####Plot points
fit_rg2
fit_rg2$neg_error <- fit_rg2$estimate - fit_rg2$low
fit_rg2$pos_error <- fit_rg2$high - fit_rg2$estimate
fig3_0 <- fig2_6 %>%
  add_markers(data = fit_rg2, 
              x = fit_rg2$X_position, 
              y = Y_location, 
              z = fit_rg2$estimate,
              type = "scatter3d", 
              mode = "markers",
              marker = list(size = 1.5, color = color_point, opacity = 1),
              showlegend = FALSE,
              error_z = list(type = "data", 
                             array = ~pos_error, 
                             arrayminus = ~neg_error, 
                             color = color_point,
                             thickness = 5,#Since this is a 3D plot, thickness has no effect here; consider adding lines manually for better control.
                             width = 5))
fig3_0

####Draw the corresponding X-axis
X_axis <- data.frame(x = c(0, cutvalue_high4),
                     y = c(Y_location, Y_location),
                     z = c(0, 0))
#Add an arrow at the end of the X-axis
arrow_x <- cutvalue_high4
arrow_y <- Y_location
arrow_z <- 0
fig3_1 <- fig3_0 %>% 
  add_trace(
    x = X_axis$x,
    y = X_axis$y,
    z = X_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3, dash = "dash"),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(1),  #Arrow points along the X-axis
    v = list(0),
    w = list(0),
    sizemode = "absolute",
    sizeref = 5,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add quantile labels on the X-axis
fig3_1 <- fig3_1 %>%
  add_trace(
    x = cutvalue_1,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_1, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = cutvalue_2,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_2, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = cutvalue_3,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_3, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig3_1

####Draw the corresponding Z-axis
Z_axis <- data.frame(x = c(0, 0),
                     y = c(Y_location, Y_location),
                     z = c(0, 2))
# Add an arrow at the end of the Z-axis
arrow_x <- 0
arrow_y <- Y_location
arrow_z <- 2
fig3_2 <- fig3_1 %>% 
  add_trace(
    x = Z_axis$x,
    y = Z_axis$y,
    z = Z_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(0),  
    v = list(0),
    w = list(1),#Arrow points along the Z-axis
    sizemode = "absolute",
    sizeref = 0.35,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add labels on the Z-axis
fig3_3 <- fig3_2 %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = 0,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 1,
    type = "scatter3d",
    mode = "text",
    text = 1,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 2,
    type = "scatter3d",
    mode = "text",
    text = 2,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) 
fig3_3


####Add a visual aid plane for the reference line
square_x <- c(0, 0, cutvalue_high4, cutvalue_high4)
square_z <- c(0, 1, 1, 0)
square_y1 <- c(Y_location, Y_location, Y_location, Y_location)
i_set <- c(0,2)
j_set <- c(1,3)
k_set <- c(2,0)
vertex_colors <- matrix(rep(color_vertex, each = 4), nrow = 4)/255
fig3_4 <- fig3_3  %>%
  add_trace(
    x = square_x, 
    z = square_z, 
    y = square_y1, 
    i = i_set, 
    j = j_set, 
    k = k_set,
    type = "mesh3d",
    vertexcolor = vertex_colors,
    opacity = 0.2,
    showscale = FALSE)
fig3_4

####Add auxiliary lines at the quantiles
line_x1 <- c(cutvalue_1, cutvalue_1)
line_z1 <- c(0, 1)
line_y1 <- c(Y_location, Y_location)
line_x2 <- c(cutvalue_2, cutvalue_2)
line_z2 <- c(0, 1)
line_y2 <- c(Y_location, Y_location)
line_x3 <- c(cutvalue_3, cutvalue_3)
line_z3 <- c(0, 1)
line_y3 <- c(Y_location, Y_location)
fig3_5 <- fig3_4  %>% 
  add_trace(
    x = line_x1,
    y = line_y1,
    z = line_z1,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE) %>% 
  add_trace(
    x = line_x2,
    y = line_y2,
    z = line_z2,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE) %>% 
  add_trace(
    x = line_x3,
    y = line_y3,
    z = line_z3,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE)
fig3_5

####Label the grouping method at the end of the X-axis
fig3_6 <- fig3_5 %>% 
  add_trace(
    x = cutvalue_high4,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = "Quartile",
    textposition = "right",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig3_6





######################Quintiles -------------------
#total_quintiles
cutvalue <- incident_followup %>%
  group_by(METh_moderate_5) %>%
  summarize(quant0 = quantile(METh_moderate, probs = 0),
            quant25 = quantile(METh_moderate, probs = 0.25),
            quant50 = quantile(METh_moderate, probs = 0.50),
            quant75 = quantile(METh_moderate, probs = 0.75),
            quant100 = quantile(METh_moderate, probs = 1.00))
cutvalue
cutvalue_low1 <- cutvalue$quant0[1]
cutvalue_low2 <- cutvalue$quant0[2]
cutvalue_low3 <- cutvalue$quant0[3]
cutvalue_low4 <- cutvalue$quant0[4]
cutvalue_low5 <- cutvalue$quant0[5]
cutvalue_high1 <- cutvalue$quant100[1]
cutvalue_high2 <- cutvalue$quant100[2]
cutvalue_high3 <- cutvalue$quant100[3]
cutvalue_high4 <- cutvalue$quant100[4]
cutvalue_high5 <- 85
cutvalue_1 <- (cutvalue_low2 + cutvalue_high1)/2
cutvalue_2 <- (cutvalue_low3 + cutvalue_high2)/2
cutvalue_3 <- (cutvalue_low4 + cutvalue_high3)/2
cutvalue_4 <- (cutvalue_low5 + cutvalue_high4)/2
cutvalue_5 <- cutvalue_high5
cutvalue_1
cutvalue_2
cutvalue_3
cutvalue_4
cutvalue_5
X_position_1 <- (0 + cutvalue_1) / 2
X_position_2 <- (cutvalue_1 + cutvalue_2) / 2
X_position_3 <- (cutvalue_2 + cutvalue_3) / 2
X_position_4 <- (cutvalue_3 + cutvalue_4) / 2
X_position_5 <- (cutvalue_5 + cutvalue_4) / 2
X_position_1
X_position_2
X_position_3
X_position_4
X_position_5

X_position <- c(X_position_1, X_position_2, X_position_3, X_position_4, X_position_5)

rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + METh_moderate_5
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
fit_rg0 <- summary(pool(rg0))
fit_rg0$low <- fit_rg0$estimate - 1.96 * fit_rg0$std.error
fit_rg0$high <- fit_rg0$estimate + 1.96 * fit_rg0$std.error
fit_rg1 <- subset(fit_rg0, select = c("term", "estimate", "low", "high"))[14:17,]
fit_rg1$estimate <- round(exp(fit_rg1$estimate), digits = 2)
fit_rg1$low <- round(exp(fit_rg1$low), digits = 2)
fit_rg1$high <- round(exp(fit_rg1$high), digits = 2)

contrast_1 <- data.frame(term = "contrast", estimate = 1, low = 1, high = 1)
fit_rg1 <- rbind(contrast_1, fit_rg1)
fit_rg1
fit_rg2 <- cbind(fit_rg1, X_position)
dput(names(fit_rg2))
fit_rg2
#Extract OR, CI, and p-value
fit_final <- subset(fit_rg0, select = c("term", "estimate", "low", "high", "p.value"))[14:17,]
fit_final$estimate <- round(exp(fit_final$estimate), digits = 2)
fit_final$low <- round(exp(fit_final$low), digits = 2)
fit_final$high <- round(exp(fit_final$high), digits = 2)
fit_final$p.value <- round(fit_final$p.value, digits = 3)
fit_final_IPH <- rbind(fit_final_IPH, fit_final)
####Define the X-axis position
Y_location <- 10

####Define colors
color_point <- "rgb(43, 106, 153)"
color <- "rgba(43, 106, 153, 1)"
color_trans <- "rgba(43, 106, 153, 0.8)"
color_vertex <- c(43, 106, 153)

####Plot points
fit_rg2
fit_rg2$neg_error <- fit_rg2$estimate - fit_rg2$low
fit_rg2$pos_error <- fit_rg2$high - fit_rg2$estimate
fig4_0 <- fig3_6 %>%
  add_markers(data = fit_rg2, 
              x = fit_rg2$X_position, 
              y = Y_location, 
              z = fit_rg2$estimate,
              type = "scatter3d", 
              mode = "markers",
              marker = list(size = 1.5, color = color_point, opacity = 1),
              showlegend = FALSE,
              error_z = list(type = "data", 
                             array = ~pos_error, 
                             arrayminus = ~neg_error, 
                             color = color_point,
                             thickness = 5,#Since this is a 3D plot, thickness has no effect here; consider adding lines manually for better control.
                             width = 5))
fig4_0

####Draw the corresponding X-axis
X_axis <- data.frame(x = c(0, cutvalue_high5),
                     y = c(Y_location, Y_location),
                     z = c(0, 0))
#Add an arrow at the end of the X-axis
arrow_x <- cutvalue_high5
arrow_y <- Y_location
arrow_z <- 0
fig4_1 <- fig4_0 %>% 
  add_trace(
    x = X_axis$x,
    y = X_axis$y,
    z = X_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3, dash = "dash"),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(1),  #Arrow points along the X-axis
    v = list(0),
    w = list(0),
    sizemode = "absolute",
    sizeref = 5,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add quantile labels on the X-axis
fig4_1 <- fig4_1 %>%
  add_trace(
    x = cutvalue_1,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_1, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = cutvalue_2,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_2, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = cutvalue_3,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_3, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = cutvalue_4,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_4, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig4_1

####Draw the corresponding Z-axis
Z_axis <- data.frame(x = c(0, 0),
                     y = c(Y_location, Y_location),
                     z = c(0, 2))
# Add an arrow at the end of the Z-axis
arrow_x <- 0
arrow_y <- Y_location
arrow_z <- 2
fig4_2 <- fig4_1 %>% 
  add_trace(
    x = Z_axis$x,
    y = Z_axis$y,
    z = Z_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(0),  
    v = list(0),
    w = list(1),#Arrow points along the Z-axis
    sizemode = "absolute",
    sizeref = 0.35,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add labels on the Z-axis
fig4_3 <- fig4_2 %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = 0,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 1,
    type = "scatter3d",
    mode = "text",
    text = 1,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 2,
    type = "scatter3d",
    mode = "text",
    text = 2,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) 


####Add a visual aid plane for the reference line
square_x <- c(0, 0, cutvalue_high5, cutvalue_high5)
square_z <- c(0, 1, 1, 0)
square_y1 <- c(Y_location, Y_location, Y_location, Y_location)
i_set <- c(0,2)
j_set <- c(1,3)
k_set <- c(2,0)
vertex_colors <- matrix(rep(color_vertex, each = 4), nrow = 4)/255
fig4_4 <- fig4_3  %>%
  add_trace(
    x = square_x, 
    z = square_z, 
    y = square_y1, 
    i = i_set, 
    j = j_set, 
    k = k_set,
    type = "mesh3d",
    vertexcolor = vertex_colors,
    opacity = 0.2,
    showscale = FALSE)
fig4_4

####Add auxiliary lines at the quantiles
line_x1 <- c(cutvalue_1, cutvalue_1)
line_z1 <- c(0, 1)
line_y1 <- c(Y_location, Y_location)
line_x2 <- c(cutvalue_2, cutvalue_2)
line_z2 <- c(0, 1)
line_y2 <- c(Y_location, Y_location)
line_x3 <- c(cutvalue_3, cutvalue_3)
line_z3 <- c(0, 1)
line_y3 <- c(Y_location, Y_location)
line_x4 <- c(cutvalue_4, cutvalue_4)
line_z4 <- c(0, 1)
line_y4 <- c(Y_location, Y_location)
fig4_5 <- fig4_4  %>% 
  add_trace(
    x = line_x1,
    y = line_y1,
    z = line_z1,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE) %>% 
  add_trace(
    x = line_x2,
    y = line_y2,
    z = line_z2,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE) %>% 
  add_trace(
    x = line_x3,
    y = line_y3,
    z = line_z3,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE) %>% 
  add_trace(
    x = line_x4,
    y = line_y4,
    z = line_z4,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE)
fig4_5

####Label the grouping method at the end of the X-axis
fig4_6 <- fig4_5 %>% 
  add_trace(
    x = cutvalue_high5,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = "Quintile",
    textposition = "right",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig4_6
fig4_7 <- fig4_6 %>% layout(
  scene = list(camera = list(eye = list(x = 0.59097400646812713, y = -1.878239943166165, z = 1.309212875732691),
                             aspectmode = "manual",
                             aspectratio = list(x = 1, y = 2, z = 1)  # Y-axis is stretched
  ))
)
###Add the PA guideline-recommended threshold
fig4_7 <- fig4_7 %>% 
  add_text(x = 15, y = 1, z = 0,
           text = '*', textposition = 'middle center',
           textfont = list(size = 20, color = 'red'))%>% 
  add_text(x = 15, y = 5, z = 0,
           text = '*', textposition = 'middle center',
           textfont = list(size = 20, color = 'red'))%>% 
  add_text(x = 15, y = 10, z = 0,
           text = '*', textposition = 'middle center',
           textfont = list(size = 20, color = 'red'))%>% 
  add_text(x = 15, y = 16, z = 0,
           text = '*', textposition = 'middle center',
           textfont = list(size = 20, color = 'red'))%>% 
  add_text(x = 15, y = -0.05, z = 0,
           text = '     15  <br>Guideline<br>Recomm\'s', textposition = 'middle center',
           textfont = list(size = 12, color = 'red')) %>%
  layout(showlegend = FALSE)

fig4_7

saveNetwork(fig4_7, "plot.html")
webshot2::webshot("plot.html", 
                  "IPH_moderate.png",
                  vwidth = 1200, vheight = 800, zoom = 3)


#Crop the PNG so it can be stitched directly later
img_imager <- load.image("IPH_moderate.png")  # Load image
plot(img_imager)
dim(img_imager)
cropped_img <- imsub(img_imager, x %in% 1050:2700, y %in% 720:1950)
plot(cropped_img)
imager::save.image(cropped_img, "IPH_moderate.png")




#####################################################Plotting: Vigorous PA------------------------------
######################Tertiles -------------------
cutvalue <- incident_followup %>%
  group_by(Vigorous_3) %>%
  summarize(quant0 = quantile(METh_vigorous, probs = 0),
            quant25 = quantile(METh_vigorous, probs = 0.25),
            quant50 = quantile(METh_vigorous, probs = 0.50),
            quant75 = quantile(METh_vigorous, probs = 0.75),
            quant100 = quantile(METh_vigorous, probs = 1.00))
cutvalue
cutvalue_low1 <- cutvalue$quant0[1]
cutvalue_low2 <- cutvalue$quant0[2]
cutvalue_low3 <- cutvalue$quant0[3]
cutvalue_high1 <- cutvalue$quant100[1]
cutvalue_high2 <- cutvalue$quant100[2]
cutvalue_high3 <- 50
cutvalue_1 <- (cutvalue_low2 + cutvalue_high1)/2
cutvalue_2 <- (cutvalue_low3 + cutvalue_high2)/2
cutvalue_3 <- cutvalue_high3
cutvalue_1
cutvalue_2
cutvalue_3
X_position_1 <- (0 + cutvalue_1) / 2
X_position_2 <- (cutvalue_1 + cutvalue_2) / 2
X_position_3 <- (cutvalue_2 + cutvalue_3) / 2
X_position_1
X_position_2
X_position_3
X_position <- c(X_position_1, X_position_2, X_position_3)

rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + Vigorous_3
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
fit_rg0 <- summary(pool(rg0))
fit_rg0$low <- fit_rg0$estimate - 1.96 * fit_rg0$std.error
fit_rg0$high <- fit_rg0$estimate + 1.96 * fit_rg0$std.error
fit_rg1 <- subset(fit_rg0, select = c("term", "estimate", "low", "high"))[14:15,]
fit_rg1$estimate <- round(exp(fit_rg1$estimate), digits = 2)
fit_rg1$low <- round(exp(fit_rg1$low), digits = 2)
fit_rg1$high <- round(exp(fit_rg1$high), digits = 2)

contrast_1 <- data.frame(term = "contrast", estimate = 1, low = 1, high = 1)
fit_rg1 <- rbind(contrast_1, fit_rg1)
fit_rg1
fit_rg2 <- cbind(fit_rg1, X_position)
dput(names(fit_rg2))
#Extract OR, CI, and p-value
fit_final <- subset(fit_rg0, select = c("term", "estimate", "low", "high", "p.value"))[14:15,]
fit_final$estimate <- round(exp(fit_final$estimate), digits = 2)
fit_final$low <- round(exp(fit_final$low), digits = 2)
fit_final$high <- round(exp(fit_final$high), digits = 2)
fit_final$p.value <- round(fit_final$p.value, digits = 3)
fit_final_IPH <- rbind(fit_final_IPH, fit_final)

####Define the X-axis position
Y_location <- 1

####Define colors
color_point <- "rgb(28, 128, 65)"
color <- "rgba(28, 128, 65, 1)"
color_trans <- "rgba(28, 128, 65, 0.8)"
color_vertex <- c(28, 128, 65)

####Plot points
fit_rg2
fit_rg2$neg_error <- fit_rg2$estimate - fit_rg2$low
fit_rg2$pos_error <- fit_rg2$high - fit_rg2$estimate
fig1_0 <- plot_ly() %>%
  add_markers(data = fit_rg2, 
              x = fit_rg2$X_position, 
              y = Y_location, 
              z = fit_rg2$estimate,
              type = "scatter3d", 
              mode = "markers",
              marker = list(size = 1.5, color = color_point, opacity = 1),
              showlegend = FALSE,
              error_z = list(type = "data", 
                             array = ~pos_error, 
                             arrayminus = ~neg_error, 
                             color = color_point,
                             thickness = 5,#Since this is a 3D plot, thickness has no effect here; consider adding lines manually for better control.
                             width = 5)) %>% 
  layout(
    scene = list(
      xaxis = list(title = "Vigorous physical activity<br>(unit: METh-week)", 
                   range = c(0, cutvalue_high3), 
                   fixedrange = TRUE, 
                   showticklabels = FALSE, 
                   showgrid = FALSE, 
                   zeroline = FALSE, 
                   showline = FALSE),
      yaxis = list(title = "", 
                   range = c(-0.5, 11.5), 
                   fixedrange = TRUE, 
                   showticklabels = FALSE, 
                   showgrid = FALSE, 
                   zeroline = FALSE, 
                   showline = FALSE), 
      zaxis = list(title = "", 
                   range = c(0, 6), 
                   fixedrange = TRUE, 
                   showticklabels = FALSE, 
                   showgrid = FALSE, 
                   zeroline = FALSE, 
                   showline = FALSE))
  ) %>%
  add_trace(
    x = c(0), y = c(0), z = c(2.2),  # Adjust Z-axis title position
    type = "scatter3d",
    mode = "text",
    text = "Odds ratio of <br> incident IPH",
    textfont = list(size = 12),
    showlegend = FALSE
  )
fig1_0

####Draw the corresponding X-axis
X_axis <- data.frame(x = c(0, cutvalue_high3),
                     y = c(Y_location, Y_location),
                     z = c(0, 0))
#Add an arrow at the end of the X-axis
arrow_x <- cutvalue_high3
arrow_y <- Y_location
arrow_z <- 0
fig2_1 <- fig1_0 %>% 
  add_trace(
    x = X_axis$x,
    y = X_axis$y,
    z = X_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3, dash = "dash"),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(1),  #Arrow points along the X-axis
    v = list(0),
    w = list(0),
    sizemode = "absolute",
    sizeref = 3,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add quantile labels on the X-axis
fig2_1 <- fig2_1 %>%
  add_trace(
    x = cutvalue_1,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_1, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = cutvalue_2,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_2, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig2_1

####Draw the corresponding Z-axis
Z_axis <- data.frame(x = c(0, 0),
                     y = c(Y_location, Y_location),
                     z = c(0, 2))
# Add an arrow at the end of the Z-axis
arrow_x <- 0
arrow_y <- Y_location
arrow_z <- 2
fig2_2 <- fig2_1 %>% 
  add_trace(
    x = Z_axis$x,
    y = Z_axis$y,
    z = Z_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(0),  
    v = list(0),
    w = list(1),#Arrow points along the Z-axis
    sizemode = "absolute",
    sizeref = 0.35,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add labels on the Z-axis
fig2_3 <- fig2_2 %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = 0,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 1,
    type = "scatter3d",
    mode = "text",
    text = 1,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 2,
    type = "scatter3d",
    mode = "text",
    text = 2,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig2_3


####Add a visual aid plane for the reference line
square_x <- c(0, 0, cutvalue_high3, cutvalue_high3)
square_z <- c(0, 1, 1, 0)
square_y1 <- c(Y_location, Y_location, Y_location, Y_location)
i_set <- c(0,2)
j_set <- c(1,3)
k_set <- c(2,0)
vertex_colors <- matrix(rep(color_vertex, each = 4), nrow = 4)/255
fig2_4 <- fig2_3  %>%
  add_trace(
    x = square_x, 
    z = square_z, 
    y = square_y1, 
    i = i_set, 
    j = j_set, 
    k = k_set,
    type = "mesh3d",
    vertexcolor = vertex_colors,
    opacity = 0.2,
    showlegend = FALSE)
fig2_4

####Add auxiliary lines at the quantiles
line_x1 <- c(cutvalue_1, cutvalue_1)
line_z1 <- c(0, 1)
line_y1 <- c(Y_location, Y_location)
line_x2 <- c(cutvalue_2, cutvalue_2)
line_z2 <- c(0, 1)
line_y2 <- c(Y_location, Y_location)
fig2_5 <- fig2_4  %>% 
  add_trace(
    x = line_x1,
    y = line_y1,
    z = line_z1,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE) %>% 
  add_trace(
    x = line_x2,
    y = line_y2,
    z = line_z2,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE)
fig2_5

####Label the grouping method at the end of the X-axis
fig2_6 <- fig2_5 %>% 
  add_trace(
    x = cutvalue_high3,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = "Tertile",
    textposition = "right",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig2_6




######################Quartiles -------------------
#total_quartiles
cutvalue <- incident_followup %>%
  group_by(Vigorous_4) %>%
  summarize(quant0 = quantile(METh_vigorous, probs = 0),
            quant25 = quantile(METh_vigorous, probs = 0.25),
            quant50 = quantile(METh_vigorous, probs = 0.50),
            quant75 = quantile(METh_vigorous, probs = 0.75),
            quant100 = quantile(METh_vigorous, probs = 1.00))
cutvalue
cutvalue_low1 <- cutvalue$quant0[1]
cutvalue_low2 <- cutvalue$quant0[2]
cutvalue_low3 <- cutvalue$quant0[3]
cutvalue_low4 <- cutvalue$quant0[4]
cutvalue_high1 <- cutvalue$quant100[1]
cutvalue_high2 <- cutvalue$quant100[2]
cutvalue_high3 <- cutvalue$quant100[3]
cutvalue_high4 <- 50
cutvalue_1 <- (cutvalue_low2 + cutvalue_high1)/2
cutvalue_2 <- (cutvalue_low3 + cutvalue_high2)/2
cutvalue_3 <- (cutvalue_low4 + cutvalue_high3)/2
cutvalue_4 <- cutvalue_high4
cutvalue_1
cutvalue_2
cutvalue_3
cutvalue_4
X_position_1 <- (0 + cutvalue_1) / 2
X_position_2 <- (cutvalue_1 + cutvalue_2) / 2
X_position_3 <- (cutvalue_2 + cutvalue_3) / 2
X_position_4 <- (cutvalue_3 + cutvalue_4) / 2
X_position_1
X_position_2
X_position_3
X_position_4
X_position <- c(X_position_1, X_position_2, X_position_3, X_position_4)

rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + Vigorous_4
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
fit_rg0 <- summary(pool(rg0))
fit_rg0$low <- fit_rg0$estimate - 1.96 * fit_rg0$std.error
fit_rg0$high <- fit_rg0$estimate + 1.96 * fit_rg0$std.error
fit_rg1 <- subset(fit_rg0, select = c("term", "estimate", "low", "high"))[14:16,]
fit_rg1$estimate <- round(exp(fit_rg1$estimate), digits = 2)
fit_rg1$low <- round(exp(fit_rg1$low), digits = 2)
fit_rg1$high <- round(exp(fit_rg1$high), digits = 2)

contrast_1 <- data.frame(term = "contrast", estimate = 1, low = 1, high = 1)
fit_rg1 <- rbind(contrast_1, fit_rg1)
fit_rg1
fit_rg2 <- cbind(fit_rg1, X_position)
dput(names(fit_rg2))
#Extract OR, CI, and p-value
fit_final <- subset(fit_rg0, select = c("term", "estimate", "low", "high", "p.value"))[14:16,]
fit_final$estimate <- round(exp(fit_final$estimate), digits = 2)
fit_final$low <- round(exp(fit_final$low), digits = 2)
fit_final$high <- round(exp(fit_final$high), digits = 2)
fit_final$p.value <- round(fit_final$p.value, digits = 3)
fit_final_IPH <- rbind(fit_final_IPH, fit_final)

####Define the X-axis position
Y_location <- 5

####Define colors
color_point <- "rgb(229, 87, 9)"
color <- "rgba(229, 87, 9, 1)"
color_trans <- "rgba(229, 87, 9, 0.8)"
color_vertex <- c(229, 87, 9)

####Plot points
fit_rg2
fit_rg2$neg_error <- fit_rg2$estimate - fit_rg2$low
fit_rg2$pos_error <- fit_rg2$high - fit_rg2$estimate
fig3_0 <- fig2_6 %>%
  add_markers(data = fit_rg2, 
              x = fit_rg2$X_position, 
              y = Y_location, 
              z = fit_rg2$estimate,
              type = "scatter3d", 
              mode = "markers",
              marker = list(size = 1.5, color = color_point, opacity = 1),
              showlegend = FALSE,
              error_z = list(type = "data", 
                             array = ~pos_error, 
                             arrayminus = ~neg_error, 
                             color = color_point,
                             thickness = 5,#Since this is a 3D plot, thickness has no effect here; consider adding lines manually for better control.
                             width = 5))
fig3_0

####Draw the corresponding X-axis
X_axis <- data.frame(x = c(0, cutvalue_high4),
                     y = c(Y_location, Y_location),
                     z = c(0, 0))
#Add an arrow at the end of the X-axis
arrow_x <- cutvalue_high4
arrow_y <- Y_location
arrow_z <- 0
fig3_1 <- fig3_0 %>% 
  add_trace(
    x = X_axis$x,
    y = X_axis$y,
    z = X_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3, dash = "dash"),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(1),  #Arrow points along the X-axis
    v = list(0),
    w = list(0),
    sizemode = "absolute",
    sizeref = 3,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add quantile labels on the X-axis
fig3_1 <- fig3_1 %>%
  add_trace(
    x = cutvalue_1,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_1, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = cutvalue_2,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_2, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = cutvalue_3,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_3, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig3_1

####Draw the corresponding Z-axis
Z_axis <- data.frame(x = c(0, 0),
                     y = c(Y_location, Y_location),
                     z = c(0, 2))
# Add an arrow at the end of the Z-axis
arrow_x <- 0
arrow_y <- Y_location
arrow_z <- 2
fig3_2 <- fig3_1 %>% 
  add_trace(
    x = Z_axis$x,
    y = Z_axis$y,
    z = Z_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(0),  
    v = list(0),
    w = list(1),#Arrow points along the Z-axis
    sizemode = "absolute",
    sizeref = 0.35,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add labels on the Z-axis
fig3_3 <- fig3_2 %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = 0,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 1,
    type = "scatter3d",
    mode = "text",
    text = 1,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 2,
    type = "scatter3d",
    mode = "text",
    text = 2,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig3_3


####Add a visual aid plane for the reference line
square_x <- c(0, 0, cutvalue_high4, cutvalue_high4)
square_z <- c(0, 1, 1, 0)
square_y1 <- c(Y_location, Y_location, Y_location, Y_location)
i_set <- c(0,2)
j_set <- c(1,3)
k_set <- c(2,0)
vertex_colors <- matrix(rep(color_vertex, each = 4), nrow = 4)/255
fig3_4 <- fig3_3  %>%
  add_trace(
    x = square_x, 
    z = square_z, 
    y = square_y1, 
    i = i_set, 
    j = j_set, 
    k = k_set,
    type = "mesh3d",
    vertexcolor = vertex_colors,
    opacity = 0.2,
    showlegend = FALSE)
fig3_4

####Add auxiliary lines at the quantiles
line_x1 <- c(cutvalue_1, cutvalue_1)
line_z1 <- c(0, 1)
line_y1 <- c(Y_location, Y_location)
line_x2 <- c(cutvalue_2, cutvalue_2)
line_z2 <- c(0, 1)
line_y2 <- c(Y_location, Y_location)
line_x3 <- c(cutvalue_3, cutvalue_3)
line_z3 <- c(0, 1)
line_y3 <- c(Y_location, Y_location)
fig3_5 <- fig3_4  %>% 
  add_trace(
    x = line_x1,
    y = line_y1,
    z = line_z1,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE) %>% 
  add_trace(
    x = line_x2,
    y = line_y2,
    z = line_z2,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE) %>% 
  add_trace(
    x = line_x3,
    y = line_y3,
    z = line_z3,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE)
fig3_5

####Label the grouping method at the end of the X-axis
fig3_6 <- fig3_5 %>% 
  add_trace(
    x = cutvalue_high4,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = "Quartile",
    textposition = "right",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig3_6





######################Quintiles -------------------
#total_quintiles
cutvalue <- incident_followup %>%
  group_by(Vigorous_5) %>%
  summarize(quant0 = quantile(METh_vigorous, probs = 0),
            quant25 = quantile(METh_vigorous, probs = 0.25),
            quant50 = quantile(METh_vigorous, probs = 0.50),
            quant75 = quantile(METh_vigorous, probs = 0.75),
            quant100 = quantile(METh_vigorous, probs = 1.00))
cutvalue
cutvalue_low1 <- cutvalue$quant0[1]
cutvalue_low2 <- cutvalue$quant0[2]
cutvalue_low3 <- cutvalue$quant0[3]
cutvalue_low4 <- cutvalue$quant0[4]
cutvalue_low5 <- cutvalue$quant0[5]
cutvalue_high1 <- cutvalue$quant100[1]
cutvalue_high2 <- cutvalue$quant100[2]
cutvalue_high3 <- cutvalue$quant100[3]
cutvalue_high4 <- cutvalue$quant100[4]
cutvalue_high5 <- 50
cutvalue_1 <- (cutvalue_low2 + cutvalue_high1)/2
cutvalue_2 <- (cutvalue_low3 + cutvalue_high2)/2
cutvalue_3 <- (cutvalue_low4 + cutvalue_high3)/2
cutvalue_4 <- (cutvalue_low5 + cutvalue_high4)/2
cutvalue_5 <- cutvalue_high5
cutvalue_1
cutvalue_2
cutvalue_3
cutvalue_4
cutvalue_5
X_position_1 <- (0 + cutvalue_1) / 2
X_position_2 <- (cutvalue_1 + cutvalue_2) / 2
X_position_3 <- (cutvalue_2 + cutvalue_3) / 2
X_position_4 <- (cutvalue_3 + cutvalue_4) / 2
X_position_5 <- (cutvalue_5 + cutvalue_4) / 2
X_position_1
X_position_2
X_position_3
X_position_4
X_position_5

X_position <- c(X_position_1, X_position_2, X_position_3, X_position_4, X_position_5)

rg0 <- with(incident_followup_1_long_3,
            geeglm(IPH
                   ~ Sex + First_age + followuptime + BMI + smokingstatus_2 
                   + rs_cohort + edu_2 + hyperchol+ hypertension + dm + Maximum
                   + Vigorous_5
                   , family = binomial, id = ergoid, corstr = "independence",  weights = IPW_stabilized))
summary(pool(rg0))
fit_rg0 <- summary(pool(rg0))
fit_rg0$low <- fit_rg0$estimate - 1.96 * fit_rg0$std.error
fit_rg0$high <- fit_rg0$estimate + 1.96 * fit_rg0$std.error
fit_rg1 <- subset(fit_rg0, select = c("term", "estimate", "low", "high"))[14:17,]
fit_rg1$estimate <- round(exp(fit_rg1$estimate), digits = 2)
fit_rg1$low <- round(exp(fit_rg1$low), digits = 2)
fit_rg1$high <- round(exp(fit_rg1$high), digits = 2)

contrast_1 <- data.frame(term = "contrast", estimate = 1, low = 1, high = 1)
fit_rg1 <- rbind(contrast_1, fit_rg1)
fit_rg1
fit_rg2 <- cbind(fit_rg1, X_position)
dput(names(fit_rg2))
fit_rg2
#Extract OR, CI, and p-value
fit_final <- subset(fit_rg0, select = c("term", "estimate", "low", "high", "p.value"))[14:17,]
fit_final$estimate <- round(exp(fit_final$estimate), digits = 2)
fit_final$low <- round(exp(fit_final$low), digits = 2)
fit_final$high <- round(exp(fit_final$high), digits = 2)
fit_final$p.value <- round(fit_final$p.value, digits = 3)
fit_final_IPH <- rbind(fit_final_IPH, fit_final)
####Define the X-axis position
Y_location <- 10

####Define colors
color_point <- "rgb(43, 106, 153)"
color <- "rgba(43, 106, 153, 1)"
color_trans <- "rgba(43, 106, 153, 0.8)"
color_vertex <- c(43, 106, 153)

####Plot points
fit_rg2
fit_rg2$neg_error <- fit_rg2$estimate - fit_rg2$low
fit_rg2$pos_error <- fit_rg2$high - fit_rg2$estimate
fig4_0 <- fig3_6 %>%
  add_markers(data = fit_rg2, 
              x = fit_rg2$X_position, 
              y = Y_location, 
              z = fit_rg2$estimate,
              type = "scatter3d", 
              mode = "markers",
              marker = list(size = 1.5, color = color_point, opacity = 1),
              showlegend = FALSE,
              error_z = list(type = "data", 
                             array = ~pos_error, 
                             arrayminus = ~neg_error, 
                             color = color_point,
                             thickness = 5,#Since this is a 3D plot, thickness has no effect here; consider adding lines manually for better control.
                             width = 5))
fig4_0

####Draw the corresponding X-axis
X_axis <- data.frame(x = c(0, cutvalue_high5),
                     y = c(Y_location, Y_location),
                     z = c(0, 0))
#Add an arrow at the end of the X-axis
arrow_x <- cutvalue_high5
arrow_y <- Y_location
arrow_z <- 0
fig4_1 <- fig4_0 %>% 
  add_trace(
    x = X_axis$x,
    y = X_axis$y,
    z = X_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3, dash = "dash"),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(1),  #Arrow points along the X-axis
    v = list(0),
    w = list(0),
    sizemode = "absolute",
    sizeref = 3,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add quantile labels on the X-axis
fig4_1 <- fig4_1 %>%
  add_trace(
    x = cutvalue_1,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_1, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = cutvalue_2,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_2, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = cutvalue_3,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_3, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = cutvalue_4,
    y = Y_location,
    z = 0.2,
    type = "scatter3d",
    mode = "text",
    text = round(cutvalue_4, 0),
    textposition = "bottom",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig4_1

####Draw the corresponding Z-axis
Z_axis <- data.frame(x = c(0, 0),
                     y = c(Y_location, Y_location),
                     z = c(0, 2))
# Add an arrow at the end of the Z-axis
arrow_x <- 0
arrow_y <- Y_location
arrow_z <- 2
fig4_2 <- fig4_1 %>% 
  add_trace(
    x = Z_axis$x,
    y = Z_axis$y,
    z = Z_axis$z,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color, width = 3),
    showlegend = FALSE) %>%
  add_trace(
    type = "cone",
    x = list(arrow_x),
    y = list(arrow_y),
    z = list(arrow_z),
    u = list(0),  
    v = list(0),
    w = list(1),#Arrow points along the Z-axis
    sizemode = "absolute",
    sizeref = 0.35,  #Adjust arrow size
    colorscale = list(c(0, color), c(1, color)),  #Color and whether solid-filled
    showscale = FALSE  #Do not display the color legend
  )
# Add labels on the Z-axis
fig4_3 <- fig4_2 %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = 0,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 1,
    type = "scatter3d",
    mode = "text",
    text = 1,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE) %>%
  add_trace(
    x = 0,
    y = Y_location,
    z = 2,
    type = "scatter3d",
    mode = "text",
    text = 2,
    textposition = "left",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig4_3


####Add a visual aid plane for the reference line
square_x <- c(0, 0, cutvalue_high5, cutvalue_high5)
square_z <- c(0, 1, 1, 0)
square_y1 <- c(Y_location, Y_location, Y_location, Y_location)
i_set <- c(0,2)
j_set <- c(1,3)
k_set <- c(2,0)
vertex_colors <- matrix(rep(color_vertex, each = 4), nrow = 4)/255
fig4_4 <- fig4_3  %>%
  add_trace(
    x = square_x, 
    z = square_z, 
    y = square_y1, 
    i = i_set, 
    j = j_set, 
    k = k_set,
    type = "mesh3d",
    vertexcolor = vertex_colors,
    opacity = 0.2,
    showlegend = FALSE)
fig4_4

####Add auxiliary lines at the quantiles
line_x1 <- c(cutvalue_1, cutvalue_1)
line_z1 <- c(0, 1)
line_y1 <- c(Y_location, Y_location)
line_x2 <- c(cutvalue_2, cutvalue_2)
line_z2 <- c(0, 1)
line_y2 <- c(Y_location, Y_location)
line_x3 <- c(cutvalue_3, cutvalue_3)
line_z3 <- c(0, 1)
line_y3 <- c(Y_location, Y_location)
line_x4 <- c(cutvalue_4, cutvalue_4)
line_z4 <- c(0, 1)
line_y4 <- c(Y_location, Y_location)
fig4_5 <- fig4_4  %>% 
  add_trace(
    x = line_x1,
    y = line_y1,
    z = line_z1,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE) %>% 
  add_trace(
    x = line_x2,
    y = line_y2,
    z = line_z2,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE) %>% 
  add_trace(
    x = line_x3,
    y = line_y3,
    z = line_z3,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE) %>% 
  add_trace(
    x = line_x4,
    y = line_y4,
    z = line_z4,
    type = "scatter3d",
    mode = "lines",
    line = list(color = color_trans, width = 3, dash = "longdash"),
    showlegend = FALSE)
fig4_5

####Label the grouping method at the end of the X-axis
fig4_6 <- fig4_5 %>% 
  add_trace(
    x = cutvalue_high5,
    y = Y_location,
    z = 0,
    type = "scatter3d",
    mode = "text",
    text = "Quintile",
    textposition = "right",
    textfont = list(color = color, size = 11),
    showlegend = FALSE)
fig4_6
fig4_7 <- fig4_6 %>% layout(
  scene = list(camera = list(eye = list(x = 0.59097400646812713, y = -1.878239943166165, z = 1.309212875732691),
                             aspectmode = "manual",
                             aspectratio = list(x = 1, y = 2, z = 1)  # Y-axis is stretched
  ))
)
###Add the PA guideline-recommended threshold
fig4_7 <- fig4_7 %>% 
  add_text(x = 7.5, y = 1, z = 0,
           text = '*', textposition = 'middle center',
           textfont = list(size = 20, color = 'red'))%>% 
  add_text(x = 7.5, y = 5, z = 0,
           text = '*', textposition = 'middle center',
           textfont = list(size = 20, color = 'red'))%>% 
  add_text(x = 7.5, y = 10, z = 0,
           text = '*', textposition = 'middle center',
           textfont = list(size = 20, color = 'red'))%>% 
  add_text(x = 7.5, y = -0.09, z = 0,
           text = '     7.5  <br>Guideline<br>Recomm\'s', textposition = 'middle center',
           textfont = list(size = 12, color = 'red')) %>%
  layout(showlegend = FALSE)

fig4_7

saveNetwork(fig4_7, "plot.html")
webshot2::webshot("plot.html", 
                  "IPH_vigorous.png",
                  vwidth = 1200, vheight = 800, zoom = 3)


#Crop the PNG so it can be stitched directly later
img_imager <- load.image("IPH_vigorous.png")  # Load image
plot(img_imager)
dim(img_imager)
cropped_img <- imsub(img_imager, x %in% 1050:2700, y %in% 720:1950)
plot(cropped_img)
imager::save.image(cropped_img, "IPH_vigorous.png")
