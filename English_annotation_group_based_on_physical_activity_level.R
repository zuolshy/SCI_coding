incident_followup <- left_join(data_total_twoscan, data_PA, by = "ergoid")
Stroke <- read_sav('Stroke_data_Luo_Oct2024.sav')
names(Stroke)
Stroke_0 <- subset(Stroke, select = c("ergoid","IC", "prevalent_stroke","Baseline_date","incident_stroke","stroke_date"))
Stroke_1 <- rename(Stroke_0, ic_Stroke = IC, prev_stroke_date = Baseline_date, prev_Stroke = prevalent_stroke, inc_stroke = incident_stroke, inc_stroke_date = stroke_date)
Stroke_1 <- subset(Stroke_1, ic_Stroke == 1)
incident_followup <- left_join(incident_followup, Stroke_1[c("ergoid", "prev_Stroke")], by = "ergoid")
table(incident_followup$prev_Stroke)

incident_followup <- subset(incident_followup, !is.na(prev_Stroke))
incident_followup <- subset(incident_followup, prev_Stroke == 0)
incident_followup <- subset(incident_followup, !is.na(total_MET))
data_followup_baselineimt <- subset(incident_followup, scan == 1)

data_followup_baselineimt <- data_followup_baselineimt[, c("ergoid", "Maximum", "stenosis","side")]
incident_followup <- left_join(incident_followup, income_edu, by="ergoid")
dput(names(incident_followup))
incident_followup <- incident_followup[,c(c("Cal", "IPH", "Lipid", "plaque", "ergoid", 
                                            "side_name", "side", "ergocar_date", "ergocarfup_date", "birth_date", 
                                            "First_age", "Follow", "Sex", "scan", "status", "status_name", 
                                            "followuptime", "total_MET", "METh_modvig", "METh_vigorous", 
                                            "METh_moderate", "prev_Stroke", "income", "income_cate", "rs_cohort", 
                                            "edu"))]
incident_followup <- left_join(incident_followup, data_followup_baselineimt, by = c("ergoid", "side"))
incident_followup <- left_join(incident_followup, energy, by="ergoid")
length(unique(incident_followup$ergoid))

Scan_covar_use <- subset(Scan_covar, select = -c(Sex,First_age))
incident_followup <- left_join(incident_followup, Scan_covar_use, by = "ergoid")
summary(incident_followup)

####group by level of physical activity----
######based on literature
incident_followup <- mutate(incident_followup, total_guideline = case_when(
  total_MET < 25 ~ 0,
  total_MET >= 25 & total_MET < 50 ~ 1,
  total_MET >= 50 ~ 2
))
incident_followup$total_guideline <- factor(incident_followup$total_guideline)

incident_followup <- mutate(incident_followup, METh_modvig_guideline = case_when(
  METh_modvig < 25 ~ 0,
  METh_modvig >= 25 & METh_modvig < 50 ~ 1,
  METh_modvig >= 50 ~ 2
))
incident_followup$METh_modvig_guideline <- factor(incident_followup$METh_modvig_guideline)

incident_followup$Vigorous <- incident_followup$METh_vigorous
incident_followup$Vigorous[incident_followup$Vigorous == 0] <- NA

######tertile
cutoff <- incident_followup %>%
  group_by(Sex) %>%
  summarise(cut_33.33 = quantile(total_MET, 0.33333),
            cut_66.66 = quantile(total_MET, 0.66666))
cutoff_men_1 <- cutoff$cut_33.33[1]
cutoff_women_1 <- cutoff$cut_33.33[2]
cutoff_men_2 <- cutoff$cut_66.66[1]
cutoff_women_2 <- cutoff$cut_66.66[2]
incident_followup <- mutate(incident_followup, total_MET_3 = case_when(
  Sex == 0 & total_MET<= cutoff_men_1 ~ 1,
  Sex == 0 & total_MET> cutoff_men_1 & total_MET <= cutoff_men_2 ~ 2,
  Sex == 0 & total_MET> cutoff_men_2 ~ 3,
  Sex == 1 & total_MET<= cutoff_women_1 ~ 1,
  Sex == 1 & total_MET> cutoff_women_1 & total_MET <= cutoff_women_2 ~ 2,
  Sex == 1 & total_MET> cutoff_women_2 ~ 3
))
incident_followup$total_MET_3 <- factor(incident_followup$total_MET_3)
table(incident_followup$total_MET_3)

cutoff <- incident_followup %>%
  group_by(Sex) %>%
  summarise(cut_33.33 = quantile(METh_modvig, 0.33333),
            cut_66.66 = quantile(METh_modvig, 0.66666))
cutoff_men_1 <- cutoff$cut_33.33[1]
cutoff_women_1 <- cutoff$cut_33.33[2]
cutoff_men_2 <- cutoff$cut_66.66[1]
cutoff_women_2 <- cutoff$cut_66.66[2]
incident_followup <- mutate(incident_followup, METh_modvig_3 = case_when(
  Sex == 0 & METh_modvig<= cutoff_men_1 ~ 1,
  Sex == 0 & METh_modvig> cutoff_men_1 & METh_modvig <= cutoff_men_2 ~ 2,
  Sex == 0 & METh_modvig> cutoff_men_2 ~ 3,
  Sex == 1 & METh_modvig<= cutoff_women_1 ~ 1,
  Sex == 1 & METh_modvig> cutoff_women_1 & METh_modvig <= cutoff_women_2 ~ 2,
  Sex == 1 & METh_modvig> cutoff_women_2 ~ 3
))
incident_followup$METh_modvig_3 <- factor(incident_followup$METh_modvig_3)
table(incident_followup$METh_modvig_3)

cutoff <- incident_followup %>%
  group_by(Sex) %>%
  summarise(cut_33.33 = quantile(METh_moderate, 0.33333),
            cut_66.66 = quantile(METh_moderate, 0.66666))
cutoff_men_1 <- cutoff$cut_33.33[1]
cutoff_women_1 <- cutoff$cut_33.33[2]
cutoff_men_2 <- cutoff$cut_66.66[1]
cutoff_women_2 <- cutoff$cut_66.66[2]
incident_followup <- mutate(incident_followup, METh_moderate_3 = case_when(
  Sex == 0 & METh_moderate<= cutoff_men_1 ~ 1,
  Sex == 0 & METh_moderate> cutoff_men_1 & METh_moderate <= cutoff_men_2 ~ 2,
  Sex == 0 & METh_moderate> cutoff_men_2 ~ 3,
  Sex == 1 & METh_moderate<= cutoff_women_1 ~ 1,
  Sex == 1 & METh_moderate> cutoff_women_1 & METh_moderate <= cutoff_women_2 ~ 2,
  Sex == 1 & METh_moderate> cutoff_women_2 ~ 3
))
incident_followup$METh_moderate_3 <- factor(incident_followup$METh_moderate_3)
table(incident_followup$METh_moderate_3)

cutoff <- incident_followup %>%
  group_by(Sex) %>%
  summarise(cut_50 = quantile(Vigorous, 0.5, na.rm = TRUE))
cutoff_men_1 <- cutoff$cut_50[1]
cutoff_women_1 <- cutoff$cut_50[2]
incident_followup <- mutate(incident_followup, Vigorous_3 = case_when(
  Sex == 0 & is.na(Vigorous) ~ 1,
  Sex == 0 & !is.na(Vigorous) & Vigorous < cutoff_men_1 ~ 2,
  Sex == 0 & !is.na(Vigorous) & Vigorous >= cutoff_men_1 ~ 3,
  Sex == 1 & is.na(Vigorous) ~ 1,
  Sex == 1 & !is.na(Vigorous) & Vigorous < cutoff_women_1 ~ 2,
  Sex == 1 & !is.na(Vigorous) & Vigorous >= cutoff_women_1 ~ 3
))
incident_followup$Vigorous_3 <- factor(incident_followup$Vigorous_3)
table(incident_followup$Vigorous_3)

######quartile
cutoff <- incident_followup %>%
  group_by(Sex) %>%
  summarise(cut_25 = quantile(total_MET, 0.25),
            cut_50 = quantile(total_MET, 0.50),
            cut_75 = quantile(total_MET, 0.75))
cutoff_men_1 <- cutoff$cut_25[1]
cutoff_women_1 <- cutoff$cut_25[2]
cutoff_men_2 <- cutoff$cut_50[1]
cutoff_women_2 <- cutoff$cut_50[2]
cutoff_men_3 <- cutoff$cut_75[1]
cutoff_women_3 <- cutoff$cut_75[2]
incident_followup <- mutate(incident_followup, total_MET_4 = case_when(
  Sex == 0 & total_MET<= cutoff_men_1 ~ 1,
  Sex == 0 & total_MET> cutoff_men_1 & total_MET <= cutoff_men_2 ~ 2,
  Sex == 0 & total_MET> cutoff_men_2 & total_MET <= cutoff_men_3 ~ 3,
  Sex == 0 & total_MET> cutoff_men_3 ~ 4,
  Sex == 1 & total_MET<= cutoff_women_1 ~ 1,
  Sex == 1 & total_MET> cutoff_women_1 & total_MET <= cutoff_women_2 ~ 2,
  Sex == 1 & total_MET> cutoff_women_2 & total_MET <= cutoff_women_3 ~ 3,
  Sex == 1 & total_MET> cutoff_women_3 ~ 4
))
incident_followup$total_MET_4 <- factor(incident_followup$total_MET_4)
table(incident_followup$total_MET_4)

cutoff <- incident_followup %>%
  group_by(Sex) %>%
  summarise(cut_25 = quantile(METh_modvig, 0.25),
            cut_50 = quantile(METh_modvig, 0.50),
            cut_75 = quantile(METh_modvig, 0.75))
cutoff_men_1 <- cutoff$cut_25[1]
cutoff_women_1 <- cutoff$cut_25[2]
cutoff_men_2 <- cutoff$cut_50[1]
cutoff_women_2 <- cutoff$cut_50[2]
cutoff_men_3 <- cutoff$cut_75[1]
cutoff_women_3 <- cutoff$cut_75[2]
incident_followup <- mutate(incident_followup, METh_modvig_4 = case_when(
  Sex == 0 & METh_modvig<= cutoff_men_1 ~ 1,
  Sex == 0 & METh_modvig> cutoff_men_1 & METh_modvig <= cutoff_men_2 ~ 2,
  Sex == 0 & METh_modvig> cutoff_men_2 & METh_modvig <= cutoff_men_3 ~ 3,
  Sex == 0 & METh_modvig> cutoff_men_3 ~ 4,
  Sex == 1 & METh_modvig<= cutoff_women_1 ~ 1,
  Sex == 1 & METh_modvig> cutoff_women_1 & METh_modvig <= cutoff_women_2 ~ 2,
  Sex == 1 & METh_modvig> cutoff_women_2 & METh_modvig <= cutoff_women_3 ~ 3,
  Sex == 1 & METh_modvig> cutoff_women_3 ~ 4
))
incident_followup$METh_modvig_4 <- factor(incident_followup$METh_modvig_4)
table(incident_followup$METh_modvig_4)

cutoff <- incident_followup %>%
  group_by(Sex) %>%
  summarise(cut_25 = quantile(METh_moderate, 0.25),
            cut_50 = quantile(METh_moderate, 0.50),
            cut_75 = quantile(METh_moderate, 0.75))
cutoff_men_1 <- cutoff$cut_25[1]
cutoff_women_1 <- cutoff$cut_25[2]
cutoff_men_2 <- cutoff$cut_50[1]
cutoff_women_2 <- cutoff$cut_50[2]
cutoff_men_3 <- cutoff$cut_75[1]
cutoff_women_3 <- cutoff$cut_75[2]
incident_followup <- mutate(incident_followup, METh_moderate_4 = case_when(
  Sex == 0 & METh_moderate<= cutoff_men_1 ~ 1,
  Sex == 0 & METh_moderate> cutoff_men_1 & METh_moderate <= cutoff_men_2 ~ 2,
  Sex == 0 & METh_moderate> cutoff_men_2 & METh_moderate <= cutoff_men_3 ~ 3,
  Sex == 0 & METh_moderate> cutoff_men_3 ~ 4,
  Sex == 1 & METh_moderate<= cutoff_women_1 ~ 1,
  Sex == 1 & METh_moderate> cutoff_women_1 & METh_moderate <= cutoff_women_2 ~ 2,
  Sex == 1 & METh_moderate> cutoff_women_2 & METh_moderate <= cutoff_women_3 ~ 3,
  Sex == 1 & METh_moderate> cutoff_women_3 ~ 4
))
incident_followup$METh_moderate_4 <- factor(incident_followup$METh_moderate_4)
table(incident_followup$METh_moderate_4)

cutoff <- incident_followup %>%
  group_by(Sex) %>%
  summarise(cut_33 = quantile(Vigorous, 0.33, na.rm = TRUE),
            cut_66 = quantile(Vigorous, 0.66, na.rm = TRUE))
cutoff_men_1 <- cutoff$cut_33[1]
cutoff_women_1 <- cutoff$cut_33[2]
cutoff_men_2 <- cutoff$cut_66[1]
cutoff_women_2 <- cutoff$cut_66[2]
incident_followup <- mutate(incident_followup, Vigorous_4 = case_when(
  Sex == 0 & is.na(Vigorous) ~ 1,
  Sex == 0 & !is.na(Vigorous) & Vigorous < cutoff_men_1 ~ 2,
  Sex == 0 & !is.na(Vigorous) & Vigorous >= cutoff_men_1 & Vigorous < cutoff_men_2 ~ 3,
  Sex == 0 & !is.na(Vigorous) & Vigorous >= cutoff_men_2 ~ 4,
  Sex == 1 & is.na(Vigorous) ~ 1,
  Sex == 1 & !is.na(Vigorous) & Vigorous < cutoff_women_1 ~ 2,
  Sex == 1 & !is.na(Vigorous) & Vigorous >= cutoff_women_1 & Vigorous < cutoff_women_2 ~ 3,
  Sex == 1 & !is.na(Vigorous) & Vigorous >= cutoff_women_2 ~ 4
))
incident_followup$Vigorous_4 <- factor(incident_followup$Vigorous_4)
table(incident_followup$Vigorous_4)


#####quintile
cutoff <- incident_followup %>%
  group_by(Sex) %>%
  summarise(cut_20 = quantile(total_MET, 0.20),
            cut_40 = quantile(total_MET, 0.40),
            cut_60 = quantile(total_MET, 0.60),
            cut_80 = quantile(total_MET, 0.82))
cutoff
cutoff_men_1 <- cutoff$cut_20[1]
cutoff_women_1 <- cutoff$cut_20[2]
cutoff_men_2 <- cutoff$cut_40[1]
cutoff_women_2 <- cutoff$cut_40[2]
cutoff_men_3 <- cutoff$cut_60[1]
cutoff_women_3 <- cutoff$cut_60[2]
cutoff_men_4 <- cutoff$cut_80[1]
cutoff_women_4 <- cutoff$cut_80[2]
incident_followup <- mutate(incident_followup, total_MET_5 = case_when(
  Sex == 0 & total_MET <= cutoff_men_1 ~ 1,
  Sex == 0 & total_MET > cutoff_men_1 & total_MET < cutoff_men_2 ~ 2,
  Sex == 0 & total_MET >= cutoff_men_2 & total_MET < cutoff_men_3 ~ 3,
  Sex == 0 & total_MET>= cutoff_men_3 & total_MET <= cutoff_men_4 ~ 4,
  Sex == 0 & total_MET> cutoff_men_4 ~ 5,
  Sex == 1 & total_MET<= cutoff_women_1 ~ 1,
  Sex == 1 & total_MET> cutoff_women_1 & total_MET <= cutoff_women_2 ~ 2,
  Sex == 1 & total_MET> cutoff_women_2 & total_MET < cutoff_women_3 ~ 3,
  Sex == 1 & total_MET>= cutoff_women_3 & total_MET <= cutoff_women_4 ~ 4,
  Sex == 1 & total_MET> cutoff_women_4 ~ 5,
))
incident_followup$total_MET_5 <- factor(incident_followup$total_MET_5)
table(incident_followup$total_MET_5)

cutoff <- incident_followup %>%
  group_by(Sex) %>%
  summarise(cut_20 = quantile(METh_modvig, 0.20),
            cut_40 = quantile(METh_modvig, 0.40),
            cut_60 = quantile(METh_modvig, 0.60),
            cut_80 = quantile(METh_modvig, 0.82))
cutoff_men_1 <- cutoff$cut_20[1]
cutoff_women_1 <- cutoff$cut_20[2]
cutoff_men_2 <- cutoff$cut_40[1]
cutoff_women_2 <- cutoff$cut_40[2]
cutoff_men_3 <- cutoff$cut_60[1]
cutoff_women_3 <- cutoff$cut_60[2]
cutoff_men_4 <- cutoff$cut_80[1]
cutoff_women_4 <- cutoff$cut_80[2]
incident_followup <- mutate(incident_followup, METh_modvig_5 = case_when(
  Sex == 0 & METh_modvig<= cutoff_men_1 ~ 1,
  Sex == 0 & METh_modvig> cutoff_men_1 & METh_modvig <= cutoff_men_2 ~ 2,
  Sex == 0 & METh_modvig> cutoff_men_2 & METh_modvig <= cutoff_men_3 ~ 3,
  Sex == 0 & METh_modvig> cutoff_men_3 & METh_modvig <= cutoff_men_4 ~ 4,
  Sex == 0 & METh_modvig> cutoff_men_4 ~ 5,
  Sex == 1 & METh_modvig<= cutoff_women_1 ~ 1,
  Sex == 1 & METh_modvig> cutoff_women_1 & METh_modvig <= cutoff_women_2 ~ 2,
  Sex == 1 & METh_modvig> cutoff_women_2 & METh_modvig <= cutoff_women_3 ~ 3,
  Sex == 1 & METh_modvig> cutoff_women_3 & METh_modvig <= cutoff_women_4 ~ 4,
  Sex == 1 & METh_modvig> cutoff_women_4 ~ 5
))
incident_followup$METh_modvig_5 <- factor(incident_followup$METh_modvig_5)
table(incident_followup$METh_modvig_5)

cutoff <- incident_followup %>%
  group_by(Sex) %>%
  summarise(cut_20 = quantile(METh_moderate, 0.20),
            cut_40 = quantile(METh_moderate, 0.40),
            cut_60 = quantile(METh_moderate, 0.60),
            cut_80 = quantile(METh_moderate, 0.82))
cutoff_men_1 <- cutoff$cut_20[1]
cutoff_women_1 <- cutoff$cut_20[2]
cutoff_men_2 <- cutoff$cut_40[1]
cutoff_women_2 <- cutoff$cut_40[2]
cutoff_men_3 <- cutoff$cut_60[1]
cutoff_women_3 <- cutoff$cut_60[2]
cutoff_men_4 <- cutoff$cut_80[1]
cutoff_women_4 <- cutoff$cut_80[2]
incident_followup <- mutate(incident_followup, METh_moderate_5 = case_when(
  Sex == 0 & METh_moderate<= cutoff_men_1 ~ 1,
  Sex == 0 & METh_moderate> cutoff_men_1 & METh_moderate <= cutoff_men_2 ~ 2,
  Sex == 0 & METh_moderate> cutoff_men_2 & METh_moderate <= cutoff_men_3 ~ 3,
  Sex == 0 & METh_moderate> cutoff_men_3 & METh_moderate <= cutoff_men_4 ~ 4,
  Sex == 0 & METh_moderate> cutoff_men_4 ~ 5,
  Sex == 1 & METh_moderate<= cutoff_women_1 ~ 1,
  Sex == 1 & METh_moderate> cutoff_women_1 & METh_moderate <= cutoff_women_2 ~ 2,
  Sex == 1 & METh_moderate> cutoff_women_2 & METh_moderate <= cutoff_women_3 ~ 3,
  Sex == 1 & METh_moderate> cutoff_women_3 & METh_moderate <= cutoff_women_4 ~ 4,
  Sex == 1 & METh_moderate> cutoff_women_4 ~ 5,
))
incident_followup$METh_moderate_5 <- factor(incident_followup$METh_moderate_5)
table(incident_followup$METh_moderate_5)

cutoff <- incident_followup %>%
  group_by(Sex) %>%
  summarise(cut_25 = quantile(Vigorous, 0.25, na.rm = TRUE),
            cut_50 = quantile(Vigorous, 0.50, na.rm = TRUE),
            cut_75 = quantile(Vigorous, 0.75, na.rm = TRUE))
cutoff_men_1 <- cutoff$cut_25[1]
cutoff_women_1 <- cutoff$cut_25[2]
cutoff_men_2 <- cutoff$cut_50[1]
cutoff_women_2 <- cutoff$cut_50[2]
cutoff_men_3 <- cutoff$cut_75[1]
cutoff_women_3 <- cutoff$cut_75[2]
incident_followup <- mutate(incident_followup, Vigorous_5 = case_when(
  Sex == 0 & is.na(Vigorous) ~ 1,
  Sex == 0 & !is.na(Vigorous) & Vigorous < cutoff_men_1 ~ 2,
  Sex == 0 & !is.na(Vigorous) & Vigorous >= cutoff_men_1 & Vigorous < cutoff_men_2 ~ 3,
  Sex == 0 & !is.na(Vigorous) & Vigorous >= cutoff_men_2 & Vigorous < cutoff_men_3 ~ 4,
  Sex == 0 & !is.na(Vigorous) & Vigorous >= cutoff_men_3 ~ 5,
  Sex == 1 & is.na(Vigorous) ~ 1,
  Sex == 1 & !is.na(Vigorous) & Vigorous < cutoff_women_1 ~ 2,
  Sex == 1 & !is.na(Vigorous) & Vigorous >= cutoff_women_1 & Vigorous < cutoff_women_2 ~ 3,
  Sex == 1 & !is.na(Vigorous) & Vigorous >= cutoff_women_2 & Vigorous < cutoff_women_3 ~ 4,
  Sex == 1 & !is.na(Vigorous) & Vigorous >= cutoff_women_3 ~ 5
))
incident_followup$Vigorous_5 <- factor(incident_followup$Vigorous_5)
table(incident_followup$Vigorous_5)

dput(names(incident_followup))
length(unique(incident_followup$ergoid))
imp0 <- mice(incident_followup, maxit=0, defaultMethod = c("pmm","logreg","polr","polyreg"))# pmm for any variable, logreg for binary, polr for ordered, polyreg for unordered
meth <- imp0$method
meth[c("plaque", "ergoid", "side_name", "side", 
       "ergocar_date", "ergocarfup_date", "birth_date",  "income_cate", "CRP", "ethnicity",
       "Follow", "scan", "status", "status_name", "METh_modvig", "income", "total_guideline", "METh_modvig_guideline",
       "stenosis",  "alcohol",  "Vigorous", 
       "total_MET_3", "METh_moderate_3", "Vigorous_3", "total_MET_4", 
       "METh_moderate_4", "Vigorous_4", "total_MET_5", "METh_moderate_5", 
       "Vigorous_5","METh_modvig_3","METh_modvig_4","METh_modvig_5")] <- "" #variables not being imputed will also not be used as predictors
# predictor matrix
pred <- imp0$predictorMatrix
pred[ ,c("plaque", "ergoid", "side_name", "side", 
         "ergocar_date", "ergocarfup_date", "birth_date",  "income_cate", "CRP", "ethnicity",
         "Follow", "scan", "status", "status_name", "METh_modvig", "income", "total_guideline", "METh_modvig_guideline",
         "stenosis", "alcohol",  "Vigorous", 
         "total_MET_3", "METh_moderate_3", "Vigorous_3", "total_MET_4", 
         "METh_moderate_4", "Vigorous_4", "total_MET_5", "METh_moderate_5", 
         "Vigorous_5","METh_modvig_3","METh_modvig_4","METh_modvig_5")] <- 0 #will not be used as predictor in any model
incident_followup_1 <- mice(incident_followup, m=30, maxit=30, meth=meth, seed=619, predictorMatrix = pred)
incident_followup_1_long <- complete(incident_followup_1, action = 'long', include = T)
incident_followup_1_check <- complete(incident_followup_1, action = 'long')
summary(incident_followup_1_check)

#cleaning the covariates
incident_followup_1_long$rs_cohort <- factor(incident_followup_1_long$rs_cohort)
#binary smoking
incident_followup_1_long <- mutate(incident_followup_1_long,
                                   smokingstatus_2 = case_when(
                                     smokingstatus == "Never_smoke" ~ 0,
                                     smokingstatus == "Past_smoke" ~ 1, 
                                     smokingstatus == "Current_smoke" ~ 1))
incident_followup_1_long$smokingstatus_2 <- factor(incident_followup_1_long$smokingstatus_2, levels = c(0,1), labels = c("Never_Past", "Current_smoke"))
#binary education
incident_followup_1_long <- mutate(incident_followup_1_long,
                                   edu_2 = case_when(
                                     edu == "primary education"  | edu == "lower/intermediate general education OR lower vocational education" | edu == "intermediate vocational education OR higher general education"  ~ 0,
                                     edu == "higher vocational education OR university"~ 1))
incident_followup_1_long$edu_2 <- factor(incident_followup_1_long$edu_2, levels = c(0,1), labels = c("low", "high"))
#hypertension
incident_followup_1_long <- mutate(incident_followup_1_long, hypertension = case_when(SBP >= 140 | DBP >= 90 | antihypertensive == "User" ~ 1,TRUE ~ 0))
incident_followup_1_long$hypertension <- factor(incident_followup_1_long$hypertension, labels = c("normal","hypertension"))
#lipid lowering
incident_followup_1_long <- mutate(incident_followup_1_long, hyperchol = case_when(TC >= 6.2 | lipidreducing == "User" ~ 1,TRUE ~ 0))
incident_followup_1_long$hyperchol <- factor(incident_followup_1_long$hyperchol, labels = c("normal","Hypercholesterolaemia"))

####baseline characteristcs by group----
data_baseline_compare <- subset(incident_followup_1_long, .imp == 1 & scan == 1)
dim(data_baseline_compare)
data_baseline_compare$Sex <- factor(data_baseline_compare$Sex)
data_baseline_compare$IPH <- factor(data_baseline_compare$IPH)
data_baseline_compare$Lipid <- factor(data_baseline_compare$Lipid)

data_baseline_compare_person <- data_baseline_compare %>% 
  distinct(ergoid, .keep_all = TRUE)
dim(data_baseline_compare_person)
dput(names(data_baseline_compare_person))

#total PA
table_1_0 <- table1(~ Sex + First_age + BMI + smokingstatus
                    + edu_2 + hyperchol + hypertension + dm | total_MET_5,
                    data=data_baseline_compare_person,
                    overall=F,
                    render.continuous=c(.="Mean (SD)"))
table_1_1 <- as.data.frame(table_1_0)
table_1_1
table_1_0 <- table1(~ Maximum + IPH + Lipid | total_MET_5,
                    data=data_baseline_compare,
                    overall=F,
                    render.continuous=c(.="Mean (SD)"))
table_1_2 <- as.data.frame(table_1_0)
table_1_2
table_totalPA <- rbind(table_1_1, table_1_2)


#modvig PA
table_1_0 <- table1(~ Sex + First_age + BMI + smokingstatus
                    + edu_2 + hyperchol + hypertension + dm | METh_modvig_5,
                    data=data_baseline_compare_person,
                    overall=F,
                    render.continuous=c(.="Mean (SD)"))
table_1_1 <- as.data.frame(table_1_0)
table_1_1
table_1_0 <- table1(~ Maximum + IPH + Lipid| METh_modvig_5,
                    data=data_baseline_compare,
                    overall=F,
                    render.continuous=c(.="Mean (SD)"))
table_1_2 <- as.data.frame(table_1_0)
table_1_2
table_modvigPA <- rbind(table_1_1, table_1_2)


#moderate PA
table_1_0 <- table1(~ Sex + First_age + BMI + smokingstatus
                    + edu_2 + hyperchol + hypertension + dm | METh_moderate_5,
                    data=data_baseline_compare_person,
                    overall=F,
                    render.continuous=c(.="Mean (SD)"))
table_1_1 <- as.data.frame(table_1_0)
table_1_1
table_1_0 <- table1(~ Maximum + IPH + Lipid | METh_moderate_5,
                    data=data_baseline_compare,
                    overall=F,
                    render.continuous=c(.="Mean (SD)"))
table_1_2 <- as.data.frame(table_1_0)
table_1_2
table_moderatePA <- rbind(table_1_1, table_1_2)


#vigorous PA
table_1_0 <- table1(~ Sex + First_age + BMI + smokingstatus
                    + edu_2 + hyperchol + hypertension + dm | Vigorous_5,
                    data=data_baseline_compare_person,
                    overall=F,
                    render.continuous=c(.="Mean (SD)"))
table_1_1 <- as.data.frame(table_1_0)
table_1_1
table_1_0 <- table1(~ Maximum + IPH + Lipid | Vigorous_5,
                    data=data_baseline_compare,
                    overall=F,
                    render.continuous=c(.="Mean (SD)"))
table_1_2 <- as.data.frame(table_1_0)
table_1_2
table_vigorousPA <- rbind(table_1_1, table_1_2)

table_totalPA_1 <- table_totalPA[c(1,4,6,8:12,15,18,21,24,25,27,30,33),]
table_modvigPA_1 <- table_modvigPA[c(1,4,6,8:12,15,18,21,24,25,27,30,33),]
table_moderatePA_1 <- table_moderatePA[c(1,4,6,8:12,15,18,21,24,25,27,30,33),]
table_vigorousPA_1 <- table_vigorousPA[c(1,4,6,8:12,15,18,21,24,25,27,30,33),]

table_PAgroup <- rbind(table_totalPA_1, table_modvigPA_1, table_moderatePA_1, table_vigorousPA_1)
table_PAgroup

####follow_ID----
incident_followup_ID <- left_join(data_total_twoscan, data_PA, by = "ergoid")
Stroke <- read_sav('Stroke_data_Luo_Oct2024.sav')
Stroke_0 <- subset(Stroke, select = c("ergoid","IC", "prevalent_stroke","Baseline_date","incident_stroke","stroke_date"))
Stroke_1 <- rename(Stroke_0, ic_Stroke = IC, prev_stroke_date = Baseline_date, prev_Stroke = prevalent_stroke, inc_stroke = incident_stroke, inc_stroke_date = stroke_date)
Stroke_1 <- subset(Stroke_1, ic_Stroke == 1)
incident_followup_ID <- left_join(incident_followup_ID, Stroke_1[c("ergoid", "prev_Stroke")], by = "ergoid")
incident_followup_ID <- subset(incident_followup_ID, !is.na(prev_Stroke)&!is.na(Cal))
incident_followup_ID <- subset(incident_followup_ID, prev_Stroke == 0)
incident_followup_ID <- subset(incident_followup_ID, !is.na(total_MET))
summary(incident_followup_ID)
ID_use <- subset(incident_followup_ID, select = c("ergoid", "side"), scan == 1)
ID_use$follow <- 1
ID_use
ID_follow <- ID_use %>%
  distinct(ergoid, .keep_all = TRUE)
ID_follow <- subset(ID_follow, select = c("ergoid", "follow"))
