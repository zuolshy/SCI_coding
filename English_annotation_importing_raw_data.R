library(readxl)
library(pubh)
library(psych)
library(plyr)
library(dplyr)
library(haven)
library(lubridate)
library(scales)
library(nnet)
library(mice)
library(grid)
library(futile.logger)
library(ggplot2)
library(questionr)
library(gnm)
library(multgee)
library(geepack)
library(broom.mixed)
library(descriptr)
library(geepack)
library(gtools)
library(labelled)
library(GPArotation)
library(sjmisc)
library(splines)
library(car)
library(table1)
library(survival)
library(JM)
library(coxme)
library(cmprsk)
library(egg)
library(patchwork)
library(contsurvplot)
library(pammtools)
library(gtable)
library(prodlim)
library(riskRegression)
library(plotly)
library(webshot2)
library(png)
library(grid)
library(imager)
library(networkD3)
rm(list=ls(all=TRUE))
options(scipen = 999)
options(digits = 5)

setwd("~")

#################data_total_twoscan----
data_baseline <- read_sav('basisbestand MRI carotid_components_10032021.sav')
dput(names(data_baseline))
data_baseline_R <- subset(data_baseline, select = c("Calc_R", "IPH_R", "LRNC_R", "Max_IMT_R", "stenosis_R", "Plaque_R", "ergoid", "Carotid_scan_new"))
data_baseline_R$side_name <- "Right"
data_baseline_R$side <- 0
data_baseline_R <- rename(data_baseline_R, Cal_1=Calc_R, Lipid_1=LRNC_R, IPH_1=IPH_R, Maximum_1=Max_IMT_R, stenosis_1=stenosis_R, plaque_1=Plaque_R)
data_baseline_L <- subset(data_baseline, select = c("Calc_L", "IPH_L", "LRNC_L", "Max_IMT_L", "stenosis_L", "Plaque_L", "ergoid", "Carotid_scan_new"))
data_baseline_L$side_name <- "Left"
data_baseline_L$side <- 1
data_baseline_L <- rename(data_baseline_L, Cal_1=Calc_L, Lipid_1=LRNC_L, IPH_1=IPH_L, Maximum_1=Max_IMT_L, stenosis_1=stenosis_L, plaque_1=Plaque_L)
data_baseline_1 <- rbind(data_baseline_R,data_baseline_L)
length(unique(data_baseline_1$ergoid))

data_secondscan <- read_excel('Second scan.xlsx')
names(data_secondscan)
data_Right <- subset(data_secondscan, select = c("ID", "Right_plaque", 
                                                 "Right_CCA_Cal", "Right_ICA_Cal",
                                                 "Right_CCA_Haem", "Right_ICA_Haem",
                                                 "Right_CCA_Lipid", "Right_ICA_Lipid",
                                                 "Right_Maximum", "Right_stenosis"))
data_Right <- rename(data_Right, 
                     Second_ID = ID,
                     plaque = Right_plaque, 
                     CCA_Cal = Right_CCA_Cal,
                     ICA_Cal = Right_ICA_Cal,
                     CCA_Haem = Right_CCA_Haem,
                     ICA_Haem = Right_ICA_Haem,
                     CCA_Lipid = Right_CCA_Lipid,
                     ICA_Lipid = Right_ICA_Lipid,
                     Maximum = Right_Maximum,
                     stenosis = Right_stenosis)
data_Right$side_name <- "Right"
data_Right$side <- 0

data_Left <- subset(data_secondscan, select = c("ID", "Left_plaque", 
                                                "Left_CCA_Cal", "Left_ICA_Cal",
                                                "Left_CCA_Haem", "Left_ICA_Haem",
                                                "Left_CCA_Lipid", "Left_ICA_Lipid",
                                                "Left_Maximum", "Left_stenosis"))
data_Left <- rename(data_Left,
                    Second_ID = ID,
                    plaque = Left_plaque , 
                    CCA_Cal = Left_CCA_Cal,
                    ICA_Cal = Left_ICA_Cal,
                    CCA_Haem = Left_CCA_Haem,
                    ICA_Haem = Left_ICA_Haem,
                    CCA_Lipid = Left_CCA_Lipid,
                    ICA_Lipid = Left_ICA_Lipid,
                    Maximum = Left_Maximum,
                    stenosis = Left_stenosis)
data_Left$side_name <- "Left"
data_Left$side <- 1

data_secondscan_0 <- rbind(data_Right, data_Left)
summary(data_secondscan_0)
names(data_secondscan_0)
data_secondscan_0$stenosis[data_secondscan_0$stenosis == 888] <- 1#recode total occlusion to 1


data_secondscan_0$Cal_Raw <- data_secondscan_0$CCA_Cal + data_secondscan_0$ICA_Cal
data_secondscan_0$Lipid_Raw <- data_secondscan_0$CCA_Lipid + data_secondscan_0$ICA_Lipid
data_secondscan_0$IPH_Raw <- data_secondscan_0$CCA_Haem + data_secondscan_0$ICA_Haem

data_secondscan_1 <- mutate(data_secondscan_0,
                            Cal = case_when(
                              Cal_Raw == 0 ~ 0,
                              Cal_Raw == 1 ~ 1,
                              Cal_Raw == 2 ~ 1,
                              TRUE ~ NA_real_
                            ))  

data_secondscan_1 <- mutate(data_secondscan_1,
                            Lipid = case_when(
                              Lipid_Raw == 0 ~ 0,
                              Lipid_Raw == 1 ~ 1,
                              Lipid_Raw == 2 ~ 1,
                              TRUE ~ NA_real_
                            )) 

data_secondscan_1 <- mutate(data_secondscan_1,
                            IPH = case_when(
                              IPH_Raw == 0 ~ 0,
                              IPH_Raw == 1 ~ 1,
                              IPH_Raw == 2 ~ 1,
                              TRUE ~ NA_real_
                            ))

#merge two date of visit
ID_0 <- read_excel('ID_Final.xlsx')
ID_R <- ID_0
ID_R$side_name <- "Right"
ID_R$side <- 0
ID_L <- ID_0
ID_L$side_name <- "Left"
ID_L$side <- 1
ID <- rbind(ID_R, ID_L)

data_secondscan_2 <- left_join(ID, data_secondscan_1, by = c("Second_ID", "side", "side_name"))
length(unique(data_secondscan_2$ergoid))
length(data_secondscan_2$ergoid)
dput(names(data_secondscan_2))
data_secondscan_3 <- subset(data_secondscan_2, select = c("Second_ID", "plaque", "Maximum", "stenosis", "side_name", 
                                                          "side", "Cal", "Lipid", "IPH","ergoid", "First_ID", "ergocar_date", "ergocarfup_date", "birth_date", 
                                                          "First_age", "Follow", "Sex"))
data_secondscan_4 <- rename(data_secondscan_3, Cal_2=Cal, Lipid_2=Lipid, IPH_2=IPH, Maximum_2=Maximum, stenosis_2=stenosis, plaque_2=plaque)
summary(data_secondscan_4)
length(unique(data_secondscan_4$ergoid))
length(data_secondscan_4$ergoid)

data_secondscan_5 <- left_join(data_baseline_1, data_secondscan_4, by = c("ergoid", "side", "side_name"))
length(unique(data_secondscan_5$ergoid))
length(data_secondscan_5$ergoid)
summary(data_secondscan_5)

data_secondscan_5 <- data_secondscan_5

data_secondscan_5 <- mutate(data_secondscan_5,
                            status_1 = case_when(
                              Cal_1 == 0 & IPH_1 == 0 & Lipid_1 == 0  ~ 8,
                              Cal_1 == 1 & IPH_1 == 0 & Lipid_1 == 0  ~ 1, 
                              Cal_1 == 0 & IPH_1 == 1 & Lipid_1 == 0  ~ 3, 
                              Cal_1 == 0 & IPH_1 == 0 & Lipid_1 == 1  ~ 2, 
                              Cal_1 == 1 & IPH_1 == 1 & Lipid_1 == 0  ~ 5, 
                              Cal_1 == 1 & IPH_1 == 0 & Lipid_1 == 1  ~ 4, 
                              Cal_1 == 0 & IPH_1 == 1 & Lipid_1 == 1  ~ 6, 
                              Cal_1 == 1 & IPH_1 == 1 & Lipid_1 == 1  ~ 7,
                              TRUE ~ NA_real_))
data_secondscan_5 <- mutate(data_secondscan_5,
                            status_2 = case_when(
                              Cal_2 == 0 & IPH_2 == 0 & Lipid_2 == 0  ~ 8,
                              Cal_2 == 1 & IPH_2 == 0 & Lipid_2 == 0  ~ 1, 
                              Cal_2 == 0 & IPH_2 == 1 & Lipid_2 == 0  ~ 3, 
                              Cal_2 == 0 & IPH_2 == 0 & Lipid_2 == 1  ~ 2, 
                              Cal_2 == 1 & IPH_2 == 1 & Lipid_2 == 0  ~ 5, 
                              Cal_2 == 1 & IPH_2 == 0 & Lipid_2 == 1  ~ 4, 
                              Cal_2 == 0 & IPH_2 == 1 & Lipid_2 == 1  ~ 6, 
                              Cal_2 == 1 & IPH_2 == 1 & Lipid_2 == 1  ~ 7,
                              TRUE ~ NA_real_))
data_secondscan_5 <- data_secondscan_5[!is.na(data_secondscan_5$status_1),]
data_secondscan_5 <- data_secondscan_5[!is.na(data_secondscan_5$status_2),]

dput(names(data_secondscan_5))
data_total_firstscan <- subset(data_secondscan_5,select = c("Cal_1", "IPH_1", "Lipid_1", "Maximum_1", "stenosis_1", "plaque_1", 
                                                            "ergoid", "side_name", "side", "ergocar_date",
                                                            "ergocarfup_date", "birth_date", "First_age", "Follow", "Sex"))
data_total_firstscan$scan <- 1
data_total_firstscan <- rename(data_total_firstscan, Cal=Cal_1, Lipid=Lipid_1, IPH=IPH_1, Maximum=Maximum_1, stenosis=stenosis_1, plaque=plaque_1)

data_total_secondscan <- subset(data_secondscan_5,select = c("Cal_2", "IPH_2", "Lipid_2", "Maximum_2", "stenosis_2", "plaque_2", 
                                                             "ergoid", "side_name", "side", "ergocar_date",
                                                             "ergocarfup_date", "birth_date", "First_age", "Follow", "Sex"))
data_total_secondscan$scan <- 2
data_total_secondscan <- rename(data_total_secondscan, Cal=Cal_2, Lipid=Lipid_2, IPH=IPH_2, Maximum=Maximum_2, stenosis=stenosis_2, plaque=plaque_2)

data_total <- rbind(data_total_firstscan,data_total_secondscan)

data_total <- mutate(data_total,
                     status = case_when(
                       Cal == 0 & IPH == 0 & Lipid == 0  ~ 7,
                       Cal == 0 & IPH == 0 & Lipid == 1  ~ 1, 
                       Cal == 0 & IPH == 1 & Lipid == 0  ~ 2, 
                       Cal == 0 & IPH == 1 & Lipid == 1  ~ 2, 
                       Cal == 1 & IPH == 0 & Lipid == 0  ~ 3, 
                       Cal == 1 & IPH == 0 & Lipid == 1  ~ 4, 
                       Cal == 1 & IPH == 1 & Lipid == 0  ~ 5,
                       Cal == 1 & IPH == 1 & Lipid == 1  ~ 6,
                       TRUE ~ NA_real_)) 

name_0 <- data.frame(status = c(1:7), status_name = c("Lipid",
                                                      "IPH ± Lipid",
                                                      "Cal",
                                                      "Cal+Lipid",
                                                      "Cal+IPH",
                                                      "Cal+Lipid&IPH",
                                                      "None"))

data_total <- left_join(data_total, name_0, by = "status")
data_total$status_name <- factor(data_total$status_name, levels = unique(data_total$status_name[order(data_total$status,decreasing = FALSE)]))
data_total <- data_total[!is.na(data_total$status_name),]
table(data_total$status_name)

data_total <- mutate(data_total,
                     followuptime = case_when(
                       scan == 1  ~ 0,
                       scan == 2  ~ time_length(difftime(data_total$ergocarfup_date, data_total$ergocar_date), "years")))

data_total_twoscan <- data_total[data_total$Follow == 1,]

length(unique(data_total_twoscan$ergoid))
summary(data_total_twoscan)
length(unique(data_total_twoscan$ergoid[!is.na(data_total_twoscan$status)]))


#########################covariates----
#########################Import ID File----
ID_0 <- read_excel("ID_Final.xlsx")
Scan_ID <- ID_0
names(Scan_ID)
Scan_ID <- distinct(Scan_ID, First_ID, .keep_all= TRUE)
Scan_ID$followup <- NA
Scan_ID$followup[Scan_ID$ergocar_date<ymd("2009-01-20")] <- "ej"
Scan_ID$followup[Scan_ID$ergocar_date>ymd("2009-01-19")] <- "e5"
Scan_ID
Scan_ID$Sex <- factor(Scan_ID$Sex, levels = c(0,1), labels = c("Men", "Women"))
summary(Scan_ID)
Scan_ID_1 <- subset(Scan_ID, select = -followup)



#########################ethnicity----
#1 for European
#2 for East-Asian
#3 for African
#4 for Admixture
#9 for Missing
Ethnicity_RS_1 <- read_sav('RS_I_GeneticEtnicity_(11-JAN-2016).sav')
Ethnicity_RS_2 <- read_sav('RS_II_GeneticEtnicity_(11-JAN-2016).sav')
Ethnicity_RS_3 <- read_sav('RS_III_GeneticEtnicity_(11-JAN-2016).sav')
names(Ethnicity_RS_1)
names(Ethnicity_RS_2)
names(Ethnicity_RS_3)
Ethnicity_RS_1<-Ethnicity_RS_1[,c("ergoid","ancestry")]
Ethnicity_RS_1<-rename(Ethnicity_RS_1, ethnicity = ancestry)
Ethnicity_RS_2<-Ethnicity_RS_2[,c("ergoid","ancestry")]
Ethnicity_RS_2<-rename(Ethnicity_RS_2, ethnicity = ancestry)
Ethnicity_RS_3<-Ethnicity_RS_3[,c("ergoid","ancestry")]
Ethnicity_RS_3<-rename(Ethnicity_RS_3, ethnicity = ancestry)
Ethnicity<-rbind(Ethnicity_RS_1, Ethnicity_RS_2)
Ethnicity<-rbind(Ethnicity, Ethnicity_RS_3)
Ethnicity$ethnicity <- factor(Ethnicity$ethnicity, levels = c(1,2,3,4,9), labels = c("European","East-Asian","African","Admixture","Missing"))
summary(Ethnicity)

#########################Diabetes----
DM_raw <- read_sav('DM.sav')
names(DM_raw)
DM_prevalence <- subset(DM_raw, Inci_DM_2015 == 8, select = c("ergoid","Inci_DM_2015"))
DM_prevalence$Inci_DM_2015 <- 1
DM_prevalence<-rename(DM_prevalence, dm = Inci_DM_2015)
DM_no <- subset(DM_raw, Inci_DM_2015 == 0, select = c("ergoid","Inci_DM_2015"))
DM_no<-rename(DM_no, dm = Inci_DM_2015)
DM_incidence <- subset(DM_raw, Inci_DM_2015 == 1, select = c("ergoid","Inci_DM_2015","Incidentdate_DM"))
DM_incidence <- left_join(DM_incidence, Scan_ID, by = "ergoid")
DM_incidence <-  DM_incidence %>% mutate(dm = ifelse(Incidentdate_DM < ergocar_date,1,0))
DM_incidence <- subset(DM_incidence, select = c("ergoid","dm"))
DM<-rbind(DM_prevalence, DM_no, DM_incidence)
DM$dm <- factor(DM$dm, levels = c(0,1), labels = c("Non_diabetes", "diabetes"))

#########################CRP----
# Because each person is only tested once for CRP, we can only use the earlier visit to impute missing
ej_CRP_raw <- read_sav('ej_(1)_AKCBASEL_(21-jun-2013).sav')
e4_CRP_raw <- read_sav('e4_(4)_LABJANSS_(09-apr-2019).sav')
ep_CRP_raw <- read_sav('ep_(1)_AKCBASEL_(11-sep-2013).sav')
e3_CRP_raw <- read_sav('e3_(3)_LABMISC_(09-apr-2019).sav')
names(ej_CRP_raw)
names(e4_CRP_raw)
names(ep_CRP_raw)
names(e3_CRP_raw)
ej_CRP<-ej_CRP_raw[,c("ergoid","ej_18452")]
ej_CRP<-rename(ej_CRP, ej_CRP = ej_18452)
e4_CRP<-e4_CRP_raw[,c("ergoid","e4_18452")]
e4_CRP<-rename(e4_CRP, e4_CRP = e4_18452)
ep_CRP<-ep_CRP_raw[,c("ergoid","ep_18452")]
ep_CRP<-rename(ep_CRP, ep_CRP = ep_18452)
e3_CRP<-e3_CRP_raw[,c("ergoid","e3_18676")]
e3_CRP<-rename(e3_CRP, e3_CRP = e3_18676)
CRP <- full_join(ej_CRP, e4_CRP, by = "ergoid")%>%
  full_join(., ep_CRP, by = "ergoid")%>%
  full_join(., e3_CRP, by = "ergoid")
CRP <- mutate(CRP, CRP = coalesce(ej_CRP, e4_CRP))
CRP <- mutate(CRP, CRP = coalesce(CRP, ep_CRP))
CRP <- mutate(CRP, CRP = coalesce(CRP, e3_CRP))
CRP_1 <- subset(CRP,select = c("ergoid", "CRP"))

#########################medication eje5----
ej_medication_raw <- read_sav('ej_MEDICATION_(22-mar-2010).sav')
e5_medication_raw <- read_sav('e5_MEDICATION_(12-MAR-2015).sav')
names(ej_medication_raw)
names(e5_medication_raw)
ej_medication<-ej_medication_raw[,c("ergoid","ej_c10", "ej_b01", "ej_a10")]
ej_medication<-rename(ej_medication, lipidreducing = ej_c10, antithrombotic = ej_b01, antidm = ej_a10)
e5_medication<-e5_medication_raw[,c("ergoid","e5_c10", "e5_b01", "e5_a10")]
e5_medication<-rename(e5_medication, lipidreducing = e5_c10, antithrombotic = e5_b01, antidm = e5_a10)
ej_medication$followup<-"ej"
e5_medication$followup<-"e5"
medication<-rbind(ej_medication, e5_medication)
#降压药选用了相对准确的基线数据
HT_raw <- read_sav('HT2018_analysisfile_(15-may-2018).sav')
names(HT_raw)
ej_HT<-HT_raw[,c("ergoid","ej_bpldrug")]
ej_HT<-rename(ej_HT, antihypertensive = ej_bpldrug)
e5_HT<-HT_raw[,c("ergoid","e5_bpldrug")]
e5_HT<-rename(e5_HT, antihypertensive = e5_bpldrug)
ej_HT$followup<-"ej"
e5_HT$followup<-"e5"
HT <- rbind(ej_HT, e5_HT)
medication <- left_join(HT, medication, by = c("ergoid","followup"))
unique(medication$antihypertensive)#9missing
unique(medication$lipidreducing)#9missing
unique(medication$antithrombotic)#9missing
unique(medication$antidm)#9missing
medication <-  medication %>% mutate(antihypertensive = replace(antihypertensive, antihypertensive==9, NA),
                                     lipidreducing = replace(lipidreducing, lipidreducing==9, NA),
                                     antithrombotic = replace(antithrombotic, antithrombotic==9, NA),
                                     antidm = replace(antidm, antithrombotic==9, NA))
unique(medication$antihypertensive)#9missing
unique(medication$lipidreducing)#9missing
unique(medication$antithrombotic)#9missing
unique(medication$antidm)#9missing
medication$antidm <- factor(medication$antidm, levels = c(0,1), labels = c("Non_user", "User"))
medication$antithrombotic <- factor(medication$antithrombotic, levels = c(0,1), labels = c("Non_user", "User"))
medication$lipidreducing <- factor(medication$lipidreducing, levels = c(0,1), labels = c("Non_user", "User"))
medication$antihypertensive <- factor(medication$antihypertensive, levels = c(0,1), labels = c("Non_user", "User"))
summary(medication)

#########################smoking eje5----
ej_smoking_raw <- read_sav('ej_intvw_SMOKING_(28-mar-2011).sav')
e5_smoking_raw <- read_sav('e5_intvw_SMOKING_(22-mar-2021).sav')
names(ej_smoking_raw)
names(e5_smoking_raw)
ej_smoking<-ej_smoking_raw[,c("ergoid","ej_yilf6","ej_yilfe")]
ej_smoking<-rename(ej_smoking, currentsmoking = ej_yilf6, pastsmoking = ej_yilfe)
e5_smoking<-e5_smoking_raw[,c("ergoid","e5_EILF6","e5_EILFE")]
e5_smoking<-rename(e5_smoking, currentsmoking = e5_EILF6, pastsmoking = e5_EILFE)
ej_smoking$followup<-"ej"
e5_smoking$followup<-"e5"
smoking<-rbind(ej_smoking, e5_smoking)
unique(smoking$currentsmoking)#9no answer，7dont know
unique(smoking$pastsmoking)#9no answer，7dont know, 8no applicable
smoking <-  smoking %>% mutate(currentsmoking = replace(currentsmoking, currentsmoking==7, NA),
                               currentsmoking = replace(currentsmoking, currentsmoking==9, NA),
                               pastsmoking = replace(pastsmoking, pastsmoking==7, NA),
                               pastsmoking = replace(pastsmoking, pastsmoking==9, NA),
                               pastsmoking = replace(pastsmoking, pastsmoking==8, 0))
unique(smoking$currentsmoking)#9no answer，7dont know
unique(smoking$pastsmoking)#9no answer，7dont know, 8no applicable

smoking <- mutate(smoking,
                  smokingstatus = case_when(
                    currentsmoking == 1 ~ 1, 
                    pastsmoking == 1 & currentsmoking == 0 ~ 2, 
                    currentsmoking == 0 & pastsmoking == 0 ~ 0,
                    TRUE ~ NA_real_))
smoking$smokingstatus <- factor(smoking$smokingstatus, levels = c(0,1,2), labels = c("Never_smoke", "Current_smoke", "Past_smoke"))
summary(smoking)

#########################alcohol eje5----
ej_alcohol_raw <- read_sav('ej_intvw_ALCOHOL_(28-mar-2011).sav')
e5_alcohol_raw <- read_sav('e5_intvw_ALCOHOL_(18-apr-2016).sav')
names(ej_alcohol_raw)
names(e5_alcohol_raw)
ej_alcohol<-ej_alcohol_raw[,c("ergoid","ej_yiaudit2")]
ej_alcohol<-rename(ej_alcohol, alcohol = ej_yiaudit2)
e5_alcohol<-e5_alcohol_raw[,c("ergoid","e5_EIAUDIT2")]
e5_alcohol<-rename(e5_alcohol, alcohol = e5_EIAUDIT2)
ej_alcohol$followup<-"ej"
e5_alcohol$followup<-"e5"
alcohol<-rbind(ej_alcohol, e5_alcohol)
alcohol$alcohol <- alcohol$alcohol + 1
unique(alcohol$alcohol)#9no answer，7dont know, 410 galss 以上，37 - 9 glass，25 - 6 glass，13 - 4 glass，01 - 2 glass
alcohol <-  alcohol %>% mutate(alcohol = replace(alcohol, alcohol==8, NA),
                               alcohol = replace(alcohol, alcohol==10, NA),
                               alcohol = replace(alcohol, alcohol==9, 0))
alcohol$alcohol <- factor(alcohol$alcohol, order = TRUE, labels = c("0", "1 - 2", "3 - 4", "5 - 6", "7 - 9", "> 10"))
summary(alcohol)

#########################bp eje5----
bp_raw <- read_sav('HT2018_analysisfile_(15-may-2018).sav')
names(bp_raw)
ej_bp<-bp_raw[,c("ergoid","ej_systolicBP","ej_diastolicBP")]
ej_bp<-rename(ej_bp, SBP = ej_systolicBP, DBP = ej_diastolicBP)
e5_bp<-bp_raw[,c("ergoid","e5_systolicBP","e5_diastolicBP")]
e5_bp<-rename(e5_bp, SBP = e5_systolicBP, DBP = e5_diastolicBP)
ej_bp$followup<-"ej"
e5_bp$followup<-"e5"
bp<-rbind(ej_bp, e5_bp)
summary(bp)#999missing
bp <-  bp %>% mutate(SBP = replace(SBP, SBP==999, NA),
                     DBP = replace(DBP, DBP==999, NA))
summary(bp)

#########################lipid + glucose eje5----
ej_lipid_raw <- read_sav('ej_(1)_LAB_(11-jun-2009)r.sav')
e5_lipid_raw <- read_sav('e5_(5)_LAB_(29-aug-2014)r.sav')
names(ej_lipid_raw)
names(e5_lipid_raw)
ej_lipid<-ej_lipid_raw[,c("ergoid","ej_3845","ej_4107","ej_15519","ej_3846")]
ej_lipid<-rename(ej_lipid, TC = ej_3845, HDLC = ej_4107, LDLC = ej_15519, GLU = ej_3846)
e5_lipid<-e5_lipid_raw[,c("ergoid","e5_3845","e5_4107","e5_15519","e5_3846")]
e5_lipid<-rename(e5_lipid, TC = e5_3845, HDLC = e5_4107, LDLC = e5_15519, GLU = e5_3846)
ej_lipid$followup<-"ej"
e5_lipid$followup<-"e5"
lipid<-rbind(ej_lipid, e5_lipid)
summary(lipid)

#########################BMI eje5----
ej_BMI_raw <- read_sav('ej_(1)_UITSCHR_(23-feb-2010)_ANTHROPO-PART.sav')
e5_BMI_raw <- read_sav('e5_(5)_ANTHROPO_(10-dec-2015).sav')
names(ej_BMI_raw)
names(e5_BMI_raw)
ej_BMI<-ej_BMI_raw[,c("ergoid","ej_229","ej_230")]
ej_BMI<-rename(ej_BMI, height_cm = ej_229, weight = ej_230)
e5_BMI<-e5_BMI_raw[,c("ergoid","e5_229","e5_230")]
e5_BMI<-rename(e5_BMI, height_cm = e5_229, weight = e5_230)
ej_BMI$followup<-"ej"
e5_BMI$followup<-"e5"
BMI<-rbind(ej_BMI, e5_BMI)
BMI$BMI <- BMI$weight / (BMI$height_cm/100)^2
summary(BMI)#999.9 = missing, 888.8 = not appropriate-default


#########################WHR eje5----
ej_WHR_raw <- read_sav('ej_(1)_LNGNPS_(18-oct-2010)_ANTHROPO-PART.sav')
e5_WHR_raw <- read_sav('e5_(5)_ANTHROPO_(10-dec-2015).sav')
names(ej_WHR_raw)
names(e5_WHR_raw)
ej_WHR<-ej_WHR_raw[,c("ergoid","ej_231","ej_232")]
ej_WHR<-rename(ej_WHR, waist_cm = ej_231, hip_cm = ej_232)
ej_WHR$WHR<-ej_WHR$waist_cm / ej_WHR$hip_cm
e5_WHR<-e5_WHR_raw[,c("ergoid","e5_231","e5_232")]
e5_WHR<-rename(e5_WHR, waist_cm = e5_231, hip_cm = e5_232)
e5_WHR$WHR<-e5_WHR$waist_cm / e5_WHR$hip_cm
ej_WHR$followup<-"ej"
e5_WHR$followup<-"e5"
WHR<-rbind(ej_WHR, e5_WHR)
summary(WHR)#999.9 = missing, 888.8 = not appropriate-default

#########################use ej impute e5----
Scan_ID
Scan_avail <- subset(Scan_ID, select = c("ergoid", "followup"))
Scan_covar <- left_join(Scan_avail, medication, by = c("ergoid","followup")) %>%
  left_join(., smoking, by = c("ergoid","followup"))%>%
  left_join(., alcohol, by = c("ergoid","followup"))%>%
  left_join(., bp, by = c("ergoid","followup"))%>%
  left_join(., lipid, by = c("ergoid","followup"))%>%
  left_join(., BMI, by = c("ergoid","followup"))%>%
  left_join(., WHR, by = c("ergoid","followup"))

dput(names(Scan_covar))
Scan_covar <- subset(Scan_covar, select = c("ergoid", "followup", "antihypertensive", "lipidreducing", "antithrombotic", "antidm", "smokingstatus", "alcohol", 
                                            "SBP", "DBP", "TC", "HDLC", "GLU", "BMI", "WHR", "waist_cm","height_cm"))

Scan_covar_complete <- Scan_covar[rowSums(is.na(Scan_covar)) == 0, ]#extract complete cases
Scan_covar_NA <- Scan_covar[rowSums(is.na(Scan_covar)) > 0, ]#extract rows with NA

Scan_covar_NA_ej <- Scan_covar_NA[Scan_covar_NA$followup == "ej",]
Scan_covar_NA_e5 <- Scan_covar_NA[Scan_covar_NA$followup == "e5",]

covar <- left_join(medication, smoking, by = c("ergoid","followup")) %>%
  left_join(., alcohol, by = c("ergoid","followup"))%>%
  left_join(., bp, by = c("ergoid","followup"))%>%
  left_join(., lipid, by = c("ergoid","followup"))%>%
  left_join(., BMI, by = c("ergoid","followup"))%>%
  left_join(., WHR, by = c("ergoid","followup"))

dput(names(covar))
covar <- subset(covar, select = c("ergoid", "antihypertensive", "followup", "lipidreducing", 
                                  "antithrombotic", "antidm", "smokingstatus", "alcohol", "SBP", "DBP", "TC", "HDLC","GLU",   
                                  "BMI", "WHR", "waist_cm","height_cm"))
covar_ej <- covar[covar$followup == "ej",]
covar_e5 <- covar[covar$followup == "e5",]

####use ej impute e5
Scan_covar_NA_e5 <- left_join(Scan_covar_NA_e5,covar_ej,by = "ergoid")
Scan_covar_NA_e5 <- mutate(Scan_covar_NA_e5, antihypertensive = coalesce(antihypertensive.x, antihypertensive.y))%>%
  mutate(Scan_covar_NA_e5, lipidreducing = coalesce(lipidreducing.x, lipidreducing.y))%>%
  mutate(Scan_covar_NA_e5, antithrombotic = coalesce(antithrombotic.x, antithrombotic.y))%>%
  mutate(Scan_covar_NA_e5, antidm = coalesce(antidm.x, antidm.y))%>%
  mutate(Scan_covar_NA_e5, smokingstatus = coalesce(smokingstatus.x, smokingstatus.y))%>%
  mutate(Scan_covar_NA_e5, BMI = coalesce(BMI.x, BMI.y))%>%
  mutate(Scan_covar_NA_e5, alcohol = coalesce(alcohol.x, alcohol.y))%>%
  mutate(Scan_covar_NA_e5, SBP = coalesce(SBP.x, SBP.y))%>%
  mutate(Scan_covar_NA_e5, DBP = coalesce(DBP.x, DBP.y))%>%
  mutate(Scan_covar_NA_e5, TC = coalesce(TC.x, TC.y))%>%
  mutate(Scan_covar_NA_e5, HDLC = coalesce(HDLC.x, HDLC.y))%>%
  mutate(Scan_covar_NA_e5, GLU = coalesce(GLU.x, GLU.y))%>%
  mutate(Scan_covar_NA_e5, waist_cm = coalesce(waist_cm.x, waist_cm.y))%>%
  mutate(Scan_covar_NA_e5, WHR = coalesce(WHR.x, WHR.y))%>%
  mutate(Scan_covar_NA_e5, height_cm = coalesce(height_cm.x, height_cm.y))
dput(names(Scan_covar_NA_e5))
Scan_covar_NA_e5 <- subset(Scan_covar_NA_e5, select = c("ergoid", "antihypertensive", "lipidreducing", "antithrombotic", 
                                                        "antidm", "smokingstatus", "BMI", "alcohol", "SBP", 
                                                        "DBP", "TC", "HDLC", "GLU", "waist_cm", "WHR","height_cm"))
Scan_covar_NA_e5$followup <- "impute_with_ej"
Scan_covar_imputed <- rbind(Scan_covar_complete,Scan_covar_NA_e5)
summary(Scan_covar_imputed)


#########################use e4 imput ej missing----
#####medication e4
e4_medication_raw <- read_sav('e4_MEDICATION_(20-jan-2010).sav')
names(e4_medication_raw)
e4_medication<-e4_medication_raw[,c("ergoid","e4_c10", "e4_b01", "e4_a10")]
e4_medication<-rename(e4_medication, lipidreducing = e4_c10, antithrombotic = e4_b01, antidm = e4_a10)
e4_medication$followup<-"e4"
names(HT_raw)
e4_HT<-HT_raw[,c("ergoid","e4_bpldrug")]
e4_HT<-rename(e4_HT, antihypertensive = e4_bpldrug)
e4_HT$followup<-"e4"
e4_medication <- left_join(e4_HT, e4_medication, by = c("ergoid","followup"))
unique(e4_medication$antihypertensive)#9missing
unique(e4_medication$lipidreducing)#9missing
unique(e4_medication$antithrombotic)#9missing
unique(e4_medication$antidm)#9missing
e4_medication <-  e4_medication %>% mutate(antihypertensive = replace(antihypertensive, antihypertensive==9, NA),
                                           lipidreducing = replace(lipidreducing, lipidreducing==9, NA),
                                           antithrombotic = replace(antithrombotic, antithrombotic==9, NA),
                                           antidm = replace(antidm, antithrombotic==9, NA)
)
unique(e4_medication$antihypertensive)#9missing
unique(e4_medication$lipidreducing)#9missing
unique(e4_medication$antithrombotic)#9missing
unique(e4_medication$antidm)#9missing
e4_medication$antidm <- factor(e4_medication$antidm, levels = c(0,1), labels = c("Non_user", "User"))
e4_medication$antithrombotic <- factor(e4_medication$antithrombotic, levels = c(0,1), labels = c("Non_user", "User"))
e4_medication$lipidreducing <- factor(e4_medication$lipidreducing, levels = c(0,1), labels = c("Non_user", "User"))
e4_medication$antihypertensive <- factor(e4_medication$antihypertensive, levels = c(0,1), labels = c("Non_user", "User"))
summary(e4_medication)

#####smoking e4
e4_smoking_raw <- read_sav('e4_intvw_SMOKING_(04-nov-2011).sav')
names(e4_smoking_raw)
e4_smoking<-e4_smoking_raw[,c("ergoid","e4_dict")]
e4_smoking<-rename(e4_smoking, smoking = e4_dict)
e4_smoking$followup<-"e4"
unique(e4_smoking$smoking)#0never smoke，1current smoke, 2&3past smoke, 7dont know, 9missing
e4_smoking <- mutate(e4_smoking,
                     smokingstatus = case_when(
                       smoking == 0 ~ 0,
                       smoking == 1 ~ 1,
                       smoking == 2 | smoking == 3 ~ 2,
                       TRUE ~ NA_real_
                     ))
e4_smoking$smokingstatus <- factor(e4_smoking$smokingstatus, levels = c(0,1,2), labels = c("Never_smoke", "Current_smoke", "Past_smoke"))
summary(e4_smoking)

#####alcohol e4
e4_alcohol_raw <- read_sav('e4_intvw_Alcoholperday_22-11-2013.sav')
names(e4_alcohol_raw)
dput(names(e4_alcohol_raw))
e4_alcohol<-e4_alcohol_raw[,c("ergoid", "e4_didbefr", "e4_didwwifr", "e4_didrwifr", "e4_didshfr", "e4_didlifr")]
e4_alcohol[is.na(e4_alcohol)] <- 0
e4_alcohol$alcohol_total <- (e4_alcohol$e4_didbefr + e4_alcohol$e4_didwwifr + e4_alcohol$e4_didrwifr + e4_alcohol$e4_didshfr + e4_alcohol$e4_didlifr) / 7
e4_alcohol <- mutate(e4_alcohol,
                     alcohol = case_when(
                       alcohol_total == 0  ~ 0,
                       alcohol_total >= 0.00001 & alcohol_total <= 2.0000 ~ 1,
                       alcohol_total >= 2.00001 & alcohol_total <= 4.0000 ~ 2,
                       alcohol_total >= 4.00001 & alcohol_total <= 6.0000 ~ 3,
                       alcohol_total >= 6.00001 & alcohol_total <= 8.0000 ~ 4,
                       alcohol_total >= 10 ~ 5,
                       TRUE ~ NA_real_
                     ))  
e4_alcohol$followup<-"e4"
e4_alcohol$alcohol <- factor(e4_alcohol$alcohol, order = TRUE, labels = c("0", "1 - 2", "3 - 4", "5 - 6", "7 - 9", "> 10"))
summary(e4_alcohol)

#####bp e4
names(bp_raw)
e4_bp<-bp_raw[,c("ergoid","e4_systolicBP","e4_diastolicBP")]
e4_bp<-rename(e4_bp, SBP = e4_systolicBP, DBP = e4_diastolicBP)
e4_bp$followup<-"e4"
e4_bp <-  e4_bp %>% mutate(SBP = replace(SBP, SBP==999, NA),
                           DBP = replace(DBP, DBP==999, NA)
)
summary(e4_bp)

#####lipid e4
e4_lipid_raw <- read_sav('e4_(4)_LAB_(10-mar-2010)b.sav')
names(e4_lipid_raw)
e4_lipid<-e4_lipid_raw[,c("ergoid","e4_3845","e4_4107","e4_3846")]
e4_lipid<-rename(e4_lipid, TC = e4_3845, HDLC = e4_4107, GLU = e4_3846)
e4_lipid$followup<-"e4"
summary(e4_lipid)


#####BMI e4
e4_BMI_raw <- read_sav('e4_(4)_UITSCHR_(06-nov-2014)_ANTHROPO-PART.sav')
names(e4_BMI_raw)
e4_BMI<-e4_BMI_raw[,c("ergoid","e4_229","e4_230")]
e4_BMI<-rename(e4_BMI, height_cm = e4_229, weight = e4_230)
e4_BMI$followup<-"e4"
e4_BMI$BMI <- e4_BMI$weight / (e4_BMI$height_cm/100)^2
summary(e4_BMI)

#####WHR e4
e4_WHR_raw <- read_sav('e4_(4)_HARTVAAT_(15-apr-2010)_ANTHROPO.sav')
names(e4_WHR_raw)
e4_WHR<-e4_WHR_raw[,c("ergoid","e4_231","e4_232")]
e4_WHR<-rename(e4_WHR, waist_cm = e4_231, hip_cm = e4_232)
e4_WHR$WHR<-e4_WHR$waist_cm / e4_WHR$hip_cm
e4_WHR$followup<-"e4"

#####merge e4
e4_covar <- left_join(e4_medication, e4_smoking, by = c("ergoid","followup")) %>%
  left_join(., e4_alcohol, by = c("ergoid","followup"))%>%
  left_join(., e4_bp, by = c("ergoid","followup"))%>%
  left_join(., e4_lipid, by = c("ergoid","followup"))%>%
  left_join(., e4_BMI, by = c("ergoid","followup"))%>%
  left_join(., e4_WHR, by = c("ergoid","followup"))

dput(names(e4_covar))
e4_covar <- subset(e4_covar, select = c("ergoid", "antihypertensive", "followup", "lipidreducing", "GLU",
                                        "antithrombotic", "antidm", "smokingstatus", "alcohol", "SBP", "DBP", "TC", "HDLC", 
                                        "BMI", "waist_cm", "WHR","height_cm"))

#####use e4 impute ej
Scan_covar_imputed_1 <- Scan_covar_NA_ej
Scan_covar_imputed_1 <- left_join(Scan_covar_imputed_1,e4_covar,by = "ergoid")
dput(names(Scan_covar_imputed_1))
Scan_covar_imputed_1 <- mutate(Scan_covar_imputed_1, antihypertensive = coalesce(antihypertensive.x, antihypertensive.y))%>%
  mutate(Scan_covar_imputed_1, lipidreducing = coalesce(lipidreducing.x, lipidreducing.y))%>%
  mutate(Scan_covar_imputed_1, antithrombotic = coalesce(antithrombotic.x, antithrombotic.y))%>%
  mutate(Scan_covar_imputed_1, antidm = coalesce(antidm.x, antidm.y))%>%
  mutate(Scan_covar_imputed_1, smokingstatus = coalesce(smokingstatus.x, smokingstatus.y))%>%
  mutate(Scan_covar_imputed_1, BMI = coalesce(BMI.x, BMI.y))%>%
  mutate(Scan_covar_imputed_1, alcohol = coalesce(alcohol.x, alcohol.y))%>%
  mutate(Scan_covar_imputed_1, SBP = coalesce(SBP.x, SBP.y))%>%
  mutate(Scan_covar_imputed_1, DBP = coalesce(DBP.x, DBP.y))%>%
  mutate(Scan_covar_imputed_1, TC = coalesce(TC.x, TC.y))%>%
  mutate(Scan_covar_imputed_1, HDLC = coalesce(HDLC.x, HDLC.y))%>%
  mutate(Scan_covar_imputed_1, GLU = coalesce(GLU.x, GLU.y))%>%
  mutate(Scan_covar_imputed_1, waist_cm = coalesce(waist_cm.x, waist_cm.y))%>%
  mutate(Scan_covar_imputed_1, WHR = coalesce(WHR.x, WHR.y))%>%
  mutate(Scan_covar_imputed_1, height_cm = coalesce(height_cm.x, height_cm.y))
names(Scan_covar_imputed_1)
names(Scan_covar_imputed)
Scan_covar_imputed_1 <- subset(Scan_covar_imputed_1, select = c("ergoid", "antihypertensive", "lipidreducing", "antithrombotic", 
                                                                "antidm", "smokingstatus", "BMI", "alcohol", "SBP", 
                                                                "DBP", "TC", "HDLC", "GLU", "waist_cm", "WHR","height_cm"))
Scan_covar_imputed_1$followup <- "impute_with_e4"
Scan_covar_imputed <- rbind(Scan_covar_imputed,Scan_covar_imputed_1)
summary(Scan_covar_imputed)

#final covariates dataset
Scan_ID
Scan_avail <- subset(Scan_ID, select = c("ergoid", "Sex", "First_age"))
Scan_covar <- left_join(Scan_avail, Scan_covar_imputed, by = c("ergoid"))
Scan_covar <- left_join(Scan_covar, CRP_1, by = "ergoid")%>%
  left_join(., Ethnicity, by = c("ergoid"))%>%
  left_join(., DM, by = c("ergoid"))

Scan_covar <- subset(Scan_covar, select = c("ergoid", "antihypertensive", "lipidreducing", "antithrombotic", "antidm", "smokingstatus",
                                            "alcohol", "SBP", "DBP", "TC", "HDLC", "GLU", "BMI", "WHR", "waist_cm", "CRP", 
                                            "ethnicity", "dm","Sex","First_age","height_cm"))

md.pattern(Scan_covar)
apply(Scan_covar, 2, function(col)sum(is.na(col))/length(col))
Scan_covar$SBP <- haven::zap_labels(Scan_covar$SBP)
Scan_covar$DBP <- haven::zap_labels(Scan_covar$DBP)
Scan_covar$TC <- haven::zap_labels(Scan_covar$TC)
Scan_covar$HDLC <- haven::zap_labels(Scan_covar$HDLC)
Scan_covar$GLU <- haven::zap_labels(Scan_covar$GLU)
Scan_covar$waist_cm <- haven::zap_labels(Scan_covar$waist_cm)
Scan_covar$CRP <- haven::zap_labels(Scan_covar$CRP)
Scan_covar$height_cm <- haven::zap_labels(Scan_covar$height_cm)
dput(names(Scan_covar))


####baseline physical activity----
total_PA_ergo5 <- read_sav('ERGO 5 (RS-I-5, II-3, III-2) - LASA (07-2020).sav')
total_PA_ergojong <- read_sav('ERGO JONG (RS-III-1) - LASA (07-2020).sav')
dput(names(total_PA_ergo5))
total_PA_ergo5_0 <- subset(total_PA_ergo5, select = c("ergoid", "e5_total_MET", "valid"))
total_PA_ergo5_0 <- total_PA_ergo5_0[total_PA_ergo5_0$valid == 1,]
total_PA_ergo5_0 <- total_PA_ergo5_0 %>%
  distinct(ergoid, .keep_all = TRUE)
dput(names(total_PA_ergojong))
total_PA_ergojong_0 <- subset(total_PA_ergojong, select = c("ergoid", "ej_total_MET", "valid"))
total_PA_ergojong_0 <- total_PA_ergojong_0[total_PA_ergojong_0$valid == 1,]
total_PA_ergojong_0 <- total_PA_ergojong_0 %>%
  distinct(ergoid, .keep_all = TRUE)
detail_PA_ergo5 <- read_sav('ERGO 5 (RS-I-5, II-3, III-2) - LASA - modvig PA (07-2020).sav')
detail_PA_ergojong <- read_sav('ERGO JONG (RS-III-1) - LASA - modvig PA (07-2020).sav')
dput(names(detail_PA_ergo5))
detail_PA_ergo5_0 <- subset(detail_PA_ergo5, select = c("ergoid" ,"METh_moderate_e5", "METh_vigorous_e5", "METh_modvig_e5"))
detail_PA_ergo5_0 <- detail_PA_ergo5_0 %>%
  distinct(ergoid, .keep_all = TRUE)
dput(names(detail_PA_ergojong))
detail_PA_ergojong_0 <- subset(detail_PA_ergojong, select = c("ergoid" ,"METh_moderate_ej", "METh_vigorous_ej", "METh_modvig_ej"))
detail_PA_ergojong_0 <- detail_PA_ergojong_0 %>%
  distinct(ergoid, .keep_all = TRUE)
total_PA <- full_join(total_PA_ergo5_0, total_PA_ergojong_0, by = "ergoid")
total_PA <- mutate(total_PA, total_MET = coalesce(ej_total_MET, e5_total_MET))
total_PA <- total_PA[!is.na(total_PA$ergoid),]
summary(total_PA)
total_PA_1 <- total_PA[,c("total_MET","ergoid")]
detail_PA <- full_join(detail_PA_ergo5_0, detail_PA_ergojong_0, by = "ergoid")
detail_PA <- mutate(detail_PA, METh_moderate = coalesce(METh_moderate_ej, METh_moderate_e5)) %>%
  mutate(detail_PA, METh_vigorous = coalesce(METh_vigorous_ej, METh_vigorous_e5)) %>%
  mutate(detail_PA, METh_modvig = coalesce(METh_modvig_ej, METh_modvig_e5))
detail_PA <- detail_PA[!is.na(detail_PA$ergoid),]
summary(detail_PA)
detail_PA_1 <- detail_PA[,c("METh_modvig","METh_vigorous","METh_moderate","ergoid")]
data_PA <- left_join(total_PA_1, detail_PA_1, by = "ergoid")
summary(data_PA)




####visit time of Physical activity----
detail_PA <- mutate(detail_PA, fo=case_when(
  !is.na(METh_moderate_ej) ~ "ej",
  is.na(METh_moderate_ej) ~ "e5"
))
detail_PA_ej <- subset(detail_PA, fo == "ej")
detail_PA_ej <- subset(detail_PA_ej, select = c("ergoid", "fo"))
detail_PA_e5 <- subset(detail_PA, fo == "e5")
detail_PA_e5 <- subset(detail_PA_e5, select = c("ergoid", "fo"))

ej_date <- read_sav('ej_(1)_RESPONS_(04-apr-2016)_excerpt.sav')
ej_date_1 <- subset(ej_date, select = c("ergoid", "ej_3493"))
ej_date_1 <- rename(ej_date_1, visit_date = ej_3493)

e5_date <- read_sav('e5_(5)_RESPONS_(03-mar-2021)_excerpt.sav')
e5_date_1 <- subset(e5_date, select = c("ergoid", "e5_3493"))
e5_date_1 <- rename(e5_date_1, visit_date = e5_3493)

e6_date <- read_sav('e6_(6)_RESPONS_(10-feb-2017)_EXCERPT.sav')
e6_date_1 <- subset(e6_date, select = c("ergoid", "e6_3493"))
e6_date_1 <- rename(e6_date_1, visit_date = e6_3493)

detail_PA_ej <- left_join(detail_PA_ej, ej_date_1)
detail_PA_e5 <- left_join(detail_PA_e5, e5_date_1)

PA_visit_date <- rbind(detail_PA_ej, detail_PA_e5)

Scan_ID <- ID_0
dput(names(Scan_ID))
Scan_ID <- subset(Scan_ID, select = c("ergoid", "ergocar_date"))

visit_time_interval <- left_join(Scan_ID, PA_visit_date, by = "ergoid")

visit_time_interval$interval <- time_length(difftime(visit_time_interval$visit_date, visit_time_interval$ergocar_date), "years")
summary(visit_time_interval$interval)


####baseline energey----
diet_data <- read_sav('ERGO5_dietary_lipids.sav')
dput(names(diet_data))
energy <- subset(diet_data, select = c("ergoid", "kcal_item_sum"))
energy <- rename(energy, energy = kcal_item_sum)


####employement----
ej_job <- read_sav('ej_intvw_SES_(23-mar-2011).sav')
dput(names(ej_job))
ej_job_1 <- subset(ej_job, select = c("ergoid","ej_yses2"))
ej_job_1 <- rename(ej_job_1, job=ej_yses2)
#0-in work，1~6not in work，7-9NA
ej_job_1 <- mutate(ej_job_1, job=case_when(
  job == 0 ~ 1,
  job == 7|job == 9 ~ NA_real_,
  TRUE ~ 0
))

e5_job <- read_sav('e5_intvw_SES_(07-feb-2014).sav')
dput(names(e5_job))
e5_job_1 <- subset(e5_job, select = c("ergoid","e5_EISES2"))
e5_job_1 <- rename(e5_job_1, job=e5_EISES2)
#0-in work，1~6not in work，7-9NA
e5_job_1 <- mutate(e5_job_1, job=case_when(
  job == 0 ~ 1,
  job == 7|job == 9 ~ NA_real_,
  TRUE ~ 0
))

Scan_ID <- ID_0
names(Scan_ID)
Scan_ID <- distinct(Scan_ID, First_ID, .keep_all= TRUE)
Scan_ID$followup <- NA
Scan_ID$followup[Scan_ID$ergocar_date<ymd("2009-01-20")] <- "ej"
Scan_ID$followup[Scan_ID$ergocar_date>ymd("2009-01-19")] <- "e5"

job_ID <- Scan_ID
job_ID <- subset(job_ID, select = c("ergoid", "followup"))

job_ID_ej <- subset(job_ID, followup == "ej")
job_ID_ej <- left_join(job_ID_ej, ej_job_1, by = "ergoid")

job_ID_e5 <- subset(job_ID, followup == "e5")
job_ID_e5 <- left_join(job_ID_e5, e5_job_1, by = "ergoid")

job_ID_final <- rbind(job_ID_ej, job_ID_e5)
job_ID_complete <- subset(job_ID_final, !is.na(job))
job_ID_NA <- subset(job_ID_final, is.na(job))
job_ID_NA <- subset(job_ID_NA, select = c("ergoid", "followup"))
job_ID_NA <- left_join(job_ID_NA, e5_job_1, by = "ergoid")
job_ID_final <- rbind(job_ID_NA, job_ID_complete)
table(job_ID_final$job)
job_ID_final <- subset(job_ID_final, select = c("ergoid", "job"))


####family disease history----
e1_his <- read_sav('e1_intvw_FAMHIST_(24-jan-2012).sav')
dput(names(e1_his))
e1_his_1 <- subset(e1_his, select = c("ergoid","e1_cvafd","e1_mifd"))
e1_his_1 <- mutate(e1_his_1, CVD_his=case_when(
  e1_cvafd > 0 | e1_mifd > 0 ~ 1,
  TRUE~0
))
e1_his_2 <- subset(e1_his_1, select = c("ergoid","CVD_his"))

ep_his <- read_sav('ep_(1)_FAMHST_1_(22-feb-2010).sav')
dput(names(ep_his))
ep_his_1 <- subset(ep_his, select = c("ergoid","ep_6604","ep_6620"))
ep_his_1$ep_6604[ep_his_1$ep_6604==8|ep_his_1$ep_6604==9] <- NA
ep_his_1$ep_6620[ep_his_1$ep_6620==8|ep_his_1$ep_6620==9] <- NA
ep_his_1 <- mutate(ep_his_1, CVD_his=case_when(
  ep_6604 == 1 | ep_6620 == 1 ~ 1,
  ep_6604 == 0 & ep_6620 == 0 ~ 0,
  TRUE~NA_real_
))
ep_his_2 <- subset(ep_his_1, select = c("ergoid","CVD_his"))

ej_his <- read_sav('ej_intvw_FAMHIST_(11-apr-2011).sav')
dput(names(ej_his))
ej_his_1 <- subset(ej_his, select = c("ergoid","ej_yihafa1","ej_yibefa1"))
ej_his_1$ej_yihafa1[ej_his_1$ej_yihafa1==8|ej_his_1$ej_yihafa1==9|ej_his_1$ej_yihafa1==7] <- NA
ej_his_1$ej_yibefa1[ej_his_1$ej_yibefa1==8|ej_his_1$ej_yibefa1==9|ej_his_1$ej_yibefa1==7] <- NA
ej_his_1 <- mutate(ej_his_1, CVD_his=case_when(
  ej_yibefa1 == 1 | ej_yibefa1 == 1 ~ 1,
  ej_yibefa1 == 0 & ej_yibefa1 == 0 ~ 0,
  TRUE~NA_real_
))
ej_his_2 <- subset(ej_his_1, select = c("ergoid","CVD_his"))

fhis <- rbind(e1_his_2,ep_his_2,ej_his_2)
fhis$CVD_his <- factor(fhis$CVD_his)

####inco_edu----
Income_RS1 <- read_sav('e1_intvw_SES_(24-jan-2012).sav')
dput(names(Income_RS1))
Income_RS1_1 <- Income_RS1[, c("ergoid", "e1_ai8eqinc")]
Income_RS1_1 <- rename(Income_RS1_1, income_1 = e1_ai8eqinc)
Income_RS1_1$income_1[Income_RS1_1$income_1 == 99.99] <- NA
Income_RS1_1 <- mutate(Income_RS1_1,
                       income = case_when(
                         income_1 < 1 ~ 1,
                         income_1 >= 1 & income_1 < 1.4 ~ 2,
                         income_1 >= 1.4 & income_1 < 1.7 ~ 3,
                         income_1 >= 1.7 & income_1 < 1.9 ~ 4,
                         income_1 >= 1.9 & income_1 < 2.1 ~ 5,
                         income_1 >= 2.1 & income_1 < 2.4 ~ 6,
                         income_1 >= 2.4 & income_1 < 2.7 ~ 7,
                         income_1 >= 2.7 & income_1 < 3.0 ~ 8,
                         income_1 >= 3.0 & income_1 < 3.5 ~ 9,
                         income_1 >= 3.5 & income_1 < 4.2 ~ 10,
                         income_1 >= 4.2 & income_1 < 5.0 ~ 11,
                         income_1 >= 5.0 & income_1 < 5.8 ~ 12,
                         income_1 >= 5.8 ~ 13,
                         TRUE ~ NA_real_
                       ))
Income_RS1_1 <- Income_RS1_1[, c("ergoid", "income")]

Income_RS2 <- read_sav('ep_intvw_SES_(27-oct-2011).sav')
dput(names(Income_RS2))
Income_RS2_1 <- Income_RS2[, c("ergoid", "ep_ses8")]
Income_RS2_1 <- rename(Income_RS2_1, income = ep_ses8)
Income_RS2_1$income[Income_RS2_1$income == 99 | Income_RS2_1$income == 77] <- NA
Income_RS2_1$income <- as.numeric(Income_RS2_1$income)

Income_RS3 <- read_sav('ej_intvw_SES_(23-mar-2011).sav')
dput(names(Income_RS3))
Income_RS3_1 <- Income_RS3[, c("ergoid", "ej_yses8")]
Income_RS3_1 <- rename(Income_RS3_1, income = ej_yses8)
Income_RS3_1$income[Income_RS3_1$income == 99 | Income_RS3_1$income == 77] <- NA
Income_RS3_1$income <- as.numeric(Income_RS3_1$income)

income <- rbind(Income_RS1_1, Income_RS2_1) %>%
  rbind(., Income_RS3_1)

income <- mutate(income,
                 income_cate = case_when(
                   income <= 5 ~ 1,
                   income > 5 & income <= 9 ~ 2,
                   income > 9 ~ 3,
                   TRUE ~ NA_real_
                 ))
income$income_cate <- factor(income$income_cate, labels = c("<=2100fl./month", "2100~3500fl./month", ">3500fl./month"))
summary(income$income_cate)

####education----
education <- read_sav('Education RS-I-II-III (UNESCO class)_(12-MAR-2015).sav')
dput(names(education))
education_1 <- education[, c("ergoid", "rs_cohort","ses_UNESCO_recoded")]
education_1 <- rename(education_1, edu = ses_UNESCO_recoded)
education_1$edu <- factor(education_1$edu, labels = c("primary education","lower/intermediate general education OR lower vocational education",
                                                      "intermediate vocational education OR higher general education",
                                                      "higher vocational education OR university"))


income_edu <- full_join(income, education_1, by = "ergoid")

