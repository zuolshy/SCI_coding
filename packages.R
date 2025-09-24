#################加载数据和包，整理原始数据----
# install.packages("descriptr")
library(readxl)
library(pubh)
library(psych)
library(plyr)
library(dplyr)
library(haven)#读spss文件
library(lubridate)#输入日期
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
library(cmprsk)#加载competing risk的包
library(egg)
library(patchwork)#组合ggplot的图
library(contsurvplot)#画连续性变量的生存曲线
library(pammtools)#画连续性变量的生存曲线需要加载的包
library(gtable)
library(prodlim)
library(riskRegression)
library(plotly)#画3d图
library(webshot2)#plotly导出的是html格式的图，需要转化为png才能进行增加坐标轴这样的操作
library(png)#读入png
library(grid)#读入png
library(imager)#对PNG图像进行裁割
library(networkD3)
rm(list=ls(all=TRUE))
options(scipen = 999)#设定不使用科学计数法
options(digits = 5)

####基线PA----
# total_PA_ergo5 <- read_sav('V:/HomeDir/069509(L_Zuo)/Lipid_Plaque/data/ERGO 5 (RS-I-5, II-3, III-2) - LASA (07-2020).sav')
# total_PA_ergojong <- read_sav('V:/HomeDir/069509(L_Zuo)/Lipid_Plaque/data/ERGO JONG (RS-III-1) - LASA (07-2020).sav')
dput(names(total_PA_ergo5))
total_PA_ergo5_0 <- subset(total_PA_ergo5, select = c("ergoid", "e5_total_MET", "valid"))
total_PA_ergo5_0 <- total_PA_ergo5_0[total_PA_ergo5_0$valid == 1,]#原始数据中，有一部分的PA是超出了正常认知范围的，所以有对应的valid的变量
total_PA_ergo5_0 <- total_PA_ergo5_0 %>%
  distinct(ergoid, .keep_all = TRUE)
dput(names(total_PA_ergojong))
total_PA_ergojong_0 <- subset(total_PA_ergojong, select = c("ergoid", "ej_total_MET", "valid"))
total_PA_ergojong_0 <- total_PA_ergojong_0[total_PA_ergojong_0$valid == 1,]
total_PA_ergojong_0 <- total_PA_ergojong_0 %>%
  distinct(ergoid, .keep_all = TRUE)
# detail_PA_ergo5 <- read_sav('V:/HomeDir/069509(L_Zuo)/Lipid_Plaque/data/ERGO 5 (RS-I-5, II-3, III-2) - LASA - modvig PA (07-2020).sav')
# detail_PA_ergojong <- read_sav('V:/HomeDir/069509(L_Zuo)/Lipid_Plaque/data/ERGO JONG (RS-III-1) - LASA - modvig PA (07-2020).sav')
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


######重复测量PA----
#rsID
ID_RS <- subset(data_baseline, select = c("ergoid", "rs_cohort"))
ID_twoPA <- subset(ID_0, select = c("ergoid", "Follow"))
ID_RS_0 <- left_join(ID_RS, ID_twoPA, by = "ergoid")
#1740人中，1有677,2有620,3有443

ID_RS_1 <- subset(ID_RS_0, rs_cohort == 1)
ID_RS_2 <- subset(ID_RS_0, rs_cohort == 2)
ID_RS_3 <- subset(ID_RS_0, rs_cohort == 3)

total_PA_ergo6 <- read_sav('ERGO 6 (RS-I-6, II-4) - LASA (07-2020).sav')
total_PA_ergo6_0 <- subset(total_PA_ergo6, select = c("ergoid", "e6_total_MET", "valid"))
detail_PA_ergo6 <- read_sav('ERGO 6 (RS-I-6, II-4) - LASA - modvig PA (07-2020).sav')
detail_PA_ergo6_0 <- subset(detail_PA_ergo6, select = c("ergoid" ,"METh_moderate_e6", "METh_vigorous_e6", "METh_modvig_e6"))


#基线测量
twoMRI_RS_1 <- left_join(ID_RS_1, total_PA_ergo5_0, by = "ergoid")
twoMRI_RS_2 <- left_join(ID_RS_2, total_PA_ergo5_0, by = "ergoid")
twoMRI_RS_3 <- left_join(ID_RS_3, total_PA_ergojong_0, by = "ergoid")
twoMRI_RS_1 <- subset(twoMRI_RS_1, !is.na(e5_total_MET))
twoMRI_RS_2 <- subset(twoMRI_RS_2, !is.na(e5_total_MET))
twoMRI_RS_3 <- subset(twoMRI_RS_3, !is.na(ej_total_MET))
twoMRI_RS_1#516
twoMRI_RS_2#507
twoMRI_RS_3#327
table(twoMRI_RS_1$Follow)#116
table(twoMRI_RS_2$Follow)#330
table(twoMRI_RS_3$Follow)#191

#两次随访测量
twoMRI_RS_1 <- left_join(ID_RS_1, total_PA_ergo5_0, by = "ergoid") %>% 
  left_join(., total_PA_ergo6_0, by = "ergoid") %>% 
  left_join(., detail_PA_ergo5_0, by = "ergoid") %>% 
  left_join(., detail_PA_ergo6_0, by = "ergoid")
twoMRI_RS_1 <- subset(twoMRI_RS_1, !is.na(e5_total_MET)&!is.na(e6_total_MET))
twoMRI_RS_1#266
twoMRI_RS_1 <- rename(twoMRI_RS_1, total_1 = e5_total_MET, total_2 = e6_total_MET, 
                      moderate_1 = METh_moderate_e5, moderate_2 = METh_moderate_e6,
                      vigorous_1 = METh_vigorous_e5, vigorous_2 = METh_vigorous_e6)
twoMRI_RS_1_final <- subset(twoMRI_RS_1, select = c("ergoid", "total_1", "total_2", "moderate_1", "moderate_2", "vigorous_1", "vigorous_2"))

twoMRI_RS_2 <- left_join(ID_RS_2, total_PA_ergo5_0, by = "ergoid") %>% 
  left_join(., total_PA_ergo6_0, by = "ergoid") %>% 
  left_join(., detail_PA_ergo5_0, by = "ergoid") %>% 
  left_join(., detail_PA_ergo6_0, by = "ergoid")
twoMRI_RS_2 <- subset(twoMRI_RS_2, !is.na(e5_total_MET)&!is.na(e6_total_MET))
twoMRI_RS_2#380
twoMRI_RS_2 <- rename(twoMRI_RS_2, total_1 = e5_total_MET, total_2 = e6_total_MET, 
                      moderate_1 = METh_moderate_e5, moderate_2 = METh_moderate_e6,
                      vigorous_1 = METh_vigorous_e5, vigorous_2 = METh_vigorous_e6)
twoMRI_RS_2_final <- subset(twoMRI_RS_2, select = c("ergoid", "total_1", "total_2", "moderate_1", "moderate_2", "vigorous_1", "vigorous_2"))

twoMRI_RS_3 <- left_join(ID_RS_3, total_PA_ergojong_0, by = "ergoid") %>% 
  left_join(., total_PA_ergo5_0, by = "ergoid") %>% 
  left_join(., detail_PA_ergojong_0, by = "ergoid") %>% 
  left_join(., detail_PA_ergo5_0, by = "ergoid")
twoMRI_RS_3 <- subset(twoMRI_RS_3, !is.na(e5_total_MET)&!is.na(ej_total_MET))
twoMRI_RS_3#233
twoMRI_RS_3 <- rename(twoMRI_RS_3, total_1 = ej_total_MET, total_2 = e5_total_MET, 
                      moderate_1 = METh_vigorous_ej, moderate_2 = METh_moderate_e5,
                      vigorous_1 = METh_modvig_ej, vigorous_2 = METh_vigorous_e5)
twoMRI_RS_3_final <- subset(twoMRI_RS_3, select = c("ergoid", "total_1", "total_2", "moderate_1", "moderate_2", "vigorous_1", "vigorous_2"))

twoPA_data <- rbind(twoMRI_RS_1_final, twoMRI_RS_2_final, twoMRI_RS_3_final)

#基线与随访的PA情况
rg1 <- lm(data=twoPA_data, total_1 ~ total_2)
summary(rg1)

rg1 <- lm(data=twoPA_data, moderate_1 ~ moderate_2)
summary(rg1)

rg1 <- lm(data=twoPA_data, vigorous_1 ~ vigorous_2)
summary(rg1)


####随访时间----
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


####基线energey----
energy <- rename(energy, energy = kcal_item_sum)#energy协变量要从膳食数据中导入


####工作情况----
ej_job <- read_sav('ej_intvw_SES_(23-mar-2011).sav')
dput(names(ej_job))
ej_job_1 <- subset(ej_job, select = c("ergoid","ej_yses2"))
ej_job_1 <- rename(ej_job_1, job=ej_yses2)
#0-在工作，1~6不工作，7-9NA
ej_job_1 <- mutate(ej_job_1, job=case_when(
  job == 0 ~ 1,
  job == 7|job == 9 ~ NA_real_,
  TRUE ~ 0
))

e5_job <- read_sav('e5_intvw_SES_(07-feb-2014).sav')
dput(names(e5_job))
e5_job_1 <- subset(e5_job, select = c("ergoid","e5_EISES2"))
e5_job_1 <- rename(e5_job_1, job=e5_EISES2)
#0-在工作，1~6不工作，7-9NA
e5_job_1 <- mutate(e5_job_1, job=case_when(
  job == 0 ~ 1,
  job == 7|job == 9 ~ NA_real_,
  TRUE ~ 0
))

Scan_ID <- ID_0
names(Scan_ID)
Scan_ID <- distinct(Scan_ID, First_ID, .keep_all= TRUE)#根据基线的ID去排除重复ID
Scan_ID$followup <- NA
Scan_ID$followup[Scan_ID$ergocar_date<ymd("2009-01-20")] <- "ej"#用mutate输入“ej”总是显示有问题mmp
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


####家族史----
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

####Geneactiv----
watch_PA_E5 <- read_sav('GeneActiv data ERGO 5 (RS-II-3, RS-III-2) - means.sav')
dput(names(watch_PA_E5))
watch_PA_E5_1 <- subset(watch_PA_E5, select = c("ergoid", "first_wear_date", "sedentary_time_mean", "light_time_mean", "moderate_time", 
                                                "vigorous_time_mean", "mod_vig_time_mean", 
                                                "awake_sedentary_time_mean", "out_bed_sedentary_time_mean"))
watch_PA_E6 <- read_sav('GeneActiv data ERGO 6 (I-6, II-4) - means.sav')
dput(names(watch_PA_E6))
watch_PA_E6_1 <- subset(watch_PA_E6, select = c("ergoid", "first_wear_date", "sedentary_time_mean", "light_time_mean", "moderate_time", 
                                                "vigorous_time_mean", "mod_vig_time_mean", 
                                                "awake_sedentary_time_mean", "out_bed_sedentary_time_mean"))
watch_PA <- rbind(watch_PA_E6_1, watch_PA_E5_1)
watch_PA_1 <- watch_PA %>%
  distinct(ergoid, .keep_all = TRUE)

watch_PA_2 <- left_join(data_baseline, watch_PA_1, by = "ergoid")
dput(names(watch_PA_2))
watch_PA_2$time <- time_length(difftime(watch_PA_2$Carotid_scan_new, watch_PA_2$first_wear_date), "years")
summary(watch_PA_2)

watch_PA_3 <- left_join(data_total_twoscan, watch_PA_1, by = "ergoid")
dput(names(watch_PA_3))
watch_PA_3$time <- time_length(difftime(watch_PA_3$ergocar_date, watch_PA_3$first_wear_date), "years")
watch_PA_3$time_1 <- time_length(difftime(watch_PA_3$ergocarfup_date, watch_PA_3$first_wear_date), "years")
watch_PA_3 <- subset(watch_PA_3, select = c("ergoid", "sedentary_time_mean", "light_time_mean", 
                                            "moderate_time", "vigorous_time_mean", "mod_vig_time_mean", "awake_sedentary_time_mean", 
                                            "out_bed_sedentary_time_mean"))
watch_PA_3 <- watch_PA_3 %>%
  distinct(ergoid, .keep_all = TRUE)

####两种检测的一致性----
PA_dif <- watch_PA_3
dput(names(PA_dif))
PA_dif <- subset(PA_dif, !is.na(sedentary_time_mean))
PA_dif <- left_join(PA_dif, data_PA, by = "ergoid")
PA_dif <- subset(PA_dif, !is.na(total_MET))
PA_dif$vigorous <- PA_dif$METh_vigorous
PA_dif$vigorous[PA_dif$vigorous == 0] <- NA

PA_dif <- PA_dif %>%
  mutate(vigorous_3 = ntile(vigorous, 2))
PA_dif$vigorous_3[is.na(PA_dif$vigorous_3)] <- 0
PA_dif$vigorous_3 <- PA_dif$vigorous_3 + 1
# PA_dif$vigorous_3 <- factor(PA_dif$vigorous_3)
quantile(PA_dif$vigorous_time_mean, c(0.66,0.83))
PA_dif <- mutate(PA_dif,
  vig_3 = case_when(
    vigorous_time_mean < 16.8 ~ 1,
    vigorous_time_mean >= 16.8 & vigorous_time_mean < 20.957 ~ 2, 
    vigorous_time_mean >= 20.957 ~3
  ))
# PA_dif$vig_3 <- factor(PA_dif$vig_3)

x=cbind(PA_dif$vigorous_3, PA_dif$vig_3)
table(PA_dif$vigorous_3, PA_dif$vig_3)
cohen.kappa(x)

summary(lm(vig_3~vigorous_3,data = PA_dif))#一致性很差
summary(lm(METh_vigorous~vigorous_time_mean,data = PA_dif))
PA_dif %>% 
  group_by(vigorous_3) %>% 
  summarise(mean = mean(vigorous_time_mean))

summary(PA_dif$vigorous_3, PA_dif$vigorous_time_mean)

####随访CVD，使用旧版CHD----
CVD_ID_data1 <- data_secondscan_5
#cvd 结局
#Stroke <- read_sav('V:/HomeDir/069509(L_Zuo)/outcome_CVD/stroke_new/Luoshiyuan_Z_Stroke.sav')
names(Stroke)
Stroke_0 <- subset(Stroke, select = c("ergoid","IC", "prevalent_stroke","Baseline_date","incident_stroke","stroke_date"))
Stroke_1 <- rename(Stroke_0, ic_Stroke = IC, prev_stroke_date = Baseline_date, prev_Stroke = prevalent_stroke, inc_stroke = incident_stroke, inc_stroke_date = stroke_date)
Stroke_1 <- Stroke_1[Stroke_1$ic_Stroke == 1,]
Stroke_2 <- mutate(Stroke_1,
                   Stroke_date = case_when(
                     prev_Stroke == 0 & inc_stroke == 1 ~ inc_stroke_date,
                     prev_Stroke == 1 ~ prev_stroke_date,
                     prev_Stroke == 0 & inc_stroke == 0 ~ as.Date("2019-12-31")))
Stroke_2 <- mutate(Stroke_2,
                   Stroke = case_when(
                     prev_Stroke == 0 & inc_stroke == 0 ~ 0,
                     prev_Stroke == 1 | inc_stroke == 1 ~ 1,
                     TRUE  ~ NA_real_))
#CHD <- read_sav('V:/HomeDir/069509(L_Zuo)/outcome_CVD/CHD_new/CHD - ZW - old prevalent merged.sav')
# names(CHD)
CHD_old <- read_sav('G:/BaiduSyncdisk/Rcoding/3PA_plaque/Coronary heart disease with dates - MGJ Leening- updated - clean(1).sav')
names(CHD_old)
CHD_0 <- subset(CHD_old, select = c("ergoid","ic_ok","fp_startdate","prev_CHD","inc_CHD","enddat_CHD"))
CHD_1 <- rename(CHD_0, ic_CHD = ic_ok, prev_CHD_date = fp_startdate, prev_CHD = prev_CHD, inc_CHD = inc_CHD, inc_CHD_date = enddat_CHD)
CHD_1 <- CHD_1[CHD_1$ic_CHD == 1,]
CHD_2 <- mutate(CHD_1,
                CHD_date = case_when(
                  prev_CHD == 0 & inc_CHD == 1 ~ inc_CHD_date,
                  prev_CHD == 1 ~ prev_CHD_date,
                  prev_CHD == 0 & inc_CHD == 0 ~ inc_CHD_date))
CHD_2 <- mutate(CHD_2,
                CHD = case_when(
                  prev_CHD == 0 & inc_CHD == 0 ~ 0,
                  prev_CHD == 1 | inc_CHD == 1 ~ 1,
                  TRUE  ~ NA_real_))
CVD_ID_data1 <- left_join(CVD_ID_data1, Stroke_2[c("ergoid", "Stroke_date","Stroke")], by = "ergoid")
CVD_ID_data1 <- left_join(CVD_ID_data1, CHD_2[c("ergoid", "CHD_date", "CHD")], by = "ergoid")
CVD_ID_data1 <- mutate(CVD_ID_data1,
                       CHD_status = case_when(
                         CHD == 1 & CHD_date <= ergocar_date ~ 1,
                         CHD == 1 & CHD_date > ergocar_date ~ 2,
                         is.na(CHD_date) == TRUE ~  NA,
                         CHD == 0 ~ 0))
CVD_ID_data1 <- mutate(CVD_ID_data1,
                       Stroke_status = case_when(
                         Stroke == 1 & Stroke_date <= ergocar_date ~ 1,
                         Stroke == 1 & Stroke_date > ergocar_date ~ 2,
                         is.na(Stroke_date) == TRUE ~  NA,
                         Stroke == 0 ~ 0))
#CVD数据
CVD_ID_data1 <- mutate(CVD_ID_data1,
                       CVD_status = case_when(
                         CHD_status == 1 | Stroke_status == 1 ~ 1,
                         CHD_status == 2 & Stroke_status == 0 ~ 2,
                         CHD_status == 0 & Stroke_status == 2 ~ 2,
                         CHD_status == 2 & Stroke_status == 2 ~ 2,
                         is.na(CHD_status) == TRUE  ~ NA_real_,
                         TRUE ~ 0))

CVD_ID_data1 <- distinct(CVD_ID_data1, ergoid, .keep_all = TRUE)
table(CVD_ID_data1$CVD_status)
table(CVD_ID_data1$Stroke_status)
table(CVD_ID_data1$CHD_status)
CVD_ID_1 <- CVD_ID_data1[,c("ergoid","CVD_status","Stroke_status","CHD_status")]
summary(CVD_ID_1)

















####随访CVD情况方法1----
CVD_ID_data1 <- data_secondscan_5
# CVD_ID_data1 <- CVD_ID_data1[!is.na(CVD_ID_data1$status_1),]
# CVD_ID_data1 <- CVD_ID_data1[!is.na(CVD_ID_data1$status_2),]
#cvd 结局
#Stroke <- read_sav('V:/HomeDir/069509(L_Zuo)/outcome_CVD/stroke_new/Luoshiyuan_Z_Stroke.sav')
names(Stroke)
Stroke_0 <- subset(Stroke, select = c("ergoid","IC", "prevalent_stroke","Baseline_date","incident_stroke","stroke_date"))
Stroke_1 <- rename(Stroke_0, ic_Stroke = IC, prev_stroke_date = Baseline_date, prev_Stroke = prevalent_stroke, inc_stroke = incident_stroke, inc_stroke_date = stroke_date)
Stroke_1 <- Stroke_1[Stroke_1$ic_Stroke == 1,]
Stroke_2 <- mutate(Stroke_1,
                   Stroke_date = case_when(
                     prev_Stroke == 0 & inc_stroke == 1 ~ inc_stroke_date,
                     prev_Stroke == 1 ~ prev_stroke_date,
                     prev_Stroke == 0 & inc_stroke == 0 ~ as.Date("2019-12-31")))
Stroke_2 <- mutate(Stroke_2,
                   Stroke = case_when(
                     prev_Stroke == 0 & inc_stroke == 0 ~ 0,
                     prev_Stroke == 1 | inc_stroke == 1 ~ 1,
                     TRUE  ~ NA_real_))
#CHD <- read_sav('V:/HomeDir/069509(L_Zuo)/outcome_CVD/CHD_new/CHD - ZW - old prevalent merged.sav')
# names(CHD)
CHD_0 <- subset(CHD, select = c("ergoid","ic_ok","startdat","Prevalent_CHD","Incident_CHD","end_of_FU_CHD"))
CHD_1 <- rename(CHD_0, ic_CHD = ic_ok, prev_CHD_date = startdat, prev_CHD = Prevalent_CHD, inc_CHD = Incident_CHD, inc_CHD_date = end_of_FU_CHD)
CHD_1 <- CHD_1[CHD_1$ic_CHD == 1,]
CHD_2 <- mutate(CHD_1,
                CHD_date = case_when(
                  prev_CHD == 0 & inc_CHD == 1 ~ inc_CHD_date,
                  prev_CHD == 1 ~ prev_CHD_date,
                  prev_CHD == 0 & inc_CHD == 0 ~ inc_CHD_date))
CHD_2 <- mutate(CHD_2,
                CHD = case_when(
                  prev_CHD == 0 & inc_CHD == 0 ~ 0,
                  prev_CHD == 1 | inc_CHD == 1 ~ 1,
                  TRUE  ~ NA_real_))
CVD_ID_data1 <- left_join(CVD_ID_data1, Stroke_2[c("ergoid", "Stroke_date","Stroke")], by = "ergoid")
CVD_ID_data1 <- left_join(CVD_ID_data1, CHD_2[c("ergoid", "CHD_date", "CHD")], by = "ergoid")
CVD_ID_data1 <- mutate(CVD_ID_data1,
               CHD_status = case_when(
                 CHD == 1 & CHD_date <= ergocar_date ~ 1,
                 CHD == 1 & CHD_date > ergocar_date ~ 2,
                 is.na(CHD_date) == TRUE ~  NA,
                 CHD == 0 ~ 0))
CVD_ID_data1 <- mutate(CVD_ID_data1,
               Stroke_status = case_when(
                 Stroke == 1 & Stroke_date <= ergocar_date ~ 1,
                 Stroke == 1 & Stroke_date > ergocar_date ~ 2,
                 is.na(Stroke_date) == TRUE ~  NA,
                 Stroke == 0 ~ 0))
#CVD数据
CVD_ID_data1 <- mutate(CVD_ID_data1,
               CVD_status = case_when(
                 CHD_status == 1 | Stroke_status == 1 ~ 1,
                 CHD_status == 2 & Stroke_status == 0 ~ 2,
                 CHD_status == 0 & Stroke_status == 2 ~ 2,
                 CHD_status == 2 & Stroke_status == 2 ~ 2,
                 is.na(CHD_status) == TRUE  ~ NA_real_,
                 TRUE ~ 0))

CVD_ID_data1 <- distinct(CVD_ID_data1, ergoid, .keep_all = TRUE)
table(CVD_ID_data1$CVD_status)
table(CVD_ID_data1$Stroke_status)
table(CVD_ID_data1$CHD_status)
CVD_ID_1 <- CVD_ID_data1[,c("ergoid","CVD_status","Stroke_status","CHD_status")]
summary(CVD_ID_1)


####随访CVD情况方法2----
CVD_ID_data2 <- data_baseline
#######cvd 结局
#Stroke <- read_sav('V:/HomeDir/069509(L_Zuo)/outcome_CVD/stroke_new/Luoshiyuan_Z_Stroke.sav')
# names(Stroke)
Stroke_0 <- subset(Stroke, select = c("ergoid","IC", "prevalent_stroke","Baseline_date","incident_stroke","stroke_date"))
Stroke_1 <- rename(Stroke_0, ic_Stroke = IC, prev_stroke_date = Baseline_date, prev_Stroke = prevalent_stroke, inc_stroke = incident_stroke, inc_stroke_date = stroke_date)
Stroke_1 <- Stroke_1[Stroke_1$ic_Stroke == 1,]

#CHD <- read_sav('V:/HomeDir/069509(L_Zuo)/outcome_CVD/CHD_new/CHD - ZW - old prevalent merged.sav')
# names(CHD)
CHD_0 <- subset(CHD, select = c("ergoid","ic_ok","startdat","Prevalent_CHD","Incident_CHD","end_of_FU_CHD"))
CHD_1 <- rename(CHD_0, ic_CHD = ic_ok, prev_CHD_date = startdat, prev_CHD = Prevalent_CHD, inc_CHD = Incident_CHD, inc_CHD_date = end_of_FU_CHD)
CHD_1 <- CHD_1[CHD_1$ic_CHD == 1,]

CVD_ID_data2 <- left_join(CVD_ID_data2, Stroke_1[c("ergoid", "prev_Stroke", "inc_stroke", "inc_stroke_date")], by = "ergoid")
CVD_ID_data2 <- left_join(CVD_ID_data2, CHD_1[c("ergoid", "prev_CHD", "inc_CHD", "inc_CHD_date")], by = "ergoid")
# CVD_ID_data2 <- CVD_ID_data2[!is.na(CVD_ID_data2$prev_Stroke),]
CVD_ID_data2$inc_stroke[is.na(CVD_ID_data2$inc_stroke)] <- 0
# CVD_ID_data2 <- CVD_ID_data2[!is.na(CVD_ID_data2$prev_CHD),]
CVD_ID_data2$inc_CHD[CVD_ID_data2$inc_CHD==8] <- 0

#并入MRI随访时间
incident_CVD_sup <- data_secondscan_5
incident_CVD_sup <- subset(incident_CVD_sup, select = c("ergoid", "ergocar_date", "ergocarfup_date", "First_age", "Sex"))
incident_CVD_sup <- incident_CVD_sup %>%
  distinct(ergoid, .keep_all = TRUE)

CVD_ID_data2 <- left_join(CVD_ID_data2, incident_CVD_sup, by = "ergoid")
names(CVD_ID_data2)
CVD_ID_data2$prev_CHD[(CVD_ID_data2$inc_CHD_date <= CVD_ID_data2$ergocar_date)&(CVD_ID_data2$inc_CHD == 1)] <- 1#代码在这里进行了修改之后，不同的cvdid整理的结果就保持了一致
CVD_ID_data2$inc_CHD[(CVD_ID_data2$inc_CHD_date <= CVD_ID_data2$ergocar_date)&(CVD_ID_data2$inc_CHD == 1)] <- 0#代码在这里进行了修改之后，不同的cvdid整理的结果就保持了一致
# CVD_ID_data2 <- CVD_ID_data2[!is.na(CVD_ID_data2$ergoid),]

CVD_ID_data2$inc_stroke_date[is.na(CVD_ID_data2$inc_stroke_date)] <- as.Date("2019-12-31")
CVD_ID_data2 <- mutate(CVD_ID_data2,
                       Stroke_status = case_when(
                         prev_Stroke == 1 & inc_stroke == 0 ~ 1,
                         prev_Stroke == 0 & inc_stroke == 1 ~ 2,
                         prev_Stroke == 0 & inc_stroke == 0 ~ 0,
                         TRUE ~ NA_real_))

CVD_ID_data2 <- mutate(CVD_ID_data2,
                       CHD_status = case_when(
                         prev_CHD == 1 & inc_CHD == 0 ~ 1,
                         prev_CHD == 0 & inc_CHD == 1 ~ 2,
                         prev_CHD == 0 & inc_CHD == 0 ~ 0,
                         TRUE ~ NA_real_))

CVD_ID_data2 <- mutate(CVD_ID_data2,
                       CVD_status = case_when(
                         CHD_status == 1 | Stroke_status == 1 ~ 1,
                         CHD_status == 2 & Stroke_status == 0 ~ 2,
                         CHD_status == 0 & Stroke_status == 2 ~ 2,
                         CHD_status == 2 & Stroke_status == 2 ~ 2,
                         is.na(CHD_status) == TRUE  ~ NA_real_,
                         TRUE ~ 0))
CVD_ID_2 <- CVD_ID_data2[,c("ergoid","CVD_status","Stroke_status","CHD_status")]
summary(CVD_ID_2)

####mortality----
death <- read_sav('G:/BaiduSyncdisk/Rcoding/3PA_plaque/fp_VitalStatus_2022-42.sav')
death_0 <- subset(death, select = c("ergoid", "fp_vitalstatus", "fp_censordate"))
death_1 <- rename(death_0, death = fp_vitalstatus, death_date = fp_censordate)
death_1$death[is.na(death_1$death)] <- 0


####找出两种方法整理出来的cvd的数据的不同之处----
colnames(CVD_ID_2)[2:4] <- paste(colnames(CVD_ID_2)[2:4],"1",sep="_")
CVD_ID_data3 <- left_join(CVD_ID_data1, CVD_ID_2, by = "ergoid")
summary(CVD_ID_data3)
CVD_ID_data3 <- mutate(CVD_ID_data3,
                       Stroke_status_dif = case_when(
                         Stroke_status == NA | Stroke_status_1 == NA ~ NA_real_,
                         TRUE ~ Stroke_status - Stroke_status_1))
CVD_ID_data3 <- mutate(CVD_ID_data3,
                       CHD_status_dif = case_when(
                         CHD_status == NA | CHD_status_1 == NA ~ NA_real_,
                         TRUE ~ CHD_status - CHD_status_1))
CVD_ID_3 <- CVD_ID_data3[,c("ergoid", "ergocar_date",  "Stroke_date", "Stroke", "CHD_date", "CHD", 
                            "CHD_status", "Stroke_status", "CVD_status", "CVD_status_1", 
                            "Stroke_status_1", "CHD_status_1", "Stroke_status_dif", "CHD_status_dif")]
CVD_ID_3$index <- 0
CVD_ID_3$index[(!(CVD_ID_3$Stroke_status_dif == 0)|!(CVD_ID_3$CHD_status_dif == 0)|(is.na(CVD_ID_3$Stroke_status_dif)|(is.na(CVD_ID_3$CHD_status_dif))))] <- 1
CVD_ID_dif <- CVD_ID_3[CVD_ID_3$index == 1,]
#！！！！！非常有意思的结果，两种不同的方法，最后得出来的结果有差别，id6461和id1531这两个人的CHD的随访截止日期在第一次mri之前，并且两个人都没有发病
#所以在第二种cvd的整理方案中，被判别错误，因为第二种判别方法单纯只用到了截止日期，并没有用到既往发病判别变量


table(CVD_ID_data2$CVD_status)
table(CVD_ID_data1$CVD_status)

table(CVD_ID_data2$Stroke_status)
table(CVD_ID_data1$Stroke_status)

table(CVD_ID_data2$CHD_status)
table(CVD_ID_data1$CHD_status)



####验证CVD发病----
CVD_verification_data <- data_baseline
#######cvd 结局
#Stroke <- read_sav('V:/HomeDir/069509(L_Zuo)/outcome_CVD/stroke_new/Luoshiyuan_Z_Stroke.sav')
# names(Stroke)
Stroke_0 <- subset(Stroke, select = c("ergoid","IC", "prevalent_stroke","Baseline_date","incident_stroke","stroke_date"))
Stroke_1 <- rename(Stroke_0, ic_Stroke = IC, prev_stroke_date = Baseline_date, prev_Stroke = prevalent_stroke, inc_stroke = incident_stroke, inc_stroke_date = stroke_date)
Stroke_1 <- Stroke_1[Stroke_1$ic_Stroke == 1,]

#CHD <- read_sav('V:/HomeDir/069509(L_Zuo)/outcome_CVD/CHD_new/CHD - ZW - old prevalent merged.sav')
# names(CHD)
CHD_0 <- subset(CHD, select = c("ergoid","ic_ok","startdat","Prevalent_CHD","Incident_CHD","end_of_FU_CHD"))
CHD_1 <- rename(CHD_0, ic_CHD = ic_ok, prev_CHD_date = startdat, prev_CHD = Prevalent_CHD, inc_CHD = Incident_CHD, inc_CHD_date = end_of_FU_CHD)
CHD_1 <- CHD_1[CHD_1$ic_CHD == 1,]

CVD_verification_data <- left_join(CVD_verification_data, Stroke_1[c("ergoid", "prev_Stroke", "inc_stroke", "inc_stroke_date")], by = "ergoid")
CVD_verification_data <- left_join(CVD_verification_data, CHD_1[c("ergoid", "prev_CHD", "inc_CHD", "inc_CHD_date")], by = "ergoid")
# CVD_verification_data <- CVD_verification_data[!is.na(CVD_verification_data$prev_Stroke),]
CVD_verification_data$inc_stroke[is.na(CVD_verification_data$inc_stroke)] <- 0
# CVD_verification_data <- CVD_verification_data[!is.na(CVD_verification_data$prev_CHD),]
CVD_verification_data$inc_CHD[CVD_verification_data$inc_CHD==8] <- 0

#并入MRI随访时间
incident_CVD_sup <- data_secondscan_5
incident_CVD_sup <- subset(incident_CVD_sup, select = c("ergoid", "ergocar_date", "ergocarfup_date", "First_age", "Sex"))
incident_CVD_sup <- incident_CVD_sup %>%
  distinct(ergoid, .keep_all = TRUE)

CVD_verification_data <- left_join(CVD_verification_data, incident_CVD_sup, by = "ergoid")
names(CVD_verification_data)
CVD_verification_data$prev_CHD[(CVD_verification_data$inc_CHD_date <= CVD_verification_data$ergocar_date)&(CVD_verification_data$inc_CHD == 1)] <- 1#代码在这里进行了修改之后，不同的cvdid整理的结果就保持了一致
CVD_verification_data$inc_CHD[(CVD_verification_data$inc_CHD_date <= CVD_verification_data$ergocar_date)&(CVD_verification_data$inc_CHD == 1)] <- 0#代码在这里进行了修改之后，不同的cvdid整理的结果就保持了一致
# CVD_verification_data <- CVD_verification_data[!is.na(CVD_verification_data$ergoid),]

CVD_verification_data$inc_stroke_date[is.na(CVD_verification_data$inc_stroke_date)] <- as.Date("2019-12-31")

#更改截止时间2015年1月1日
CVD_verification_data$inc_stroke[(CVD_verification_data$inc_stroke_date > as.Date("2015-01-01")) & (CVD_verification_data$inc_stroke == 1)] <- 0
CVD_verification_data$inc_CHD[(CVD_verification_data$inc_CHD_date > as.Date("2015-01-01")) & (CVD_verification_data$inc_CHD == 1)] <- 0

CVD_verification_data <- mutate(CVD_verification_data,
                       Stroke_status = case_when(
                         prev_Stroke == 1 & inc_stroke == 0 ~ 1,
                         prev_Stroke == 0 & inc_stroke == 1 ~ 2,
                         prev_Stroke == 0 & inc_stroke == 0 ~ 0,
                         TRUE ~ NA_real_))

CVD_verification_data <- mutate(CVD_verification_data,
                       CHD_status = case_when(
                         prev_CHD == 1 & inc_CHD == 0 ~ 1,
                         prev_CHD == 0 & inc_CHD == 1 ~ 2,
                         prev_CHD == 0 & inc_CHD == 0 ~ 0,
                         TRUE ~ NA_real_))

CVD_verification_data <- mutate(CVD_verification_data,
                       CVD_status = case_when(
                         CHD_status == 1 | Stroke_status == 1 ~ 1,
                         CHD_status == 2 & Stroke_status == 0 ~ 2,
                         CHD_status == 0 & Stroke_status == 2 ~ 2,
                         CHD_status == 2 & Stroke_status == 2 ~ 2,
                         is.na(CHD_status) == TRUE  ~ NA_real_,
                         TRUE ~ 0))
table(CVD_verification_data$CVD_status)
CVD_verification_data_1 <- subset(CVD_verification_data, CVD_status == 0 | CVD_status == 2)
table(CVD_verification_data_1$Stroke_status)
table(CVD_verification_data_1$CHD_status)
table(CVD_verification_data_1$CVD_status)
dput(names(CVD_verification_data_1))
CVD_verification_ID <- subset(CVD_verification_data_1, select = c("ergoid", "CHD_status", "inc_CHD_date"))#,"inc_stroke_date", "Stroke_status"
# CVD_verification_ID <- subset(CVD_verification_data, select = c("ergoid", "CHD_status", "inc_CHD_date"))#,"inc_stroke_date", "Stroke_status"


#JACC文章的原始数据
CVD_JACC <- read_sav('G:/BaiduSyncdisk/Rcoding/3PA_plaque/MRI_Plaque_strokechd_v8_with plaque.sav')
dput(names(CVD_JACC))
CVD_JACC_0 <- subset(CVD_JACC, select = c("ergoid", "Carotid_scan_new", "inc_CHD", "enddat_CHD", "inc_stroke2015", "enddate_obs2015"))
table(CVD_JACC_0$inc_CHD)
#CHD的原始数据
CHD_2 <- rename(CHD_1, inc_CHD_0 = inc_CHD)
CHD_2$inc_CHD_0[(CHD_2$inc_CHD_date > as.Date("2015-01-01")) & (CHD_2$inc_CHD_0 == 1)] <- 0
#Storkd的原始数据
Stroke_2 <- rename(Stroke_1, inc_stroke_0 = inc_stroke)
Stroke_2$inc_stroke_0[(Stroke_2$inc_stroke_date > as.Date("2015-01-01")) & (Stroke_2$inc_stroke_0 == 1)] <- 0

CVD_dif <- left_join(CVD_JACC_0,CHD_2,by="ergoid")
CVD_dif <- left_join(CVD_dif,Stroke_2,by="ergoid")
table(CVD_dif$inc_CHD, CVD_dif$inc_CHD_0)
table(CVD_dif$inc_stroke2015, CVD_dif$inc_stroke_0)

CVD_dif$CHD_dif <- CVD_dif$inc_CHD_0 - CVD_dif$inc_CHD
CVD_dif_1 <- subset(CVD_dif, !CHD_dif == 0)

#JACC的原始数据和我整理的数据里面的数据虽然都是1349，但是人是有区别的
#提取1740人的ID
ID_total <- data_baseline[, "ergoid"]

#合并JACC文章的数据
CVD_dif_2 <- left_join(ID_total, CVD_JACC_0, by = "ergoid")

#合并最新的数据
CVD_dif_3 <- left_join(CVD_dif_2, CVD_verification_data, by = "ergoid")

#先看MRI日期是否相同
CVD_verification_data$date_check <- ifelse(CVD_verification_data$ergocar_date == CVD_verification_data$Carotid_scan_new, 0, 1)
#结果发现日期确实是一样的，没有错

#只留下必要的变量，让data看起来更简单
CVD_dif_4 <- subset(CVD_dif_3, select = c("ergoid", "inc_CHD.x", "enddat_CHD", "inc_stroke2015", "enddate_obs2015",
                                          "prev_Stroke", "inc_stroke", "inc_stroke_date", "prev_CHD", "inc_CHD.y", 
                                          "inc_CHD_date", "ergocar_date", "ergocarfup_date", "Stroke_status", "CHD_status", "CVD_status"))
#根据JACC数据生成基线有CVD的指示变量
CVD_dif_4$CVD_status_JACC <- 0
CVD_dif_4$CVD_status_JACC[is.na(CVD_dif_4$inc_CHD.x)] <- 1
#根据自己的数据生成CVD的指示变量，主要是将发病的数据给并到基线无CVD的分类里面去
CVD_dif_4$CVD_status_Luo <- CVD_dif_4$CVD_status
CVD_dif_4$CVD_status_Luo[CVD_dif_4$CVD_status_Luo == 2] <- 0
#做表格比较
table_CVD_check <- table(CVD_dif_4$CVD_status_JACC, CVD_dif_4$CVD_status_Luo)
names(dimnames(table_CVD_check)) <- c("CVD_status_JACC", "CVD_status_Luo")
table_CVD_check
#提取出基线CVD判定不同的人
CVD_dif_4$date_check_baseline_CVD <- ifelse(CVD_dif_4$CVD_status_JACC == CVD_dif_4$CVD_status_Luo, 0, 1)
data_check_CVD <- subset(CVD_dif_4, date_check_baseline_CVD == 1)
dput(names(data_check_CVD))
data_check_CVD_ID <- data_check_CVD[,c("ergoid","CVD_status_JACC","CVD_status_Luo")]



#CHD数据对比
#2015年版chd数据
CHD_2015 <- read_sav('G:/BaiduSyncdisk/Rcoding/3PA_plaque/Coronary heart disease with dates - MGJ Leening- updated - clean(1).sav')
CHD_2015_0 <- subset(CHD_2015, select = c("ergoid", "fp_startdate", "prev_CHD", "inc_CHD", "enddat_CHD"))
CHD_2015_1 <- rename(CHD_2015_0, prev_CHD_date_2015 = fp_startdate, prev_CHD_2015 = prev_CHD, inc_CHD_2015 = inc_CHD, inc_CHD_date_2015 = enddat_CHD)
#2020年版chd数据
CHD_2020 <- read_sav('G:/BaiduSyncdisk/Rcoding/3PA_plaque/CHD - ZW - old prevalent merged_Luoshiyuan Zuo - 10102023.sav')
CHD_2020_0 <- subset(CHD_2020, select = c("ergoid","startdat","Prevalent_CHD","Incident_CHD","end_of_FU_CHD"))
CHD_2020_1 <- rename(CHD_2020_0, prev_CHD_date_2020 = startdat, prev_CHD_2020 = Prevalent_CHD, inc_CHD_2020 = Incident_CHD, inc_CHD_date_2020 = end_of_FU_CHD)
# 
data_check_CHD_1 <- left_join(CHD_2020_1, CHD_2015_1, by="ergoid")
# data_check_CHD_1 <- left_join(ID_total, CHD_2015_1, by="ergoid") %>% 
#   left_join(., CHD_2020_1, by="ergoid")
# data_check_CHD_1$inc_CHD_2020[(data_check_CHD_1$inc_CHD_date_2020 > as.Date("2015-01-01")) & (data_check_CHD_1$inc_CHD_2020 == 1)] <- 0

data_check_CHD_2 <- subset(data_check_CHD_1, prev_CHD_2015 == 1 | prev_CHD_2015 == 0) %>% 
  subset(., prev_CHD_2020 == 1 | prev_CHD_2020 == 0)
data_check_CHD_2$date_check_baseline_CHD <- ifelse(data_check_CHD_2$prev_CHD_2015 == data_check_CHD_2$prev_CHD_2020, 0, 1)
data_check_CHD_ID <- subset(data_check_CHD_2, date_check_baseline_CHD == 1)
dput(names(data_check_CHD_ID))
data_check_CHD_ID <- subset(data_check_CHD_ID, select = c("ergoid", "prev_CHD_2015", "prev_CHD_2020"))
write_sav(data_check_CHD_ID, "G:/BaiduSyncdisk/Rcoding/3PA_plaque/prevalent_CHD_check.sav")

data_check_CHD_3 <- subset(data_check_CHD_1, inc_CHD_2015 == 1 & inc_CHD_2020 == 0)
data_check_CHD_4 <- subset(data_check_CHD_3, select = c("ergoid", "inc_CHD_2015", "inc_CHD_date_2015", "inc_CHD_2020", "inc_CHD_date_2020"))
write_sav(data_check_CHD_4, "G:/BaiduSyncdisk/Rcoding/3PA_plaque/incident_CHD_check.sav")



#stroke数据对比
#2015年版chd数据
stroke_2015 <- read_sav('G:/BaiduSyncdisk/Rcoding/3PA_plaque/Strokes RSI-III 01-01-2016 (28-08-2017)(1).sav')
dput(names(stroke_2015))
stroke_2015_0 <- subset(stroke_2015, select = c("ergoid", "startdat", "prev_stroke_2016", "incid_stroke_2016", "eventdate_2016"))
stroke_2015_1 <- rename(stroke_2015_0, prev_stroke_date_2015 = startdat, prev_stroke_2015 = prev_stroke_2016, inc_stroke_2015 = incid_stroke_2016, inc_stroke_date_2015 = eventdate_2016)

#2020年版stroke数据
stroke_2020 <- read_sav('G:/BaiduSyncdisk/Rcoding/3PA_plaque/Luoshiyuan_Z_Stroke.sav')
dput(names(stroke_2020))
stroke_2020_0 <- subset(stroke_2020, select = c("ergoid", "prevalent_stroke","Baseline_date","incident_stroke","stroke_date"))
stroke_2020_1 <- rename(stroke_2020_0, prev_stroke_date_2020 = Baseline_date, prev_stroke_2020 = prevalent_stroke, inc_stroke_2020 = incident_stroke, inc_stroke_date_2020 = stroke_date)

incident_CVD_sup_1 <- subset(incident_CVD_sup, select = c("ergoid", "ergocar_date"))

data_check_stroke_1 <- left_join(ID_total, stroke_2015_1, by="ergoid") %>% 
  left_join(., stroke_2020_1, by="ergoid") %>% 
  left_join(., incident_CVD_sup_1, by="ergoid")

data_check_stroke_1$inc_stroke_2020[(data_check_stroke_1$inc_stroke_date_2020 > as.Date("2016-01-01")) & (data_check_stroke_1$inc_stroke_2020 == 1)] <- 0

data_check_stroke_1$prev_stroke_2015[(data_check_stroke_1$inc_stroke_date_2015 <= data_check_stroke_1$ergocar_date)&(data_check_stroke_1$inc_stroke_2015 == 1)] <- 1

table_stroke_check <- table(data_check_stroke_1$prev_stroke_2015, data_check_stroke_1$prev_stroke_2020)
names(dimnames(table_stroke_check)) <- c("prev_stroke_2015", "prev_stroke_2020")
table_stroke_check


table_stroke_check <- table(data_check_stroke_1$inc_stroke_2015, data_check_stroke_1$inc_stroke_2020)
names(dimnames(table_stroke_check)) <- c("inc_stroke_2015", "inc_stroke_2020")
table_stroke_check


