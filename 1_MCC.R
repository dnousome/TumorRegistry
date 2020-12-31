###############################
#
#MCC ACTUR-SEER PRoject
#By: Darryl Nousome
#
###############################
library(tidyverse)
library(haven)
setwd("~/Documents/G/Work_Bench/Epidemiology/Projects/MCC_CPDR SEERACTUR 2020/")

dt=read_sas("~/Documents/G/Work_Bench/Epidemiology/Projects/MCC_CPDR SEERACTUR 2020/Data/prostate_delivery_pwd.sas7bdat")
#dt=read_csv("~/Documents/G/Work_Bench/Epidemiology/Projects/MCC_CPDR SEERACTUR 2020/Data/dt.csv")

##Summarize the variables and cohort
#Age group- 18‐39 y,40‐54 y,55‐64 y,≥65 y
#Race-White, Black, Asian or Pacific Islander, Other
#?Ethnicity
#Year of diagnosis, 1987‐1989, 1990‐1994, 1995‐1999, 2000‐2004, 2005‐2009,2010‐2013
#Tumor grade, Well differentiated, grade1, Moderately differentiated, grade 2,Poorly differentiated, grade 3, Undifferentiated, grade 4

#Tumor location
#Supratentorial
#Infratentorial
#Other locations
#Cancer‐directed surgery
#No
#Yes
#Unknown or missing
#Radiation,No, Yes, Unknown or missing

#CS Site-Specific Factor 1-ER 
##PSA997	Test ordered, results not in chart
#998	Test not done (test not ordered and not performed)
#999	Unknown or no information,Not documented in patient record
##Clean new radiation

#2994-PSA not mising

#1,2_psa
#6,8 gleason
#10- RPGleason
datefx=function(x)as.Date(x,"%m-%d-%Y")

dt_final=dt %>% 
  mutate(radiation_new=ifelse(dataset=="SEER" & Radiation %in% c("Beam radiation",
      "Combination of beam with implants or isotopes","Radiation, NOS  method or source not specified",
      "Radioactive implants","Radioisotopes"),"Yes",
      ifelse(dataset=="SEER" & Radiation %in% c("None","Refused"),"No",
          ifelse(dataset=="SEER" & Radiation %in% c("Unknown","Recommended, unknown if administered"),"Unknown/Missing",
                 ifelse(dataset=="ACTUR" & RX_Summ__Radiation %in% as.character(c(1:6)),"Yes",
                        ifelse(dataset=="ACTUR" & RX_Summ__Radiation %in% c("8","9",""),"Unknown/Missing",
                              ifelse(dataset=="ACTUR" & RX_Summ__Radiation %in% c("0","7"),"No",Radiation))))))) %>%
  mutate(PSA=ifelse(dataset=="ACTUR" & CS_Site_Specific_Factor_1 %in% c('997','998','999'),NA,
                  ifelse(dataset=="ACTUR",as.numeric(CS_Site_Specific_Factor_1)/10,
                    ifelse(dataset=="SEER" & CSsitespecificfactor12004varyi=="Blank(s)",NA,
                           ifelse(dataset=="SEER" & CSsitespecificfactor12004varyi %in% c('997','998','999'),NA,
                                  ifelse(dataset=="SEER", as.numeric(CSsitespecificfactor12004varyi)/10,NA)))))) %>%
  mutate(PSA_category=ifelse(dataset=="ACTUR" & CS_Site_Specific_Factor_2 %in% c('997','998','999'),NA,
                    ifelse(dataset=="ACTUR" & CS_Site_Specific_Factor_2=='020',"Negative",
                           ifelse(dataset=="ACTUR" & CS_Site_Specific_Factor_2=='010',"Positive",
                                  ifelse(dataset=="ACTUR" & CS_Site_Specific_Factor_2=='030',NA,
                           ifelse(dataset=="SEER" & CSsitespecificfactor22004varyi %in% c("Blank(s)","997","998","999"),NA,
                                  ifelse(dataset=="SEER" & CSsitespecificfactor22004varyi=='020',"Negative",
                                         ifelse(dataset=="SEER" & CSsitespecificfactor22004varyi=='010',"Positive",
                                                ifelse(dataset=="SEER" & CSsitespecificfactor22004varyi=='030',NA,NA))))))))) %>%
  mutate(Gleason_old=ifelse(dataset=="ACTUR" & CS_Site_Specific_Factor_6 %in% c('988','998','999'),NA,
                             ifelse(dataset=="ACTUR" & CS_Site_Specific_Factor_6 %in% c('000','002','004','005'),NA,
                                    ifelse(dataset=="ACTUR", as.numeric(CS_Site_Specific_Factor_6),
                                                  ifelse(dataset=="SEER" & CSsitespecificfactor62004varyi %in% c("Blank(s)","988","999"),NA,
                                                         ifelse(dataset=="SEER" & CSsitespecificfactor62004varyi %in% c('000','002','004','005'),NA,
                                                         ifelse(dataset=="SEER", as.numeric(CSsitespecificfactor62004varyi),NA))))))) %>%
  mutate(Gleason_biopsy=ifelse(dataset=="ACTUR" & CS_Site_Specific_Factor_8 %in% c('988','998','999'),NA,
                            ifelse(dataset=="ACTUR" & CS_Site_Specific_Factor_8 %in% c('000','002','003','004','005'),NA,
                                   ifelse(dataset=="ACTUR", as.numeric(CS_Site_Specific_Factor_8),
                                          ifelse(dataset=="SEER" & CSsitespecificfactor82004varyi %in% c("Blank(s)","988","998","999"),NA,
                                                 ifelse(dataset=="SEER" & CSsitespecificfactor82004varyi %in% c('000','002','003','004','005'),NA,
                                                        ifelse(dataset=="SEER", as.numeric(CSsitespecificfactor82004varyi),NA))))))) %>%
  mutate(Gleason_prostatectomy=ifelse(dataset=="ACTUR" & CS_Site_Specific_Factor10 %in% c('988','998','999'),NA,
                        ifelse(dataset=="ACTUR" & CS_Site_Specific_Factor10 %in% c('000','002','003','004','005'),NA,
                               ifelse(dataset=="ACTUR", as.numeric(CS_Site_Specific_Factor10),
                                      ifelse(dataset=="SEER" & CSsitespecificfactor102004vary %in% c("Blank(s)","988","998","999"),NA,
                                             ifelse(dataset=="SEER" & CSsitespecificfactor102004vary %in% c('000','002','003','004','005'),NA,
                                                    ifelse(dataset=="SEER", as.numeric(CSsitespecificfactor102004vary),NA))))))) %>%
  mutate(Surgery=ifelse(dataset=="ACTUR" & RX_Summ__Surg_Prim_Site %in% c('00'),"No",
                           ifelse(dataset=="ACTUR" & RX_Summ__Surg_Prim_Site %in% c('99',""),"Unknown/Missing",
                                ifelse(dataset=="ACTUR" & RX_Summ__Surg_Prim_Site %in% c('14','19','21','22','23','30','50','70','80','90'),"Yes",
                                       ifelse((dataset=="SEER" & Sitespecificsurgery19731997var %in% c('00','01','02','03','05','06','07')) | 
                                                (dataset=="SEER" & RX_Summ_Surg_Prim_Site_1998 %in% c("00")),"No",
                                              ifelse((dataset=="SEER" & Sitespecificsurgery19731997var %in% c('10','20','30','40','50','60','70','80',"90")) | 
                                                      (dataset=="SEER" &   RX_Summ_Surg_Prim_Site_1998 %in% c('10','14','15','16','17','18','19','20','21','22','23','24','25','26','30','50','70','80','90')),"Yes",
                                                        ifelse((dataset=="SEER" & Sitespecificsurgery19731997var %in% c('09','58','99','Blank(s)')) | 
                                                            (dataset=="SEER" &  RX_Summ_Surg_Prim_Site_1998 %in% c("99",'Blank(s)')),"Unknown/Missing",NA
                                                                 ))))))) %>%
  mutate(Gleason_any=ifelse(Gleason_biopsy %in% c(6:10)|Gleason_prostatectomy %in% c(6:10),"Yes",
                            ifelse(is.na(Gleason_biopsy)|is.na(Gleason_prostatectomy),"No",NA))) %>%
  mutate(Gleason_final=Gleason_prostatectomy) %>%
  mutate(Gleason_final=ifelse(is.na(Gleason_final) & !is.na(Gleason_biopsy),Gleason_biopsy,Gleason_final)) %>%
  mutate(Gleason_final_method=ifelse(!is.na(Gleason_prostatectomy),"Prostatectomy",
                                     ifelse(!is.na(Gleason_biopsy),"Biopsy",NA))) %>%

  mutate(tumor_grade=ifelse(tumor_grade==9,"Unknown/Missing",tumor_grade))   %>%
  mutate(stage_final1=as.factor(stage_final)) %>%
  mutate(PSA_cat=cut(PSA,breaks = c(0,4,10,20,96),include.lowest = F,right = F,labels=c("0-<4","4-<10","10-<20","\u226520"))) %>%
  mutate(Hispanic=recode(hispanic,`0`="Non-Hispanic",`1`="Hispanic",`9`="Unknown")) %>%
  #mutate(Surgery=ifelse(Reasonnocancerdirected_surgery=="Surgery performed" & Surgery=="No","Yes",Surgery)) %>%
  mutate(PSA_cat_rem=ifelse(!(year_cat %in% '2010-2013') & !is.na(PSA_cat),NA,PSA_cat))  %>%
  mutate(Gleason_final_rem=ifelse(!(year_cat %in% '2010-2013') & !is.na(Gleason_final),NA,Gleason_final)) %>%
  
  mutate(age_cat1=cut(age,breaks = c(0,50,65,100),labels=c("<50","50-<65","\u226565"),right = F)) %>%
  mutate(Dxyr_PSA_all=ifelse(year_cat=="2010-2013" & is.na(PSA_cat),8,
                             ifelse(!year_cat=="2010-2013" & is.na(PSA_cat),9,
                         ifelse(PSA_cat=="0-<4",0,
                                ifelse(PSA_cat=="4-<10",1,
                                       ifelse(PSA_cat=="10-<20",2,
                                              ifelse(PSA_cat=="\u226520",3,9
                                                     ))))))) %>%
  
  
  mutate(Dxyr_PSA=ifelse(year_cat=="2010-2013" & is.na(PSA_cat),8,
                                ifelse(year_cat=="2010-2013" & PSA_cat=="0-<4",0,
                         ifelse(year_cat=="2010-2013" & PSA_cat=="4-<10",1,
                                ifelse(year_cat=="2010-2013" &PSA_cat=="10-<20",2,
                                       ifelse(year_cat=="2010-2013" & PSA_cat=="\u226520",3,
                                              
                                                     ifelse(!year_cat=="2010-2013" & is.na(PSA_cat),9,9
                                                     ))))))) %>%
  
  ##Use Gleason first for grade
  mutate(tumor_grade_gleason=ifelse(is.na(Gleason_final_rem),tumor_grade,
                                    ifelse(Gleason_final_rem=='6','1',
                                     ifelse(Gleason_final_rem=='7','2',
                                             ifelse(Gleason_final_rem=='8-10','3',tumor_grade))))) %>%

#  mutate(tumor_grade_gleason=ifelse(tumor_grade_gleason=='4','3',tumor_grade_gleason)) %>%
#  mutate(tumor_grade_gleason=ifelse(is.na(tumor_grade_gleason) | tumor_grade_gleason=='5',"Unknown/Missing",tumor_grade_gleason)) %>%
  mutate(tumor_grade_gleason=factor(tumor_grade_gleason,levels=c("1","2","3","4","Unknown/Missing"),labels=c("1","2","3","4","Unknown/Missing")))


###Check the surgery variable
#dt1=dt %>% 
#  filter(grepl("recommended",Reasonnocancerdirected_surgery)|grepl("Recommended",Reasonnocancerdirected_surgery)|
#                                                                     grepl("Unknown",Reasonnocancerdirected_surgery)) %>% 
dt1=dt_final %>% filter(Reasonnocancerdirected_surgery=="Surgery performed")
table(dt_final$Surgery,dt_final$RX_Summ_Surg_Prim_Site_1998)
table(dt_final$Surgery,dt_final$Sitespecificsurgery19731997var)

dt1=dt_final %>% filter(Reasonnocancerdirected_surgery=="Surgery performed" & dt_final$Surgery=="No")


##Create death_new which is 
dt_time=dt_final %>% 
  mutate(death_new=ifelse(dataset=="ACTUR" & is.na(datefx(DEERS_Death_Date)) & death %in% 1,1,
                        ifelse(dataset=="ACTUR" & datefx(DEERS_Death_Date)>studyend & death %in% 1,0,death))) %>%
  mutate(actur_FU=if_else(last_contact_date>studyend,studyend,last_contact_date)) %>%
  mutate(time_dx_death=ifelse(dataset=="ACTUR" & death_new==1,round((datefx(DEERS_Death_Date)-dx_date_final)/30.25),
                              ifelse(dataset=="ACTUR" & death_new==0,round((actur_FU-dx_date_final)/30.25),surtime_month))) %>%
 # mutate(time_dx_death=ifelse(dataset=="ACTUR" & is.na(time_death_death) & death_new==0))
  mutate(time_dx_death=ifelse(time_dx_death<1,0,time_dx_death))  

##162 deaths occur after 2013 and are censored
dt_time %>% filter(datefx(DEERS_Death_Date)>datefx("12-31-2013"))
##176 no death dates     

no_dod=dt_time %>% filter(death_new==1&is.na(time_dx_death)) %>% select(dataset,pair,dx_date_final,DEERS_Death_Date,last_contact_date,death,death_new,time_dx_death)
###2 with deaths before dx, ACTUR
dt_time %>% filter(datefx(DEERS_Death_Date)<dx_date_final) %>% select(dataset,pair,dx_date_final,DEERS_Death_Date,last_contact_date,death,death_new)


##Check to make sure fix
#dt_time_fix %>% filter(death_new==1&is.na(time_dx_death)) %>% select(dataset,pair,dx_date_final,DEERS_Death_Date,last_contact_date,death,death_new,time_dx_death)
##Set last contact date as death 
#3703, 5212 drop
fix1=dt_time %>% filter(is.na(time_dx_death)) %>% select(time_dx_death,last_contact_date,DEERS_Death_Date,dx_date_final,death,death_new)

dt_time_fix=dt_time %>% 
  filter(!pair %in% c("3703", "5212")) %>%
  mutate(time_dx_death=ifelse(death_new==1 & is.na(time_dx_death),last_contact_date-dx_date_final,time_dx_death)) %>%
  mutate(time_dx_death=ifelse(time_dx_death<1,0,time_dx_death))  %>%
 
  mutate(fiveyear_FU=ifelse(dataset=="SEER" &death_new==1 & time_dx_death<=60,time_dx_death,
                            ifelse(dataset=="SEER" & death_new==1 & time_dx_death>60,60,
                                   ifelse(dataset=="SEER" & death_new==0 & time_dx_death<=60,time_dx_death,
                                          ifelse(dataset=="SEER" & death_new==0 & time_dx_death>60,60,
                                                 ifelse(dataset=="ACTUR" &death_new==1 & time_dx_death<=60,time_dx_death,
                                                        ifelse(dataset=="ACTUR" & death_new==1 & time_dx_death>60,60,
                                                               ifelse(dataset=="ACTUR" & death_new==0 & time_dx_death<=60,time_dx_death,
                                                                      ifelse(dataset=="ACTUR" & death_new==0 & time_dx_death>60,60,NA))))))))) %>%
  mutate(fiveyear_death=ifelse(dataset=="SEER" &death_new==1 & time_dx_death<=60,1,
                            ifelse(dataset=="SEER" & death_new==1 & time_dx_death>60,0,
                                   ifelse(dataset=="SEER" & death_new==0 & time_dx_death<=60,0,
                                          ifelse(dataset=="SEER" & death_new==0 & time_dx_death>60,0,
                                                 ifelse(dataset=="ACTUR" &death_new==1 & time_dx_death<=60,1,
                                                        ifelse(dataset=="ACTUR" & death_new==1 & time_dx_death>60,0,
                                                               ifelse(dataset=="ACTUR" & death_new==0 & time_dx_death<=60,0,
                                                                      ifelse(dataset=="ACTUR" & death_new==0 & time_dx_death>60,0,NA)))))))))



table(dt_time$death_new,dt_time$dataset)
group_by(dt_time,dataset) %>% summarize(mean(actur_FU,na.rm=T))  
group_by(dt_time_fix,dataset) %>% summarize(mean(time_dx_death,na.rm=T),median(time_dx_death,na.rm=T))  

group_by(dt_time_fix,dataset) %>% summarize(mean(actur_FU,na.rm=T))  
group_by(dt_time_fix,dataset) %>% summarize(mean(fiveyear_FU,na.rm=T),mean(fiveyear_death,na.rm=T))  



  
dt1=dt_time %>% 
  select(region,pair,age,age_cat1,race_new,Hispanic,year_cat,Year_of_Diagnosis,age_cat5,radiation_new,
         dataset,PSA,PSA_category,PSA_cat,PSA_cat_rem,
         Gleason_old,Gleason_biopsy,Gleason_prostatectomy,Gleason_any,Gleason_final,Gleason_final_rem,
         Gleason_final_method,tumor_grade_gleason,Dxyr_PSA,Dxyr_PSA_all,tumor_grade,
         Surgery,RX_Summ__Surg_Prim_Site,Sitespecificsurgery19731997var,RX_Summ_Surg_Prim_Site_1998,
         death_new,time_dx_death,stage_final1,Radiation,RX_Summ__Radiation) 


##Dt1 is the full cohort
###dt_time fix removes the 2 matched set and fixes the followup time
write_rds(dt1,"Data/mcc_actur.rds")
write_rds(dt_time_fix,"Data/mcc_actur_surv.rds")















library(survminer)
library(survival)
fit <- survfit(Surv(time_dx_death, death_new) ~ dataset,
               data = dt1)
ggsurvplot(fit, data = dt1, risk.table = TRUE)

summary(dt_final$time_dx_death)
s=filter(dt_final,time_dx_death<0)
datefx(DEERS_Death_Date)-dx_date_final


Studyend: ACTUR study end date.
DEERS_Death_Date: ACTUR date of Death variable
Dx_Date_final: ACTUR date of diagnosis
Last_contact_date: ACTUR date of last contact







lapply(names(dt),function(x){
  print(x)
  s=sum(dt[[x]]==dt1[[x]],na.rm=T)
  print(s)
  })
