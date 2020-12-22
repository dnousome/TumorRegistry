dt_surv=dt_time_fix
dt_surv_check=dt_surv %>% filter(DEERS_Death_Date=="" & death_new==1 &dataset=="ACTUR")
dt_surv_check=dt_surv %>% filter(!pair %in% dt_surv_check$pair) %>% 
  mutate(age_cat5=factor(age_cat5)) %>%
  mutate(age_cat4=recode(age_cat5,
      `35-39`="35-49",`40-44`="35-49",`45-49`='35-49',
      `50-54`="50-59", `55-59`="50-59", 
      `60-64`="60-69",`65-69`="60-69", 
        `70-74`="70+", `75-79`="70+", `80-84`="70+", `85-89`="70+", `90-94`="70+", `95-99`="70+")) 

#HRs were further adjusted for age (as a continuous variable), Hispanic origin, tumor grade, tumor location, surgery, and radiation. 
#In the stratified analysis, all variables adjusted in the overall analysis were adjusted except for the stratified variable itself.

##Across all
library(survival)
summary(coxph(Surv(fiveyear_FU, fiveyear_death) ~dataset+stage_final1+Surgery+radiation_new+age+as.factor(hispanic)+strata(pair),data=dt_surv_check))


##Filter only to t3
t3=dt_surv_check %>% filter(stage_final1==3)
summary(coxph(Surv(fiveyear_FU, fiveyear_death) ~dataset+Surgery+radiation_new+age+as.factor(hispanic),data=t3))

##t4
t4=dt_surv_check %>% filter(stage_final1==4)
summary(coxph(Surv(fiveyear_FU, fiveyear_death) ~dataset+Surgery+radiation_new+age+as.factor(hispanic),data=t4))


lapply(unique(dt_surv_check$age_cat4),function(x){
  print(x)
  temp_dt=dt_surv_check %>% filter(age_cat4 %in% x)
  summary(coxph(Surv(fiveyear_FU, fiveyear_death) ~dataset+,data=temp_dt))$coefficients
})



fit1=survfit(Surv(fiveyear_FU, fiveyear_death)~ dataset+stage_final1+Surgery+radiation_new+age+as.factor(hispanic),data=dt_surv)
ggsurvplot(fit1, data = dt_surv, fun = "cloglog",title="Log-Log plot for Assumptions")

plot(fit1, fun='cloglog')

fit1=coxph(Surv(fiveyear_FU, fiveyear_death) ~dataset+strata(stage_final1)+Surgery+radiation_new+age+as.factor(hispanic)+strata(pair),data=dt_surv)
cox.zph(fit1)
ggcoxdiagnostics(fit1, type = "schoenfeld")
