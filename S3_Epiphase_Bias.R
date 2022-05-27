
##########################################################################################
#                                                                                        #
#   sensitivity ANALYSES: S1: Epidemic phase bias effect on age-group specific HRs       #
#                                                                                        #
##########################################################################################

#Epidemic phase bias sensitivity analysis: Shift start of follow up, dropping individuals who fall outside the study time interval, shift those with any outcome (hosp, attend, died_cod)

voc_startFolupSens <- voc %>%
  mutate(
    inclIntervalLower = date.start.inclusion,
    inclIntervalUpper = date.data,
    inclIntervalLowerPLUS1 = inclIntervalLower+1,
    inclIntervalLowerPLUS2 = inclIntervalLower+2,
    inclIntervalLowerPLUS3 = inclIntervalLower+3,
    inclIntervalLowerPLUS4 = inclIntervalLower+4,
    
    specimen_date_PLUS1outc=if_else(HospAdm2dLoS=="Yes" | AttendanceHosp=="Yes" | died_28codyn=="Yes",specimen_date+1,specimen_date),
    specimen_date_PLUS2outc=if_else(HospAdm2dLoS=="Yes" | AttendanceHosp=="Yes" | died_28codyn=="Yes",specimen_date+2,specimen_date),
    specimen_date_PLUS3outc=if_else(HospAdm2dLoS=="Yes" | AttendanceHosp=="Yes" | died_28codyn=="Yes",specimen_date+3,specimen_date),
    specimen_date_PLUS4outc=if_else(HospAdm2dLoS=="Yes" | AttendanceHosp=="Yes" | died_28codyn=="Yes",specimen_date+4,specimen_date),
    week_PLUS1outc=isoweek(specimen_date_PLUS1outc),
    week_PLUS2outc=isoweek(specimen_date_PLUS2outc),
    week_PLUS3outc=isoweek(specimen_date_PLUS3outc),
    week_PLUS4outc=isoweek(specimen_date_PLUS4outc)
    
  )


#C0 - hospital admission

cxmod <- coxph( Surv(timeHospAdm2dLoS, HospAdm2dLoS=="Yes") ~ ageGrp10yr.Children5yo_10yr:I(1*(variantConfProbSGTF=="BA2")) + 
                  strata(specimen_date_DAY) + 
                  strata(ageGrp10yr.Children5yo_10yr) + 
                  strata(vaccGp_6grp) +  
                  strata(nhser_name) + 
                  ageGrp10yr.Children5yo_10yr:age +
                  reinfection_status + 
                  sexHC + 
                  ethGrp4InclUnkn + 
                  imd_quintileGp + 
                  imd_quintileGp:imd_rankNUM, 
                data = voc )
hosp0 <- cbind(exp(coef(cxmod)),exp(confint(cxmod)))
hosp0 <- hosp0[10:19, ]
hosp0 <- as.data.frame(hosp0) %>% 
  rename(HR_0 = V1, 
         `2.5 %_0` = `2.5 %`,
         `97.5 %_0` = `97.5 %`)
hosp0 <- tibble::rownames_to_column(hosp0)

cxmod_all <- coxph(Surv(timeHospAdm2dLoS, HospAdm2dLoS =="Yes") ~ variantConfProbSGTF + 
                     strata(specimen_date_DAY) + 
                     strata(ageGrp10yr.Children5yo_10yr) + 
                     strata(vaccGp_6grp) +  
                     strata(nhser_name) + 
                     ageGrp10yr.Children5yo_10yr:age +
                     reinfection_status + 
                     sexHC + 
                     ethGrp4InclUnkn + 
                     imd_quintileGp + 
                     imd_quintileGp:imd_rankNUM, 
                   data = voc)
hosp0_all <- cbind(exp(coef(cxmod_all)),exp(confint(cxmod_all)))
hosp0_all <- as.data.frame(hosp0_all) %>%
  rename(HR_0 = V1, 
         `2.5 %_0` = `2.5 %`,
         `97.5 %_0` = `97.5 %`) 
hosp0_all <- tibble::rownames_to_column(hosp0_all) %>%
  dplyr::filter(rowname == "variantConfProbSGTFBA2") %>%
  rbind(hosp0)

rm(hosp0)
hosp0 <- hosp0_all
rm(hosp0_all)

#C1 - admission  hospital admission

cxmod <- coxph( Surv(timeHospAdm2dLoS, HospAdm2dLoS=="Yes") ~ ageGrp10yr.Children5yo_10yr:I(1*(variantConfProbSGTF=="BA2")) + 
                  strata(specimen_date_PLUS1outc) + 
                  strata(ageGrp10yr.Children5yo_10yr) + 
                  strata(vaccGp_6grp) +  
                  strata(nhser_name) + 
                  ageGrp10yr.Children5yo_10yr:age +
                  reinfection_status + 
                  sexHC + 
                  ethGrp4InclUnkn + 
                  imd_quintileGp + 
                  imd_quintileGp:imd_rankNUM, 
                data = voc_startFolupSens %>% 
                  filter(inclIntervalLowerPLUS1 <= specimen_date_PLUS1outc & specimen_date_PLUS1outc <= inclIntervalUpper) 
)
hosp1 <- cbind(exp(coef(cxmod)),exp(confint(cxmod)))
hosp1 <- hosp1[10:19, ]
hosp1 <- as.data.frame(hosp1) %>% 
  rename(HR_1 = V1, 
         `2.5 %_1` = `2.5 %`,
         `97.5 %_1` = `97.5 %`)
hosp1 <- tibble::rownames_to_column(hosp1)

cxmod_all <- coxph(Surv(timeHospAdm2dLoS, HospAdm2dLoS =="Yes") ~ variantConfProbSGTF + 
                     strata(specimen_date_PLUS1outc) + 
                     strata(ageGrp10yr.Children5yo_10yr) + 
                     strata(vaccGp_6grp) +  
                     strata(nhser_name) + 
                     ageGrp10yr.Children5yo_10yr:age +
                     reinfection_status + 
                     sexHC + 
                     ethGrp4InclUnkn + 
                     imd_quintileGp + 
                     imd_quintileGp:imd_rankNUM, 
                   data = voc_startFolupSens %>% 
                     filter(inclIntervalLowerPLUS1 <= specimen_date_PLUS1outc & specimen_date_PLUS1outc <= inclIntervalUpper)
)
hosp1_all <- cbind(exp(coef(cxmod_all)),exp(confint(cxmod_all)))
hosp1_all <- as.data.frame(hosp1_all) %>%
  rename(HR_1 = V1, 
         `2.5 %_1` = `2.5 %`,
         `97.5 %_1` = `97.5 %`) 
hosp1_all <- tibble::rownames_to_column(hosp1_all) %>%
  dplyr::filter(rowname == "variantConfProbSGTFBA2") %>%
  rbind(hosp1)

rm(hosp1)
hosp1 <- hosp1_all
rm(hosp1_all)

#C2 -  hospital admission

cxmod <- coxph( Surv(timeHospAdm2dLoS, HospAdm2dLoS=="Yes") ~ ageGrp10yr.Children5yo_10yr:I(1*(variantConfProbSGTF=="BA2")) + 
                  strata(specimen_date_PLUS2outc) + 
                  strata(ageGrp10yr.Children5yo_10yr) + 
                  strata(vaccGp_6grp) +  
                  strata(nhser_name) + 
                  ageGrp10yr.Children5yo_10yr:age +
                  reinfection_status + 
                  sexHC + 
                  ethGrp4InclUnkn + 
                  imd_quintileGp + 
                  imd_quintileGp:imd_rankNUM, 
                data = voc_startFolupSens %>% 
                  filter(inclIntervalLowerPLUS2 <= specimen_date_PLUS2outc & specimen_date_PLUS2outc <= inclIntervalUpper) 
)
hosp2 <- cbind(exp(coef(cxmod)),exp(confint(cxmod)))
hosp2 <- hosp2[10:19, ]
hosp2 <- as.data.frame(hosp2) %>% 
  rename(HR_2 = V1, 
         `2.5 %_2` = `2.5 %`,
         `97.5 %_2` = `97.5 %`)
hosp2 <- tibble::rownames_to_column(hosp2)

cxmod_all <- coxph(Surv(timeHospAdm2dLoS, HospAdm2dLoS =="Yes") ~ variantConfProbSGTF + 
                     strata(specimen_date_PLUS2outc) + 
                     strata(ageGrp10yr.Children5yo_10yr) + 
                     strata(vaccGp_6grp) +  
                     strata(nhser_name) + 
                     ageGrp10yr.Children5yo_10yr:age +
                     reinfection_status + 
                     sexHC + 
                     ethGrp4InclUnkn + 
                     imd_quintileGp + 
                     imd_quintileGp:imd_rankNUM, 
                   data = voc_startFolupSens %>% 
                     filter(inclIntervalLowerPLUS2 <= specimen_date_PLUS2outc & specimen_date_PLUS2outc <= inclIntervalUpper)
)
hosp2_all <- cbind(exp(coef(cxmod_all)),exp(confint(cxmod_all)))
hosp2_all <- as.data.frame(hosp2_all) %>%
  rename(HR_2 = V1, 
         `2.5 %_2` = `2.5 %`,
         `97.5 %_2` = `97.5 %`) 
hosp2_all <- tibble::rownames_to_column(hosp2_all) %>%
  dplyr::filter(rowname == "variantConfProbSGTFBA2") %>%
  rbind(hosp2)

rm(hosp2)
hosp2 <- hosp2_all
rm(hosp2_all)

#C3 -  hospital admission

cxmod <- coxph( Surv(timeHospAdm2dLoS, HospAdm2dLoS=="Yes") ~ ageGrp10yr.Children5yo_10yr:I(1*(variantConfProbSGTF=="BA2")) + 
                  strata(specimen_date_PLUS3outc) + 
                  strata(ageGrp10yr.Children5yo_10yr) + 
                  strata(vaccGp_6grp) +  
                  strata(nhser_name) + 
                  ageGrp10yr.Children5yo_10yr:age +
                  reinfection_status + 
                  sexHC + 
                  ethGrp4InclUnkn + 
                  imd_quintileGp + 
                  imd_quintileGp:imd_rankNUM, 
                data = voc_startFolupSens %>% 
                  filter(inclIntervalLowerPLUS3 <= specimen_date_PLUS3outc & specimen_date_PLUS3outc <= inclIntervalUpper) 
)
hosp3 <- cbind(exp(coef(cxmod)),exp(confint(cxmod)))
hosp3 <- hosp3[10:19, ]
hosp3 <- as.data.frame(hosp3) %>% 
  rename(HR_3 = V1, 
         `2.5 %_3` = `2.5 %`,
         `97.5 %_3` = `97.5 %`)
hosp3 <- tibble::rownames_to_column(hosp3)

cxmod_all <- coxph(Surv(timeHospAdm2dLoS, HospAdm2dLoS =="Yes") ~ variantConfProbSGTF + 
                     strata(specimen_date_PLUS3outc) + 
                     strata(ageGrp10yr.Children5yo_10yr) + 
                     strata(vaccGp_6grp) +  
                     strata(nhser_name) + 
                     ageGrp10yr.Children5yo_10yr:age +
                     reinfection_status + 
                     sexHC + 
                     ethGrp4InclUnkn + 
                     imd_quintileGp + 
                     imd_quintileGp:imd_rankNUM, 
                   data = voc_startFolupSens %>% 
                     filter(inclIntervalLowerPLUS3 <= specimen_date_PLUS3outc & specimen_date_PLUS3outc <= inclIntervalUpper)
)
hosp3_all <- cbind(exp(coef(cxmod_all)),exp(confint(cxmod_all)))
hosp3_all <- as.data.frame(hosp3_all) %>%
  rename(HR_3 = V1, 
         `2.5 %_3` = `2.5 %`,
         `97.5 %_3` = `97.5 %`) 
hosp3_all <- tibble::rownames_to_column(hosp3_all) %>%
  dplyr::filter(rowname == "variantConfProbSGTFBA2") %>%
  rbind(hosp3)

rm(hosp3)
hosp3 <- hosp3_all
rm(hosp3_all)

#C4 -  hospital admission

cxmod <- coxph( Surv(timeHospAdm2dLoS, HospAdm2dLoS=="Yes") ~ ageGrp10yr.Children5yo_10yr:I(1*(variantConfProbSGTF=="BA2")) + 
                  strata(specimen_date_PLUS4outc) + 
                  strata(ageGrp10yr.Children5yo_10yr) + 
                  strata(vaccGp_6grp) +  
                  strata(nhser_name) + 
                  ageGrp10yr.Children5yo_10yr:age +
                  reinfection_status + 
                  sexHC + 
                  ethGrp4InclUnkn + 
                  imd_quintileGp + 
                  imd_quintileGp:imd_rankNUM, 
                data = voc_startFolupSens %>% filter(inclIntervalLowerPLUS4 <= specimen_date_PLUS4outc & specimen_date_PLUS4outc <= inclIntervalUpper) 
)
hosp4 <- cbind(exp(coef(cxmod)),exp(confint(cxmod)))
hosp4 <- hosp4[10:19, ]
hosp4 <- as.data.frame(hosp4) %>% 
  rename(HR_4 = V1, 
         `2.5 %_4` = `2.5 %`,
         `97.5 %_4` = `97.5 %`)
hosp4 <- tibble::rownames_to_column(hosp4)

cxmod_all <- coxph(Surv(timeHospAdm2dLoS, HospAdm2dLoS =="Yes") ~ variantConfProbSGTF + 
                     strata(specimen_date_PLUS4outc) + 
                     strata(ageGrp10yr.Children5yo_10yr) + 
                     strata(vaccGp_6grp) +  
                     strata(nhser_name) + 
                     ageGrp10yr.Children5yo_10yr:age +
                     reinfection_status + 
                     sexHC + 
                     ethGrp4InclUnkn + 
                     imd_quintileGp + 
                     imd_quintileGp:imd_rankNUM, 
                   data = voc_startFolupSens %>% 
                     filter(inclIntervalLowerPLUS4 <= specimen_date_PLUS4outc & specimen_date_PLUS4outc <= inclIntervalUpper)
)
hosp4_all <- cbind(exp(coef(cxmod_all)),exp(confint(cxmod_all)))
hosp4_all <- as.data.frame(hosp4_all) %>%
  rename(HR_4 = V1, 
         `2.5 %_4` = `2.5 %`,
         `97.5 %_4` = `97.5 %`) 
hosp4_all <- tibble::rownames_to_column(hosp4_all) %>%
  dplyr::filter(rowname == "variantConfProbSGTFBA2") %>%
  rbind(hosp4)

rm(hosp4)
hosp4 <- hosp4_all
rm(hosp4_all)


epiphase_adm <-left_join(hosp0, hosp1, by = "rowname") 
epiphase_adm <-left_join(epiphase_adm, hosp2, by = "rowname") 
epiphase_adm <-left_join(epiphase_adm, hosp3, by = "rowname") 
epiphase_adm <-left_join(epiphase_adm, hosp4, by = "rowname") %>%
  mutate ( outcome="admission")

rm(hosp0)
rm(hosp1)
rm(hosp2)
rm(hosp3)
rm(hosp4)
rm(cxmod)
rm(cxmod_all)


# C0 - attendance

cxmod <- coxph( Surv(timeAttendanceHosp, AttendanceHosp=="Yes") ~ ageGrp10yr.Children5yo_10yr:I(1*(variantConfProbSGTF=="BA2")) + 
                  strata(specimen_date_DAY) + 
                  strata(ageGrp10yr.Children5yo_10yr) + 
                  strata(vaccGp_6grp) +  
                  strata(nhser_name) + 
                  ageGrp10yr.Children5yo_10yr:age +
                  reinfection_status + 
                  sexHC + 
                  ethGrp4InclUnkn + 
                  imd_quintileGp + 
                  imd_quintileGp:imd_rankNUM, 
                data = voc )

att0 <- cbind(exp(coef(cxmod)),exp(confint(cxmod)))
att0 <- att0[10:19, ]
att0 <- as.data.frame(att0) %>% 
  rename(HR_0 = V1, 
         `2.5 %_0` = `2.5 %`,
         `97.5 %_0` = `97.5 %`)
att0 <- tibble::rownames_to_column(att0)

cxmod_all <- coxph( Surv(timeAttendanceHosp, AttendanceHosp=="Yes") ~ variantConfProbSGTF + 
                      strata(specimen_date_DAY) + 
                      strata(ageGrp10yr.Children5yo_10yr) + 
                      strata(vaccGp_6grp) +  
                      strata(nhser_name) + 
                      ageGrp10yr.Children5yo_10yr:age +
                      reinfection_status + 
                      sexHC + 
                      ethGrp4InclUnkn + 
                      imd_quintileGp + 
                      imd_quintileGp:imd_rankNUM, 
                    data = voc)

att0_all <- cbind(exp(coef(cxmod_all)),exp(confint(cxmod_all)))
att0_all <- as.data.frame(att0_all) %>%
  rename(HR_0 = V1, 
         `2.5 %_0` = `2.5 %`,
         `97.5 %_0` = `97.5 %`) 
att0_all <- tibble::rownames_to_column(att0_all) %>%
  dplyr::filter(rowname == "variantConfProbSGTFBA2") %>%
  rbind(att0)

rm(att0)
att0 <- att0_all
rm(att0_all)


#C1  - attendance

cxmod <- coxph( Surv(timeAttendanceHosp, AttendanceHosp=="Yes") ~ ageGrp10yr.Children5yo_10yr:I(1*(variantConfProbSGTF=="BA2")) + 
                  strata(specimen_date_PLUS1outc) + 
                  strata(ageGrp10yr.Children5yo_10yr) + 
                  strata(vaccGp_6grp) +  
                  strata(nhser_name) + 
                  ageGrp10yr.Children5yo_10yr:age +
                  reinfection_status + 
                  sexHC + 
                  ethGrp4InclUnkn + 
                  imd_quintileGp + 
                  imd_quintileGp:imd_rankNUM, 
                data = voc_startFolupSens %>% 
                  filter(inclIntervalLowerPLUS1 <= specimen_date_PLUS1outc & specimen_date_PLUS1outc <= inclIntervalUpper) 
)
att1 <- cbind(exp(coef(cxmod)),exp(confint(cxmod)))
att1 <- att1[10:19, ]
att1 <- as.data.frame(att1) %>% 
  rename(HR_1 = V1, 
         `2.5 %_1` = `2.5 %`,
         `97.5 %_1` = `97.5 %`)
att1 <- tibble::rownames_to_column(att1)

cxmod_all <- coxph( Surv(timeAttendanceHosp, AttendanceHosp=="Yes") ~ variantConfProbSGTF + 
                      strata(specimen_date_PLUS1outc) + 
                      strata(ageGrp10yr.Children5yo_10yr) + 
                      strata(vaccGp_6grp) +  
                      strata(nhser_name) + 
                      ageGrp10yr.Children5yo_10yr:age +
                      reinfection_status + 
                      sexHC + 
                      ethGrp4InclUnkn + 
                      imd_quintileGp + 
                      imd_quintileGp:imd_rankNUM, 
                    data = voc_startFolupSens %>% 
                      filter(inclIntervalLowerPLUS1 <= specimen_date_PLUS1outc & specimen_date_PLUS1outc <= inclIntervalUpper)
)

att1_all <- cbind(exp(coef(cxmod_all)),exp(confint(cxmod_all)))
att1_all <- as.data.frame(att1_all) %>%
  rename(HR_1 = V1, 
         `2.5 %_1` = `2.5 %`,
         `97.5 %_1` = `97.5 %`) 
att1_all <- tibble::rownames_to_column(att1_all) %>%
  dplyr::filter(rowname == "variantConfProbSGTFBA2") %>%
  rbind(att1)

rm(att1)
att1 <- att1_all
rm(att1_all)

#C2
cxmod <- coxph( Surv(timeAttendanceHosp, AttendanceHosp=="Yes") ~ ageGrp10yr.Children5yo_10yr:I(1*(variantConfProbSGTF=="BA2")) + 
                  strata(specimen_date_PLUS2outc) + 
                  strata(ageGrp10yr.Children5yo_10yr) + 
                  strata(vaccGp_6grp) +  
                  strata(nhser_name) + 
                  ageGrp10yr.Children5yo_10yr:age +
                  reinfection_status + 
                  sexHC + 
                  ethGrp4InclUnkn + 
                  imd_quintileGp + 
                  imd_quintileGp:imd_rankNUM, 
                data = voc_startFolupSens %>% 
                  filter(inclIntervalLowerPLUS2 <= specimen_date_PLUS2outc & specimen_date_PLUS2outc <= inclIntervalUpper) 
)
att2 <- cbind(exp(coef(cxmod)),exp(confint(cxmod)))
att2 <- att2[10:19, ]
att2 <- as.data.frame(att2) %>% 
  rename(HR_2 = V1, 
         `2.5 %_2` = `2.5 %`,
         `97.5 %_2` = `97.5 %`)
att2 <- tibble::rownames_to_column(att2)

cxmod_all <- coxph( Surv(timeAttendanceHosp, AttendanceHosp=="Yes") ~ variantConfProbSGTF  + 
                      strata(specimen_date_PLUS2outc) + 
                      strata(ageGrp10yr.Children5yo_10yr) + 
                      strata(vaccGp_6grp) +  
                      strata(nhser_name) + 
                      ageGrp10yr.Children5yo_10yr:age +
                      reinfection_status + 
                      sexHC + 
                      ethGrp4InclUnkn + 
                      imd_quintileGp + 
                      imd_quintileGp:imd_rankNUM, 
                    data = voc_startFolupSens %>% 
                      filter(inclIntervalLowerPLUS2 <= specimen_date_PLUS2outc & specimen_date_PLUS2outc <= inclIntervalUpper) 
)

att2_all <- cbind(exp(coef(cxmod_all)),exp(confint(cxmod_all)))
att2_all <- as.data.frame(att2_all) %>%
  rename(HR_2 = V1, 
         `2.5 %_2` = `2.5 %`,
         `97.5 %_2` = `97.5 %`) 
att2_all <- tibble::rownames_to_column(att2_all) %>%
  dplyr::filter(rowname == "variantConfProbSGTFBA2") %>%
  rbind(att2)

rm(att2)
att2 <- att2_all
rm(att2_all)



#C3  - attendance
cxmod <- coxph( Surv(timeAttendanceHosp, AttendanceHosp=="Yes") ~ ageGrp10yr.Children5yo_10yr:I(1*(variantConfProbSGTF=="BA2")) + 
                  strata(specimen_date_PLUS3outc) + 
                  strata(ageGrp10yr.Children5yo_10yr) + 
                  strata(vaccGp_6grp) +  
                  strata(nhser_name) + 
                  ageGrp10yr.Children5yo_10yr:age +
                  reinfection_status + 
                  sexHC + 
                  ethGrp4InclUnkn + 
                  imd_quintileGp + 
                  imd_quintileGp:imd_rankNUM, 
                data = voc_startFolupSens %>% 
                  filter(inclIntervalLowerPLUS3 <= specimen_date_PLUS3outc & specimen_date_PLUS3outc <= inclIntervalUpper) 
)
att3 <- cbind(exp(coef(cxmod)),exp(confint(cxmod)))
att3 <- att3[10:19, ]
att3 <- as.data.frame(att3) %>% 
  rename(HR_3 = V1, 
         `2.5 %_3` = `2.5 %`,
         `97.5 %_3` = `97.5 %`)
att3 <- tibble::rownames_to_column(att3)


cxmod_all <- coxph( Surv(timeAttendanceHosp, AttendanceHosp=="Yes") ~ variantConfProbSGTF  + 
                      strata(specimen_date_PLUS3outc) + 
                      strata(ageGrp10yr.Children5yo_10yr) + 
                      strata(vaccGp_6grp) +  
                      strata(nhser_name) + 
                      ageGrp10yr.Children5yo_10yr:age +
                      reinfection_status + 
                      sexHC + 
                      ethGrp4InclUnkn + 
                      imd_quintileGp + 
                      imd_quintileGp:imd_rankNUM, 
                    data = voc_startFolupSens %>% 
                      filter(inclIntervalLowerPLUS3 <= specimen_date_PLUS3outc & specimen_date_PLUS3outc <= inclIntervalUpper) 
)

att3_all <- cbind(exp(coef(cxmod_all)),exp(confint(cxmod_all)))
att3_all <- as.data.frame(att3_all) %>%
  rename(HR_3 = V1, 
         `2.5 %_3` = `2.5 %`,
         `97.5 %_3` = `97.5 %`) 
att3_all <- tibble::rownames_to_column(att3_all) %>%
  dplyr::filter(rowname == "variantConfProbSGTFBA2") %>%
  rbind(att3)

rm(att3)
att3 <- att3_all
rm(att3_all)


#C4  - attendance
cxmod <- coxph( Surv(timeAttendanceHosp, AttendanceHosp=="Yes") ~ ageGrp10yr.Children5yo_10yr:I(1*(variantConfProbSGTF=="BA2")) + 
                  strata(specimen_date_PLUS4outc) + 
                  strata(ageGrp10yr.Children5yo_10yr) + 
                  strata(vaccGp_6grp) +  
                  strata(nhser_name) + 
                  ageGrp10yr.Children5yo_10yr:age +
                  reinfection_status + 
                  sexHC + 
                  ethGrp4InclUnkn + 
                  imd_quintileGp + 
                  imd_quintileGp:imd_rankNUM, 
                data = voc_startFolupSens %>% filter(inclIntervalLowerPLUS4 <= specimen_date_PLUS4outc & specimen_date_PLUS4outc <= inclIntervalUpper) 
)
att4 <- cbind(exp(coef(cxmod)),exp(confint(cxmod)))
att4 <- att4[10:19, ]
att4 <- as.data.frame(att4) %>% 
  rename(HR_4 = V1, 
         `2.5 %_4` = `2.5 %`,
         `97.5 %_4` = `97.5 %`)
att4 <- tibble::rownames_to_column(att4)

cxmod_all <- coxph( Surv(timeAttendanceHosp, AttendanceHosp=="Yes") ~ variantConfProbSGTF  + 
                      strata(specimen_date_PLUS4outc) + 
                      strata(ageGrp10yr.Children5yo_10yr) + 
                      strata(vaccGp_6grp) +  
                      strata(nhser_name) + 
                      ageGrp10yr.Children5yo_10yr:age +
                      reinfection_status + 
                      sexHC + 
                      ethGrp4InclUnkn + 
                      imd_quintileGp + 
                      imd_quintileGp:imd_rankNUM, 
                    data = voc_startFolupSens %>% filter(inclIntervalLowerPLUS4 <= specimen_date_PLUS4outc & specimen_date_PLUS4outc <= inclIntervalUpper) 
)

att4_all <- cbind(exp(coef(cxmod_all)),exp(confint(cxmod_all)))
att4_all <- as.data.frame(att4_all) %>%
  rename(HR_4 = V1, 
         `2.5 %_4` = `2.5 %`,
         `97.5 %_4` = `97.5 %`) 
att4_all <- tibble::rownames_to_column(att4_all) %>%
  dplyr::filter(rowname == "variantConfProbSGTFBA2") %>%
  rbind(att4)

rm(att4)
att4 <- att4_all
rm(att4_all)

epiphase_att <-left_join(att0, att1, by = "rowname") 
epiphase_att <-left_join(epiphase_att, att2, by = "rowname") 
epiphase_att <-left_join(epiphase_att, att3, by = "rowname") 
epiphase_att <-left_join(epiphase_att, att4, by = "rowname") %>%
  mutate ( outcome="attendance")


rm(att0)
rm(att1)
rm(att2)
rm(att3)
rm(att4)
rm(cxmod)



# C0 - death

cxmod <- coxph( Surv(timeDied28, died_28yn == "Yes") ~ ageGrp10yr:I(1*(variantConfProbSGTF=="BA2")) + 
                  strata(specimen_date_DAY) + 
                  strata(ageGrp10yr) + 
                  strata(vaccGp_6grp) +  
                  strata(nhser_name) + 
                  ageGrp10yr:age +
                  reinfection_status + 
                  sexHC + 
                  ethGrp4InclUnkn + 
                  imd_quintileGp + 
                  imd_quintileGp:imd_rankNUM, 
                data = voc )

death0 <- cbind(exp(coef(cxmod)),exp(confint(cxmod)))
death0 <- death0[10:19, ]
death0 <- as.data.frame(death0) %>% 
  rename(HR_0 = V1, 
         `2.5 %_0` = `2.5 %`,
         `97.5 %_0` = `97.5 %`)
death0 <- tibble::rownames_to_column(death0)

cxmod_all <- coxph( Surv(timeDied28, died_28yn == "Yes") ~ variantConfProbSGTF + 
                      strata(specimen_date_DAY) + 
                      strata(ageGrp10yr) + 
                      strata(vaccGp_6grp) +  
                      strata(nhser_name) + 
                      ageGrp10yr:age +
                      reinfection_status + 
                      sexHC + 
                      ethGrp4InclUnkn + 
                      imd_quintileGp + 
                      imd_quintileGp:imd_rankNUM, 
                    data = voc )

death0_all <- cbind(exp(coef(cxmod_all)),exp(confint(cxmod_all)))
death0_all <- as.data.frame(death0_all) %>%
  rename(HR_0 = V1, 
         `2.5 %_0` = `2.5 %`,
         `97.5 %_0` = `97.5 %`) 
death0_all <- tibble::rownames_to_column(death0_all) %>%
  dplyr::filter(rowname == "variantConfProbSGTFBA2") %>%
  rbind(death0)

rm(death0)
death0 <- death0_all
rm(death0_all)

#C1  - death

cxmod <- coxph( Surv(timeDied28, died_28yn == "Yes") ~ ageGrp10yr:I(1*(variantConfProbSGTF=="BA2")) + 
                  strata(specimen_date_PLUS1outc) + 
                  strata(ageGrp10yr) + 
                  strata(vaccGp_6grp) +  
                  strata(nhser_name) + 
                  ageGrp10yr:age +
                  reinfection_status + 
                  sexHC + 
                  ethGrp4InclUnkn + 
                  imd_quintileGp + 
                  imd_quintileGp:imd_rankNUM, 
                data = voc_startFolupSens %>% 
                  filter(inclIntervalLowerPLUS1 <= specimen_date_PLUS1outc & specimen_date_PLUS1outc <= inclIntervalUpper) 
)
death1 <- cbind(exp(coef(cxmod)),exp(confint(cxmod)))
death1 <- death1[10:19, ]
death1 <- as.data.frame(death1) %>% 
  rename(HR_1 = V1, 
         `2.5 %_1` = `2.5 %`,
         `97.5 %_1` = `97.5 %`)
death1 <- tibble::rownames_to_column(death1)

cxmod_all <- coxph( Surv(timeDied28, died_28yn == "Yes") ~ variantConfProbSGTF + 
                      strata(specimen_date_PLUS1outc) + 
                      strata(ageGrp10yr) + 
                      strata(vaccGp_6grp) +  
                      strata(nhser_name) + 
                      ageGrp10yr:age +
                      reinfection_status + 
                      sexHC + 
                      ethGrp4InclUnkn + 
                      imd_quintileGp + 
                      imd_quintileGp:imd_rankNUM, 
                    data = voc_startFolupSens %>% 
                      filter(inclIntervalLowerPLUS1 <= specimen_date_PLUS1outc & specimen_date_PLUS1outc <= inclIntervalUpper) 
)

death1_all <- cbind(exp(coef(cxmod_all)),exp(confint(cxmod_all)))
death1_all <- as.data.frame(death1_all) %>%
  rename(HR_1 = V1, 
         `2.5 %_1` = `2.5 %`,
         `97.5 %_1` = `97.5 %`) 
death1_all <- tibble::rownames_to_column(death1_all) %>%
  dplyr::filter(rowname == "variantConfProbSGTFBA2") %>%
  rbind(death1)

rm(death1)
death1 <- death1_all
rm(death1_all)

#C2
cxmod <- coxph( Surv(timeDied28, died_28yn == "Yes") ~ ageGrp10yr:I(1*(variantConfProbSGTF=="BA2")) + 
                  strata(specimen_date_PLUS2outc) + 
                  strata(ageGrp10yr) + 
                  strata(vaccGp_6grp) +  
                  strata(nhser_name) + 
                  ageGrp10yr:age +
                  reinfection_status + 
                  sexHC + 
                  ethGrp4InclUnkn + 
                  imd_quintileGp + 
                  imd_quintileGp:imd_rankNUM, 
                data = voc_startFolupSens %>% 
                  filter(inclIntervalLowerPLUS2 <= specimen_date_PLUS2outc & specimen_date_PLUS2outc <= inclIntervalUpper) 
)
death2 <- cbind(exp(coef(cxmod)),exp(confint(cxmod)))
death2 <- death2[10:19, ]
death2 <- as.data.frame(death2) %>% 
  rename(HR_2 = V1, 
         `2.5 %_2` = `2.5 %`,
         `97.5 %_2` = `97.5 %`)
death2 <- tibble::rownames_to_column(death2)

cxmod_all <- coxph( Surv(timeDied28, died_28yn == "Yes") ~ variantConfProbSGTF  + 
                      strata(specimen_date_PLUS2outc) + 
                      strata(ageGrp10yr) + 
                      strata(vaccGp_6grp) +  
                      strata(nhser_name) + 
                      ageGrp10yr:age +
                      reinfection_status + 
                      sexHC + 
                      ethGrp4InclUnkn + 
                      imd_quintileGp + 
                      imd_quintileGp:imd_rankNUM, 
                    data = voc_startFolupSens %>% 
                      filter(inclIntervalLowerPLUS2 <= specimen_date_PLUS2outc & specimen_date_PLUS2outc <= inclIntervalUpper) 
)

death2_all <- cbind(exp(coef(cxmod_all)),exp(confint(cxmod_all)))
death2_all <- as.data.frame(death2_all) %>%
  rename(HR_2 = V1, 
         `2.5 %_2` = `2.5 %`,
         `97.5 %_2` = `97.5 %`) 
death2_all <- tibble::rownames_to_column(death2_all) %>%
  dplyr::filter(rowname == "variantConfProbSGTFBA2") %>%
  rbind(death2)

rm(death2)
death2 <- death2_all
rm(death2_all)

#C3  - death
cxmod <- coxph( Surv(timeDied28, died_28yn == "Yes") ~ ageGrp10yr:I(1*(variantConfProbSGTF=="BA2")) + 
                  strata(specimen_date_PLUS3outc) + 
                  strata(ageGrp10yr) + 
                  strata(vaccGp_6grp) +  
                  strata(nhser_name) + 
                  ageGrp10yr:age +
                  reinfection_status + 
                  sexHC + 
                  ethGrp4InclUnkn + 
                  imd_quintileGp + 
                  imd_quintileGp:imd_rankNUM, 
                data = voc_startFolupSens %>% 
                  filter(inclIntervalLowerPLUS3 <= specimen_date_PLUS3outc & specimen_date_PLUS3outc <= inclIntervalUpper) 
)
death3 <- cbind(exp(coef(cxmod)),exp(confint(cxmod)))
death3 <- death3[10:19, ]
death3 <- as.data.frame(death3) %>% 
  rename(HR_3 = V1, 
         `2.5 %_3` = `2.5 %`,
         `97.5 %_3` = `97.5 %`)
death3 <- tibble::rownames_to_column(death3)


cxmod_all <- coxph( Surv(timeDied28, died_28yn == "Yes") ~ variantConfProbSGTF  + 
                      strata(specimen_date_PLUS3outc) + 
                      strata(ageGrp10yr) + 
                      strata(vaccGp_6grp) +  
                      strata(nhser_name) + 
                      ageGrp10yr:age +
                      reinfection_status + 
                      sexHC + 
                      ethGrp4InclUnkn + 
                      imd_quintileGp + 
                      imd_quintileGp:imd_rankNUM, 
                    data = voc_startFolupSens %>% 
                      filter(inclIntervalLowerPLUS3 <= specimen_date_PLUS3outc & specimen_date_PLUS3outc <= inclIntervalUpper) 
)

death3_all <- cbind(exp(coef(cxmod_all)),exp(confint(cxmod_all)))
death3_all <- as.data.frame(death3_all) %>%
  rename(HR_3 = V1, 
         `2.5 %_3` = `2.5 %`,
         `97.5 %_3` = `97.5 %`) 
death3_all <- tibble::rownames_to_column(death3_all) %>%
  dplyr::filter(rowname == "variantConfProbSGTFBA2") %>%
  rbind(death3)

rm(death3)
death3 <- death3_all
rm(death3_all)

#C4  - death
cxmod <- coxph( Surv(timeDied28, died_28yn == "Yes") ~ ageGrp10yr:I(1*(variantConfProbSGTF=="BA2")) + 
                  strata(specimen_date_PLUS4outc) + 
                  strata(ageGrp10yr) + 
                  strata(vaccGp_6grp) +  
                  strata(nhser_name) + 
                  ageGrp10yr:age +
                  reinfection_status + 
                  sexHC + 
                  ethGrp4InclUnkn + 
                  imd_quintileGp + 
                  imd_quintileGp:imd_rankNUM, 
                data = voc_startFolupSens %>% filter(inclIntervalLowerPLUS4 <= specimen_date_PLUS4outc & specimen_date_PLUS4outc <= inclIntervalUpper) 
)
death4 <- cbind(exp(coef(cxmod)),exp(confint(cxmod)))
death4 <- death4[10:19, ]
death4 <- as.data.frame(death4) %>% 
  rename(HR_4 = V1, 
         `2.5 %_4` = `2.5 %`,
         `97.5 %_4` = `97.5 %`)
death4 <- tibble::rownames_to_column(death4)

cxmod_all <- coxph( Surv(timeDied28, died_28yn == "Yes") ~ variantConfProbSGTF  + 
                      strata(specimen_date_PLUS4outc) + 
                      strata(ageGrp10yr) + 
                      strata(vaccGp_6grp) +  
                      strata(nhser_name) + 
                      ageGrp10yr:age +
                      reinfection_status + 
                      sexHC + 
                      ethGrp4InclUnkn + 
                      imd_quintileGp + 
                      imd_quintileGp:imd_rankNUM, 
                    data = voc_startFolupSens %>% 
                      filter(inclIntervalLowerPLUS4 <= specimen_date_PLUS4outc & specimen_date_PLUS4outc <= inclIntervalUpper) 
)

death4_all <- cbind(exp(coef(cxmod_all)),exp(confint(cxmod_all)))
death4_all <- as.data.frame(death4_all) %>%
  rename(HR_4 = V1, 
         `2.5 %_4` = `2.5 %`,
         `97.5 %_4` = `97.5 %`) 
death4_all <- tibble::rownames_to_column(death4_all) %>%
  dplyr::filter(rowname == "variantConfProbSGTFBA2") %>%
  rbind(death4)

rm(death4)
death4 <- death4_all
rm(death4_all)

epiphase_death <-left_join(death0, death1, by = "rowname") 
epiphase_death <-left_join(epiphase_death, death2, by = "rowname") 
epiphase_death <-left_join(epiphase_death, death3, by = "rowname") 
epiphase_death <-left_join(epiphase_death, death4, by = "rowname") %>%
  mutate ( outcome="death28")

rm(death0)
rm(death1)
rm(death2)
rm(death3)
rm(death4)
rm(cxmod)

epiphase <- rbind(epiphase_adm,
                  epiphase_att,
                  epiphase_death)
epiphase <- as_grouped_data(x = epiphase, groups = "outcome")

epiphase <- epiphase %>%
  mutate(rowname = case_when(
    (grepl("variantConfProbSGTFBA2",rowname) ~ "All ages"),
    (grepl("0,5)",rowname) ~ "<5"),
    (grepl("5,10",rowname) ~ "5-9"),
    (grepl("10,20",rowname) ~ "10-19"),
    (grepl("20,30",rowname) ~ "20-29"),
    (grepl("30,40",rowname) ~ "30-39"),
    (grepl("40,50",rowname) ~ "40-49"),
    (grepl("50,60",rowname) ~ "50-59"),
    (grepl("60,70",rowname) ~ "60-69"),
    (grepl("70,80",rowname) ~ "70-79"),
    (grepl("80,Inf",rowname) ~ "\u2265 80")
  ))


epiphase_adm_ft <- flextable(epiphase, c("outcome","rowname", "c0","c1","c2","c3","c4")) %>%
  
  colformat_double(i=2, digits = 3) %>%
  
  set_caption(caption = "S3: Effect of epidemic phase bias on ascertainment of risk of hospital admission, by age group") %>%
  
  compose(j = "c0",
          value = as_paragraph(as_chunk(sprintf("%.2f",`HR_0`)), as_chunk(" ("), as_chunk(sprintf("%.2f",`2.5 %_0`)), as_chunk(" - "), as_chunk(sprintf("%.2f",`97.5 %_0`)), as_chunk(")")
          )) %>%
  
  compose(j = "c1",
          value = as_paragraph(as_chunk(sprintf("%.2f",`HR_1`)), as_chunk(" ("), as_chunk(sprintf("%.2f",`2.5 %_1`)), as_chunk(" - "), as_chunk(sprintf("%.2f",`97.5 %_1`)), as_chunk(")")
          )) %>%
  
  compose(j = "c2",
          value = as_paragraph(as_chunk(sprintf("%.2f",`HR_2`)), as_chunk(" ("), as_chunk(sprintf("%.2f",`2.5 %_2`)), as_chunk(" - "), as_chunk(sprintf("%.2f",`97.5 %_2`)), as_chunk(")")
          ))  %>%
  compose(j = "c3",
          value = as_paragraph(as_chunk(sprintf("%.2f",`HR_3`)), as_chunk(" ("), as_chunk(sprintf("%.2f",`2.5 %_3`)), as_chunk(" - "), as_chunk(sprintf("%.2f",`97.5 %_3`)), as_chunk(")")
          ))  %>%
  compose(j = "c4",
          value = as_paragraph(as_chunk(sprintf("%.2f",`HR_4`)), as_chunk(" ("), as_chunk(sprintf("%.2f",`2.5 %_4`)), as_chunk(" - "), as_chunk(sprintf("%.2f",`97.5 %_4`)), as_chunk(")")
          ))  %>%
  
  set_header_labels("outcome" = "Outcome",
                    "rowname" = "Age Group",
                    "c0" = "0 days (primary analysis) ",
                    "c1"= "1 day", 
                    "c2" = "2 days", 
                    "c3" = "3 days",
                    "c4"= "4 days" ) %>%
  
  set_table_properties(layout = "autofit") %>%
  
  add_header_row(colwidths = c(1,6), 
                 values = c("","HR (CI 95%) Assumed difference in mean time from infection to positive test between cases without outcome and cases with outcome")) %>%
  merge_at(i = 1:1, j = 1:7) %>%
  merge_at(i = 13:13, j = 1:7) %>%
  merge_at(i = 25:25, j = 1:7) 

print(epiphase_adm_ft, preview = "docx")


save_as_docx(epiphase_adm_ft, path = glue("{project}/Tasks/Alternate_Hosp_Rerun/Outputs/S3_epiphase_{Sys.Date()}.docx"))
