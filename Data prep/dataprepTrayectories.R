# dataprep

library(lubridate)
library(tidyverse)
library(tableone)
library(kableExtra)

fun.age<-function(date1,date2){
  #date1 birth
  output<-year(date2)-year(date1)
  ifelse(month(date2)<month(date1)|(month(date2)==month(date1)&day(date2)<day(date1)),output-1,
         output)
}


pobl <- readRDS("//epofs/apistillo/EPOFS/SIDIAP/20180514 OBECAN - Talita/Bases de datos/poblacio.rds")
var <- readRDS("//epofs/apistillo/EPOFS/SIDIAP/20180514 OBECAN - Talita/Bases de datos/variables_cliniques.rds")
socio <- readRDS("//epofs/apistillo/EPOFS/SIDIAP/20180514 OBECAN - Talita/Bases de datos/variables_socioeconomiques.rds")
probl<- readRDS("//epofs/apistillo/EPOFS/SIDIAP/20180514 OBECAN - Talita/Bases de datos/problemes.rds")


n0 <- nrow(pobl %>% distinct(id)) 


#poblacion que entra y se va el mismo dia
pop <- pobl%>%
  filter(!sortida==entrada)

# index.date 1 year after entering

pop <- pop %>% mutate(index.date=entrada+365)

# adults at index date
pop <- pop %>% mutate(age.index=fun.age(dnaix, index.date))
pop <- pop %>% filter(age.index>=18)

# cancer -----
all_cancer <- readRDS("E:/Andrea/Trajectories OBECAN/cancer_sidiap.rds")
all_cancer_cmbd <- readRDS("E:/Andrea/Trajectories OBECAN/cancer_cmbd.rds")
cancer <- all_cancer %>% select(id, dat, group, general) %>% mutate(sidiap=1) %>% 
  rbind(all_cancer_cmbd %>% select(id, dat=date_cancer_cmbd, group=group_cmbd, general=general_cmbd) %>% mutate(sidiap=2))

cancer1 <- cancer %>% arrange(dat) %>% distinct(id, .keep_all=TRUE)
# cancer_test <-  all_cancer %>% select(id, dat, group, general)%>% 
#   left_join((all_cancer_cmbd %>% select(id, date_cancer_cmbd, group=group_cmbd, general_cmbd) )) %>% 
#   filter(!is.na(date_cancer_cmbd)) %>% 
#   mutate(diff=as.numeric(date_cancer_cmbd-dat))
# hist(cancer_test$diff)

# End of follow-up: exit from SIDIAP, death, cancer Dx (any but C44) or end of study (31st dec 2018).

pop <- pop %>% left_join(cancer1 %>% select(-general, -sidiap))
pop <- pop %>% mutate(dat.fake=if_else(is.na(dat), ymd("2018/12/31"), dat),
                      end = pmin(sortida, dat.fake, ymd("2018/12/31")))
summary(pop)

# add variables ----

# cvd -----

# wholeCVD <- all_CVD %>% select(dat, id, group) %>% 
#   left_join(all_CVD_cmbd %>% select(dat, id, group), by="id")
#wholeCVD <- wholeCVD %>% filter(dat.x>=dat.y) %>% mutate(diff=as.numeric(dat.y-dat.x))

wholeCVD <- all_CVD %>% select(dat, id, group) %>% 
  rbind(all_CVD_cmbd %>% select(dat, id, group)) %>% 
  arrange(dat) %>% distinct(id, .keep_all=TRUE)
wholeCVD <- wholeCVD %>% transmute(id, date.cvd=dat, cvd=1)

pop <- pop %>% left_join(wholeCVD)
pop <- pop %>% mutate(cvd=ifelse(is.na(cvd), 0, cvd))

# hta -------
hta <- HTA %>% transmute(id, date.hta=dat, hta=1)
pop <- pop %>% left_join(hta)
pop <- pop %>% mutate(hta=ifelse(is.na(hta), 0, hta))

# dm --------
dm <- DM %>% transmute(id, date.dm=dat, dm=1)
pop <- pop %>% left_join(dm)
pop <- pop %>% mutate(dm=ifelse(is.na(dm), 0, dm))

# smoking ------- 

# fumador Fecha más cercana a la fecha índice: cualquier info previa, 1 año despúes fecha Índice
smoking <- readRDS("//epofs/apistillo/EPOFS/SIDIAP/20180514 OBECAN - Talita/Bases de datos/tabac.rds")
index.date2 <- ymd("2009-01-01")
smoking <- smoking %>% left_join(pop %>% select(id, index.date)) %>% 
  #filter(dat>=index.date-years(5),dat<=index.date+years(1)) %>% 
  arrange(abs(dat-index.date2))%>%
  distinct(id, .keep_all = TRUE) %>% mutate(smoking=ifelse(val==0, "Never smoker", 
                                                           ifelse(val==1, "Current smoker",
                                                                  ifelse(val==2, "Previous smoker", "Missing"))))

smoking <- smoking %>% select(-dbaixa,-dat,-val)
pop <- pop %>% left_join(smoking)

rm(smoking)

# alcohol ---------

# alcohol Fecha más cercana a la fecha índice: cualquier info previa, 1 año despúes fecha Índice

clini <- var%>%filter(agr=="ALCOHOL")

View(head(clini))

clini <- clini %>% left_join(pop %>% select(id, index.date)) %>% 
  arrange(abs(dat-index.date2))%>%
  distinct(id, .keep_all = TRUE) %>%  mutate(alcohol=ifelse(val==0, "No alcohol", 
                                                            ifelse(val==1, "Alcohol small risk",
                                                                   ifelse(val==2, "Alcohol high risk", "Missing"))))%>%
  select(id,alcohol)%>%
  mutate(alcohol=as.factor(alcohol))

pop<-pop%>%left_join(clini)
rm(clini)

# charlson ---------

charl <- readRDS("//epofs/apistillo/EPOFS/SIDIAP/20180514 OBECAN - Talita/Bases de datos/charlson.rds")

charl <- charl%>% left_join(pop %>% select(id, index.date)) %>% 
  arrange(abs(dat-index.date2))%>%
  distinct(id, .keep_all = TRUE)%>%
  select(id, index_ch) %>% 
  mutate(index_ch=ifelse(index_ch >=3, "3+", as.character(index_ch)),
         index_ch=as.factor(index_ch)) 

pop<-pop%>%left_join(charl)



#  nationality variable ------


pop<-pop%>%mutate(nationality=ifelse(agr_pais==1, "Spanish", 
                                     ifelse(agr_pais==2|agr_pais==3|agr_pais==4|agr_pais==5, "Global North",
                                            ifelse(agr_pais==6|agr_pais==7|agr_pais==8|agr_pais==9, "Global South", NA))))


# We add MEDEA --------- 

ses_var <- readRDS("//epofs/apistillo/EPOFS/SIDIAP/20180514 OBECAN - Talita/Bases de datos/variables_socioeconomiques.rds")


pop <- pop %>% left_join(ses_var %>% select(-ind_aquas)%>%distinct())
rm(ses_var)

pop<-pop%>%mutate(qmedea=ifelse(qmedea=="", "missing", qmedea))




# years of follow up ------

pop<-pop%>%mutate(follow_up=as.numeric((end-index.date)/365.25))

# at least 1 year
pop <- pop %>% filter(follow_up>=1)


saveRDS(pop,"popSinBMI.rds")

# add all BMI ----

#mini exploracio, gent amb bmi te pes i talla i contrari
bmi<-var %>% filter(agr=="IMC") %>% transmute(id, bmi=1) %>% distinct()
pes<-var %>% filter(agr=="PES") %>% transmute(id, pes=1) %>% distinct()
talla<-var %>% filter(agr=="TALLA") %>% transmute(id, talla=1) %>% distinct()
expl <- pobl %>% 
  left_join(bmi) %>% 
  left_join(pes) %>% 
  left_join(talla)
expl <- expl %>% mutate(bmi=ifelse(is.na(bmi), 0 ,bmi),
                        pes=ifelse(is.na(pes), 0 ,pes),
                        talla=ifelse(is.na(talla), 0 ,talla))
round(prop.table(table(expl$bmi, expl$pes, dnn = c("BMI", "Peso")), margin=2)*100,2)
round(prop.table(table(expl$bmi, expl$pes, dnn = c("BMI", "Peso")), margin=1)*100,2)

table(expl$bmi)
prop.table(table((expl %>% filter(bmi==1))$pes))

expl <- expl %>% mutate(age=fun.age(dnaix, as.Date("2018-12-31")))

expl <- expl %>% 
  mutate(age.gr=ifelse(age<=39, "-39", 
                      
                              ifelse(age>39&age<=59, "40-59",
                                     ifelse(age>59&age<=79, "60-79",
                                            ifelse(age>=80,"80+", "Missing"))))) %>% 
  mutate(age.gr=as.factor(age.gr))

round(prop.table(table(expl$age.gr, expl$pes, expl$bmi,dnn = c("Edad", "Peso","BMI")), margin=c(1,3))*100,2)
table(expl$age.gr, expl$pes, expl$bmi,dnn = c("Edad", "Peso","BMI"))

#end of mini exploracio
rm(talla)
rm(pes)

bmi<-var %>% filter(agr=="IMC")%>%
  select(-cod, -agr)%>%rename(date.bmi=dat,bmi=val)

bmi<-bmi%>%filter(bmi>=15)
bmi <- bmi %>% filter(year(date.bmi)>2005)

# also, drop any BMI records from when someone was under 18
bmi<-bmi %>% 
  left_join(pop %>% select(id, dnaix)) %>% 
  mutate(age.bmi=fun.age(dnaix, date.bmi)) %>% 
  filter(age.bmi>=18) 


data.bmi<-pop %>% 
  left_join(bmi %>% select(id, date.bmi, bmi, age.bmi))

saveRDS(data.bmi, "dataBMI_indexdate.rds")

# new ---------

# new ---------
dataBMI_indexdate1<-readRDS("dataBMI_indexdate.rds")
# add rurality to medea
dataBMI_indexdate <- dataBMI_indexdate %>% mutate(qmedea=ifelse(qmedea=="missing"&ruralitat=="R", "Rural", qmedea))

#drom bmi during pregnancy
emb <- readRDS("//epofs/apistillo/EPOFS/SIDIAP/20180514 OBECAN - Talita/Bases de datos/OBECAN_EMB.rds")
emb <- emb %>% group_by(id) %>% arrange(dur) %>% mutate(n.emb=row_number())
emb <- emb %>% ungroup()
#test <-  dataBMI_indexdate 
for (i in 1:max(emb$n.emb)){
  dataBMI_indexdate <- dataBMI_indexdate%>% left_join(emb %>% filter(n.emb==i))
  dataBMI_indexdate <- dataBMI_indexdate %>%
    mutate(bmi=ifelse((!is.na(dur))&date.bmi>=dur+months(3)&date.bmi<=datfi+months(2), NA, bmi))
  dataBMI_indexdate <- dataBMI_indexdate %>% select(-dur, -datfi,-n.emb)
}

# nrow(dataBMI_indexdate %>% filter(is.na(bmi)) %>% distinct(id))/nrow(test %>% distinct(id))
# nrow(test %>% filter(is.na(bmi)) %>% distinct(id))/nrow(test %>% distinct(id))


# add bariatric surgery ------
bar<- readRDS("//epofs/apistillo/EPOFS/SIDIAP/20180514 OBECAN - Talita/Bases de datos/cmbdah_dx.rds")
bariatric <- bar %>% filter(agr=="Bariatric_surgery")
bariatric <- bariatric %>% mutate(date.bariatric=ymd(dat)) %>% arrange(date.bariatric) %>% 
  distinct(idp, .keep_all=TRUE) %>% 
  select(idp, date.bariatric)

dataBMI_indexdate <- dataBMI_indexdate %>% left_join(bariatric) %>% 
  mutate(bariatric=ifelse(is.na(date.bariatric), 0, 1))

dataBMI_indexdate <- dataBMI_indexdate %>% 
  mutate(age.bmi=ifelse(is.na(bmi), NA, age.bmi))

saveRDS(dataBMI_indexdate, "dataBMI_indexdate_new.rds")

#extra -------

# n of visits per year -------

visits <- readRDS("//epofs/apistillo/EPOFS/SIDIAP/20180514 OBECAN - Talita/Bases de datos/visites.rds")

visits1 <- visits %>% 
  filter(id%in%pop$id) 

visits1 <- visits1%>%
  select(id,dat) %>% distinct()

visits2 <- visits1 %>%
  group_by(id) %>% 
  mutate(n.visits=n()) 


pop.new<-pop.new%>%left_join(visits2, by="idp")
rm(visits1,visits2)

pop.new<-pop.new%>%mutate(n.visits.prior = replace_na(n.visits.prior, 0))
pop.new<-mutate(pop.new, n_visits_index= n.visits.prior)
pop.new<-pop.new%>%select(-n.visits.prior)

rm(visits)

saveRDS(pop.new,"//EPOFS/mrecalde/Adcocans/study_population.rds")

# Age at index date -----

pop.new <- pop.new %>% mutate(pop.new, age_index= age.start)
pop.new <- pop.new %>% select(-age.start)

# change end of follow up -------
# pon situacion al final

pop.new <- pop.new %>%
  mutate(situacio=if_else(!is.na(date_cancer_obesity) & date_cancer_obesity <= exit_date, "Obesity cancer", as.character(situacio)),
         situacio=as.factor(situacio))

pop.new <- pop.new %>%
  mutate(situacio=if_else(!is.na(date_cancer_non_obesity) & date_cancer_non_obesity <= exit_date, "Non-Obesity cancer", as.character(situacio)),
         situacio=as.factor(situacio))

# do the same we did for situacio, but with exit date (exit date must be corrected if the cause of exit is cancer)
pop.new <- pop.new %>%
  mutate(exit_date=if_else(!is.na(date_cancer_obesity)&date_cancer_obesity <= exit_date, date_cancer_obesity, as.Date(exit_date)),
         exit_date=as.Date(exit_date))

pop.new <- pop.new %>%
  mutate(exit_date=if_else(!is.na(date_cancer_non_obesity)&date_cancer_non_obesity <= exit_date, date_cancer_non_obesity, as.Date(exit_date)),
         exit_date=as.Date(exit_date))