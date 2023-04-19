## Menopausia
library(lubridate)
library(tidyverse)
library(tableone)
library(kableExtra)

pobl <- readRDS("//epofs/apistillo/EPOFS/SIDIAP/20180514 OBECAN - Talita/Bases de datos/poblacio.rds")
probl1<-  readRDS("//epofs/apistillo/EPOFS/SIDIAP/20180514 OBECAN - Talita/Bases de datos/problemes.rds")
var_assir <- readRDS("//epofs/apistillo/EPOFS/SIDIAP/20180514 OBECAN - Talita/Bases de datos/variables_assir.rds")
var_clin <- readRDS("//epofs/apistillo/EPOFS/SIDIAP/20180514 OBECAN - Talita/Bases de datos/variables_cliniques.rds")

#we only keep women

pobl <-pobl %>% filter(sexe=="D")


#we create a function of age
fun.age<-function(date1,date2){
  #date1 birth
  output<-year(date2)-year(date1)
  ifelse(month(date2)<month(date1)|(month(date2)==month(date1)&day(date2)<day(date1)),output-1,
         output)
}



################################
#####  First criteria: a Diagnostis of Menopause
################################

# filtramos por el agrupador de manopausia y obtenemos informacion del a?o de nacimiento
probl1 <- probl1%>%filter(str_detect(agr, "MNP"))


probl1 <- probl1%>% left_join(pobl %>% select(id, dnaix),by="id")


# Calculamos la edad del diagn?stico
# y nos quedamos con aquellas mujeres que tengan diagn?sticos entre los 45 y 55

probl1<-probl1%>%
  mutate(age=fun.age(dnaix, dat)) %>% 
  filter(age>=45&age<=55)

## Nos quedamos con el primer registro para cada mujer
n0_meno_dx <- nrow(probl1 %>% distinct(id))

probl1 <- probl1 %>% group_by(id) %>% filter(dat == min(dat)) %>%
  ungroup()
n1_meno_dx <- nrow(probl1 %>% distinct(id))


probl1 <- probl1 %>%
  group_by(id) %>%
  mutate(n_diagnoses = 1:n())

probl1 <- probl1 %>% group_by(id) %>% filter(n_diagnoses<=1) %>%
  ungroup()

n2_meno_dx <- nrow(probl1 %>% distinct(id))

#nos aseguramos que el a?o de diagnostico no sea luego de la salida del sidiap

probl1 <- probl1%>% left_join(pobl %>% select(id, sortida),by="id")

probl1 <- probl1 %>% filter(dat<=sortida)


################################
#####   Age of menopause: data comes from assir and clinical variables ----
################################

var_assir <- var_assir%>%filter(str_detect(agr, "MENOPAUSA_ASSIR"))

age_meno <- var_assir%>%filter(str_detect(cod, "PAFG003"))

library(Hmisc)
describe(age_meno$val)


## Nos quedamos con los valores entre 45 y 55
age_meno<-age_meno%>%
  filter(val>=45&val<=55)



## Nos quedamos con el primer registro para cada mujer
n0_meno_agemeno <- nrow(age_meno %>% distinct(id))

age_meno <- age_meno %>% group_by(id) %>% filter(dat == min(dat)) %>%
  ungroup()
n1_meno_agemeno <- nrow(age_meno %>% distinct(id))


age_meno <- age_meno %>%
  group_by(id) %>%
  mutate(n_diagnoses = 1:n())

age_meno <- age_meno %>% group_by(id) %>% filter(n_diagnoses<=1) %>%
  ungroup()

n2_meno_agemeno <- nrow(age_meno %>% distinct(id))


## Primero debemos a?adir la fecha de nacimiento y salida
## Calularemos que no sea algun a?o posterior a su salida.


age_meno <- age_meno%>% left_join(pobl %>% select(id, dnaix, sortida),by="id")

age_meno <- age_meno%>% mutate(year_birth=as.numeric(format(dnaix,'%Y')))
                                 
age_meno <-age_meno%>% mutate(year_meno= year_birth+val)


## Comprobamos que ningun a?o sea posterior al a?o de salida

age_meno <- age_meno%>% mutate(year_exit=as.numeric(format(sortida,'%Y')))

age_meno <- age_meno %>% filter(year_meno<=year_exit)

n3_meno_agemeno <- nrow(age_meno %>% distinct(id))


age_meno <- age_meno %>% mutate(month_meno="01")

age_meno <- age_meno %>% mutate(day_meno="01")

age_meno$date_meno <- as.Date(with(age_meno, paste(year_meno, month_meno, day_meno)), "%Y%m%d")
                                                  
################################
## hacemos lo mismo para la variable clinica -------


var_clin <- var_clin%>%filter(str_detect(agr, "edat_menopausia"))
describe(var_clin$val)


## Nos quedamos con los valores entre 45 y 55
age_meno2<-var_clin%>%
  filter(val>=45&val<=55)

## Nos quedamos con el primer registro para cada mujer
n0_agemeno2 <- nrow(age_meno2 %>% distinct(id))

age_meno2 <- age_meno2 %>% group_by(id) %>% filter(dat == min(dat)) %>%
  ungroup()
n1_agemeno2 <- nrow(age_meno2 %>% distinct(id))


age_meno2 <- age_meno2 %>%
  group_by(id) %>%
  mutate(n_diagnoses = 1:n())

age_meno2 <- age_meno2 %>% group_by(id) %>% filter(n_diagnoses<=1) %>%
  ungroup()

n2_agemeno2 <- nrow(age_meno2 %>% distinct(id))



## Primero debemos a?adir la fecha de nacimiento y salida
## Calularemos que no sea algun a?o posterior a su salida.


age_meno2 <- age_meno2%>% left_join(pobl %>% select(id, dnaix, sortida),by="id")


age_meno2 <- age_meno2%>% mutate(year_birth=as.numeric(format(dnaix,'%Y')))


age_meno2 <-age_meno2%>% mutate(year_meno= year_birth+val)



## Comprobamos que ningun a?o sea posterior al a?o de salida

age_meno2 <- age_meno2%>% mutate(year_exit=as.numeric(format(sortida,'%Y')))

age_meno2 <- age_meno2 %>% filter(year_meno<=year_exit)

n3_agemeno2 <- nrow(age_meno2 %>% distinct(id))


age_meno2 <- age_meno2 %>% mutate(month_meno="01")

age_meno2 <- age_meno2 %>% mutate(day_meno="01")

age_meno2$date_meno <- as.Date(with(age_meno2, paste(year_meno, month_meno, day_meno)), "%Y%m%d")


################################
#####   COD PACI017 -- last date of period or menopause ------
################################

####### Variable Menopausia


last_period <- var_assir%>%filter(str_detect(cod, "PACI017"))

## a?adimos a fecha de nacimineto


last_period <- last_period %>% left_join(pobl %>% select(id, dnaix, sortida),by="id")

last_period<-last_period%>%
  mutate(age=fun.age(dnaix, dat)) %>% 
  filter(age>=45&age<=55)


##Nos quedamos con el primer registro
n0_last_period <- nrow(last_period %>% distinct(id))

last_period <- last_period %>% group_by(id) %>% filter(dat == min(dat)) %>%
  ungroup()
n1_last_period <- nrow(last_period %>% distinct(id))


last_period <- last_period %>%
  group_by(id) %>%
  mutate(n_diagnoses = 1:n())

last_period <- last_period %>% group_by(id) %>% filter(n_diagnoses<=1) %>%
  ungroup()

n2_last_period <- nrow(last_period %>% distinct(id))


## Comprobamos que ningun a?o sea posterior al a?o de salida

last_period <- last_period %>% filter(dat<=sortida)

n3_last_period <- nrow(last_period %>% distinct(id))

last_period  <- last_period  %>% rename(date_meno=dat)


################################
#####   COD PAED007 -- Year of menopause -----
################################



year_meno <- var_assir%>%filter(str_detect(cod, "PAED007"))


## a?adimos a fecha de nacimineto


year_meno <- year_meno %>% left_join(pobl %>% select(id, dnaix, sortida),by="id")


year_meno<-year_meno%>%
  mutate(age=fun.age(dnaix, dat)) %>% 
  filter(age>=45&age<=55)


##Nos quedamos con el primer registro
n0_year_meno <- nrow(year_meno %>% distinct(id))

year_meno <- year_meno %>% group_by(id) %>% filter(dat == min(dat)) %>%
  ungroup()
n1_year_meno <- nrow(year_meno %>% distinct(id))


year_meno <- year_meno %>%
  group_by(id) %>%
  mutate(n_diagnoses = 1:n())

year_meno <- year_meno %>% group_by(id) %>% filter(n_diagnoses<=1) %>%
  ungroup()

n2_year_meno <- nrow(year_meno %>% distinct(id))


## Comprobamos que ningun a?o sea posterior al a?o de salida

year_meno <- year_meno %>% filter(dat<=sortida)

n3_year_meno <- nrow(year_meno %>% distinct(id))

year_meno  <- year_meno  %>% rename(date_meno=dat)

###### we save the datasets


saveRDS(probl1,"meno1_dx.rds")

saveRDS(age_meno,"meno2_age.rds")

saveRDS(age_meno2,"meno3_age.rds")

saveRDS(last_period,"meno4_lastperiod.rds")

saveRDS(year_meno,"meno5_year.rds")
