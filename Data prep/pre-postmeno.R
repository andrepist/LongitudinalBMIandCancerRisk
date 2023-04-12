#We add menopause
library(tidyverse, lubridate)

fun.age<-function(date1,date2){
  #date1 birth
  output<-year(date2)-year(date1)
  ifelse(month(date2)<month(date1)|(month(date2)==month(date1)&day(date2)<day(date1)),output-1,
         output)
}

pop.new <- readRDS("E:/Andrea/Trajectories OBECAN/populationCox.rds")

meno1_dx <- readRDS("E:/Andrea/Trajectories OBECAN/meno1_dx.rds")

meno2_age <- readRDS("E:/Andrea/Trajectories OBECAN/meno2_age.rds")

meno3_age <- readRDS("E:/Andrea/Trajectories OBECAN/meno3_age.rds")

meno4_lastperiod <- readRDS("E:/Andrea/Trajectories OBECAN/meno4_lastperiod.rds")

meno5_year <- readRDS("E:/Andrea/Trajectories OBECAN/meno5_year.rds")

# 1. we add diagnosis of meno

pop.new <- pop %>% ungroup() %>% rename(date_cancer=dat) %>% filter(.imp==1) %>% select(-.imp) %>% 
  left_join(pop.gen %>% select(id, dnaix))
pop.new<-pop.new %>% left_join(meno1_dx %>% select(id, dat),by="id")

pop.new<-pop.new%>%rename(date_menopause=dat)


# 2. we add age at menopause 1

#hacemos una bbdd de personas sin meno

pop.new.nomeno<-pop.new %>% filter(is.na(date_menopause))
n_nomeno1 <- nrow(pop.new.nomeno %>% distinct(id))

# We add the menopause varible to the dataset of people with no menopause

pop.new.nomeno <- pop.new.nomeno %>% left_join(meno2_age %>% select(id, date_meno),by="id")
n_nomeno2 <- nrow(pop.new.nomeno %>% distinct(id))

pop.new.nomeno<-pop.new.nomeno %>% filter(!is.na(date_meno))
#n_nomeno3 <- nrow(pop.new.nomeno %>% distinct(id))

# We add the date of meno varible to the main dataset

pop.new <- pop.new %>% left_join(pop.new.nomeno %>% select(id, date_meno),by="id")

#we add the date of the 2nd menopause variable to the main dataset

pop.new <- mutate(pop.new, date_menopause= if_else(is.na(date_menopause), date_meno, date_menopause))

pop.new <- select (pop.new, -date_meno)

# 3. we add age at menopause 2

#hacemos una bbdd de personas sin meno

pop.new.nomeno<-pop.new %>% filter(is.na(date_menopause))
#n_nomeno1 <- nrow(pop.new.nomeno %>% distinct(id))

# We add the menopause varible to the dataset of people with no menopause

pop.new.nomeno <- pop.new.nomeno %>% left_join(meno3_age %>% select(id, date_meno),by="id")
#n_nomeno2 <- nrow(pop.new.nomeno %>% distinct(id))

pop.new.nomeno<-pop.new.nomeno %>% filter(!is.na(date_meno))
#n_nomeno3 <- nrow(pop.new.nomeno %>% distinct(id))

# We add the date of meno varible to the main dataset

pop.new <- pop.new %>% left_join(pop.new.nomeno %>% select(id, date_meno),by="id")

#we add the date of the 2nd menopause variable to the main dataset

pop.new <- mutate(pop.new, date_menopause= if_else(is.na(date_menopause), date_meno, date_menopause))

pop.new <- select (pop.new, -date_meno)


# 4. we add last period information

#hacemos una bbdd de personas sin meno

pop.new.nomeno<-pop.new %>% filter(is.na(date_menopause))
#n_nomeno1 <- nrow(pop.new.nomeno %>% distinct(id))

# We add the menopause varible to the dataset of people with no menopause

pop.new.nomeno <- pop.new.nomeno %>% left_join(meno4_lastperiod %>% select(id, date_meno),by="id")
#n_nomeno2 <- nrow(pop.new.nomeno %>% distinct(id))

pop.new.nomeno<-pop.new.nomeno %>% filter(!is.na(date_meno))
#n_nomeno3 <- nrow(pop.new.nomeno %>% distinct(id))

# We add the date of meno varible to the main dataset

pop.new <- pop.new %>% left_join(pop.new.nomeno %>% select(id, date_meno),by="id")

#we add the date of the 2nd menopause variable to the main dataset

pop.new <- mutate(pop.new, date_menopause= if_else(is.na(date_menopause), date_meno, date_menopause))

pop.new <- select (pop.new, -date_meno)


# 5.  we add year of menopause

#hacemos una bbdd de personas sin meno

pop.new.nomeno<-pop.new %>% filter(is.na(date_menopause))
#n_nomeno1 <- nrow(pop.new.nomeno %>% distinct(id))

# We add the menopause varible to the dataset of people with no menopause

pop.new.nomeno <- pop.new.nomeno %>% left_join(meno5_year %>% select(id, date_meno),by="id")
#n_nomeno2 <- nrow(pop.new.nomeno %>% distinct(id))

pop.new.nomeno<-pop.new.nomeno %>% filter(!is.na(date_meno))
#n_nomeno3 <- nrow(pop.new.nomeno %>% distinct(id))

# We add the date of meno varible to the main dataset

pop.new <- pop.new %>% left_join(pop.new.nomeno %>% select(id, date_meno),by="id")

#we add the date of the 2nd menopause variable to the main dataset

pop.new <- mutate(pop.new, date_menopause= if_else(is.na(date_menopause), date_meno, date_menopause))

pop.new <- select (pop.new, -date_meno)


# The variable we constructed, date_menopause, will be used in combination with date of breast cancer diagnosis
#to stratify that cancer into pre and post menopause

# here we create a menopausal variable at index date

#we see which women are 50 or more at index and have a date of menopause prior to index date

pop.new<-pop.new%>%mutate(fifty_index=ifelse(sexe=="D" & age.index>=50, 1,
                                             ifelse(sexe=="D"& age.index<50, 0, NA)))

index.date_date <- ymd("2009-01-01")
pop.new<-pop.new%>%mutate(meno_index=ifelse(sexe=="D" & date_menopause<=index.date_date & (!is.na(date_menopause)), 1,
                                            ifelse(sexe=="D"& ((is.na(date_menopause))| date_menopause>index.date_date), 0, NA)))


pop.new<-pop.new%>%mutate(meno_status_index=ifelse(sexe=="D" & meno_index==1 | fifty_index==1, 1,
                                                   ifelse(sexe=="D"& (meno_index==0 & fifty_index==0), 0, NA)))

pop.new %>% group_by(sexe) %>% count(meno_status_index) %>% ungroup()

# meno status index para discriminar pre post meno at index

#____________________________________________________________________________________________________________________
#We create obesity related cancers 


# First, we divide breast cancer into pre and post menopausal


pop.new<-pop.new%>%
  mutate(age.end=fun.age(dnaix, end)) 


pop.new<-pop.new%>%mutate(fifty_breast_can=ifelse(group=="Breast" & age.end>=50, 1,
                                                  ifelse(group=="Breast" & age.end<50, 0, NA)))


pop.new<-pop.new%>%mutate(meno_breast_can=ifelse(group=="Breast" & date_menopause<=date_cancer & (!is.na(date_menopause)), 1,
                                                 ifelse(group=="Breast"& ((is.na(date_menopause))| date_menopause>date_cancer), 0, NA)))


pop.new<-pop.new%>%mutate(meno_status_breastcan=ifelse(group=="Breast" & meno_breast_can==1 | fifty_breast_can==1, 1,
                                                       ifelse(group=="Breast" & (meno_breast_can==0 & fifty_breast_can==0), 0, NA)))


pop.new<-pop.new%>%mutate(group=ifelse(group=="Breast" & meno_status_breastcan==1, "Breast post",
                                              ifelse(group=="Breast" & meno_status_breastcan==0, "Breast pre", group)))

pop.new<-pop.new%>% select (-fifty_index, -meno_index, -fifty_breast_can, -meno_breast_can)

#first we check the sexe distribution of cancers (we make sure no women cancers happen in men nor viceversa)

pop.new %>% group_by(sexe) %>% count(group) %>%print (n=70) %>% ungroup()




saveRDS(pop.new,"E:/Andrea/Trajectories OBECAN/populationMenopause.rds")
