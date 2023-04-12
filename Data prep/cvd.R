#CVD definition


## CVD definition

angina <- probl%>%filter(str_detect(cod,"I20"))%>%
  mutate(group="Angina", general="CVD")

MI <- probl%>%filter(str_detect(cod,"I21|I22"))%>%
  mutate(group="Myocardial infarction", general="CVD")

other_CHD <- probl%>%filter(str_detect(cod,"I23|I24|I25"))%>%
  mutate(group="Other CHD", general="CVD")


hemo_stroke <- probl%>%filter(str_detect(cod,"I60|I61"))%>%
  mutate(group="Haemorrhagic stroke", general="CVD")

ische_stroke <- probl%>%filter(str_detect(cod,"I63"))%>%
  mutate(group="Cerebral infarction", general="CVD")

unclass_stroke <- probl%>%filter(str_detect(cod,"I64"))%>%
  mutate(group="Unclassified stroke", general="CVD")

oth_acute_cere <- probl%>%filter(str_detect(cod,"I62|I65|I66|I67|I68|I69"))%>%
  mutate(group="Other acute cerebrovascular events", general="CVD")


#We create a CVD variable
all_CVD<-rbind(angina, MI, other_CHD, hemo_stroke, ische_stroke, unclass_stroke, oth_acute_cere) %>% 
  group_by(id) 

#For people with several CVD diagnoses, we keep the first one (although this also keeps 2 dx on the same date)
all_CVD <- all_CVD %>% group_by(id) %>% filter(dat == min(dat)) %>%
  ungroup()

n1_CVD <- nrow(all_CVD %>% distinct(id))
#403,697 people 

#we only keep the first dx per person
all_CVD <- all_CVD %>%
  group_by(id) %>%
  mutate(n_diagnoses = 1:n())

all_CVD <- all_CVD %>% group_by(id) %>% filter(n_diagnoses<=1) %>%
  ungroup()

n2_CVD <- nrow(all_CVD %>% distinct(id))
#403,697 people

#n1_CVD and n2_CVD now have the same N as the number of obs in the all_CVD dataframe

rm (angina, MI, other_CHD, hemo_stroke, ische_stroke, unclass_stroke, oth_acute_cere)

saveRDS(all_CVD,"CVD_sidiap.rds")


##______________________________________________________________________
# definition in cmbd


## CVD definition in CMBD

angina <- probl_cmbd%>%filter(str_detect(cod,"^413"))%>% 
  mutate(group="Angina", general="CVD")

MI <- probl_cmbd%>%filter(str_detect(cod,"^410"))%>%
  mutate(group="Myocardial infarction", general="CVD")

other_CHD <- probl_cmbd%>%filter(str_detect(cod,"^411|^412|^414"))%>%
  mutate(group="Other CHD", general="CVD")

hemo_stroke <- probl_cmbd%>%filter(str_detect(cod,"^430|^431"))%>%
  mutate(group="Haemorrhagic stroke", general="CVD")

#ische_stroke --> no ICD 9 codes

#unclass_stroke --> no ICD 9 codes

oth_acute_cere <- probl_cmbd%>%filter(str_detect(cod,"^432|^433|^434|^436|^437|^438"))%>%
  mutate(group="Other acute cerebrovascular events", general="CVD")

#We create a CVD variable
all_CVD_cmbd<-rbind(angina, MI, other_CHD, hemo_stroke, oth_acute_cere) %>% 
  group_by(id) 

#For people with several CVD diagnoses, we keep the first one (although this also keeps 2 dx on the same date)
all_CVD_cmbd <- all_CVD_cmbd %>% group_by(id) %>% filter(dat == min(dat)) %>%
  ungroup()

n1_CVD_cmbd <- nrow(all_CVD_cmbd %>% distinct(id))
#265,151 people 

#we only keep the first dx per person
all_CVD_cmbd <- all_CVD_cmbd %>%
  group_by(id) %>%
  mutate(n_diagnoses = 1:n())

all_CVD_cmbd <- all_CVD_cmbd %>% group_by(id) %>% filter(n_diagnoses<=1) %>%
  ungroup()

n2_CVD_cmbd <- nrow(all_CVD_cmbd %>% distinct(id))
#265,151 people

#n1_CVD and n2_CVD now have the same N as the number of obs in the all_CVD dataframe

rm (angina, MI, other_CHD, hemo_stroke, oth_acute_cere)

saveRDS(all_CVD_cmbd,"CVD_cmbd.rds")