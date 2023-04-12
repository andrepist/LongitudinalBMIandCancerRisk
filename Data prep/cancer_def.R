library(lubridate)
library(tidyverse)
library(tableone)
library(kableExtra)

probl<- readRDS("//epofs/apistillo/EPOFS/SIDIAP/20180514 OBECAN - Talita/Bases de datos/problemes.rds")

pobl <- readRDS("//epofs/apistillo/EPOFS/SIDIAP/20180514 OBECAN - Talita/Bases de datos/poblacio.rds")

## cancer definition 

# We define cancers in SIDIAP

head_can <- probl%>%filter(str_detect(cod,"C00|C01|C02|C03|C04|C05|C06|C07|C08|C09|C10|C11|C12|C13|C14"))%>%
  mutate(group="Head and Neck", general="Cancer")

esopha_can <- probl%>%filter(str_detect(cod,"C15"))%>%
  mutate(group="Esophagus", general="Cancer")

stomach_can <- probl%>%filter(str_detect(cod,"C16"))%>%
  mutate(group="Stomach", general="Cancer")

colorectal_can <- probl%>%filter(str_detect(cod,"C18|C19|C20|C21"))%>%
  mutate(group="Colorectal", general="Cancer")

liver_can <- probl%>%filter(str_detect(cod,"C22"))%>%
  mutate(group="Liver", general="Cancer")

gall_can <- probl%>%filter(str_detect(cod,"C23|C24"))%>%
  mutate(group="Gallbladder & biliary tract", general="Cancer")

pancreas_can <- probl%>%filter(str_detect(cod,"C25"))%>%
  mutate(group="Pancreas", general="Cancer")

larynx_can <- probl%>%filter(str_detect(cod,"C32"))%>%
  mutate(group="Larynx", general="Cancer")

lung_can <- probl%>%filter(str_detect(cod,"C33|C34"))%>%
  mutate(group="Trachea, bronchus & Lung", general="Cancer")

bone_can <- probl%>%filter(str_detect(cod,"C40|C41"))%>%
  mutate(group="Bone and articular cartilage", general="Cancer")

mms_can <- probl%>%filter(str_detect(cod,"C43"))%>%
  mutate(group="Malignant melanoma of skin", general="Cancer")

softtissue_can <- probl%>%filter(str_detect(cod,"C47|C49"))%>%
  mutate(group="Connective and soft tissue", general="Cancer")

breast_can <- probl%>%filter(str_detect(cod,"C50"))%>%
  mutate(group="Breast", general="Cancer")
breast_can <- breast_can%>% left_join(pobl %>% select(id, sexe),by="id")
breast_can <- breast_can%>%filter(sexe=="D")
breast_can <- breast_can%>%select(-sexe)

cervix_can <- probl%>%filter(str_detect(cod,"C53"))%>%
  mutate(group="Cervix Uteri", general="Cancer")
cervix_can <- cervix_can%>% left_join(pobl %>% select(id, sexe),by="id")
cervix_can <- cervix_can%>%filter(sexe=="D")
cervix_can <- cervix_can%>%select(-sexe)


uterus_can <- probl%>%filter(str_detect(cod,"C54|C55"))%>%
  mutate(group="Corpus Uteri", general="Cancer")
uterus_can <- uterus_can%>% left_join(pobl %>% select(id, sexe),by="id")
uterus_can <- uterus_can%>%filter(sexe=="D")
uterus_can <- uterus_can%>%select(-sexe)


ovary_can <- probl%>%filter(str_detect(cod,"C56"))%>%
  mutate(group="Ovary", general="Cancer")
ovary_can <- ovary_can%>% left_join(pobl %>% select(id, sexe),by="id")
ovary_can <- ovary_can%>%filter(sexe=="D")
ovary_can <- ovary_can%>%select(-sexe)


penis_can <- probl%>%filter(str_detect(cod,"C60"))%>%
  mutate(group="Penis", general="Cancer")
penis_can <- penis_can%>% left_join(pobl %>% select(id, sexe),by="id")
penis_can <- penis_can%>%filter(sexe=="H")
penis_can <- penis_can%>%select(-sexe)

prostate_can <- probl%>%filter(str_detect(cod,"C61"))%>%
  mutate(group="Prostate", general="Cancer")
prostate_can <- prostate_can%>% left_join(pobl %>% select(id, sexe),by="id")
prostate_can <- prostate_can%>%filter(sexe=="H")
prostate_can <- prostate_can%>%select(-sexe)

testis_can <- probl%>%filter(str_detect(cod,"C62"))%>%
  mutate(group="Testis", general="Cancer")
testis_can <- testis_can%>% left_join(pobl %>% select(id, sexe),by="id")
testis_can <- testis_can%>%filter(sexe=="H")
testis_can <- testis_can%>%select(-sexe)

kidney_can <- probl%>%filter(str_detect(cod,"C64"))%>%
  mutate(group="Kidney", general="Cancer")

urinary_can <- probl%>%filter(str_detect(cod,"C65|C66|C68"))%>%
  mutate(group="Urinary tract", general="Cancer")

bladder_can <- probl%>%filter(str_detect(cod,"C67"))%>%
  mutate(group="Bladder", general="Cancer")

brain_can <- probl%>%filter(str_detect(cod,"C70|C71|C72|C75.1|C75.2|C75.3"))%>%
  mutate(group="Brain and CNS", general="Cancer")

thyroid_can <- probl%>%filter(str_detect(cod,"C73"))%>%
  mutate(group="Thyroid", general="Cancer")

hl_can <- probl%>%filter(str_detect(cod,"C81"))%>%
  mutate(group="Hodgkin lymphoma", general="Cancer")

nhl_can <- probl%>%filter(str_detect(cod,"C82|C83|C84|C85|C86|C96"))%>%
  mutate(group="Non-Hodgkin Lymphoma", general="Cancer")

mmyeloma_can <- probl%>%filter(str_detect(cod,"C90"))%>%
  mutate(group="Multiple myeloma", general="Cancer")

leukemia_can <- probl%>%filter(str_detect(cod,"C91|C92|C93|C94|C95"))%>%
  mutate(group="Leukemia", general="Cancer")

others_can <- probl%>%filter(str_detect(cod,"C17|C26|C30|C31|C37|C38|C39|C4A|C45|C46|C48|C51|
                                        C52|C57|C58|C63|C69|C74|C75.0|C75.4|C75.5|C75.8|C75.9|
                                        C7A|C76|C80|C88|C97"))%>%
  mutate(group="Others and non-specific", general="Cancer")


#We create a cancer variable
all_cancer<-rbind(head_can, 
                  esopha_can, stomach_can, colorectal_can, liver_can, gall_can,
                  pancreas_can, larynx_can, lung_can, bone_can, mms_can, softtissue_can, 
                  breast_can, cervix_can, uterus_can, ovary_can, penis_can, prostate_can, 
                  testis_can, kidney_can, urinary_can, bladder_can, brain_can, thyroid_can, 
                  hl_can, nhl_can, mmyeloma_can, leukemia_can, others_can) %>% 
  group_by(id) 

#For people with several cancer diagnoses, we keep the first one (although this also keeps 2 dx on the same date)
all_cancer <- all_cancer %>% group_by(id) %>% filter(dat == min(dat)) %>%
  ungroup()
n1_canc <- nrow(all_cancer %>% distinct(id))


#we only keep the first dx per person
all_cancer <- all_cancer %>%
  group_by(id) %>%
  mutate(n_diagnoses = 1:n())

all_cancer <- all_cancer %>% group_by(id) %>% filter(n_diagnoses<=1) %>%
  ungroup()

n2_canc <- nrow(all_cancer %>% distinct(id))

#n1_can and n2_can now have the same N as the number of obs in the all_cancer dataframe

rm (head_can, 
    esopha_can, stomach_can, colorectal_can, liver_can, gall_can,
    pancreas_can, larynx_can, lung_can, bone_can, mms_can, softtissue_can, 
    breast_can, cervix_can, uterus_can, ovary_can, penis_can, prostate_can, 
    testis_can, kidney_can, urinary_can, bladder_can, brain_can, thyroid_can, 
    hl_can, nhl_can, mmyeloma_can, leukemia_can, others_can)

saveRDS(all_cancer,"cancer_sidiap.rds")

##______________________________________________________________________________________________

## We define cancers in CMBD


probl_cmbd<- readRDS("//epofs/apistillo/EPOFS/SIDIAP/20180514 OBECAN - Talita/Bases de datos/cmbd_dx_padris.rds")


head_can <- probl_cmbd%>%filter(str_detect(cod,"^140|^141|^142|^143|^144|^145|^146|^147|^148|^149"))%>%
  mutate(group="Head and Neck", general="Cancer")

esopha_can <- probl_cmbd%>%filter(str_detect(cod,"^150"))%>%
  mutate(group="Esophagus", general="Cancer")

stomach_can <- probl_cmbd%>%filter(str_detect(cod,"^151"))%>%
  mutate(group="Stomach", general="Cancer")

colorectal_can <- probl_cmbd%>%filter(str_detect(cod,"^153|^154"))%>%
  mutate(group="Colorectal", general="Cancer")

liver_can <- probl_cmbd%>%filter(str_detect(cod,"^155"))%>%
  mutate(group="Liver", general="Cancer")

gall_can <- probl_cmbd%>%filter(str_detect(cod,"^156"))%>%
  mutate(group="Gallbladder & biliary tract", general="Cancer")

pancreas_can <- probl_cmbd%>%filter(str_detect(cod,"^157"))%>%
  mutate(group="Pancreas", general="Cancer")

larynx_can <- probl_cmbd%>%filter(str_detect(cod,"^161"))%>%
  mutate(group="Larynx", general="Cancer")

lung_can <- probl_cmbd%>%filter(str_detect(cod,"^162"))%>%
  mutate(group="Trachea, bronchus & Lung", general="Cancer")

bone_can <- probl_cmbd%>%filter(str_detect(cod,"^170"))%>%
  mutate(group="Bone and articular cartilage", general="Cancer")

mms_can <- probl_cmbd%>%filter(str_detect(cod,"^172"))%>%
  mutate(group="Malignant melanoma of skin", general="Cancer")

softtissue_can <- probl_cmbd%>%filter(str_detect(cod,"^171"))%>%
  mutate(group="Connective and soft tissue", general="Cancer")

breast_can <- probl_cmbd%>%filter(str_detect(cod,"^174|^175"))%>%
  mutate(group="Breast", general="Cancer")
breast_can <- breast_can%>% left_join(pobl %>% select(id, sexe),by="id")
breast_can <- breast_can%>%filter(sexe=="D")
breast_can <- breast_can%>%select(-sexe)

cervix_can <- probl_cmbd%>%filter(str_detect(cod,"^180"))%>%
  mutate(group="Cervix Uteri", general="Cancer")
cervix_can <- cervix_can%>% left_join(pobl %>% select(id, sexe),by="id")
cervix_can <- cervix_can%>%filter(sexe=="D")
cervix_can <- cervix_can%>%select(-sexe)


uterus_can <- probl_cmbd%>%filter(str_detect(cod,"^179|^182"))%>%
  mutate(group="Corpus Uteri", general="Cancer")
uterus_can <- uterus_can%>% left_join(pobl %>% select(id, sexe),by="id")
uterus_can <- uterus_can%>%filter(sexe=="D")
uterus_can <- uterus_can%>%select(-sexe)

ovary_can <- probl_cmbd%>%filter(str_detect(cod,"^183"))%>%
  mutate(group="Ovary", general="Cancer")
ovary_can <- ovary_can%>% left_join(pobl %>% select(id, sexe),by="id")
ovary_can <- ovary_can%>%filter(sexe=="D")
ovary_can <- ovary_can%>%select(-sexe)


penis_can <- probl_cmbd%>%filter(str_detect(cod,"^1871|^1872|^1873|^1874"))%>%
  mutate(group="Penis", general="Cancer")
penis_can <- penis_can%>% left_join(pobl %>% select(id, sexe),by="id")
penis_can <- penis_can%>%filter(sexe=="H")
penis_can <- penis_can%>%select(-sexe)

prostate_can <- probl_cmbd%>%filter(str_detect(cod,"^185"))%>%
  mutate(group="Prostate", general="Cancer")
prostate_can <- prostate_can%>% left_join(pobl %>% select(id, sexe),by="id")
prostate_can <- prostate_can%>%filter(sexe=="H")
prostate_can <- prostate_can%>%select(-sexe)

testis_can <- probl_cmbd%>%filter(str_detect(cod,"^186"))%>%
  mutate(group="Testis", general="Cancer")
testis_can <- testis_can%>% left_join(pobl %>% select(id, sexe),by="id")
testis_can <- testis_can%>%filter(sexe=="H")
testis_can <- testis_can%>%select(-sexe)

kidney_can <- probl_cmbd%>%filter(str_detect(cod,"^1890"))%>%
  mutate(group="Kidney", general="Cancer")

urinary_can <- probl_cmbd%>%filter(str_detect(cod,"1891|^1892|^1893|^1894|^1898|^1899"))%>%
  mutate(group="Urinary tract", general="Cancer")

bladder_can <- probl_cmbd%>%filter(str_detect(cod,"^188"))%>%
  mutate(group="Bladder", general="Cancer")

brain_can <- probl_cmbd%>%filter(str_detect(cod,"^191|^192|^1943|^1944"))%>%
  mutate(group="Brain and CNS", general="Cancer")

thyroid_can <- probl_cmbd%>%filter(str_detect(cod,"^193"))%>%
  mutate(group="Thyroid", general="Cancer")

hl_can <- probl_cmbd%>%filter(str_detect(cod,"^201"))%>%
  mutate(group="Hodgkin lymphoma", general="Cancer")

nhl_can <- probl_cmbd%>%filter(str_detect(cod,"^200|^202"))%>%
  mutate(group="Non-Hodgkin Lymphoma", general="Cancer")

mmyeloma_can <- probl_cmbd%>%filter(str_detect(cod,"^203"))%>%
  mutate(group="Multiple myeloma", general="Cancer")

leukemia_can <- probl_cmbd%>%filter(str_detect(cod,"204|^205|^206|^207|^208"))%>%
  mutate(group="Leukemia", general="Cancer")

others_can <- probl_cmbd%>%filter(str_detect(cod,"^152|^158|^159|^160|^163|^164|^165|^176|^181|^184|
                                        ^1875|^1876|^1877|^1878|^1879|^190|^1940|^1941|^1945|^1946|^1948|
                                        ^1949|^195|^199|^2091|^2092|^2093|^2733|^2795"))%>%
  mutate(group="Others and non-specific", general="Cancer")



#We create a cancer variable
all_cancer_cmbd<-rbind(head_can, 
                       esopha_can, stomach_can, colorectal_can, liver_can, gall_can,
                       pancreas_can, larynx_can, lung_can, bone_can, mms_can, softtissue_can, 
                       breast_can, cervix_can, uterus_can, ovary_can, penis_can, prostate_can, 
                       testis_can, kidney_can, urinary_can, bladder_can, brain_can, thyroid_can, 
                       hl_can, nhl_can, mmyeloma_can, leukemia_can, others_can) %>% 
  group_by(id) 

n0_canc_cmbd <- nrow(all_cancer_cmbd %>% distinct(id))

all_cancer_cmbd <- all_cancer_cmbd %>% group_by(id) %>% filter(dat == min(dat)) %>%
  ungroup()

n1_canc_cmbd <- nrow(all_cancer_cmbd %>% distinct(id))

all_cancer_cmbd <- all_cancer_cmbd %>%
  group_by(id) %>%
  mutate(n_diagnoses = 1:n())

all_cancer_cmbd <- all_cancer_cmbd %>% group_by(id) %>% filter(n_diagnoses<=1) 

n2_canc <- nrow(all_cancer_cmbd %>% distinct(id))

all_cancer_cmbd <- all_cancer_cmbd%>%rename(date_cancer_cmbd=dat)
all_cancer_cmbd <- all_cancer_cmbd%>%rename(group_cmbd=group)
all_cancer_cmbd <- all_cancer_cmbd%>%rename(general_cmbd=general)

rm (head_can, 
    esopha_can, stomach_can, colorectal_can, liver_can, gall_can,
    pancreas_can, larynx_can, lung_can, bone_can, mms_can, softtissue_can, 
    breast_can, cervix_can, uterus_can, ovary_can, penis_can, prostate_can, 
    testis_can, kidney_can, urinary_can, bladder_can, brain_can, thyroid_can, 
    hl_can, nhl_can, mmyeloma_can, leukemia_can, others_can)

saveRDS(all_cancer_cmbd,"cancer_cmbd.rds")