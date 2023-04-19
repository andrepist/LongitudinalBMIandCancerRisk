library(mice)
library(tidyverse)
library(rms)

#define cancer variables
pop.analysis <- readRDS("E:/Andrea/Trajectories OBECAN/populationCox.rds")
pop.analysis <- pop.analysis %>% ungroup()
meno<-readRDS("E:/Andrea/Trajectories OBECAN/populationMenopause.rds")
#meno<-pop.new
pop.analysis <- pop.analysis %>% select(-group) %>% left_join(meno %>% select(id, "dnaix", "date_menopause", "meno_status_index",    
                                                                              "age.end","meno_status_breastcan",
                                                                              group))
pop.analysis <- pop.analysis %>% pivot_wider(names_from = group, values_from=group,names_glue = "cancer_{group}")
pop.analysis <- pop.analysis %>% mutate_at(vars(starts_with("cancer_")), funs(ifelse(is.na(.),0,1)))
pop.analysis <- pop.analysis %>% mutate(age_gr= ifelse(age.index<=44, "Aged 40 to 44",
                                                       ifelse(age.index<=49, "Aged 45 to 49",
                                                              ifelse(age.index<=54, "Aged 50 to 54",
                                                                     ifelse(age.index<=59, "Aged 55 to 59",
                                                                            ifelse(age.index<=64, "Aged 60 to 64",
                                                                                   ifelse(age.index<=69, "Aged 65 to 69",
                                                                                          ifelse(age.index<=74, "Aged 70 to 74",
                                                                                                 ifelse(age.index<=79, "Aged 75 to 79",
                                                                                                        ifelse(age.index<=84, "Aged 80 to 84",
                                                                                                               ifelse(age.index<=89, "Aged 85 to 89",
                                                                                                                      "Aged 90 or above" )))))))))))

bmiindex <- readRDS("E:/Andrea/Trajectories OBECAN/bmiindex.rds")
pop.analysis <- pop.analysis %>% left_join(bmiindex)

pop.analysis <- pop.analysis %>% select(-c("bmi_18", "bmi_30" ,                             "bmi_40"                            ,  "bmi_55"                    ,         
                                           "bmi_70"        ,                      "bmi_118",                             "n"                                  ,
                                           "bmi_19"      ,                        "bmi_20" ,                             "bmi_21"                             ,
                                           "bmi_22"     ,                         "bmi_23"  ,                            "bmi_24"                             ,
                                           "bmi_25"    ,                          "bmi_26"   ,                           "bmi_27"                             ,
                                           "bmi_28"   ,                           "bmi_29"    ,                          "bmi_31"                             ,
                                           "bmi_32"  ,                            "bmi_33"     ,                         "bmi_34"                             ,
                                           "bmi_35" ,                             "bmi_36"      ,                        "bmi_37"                            , 
                                           "bmi_38",                              "bmi_39" , 
                                           cancer_NA, n))

exposures <- c("yearsBMI25more",
               "yearsBMI30more",
               "ageOnsetOver",
               "ageOnsetObese" ,
               "cumOver" ,
               "cumObese",
               "bmi_atindex"
)
canceres <- c("cancer_Non-Hodgkin Lymphoma" ,
              "cancer_Colorectal" ,
              "cancer_Leukemia"    ,
              "cancer_Liver" ,                      
              "cancer_Brain and CNS"  ,
              "cancer_Pancreas"       ,
              "cancer_Trachea, bronchus & Lung"   , 
              "cancer_Ovary"            ,
              "cancer_Head and Neck"     ,
              "cancer_Others and non-specific"  ,   
              "cancer_Cervix Uteri" ,
              "cancer_Kidney"       ,
              "cancer_Malignant melanoma of skin" , 
              "cancer_Bone and articular cartilage" ,
              "cancer_Testis"        ,
              "cancer_Stomach"     ,                
              "cancer_Thyroid"   ,
              "cancer_Corpus Uteri"    ,
              "cancer_Esophagus"   ,                
              "cancer_Prostate"  ,
              "cancer_Bladder"          ,
              "cancer_Penis"                ,       
              "cancer_Hodgkin lymphoma"    ,
              "cancer_Larynx"                 ,
              "cancer_Connective and soft tissue"  ,
              "cancer_Multiple myeloma"  ,
              "cancer_Gallbladder & biliary tract",
              "cancer_Urinary tract"  ,
              "cancer_Breast pre"  ,
              "cancer_Breast post"  
)


#sensibility - sex -----


out <- list()
harrell<-list()
i = 0
current.pop<-pop.analysis%>%filter(sexe=="D") %>%  group_by(.imp)
for (cancer in canceres){
  print(paste0("Working on cancer", cancer))
  
  # if(cancer =="cancer_Breast post"){
  #   current.pop<-pop.analysis %>% filter(meno_status_index==1) %>% group_by(.imp)
  #   print("breast post")
  # } else {
  #   if(cancer=="cancer_Breast pre"){
  #     current.pop<-pop.analysis %>% filter(meno_status_index==0)%>% group_by(.imp)
  #     print("breast pre")
  #   }
  #   else{print("normal")}
  # }
  
  
  for (exposure in exposures){
    print(paste0("Working on exposure", exposure))  
    
    print(paste0("Working on fully adjusted")) 
    if(cancer %in%c("cancer_Breast post","cancer_Breast pre" )){
      test <- current.pop %>%  
        do(model = coxph(Surv(age.index, age.index + follow_up, get(cancer)) ~
                           get(exposure) + qmedea + nationality + smoking +
                           alcohol +strata(age_gr),  data = .))   
    }else{
      test <- current.pop %>%  
        do(model = coxph(Surv(age.index, age.index + follow_up, get(cancer)) ~
                           get(exposure) + qmedea + nationality + smoking +
                           alcohol +strata(age_gr),  data = .))   
    }
    aic<-sapply(test$model, AIC)
    bic<-sapply(test$model, BIC)
    test <- test%>%
      as.list()
    
    test.pool <- test %>%
      .[[-1]]%>% 
      pool()
    tidy <- tidy(test.pool, conf.int = TRUE, conf.level=0.95)
    i <- i+1
    out[[i]]<-tidy[1,] %>%
      mutate(exposure=exposure, cancer=cancer, group="no",
             model="Fully adjusted", 
             aic1=aic[1], aic2=aic[2], aic3=aic[3], aic4=aic[4],aic5=aic[5],
             bic1=bic[1], bic2=bic[2], bic3=bic[3], bic4=bic[4],bic5=bic[5])
    
    c_index<-tibble()
    for (j in 1:5){
      c_index <- rbind(c_index,(test$model[[j]])$concordance %>% bind_rows()) 
    }
  }
  
  harrell[[i]]<-summarize_if(c_index,is.numeric, mean) %>% mutate(exposure=exposure, cancer=cancer, group="no",
                                                                  model="Male")
  
}

df <- plyr::ldply (out, tibble)
#char<-plyr::ldply (harrell, tibble)

saveRDS(df, "out.females.rds")
#saveRDS(char, "charlson.rds")

#plot ------

#women
df<-readRDS("out.females.rds") 

pop.never <- pop.analysis %>% filter(sexe=="D")
Ncasos <- pop.never %>% select(starts_with("cancer")) %>% summarise(across(where(is.numeric), sum)) %>% 
  mutate_all(funs(round(./5)))

toplotSD <- df %>% 
  #filter(!exposure%in%c("ageOnsetOver", "ageOnsetObese")) %>% 
  # rbind(out1) %>% 
  mutate(expo=recode(exposure,
                     "yearsBMI25more"="Duration overweight",
                     "yearsBMI30more"= "Duration obesity",
                     "cumOver"="Cumulative overweight",
                     "cumObese" ="Cumulative obesity" ,
                     "ageOnsetOver"="Age onset overweight" ,
                     "ageOnsetObese"="Age onset obesity" ,
                     "bmi_atindex"="BMI at index"
                     
  )) %>% mutate(expo=ifelse(model!="Fully adjusted", model, expo))

toplotSD <-toplotSD %>% mutate(
  estimate=ifelse(expo%in%c("Age onset overweight" , "Age onset obesity"), -estimate, estimate),
  conf.low=ifelse(expo%in%c("Age onset overweight" , "Age onset obesity"), -conf.low, conf.low),
  conf.high=ifelse(expo%in%c("Age onset overweight" , "Age onset obesity"), -conf.high, conf.high))
Ncasos <- Ncasos %>% pivot_longer(cols=2:length(Ncasos)) %>% select(-(1))
toplotSD <-toplotSD %>%left_join(Ncasos %>% rename(cancer=name, Nf=value))

toplotSD <-toplotSD %>% separate(cancer, into=c(NA, "cancer"), sep="_") %>% 
  mutate(HR=ifelse(expo=="Duration overweight", exp(estimate*sd(pop.never$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(estimate*sd(pop.never$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(estimate*sd(pop.never$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(estimate*sd(pop.never$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(estimate*sd(pop.never$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(estimate*sd(pop.never$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(estimate*sd(pop.never$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(estimate*sd(pop.never$bmi_atindex)),
                                                                    "Missing")))))))),
         low=ifelse(expo=="Duration overweight", exp(conf.low*sd(pop.never$yearsBMI25more)),
                    ifelse(expo=="Duration obesity", exp(conf.low*sd(pop.never$yearsBMI30more)),
                           ifelse(expo=="Cumulative overweight", exp(conf.low*sd(pop.never$cumOver)), 
                                  ifelse(expo=="Cumulative obesity", exp(conf.low*sd(pop.never$cumObese)), 
                                         ifelse(expo=="Age onset overweight", exp(conf.low*sd(pop.never$ageOnsetOver, na.rm = T)), 
                                                ifelse(expo=="Age onset obesity", exp(conf.low*sd(pop.never$ageOnsetObese, na.rm = T)), 
                                                       ifelse(expo=="Yearly BMI change", exp(conf.low*sd(pop.never$trajectorySlope)),
                                                              ifelse(expo=="BMI at index", exp(conf.low*sd(pop.never$bmi_atindex)),"Missing")))))))),
         up=ifelse(expo=="Duration overweight", exp(conf.high*sd(pop.never$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(conf.high*sd(pop.never$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(conf.high*sd(pop.never$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(conf.high*sd(pop.never$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(conf.high*sd(pop.never$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(conf.high*sd(pop.never$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(conf.high*sd(pop.never$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(conf.high*sd(pop.never$bmi_atindex)),
                                                                    "Missing")))))))))%>%  
  mutate(expo2=ifelse(exposure%in%c("yearsBMI25more", "yearsBMI30more"), "Duration",
                      ifelse(exposure%in%c( "cumOver","cumObese" ), "Cumulative",
                             ifelse(exposure%in%c( "ageOnsetOver","ageOnsetObese" ), "Age of onset",
                                    "Slope"))),
         Type=ifelse(exposure%in%c("yearsBMI25more", "cumOver","ageOnsetOver"), "Overweight/Obese",
                     ifelse(exposure%in%c("yearsBMI30more", "cumObese","ageOnsetObese"), "Obese",
                            "Overall")))

# to test
women <-toplotSD %>% filter(cancer%in%c( "Cervix Uteri","Corpus Uteri","Ovary")) %>% 
  mutate(HR=ifelse(expo=="Duration overweight", exp(estimate*sd((pop.never%>%filter(sexe=="D"))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(estimate*sd((pop.never%>%filter(sexe=="D"))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(estimate*sd((pop.never%>%filter(sexe=="D"))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(estimate*sd((pop.never%>%filter(sexe=="D"))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(estimate*sd((pop.never%>%filter(sexe=="D"))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(estimate*sd((pop.never%>%filter(sexe=="D"))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(estimate*sd((pop.never%>%filter(sexe=="D"))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(estimate*sd((pop.never%>%filter(sexe=="D"))$bmi_atindex)),
                                                                    "Missing")))))))),
         low=ifelse(expo=="Duration overweight", exp(conf.low*sd((pop.never%>%filter(sexe=="D"))$yearsBMI25more)),
                    ifelse(expo=="Duration obesity", exp(conf.low*sd((pop.never%>%filter(sexe=="D"))$yearsBMI30more)),
                           ifelse(expo=="Cumulative overweight", exp(conf.low*sd((pop.never%>%filter(sexe=="D"))$cumOver)), 
                                  ifelse(expo=="Cumulative obesity", exp(conf.low*sd((pop.never%>%filter(sexe=="D"))$cumObese)), 
                                         ifelse(expo=="Age onset overweight", exp(conf.low*sd((pop.never%>%filter(sexe=="D"))$ageOnsetOver, na.rm = T)), 
                                                ifelse(expo=="Age onset obesity", exp(conf.low*sd((pop.never%>%filter(sexe=="D"))$ageOnsetObese, na.rm = T)), 
                                                       ifelse(expo=="Yearly BMI change", exp(conf.low*sd((pop.never%>%filter(sexe=="D"))$trajectorySlope)),
                                                              ifelse(expo=="BMI at index", exp(conf.low*sd((pop.never%>%filter(sexe=="D"))$bmi_atindex)),"Missing")))))))),
         up=ifelse(expo=="Duration overweight", exp(conf.high*sd((pop.never%>%filter(sexe=="D"))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(conf.high*sd((pop.never%>%filter(sexe=="D"))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(conf.high*sd((pop.never%>%filter(sexe=="D"))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(conf.high*sd((pop.never%>%filter(sexe=="D"))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(conf.high*sd((pop.never%>%filter(sexe=="D"))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(conf.high*sd((pop.never%>%filter(sexe=="D"))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(conf.high*sd((pop.never%>%filter(sexe=="D"))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(conf.high*sd((pop.never%>%filter(sexe=="D"))$bmi_atindex)),
                                                                    "Missing")))))))))

men <- toplotSD %>% filter(cancer%in%c("Prostate", "Testis")) %>% 
  mutate(HR=ifelse(expo=="Duration overweight", exp(estimate*sd((pop.never%>%filter(sexe=="H"))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(estimate*sd((pop.never%>%filter(sexe=="H"))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(estimate*sd((pop.never%>%filter(sexe=="H"))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(estimate*sd((pop.never%>%filter(sexe=="H"))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(estimate*sd((pop.never%>%filter(sexe=="H"))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(estimate*sd((pop.never%>%filter(sexe=="H"))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(estimate*sd((pop.never%>%filter(sexe=="H"))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(estimate*sd((pop.never%>%filter(sexe=="H"))$bmi_atindex)),
                                                                    "Missing")))))))),
         low=ifelse(expo=="Duration overweight", exp(conf.low*sd((pop.never%>%filter(sexe=="H"))$yearsBMI25more)),
                    ifelse(expo=="Duration obesity", exp(conf.low*sd((pop.never%>%filter(sexe=="H"))$yearsBMI30more)),
                           ifelse(expo=="Cumulative overweight", exp(conf.low*sd((pop.never%>%filter(sexe=="H"))$cumOver)), 
                                  ifelse(expo=="Cumulative obesity", exp(conf.low*sd((pop.never%>%filter(sexe=="H"))$cumObese)), 
                                         ifelse(expo=="Age onset overweight", exp(conf.low*sd((pop.never%>%filter(sexe=="H"))$ageOnsetOver, na.rm = T)), 
                                                ifelse(expo=="Age onset obesity", exp(conf.low*sd((pop.never%>%filter(sexe=="H"))$ageOnsetObese, na.rm = T)), 
                                                       ifelse(expo=="Yearly BMI change", exp(conf.low*sd((pop.never%>%filter(sexe=="H"))$trajectorySlope)),
                                                              ifelse(expo=="BMI at index", exp(conf.low*sd((pop.never%>%filter(sexe=="H"))$bmi_atindex)),"Missing")))))))),
         up=ifelse(expo=="Duration overweight", exp(conf.high*sd((pop.never%>%filter(sexe=="H"))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(conf.high*sd((pop.never%>%filter(sexe=="H"))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(conf.high*sd((pop.never%>%filter(sexe=="H"))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(conf.high*sd((pop.never%>%filter(sexe=="H"))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(conf.high*sd((pop.never%>%filter(sexe=="H"))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(conf.high*sd((pop.never%>%filter(sexe=="H"))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(conf.high*sd((pop.never%>%filter(sexe=="H"))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(conf.high*sd((pop.never%>%filter(sexe=="H"))$bmi_atindex)),
                                                                    "Missing")))))))))

breast_pre <-toplotSD %>% filter(cancer%in%c( "Breast pre")) %>% 
  mutate(HR=ifelse(expo=="Duration overweight", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$bmi_atindex)),
                                                                    "Missing")))))))),
         low=ifelse(expo=="Duration overweight", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$yearsBMI25more)),
                    ifelse(expo=="Duration obesity", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$yearsBMI30more)),
                           ifelse(expo=="Cumulative overweight", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$cumOver)), 
                                  ifelse(expo=="Cumulative obesity", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$cumObese)), 
                                         ifelse(expo=="Age onset overweight", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$ageOnsetOver, na.rm = T)), 
                                                ifelse(expo=="Age onset obesity", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$ageOnsetObese, na.rm = T)), 
                                                       ifelse(expo=="Yearly BMI change", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$trajectorySlope)),
                                                              ifelse(expo=="BMI at index", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$bmi_atindex)),"Missing")))))))),
         up=ifelse(expo=="Duration overweight", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$bmi_atindex)),
                                                                    "Missing")))))))))

breast_post <-toplotSD %>% filter(cancer%in%c( "Breast post")) %>% 
  mutate(HR=ifelse(expo=="Duration overweight", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$bmi_atindex)),
                                                                    "Missing")))))))),
         low=ifelse(expo=="Duration overweight", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$yearsBMI25more)),
                    ifelse(expo=="Duration obesity", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$yearsBMI30more)),
                           ifelse(expo=="Cumulative overweight", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$cumOver)), 
                                  ifelse(expo=="Cumulative obesity", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$cumObese)), 
                                         ifelse(expo=="Age onset overweight", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$ageOnsetOver, na.rm = T)), 
                                                ifelse(expo=="Age onset obesity", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$ageOnsetObese, na.rm = T)), 
                                                       ifelse(expo=="Yearly BMI change", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$trajectorySlope)),
                                                              ifelse(expo=="BMI at index", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$bmi_atindex)),"Missing")))))))),
         up=ifelse(expo=="Duration overweight", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$bmi_atindex)),
                                                                    "Missing")))))))))

#toplot <- toplot %>% mutate(Type=factor(Type, levels=rev(c("Overweight/Obese", "Overall", "Obese"))))
woman <- toplotSD %>% filter(!cancer%in%c( "Cervix Uteri","Corpus Uteri","Ovary","Prostate", "Testis", "Breast pre", "Breast post")) %>% 
  rbind(women) %>% 
  rbind(men) %>% 
  rbind(breast_pre) %>% 
  rbind(breast_post)

#men
df<-readRDS("out.males.rds") 

pop.never <- pop.analysis %>% filter(sexe=="H")
Ncasos <- pop.never %>% select(starts_with("cancer")) %>% summarise(across(where(is.numeric), sum)) %>% 
  mutate_all(funs(round(./5)))

toplotSD <- df %>% 
  #filter(!exposure%in%c("ageOnsetOver", "ageOnsetObese")) %>% 
  # rbind(out1) %>% 
  mutate(expo=recode(exposure,
                     "yearsBMI25more"="Duration overweight",
                     "yearsBMI30more"= "Duration obesity",
                     "cumOver"="Cumulative overweight",
                     "cumObese" ="Cumulative obesity" ,
                     "ageOnsetOver"="Age onset overweight" ,
                     "ageOnsetObese"="Age onset obesity" ,
                     "bmi_atindex"="BMI at index"
                     
  )) %>% mutate(expo=ifelse(model!="Fully adjusted", model, expo))

toplotSD <-toplotSD %>% mutate(
  estimate=ifelse(expo%in%c("Age onset overweight" , "Age onset obesity"), -estimate, estimate),
  conf.low=ifelse(expo%in%c("Age onset overweight" , "Age onset obesity"), -conf.low, conf.low),
  conf.high=ifelse(expo%in%c("Age onset overweight" , "Age onset obesity"), -conf.high, conf.high))

Ncasos <- Ncasos %>% pivot_longer(cols=2:length(Ncasos)) %>% select(-(1))
toplotSD <-toplotSD %>%left_join(Ncasos %>% rename(cancer=name, Nm=value))

toplotSD <-toplotSD %>% separate(cancer, into=c(NA, "cancer"), sep="_") %>% 
  mutate(HR=ifelse(expo=="Duration overweight", exp(estimate*sd(pop.never$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(estimate*sd(pop.never$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(estimate*sd(pop.never$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(estimate*sd(pop.never$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(estimate*sd(pop.never$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(estimate*sd(pop.never$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(estimate*sd(pop.never$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(estimate*sd(pop.never$bmi_atindex)),
                                                                    "Missing")))))))),
         low=ifelse(expo=="Duration overweight", exp(conf.low*sd(pop.never$yearsBMI25more)),
                    ifelse(expo=="Duration obesity", exp(conf.low*sd(pop.never$yearsBMI30more)),
                           ifelse(expo=="Cumulative overweight", exp(conf.low*sd(pop.never$cumOver)), 
                                  ifelse(expo=="Cumulative obesity", exp(conf.low*sd(pop.never$cumObese)), 
                                         ifelse(expo=="Age onset overweight", exp(conf.low*sd(pop.never$ageOnsetOver, na.rm = T)), 
                                                ifelse(expo=="Age onset obesity", exp(conf.low*sd(pop.never$ageOnsetObese, na.rm = T)), 
                                                       ifelse(expo=="Yearly BMI change", exp(conf.low*sd(pop.never$trajectorySlope)),
                                                              ifelse(expo=="BMI at index", exp(conf.low*sd(pop.never$bmi_atindex)),"Missing")))))))),
         up=ifelse(expo=="Duration overweight", exp(conf.high*sd(pop.never$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(conf.high*sd(pop.never$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(conf.high*sd(pop.never$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(conf.high*sd(pop.never$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(conf.high*sd(pop.never$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(conf.high*sd(pop.never$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(conf.high*sd(pop.never$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(conf.high*sd(pop.never$bmi_atindex)),
                                                                    "Missing")))))))))%>%  
  mutate(expo2=ifelse(exposure%in%c("yearsBMI25more", "yearsBMI30more"), "Duration",
                      ifelse(exposure%in%c( "cumOver","cumObese" ), "Cumulative",
                             ifelse(exposure%in%c( "ageOnsetOver","ageOnsetObese" ), "Age of onset",
                                    "Slope"))),
         Type=ifelse(exposure%in%c("yearsBMI25more", "cumOver","ageOnsetOver"), "Overweight/Obese",
                     ifelse(exposure%in%c("yearsBMI30more", "cumObese","ageOnsetObese"), "Obese",
                            "Overall")))

# to test
women <-toplotSD %>% filter(cancer%in%c( "Cervix Uteri","Corpus Uteri","Ovary")) %>% 
  mutate(HR=ifelse(expo=="Duration overweight", exp(estimate*sd((pop.never%>%filter(sexe=="D"))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(estimate*sd((pop.never%>%filter(sexe=="D"))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(estimate*sd((pop.never%>%filter(sexe=="D"))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(estimate*sd((pop.never%>%filter(sexe=="D"))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(estimate*sd((pop.never%>%filter(sexe=="D"))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(estimate*sd((pop.never%>%filter(sexe=="D"))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(estimate*sd((pop.never%>%filter(sexe=="D"))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(estimate*sd((pop.never%>%filter(sexe=="D"))$bmi_atindex)),
                                                                    "Missing")))))))),
         low=ifelse(expo=="Duration overweight", exp(conf.low*sd((pop.never%>%filter(sexe=="D"))$yearsBMI25more)),
                    ifelse(expo=="Duration obesity", exp(conf.low*sd((pop.never%>%filter(sexe=="D"))$yearsBMI30more)),
                           ifelse(expo=="Cumulative overweight", exp(conf.low*sd((pop.never%>%filter(sexe=="D"))$cumOver)), 
                                  ifelse(expo=="Cumulative obesity", exp(conf.low*sd((pop.never%>%filter(sexe=="D"))$cumObese)), 
                                         ifelse(expo=="Age onset overweight", exp(conf.low*sd((pop.never%>%filter(sexe=="D"))$ageOnsetOver, na.rm = T)), 
                                                ifelse(expo=="Age onset obesity", exp(conf.low*sd((pop.never%>%filter(sexe=="D"))$ageOnsetObese, na.rm = T)), 
                                                       ifelse(expo=="Yearly BMI change", exp(conf.low*sd((pop.never%>%filter(sexe=="D"))$trajectorySlope)),
                                                              ifelse(expo=="BMI at index", exp(conf.low*sd((pop.never%>%filter(sexe=="D"))$bmi_atindex)),"Missing")))))))),
         up=ifelse(expo=="Duration overweight", exp(conf.high*sd((pop.never%>%filter(sexe=="D"))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(conf.high*sd((pop.never%>%filter(sexe=="D"))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(conf.high*sd((pop.never%>%filter(sexe=="D"))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(conf.high*sd((pop.never%>%filter(sexe=="D"))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(conf.high*sd((pop.never%>%filter(sexe=="D"))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(conf.high*sd((pop.never%>%filter(sexe=="D"))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(conf.high*sd((pop.never%>%filter(sexe=="D"))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(conf.high*sd((pop.never%>%filter(sexe=="D"))$bmi_atindex)),
                                                                    "Missing")))))))))

men <- toplotSD %>% filter(cancer%in%c("Prostate", "Testis")) %>% 
  mutate(HR=ifelse(expo=="Duration overweight", exp(estimate*sd((pop.never%>%filter(sexe=="H"))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(estimate*sd((pop.never%>%filter(sexe=="H"))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(estimate*sd((pop.never%>%filter(sexe=="H"))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(estimate*sd((pop.never%>%filter(sexe=="H"))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(estimate*sd((pop.never%>%filter(sexe=="H"))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(estimate*sd((pop.never%>%filter(sexe=="H"))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(estimate*sd((pop.never%>%filter(sexe=="H"))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(estimate*sd((pop.never%>%filter(sexe=="H"))$bmi_atindex)),
                                                                    "Missing")))))))),
         low=ifelse(expo=="Duration overweight", exp(conf.low*sd((pop.never%>%filter(sexe=="H"))$yearsBMI25more)),
                    ifelse(expo=="Duration obesity", exp(conf.low*sd((pop.never%>%filter(sexe=="H"))$yearsBMI30more)),
                           ifelse(expo=="Cumulative overweight", exp(conf.low*sd((pop.never%>%filter(sexe=="H"))$cumOver)), 
                                  ifelse(expo=="Cumulative obesity", exp(conf.low*sd((pop.never%>%filter(sexe=="H"))$cumObese)), 
                                         ifelse(expo=="Age onset overweight", exp(conf.low*sd((pop.never%>%filter(sexe=="H"))$ageOnsetOver, na.rm = T)), 
                                                ifelse(expo=="Age onset obesity", exp(conf.low*sd((pop.never%>%filter(sexe=="H"))$ageOnsetObese, na.rm = T)), 
                                                       ifelse(expo=="Yearly BMI change", exp(conf.low*sd((pop.never%>%filter(sexe=="H"))$trajectorySlope)),
                                                              ifelse(expo=="BMI at index", exp(conf.low*sd((pop.never%>%filter(sexe=="H"))$bmi_atindex)),"Missing")))))))),
         up=ifelse(expo=="Duration overweight", exp(conf.high*sd((pop.never%>%filter(sexe=="H"))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(conf.high*sd((pop.never%>%filter(sexe=="H"))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(conf.high*sd((pop.never%>%filter(sexe=="H"))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(conf.high*sd((pop.never%>%filter(sexe=="H"))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(conf.high*sd((pop.never%>%filter(sexe=="H"))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(conf.high*sd((pop.never%>%filter(sexe=="H"))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(conf.high*sd((pop.never%>%filter(sexe=="H"))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(conf.high*sd((pop.never%>%filter(sexe=="H"))$bmi_atindex)),
                                                                    "Missing")))))))))

breast_pre <-toplotSD %>% filter(cancer%in%c( "Breast pre")) %>% 
  mutate(HR=ifelse(expo=="Duration overweight", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$bmi_atindex)),
                                                                    "Missing")))))))),
         low=ifelse(expo=="Duration overweight", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$yearsBMI25more)),
                    ifelse(expo=="Duration obesity", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$yearsBMI30more)),
                           ifelse(expo=="Cumulative overweight", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$cumOver)), 
                                  ifelse(expo=="Cumulative obesity", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$cumObese)), 
                                         ifelse(expo=="Age onset overweight", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$ageOnsetOver, na.rm = T)), 
                                                ifelse(expo=="Age onset obesity", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$ageOnsetObese, na.rm = T)), 
                                                       ifelse(expo=="Yearly BMI change", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$trajectorySlope)),
                                                              ifelse(expo=="BMI at index", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$bmi_atindex)),"Missing")))))))),
         up=ifelse(expo=="Duration overweight", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==0))$bmi_atindex)),
                                                                    "Missing")))))))))

breast_post <-toplotSD %>% filter(cancer%in%c( "Breast post")) %>% 
  mutate(HR=ifelse(expo=="Duration overweight", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(estimate*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$bmi_atindex)),
                                                                    "Missing")))))))),
         low=ifelse(expo=="Duration overweight", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$yearsBMI25more)),
                    ifelse(expo=="Duration obesity", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$yearsBMI30more)),
                           ifelse(expo=="Cumulative overweight", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$cumOver)), 
                                  ifelse(expo=="Cumulative obesity", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$cumObese)), 
                                         ifelse(expo=="Age onset overweight", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$ageOnsetOver, na.rm = T)), 
                                                ifelse(expo=="Age onset obesity", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$ageOnsetObese, na.rm = T)), 
                                                       ifelse(expo=="Yearly BMI change", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$trajectorySlope)),
                                                              ifelse(expo=="BMI at index", exp(conf.low*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$bmi_atindex)),"Missing")))))))),
         up=ifelse(expo=="Duration overweight", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(conf.high*sd((pop.never%>%filter(sexe=="D",meno_status_index==1))$bmi_atindex)),
                                                                    "Missing")))))))))

#toplot <- toplot %>% mutate(Type=factor(Type, levels=rev(c("Overweight/Obese", "Overall", "Obese"))))
man <- toplotSD %>% filter(!cancer%in%c( "Cervix Uteri","Corpus Uteri","Ovary","Prostate", "Testis", "Breast pre", "Breast post")) %>% 
  rbind(women) %>% 
  rbind(men) %>% 
  rbind(breast_pre) %>% 
  rbind(breast_post)

toplotSD<- rbind(man %>% mutate(strata="M") %>% rename(N=Nm), 
                 woman %>% mutate(strata="F")%>% rename(N=Nf))

order.cancer <- c("Head and Neck",
                  "Esophagus",
                  "Stomach",
                  "Colorectal",
                  "Liver",
                  "Gallbladder & \nbiliary tract",
                  "Pancreas",
                  "Larynx",
                  "Trachea, \nbronchus & Lung",
                  "Bone and \narticular cartilage",
                  "Malignant melanoma \nof skin",
                  #"Connective and soft tissue",
                  "Breast \npremenopausal",
                  "Breast \npostmenopausal",
                  "Breast",
                  "Cervix Uteri",
                  "Corpus Uteri",
                  "Ovary",
                  #"Penis",
                  "Prostate",
                  "Testis",
                  "Kidney",
                  #"Urinary Tract",
                  "Bladder",
                  "Brain and CNS",
                  "Thyroid",
                  "Hodgkin lymphoma",
                  "Non-Hodgkin Lymph.",
                  "Multiple myeloma",
                  "Leukemia"
)

toplotSD <- toplotSD %>% mutate(cancer=recode(cancer, 
                                              "Gallbladder & biliary tract"= "Gallbladder & \nbiliary tract",
                                              "Trachea, bronchus & Lung" ="Trachea, \nbronchus & Lung" ,
                                              "Bone and articular cartilage"="Bone and \narticular cartilage",
                                              "Malignant melanoma of skin" ="Malignant melanoma \nof skin",
                                              "Breast pre"="Breast \npremenopausal",
                                              "Breast post"="Breast \npostmenopausal",
                                              "Non-Hodgkin Lymphoma"="Non-Hodgkin Lymph."
)) %>% 
  filter(cancer%in%order.cancer) 

toplotSD <- toplotSD %>%mutate(HR1=as.numeric(HR),
                               low=as.numeric(low),
                               up=as.numeric(up))
toplotSD <- toplotSD %>% mutate(low1=ifelse(low>up, up, low),
                                up1=ifelse(up<low, low, up)) %>% 
  mutate(low=low1, up=up1)
toplotSD <- toplotSD %>%
  mutate(expo=factor(expo, levels=rev(c("BMI at index","Duration overweight", "Duration obesity",
                                        "Cumulative overweight", "Cumulative obesity",
                                        "Age onset overweight" , "Age onset obesity"
  ))),
  Type=recode(Type, "Obese"="Obesity",
              "Overweight/Obese"= "Overweight/Obesity"))

toplotSD <- toplotSD %>% mutate(expo3=ifelse(expo=="BMI at index", "BMI at index", expo2))

toplotSD <- toplotSD %>% mutate(expo=recode(expo,"BMI at index"="BMI at index date",
                           "Duration overweight"= "Duration BMI \u226525",
                           "Duration obesity"= "Duration BMI \u226530",
                           "Cumulative overweight"= "Cum. exposure BMI \u226525",
                           "Cumulative obesity"=  "Cum. exposure BMI \u226530",
                            "Age onset overweight"="Age of onset BMI \u226525",
                           "Age onset obesity"= "Age of onset BMI \u226530"
))

order.cancer1<-readRDS("order.cancer.rds")
#View((toplotSD%>%filter(group=="no", Type!="Overall", expo=="BMI at index")%>%arrange(HR1) ))
#saveRDS(order.cancer1, "order.cancer.rds")
toplotSD<-toplotSD %>% 
  mutate(cancer=factor(cancer, levels=rev(order.cancer1)))

toplotSD<-toplotSD %>% 
  mutate(expostrata=paste0(expo, " (", strata, ")"),
         expostrata= factor(expostrata, levels=rev(c("BMI at index date (M)","BMI at index date (F)",
                                         "Duration BMI \u226525 (M)","Duration BMI \u226525 (F)", 
                                         "Duration BMI \u226530 (M)", "Duration BMI \u226530 (F)",
                                   "Cum. exposure BMI \u226525 (M)", "Cum. exposure BMI \u226525 (F)",
                                   "Cum. exposure BMI \u226530 (M)", "Cum. exposure BMI \u226530 (F)",
                                   "Age of onset BMI \u226525 (M)" ,"Age of onset BMI \u226525 (F)" ,
                                   "Age of onset BMI \u226530 (M)", "Age of onset BMI \u226530 (F)"
         ))) )

Nm <- toplotSD %>% filter(strata=="M") %>% rename(Nm=N) %>% select(Nm,cancer)
Nf <- toplotSD %>% filter(strata=="F") %>% rename(Nf=N)%>% select(Nf, cancer)

toplotSD <- toplotSD %>% left_join(Nm) %>% left_join(Nf) %>% distinct()

toplotSD <- toplotSD %>% 
  mutate(label=trimws(paste0(cancer, "\n(N M = ", trimws(format(Nm, big.mark=",")),"; N F = " ,
                             trimws(format(Nf, big.mark=",")),")")))

ggplot(toplotSD%>% filter(group=="no",  expo!="Missing") )+
  geom_errorbar(aes(xmin=low, xmax=up, y=expostrata, group=expo),  width = 0.5, colour="grey30")+
  geom_point(aes(HR1, expostrata, col=expo, shape=expo3), 
             size=1.8)+
  facet_wrap(cancer~., ncol=5,labeller=as_labeller(setNames(toplotSD$label,toplotSD$cancer)))+
  geom_vline(xintercept=1, linetype="dashed", colour="grey20")+
  theme_bw() +
  theme(panel.background = element_blank(),
        # legend.title = element_blank(),
        legend.position = "null",
        #panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        # strip.text.y.left = element_text(angle = 0),
        axis.title.y.left = element_blank(),
        #axis.ticks.y= element_blank(),
        #axis.text.y.left= element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.background = element_rect( fill="grey92"))+ 
  
  labs(colour="Legend")+
  # scale_y_discrete(labels = rev(c("BMI at index date",
  #                                 expression("Duration BMI">="25"),
  #                                 expression("Duration BMI">="30"),
  #                                 expression("Cum. exposure BMI">="25"),
  #                                 expression("Cum. exposure BMI">="30"),
  #                                 expression("Age onset BMI">="25"),
  #                                 expression("Age onset BMI">="30"))))+
  scale_color_manual(values=rev(c("black","#03992A", "#FF9F2F", "#32A5DD","#FF5B2F","#C51A95","#AA6C39")))+
  scale_shape_manual(values=c(15,16,17,8))+
  xlab("HR (95% CI)")+
  coord_cartesian(xlim=c(0.7, 1.6))

table<-toplotSD%>% filter(group=="no", Type!="Overall", expo!="Missing")
htmlTable::htmlTable(table %>% select(strata, cancer, expo, HR1,low, up, N))

ggsave("S7B.png",
       dpi=300,
       width = 10, height = 12)

# table  -----------
table<-toplotSD%>% filter(group=="no",  expo!="Missing")
table <- table %>% select(HR1, cancer, expostrata, low, up)
table <- table %>% 
  mutate(HR=paste0(format(round(HR1,2),nsmall=2),
                   " (", 
                   format(round(low, 2),nsmall=2), "-", 
                   format(round(up,2),nsmall=2), ")")) %>% 
  select(-HR1, -low,-up)



table <- table %>% pivot_wider(names_from = "expostrata", values_from="HR")
table[match(rev(order.cancer1), table$cancer),] %>% htmlTable::htmlTable()

# table casos --------
Ncasos <- pop.analysis %>% filter(sexe=="H",everObese==1) %>% 
  select(starts_with("cancer")) %>% summarise(across(where(is.numeric), sum)) %>% 
  mutate_all(funs(round(./5)))

Ncasos <- Ncasos %>% pivot_longer(cols=2:length(Ncasos)) %>% select(-(1)) 
Ncasos <- Ncasos%>%
  separate(name, into=c(NA, "cancer"), sep="_")

Ncasos <- Ncasos %>% mutate(cancer=recode(cancer, 
                                          "Gallbladder & biliary tract"= "Gallbladder & \nbiliary tract",
                                          "Trachea, bronchus & Lung" ="Trachea, \nbronchus & Lung" ,
                                          "Bone and articular cartilage"="Bone and \narticular cartilage",
                                          "Malignant melanoma of skin" ="Malignant melanoma \nof skin",
                                          "Breast pre"="Breast \npremenopausal",
                                          "Breast post"="Breast \npostmenopausal",
                                          "Non-Hodgkin Lymphoma"="Non-Hodgkin Lymph."
)) 
table <- Ncasos %>% mutate(value=format(value, big.mark = ","))

table[match(rev(order.cancer1), table$cancer),] %>% htmlTable::htmlTable()


