#analysis
library(mice)
library(tidyverse)
library(rms)
library(lubridate)

check <- Ncasos %>% select(-c(cancer.outcome, 
                              "cancer_Others and non-specific"   ,
                              "cancer_Penis",
                              "cancer_Connective and soft tissue" ,
                              "cancer_Urinary tract"  ))
sum(check)
#define cancer variables
pop.analysis <- readRDS("E:/Andrea/Trajectories OBECAN/populationCox.rds")
pop.analysis <- pop.analysis %>% ungroup()
#pop.analysis<-pop
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
bmiindex <- bmiindex %>% ungroup()
pop.analysis <- pop.analysis %>% left_join(bmiindex)
# linear ----
exposures <- c("yearsBMI25more",
               "yearsBMI30more",
               "ageOnsetOver",
               "ageOnsetObese" ,
               "cumOver" ,
               "cumObese" 
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


exposures<-c("yearsBMI25more",
             #"yearsBMI30more",
             "ageOnsetOver",
             #"ageOnsetObese" ,
             "cumOver")
canceres <- c("cancer_Corpus Uteri")

pop.analysis<-pop.analysis %>% mutate(bmi.gr18=ifelse(bmi_18>=25&bmi_18<30, "Overweight", 
                                                      ifelse(bmi_18>=30, "Obese", "Normal")))

#change now its restricted to overweight
out <- list()
harrell<-list()
i = 0
current.pop<-pop.analysis%>% group_by(.imp)
for (cancer in canceres){
  print(paste0("Working on cancer", cancer))
  
  if(cancer =="cancer_Breast post"){
    current.pop<-pop.analysis %>% filter(meno_status_index==1) %>% group_by(.imp)
    print("breast post")
  } else {
    if(cancer=="cancer_Breast pre"){
      current.pop<-pop.analysis %>% filter(meno_status_index==0)%>% group_by(.imp)
      print("breast pre")
    }
    else{print("normal")}
  }
  
  
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
                           get(exposure) +sexe+ qmedea + nationality + smoking +
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
    harrell[[i]]<-summarize_if(c_index,is.numeric, mean) %>% mutate(exposure=exposure, cancer=cancer, group="no",
                                                                    model="Fully adjusted")
    
    
    # print(paste0("Working on stratify the slope")) 
    # for (bmigr in c("Normal","Overweight", "Obese" )){
    #   
    #   if(cancer %in%c("cancer_Breast post","cancer_Breast pre" )){
    #     test <- current.pop %>% 
    #       do(model = coxph(Surv(age.index, age.index + follow_up, get(cancer)) ~
    #                          get(exposure)  + qmedea + nationality + smoking +
    #                          alcohol+strata(age_gr), subset=(bmi.gr18==bmigr), data = .))  
    #   }else{
    #     test <- current.pop %>% 
    #       do(model = coxph(Surv(age.index, age.index + follow_up, get(cancer)) ~
    #                          get(exposure) +sexe + qmedea + nationality + smoking +
    #                          alcohol+strata(age_gr), subset=(bmi.gr18==bmigr), data = .))   
    #   }
    #   
    #   
    #   aic<-sapply(test$model, AIC)
    #   bic<-sapply(test$model, BIC)
    #   test <- test%>%
    #     as.list()
    #   test.pool <- test %>%
    #     .[[-1]]%>%
    #     pool()
    #   
    #   tidy <- tidy(test.pool, conf.int = TRUE, conf.level=0.95)
    #   i <- i+1
    #   out[[i]]<-tidy[1,] %>% 
    #     mutate(exposure=exposure, cancer=cancer, model="Slope stratified",group=bmigr, 
    #            aic1=aic[1], aic2=aic[2], aic3=aic[3], aic4=aic[4],aic5=aic[5],
    #            bic1=bic[1], bic2=bic[2], bic3=bic[3], bic4=bic[4],bic5=bic[5])
    # }
    
  }
  
  print(paste0("Working on exposure BMI at index")) 
  
  if(cancer %in%c("cancer_Breast post","cancer_Breast pre" )){
    test <- current.pop %>%  
      do(model = coxph(Surv(age.index, age.index + follow_up, get(cancer)) ~
                         bmi_atindex + qmedea + nationality + smoking +
                         alcohol +strata(age_gr),  data = .))   
  }else{
    test <- current.pop %>%  
      do(model = coxph(Surv(age.index, age.index + follow_up, get(cancer)) ~
                         bmi_atindex +sexe + qmedea + nationality + smoking +
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
    mutate(exposure=exposure, cancer=cancer, group="no", model="BMI at index", 
           aic1=aic[1], aic2=aic[2], aic3=aic[3], aic4=aic[4],aic5=aic[5],
           bic1=bic[1], bic2=bic[2], bic3=bic[3], bic4=bic[4],bic5=bic[5])
  
  c_index<-tibble()
  for (j in 1:5){
    c_index <- rbind(c_index,(test$model[[j]])$concordance %>% bind_rows()) 
  }
  harrell[[i]]<-summarize_if(c_index,is.numeric, mean) %>%
    mutate(exposure=exposure, cancer=cancer, group="no", model="BMI at index")
  
  print(paste0("Working on null model")) 
  
  if(cancer %in%c("cancer_Breast post","cancer_Breast pre" )){
    test <- current.pop %>%  
      do(model = coxph(Surv(age.index, age.index + follow_up, get(cancer)) ~
                         qmedea + nationality + smoking +
                         alcohol +strata(age_gr),  data = .))   
  }else{
    
    test <- current.pop %>%  
      do(model = coxph(Surv(age.index, age.index + follow_up, get(cancer)) ~
                         sexe + qmedea + nationality + smoking +
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
    mutate(exposure=exposure, cancer=cancer, group="no", model="Null model", 
           aic1=aic[1], aic2=aic[2], aic3=aic[3], aic4=aic[4],aic5=aic[5],
           bic1=bic[1], bic2=bic[2], bic3=bic[3], bic4=bic[4],bic5=bic[5])
  c_index<-tibble()
  for (j in 1:5){
    c_index <- rbind(c_index,(test$model[[j]])$concordance %>% bind_rows()) 
  }
  harrell[[i]]<-summarize_if(c_index,is.numeric, mean) %>%
    mutate(exposure=exposure, cancer=cancer, group="no", model="Null model")
  
  
}

df <- plyr::ldply (out, tibble)
char<-plyr::ldply (harrell, tibble)

saveRDS(df, "out.951Y.rds")
saveRDS(char, "charlson.rds")
# plot ----
out <- df
toplot <- out %>% 
  mutate(expo=ifelse(exposure%in%c("yearsBMI25more", "yearsBMI30more"), "Duration",
                     ifelse(exposure%in%c( "cumOver","cumObese" ), "Cumulative",
                            ifelse(exposure%in%c( "ageOnsetOver","ageOnsetObese" ), "Age of onset",
                                   "Slope"))),
         Type=ifelse(exposure%in%c("yearsBMI25more", "cumOver","ageOnsetOver"), "Overweight/Obese",
                     ifelse(exposure%in%c("yearsBMI30more", "cumObese","ageOnsetObese"), "Obese",
                            "Overall")))

toplot <- toplot %>% mutate(HR=ifelse(expo=="Duration", exp(estimate*10),
                                      ifelse(expo=="Cumulative", exp(estimate*100), 
                                             ifelse(expo=="Slope", exp(estimate*0.1),
                                                    ifelse(expo=="Age of onset", exp(estimate*1),exp(estimate))))),
                            low=ifelse(expo=="Duration", exp(conf.low*10),
                                       ifelse(expo=="Cumulative", exp(conf.low*100), 
                                              ifelse(expo=="Slope", exp(conf.low*0.1),
                                                     ifelse(expo=="Age of onset", exp(conf.low*1), exp(conf.low))))),
                            up=ifelse(expo=="Duration", exp(conf.high*10),
                                      ifelse(expo=="Cumulative", exp(conf.high*100), 
                                             ifelse(expo=="Slope", exp(conf.high*0.1),
                                                    ifelse(expo=="Age of onset", exp(conf.high*1), exp(conf.high))))))


toplot <- toplot %>% separate(cancer, into=c(NA, "cancer"), sep="_")


toplot <- toplot %>% mutate(Type=factor(Type, levels=rev(c("Overweight/Obese", "Overall", "Obese"))))

toplot <- toplot %>% mutate(cancer=recode(cancer, 
                                          
                                          "Breast pre"="Breast premenopausal",
                                          "Breast post"="Breast postmenopausal"))

cancer.include <- c("Head and Neck",
                    "Esophagus",
                    "Stomach",
                    "Colorectal",
                    "Liver",
                    "Gallbladder & biliary tract",
                    "Pancreas",
                    "Larynx",
                    "Trachea, bronchus & Lung",
                    "Bone and articular cartilage",
                    "Malignant melanoma of skin",
                    #"Connective and soft tissue",
                    "Breast premenopausal",
                    "Breast postmenopausal",
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
                    "Non-Hodgkin Lymphoma",
                    "Multiple myeloma",
                    "Leukemia"
)
toplot1 <- toplot %>% filter(cancer%in%cancer.include) %>% 
  mutate(cancer=factor(cancer, levels=cancer.include))
toplot1 <- toplot1 %>% mutate(expo=recode(expo, 
                                          
                                          "Duration"="Duration \n(per 10-year increment)",
                                          "Cumulative"="Cumulative exposure \n(per 100-unit increment)",
                                          "Age of onset"="Age of onset \n(per 1-year increment)",
                                          "Slope"="Yearly changes in BMI levels \n(per 0.1-kg/m2 increment)"))

toplot1<-toplot1 %>% 
  mutate(Type=ifelse(group=="Overweight"&exposure=="trajectorySlope", "Overweight/Obese",
                     ifelse(group=="Obese"&exposure=="trajectorySlope", "Obese",
                            as.character(Type)))) %>% 
  mutate(group=ifelse(group=="Overweight"&exposure=="trajectorySlope", "no",
                      ifelse(group=="Obese"&exposure=="trajectorySlope", "no",
                             group)))



ggplot(toplot1 %>% filter(group=="no", Type!="Overall"))+
  geom_errorbar(aes(xmin=low, xmax=up, y=Type),  width = 0.5, colour="grey30")+
  geom_point(aes(HR, Type, col=factor(Type, levels=c("Overweight/Obese", "Overall", "Obese"))),
             size=1.8)+
  
  facet_grid(cancer~factor(expo, levels=c("Duration \n(per 10-year increment)",
                                          "Cumulative exposure \n(per 100-unit increment)",
                                          "Age of onset \n(per 1-year increment)",
                                          "Yearly changes in BMI levels \n(per 0.1-kg/m2 increment)")),
             switch="y" , scales="free")+
  geom_vline(xintercept=1, linetype="dashed", colour="grey20")+
  theme_bw() +
  theme(panel.background = element_blank(),
        # legend.title = element_blank(),
        legend.position = "bottom",
        #panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.text.y.left = element_text(angle = 0),
        axis.title.y.left = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y.left= element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.background = element_rect( fill="grey92"))+ 
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  labs(colour="Legend")






ggplot(toplot1)+
  geom_errorbar(aes(xmin=low, xmax=up, y=group),  width = 0.5, colour="grey30")+
  geom_point(aes(HR, group, col=group),
             size=1.8)+
  
  facet_grid(cancer~factor(expo, levels=c("Duration","Cumulative", "Age of onset", "Slope")),
             switch="y" , scales="free")+
  geom_vline(xintercept=1, linetype="dashed", colour="grey20")+
  theme_bw() +
  theme(panel.background = element_blank(),
        # legend.title = element_blank(),
        legend.position = "bottom",
        #panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.text.y.left = element_text(angle = 0),
        axis.title.y.left = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y.left= element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.background = element_rect( fill="grey92"))+ 
  labs(colour="Legend")

ggsave("cancerLinearNEW.png",
       dpi=300,
       width = 8.5, height = 10)

# table ----
table <- toplot1 %>% select("exposure","cancer"  ,  "expo"   ,   "Type"  ,    "HR"     ,   "low" ,     "up")
table <- table %>% mutate(HR= paste0(format(round(HR,2), nsmall=2), " (",
                                     format(round(low, 2), nsmall=2),"-",format(round(up,2), nsmall=2),
                                     ")")) %>% select(-low, -up, -Type, -expo)
table <- table %>% pivot_wider(names_from=exposure, values_from=HR)
table[match(cancer.include, table$cancer),] %>% htmlTable::htmlTable()
# increment per standard deviation ------------
# add bmi model
outBMI <- tibble()
for (cancer in canceres){
  print(paste0("Working on cancer", cancer))
  test <- pop.analysis %>% filter(.imp!=0) %>% group_by(.imp) %>%
    do(model = coxph(Surv(age.index, age.index + follow_up, get(cancer)) ~
                       bmi.index +sexe + qmedea + nationality + smoking +
                       alcohol+strata(age_gr),  data = .)) %>%
    as.list()
  
  test.pool <- test %>%
    .[[-1]]%>% 
    pool()
  
  tidy <- tidy(test.pool, conf.int = TRUE, conf.level=0.95)
  
  outBMI <- rbind(outBMI, tidy[1,] %>% mutate(exposure="BMI", cancer=cancer))
}

# saca sd de cada exposicion

toplotSD <- out %>% #rbind(outBMI) %>% 
  mutate(expo=recode(exposure,
                     "yearsBMI25more"="Duration overweight",
                     "yearsBMI30more"= "Duration obesity",
                     "cumOver"="Cumulative overweight",
                     "cumObese" ="Cumulative obesity" ,
                     "ageOnsetOver"="Age onset overweight" ,
                     "ageOnsetObese"="Age onset obesity" ,
                     "trajectorySlope"="Yearly BMI change" 
                     
  ))
toplotSD <-toplotSD %>% separate(cancer, into=c(NA, "cancer"), sep="_") %>% 
  mutate(HR=ifelse(expo=="Duration overweight", exp(estimate*sd(pop.analysis$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(estimate*sd(pop.analysis$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(estimate*sd(pop.analysis$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(estimate*sd(pop.analysis$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(estimate*sd(pop.analysis$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(estimate*sd(pop.analysis$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(estimate*sd(pop.analysis$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(estimate*sd(pop.analysis$bmi_atindex)),
                                                                    "Missing")))))))),
         low=ifelse(expo=="Duration overweight", exp(conf.low*sd(pop.analysis$yearsBMI25more)),
                    ifelse(expo=="Duration obesity", exp(conf.low*sd(pop.analysis$yearsBMI30more)),
                           ifelse(expo=="Cumulative overweight", exp(conf.low*sd(pop.analysis$cumOver)), 
                                  ifelse(expo=="Cumulative obesity", exp(conf.low*sd(pop.analysis$cumObese)), 
                                         ifelse(expo=="Age onset overweight", exp(conf.low*sd(pop.analysis$ageOnsetOver, na.rm = T)), 
                                                ifelse(expo=="Age onset obesity", exp(conf.low*sd(pop.analysis$ageOnsetObese, na.rm = T)), 
                                                       ifelse(expo=="Yearly BMI change", exp(conf.low*sd(pop.analysis$trajectorySlope)),
                                                              ifelse(expo=="BMI at index", exp(conf.low*sd(pop.analysis$bmi_atindex)),"Missing")))))))),
         up=ifelse(expo=="Duration overweight", exp(conf.high*sd(pop.analysis$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(conf.high*sd(pop.analysis$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(conf.high*sd(pop.analysis$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(conf.high*sd(pop.analysis$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(conf.high*sd(pop.analysis$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(conf.high*sd(pop.analysis$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(conf.high*sd(pop.analysis$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(conf.high*sd(pop.analysis$bmi_atindex)),
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
  mutate(HR=ifelse(expo=="Duration overweight", exp(estimate*sd((pop.analysis%>%filter(sexe=="D"))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(estimate*sd((pop.analysis%>%filter(sexe=="D"))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(estimate*sd((pop.analysis%>%filter(sexe=="D"))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(estimate*sd((pop.analysis%>%filter(sexe=="D"))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(estimate*sd((pop.analysis%>%filter(sexe=="D"))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(estimate*sd((pop.analysis%>%filter(sexe=="D"))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(estimate*sd((pop.analysis%>%filter(sexe=="D"))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(estimate*sd((pop.analysis%>%filter(sexe=="D"))$bmi_atindex)),
                                                                    "Missing")))))))),
         low=ifelse(expo=="Duration overweight", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D"))$yearsBMI25more)),
                    ifelse(expo=="Duration obesity", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D"))$yearsBMI30more)),
                           ifelse(expo=="Cumulative overweight", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D"))$cumOver)), 
                                  ifelse(expo=="Cumulative obesity", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D"))$cumObese)), 
                                         ifelse(expo=="Age onset overweight", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D"))$ageOnsetOver, na.rm = T)), 
                                                ifelse(expo=="Age onset obesity", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D"))$ageOnsetObese, na.rm = T)), 
                                                       ifelse(expo=="Yearly BMI change", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D"))$trajectorySlope)),
                                                              ifelse(expo=="BMI at index", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D"))$bmi_atindex)),"Missing")))))))),
         up=ifelse(expo=="Duration overweight", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D"))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D"))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D"))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D"))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D"))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D"))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D"))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D"))$bmi_atindex)),
                                                                    "Missing")))))))))

men <- toplotSD %>% filter(cancer%in%c("Prostate", "Testis")) %>% 
  mutate(HR=ifelse(expo=="Duration overweight", exp(estimate*sd((pop.analysis%>%filter(sexe=="H"))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(estimate*sd((pop.analysis%>%filter(sexe=="H"))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(estimate*sd((pop.analysis%>%filter(sexe=="H"))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(estimate*sd((pop.analysis%>%filter(sexe=="H"))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(estimate*sd((pop.analysis%>%filter(sexe=="H"))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(estimate*sd((pop.analysis%>%filter(sexe=="H"))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(estimate*sd((pop.analysis%>%filter(sexe=="H"))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(estimate*sd((pop.analysis%>%filter(sexe=="H"))$bmi_atindex)),
                                                                    "Missing")))))))),
         low=ifelse(expo=="Duration overweight", exp(conf.low*sd((pop.analysis%>%filter(sexe=="H"))$yearsBMI25more)),
                    ifelse(expo=="Duration obesity", exp(conf.low*sd((pop.analysis%>%filter(sexe=="H"))$yearsBMI30more)),
                           ifelse(expo=="Cumulative overweight", exp(conf.low*sd((pop.analysis%>%filter(sexe=="H"))$cumOver)), 
                                  ifelse(expo=="Cumulative obesity", exp(conf.low*sd((pop.analysis%>%filter(sexe=="H"))$cumObese)), 
                                         ifelse(expo=="Age onset overweight", exp(conf.low*sd((pop.analysis%>%filter(sexe=="H"))$ageOnsetOver, na.rm = T)), 
                                                ifelse(expo=="Age onset obesity", exp(conf.low*sd((pop.analysis%>%filter(sexe=="H"))$ageOnsetObese, na.rm = T)), 
                                                       ifelse(expo=="Yearly BMI change", exp(conf.low*sd((pop.analysis%>%filter(sexe=="H"))$trajectorySlope)),
                                                              ifelse(expo=="BMI at index", exp(conf.low*sd((pop.analysis%>%filter(sexe=="H"))$bmi_atindex)),"Missing")))))))),
         up=ifelse(expo=="Duration overweight", exp(conf.high*sd((pop.analysis%>%filter(sexe=="H"))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(conf.high*sd((pop.analysis%>%filter(sexe=="H"))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(conf.high*sd((pop.analysis%>%filter(sexe=="H"))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(conf.high*sd((pop.analysis%>%filter(sexe=="H"))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(conf.high*sd((pop.analysis%>%filter(sexe=="H"))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(conf.high*sd((pop.analysis%>%filter(sexe=="H"))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(conf.high*sd((pop.analysis%>%filter(sexe=="H"))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(conf.high*sd((pop.analysis%>%filter(sexe=="H"))$bmi_atindex)),
                                                                    "Missing")))))))))


#toplot <- toplot %>% mutate(Type=factor(Type, levels=rev(c("Overweight/Obese", "Overall", "Obese"))))
toplotSD <- toplotSD %>% filter(!cancer%in%c( "Cervix Uteri","Corpus Uteri","Ovary","Prostate", "Testis")) %>% 
  rbind(women) %>% 
  rbind(men)

toplotSD <- toplotSD %>% filter(cancer%in%cancer.include) %>% 
  mutate(cancer=factor(cancer, levels=cancer.include))

toplotSD <- toplotSD %>%mutate(HR1=as.numeric(HR),
                               low=as.numeric(low),
                               up=as.numeric(up))

toplotSD <- toplotSD %>%
  mutate(expo=factor(expo, levels=rev(c("BMI","Duration overweight", "Duration obesity",
                                        "Cumulative overweight", "Cumulative obesity",
                                        "Age onset overweight" , "Age onset obesity", 
                                        "Yearly BMI change" ))))

ggplot(toplotSD%>% filter(group=="no", Type!="Overall"))+
  geom_errorbar(aes(xmin=low, xmax=up, y=expo),  width = 0.5, colour="grey30")+
  geom_point(aes(HR1, expo, col=factor(Type, levels=c("Overweight/Obese", "Overall", "Obese"))),
             size=1.8)+
  facet_wrap(cancer~.)+
  geom_vline(xintercept=1, linetype="dashed", colour="grey20")+
  theme_bw() +
  theme(panel.background = element_blank(),
        # legend.title = element_blank(),
        legend.position = "bottom",
        #panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        # strip.text.y.left = element_text(angle = 0),
        axis.title.y.left = element_blank(),
        #axis.ticks.y= element_blank(),
        #axis.text.y.left= element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.background = element_rect( fill="grey92"))+ 
  coord_cartesian(xlim=c(0.63, 1.6))+
  labs(colour="Legend")+
  xlab("HR")

ggsave("exposure1SDNEW.png",
       dpi=300,
       width = 9, height = 9)

# table ----
table <- toplotSD %>% select("exposure","cancer"  ,  "expo"   ,   "Type"  ,    "HR1"     ,   "low" ,     "up")
table <- table %>% mutate(HR1= paste0(format(round(HR1,2), nsmall=2), " (",
                                      format(round(low, 2), nsmall=2),"-",format(round(up,2), nsmall=2),
                                      ")")) %>% select(-low, -up, -Type, -expo)
table <- table %>% pivot_wider(names_from=exposure, values_from=HR1)
table[match(cancer.include, table$cancer),] %>% htmlTable::htmlTable()



test <- test %>% mutate(cancer=recode(cancer, 
                                      "Gallbladder & biliary tract"= "Gallbladder & \nbiliary tract",
                                      "Trachea, bronchus & Lung" ="Trachea, \nbronchus & Lung" ,
                                      "Bone and articular cartilage"="Bone and \narticular cartilage",
                                      "Malignant melanoma of skin" ="Malignant melanoma \nof skin" 
))
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
                  "Breast premenopausal",
                  "Breast postmenopausal",
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
                  "Non-Hodgkin Lymphoma",
                  "Multiple myeloma",
                  "Leukemia"
)

ggplot()+
  geom_hline(yintercept = 1, colour = "#000000",
             linetype=2, size=0.8)+
  geom_line(data=test %>% filter(expo=="Cumulative"),aes(value, HR, colour=expo), size=2) +
  geom_ribbon(data=test %>% filter(expo=="Cumulative"),aes(x=value, y=HR,ymin=low, ymax=up, fill=expo), 
              alpha=0.4, linetype=0,show.legend=F)+
  
  geom_line(data=test %>% filter(!expo=="Cumulative"),aes(value*coef, HR, colour=expo), size=2) +
  geom_ribbon(data=test %>% filter(!expo=="Cumulative"),aes(x=value*coef, y=HR,ymin=low, ymax=up, fill=expo), 
              alpha=0.4, linetype=0,show.legend=F)+
  facet_wrap(factor(cancer, levels=order.cancer)~.)+#, scales="free_y")+
  theme_bw() +
  theme(panel.spacing = unit(0, "lines"),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text =element_text(size=12),
        panel.grid.major.y= element_line(colour = "grey", linetype="dashed"),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        #axis.text.x=element_text(size=12),
        #axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size=9, face="bold"),
        strip.background = element_rect( fill="#f7f7f7"),
        strip.text.y.left = element_text(angle = 0)) +
  ylab("Relative\nhazard ratio\n")+
  #xlab("\nYears lived with a BMI greater than 25")+
  scale_x_continuous( name = "Cumulative overweight exposure",
                      sec.axis = sec_axis(~./coef, name="Duration of overweight in years"),
                      breaks = c(0,300,600)) +
  coord_cartesian(ylim = c(0.5, 2))

ggsave("cancerSplineContr.png",
       dpi=300,
       width = 8.7, height = 10)

# las otras 2 exposiciones
test <- test %>% filter(exposure%in%c("ageOnsetOver", "trajectorySlope"))
summary((test %>% filter(exposure=="ageOnsetOver"))$value)
summary((test %>% filter(!exposure=="ageOnsetOver"))$value)
coef<-8.9
term<-27.5

ggplot()+
  geom_hline(yintercept = 1, colour = "#000000",
             linetype=2, size=0.8)+
  geom_line(data=test %>% filter(exposure=="ageOnsetOver"),aes(value, HR, colour=expo), size=2) +
  geom_ribbon(data=test %>% filter(exposure=="ageOnsetOver"),aes(x=value, y=HR,ymin=low, ymax=up, fill=expo), 
              alpha=0.4, linetype=0,show.legend=F)+
  
  geom_line(data=test %>% filter(!exposure=="ageOnsetOver"),aes(value*coef+term, HR, colour=expo), size=2) +
  geom_ribbon(data=test %>% filter(!exposure=="ageOnsetOver"),aes(x=value*coef+term, y=HR,ymin=low, ymax=up, fill=expo), 
              alpha=0.4, linetype=0,show.legend=F)+
  facet_wrap(factor(cancer, levels=order.cancer)~.)+#, scales="free_y")+
  theme_bw() +
  theme(panel.spacing = unit(0, "lines"),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text =element_text(size=12),
        panel.grid.major.y= element_line(colour = "grey", linetype="dashed"),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        #axis.text.x=element_text(size=12),
        #axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size=9, face="bold"),
        strip.background = element_rect( fill="#f7f7f7"),
        strip.text.y.left = element_text(angle = 0)) +
  ylab("Relative\nhazard ratio\n")+
  #xlab("\nYears lived with a BMI greater than 25")+
  scale_x_continuous( name = "Age of onset of overweight",
                      sec.axis = sec_axis(~./coef-term, name="Yearly change in BMI",
                                          labels=round(seq(from=min((test %>% filter(!exposure=="ageOnsetOver"))$value),
                                                           to =max((test %>% filter(!exposure=="ageOnsetOver"))$value),
                                                           length.out = 6),1))
                      #    , breaks = c(0,300,600)
  ) +
  coord_cartesian(ylim = c(0.5, 2))

ggsave("cancerSplineContr2.png",
       dpi=300,
       width = 8.7, height = 10)

# las otras 2 exposiciones
test <- test %>% filter(exposure%in%c("ageOnsetObese", "trajectorySlope"))
summary((test %>% filter(exposure=="ageOnsetObese"))$value)
summary((test %>% filter(!exposure=="ageOnsetObese"))$value)
coef<-8.9
term<-27.5

ggplot()+
  geom_hline(yintercept = 1, colour = "#000000",
             linetype=2, size=0.8)+
  geom_line(data=test %>% filter(exposure=="ageOnsetObese"),aes(value, HR, colour=expo), size=2) +
  geom_ribbon(data=test %>% filter(exposure=="ageOnsetObese"),aes(x=value, y=HR,ymin=low, ymax=up, fill=expo), 
              alpha=0.4, linetype=0,show.legend=F)+
  
  geom_line(data=test %>% filter(!exposure=="ageOnsetObese"),aes(value*coef+term, HR, colour=expo), size=2) +
  geom_ribbon(data=test %>% filter(!exposure=="ageOnsetObese"),aes(x=value*coef+term, y=HR,ymin=low, ymax=up, fill=expo), 
              alpha=0.4, linetype=0,show.legend=F)+
  facet_wrap(factor(cancer, levels=order.cancer)~.)+#, scales="free_y")+
  theme_bw() +
  theme(panel.spacing = unit(0, "lines"),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text =element_text(size=12),
        panel.grid.major.y= element_line(colour = "grey", linetype="dashed"),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        #axis.text.x=element_text(size=12),
        #axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size=9, face="bold"),
        strip.background = element_rect( fill="#f7f7f7"),
        strip.text.y.left = element_text(angle = 0)) +
  ylab("Relative\nhazard ratio\n")+
  #xlab("\nYears lived with a BMI greater than 25")+
  scale_x_continuous( name = "Age of onset of overweight",
                      sec.axis = sec_axis(~./coef-term, name="Yearly change in BMI",
                                          labels=round(seq(from=min((test %>% filter(!exposure=="ageOnsetObese"))$value),
                                                           to =max((test %>% filter(!exposure=="ageOnsetObese"))$value),
                                                           length.out = 6),2))
                      #    , breaks = c(0,300,600)
  ) +
  coord_cartesian(ylim = c(0.5, 2))


# linear , slope stratified for increasing.. -------
pop.analysis <- pop.analysis %>% 
  mutate(splope.gr=ifelse(trajectorySlope<=(-0.1), "decreasing",
                          ifelse(trajectorySlope>=0.1, "increasing", "stable")
  ))

pop.analysis <- pop.analysis %>% 
  mutate(splope.gr=factor(splope.gr, levels=c("stable", "decreasing", "increasing")))

out <- list()
i = 0
for (cancer in canceres){
  print(paste0("Working on cancer", cancer))
  
  if(cancer =="cancer_Breast post"){
    current.pop<-pop.analysis %>% filter(meno_status_index==1) %>% group_by(.imp)
    print("breast post")
  } else {
    if(cancer=="cancer_Breast pre"){
      current.pop<-pop.analysis %>% filter(meno_status_index==0)%>% group_by(.imp)
      print("breast pre")
    }
    else{current.pop<-pop.analysis %>% group_by(.imp)
    print("normal")}
  }
  
  
  for(bmigr in c("Overweight", "Obese")){
    print(paste0("working on ", bmigr))
    
    if(cancer %in%c("cancer_Breast post","cancer_Breast pre" )){
      test <- current.pop %>% 
        do(model = coxph(Surv(age.index, age.index + follow_up, get(cancer)) ~
                           splope.gr  + qmedea + nationality + smoking +
                           alcohol+strata(age_gr), subset=(bmi.gr18==bmigr), data = .))  
    }else{
      test <- current.pop %>% 
        do(model = coxph(Surv(age.index, age.index + follow_up, get(cancer)) ~
                           splope.gr +sexe + qmedea + nationality + smoking +
                           alcohol+strata(age_gr), subset=(bmi.gr18==bmigr), data = .))   
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
    out[[i]]<-tidy[1:2,] %>% 
      mutate(exposure=exposure, cancer=cancer, model="Slope stratified",group=bmigr, 
             aic1=aic[1], aic2=aic[2], aic3=aic[3], aic4=aic[4],aic5=aic[5],
             bic1=bic[1], bic2=bic[2], bic3=bic[3], bic4=bic[4],bic5=bic[5])
  }
  
}


df <- plyr::ldply(out, tibble)
saveRDS(df, "splopeCatogories.rds")

out <- df
toplot <- out %>% 
  mutate(expo=ifelse(exposure%in%c("yearsBMI25more", "yearsBMI30more"), "Duration",
                     ifelse(exposure%in%c( "cumOver","cumObese" ), "Cumulative",
                            ifelse(exposure%in%c( "ageOnsetOver","ageOnsetObese" ), "Age of onset",
                                   "Slope"))),
         Type=ifelse(exposure%in%c("yearsBMI25more", "cumOver","ageOnsetOver"), "Overweight/Obese",
                     ifelse(exposure%in%c("yearsBMI30more", "cumObese","ageOnsetObese"), "Obese",
                            "Overall")))

toplot <- toplot %>% mutate(HR=ifelse(expo=="Duration", exp(estimate*10),
                                      ifelse(expo=="Cumulative", exp(estimate*100), 
                                             ifelse(expo=="Slope", exp(estimate*0.1),
                                                    ifelse(expo=="Age of onset", exp(estimate*1),exp(estimate))))),
                            low=ifelse(expo=="Duration", exp(conf.low*10),
                                       ifelse(expo=="Cumulative", exp(conf.low*100), 
                                              ifelse(expo=="Slope", exp(conf.low*0.1),
                                                     ifelse(expo=="Age of onset", exp(conf.low*1), exp(conf.low))))),
                            up=ifelse(expo=="Duration", exp(conf.high*10),
                                      ifelse(expo=="Cumulative", exp(conf.high*100), 
                                             ifelse(expo=="Slope", exp(conf.high*0.1),
                                                    ifelse(expo=="Age of onset", exp(conf.high*1), exp(conf.high))))))


toplot <- toplot %>% separate(cancer, into=c(NA, "cancer"), sep="_")


toplot <- toplot %>% mutate(Type=factor(Type, levels=rev(c("Overweight/Obese", "Overall", "Obese"))))

toplot <- toplot %>% mutate(cancer=recode(cancer, 
                                          
                                          "Breast pre"="Breast premenopausal",
                                          "Breast post"="Breast postmenopausal"))


toplot1 <- toplot %>% filter(cancer%in%cancer.include) %>% 
  mutate(cancer=factor(cancer, levels=cancer.include))
toplot1 <- toplot1 %>% mutate(expo=recode(expo, 
                                          
                                          "Duration"="Duration \n(per 10-year increment)",
                                          "Cumulative"="Cumulative exposure \n(per 100-unit increment)",
                                          "Age of onset"="Age of onset \n(per 1-year increment)",
                                          "Slope"="Yearly changes in BMI levels \n(per 0.1-kg/m2 increment)"))



ggplot(toplot1 %>% filter(term=="splope.grdecreasing"))+
  geom_errorbar(aes(xmin=low, xmax=up, y=group),  width = 0.5, colour="grey30")+
  geom_point(aes(HR, group, col=group),
             size=1.8)+
  
  facet_grid(cancer~factor(expo, levels=c("Duration","Cumulative", "Age of onset", "Yearly changes in BMI levels \n(per 0.1-kg/m2 increment)")),
             switch="y" , scales="free")+
  geom_vline(xintercept=1, linetype="dashed", colour="grey20")+
  theme_bw() +
  theme(panel.background = element_blank(),
        # legend.title = element_blank(),
        legend.position = "bottom",
        #panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.text.y.left = element_text(angle = 0),
        axis.title.y.left = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y.left= element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.background = element_rect( fill="grey92"))+ 
  xlim(c(0.2, 1.5))+
  labs(colour="Legend")


saveRDS(df, "out.coxSlopeGR.rds")
# Harrell index -------
install.packages("dynpred")
library(dynpred)
?cindex 

m1 <- coxph(Surv(age.index, age.index + follow_up, get(cancer)) ~
              get(exposure) +sexe + qmedea + nationality + smoking +
              alcohol+strata(age_gr),  data = pop.analysis %>% group_by(.imp))
# maybe likelihood test?
m0 <- coxph(Surv(age.index, age.index + follow_up, get(cancer)) ~
              sexe + qmedea + nationality + smoking +
              alcohol+strata(age_gr),  data = pop.analysis %>% group_by(.imp))
anova(m0, m1) # if <0.05, H0 rejected so we need to include the variable

m1 <- coxph(Surv(age.index, age.index + follow_up, get(cancer)) ~
              bmi.index+sexe + qmedea + nationality + smoking +
              alcohol+strata(age_gr),  data = pop.analysis %>% filter(.imp==1))


sum.surv <- summary(m1)

c_index<-tibble()
for (i in 1:5){
  c_index <- rbind(c_index,(test$model[[i]])$concordance)
}

# aic ------

aic <- readRDS("out.coxAllAdjustNEW2.rds")

aic<-aic %>% filter(model%in%c("Null model","BMI at index","Fully adjusted" ))

aic<-aic %>% mutate(BIC=(bic1+ bic2+bic3+bic4+bic5)/5, 
                    AIC=(aic1+ aic2+aic3+aic4+aic5)/5)

table <- aic %>% select(cancer, model, BIC, exposure)
table <- table %>% mutate(exposure=ifelse(model!="Fully adjusted", model, exposure)) 
table<-table %>% select(-model) %>% distinct()
table <- table%>% 
  pivot_wider(names_from=exposure, values_from=BIC)
table <- table %>% separate(cancer, into=c(NA, "cancer"), sep="_")
table[match(cancer.include, table$cancer),] %>% htmlTable::htmlTable()

table <- aic %>% select(cancer, model, AIC, exposure)
table <- table %>% mutate(exposure=ifelse(model!="Fully adjusted", model, exposure)) 
table<-table %>% select(-model) %>% distinct()
table <- table%>% 
  pivot_wider(names_from=exposure, values_from=AIC)
table <- table %>% separate(cancer, into=c(NA, "cancer"), sep="_")
table[match(cancer.include, table$cancer),] %>% htmlTable::htmlTable()

# harrell
cancer.include <- c( "Corpus Uteri",
                     "Kidney",
                     "Gallbladder & biliary tract",
                     "Thyroid",
                     "Breast postmenopausal",
                     "Leukemia",
                     "Multiple myeloma",
                     "Testis",
                     "Brain and CNS",
                     "Colorectal",
                     "Hodgkin lymphoma",
                     "Liver",
                     "Ovary",
                     "Non-Hodgkin Lymphoma",
                     "Malignant melanoma of skin",
                     "Cervix Uteri",
                     "Prostate",
                     "Bladder",
                     "Bone and articular cartilage",
                     "Pancreas",
                     "Breast premenopausal",
                     "Stomach",
                     "Head and Neck",
                     "Trachea, bronchus & Lung",
                     "Esophagus",
                     "Larynx"
                     
)
char <- readRDS("E:/Andrea/Trajectories OBECAN/charlson.rds")
h<-char %>% filter(model%in%c("Null model","BMI at index","Fully adjusted" ))

table <- h %>% select(cancer, model, concordance, exposure, std)
table <- table %>% mutate(exposure=ifelse(model!="Fully adjusted", model, exposure)) 
table<-table %>% select(-model) %>% distinct()
table <- table%>% 
  pivot_wider(names_from=exposure, values_from=c(concordance, std))
table <- table %>% separate(cancer, into=c(NA, "cancer"), sep="_")%>% 
  mutate(cancer=recode(cancer, 
                       "Breast pre"="Breast premenopausal",
                       "Breast post"="Breast postmenopausal"
  ))
table[match(cancer.include, table$cancer),] %>% htmlTable::htmlTable()

library(kableExtra)

fit<-table
#res<-matrix(nrow=30, ncol=10)
fit <- fit %>%
  mutate_at(2:ncol(fit), funs(as.numeric(str_replace_all(., ",", "")))) %>% 
  mutate_at(2:ncol(fit), funs(round(.,5)))
# for (i in 2:nrow(fit)) {
#   val.min <- min(fit[i,2:ncol(fit)])
#   j<-which(fit[i, 2:ncol(fit)]==val.min)
#   res[i,j]<- cell_spec(fit[i,j], background = "orange")
# }
# 
# kable(escape = F, res) %>% 
#   kable_styling(bootstrap_options = c("striped", "bordered"))



df <- pivot_longer(fit, cols=2:length(fit),names_to=c("tipo","variable"), values_to=c("value"), names_sep="_") 
df<-df %>% pivot_wider( names_from = "tipo", values_from="value")
df <- df %>% rename(value=concordance, sd=std)
df <- df %>% group_by(cancer) %>% arrange(value) %>% mutate(val_norm=row_number())%>%
  ungroup()  %>% 
  filter(cancer%in%cancer.include) %>% 
  mutate(cancer=factor(cancer, levels=rev(cancer.include))) 

df<-df %>% mutate(char=paste0(format(round(value,3), nsmall=3),
                              "\n(", 
                              format(round(value-1.96*sd, 3), nsmall=3),
                              "-", 
                              format(round(value+1.96*sd, 3), nsmall=3),
                              ")"))

ggplot(data = df ,#%>% filter(!variable%in%c( "ageOnsetObese"))  , 
       aes(x = factor(variable, levels=c("Null model", "BMI at index",
                                         "yearsBMI25more",  "yearsBMI30more",
                                         "cumOver" , "cumObese"  ,
                                         "ageOnsetOver" ,   "ageOnsetObese", 
                                         "trajectorySlope")), y = cancer,
           group=cancer),
       position_dodge()) +
  geom_tile(aes(fill = val_norm), colour = "white")+
  scale_fill_gradient2(low="red", mid="orange",high="green", midpoint=4, na.value="black") +     
  geom_text(aes(label=char), lineheight = .7)+
  theme_minimal()+
  theme(legend.position="none", 
        text=element_text(size = 16),
        axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 8))+
  xlab("")+
  ylab("")+
  scale_x_discrete(labels = (c("Null model","BMI at index date",
                               expression("Duration BMI">="25"),
                               expression("Duration BMI">="30"),
                               expression("Cum. exposure BMI">="25"),
                               expression("Cum. exposure BMI">="30"),
                               expression("Age onset BMI">="25"),
                               expression("Age onset BMI">="30"))))

ggsave("harrellindex.png",
       dpi=300,
       width = 14, height = 11)

# final - increase in SD -----
Ncasos <- pop.analysis %>% select(starts_with("cancer")) %>% summarise(across(where(is.numeric), sum)) %>% 
  mutate_all(funs(./5))

out<-readRDS("out.95.rds")
#out1<-readRDS("out.95_duration.rds")

toplotSD <- out %>% 
  #filter(!exposure%in%c("ageOnsetOver", "ageOnsetObese")) %>% 
  # rbind(out1) %>% 
  mutate(expo=recode(exposure,
                     "yearsBMI25more"="Duration overweight",
                     "yearsBMI30more"= "Duration obesity",
                     "cumOver"="Cumulative overweight",
                     "cumObese" ="Cumulative obesity" ,
                     "ageOnsetOver"="Age onset overweight" ,
                     "ageOnsetObese"="Age onset obesity" 
                     
  )) %>% mutate(expo=ifelse(model!="Fully adjusted", model, expo))

toplotSD <-toplotSD %>% mutate(
  estimate=ifelse(expo%in%c("Age onset overweight" , "Age onset obesity"), -estimate, estimate),
  conf.low=ifelse(expo%in%c("Age onset overweight" , "Age onset obesity"), -conf.low, conf.low),
  conf.high=ifelse(expo%in%c("Age onset overweight" , "Age onset obesity"), -conf.high, conf.high))

Ncasos1 <- Ncasos %>% pivot_longer(cols=3:length(Ncasos)) %>% select(-(1:2))

toplotSD <-toplotSD %>%left_join(Ncasos1 %>% rename(cancer=name, N=value))

toplotSD <-toplotSD %>% separate(cancer, into=c(NA, "cancer"), sep="_") %>% 
  mutate(HR=ifelse(expo=="Duration overweight", exp(estimate*sd(pop.analysis$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(estimate*sd(pop.analysis$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(estimate*sd(pop.analysis$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(estimate*sd(pop.analysis$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(estimate*sd(pop.analysis$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(estimate*sd(pop.analysis$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(estimate*sd(pop.analysis$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(estimate*sd(pop.analysis$bmi_atindex)),
                                                                    "Missing")))))))),
         low=ifelse(expo=="Duration overweight", exp(conf.low*sd(pop.analysis$yearsBMI25more)),
                    ifelse(expo=="Duration obesity", exp(conf.low*sd(pop.analysis$yearsBMI30more)),
                           ifelse(expo=="Cumulative overweight", exp(conf.low*sd(pop.analysis$cumOver)), 
                                  ifelse(expo=="Cumulative obesity", exp(conf.low*sd(pop.analysis$cumObese)), 
                                         ifelse(expo=="Age onset overweight", exp(conf.low*sd(pop.analysis$ageOnsetOver, na.rm = T)), 
                                                ifelse(expo=="Age onset obesity", exp(conf.low*sd(pop.analysis$ageOnsetObese, na.rm = T)), 
                                                       ifelse(expo=="Yearly BMI change", exp(conf.low*sd(pop.analysis$trajectorySlope)),
                                                              ifelse(expo=="BMI at index", exp(conf.low*sd(pop.analysis$bmi_atindex)),"Missing")))))))),
         up=ifelse(expo=="Duration overweight", exp(conf.high*sd(pop.analysis$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(conf.high*sd(pop.analysis$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(conf.high*sd(pop.analysis$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(conf.high*sd(pop.analysis$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(conf.high*sd(pop.analysis$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(conf.high*sd(pop.analysis$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(conf.high*sd(pop.analysis$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(conf.high*sd(pop.analysis$bmi_atindex)),
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
  mutate(HR=ifelse(expo=="Duration overweight", exp(estimate*sd((pop.analysis%>%filter(sexe=="D"))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(estimate*sd((pop.analysis%>%filter(sexe=="D"))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(estimate*sd((pop.analysis%>%filter(sexe=="D"))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(estimate*sd((pop.analysis%>%filter(sexe=="D"))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(estimate*sd((pop.analysis%>%filter(sexe=="D"))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(estimate*sd((pop.analysis%>%filter(sexe=="D"))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(estimate*sd((pop.analysis%>%filter(sexe=="D"))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(estimate*sd((pop.analysis%>%filter(sexe=="D"))$bmi_atindex)),
                                                                    "Missing")))))))),
         low=ifelse(expo=="Duration overweight", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D"))$yearsBMI25more)),
                    ifelse(expo=="Duration obesity", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D"))$yearsBMI30more)),
                           ifelse(expo=="Cumulative overweight", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D"))$cumOver)), 
                                  ifelse(expo=="Cumulative obesity", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D"))$cumObese)), 
                                         ifelse(expo=="Age onset overweight", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D"))$ageOnsetOver, na.rm = T)), 
                                                ifelse(expo=="Age onset obesity", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D"))$ageOnsetObese, na.rm = T)), 
                                                       ifelse(expo=="Yearly BMI change", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D"))$trajectorySlope)),
                                                              ifelse(expo=="BMI at index", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D"))$bmi_atindex)),"Missing")))))))),
         up=ifelse(expo=="Duration overweight", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D"))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D"))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D"))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D"))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D"))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D"))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D"))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D"))$bmi_atindex)),
                                                                    "Missing")))))))))

men <- toplotSD %>% filter(cancer%in%c("Prostate", "Testis")) %>% 
  mutate(HR=ifelse(expo=="Duration overweight", exp(estimate*sd((pop.analysis%>%filter(sexe=="H"))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(estimate*sd((pop.analysis%>%filter(sexe=="H"))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(estimate*sd((pop.analysis%>%filter(sexe=="H"))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(estimate*sd((pop.analysis%>%filter(sexe=="H"))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(estimate*sd((pop.analysis%>%filter(sexe=="H"))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(estimate*sd((pop.analysis%>%filter(sexe=="H"))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(estimate*sd((pop.analysis%>%filter(sexe=="H"))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(estimate*sd((pop.analysis%>%filter(sexe=="H"))$bmi_atindex)),
                                                                    "Missing")))))))),
         low=ifelse(expo=="Duration overweight", exp(conf.low*sd((pop.analysis%>%filter(sexe=="H"))$yearsBMI25more)),
                    ifelse(expo=="Duration obesity", exp(conf.low*sd((pop.analysis%>%filter(sexe=="H"))$yearsBMI30more)),
                           ifelse(expo=="Cumulative overweight", exp(conf.low*sd((pop.analysis%>%filter(sexe=="H"))$cumOver)), 
                                  ifelse(expo=="Cumulative obesity", exp(conf.low*sd((pop.analysis%>%filter(sexe=="H"))$cumObese)), 
                                         ifelse(expo=="Age onset overweight", exp(conf.low*sd((pop.analysis%>%filter(sexe=="H"))$ageOnsetOver, na.rm = T)), 
                                                ifelse(expo=="Age onset obesity", exp(conf.low*sd((pop.analysis%>%filter(sexe=="H"))$ageOnsetObese, na.rm = T)), 
                                                       ifelse(expo=="Yearly BMI change", exp(conf.low*sd((pop.analysis%>%filter(sexe=="H"))$trajectorySlope)),
                                                              ifelse(expo=="BMI at index", exp(conf.low*sd((pop.analysis%>%filter(sexe=="H"))$bmi_atindex)),"Missing")))))))),
         up=ifelse(expo=="Duration overweight", exp(conf.high*sd((pop.analysis%>%filter(sexe=="H"))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(conf.high*sd((pop.analysis%>%filter(sexe=="H"))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(conf.high*sd((pop.analysis%>%filter(sexe=="H"))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(conf.high*sd((pop.analysis%>%filter(sexe=="H"))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(conf.high*sd((pop.analysis%>%filter(sexe=="H"))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(conf.high*sd((pop.analysis%>%filter(sexe=="H"))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(conf.high*sd((pop.analysis%>%filter(sexe=="H"))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(conf.high*sd((pop.analysis%>%filter(sexe=="H"))$bmi_atindex)),
                                                                    "Missing")))))))))

breast_pre <-toplotSD %>% filter(cancer%in%c( "Breast pre")) %>% 
  mutate(HR=ifelse(expo=="Duration overweight", exp(estimate*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(estimate*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(estimate*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(estimate*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(estimate*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(estimate*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(estimate*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(estimate*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$bmi_atindex)),
                                                                    "Missing")))))))),
         low=ifelse(expo=="Duration overweight", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$yearsBMI25more)),
                    ifelse(expo=="Duration obesity", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$yearsBMI30more)),
                           ifelse(expo=="Cumulative overweight", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$cumOver)), 
                                  ifelse(expo=="Cumulative obesity", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$cumObese)), 
                                         ifelse(expo=="Age onset overweight", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$ageOnsetOver, na.rm = T)), 
                                                ifelse(expo=="Age onset obesity", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$ageOnsetObese, na.rm = T)), 
                                                       ifelse(expo=="Yearly BMI change", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$trajectorySlope)),
                                                              ifelse(expo=="BMI at index", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$bmi_atindex)),"Missing")))))))),
         up=ifelse(expo=="Duration overweight", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$bmi_atindex)),
                                                                    "Missing")))))))))

breast_post <-toplotSD %>% filter(cancer%in%c( "Breast post")) %>% 
  mutate(HR=ifelse(expo=="Duration overweight", exp(estimate*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(estimate*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(estimate*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(estimate*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(estimate*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(estimate*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(estimate*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(estimate*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$bmi_atindex)),
                                                                    "Missing")))))))),
         low=ifelse(expo=="Duration overweight", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$yearsBMI25more)),
                    ifelse(expo=="Duration obesity", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$yearsBMI30more)),
                           ifelse(expo=="Cumulative overweight", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$cumOver)), 
                                  ifelse(expo=="Cumulative obesity", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$cumObese)), 
                                         ifelse(expo=="Age onset overweight", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$ageOnsetOver, na.rm = T)), 
                                                ifelse(expo=="Age onset obesity", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$ageOnsetObese, na.rm = T)), 
                                                       ifelse(expo=="Yearly BMI change", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$trajectorySlope)),
                                                              ifelse(expo=="BMI at index", exp(conf.low*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$bmi_atindex)),"Missing")))))))),
         up=ifelse(expo=="Duration overweight", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$yearsBMI25more)),
                   ifelse(expo=="Duration obesity", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$yearsBMI30more)),
                          ifelse(expo=="Cumulative overweight", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$cumOver)), 
                                 ifelse(expo=="Cumulative obesity", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$cumObese)), 
                                        ifelse(expo=="Age onset overweight", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$ageOnsetOver, na.rm = T)), 
                                               ifelse(expo=="Age onset obesity", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$ageOnsetObese, na.rm = T)), 
                                                      ifelse(expo=="Yearly BMI change", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$trajectorySlope)),
                                                             ifelse(expo=="BMI at index", exp(conf.high*sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$bmi_atindex)),
                                                                    "Missing")))))))))

#toplot <- toplot %>% mutate(Type=factor(Type, levels=rev(c("Overweight/Obese", "Overall", "Obese"))))
toplotSD <- toplotSD %>% filter(!cancer%in%c( "Cervix Uteri","Corpus Uteri","Ovary","Prostate", "Testis", "Breast pre", "Breast post")) %>% 
  rbind(women) %>% 
  rbind(men) %>% 
  rbind(breast_pre) %>% 
  rbind(breast_post)

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

# añadir valor sd en nombres
# paste0("Duration overweight", round(sd(pop.analysis$yearsBMI25more),2)),
# ifelse(expo=="Duration obesity", exp(conf.high*sd(pop.analysis$yearsBMI30more)),
#        ifelse(expo=="Cumulative overweight", exp(conf.high*sd(pop.analysis$cumOver)), 
#               ifelse(expo=="Cumulative obesity", exp(conf.high*sd(pop.analysis$cumObese)), 
#                      ifelse(expo=="Age onset overweight", exp(conf.high*sd(pop.analysis$ageOnsetOver, na.rm = T)), 
#                             ifelse(expo=="Age onset obesity", exp(conf.high*sd(pop.analysis$ageOnsetObese, na.rm = T)), 
#                                    ifelse(expo=="Yearly BMI change", exp(conf.high*sd(pop.analysis$trajectorySlope)),
#                                           ifelse(expo=="BMI at index", exp(conf.high*sd(pop.analysis$bmi_atindex)
# 
#                                           )))))))))


order.cancer1<-(toplotSD%>%filter(group=="no", Type!="Overall", expo=="BMI at index")%>%arrange(HR1) %>% distinct(cancer))$cancer
#View((toplotSD%>%filter(group=="no", Type!="Overall", expo=="BMI at index")%>%arrange(HR1) ))
#saveRDS(order.cancer1, "order.cancer.rds")
toplotSD<-toplotSD %>% 
  mutate(cancer=factor(cancer, levels=rev(order.cancer1)))



toplotSD <- toplotSD %>% 
  mutate(label=paste0(cancer, "\n(N =", format(N, big.mark=","), ")"))

ggplot(toplotSD%>% filter(group=="no", Type!="Overall", expo!="Missing"))+
  geom_errorbar(aes(xmin=low, xmax=up, y=expo),  width = 0.5, colour="grey30")+
  geom_point(aes(HR1, expo, col=expo, shape=expo3),
             size=1.5)+
  facet_wrap(cancer~., ncol=5, labeller=as_labeller(setNames(toplotSD$label,toplotSD$cancer)))+
  geom_vline(xintercept=1, linetype="dashed", colour="grey20")+
  theme_bw() +
  theme(panel.background = element_blank(),
        axis.text=element_text(size=7),
        # legend.title = element_blank(),
        legend.position = "null",
        #panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        # strip.text.y.left = element_text(angle = 0),
        axis.title.y.left = element_blank(),
        #axis.ticks.y= element_blank(),
        #axis.text.y.left= element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.background = element_rect( fill="grey92"),
        strip.text = element_text(size=7))+ 
  
  labs(colour="Legend")+
  scale_y_discrete(labels = rev(c("BMI at index date",
                                  expression("Duration BMI">="25"),
                                  expression("Duration BMI">="30"),
                                  expression("Cum. exposure BMI">="25"),
                                  expression("Cum. exposure BMI">="30"),
                                  expression("Age onset BMI">="25"),
                                  expression("Age onset BMI">="30"))))+
  scale_color_manual(values=rev(c("black","#03992A", "#FF9F2F", "#32A5DD","#FF5B2F","#C51A95","#AA6C39")))+
  scale_shape_manual(values=c(15,16,17,8))+
  xlab("HR (95% CI)")+
  coord_cartesian(xlim=c(0.7, 1.6))


# ggsave("exposureSD_lastN2_correct_rev.png",
#        dpi=300,
#        width = 9, height = 9)

ggsave("editor/figura1.pdf",
       dpi=300,
       width = 180, height = 185,
       units="mm")

ggsave("editor/figura1.png",
       dpi=300,
       width = 180, height = 185,
       units="mm")


#table hR----
table<-toplotSD%>% filter(group=="no", Type!="Overall", expo!="Missing")

htmlTable::htmlTable(table %>% select(HR1, cancer, expo, low, up, N))

table <- table %>% select(HR1, cancer, expo, low, up)
table <- table %>% 
  mutate(HR=paste0(format(round(HR1,2),nsmall=2),
                   " (", 
                   format(round(low, 2),nsmall=2), "-", 
                   format(round(up,2),nsmall=2), ")")) %>% 
  select(-HR1, -low,-up)



table <- table %>% pivot_wider(names_from = "expo", values_from="HR")
table[match(rev(order.cancer1), table$cancer),] %>% htmlTable::htmlTable()

# table sd ------
sd<-list()
sd[["Duration overweight"]]<-sd(pop.analysis$yearsBMI25more)
sd[["Duration obesity"]]<-sd(pop.analysis$yearsBMI30more)
sd[["Cumulative overweight"]]<-sd(pop.analysis$cumOver)
sd[["Cumulative obesity"]]<-sd(pop.analysis$cumObese) 
sd[["Age onset overweight"]]<-sd(pop.analysis$ageOnsetOver, na.rm = T)
sd[["Age onset obesity"]]<-sd(pop.analysis$ageOnsetObese, na.rm = T)
sd[["BMI at index"]]<-sd(pop.analysis$bmi_atindex)

sd[["Duration overweight- men"]]<-sd((pop.analysis%>%filter(sexe=="H"))$yearsBMI25more)
sd[["Duration obesity- men"]]<-sd((pop.analysis%>%filter(sexe=="H"))$yearsBMI30more)
sd[["Cumulative overweight- men"]]<-sd((pop.analysis%>%filter(sexe=="H"))$cumOver)
sd[["Cumulative obesity- men"]]<-sd((pop.analysis%>%filter(sexe=="H"))$cumObese) 
sd[["Age onset overweight- men"]]<-sd((pop.analysis%>%filter(sexe=="H"))$ageOnsetOver, na.rm = T)
sd[["Age onset obesity- men"]]<-sd((pop.analysis%>%filter(sexe=="H"))$ageOnsetObese, na.rm = T)
sd[["BMI at index- men"]]<-sd((pop.analysis%>%filter(sexe=="H"))$bmi_atindex)

sd[["Duration overweight- women"]]<-sd((pop.analysis%>%filter(sexe=="D"))$yearsBMI25more)
sd[["Duration obesity- women"]]<-sd((pop.analysis%>%filter(sexe=="D"))$yearsBMI30more)
sd[["Cumulative overweight- women"]]<-sd((pop.analysis%>%filter(sexe=="D"))$cumOver)
sd[["Cumulative obesity- women"]]<-sd((pop.analysis%>%filter(sexe=="D"))$cumObese) 
sd[["Age onset overweight- women"]]<-sd((pop.analysis%>%filter(sexe=="D"))$ageOnsetOver, na.rm = T)
sd[["Age onset obesity- women"]]<-sd((pop.analysis%>%filter(sexe=="D"))$ageOnsetObese, na.rm = T)
sd[["BMI at index- women"]]<-sd((pop.analysis%>%filter(sexe=="D"))$bmi_atindex)

sd[["Duration overweight- breast pre"]]<-sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$yearsBMI25more)
sd[["Duration obesity- breast pre"]]<-sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$yearsBMI30more)
sd[["Cumulative overweight- breast pre"]]<-sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$cumOver)
sd[["Cumulative obesity- breast pre"]]<-sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$cumObese) 
sd[["Age onset overweight- breast pre"]]<-sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$ageOnsetOver, na.rm = T)
sd[["Age onset obesity- breast pre"]]<-sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$ageOnsetObese, na.rm = T)
sd[["BMI at index- breast pre"]]<-sd((pop.analysis%>%filter(sexe=="D",meno_status_index==0))$bmi_atindex)

sd[["Duration overweight- breast post"]]<-sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$yearsBMI25more)
sd[["Duration obesity- breast post"]]<-sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$yearsBMI30more)
sd[["Cumulative overweight- breast post"]]<-sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$cumOver)
sd[["Cumulative obesity- breast post"]]<-sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$cumObese) 
sd[["Age onset overweight- breast post"]]<-sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$ageOnsetOver, na.rm = T)
sd[["Age onset obesity- breast post"]]<-sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$ageOnsetObese, na.rm = T)
sd[["BMI at index- breast post"]]<-sd((pop.analysis%>%filter(sexe=="D",meno_status_index==1))$bmi_atindex)


tableSD <- plyr::ldply(sd)
tableSD <- tableSD %>% mutate(V1=round(V1, 2))
tableSD <- tableSD %>% separate(col=".id", into=c(".id", "group"), sep="-")
tableSD <- tableSD %>% pivot_wider(names_from="group", values_from="V1")
tableSD %>% htmlTable::htmlTable()

# ajustar age of onset por duration ----
exposures <- c("ageOnsetOver")
,
"ageOnsetObese" 
)
canceres <- c(#"cancer_Non-Hodgkin Lymphoma" ,
  # "cancer_Colorectal"#,
  #  "cancer_Leukemia"    ,
  #  "cancer_Liver" ,
  #  "cancer_Brain and CNS"  ,
  #  "cancer_Pancreas"       ,
  #  "cancer_Trachea, bronchus & Lung"   ,
  #  "cancer_Ovary"            ,
  #  "cancer_Head and Neck"     ,
  #  #"cancer_Others and non-specific"  ,
  #  "cancer_Cervix Uteri" ,
  #  "cancer_Kidney"       ,
  #  "cancer_Malignant melanoma of skin" ,
  #  "cancer_Bone and articular cartilage" ,
  #  "cancer_Testis"        ,
  #  "cancer_Stomach"     ,
  #  "cancer_Thyroid"   ,
  "cancer_Corpus Uteri"#,
  #  "cancer_Esophagus"   ,
  #  "cancer_Prostate"  ,
  #  "cancer_Bladder"          ,
  # # "cancer_Penis"                ,
  #  "cancer_Hodgkin lymphoma"    ,
  #  "cancer_Larynx"                 ,
  #  #"cancer_Connective and soft tissue"  ,
  #  "cancer_Multiple myeloma"  ,
  #  "cancer_Gallbladder & biliary tract",
  #  #"cancer_Urinary tract"  ,
  #  "cancer_Breast pre"  ,
  #  "cancer_Breast post"
)



pop.analysis<-pop.analysis %>% mutate(durationOver.gr=ifelse(yearsBMI25more==23, 1, 0),
                                      durationOver.gr=as.factor(durationOver.gr))

#change now its restricted to overweight
out <- list()
harrell<-list()
i = 0
current.pop<-pop.analysis%>% group_by(.imp)
for (cancer in canceres){
  print(paste0("Working on cancer", cancer))
  
  if(cancer =="cancer_Breast post"){
    current.pop<-pop.analysis %>% filter(meno_status_index==1) %>% group_by(.imp)
    print("breast post")
  } else {
    if(cancer=="cancer_Breast pre"){
      current.pop<-pop.analysis %>% filter(meno_status_index==0)%>% group_by(.imp)
      print("breast pre")
    }
    else{print("normal")}
  }
  
  
  for (exposure in exposures){
    print(paste0("Working on exposure", exposure))  
    
    print(paste0("Working on fully adjusted")) 
    
    if(cancer %in%c("cancer_Breast post","cancer_Breast pre" )){
      if(exposure == "ageOnsetOver"){
        test <- current.pop %>%  
          do(model = coxph(Surv(age.index, age.index + follow_up, get(cancer)) ~
                             get(exposure) + qmedea + nationality + smoking +
                             alcohol+yearsBMI25more +strata(age_gr),  data = .))  
      }else{
        if(exposure == "ageOnsetObese"){
          test <- current.pop %>%  
            do(model = coxph(Surv(age.index, age.index + follow_up, get(cancer)) ~
                               get(exposure) + qmedea + nationality + smoking +
                               alcohol+yearsBMI30more +strata(age_gr),  data = .))  
        }
        else{
          test <- current.pop %>%  
            do(model = coxph(Surv(age.index, age.index + follow_up, get(cancer)) ~
                               get(exposure) + qmedea + nationality + smoking +
                               alcohol +strata(age_gr),  data = .))  
        }
      }
      
    }else{
      if(exposure == "ageOnsetOver"){
        test <- current.pop %>%  
          do(model = coxph(Surv(age.index, age.index + follow_up, get(cancer)) ~
                             get(exposure) +sexe+ qmedea + nationality + smoking +
                             alcohol+rcs(yearsBMI25more,3) +strata(age_gr),  data = .))  
      }else{
        if(exposure == "ageOnsetObese"){
          test <- current.pop %>%  
            do(model = coxph(Surv(age.index, age.index + follow_up, get(cancer)) ~
                               get(exposure) +sexe+ qmedea + nationality + smoking +
                               alcohol +yearsBMI30more+strata(age_gr),  data = .))  
        }
        else{
          test <- current.pop %>%  
            do(model = coxph(Surv(age.index, age.index + follow_up, get(cancer)) ~
                               get(exposure) +sexe+ qmedea + nationality + smoking +
                               alcohol +strata(age_gr),  data = .))  
        }
      }  
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
    harrell[[i]]<-summarize_if(c_index,is.numeric, mean) %>% mutate(exposure=exposure, cancer=cancer, group="no",
                                                                    model="Fully adjusted")
    
    
    
    
  }
  
  
  
}

df <- plyr::ldply (out, tibble)
char<-plyr::ldply (harrell, tibble)


saveRDS(char, "charlson_duration.rds")

saveRDS(df, "out.95_duration.rds")


# histograma duration --------
ggplot(pop.analysis %>% filter(everOver==1))+
  geom_histogram(aes(yearsBMI25more), bins=23)

ggplot(pop.analysis %>% filter(everObese==1))+
  geom_histogram(aes(yearsBMI30more), bins=23)


# harrell index with bmi + exposures --------
exposures <- c("yearsBMI25more",
               "yearsBMI30more",
               "ageOnsetOver",
               "ageOnsetObese" ,
               "cumOver" ,
               "cumObese" 
)
out <- list()
harrell<-list()
i = 0
current.pop<-pop.analysis%>% group_by(.imp)
for (cancer in canceres){
  print(paste0("Working on cancer", cancer))
  
  if(cancer =="cancer_Breast post"){
    current.pop<-pop.analysis %>% filter(meno_status_index==1) %>% group_by(.imp)
    print("breast post")
  } else {
    if(cancer=="cancer_Breast pre"){
      current.pop<-pop.analysis %>% filter(meno_status_index==0)%>% group_by(.imp)
      print("breast pre")
    }
    else{print("normal")}
  }
  
  
  for (exposure in exposures){
    print(paste0("Working on exposure", exposure))  
    
    print(paste0("Working on fully adjusted")) 
    
    if(cancer %in%c("cancer_Breast post","cancer_Breast pre" )){
      
      test <- current.pop %>%  
        do(model = coxph(Surv(age.index, age.index + follow_up, get(cancer)) ~
                           get(exposure)+bmi_atindex + qmedea + nationality + smoking +
                           alcohol +strata(age_gr),  data = .))  
      
      
      
    }else{
      
      test <- current.pop %>%  
        do(model = coxph(Surv(age.index, age.index + follow_up, get(cancer)) ~
                           get(exposure)+bmi_atindex +sexe+ qmedea + nationality + smoking +
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
    harrell[[i]]<-summarize_if(c_index,is.numeric, mean) %>% mutate(exposure=exposure, cancer=cancer, group="no",
                                                                    model="Fully adjusted")
    
  }
  
 
  print("working on BMI at index model")
  if(cancer %in%c("cancer_Breast post","cancer_Breast pre" )){
    
    test <- current.pop %>%  
      do(model = coxph(Surv(age.index, age.index + follow_up, get(cancer)) ~
                         bmi_atindex + qmedea + nationality + smoking +
                         alcohol +strata(age_gr),  data = .))  
  }else{
    
    test <- current.pop %>%  
      do(model = coxph(Surv(age.index, age.index + follow_up, get(cancer)) ~
                         bmi_atindex +sexe+ qmedea + nationality + smoking +
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
           model="BMI at index", 
           aic1=aic[1], aic2=aic[2], aic3=aic[3], aic4=aic[4],aic5=aic[5],
           bic1=bic[1], bic2=bic[2], bic3=bic[3], bic4=bic[4],bic5=bic[5])
  
  c_index<-tibble()
  for (j in 1:5){
    c_index <- rbind(c_index,(test$model[[j]])$concordance %>% bind_rows()) 
  }
  harrell[[i]]<-summarize_if(c_index,is.numeric, mean) %>% mutate(exposure=exposure, cancer=cancer, group="no",
                                                                  model="BMI at index")
  
  
  
}

df <- plyr::ldply (out, tibble)
char<-plyr::ldply (harrell, tibble)

#saveRDS(df, "outPlusBMI.rds")
saveRDS(char, "HarellPlusBMI_last.rds")

# plot charlson -------
char <- readRDS("HarellPlusBMI.rds")
h<-char %>% filter(model%in%c("Null model","BMI at index","Fully adjusted" ))

table <- h %>% select(cancer, model, concordance, exposure, std)
table <- table %>% mutate(exposure=ifelse(model!="Fully adjusted", model, exposure)) 
table<-table %>% select(-model) %>% distinct()
table <- table%>% 
  pivot_wider(names_from=exposure, values_from=c(concordance, std))
table <- table %>% separate(cancer, into=c(NA, "cancer"), sep="_")%>% 
  mutate(cancer=recode(cancer, 
                       "Breast pre"="Breast premenopausal",
                       "Breast post"="Breast postmenopausal"
  ))
table[match(order.cancer1, table$cancer),] %>% htmlTable::htmlTable()

library(kableExtra)

fit<-table
#res<-matrix(nrow=30, ncol=10)
fit <- fit %>%
  mutate_at(2:ncol(fit), funs(as.numeric(str_replace_all(., ",", "")))) %>% 
  mutate_at(2:ncol(fit), funs(round(.,5)))
# for (i in 2:nrow(fit)) {
#   val.min <- min(fit[i,2:ncol(fit)])
#   j<-which(fit[i, 2:ncol(fit)]==val.min)
#   res[i,j]<- cell_spec(fit[i,j], background = "orange")
# }
# 
# kable(escape = F, res) %>% 
#   kable_styling(bootstrap_options = c("striped", "bordered"))



df <- pivot_longer(fit, cols=2:length(fit),names_to=c("tipo","variable"), values_to=c("value"), names_sep="_") 
df<-df %>% pivot_wider( names_from = "tipo", values_from="value")
df <- df %>% rename(value=concordance, sd=std)
df <- df %>% left_join(df %>% filter(variable=="BMI at index") %>% select(bmi.index=value, cancer) )
df <- df %>% mutate(diff=value-bmi.index)
# df <- df %>% group_by(cancer) %>% arrange(value) %>% mutate(val_norm=row_number())%>%
#   ungroup()  %>% 
#   filter(cancer%in%cancer.include) %>% 
#   mutate(cancer=factor(cancer, levels=rev(cancer.include))) 

df<-df %>% mutate(char=paste0(format(round(value,3), nsmall=3),
                              "\n(", 
                              format(round(value-1.96*sd, 3), nsmall=3),
                              "-", 
                              format(round(value+1.96*sd, 3), nsmall=3),
                              ")")) 
%>% 
  mutate(val_norm=ifelse(variable%in%c("Null model", "BMI at index"), NA, val_norm))

order.cancer2<-str_replace_all(order.cancer1, "\n", "")

ggplot(data = df %>% filter(!variable%in%c("Null model"))%>% 
         mutate(diff=ifelse(variable%in%c("Null model", "BMI at index"), NA, diff),
                cancer=ifelse(cancer=="Non-Hodgkin Lymphoma" , "Non-Hodgkin Lymph.", cancer ))     , 
       aes(x = factor(variable, levels=c("BMI at index",
                                         "yearsBMI25more",  "yearsBMI30more",
                                         "cumOver" , "cumObese"  ,
                                         "ageOnsetOver" ,   "ageOnsetObese", 
                                         "trajectorySlope")), 
           y = factor(cancer, levels=order.cancer2),
           group=cancer),
       position_dodge()) +
  geom_tile(aes(fill = diff), colour = "white")+ 
  scale_fill_gradient2(low="tomato", high="darkblue",na.value="white", name="Difference") +     
  geom_text(aes(label=char), lineheight = .7)+
  theme_minimal()+
  theme(legend.position="right", 
        text=element_text(size = 16),
        axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 8),
        panel.grid = element_blank())+
  xlab("")+
  ylab("")+
  scale_x_discrete(labels = (c("BMI at index date",
                               expression("+duration BMI">="25"),
                               expression("+duration BMI">="30"),
                               expression("+cum. exposure BMI">="25"),
                               expression("+cum. exposure BMI">="30"),
                               expression("+age onset BMI">="25"),
                               expression("+age onset BMI">="30"))))

ggsave("harrellindexBMITest2_colorblind.png",
       dpi=300,
       width = 16, height = 11)

# sensitivity - follow up starting 1 y after ----
exposures <- c("yearsBMI25more",
               "yearsBMI30more",
               "ageOnsetOver",
               "ageOnsetObese" ,
               "cumOver" ,
               "cumObese" 
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

#change age index, 1 year later
index.date <- ymd("2010-01-01")

# adults at index date
pop.analysis <- pop.analysis %>% mutate(age.index=fun.age(dnaix, index.date))

#change now its restricted to overweight
out <- list()
harrell<-list()
i = 0
current.pop<-pop.analysis%>% group_by(.imp)
for (cancer in canceres){
  print(paste0("Working on cancer", cancer))
  
  if(cancer =="cancer_Breast post"){
    current.pop<-pop.analysis %>% filter(meno_status_index==1) %>% group_by(.imp)
    print("breast post")
  } else {
    if(cancer=="cancer_Breast pre"){
      current.pop<-pop.analysis %>% filter(meno_status_index==0)%>% group_by(.imp)
      print("breast pre")
    }
    else{print("normal")}
  }
  
  
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
                           get(exposure) +sexe+ qmedea + nationality + smoking +
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
    harrell[[i]]<-summarize_if(c_index,is.numeric, mean) %>% mutate(exposure=exposure, cancer=cancer, group="no",
                                                                    model="Fully adjusted")
    
    
  
    
  }
  
  print(paste0("Working on exposure BMI at index")) 
  
  if(cancer %in%c("cancer_Breast post","cancer_Breast pre" )){
    test <- current.pop %>%  
      do(model = coxph(Surv(age.index, age.index + follow_up, get(cancer)) ~
                         bmi_atindex + qmedea + nationality + smoking +
                         alcohol +strata(age_gr),  data = .))   
  }else{
    test <- current.pop %>%  
      do(model = coxph(Surv(age.index, age.index + follow_up, get(cancer)) ~
                         bmi_atindex +sexe + qmedea + nationality + smoking +
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
    mutate(exposure=exposure, cancer=cancer, group="no", model="BMI at index", 
           aic1=aic[1], aic2=aic[2], aic3=aic[3], aic4=aic[4],aic5=aic[5],
           bic1=bic[1], bic2=bic[2], bic3=bic[3], bic4=bic[4],bic5=bic[5])
  
  c_index<-tibble()
  for (j in 1:5){
    c_index <- rbind(c_index,(test$model[[j]])$concordance %>% bind_rows()) 
  }
  harrell[[i]]<-summarize_if(c_index,is.numeric, mean) %>%
    mutate(exposure=exposure, cancer=cancer, group="no", model="BMI at index")
  
  
}

df <- plyr::ldply (out, tibble)
char<-plyr::ldply (harrell, tibble)

saveRDS(df, "out.951Y.rds")


