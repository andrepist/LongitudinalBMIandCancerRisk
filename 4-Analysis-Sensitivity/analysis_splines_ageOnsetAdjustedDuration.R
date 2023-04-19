#analysis
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
pop.analysis<-pop.analysis %>% left_join(bmiindex)
# age of onset with duration and bmi at index  ---------
exposures <- c("ageOnsetOver",
               "ageOnsetObese" ,
               "bmi_atindex"
               
)
canceres <- c("cancer_Non-Hodgkin Lymphoma" ,
              "cancer_Colorectal",
              "cancer_Leukemia"    ,
              "cancer_Liver" ,
              "cancer_Brain and CNS"  ,
              "cancer_Pancreas"       ,
              "cancer_Trachea, bronchus & Lung"   ,
              "cancer_Ovary"            ,
              "cancer_Head and Neck"     ,
              #"cancer_Others and non-specific"  ,
              "cancer_Cervix Uteri" ,
              "cancer_Kidney"       ,
              "cancer_Malignant melanoma of skin" ,
              "cancer_Bone and articular cartilage" ,
              "cancer_Testis"        ,
              "cancer_Stomach"     ,
              "cancer_Thyroid"   ,
              "cancer_Corpus Uteri",
              "cancer_Esophagus"   ,
              "cancer_Prostate"  ,
              "cancer_Bladder"          ,
              # "cancer_Penis"                ,
              "cancer_Hodgkin lymphoma"    ,
              "cancer_Larynx"                 ,
              #"cancer_Connective and soft tissue"  ,
              "cancer_Multiple myeloma"  ,
              "cancer_Gallbladder & biliary tract",
              #"cancer_Urinary tract"  ,
              "cancer_Breast pre"  ,
              "cancer_Breast post"
)


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
                    "Breast pre",
                    "Breast post",
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


# splines --------
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
# BREAST HA DE SER ULTIMA

working.data<-pop.analysis
dd<<-datadist(working.data); options(datadist = "dd" )
output <- tibble()
for (cancer in canceres){
  if(cancer =="cancer_Breast post"){
    working.data<-pop.analysis %>% filter(meno_status_index==1) 
    dd<<-datadist(working.data); options(datadist = "dd" )
    print("breast post")
  } else {
    if(cancer=="cancer_Breast pre"){
      working.data<-pop.analysis %>% filter(meno_status_index==0)
      dd<<-datadist(working.data); options(datadist = "dd" )
      print("breast pre")
    }
    else{
      print("normal")}
  }
  
  print(paste0("Working on cancer ", cancer))
  exp=tibble()
  for (exposure in exposures){
    print(paste0("Working on exposure ", exposure)) 
    sequence <-  seq(from=min(working.data %>% select(!!!exposure), na.rm = TRUE),
                     to=max(working.data %>% select(!!!exposure), na.rm = TRUE) ,length.out=15 )
    imput=tibble()
    for (i in 1:5){
      print(paste0("Working on imputation ", i))
      if(cancer %in%c("cancer_Breast post","cancer_Breast pre" )){
        if(exposure == "ageOnsetOver"){
          iformula <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,"`)~rcs(",exposure,
                                     ",3)+yearsBMI25more+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
        }else{
          if(exposure == "ageOnsetObese"){
            iformula <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,"`)~rcs(",exposure,
                                       ",3)+yearsBMI30more+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
          }else{
            iformula <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,"`)~rcs(",exposure,
                                       ",3)+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
          
          }
        }
       
      }else{
        if(exposure == "ageOnsetOver"){
          iformula <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,"`)~rcs(",exposure,
                                     ",3)+sexe+yearsBMI25more+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
        }else{
          if(exposure == "ageOnsetObese"){
            iformula <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,"`)~rcs(",exposure,
                                       ",3)+sexe+yearsBMI30more+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
          }
          else{
            iformula <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,"`)~rcs(",exposure,
                                       ",3)+sexe+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
            
        }
        }
      }
      
      fit<- cph(iformula,  subset=(.imp==i),
                surv=TRUE,x=TRUE,y=TRUE, data = working.data)
      a <- tibble()
      for (val in sequence){
        summ <- paste0("as.data.frame(summary(fit,", exposure, "=c(0,val), antilog=FALSE))")
        a <- rbind(a,head(eval(parse(text=summ)),1) %>% 
                     rename(estimate=Effect) %>% mutate(value=val))
      }
      imput <- rbind(imput, a %>% mutate(imp=i))
    }
    exp <- rbind(exp, imput %>% mutate(expo=exposure))
  }
  output <- rbind(output, exp %>% mutate(cancer=cancer))
}
#In regularize.values(x, y, ties, missing(ties), na.rm = na.rm) :collapsing to unique 'x' values
saveRDS(output, "outSplines_ageOnsetDurationANDBmiindex.rds")
