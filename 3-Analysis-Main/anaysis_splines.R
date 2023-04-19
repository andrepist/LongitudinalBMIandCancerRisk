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
# missing categories medea and smoking ----- 
# CHANGE THIS -----
#pop.analysis <- pop.analysis %>% mutate(smoking=ifelse(is.na(smoking), "missing", smoking))

# colorectal ----
exposures <- c("yearsBMI25more",
               "yearsBMI30more",
               "ageOnsetOver",
               "ageOnsetObese" ,
               "cumOver" ,
               "cumObese" ,
               "trajectorySlope"
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
              "cancer_Urinary tract" ,
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

# BREAST HA DE SER ULTIMA
working.data<-pop.analysis
dd<<-datadist(working.data); options(datadist = "dd" )
output <- tibble()
for (cancer in canceres){
  if(cancer =="cancer_Breast post"){
    working.data<-pop.analysis %>% filter(meno_status_index==1) 
    print("breast post")
  } else {
    if(cancer=="cancer_Breast pre"){
      working.data<-pop.analysis %>% filter(meno_status_index==0)
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
        iformula <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,"`)~rcs(",exposure,",3)+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
      }else{
        iformula <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,"`)~rcs(",exposure,",3)+sexe+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
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
saveRDS(output, "outSpline_FINAL.rds")

output <- output %>% rbind(output1)



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

#age of onset only -----
output1 <- tibble()
working.data<-pop.analysis

dd<<-datadist(working.data); options(datadist = "dd" )
for (cancer in canceres){
 
  print(paste0("Working on cancer ", cancer))
  exp=tibble()
  for (exposure in c("ageOnsetOver" ,   "ageOnsetObese")){
    print(paste0("Working on exposure ", exposure)) 
    sequence <-  seq(from=min(working.data %>% select(!!!exposure), na.rm = TRUE),
                     to=max(working.data %>% select(!!!exposure), na.rm = TRUE) ,length.out=15 )
    imput=tibble()
    for (i in 1:5){
      print(paste0("Working on imputation ", i))
      if(cancer %in%c("cancer_Breast post","cancer_Breast pre" )){
        iformula <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,"`)~rcs(",exposure,",3)+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
      }else{
        iformula <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,"`)~rcs(",exposure,",3)+sexe+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
      }
      
      fit<- cph(iformula,  subset=(.imp==i),
                surv=TRUE,x=TRUE,y=TRUE, data = working.data)
      a <- tibble()
      for (val in sequence){
        summ <- paste0("as.data.frame(summary(fit,", exposure, "=c(40,val), antilog=FALSE))")
        a <- rbind(a,head(eval(parse(text=summ)),1) %>% 
                     rename(estimate=Effect) %>% mutate(value=val))
      }
      imput <- rbind(imput, a %>% mutate(imp=i))
    }
    exp <- rbind(exp, imput %>% mutate(expo=exposure))
  }
  output1 <- rbind(output1, exp %>% mutate(cancer=cancer))
}
#In regularize.values(x, y, ties, missing(ties), na.rm = na.rm) :collapsing to unique 'x' values
saveRDS(output1, "outSpline_ageOnset.rds")

# spline stratified -------

working.data<-pop.analysis
dd<<-datadist(working.data); options(datadist = "dd" )
output <- tibble()
if(canceres[length(canceres)]=="cancer_Breast post"){print("Ok, continue")}
for (cancer in canceres){
  if(cancer =="cancer_Breast post"){
    working.data<-pop.analysis %>% filter(meno_status_index==1) 
    print("breast post")
  } else {
    if(cancer=="cancer_Breast pre"){
      working.data<-pop.analysis %>% filter(meno_status_index==0)
      print("breast pre")
    }
    else{
      print("normal")}
  }
  
  print(paste0("Working on cancer ", cancer))
  exp=tibble()
  for (exposure in c("trajectorySlope")){
    print(paste0("Working on exposure ", exposure)) 
    sequence <-  seq(from=min(working.data %>% select(!!!exposure), na.rm = TRUE),
                     to=max(working.data %>% select(!!!exposure), na.rm = TRUE) ,length.out=15 )
    imput=tibble()
    for (i in 1:5){
      print(paste0("Working on imputation ", i))
      
      strati <-tibble()
      for (bmigr in c("Overweight", "Obese" )){
        
      if(cancer %in%c("cancer_Breast post","cancer_Breast pre" )){
        iformula <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,"`)~rcs(",exposure,",3)+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
      }else{
        iformula <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,"`)~rcs(",exposure,",3)+sexe+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
      }
      
      fit<- cph(iformula,  subset=(.imp==i&bmi.gr18==bmigr),
                surv=TRUE,x=TRUE,y=TRUE, data = working.data)
      
      a <- tibble()
      for (val in sequence){
        summ <- paste0("as.data.frame(summary(fit,", exposure, "=c(0,val), antilog=FALSE))")
        a <- rbind(a,head(eval(parse(text=summ)),1) %>% 
                     rename(estimate=Effect) %>% mutate(value=val))
      }
      strati <- rbind(strati, a %>% mutate(gr=bmigr))
      }
      imput <- rbind(imput, strati %>% mutate(imp=i))
    }
    exp <- rbind(exp, imput %>% mutate(expo=exposure))
  }
  output <- rbind(output, exp %>% mutate(cancer=cancer))
}
#In regularize.values(x, y, ties, missing(ties), na.rm = na.rm) :collapsing to unique 'x' values
saveRDS(output, "outSplineStrat_FINAL.rds")

test <- output %>% group_by(value, cancer, expo, gr) %>% mutate(HR=exp(mean(estimate)),
                                                            low=exp(mean(`Lower 0.95`)),
                                                            up=exp(mean(`Upper 0.95`))) %>% distinct(HR, expo, value, low, up)

test <- test %>% separate(cancer, into=c(NA, "cancer"), sep="_")

test <- test %>% mutate(exposure=expo) %>% 
  mutate(expo=ifelse(exposure%in%c("yearsBMI25more", "yearsBMI30more"), "Duration",
                     ifelse(exposure%in%c( "cumOver","cumObese" ), "Cumulative",
                            ifelse(exposure%in%c( "ageOnsetOver","ageOnsetObese" ), "Age of onset",
                                   "Slope"))),
         Type=ifelse(exposure%in%c("yearsBMI25more", "cumOver","ageOnsetOver"), "Overweight/Obese",
                     ifelse(exposure%in%c("yearsBMI30more", "cumObese","ageOnsetObese"), "Obese",
                            "Overall")))


test <- test %>% filter(cancer%in%cancer.include)

ggplot(data=test)+
  geom_hline(yintercept = 1, colour = "#000000",
             linetype=2, size=0.8)+
  
  geom_line(aes(value, HR, colour=gr), size=2) +
  geom_ribbon(aes(x=value, y=HR,ymin=low, ymax=up, fill=gr), 
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
  ylab("Relative\nhazard ratio\n") +
  coord_cartesian(ylim = c(0.5, 2))


# p non linearity -------
exposures <- c(
  # "yearsBMI25more",
  #              "yearsBMI30more",
               "ageOnsetOver",
               "ageOnsetObese") #,
               # "cumOver" ,
               # "cumObese" )
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
    
    out<-tibble()
    p <- tibble()
    for (i in 1:5){
      print(paste0("Working on imputation ", i))
      
      if(cancer %in%c("cancer_Breast post","cancer_Breast pre" )){
        print("breast")
       
            formulaNL <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,
                                        "`)~rcs(",exposure,",3)+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
            formulaL <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,
                                       "`)~",exposure,"+qmedea+nationality+smoking+alcohol+strat(age_gr) ")) 
          }else{
        
            print("normal")
            formulaNL <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,
                                        "`)~rcs(",exposure,",3)+sexe+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
            formulaL <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,
                                       "`)~",exposure,"+sexe+qmedea+nationality+smoking+alcohol+strat(age_gr) ")) 
          }
        
      a<-coxph(formulaNL,  subset=(.imp==i), data = working.data)
      b<-coxph(formulaL,  subset=(.imp==i),   data = working.data)
      c<-anova(a,b)
      p<-rbind(p,c %>% select( "P(>|Chi|)") %>% filter(!`P(>|Chi|)`=="") %>% mutate(imp=i))
      
      
    }
    q<-pivot_wider(p, names_from=imp, values_from=`P(>|Chi|)`)
    out <- rbind(out, q)
    
    exp <- rbind(exp, out %>% mutate(expo=exposure))
  }
  output <- rbind(output, exp %>% mutate(cancer=cancer))
}
saveRDS(output, "p_nonlinearity_ageOnsetNoAdj.rds")

# old adjusting age of onset for duration
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
  
    out<-tibble()
    p <- tibble()
    for (i in 1:5){
      print(paste0("Working on imputation ", i))
      
      if(cancer %in%c("cancer_Breast post","cancer_Breast pre" )){
        print("breast")
        if(exposure == "ageOnsetOver"){
          formulaNL <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,
                                      "`)~rcs(",exposure,",3)+yearsBMI25more+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
          formulaL <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,
                                     "`)~",exposure,"+yearsBMI25more+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
          }else{
          if(exposure == "ageOnsetObese"){
            formulaNL <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,
                                        "`)~rcs(",exposure,",3)+yearsBMI30more+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
            formulaL <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,
                                       "`)~",exposure,"+yearsBMI30more+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
          }
          else{
            formulaNL <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,
                                        "`)~rcs(",exposure,",3)+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
            formulaL <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,
                                       "`)~",exposure,"+qmedea+nationality+smoking+alcohol+strat(age_gr) ")) 
          }
        }
      
        
      }else{
        if(exposure == "ageOnsetOver"){
          print("onset over")
          formulaNL <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,
                                      "`)~rcs(",exposure,",3)+sexe+yearsBMI25more+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
          formulaL <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,
                                     "`)~",exposure,"+sexe+yearsBMI25more+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
        }else{
          if(exposure == "ageOnsetObese"){
            print("onset obese")
            formulaNL <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,
                                        "`)~rcs(",exposure,",3)+sexe+yearsBMI30more+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
            formulaL <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,
                                       "`)~",exposure,"+sexe+yearsBMI30more+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
          }
          else{
            print("normal")
            formulaNL <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,
                                        "`)~rcs(",exposure,",3)+sexe+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
            formulaL <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,
                                       "`)~",exposure,"+sexe+qmedea+nationality+smoking+alcohol+strat(age_gr) ")) 
          }
        }
    
      }
      
      a<-coxph(formulaNL,  subset=(.imp==i), data = working.data)
      b<-coxph(formulaL,  subset=(.imp==i),   data = working.data)
      c<-anova(a,b)
      p<-rbind(p,c %>% select( "P(>|Chi|)") %>% filter(!`P(>|Chi|)`=="") %>% mutate(imp=i))
    
    
    }
    q<-pivot_wider(p, names_from=imp, values_from=`P(>|Chi|)`)
    out <- rbind(out, q)
    
    exp <- rbind(exp, out %>% mutate(expo=exposure))
  }
  output <- rbind(output, exp %>% mutate(cancer=cancer))
}

# p non linearity bmi at index-------
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
  
    print(paste0("Working on bmi at index ", exposure)) 
    
    out<-tibble()
    p <- tibble()
    for (i in 1:5){
      print(paste0("Working on imputation ", i))
      
      if(cancer %in%c("cancer_Breast post","cancer_Breast pre" )){
        print("breast")
       
            formulaNL <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,
                                        "`)~rcs(bmi_atindex,3)+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
            formulaL <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,
                                       "`)~bmi_atindex+qmedea+nationality+smoking+alcohol+strat(age_gr) ")) 
          
      }else{
        
            print("normal")
            formulaNL <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,
                                        "`)~rcs(bmi_atindex,3)+sexe+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
            formulaL <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,
                                       "`)~bmi_atindex+sexe+qmedea+nationality+smoking+alcohol+strat(age_gr) ")) 
          
        
      }
      
      a<-coxph(formulaNL,  subset=(.imp==i), data = working.data)
      b<-coxph(formulaL,  subset=(.imp==i),   data = working.data)
      c<-anova(a,b)
      p<-rbind(p,c %>% select( "P(>|Chi|)") %>% filter(!`P(>|Chi|)`=="") %>% mutate(imp=i))
      
      
    }
    q<-pivot_wider(p, names_from=imp, values_from=`P(>|Chi|)`)
    out <- rbind(out, q)
    
    exp <- rbind(exp, out %>% mutate(expo=exposure))
  
  output <- rbind(output, exp %>% mutate(cancer=cancer))
}


#In regularize.values(x, y, ties, missing(ties), na.rm = na.rm) :collapsing to unique 'x' values
saveRDS(output, "p_nonlinearity_bmiindex.rds")



output<-readRDS("p_nonlinearity.rds")
output<-output %>% mutate_if(~is.numeric(.), ~round(.,15))
output<-output %>% mutate(p=(`1`+`2`+`3`+`4`+`5`)/5)
output<-output %>% mutate(sig=ifelse(p<0.05, "*", "0"))

output %>% htmlTable::htmlTable()

# bmi at index reference 22 --------
output1 <- tibble()
working.data<-pop.analysis

dd<<-datadist(working.data); options(datadist = "dd" )
for (cancer in canceres){
  
  print(paste0("Working on cancer ", cancer))
  exp=tibble()
  for (exposure in "bmi_atindex"){
    print(paste0("Working on exposure ", exposure)) 
    sequence <-  seq(from=min(working.data %>% select(!!!exposure), na.rm = TRUE),
                     to=max(working.data %>% select(!!!exposure), na.rm = TRUE) ,length.out=15 )
    imput=tibble()
    for (i in 1:5){
      print(paste0("Working on imputation ", i))
      if(cancer %in%c("cancer_Breast post","cancer_Breast pre" )){
        iformula <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,"`)~rcs(",exposure,",3)+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
      }else{
        iformula <- formula(paste0("Surv(age.index,age.index+follow_up,`", cancer,"`)~rcs(",exposure,",3)+sexe+qmedea+nationality+smoking+alcohol+strat(age_gr) "))
      }
      
      fit<- cph(iformula,  subset=(.imp==i),
                surv=TRUE,x=TRUE,y=TRUE, data = working.data)
      a <- tibble()
      for (val in sequence){
        summ <- paste0("as.data.frame(summary(fit,", exposure, "=c(22,val), antilog=FALSE))")
        a <- rbind(a,head(eval(parse(text=summ)),1) %>% 
                     rename(estimate=Effect) %>% mutate(value=val))
      }
      imput <- rbind(imput, a %>% mutate(imp=i))
    }
    exp <- rbind(exp, imput %>% mutate(expo=exposure))
  }
  output1 <- rbind(output1, exp %>% mutate(cancer=cancer))
}
#In regularize.values(x, y, ties, missing(ties), na.rm = na.rm) :collapsing to unique 'x' values
saveRDS(output1, "outSpline_bmiAtIndex.rds")
