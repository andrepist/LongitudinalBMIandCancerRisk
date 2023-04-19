library(tidyverse)
library(lubridate)
library(tableone)


# pop.exp <- readRDS("popExposiciones.rds")
# 
# pop.exp <- pop.exp %>% select(-c("follow_up", "age.index", "group"))
# 
# pop <- pop.exp %>% 
#   inner_join(pop.gen, by = c("sexe", "dnaix", "entrada", "end", "dat.fake",
#                             "situacio", "sortida", "agr_pais", "id" ))

pop.gen <- readRDS("analysis.pop.rds")
pop <- readRDS("popExposicionesFINAL.rds")


pop <- test4
rm(test, test1, test2, test4, wide)
gc()
dataBMI_indexdate <- readRDS("E:/Andrea/Trajectories OBECAN/dataBMI_indexdate_new.rds")
dataBMI_indexdate <- dataBMI_indexdate %>% distinct(id,index_ch, ruralitat, nationality)
pop <- pop %>% left_join(dataBMI_indexdate)
pop <- pop %>% left_join(pop.gen %>% 
                           select(id, follow_up, group, end, dat,sortida, situacio, age.index, sexe, dnaix ))

pop <- pop %>% mutate(cause.exit=ifelse(end==sortida, situacio, "Any cancer"))
pop <- pop %>% mutate(cancer.outcome=ifelse(is.na(group)|group%in%c("Urinary tract",
                                                        "Penis",
                                                        "Connective and soft tissue",
                                                        "Others and non-specific"), 0, 1 ))

bmi.index <- pop %>% pivot_longer(cols=paste0("bmi_", c(40,55,70,118)))
bmi.index1 <- bmi.index %>% mutate(bmi.age=as.numeric(substr(name, 5,7)))
bmi.index1 <- bmi.index1 %>% filter(bmi.age<=age.index) #CHANGE THIS -----
bmi.index1 <- bmi.index1 %>% arrange(abs(bmi.age-age.index)) %>% 
  distinct(id, .imp, .keep_all=TRUE)
bmi.index1 <- bmi.index1 %>% transmute(id, .imp, bmi.index=value)

pop <- pop %>% left_join(bmi.index1)
rm(bmi.index1, bmi.index)

pop <- pop %>%select(-dnaix) %>%  rbind(populationCox)

saveRDS(pop, "populationCoxReviewer.rds")

pop<-pop %>% ungroup() 
rm(pp)
# tabla descriptiva ------

# Tabla 1 - Baseline characteristics of the study population by overweight and obesity during adulthood, after multiple imputations --------

vars<-c("follow_up",
        "yearsBMI25more",
        "yearsBMI30more",
        "cumOver",
        "cumObese",
        "trajectorySlope",
        "ageOnsetOver",
        "ageOnsetObese",
        "bmi.index",
        "age.index",
        "sexe",
        "nationality",
        "qmedea",
        "smoking",
        "alcohol",
        "index_ch",
        "cause.exit",
        "cancer.outcome"
        
)
factor.vars<- c( "sexe",
                 "nationality",
                 "qmedea",
                 "smoking",
                 "alcohol",
                 "index_ch",
                 "cause.exit",
                 "cancer.outcome"
                 )

summary<-cbind(
  
  print(CreateTableOne(
    vars =  vars,
    factorVars = factor.vars,
    includeNA=T,
    data = pop ,
    #strata = "",
    test = F), 
    showAllLevels=F,smd=F,
    nonnormal = vars, 
    noSpaces = TRUE,
    contDigits = 1,
    printToggle=FALSE),
  
  print(CreateTableOne(
    vars =  vars,
    factorVars = factor.vars,
    includeNA=T,
    data = pop ,
    strata = "everOver",
    test = F), 
    showAllLevels=F,smd=F,
    nonnormal = vars, 
    noSpaces = TRUE,
    contDigits = 1,
    printToggle=FALSE),

print(CreateTableOne(
    vars =  vars,
    factorVars = factor.vars,
    includeNA=T,
    data = pop ,
    strata = "everObese",
    test = F), 
    showAllLevels=F,smd=F,
    nonnormal = vars, 
    noSpaces = TRUE,
    contDigits = 1,
    printToggle=FALSE)
)


summary.characteristics <- summary
# divide for 5 
nice.num<-function(x){
  prettyNum(x, big.mark=",", nsmall = 0, digits=0,scientific = FALSE)}
for(i in 1:ncol(summary.characteristics)) {
  # tidy up 
  cur_column <- summary.characteristics[, i]
  cur_column <- str_extract(cur_column, '[0-9.]+\\b') %>% 
    as.numeric() 
  cur_column <-nice.num(cur_column/5)
  # add back in
  summary.characteristics[, i] <- str_replace(string=summary.characteristics[, i], 
                                              pattern='[0-9.]+\\b', 
                                              replacement=cur_column)    
}

library(htmlTable)


summary.characteristics %>% htmlTable()
summary %>% htmlTable()

# extra --------
# categories of overweight duration

pop <- pop %>% mutate(overDuration=ifelse(yearsBMI25more==0 , "0", 
                                          ifelse(yearsBMI25more!=0&yearsBMI25more<=11, "<=11",
                                                 ifelse(yearsBMI25more!=0&yearsBMI25more>11, ">11",
                                                        "Missing"))),
                      obeseDuration=ifelse(yearsBMI30more==0 , "0", 
                                           ifelse(yearsBMI30more!=0&yearsBMI30more<=11, "<=11",
                                                  ifelse(yearsBMI30more!=0&yearsBMI30more>11, ">11",
                                                         "Missing"))))
# mira distrbucion de exposiciones
# Exposure 2 (cumulative exposure),  0, <100 units, >100 units 
summary(pop$cumOver)
summary(pop$cumObese)
pop <- pop %>% mutate(cumOver.gr=ifelse(cumOver==0, "0",
                                        ifelse(cumOver!=0&cumOver<=100,"<=100",
                                               ifelse(cumOver!=0&cumOver>100,">100",
                                                      "Missing"))),
                      cumObese.gr=ifelse(cumObese==0, "0",
                                         ifelse(cumObese!=0&cumObese<=100,"<=100",
                                                ifelse(cumObese!=0&cumObese>100,">100",
                                                       "Missing")))
)

# exposure 3
summary(pop$trajectorySlope)
pop <- pop %>% mutate(slope.gr=ifelse(trajectorySlope<=0.015, "<=0.015", 
                                      ifelse(trajectorySlope>0.015, "0.015",
                                             "Missing")))

#exposure 4
summary(pop$ageOnsetObese)
summary(pop$ageOnsetOver)

pop <- pop %>% mutate(ageOnsetOver.gr=ifelse(!is.na(ageOnsetOver)&ageOnsetOver<=30, "<=30", 
                                             ifelse(!is.na(ageOnsetOver)&ageOnsetOver>30, "30",
                                                    "Missing")),
                      ageOnsetObese.gr=ifelse(ageOnsetObese<=30, "<=30", 
                                              ifelse(ageOnsetObese>30, "30",
                                                     "Missing"))
)

# Table S2. Baseline characteristics of the study population by type or lack of BMI assessment, after multiple imputations ------
mi.pop <- dataBMI_indexdate %>% filter(id%in%pop$id)
table(mi.pop$qmedea, useNA="always")
mi.pop <- mi.pop %>%  mutate(qmedea=ifelse(qmedea=="missing", NA_character_, qmedea))
index.missingCov <- mi.pop %>% filter(is.na(qmedea)|is.na(smoking)|is.na(alcohol))
index.oneBMI <- mi.pop %>% filter(!is.na(bmi))

vars<-c("follow_up",
        "yearsBMI25more",
        "yearsBMI30more",
        "cumOver",
        "cumObese",
        "trajectorySlope",
        "ageOnsetOver",
        "ageOnsetObese",
        "bmi.index",
        "age.index",
        "sexe",
        "nationality",
        "qmedea",
        "smoking",
        "alcohol",
        "index_ch",
        "cause.exit",
        "cancer.outcome"
        
)
factor.vars<- c( "sexe",
                 "nationality",
                 "qmedea",
                 "smoking",
                 "alcohol",
                 "index_ch",
                 "cause.exit",
                 "cancer.outcome"
)

summary<-cbind(
  
  print(CreateTableOne(
    vars =  vars,
    factorVars = factor.vars,
    includeNA=T,
    data = pop ,
    test = F), 
    showAllLevels=F,smd=F,
    nonnormal = vars, 
    noSpaces = TRUE,
    contDigits = 1,
    printToggle=FALSE),
  
  print(CreateTableOne(
    vars =  vars,
    factorVars = factor.vars,
    includeNA=T,
    data = pop %>% filter(!id%in%index.oneBMI$id) ,
    test = F), 
    showAllLevels=F,smd=F,
    nonnormal = vars, 
    noSpaces = TRUE,
    contDigits = 1,
    printToggle=FALSE),
  
  print(CreateTableOne(
    vars =  vars,
    factorVars = factor.vars,
    includeNA=T,
    data = pop %>% filter(id%in%index.oneBMI$id) ,
    test = F), 
    showAllLevels=F,smd=F,
    nonnormal = vars, 
    noSpaces = TRUE,
    contDigits = 1,
    printToggle=FALSE)
  # ,
  # 
  # print(CreateTableOne(
  #   vars =  vars,
  #   factorVars = factor.vars,
  #   includeNA=T,
  #   data = pop %>% filter(!id%in%index.oneBMI$id) ,
  #   test = F), 
  #   showAllLevels=F,smd=F,
  #   nonnormal = vars, 
  #   noSpaces = TRUE,
  #   contDigits = 1,
  #   printToggle=FALSE)
)

summary.characteristics <- summary
# divide for 5 
nice.num<-function(x){
  prettyNum(x/5, big.mark=",", nsmall = 0, digits=0,scientific = FALSE)}
for(i in 1:ncol(summary.characteristics)) {
  # tidy up 
  cur_column <- summary.characteristics[, i]
  cur_column <- str_extract(cur_column, '[0-9.]+\\b') %>% 
    as.numeric() 
  cur_column <-nice.num(cur_column)
  # add back in
  summary.characteristics[, i] <- str_replace(string=summary.characteristics[, i], 
                                              pattern='[0-9.]+\\b', 
                                              replacement=cur_column)    
}



summary.characteristics %>% htmlTable()
summary %>% htmlTable()

# Table S3. Baseline characteristics of the study population by different exposures, after multiple imputations --------

vars<-c("follow_up",
        "yearsBMI25more",
        "yearsBMI30more",
        "cumOver",
        "cumObese",
        "trajectorySlope",
        "ageOnsetOver",
        "ageOnsetObese",
        "bmi.index",
        "age.index",
        "sexe",
        "nationality",
        "qmedea",
        "smoking",
        "alcohol",
        "index_ch",
        "cause.exit",
        "cancer.outcome"
        
)
factor.vars<- c( "sexe",
                 "nationality",
                 "qmedea",
                 "smoking",
                 "alcohol",
                 "index_ch",
                 "cause.exit",
                 "cancer.outcome"
)

summary<-cbind(
  
  print(CreateTableOne(
    vars =  vars,
    factorVars = factor.vars,
    includeNA=T,
    data = pop ,
    strata = "overDuration",
    test = F), 
    showAllLevels=F,smd=F,
    nonnormal = vars, 
    noSpaces = TRUE,
    contDigits = 1,
    printToggle=FALSE),
  
  print(CreateTableOne(
    vars =  vars,
    factorVars = factor.vars,
    includeNA=T,
    data = pop ,
    strata = "cumOver.gr",
    test = F), 
    showAllLevels=F,smd=F,
    nonnormal = vars, 
    noSpaces = TRUE,
    contDigits = 1,
    printToggle=FALSE),
  
  print(CreateTableOne(
    vars =  vars,
    factorVars = factor.vars,
    includeNA=T,
    data = pop ,
    strata = "ageOnsetOver.gr",
    test = F), 
    showAllLevels=F,smd=F,
    nonnormal = vars, 
    noSpaces = TRUE,
    contDigits = 1,
    printToggle=FALSE),
  
  print(CreateTableOne(
    vars =  vars,
    factorVars = factor.vars,
    includeNA=T,
    data = pop ,
    strata = "slope.gr",
    test = F), 
    showAllLevels=F,smd=F,
    nonnormal = vars, 
    noSpaces = TRUE,
    contDigits = 1,
    printToggle=FALSE)
)


summary.characteristics <- summary
# divide for 5 
nice.num<-function(x){
  prettyNum(x/5, big.mark=",", nsmall = 0, digits=0,scientific = FALSE)}
for(i in 1:ncol(summary.characteristics)) {
  # tidy up 
  cur_column <- summary.characteristics[, i]
  cur_column <- str_extract(cur_column, '[0-9.]+\\b') %>% 
    as.numeric() 
  cur_column <-nice.num(cur_column)
  # add back in
  summary.characteristics[, i] <- str_replace(string=summary.characteristics[, i], 
                                              pattern='[0-9.]+\\b', 
                                              replacement=cur_column)    
}



summary.characteristics %>% htmlTable()
summary %>% htmlTable()

# Histograms exposures -----------

hist(populationCox$trajectorySlope, main="Slope")
hist(pop$ageOnsetOver, main="Age onset overweight")
hist(pop$ageOnsetObese, main="Age onset obesity")
hist(populationCox$cumOver, main="Cumulative overweight" )
hist(populationCox$cumObese, main="Cumulative obesity" )
hist(populationCox$yearsBMI25more, main="Duration overweight" )
hist(populationCox$yearsBMI30more, main="Duration obesity" )

imp.set <- complete(imp, "long", include=FALSE)

summary.characteristics %>% htmlTable()
summary %>% htmlTable()

populationCox <- populationCox %>% mutate(when=ifelse(bmi.index==bmi_40, 40,
                                                      ifelse(bmi.index==bmi_55, 55,
                                                             ifelse(bmi.index==bmi_70, 70,
                                                                    ifelse(bmi.index==bmi_95, 95,NA)))))
prop.table(table(populationCox$when))
