library(tidyverse)
library(lubridate)
library(mice)

fun.age<-function(date1,date2){
  #date1 birth
  output<-year(date2)-year(date1)
  ifelse(month(date2)<month(date1)|(month(date2)==month(date1)&day(date2)<day(date1)),output-1,
         output)
}

# imp <- readRDS("imp_500000id.rds")
# 
# data <- complete(imp, "long", include=FALSE)
# nrow(data %>% distinct(id))
# study population -----
# people at least 40 y old at index date
index.date <- ymd("2009-01-01")

# pop <- data %>% mutate(age.index=fun.age(dnaix, index.date))
# nrow(pop %>% distinct(id))
# pop <- pop %>% filter(age.index>=40)
# nrow(pop %>% distinct(id)) # 289,019
# # *5 = 1,445,095 

analysis.pop <- readRDS("E:/Andrea/Trajectories OBECAN/analysis.pop.rds")
trajectories_nosmokingmedea <- readRDS("E:/Andrea/Trajectories OBECAN/trajectories_nosmokingmedea.rds")
pop <- trajectories_nosmokingmedea %>% filter(id%in%analysis.pop$id)

rm(analysis.pop)
rm(trajectories_nosmokingmedea)

#summary(pop$age.index)
# cambio de formato de la base de datos
# me quedo con smoking, medea más cercanos a index date
# wide <- pop %>% select(-c(idp, idup, x1, x2, x3, x4, x5, x6, x7, x8, .id,  "i.cvd","i.htn" ,"i.dm" ,"i.cancer1yprior" , year.bmi, date.bmi)) %>% 
#   filter(fake=="yes") %>% 
#   arrange(abs(age.index-age.bmi)) %>% 
#   distinct(id, .imp, age.bmi, .keep_all=TRUE)

# 15 60 limits
wide <- pop %>% filter(fake=="yes") %>% 
  mutate(bmi=ifelse(bmi<15, 15, 
                    ifelse(bmi>60, 60, bmi)))

traj <- readRDS("imputedforstudypop.rds")

wide <- traj %>% 
  mutate(bmi=ifelse(bmi<15, 15, 
                    ifelse(bmi>60, 60, bmi)))
rm(traj)

# I select the first medea, alcohol, smoking
wide1 <- wide %>% select(-smoking, -alcohol, -qmedea) %>% distinct() %>% 
  left_join(wide %>% arrange(age.bmi) %>% distinct(id, .imp, .keep_all=TRUE) %>% select(-bmi, -age.bmi))

timepoints <- distinct(wide, age.bmi) 

wide2 <- wide1 %>%  
  pivot_wider(names_from = c("age.bmi"), values_from="bmi", 
              names_glue = "bmi_{age.bmi}")

wide <- wide2 %>% group_by(id, .imp) %>% mutate(n=row_number())

rm(dataBMI_indexdate, id1, ids, imp2 ,pop, pop.analysis, pop.gen, summary, traj, wide2, wide1)

  timepoints <- c(18,25,30,35,40)#,55,70,95)  # until 40 years
timepoints <- c(18,30,40)
test <- data.frame(wide$id, wide$.imp)
i <- 1
for (time in timepoints){
  if (! time==timepoints[length(timepoints)]){
    x1 <- timepoints[i]
    x2 <- timepoints[i+1]
    y1 <- paste0("bmi_", x1)
    y2 <- paste0("bmi_", x2)
    line <- function(t, x1, x2, y1, y2){
      y <- y1+(y2-y1)*(t-x1)/(x2-x1)
      y
    }
    for (t in seq(x1, x2, by=1)){
      name <- paste0("bmi_", t)
      test[name] <- line(t, x1, x2, wide[y1],wide[y2])
    }
    i <- i+1
  }
}

test <- test %>% rename(id=wide.id, `.imp`=wide..imp) 
# %>% 
  # left_join(wide %>% distinct(dnaix, id))
# test <- test %>% mutate(age.index=fun.age(dnaix,index.date ))
# test <- test %>% mutate(bmi.index=ifelse(paste0("bmi_",age.index))
# test <- test %>% select()


test1 <- wide %>% left_join(test %>% select(-paste0("bmi_", timepoints)),
                            by=c("id", ".imp"))

test2 <- test1 %>% pivot_longer(cols=paste0("bmi_", 18:40))
test2 <- test2 %>% mutate(value=round(value, 1))
test2 <- test2 %>% mutate(overweight=ifelse(value>=25,1, 0),
                          obese=ifelse(value>=30,1, 0),
                          ageOver=ifelse(value>=25,as.numeric(substr(name, 5,6)), 999),
                          ageObese=ifelse(value>=30,as.numeric(substr(name, 5,6)), 999),
                          deltaOver=ifelse(value>24.9,value-24.9, 0),
                          deltaObese=ifelse(value>29.9,value-29.9, 0))
test2 <- test2 %>% group_by(id, .imp) %>% mutate(yearsBMI25more=sum(overweight),
                                                 yearsBMI30more=sum(obese),
                                                 ageOnsetOver=min(ageOver),
                                                 ageOnsetObese=min(ageObese),
                                                 cumOver=sum(deltaOver),
                                                 cumObese=sum(deltaObese))
test2 <- test2 %>% mutate(everOver=ifelse(yearsBMI25more==0, 0, 1),
                          everObese=ifelse(yearsBMI30more==0, 0, 1))
test2 <- test2 %>% ungroup()
test2 <- test2 %>% mutate(ageOnsetOver=ifelse(ageOnsetOver==999, NA, ageOnsetOver),
                          ageOnsetObese=ifelse(ageOnsetObese==999, NA, ageOnsetObese))

#saveRDS(test2, "exposicionesSinTrabajar.rds")
# View(head(test2))
# summary(test2$test)
test4 <- test2 %>% distinct(.imp, id, 
                            yearsBMI25more, yearsBMI30more,
                            cumOver, cumObese, everOver, everObese,
                            ageOnsetOver, ageOnsetObese)

test4 <- test1 %>% left_join(test4)
test4 <- test4 %>% mutate(trajectorySlope=(bmi_40-bmi_18)/(40-18+1))# CAMBIA ESTO ----

test4 <- test4 %>% ungroup()
saveRDS(test4, "popExposicionesOmitted.rds")
#deja los varios bmi y guarda
saveRDS(test4, "popExposicionesFINAL.rds")

# test3 <- test2 %>% ungroup() %>% select(-c(X1.dim.wide..1., name, value, greater25)) %>% distinct()
# 
# saveRDS(test3, "popWithAgeBMIgreater25.rds")




