library(tidyverse)
library(lubridate)
library(mice)

fun.age<-function(date1,date2){
  #date1 birth
  output<-year(date2)-year(date1)
  ifelse(month(date2)<month(date1)|(month(date2)==month(date1)&day(date2)<day(date1)),output-1,
         output)
}

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
timepoints <- c(40,   55,70,118)
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

# add age at index
pop <- readRDS("E:/Andrea/Trajectories OBECAN/populationCox.rds")
pop <- pop %>% ungroup() %>% select(id, age.index) %>% distinct()

test1<-test1 %>% left_join(pop)
rm(pop, wide1, wide2)

#change manually .imp
test2<-test1 %>% filter(.imp==5) %>% select(-c(bmi_18, bmi_30, paste0("bmi_", 110:118), n, qmedea, smoking, alcohol))%>%
  pivot_longer(cols=paste0("bmi_", 40:109))

test2<-test2 %>% mutate(match=as.numeric(substr(name, 5,7)))

test3<-rbind(test3,
             test2 %>% filter(age.index==match) %>% rename(bmi_atindex=value))

saveRDS(test3, "bmiindex.rds")
