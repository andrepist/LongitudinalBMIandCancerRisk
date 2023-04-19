library(mice)
library(tidyverse)
library(tictoc)
library(lubridate)

imp <- complete(imp_single3 , "long", include=FALSE)
saveRDS(imp, "imp/imp3.rds")
rm(imp1,imp_single1)
gc()

imp2 <- complete(imp2, "long", include=FALSE)
imp2<-imp2  %>% select(-c(date.cancer, group, dat.fake, end, ruralitat, date.cvd, date.hta, date.dm, year.bmi, 
                                          bariatric, date.bariatric, age.bariatric, age.cvd, age.hta, age.dm, age.cancer))
saveRDS(imp2, "imp/imp2.rds")
rm(imp2)

imp3 <- complete(imp3,"long", include=FALSE)
saveRDS(imp3, "imp/imp3.rds")
rm(imp3)

imp4 <- complete(imp, 4, include=FALSE)
saveRDS(imp4, "imp4.rds")
rm(imp4)
imp5 <- complete(imp, 5, include=FALSE)
saveRDS(imp5, "imp5.rds")
rm(imp5)