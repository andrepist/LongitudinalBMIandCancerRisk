# selecciona muestra m?s peque?a del dataset, random-----------
library(splines)
library(mice)
library(tidyverse)
library(tictoc)
library(lubridate)

fun.age<-function(date1,date2){
  #date1 birth
  output<-year(date2)-year(date1)
  ifelse(month(date2)<month(date1)|(month(date2)==month(date1)&day(date2)<day(date1)),output-1,
         output)
}

data.bmi1 <- readRDS("dataBMI_indexdate_new.rds")
data.bmi1<-dataBMI_indexdate_new
# variables como factores!!!
data.bmi <- data.bmi1 %>% rename(date.cancer=dat) %>% 
  mutate(qmedea=ifelse(qmedea=="missing", NA_character_, qmedea))
data.bmi <- data.bmi %>% mutate(qmedea=as.factor(qmedea),
                                smoking=as.factor(smoking))
data.bmi <- data.bmi %>% 
  mutate(date.bmi=if_else(is.na(bmi), NA_Date_, date.bmi))

 # set.seed(5)
 # sample.id <- sample(data.bmi$id, 5000)
 # data <- data.bmi %>% filter(id%in%sample.id)

data <- data.bmi


# sobre esta, imputamos BMI al 01/01/2010, 2013, 2016
# 18 25 30 35 40  55 70 95
minage <- min(data$age.bmi, na.rm = T)
age2 <- 25
age3 <- 30
age4 <- 35
age5 <- 40
age6 <- 55
age7 <- 70
age8 <- 95
maxage <- max(data$age.bmi, na.rm = T)
timepoints <- c(minage,age2,age3,age4,age5,age6,age7,age8, maxage)

timepoints <- c(minage,age3,age5,age6,age7, maxage)

#pongo variable year que es el a?o de medida del bmi

data <- data %>% mutate(year.bmi=year(date.bmi))

# si mas medidas en un a?o, la mas cercana a la media

data <- data %>%
  group_by(id,year.bmi) %>% 
  mutate(media=mean(bmi))%>%
  ungroup() %>% 
  arrange(abs(media-bmi), year.bmi) %>% 
  distinct(id, year.bmi,.keep_all = TRUE) %>% 
  select(-media)

#a?ado filas "falsas" a las edades que me interesan, cogiendo informacion previa y borrando bmi
data <- data %>% mutate(fake="no")
# hago subset de info mas cercana a esto
for (j in timepoints){
  fila <- data %>% arrange(abs(age.bmi-j)) %>% 
    distinct(id, .keep_all = TRUE) %>% 
    mutate(bmi=NA, age.bmi=j)%>% 
    mutate(fake="yes")
  data <- rbind(data, fila)
}

data <- data %>% 
  filter(!is.na(age.bmi)) 

nudos <- timepoints
k <- length(nudos)
### calculate B-spline
X <- bs(data$age.bmi,
        knots = nudos,
        Boundary.knots = c(nudos[1], nudos[k]+1),
        degree = 1)
X <- X[,-(k + 1)]
dimnames(X)[[2]] <- paste("x", 1:ncol(X), sep = "")
#a?ado como columna a los datos
data <- cbind(data, X)

#data <- data %>% rename(date.cancer=dat)

# add indicators of CVD, HTN, DM 
data <- data %>% mutate(age.cvd=fun.age(dnaix, date.cvd),
                        age.hta=fun.age(dnaix, date.hta),
                        age.dm=fun.age(dnaix, date.dm))

data <- data %>% mutate(i.cvd=ifelse(is.na(age.cvd)|(!is.na(age.bmi)&age.bmi<age.cvd)|age.bmi<age.cvd, 0, 1),
                        i.htn=ifelse(is.na(age.hta)|(!is.na(age.bmi)&age.bmi<age.hta)|age.bmi<age.hta, 0, 1),
                        i.dm=ifelse(is.na(age.dm)|(!is.na(age.bmi)&age.bmi<age.dm)|age.bmi<age.dm, 0, 1))

data <- data %>% mutate(index_ch=ifelse(is.na(index_ch), 0 , index_ch))
# age cancer
data <- data %>% mutate(age.cancer=fun.age(dnaix, date.cancer))
data <- data %>% mutate(i.cancer1yprior=ifelse(is.na(age.cancer)|(!is.na(age.cancer)&age.bmi<(age.cancer-1))|age.bmi<(age.cancer-1), 0, 1))
data <- data %>% mutate(cancer=ifelse(!is.na(age.cancer), 1, 0))
data <- data %>% mutate(cancer.group=ifelse(is.na(group), "No cancer", group))

# indicator bariatric surgery
data <- data %>% mutate(age.bariatric=fun.age(dnaix, date.bariatric))
data <- data %>% mutate(i.bariatric=ifelse(is.na(age.bariatric)|(!is.na(age.bariatric)&age.bmi<(age.bariatric))|age.bmi<(age.bariatric), 0, 1))


# light the dataset, remove unnecessary columns

data <- data %>% select(-c(date.cancer, group, dat.fake, end, ruralitat, date.cvd, date.hta, date.dm, year.bmi, 
                            bariatric, date.bariatric, age.bariatric, age.cvd, age.hta, age.dm, age.cancer))
saveRDS(data, "data.beforeImputation.rds")
saveRDS(data, "data.imp_new.rds")

Y <- c("smoking","alcohol","qmedea","bmi")
metodo <- make.method(data)
metodo[1:length(metodo)] <- ""
metodo[Y] <- c("pmm","pmm", "pmm","2l.pan")

#set up predictor matrix,  It specifies the target variable or block in the rows, and the predictor variables on the columns. 
# An entry of 0 means that the column variable is NOT used to impute the row variable or block. A nonzero value indicates that it is used.
pred <- make.predictorMatrix(data)
pred[1:nrow(pred), 1:ncol(pred)] <- 0
pred[Y, "id"] <- (-2)
pred["bmi", paste("x", 1:6, sep = "")] <- 2 # change
pred[Y, "sexe"] <- 1
pred[Y, "nationality"] <- 1
pred[Y, "follow_up"] <- 1
pred[Y, "index_ch"] <- 1
pred[Y, "cancer.group"] <- 1
pred[Y[1:3], "age.index"] <- 1
pred["smoking","qmedea" ] <- 1
pred["smoking","alcohol" ] <- 1
pred["qmedea","smoking" ] <- 1
pred["qmedea","alcohol" ] <- 1
pred["alcohol","smoking" ] <- 1
pred["alcohol","qmedea" ] <- 1
pred["bmi","smoking" ] <- 1
pred["bmi","alcohol" ] <- 1
pred["bmi","qmedea" ] <- 1
pred["bmi", "i.cvd"] <- 2
pred["bmi", "i.htn"] <- 2
pred["bmi", "i.dm"] <- 2
pred["bmi", "i.cancer1yprior"] <- 2
pred["bmi", "i.bariatric"] <- 2

saveRDS(pred, "pred_new2.rds")
saveRDS(metodo, "metodo2.rds")
#free space
#rm(X, fila, data.bmi, data1, data.bmi1)
rm(imp.set, data.plot, imp)
# library(micemd)
# tic()
# imp.parall <- mice.par(data, m = 5, method = metodo, predictorMatrix = pred, maxit = 10,
#                        seed = 1, nnodes = 5)
# toc()
# 

# test <- list()
# for (mi in 6:7){
# data <- readRDS("data.imp_new.rds")
# pred<-readRDS("pred_new.rds")
# metodo<-readRDS("metodo.rds")
# library(splines)
# library(mice)
# library(tidyverse)
# library(tictoc)
# library(lubridate)
gc()
tic()
#data <- data %>% as.data.frame()
mi=4
# imp <- mice.par(data, m = 1, method = metodo, predictorMatrix = pred, maxit = 10,
#                 seed = mi, nnodes = 3, print = TRUE)
imp <- mice(data, method = metodo, pred = pred, m = 1,
            maxit = 10, seed = mi, print = TRUE)
toc()
saveRDS(imp, paste0("imp/imp_single", mi, ".rds"))
# rm(imp)
# }

t <- rbind( complete(test[[1]], "long", include=FALSE),
            complete(test[[2]], "long", include=FALSE), 
            complete(test[[3]], "long", include=FALSE),
            complete(test[[4]], "long", include=FALSE),
            complete(test[[5]], "long", include=FALSE))


plot(test[[3]], layout=c(2, 4))
plot(test[[5]], layout=c(2, 4))

saveRDS(imp, "imp_total_new.rds")


# diagnostic --------
plot(imp, layout=c(2, 4))
plot(imp.parall, layout=c(2, 4))
densityplot(imp)

imp<-readRDS("//epofs/RWCancerEpi/Multilevel raster imputation/Alicia/imp_single5.rds")

imp.set <- complete(imp, "long", include=FALSE)
 imp2 <- complete(imp, 2, include=FALSE)
saveRDS(imp2, "imp2.rds")
rm(imp2)
imp3 <- complete(imp, 3, include=FALSE)
saveRDS(imp3, "imp3.rds")
rm(imp3)
imp4 <- complete(imp, 4, include=FALSE)
saveRDS(imp.set, "imp/imp4.rds")
rm(imp4)
imp5 <- complete(imp, 5, include=FALSE)
saveRDS(imp5, "imp5.rds")
rm(imp5)

imp.set <- imp.set%>% mutate(age.cancer=fun.age(dnaix, date.cancer))

imp.set <- imp.set %>% left_join(data.imp_new)

trajectories <- rbind(imp1 %>% filter(fake=="yes")%>% select(id, bmi, age.bmi, qmedea, smoking, alcohol) %>% mutate(`.imp`=1),
                      imp2 %>% filter(fake=="yes")%>% select(id, bmi, age.bmi, qmedea, smoking, alcohol) %>% mutate(`.imp`=2),
                      imp3 %>% filter(fake=="yes")%>% select(id, bmi, age.bmi, qmedea, smoking, alcohol) %>% mutate(`.imp`=3),
                      imp4 %>% filter(fake=="yes")%>% select(id, bmi, age.bmi, qmedea, smoking, alcohol) %>% mutate(`.imp`=4),
                      imp5 %>% filter(fake=="yes")%>% select(id, bmi, age.bmi, qmedea, smoking, alcohol) %>% mutate(`.imp`=5))
saveRDS(trajectories, "trajectories_nosmokingmedea.rds")

trajectories <- rbind(imp1 %>% select(fake, entrada, sortida, dnaix, id, bmi, age.bmi, qmedea, smoking, alcohol) %>% mutate(`.imp`=1),
                      imp2 %>% select(fake,entrada, sortida, dnaix,id, bmi, age.bmi, qmedea, smoking, alcohol) %>% mutate(`.imp`=2),
                      imp3 %>% select(fake,entrada, sortida, dnaix,id, bmi, age.bmi, qmedea, smoking, alcohol) %>% mutate(`.imp`=3),
                      imp4 %>% select(fake,entrada, sortida, dnaix,id, bmi, age.bmi, qmedea, smoking, alcohol) %>% mutate(`.imp`=4),
                      imp5 %>% select(fake,entrada, sortida, dnaix,id, bmi, age.bmi, qmedea, smoking, alcohol) %>% mutate(`.imp`=5))

set.seed(5)
sample.id2 <- sample(imp.set$id,20)


data.plot <- rbind(imp1 %>% filter(id%in%sample.id2)%>% mutate(`.imp`=1),
                      imp2 %>% filter(id%in%sample.id2) %>% mutate(`.imp`=2),
                      imp3 %>% filter(id%in%sample.id2) %>% mutate(`.imp`=3),
                      imp4 %>% filter(id%in%sample.id2) %>% mutate(`.imp`=4),
                      imp5 %>% filter(id%in%sample.id2) %>% mutate(`.imp`=5))

data.plot <- imp.set %>% filter(id%in%sample.id2)

ggplot()+
  geom_line(data=data.plot %>% filter(fake=="yes"),aes(age.bmi,bmi, group=.imp),colour="dodgerblue1")+
  geom_point(data=data.plot %>% filter(fake=="no"),aes(age.bmi,bmi))+
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        panel.border = element_rect(color = "grey", fill = NA))+
  ylim(c(10,60))+
  #geom_vline(xintercept = nudos, linetype="dashed", color = "grey")+
  facet_wrap(.~id, drop = FALSE )+
  geom_vline(aes(xintercept=age.cvd, colour="CVD"), data=data.plot ) +
  geom_vline(aes(xintercept=age.hta, colour="HTA"), data=data.plot ) +
  geom_vline(aes(xintercept=age.dm, colour="DM"), data=data.plot ) +
  geom_vline(aes(xintercept=age.cancer, colour="Cancer"), data=data.plot) +
  geom_vline(aes(xintercept=age.bariatric, colour="Bariatric surgery"), data=data.plot) +
 # scale_x_continuous(breaks = nudos, labels=nudos)+
  ylab("BMI")+
  xlab("Age of measurement")+
  labs(colour="Condition")


# extra -------

# check N of BMI measurements per age group
# data <- data %>% 
#   mutate(agebmi.gr=
#            ifelse(age.bmi<=age2, paste(minage, "-",age2), 
#                   ifelse(age.bmi>age2&age.bmi<=age3, paste0(age2+1, "-", age3), 
#                          ifelse(age.bmi>age3&age.bmi<=age4, paste0(age3+1, "-", age4),
#                                 ifelse(age.bmi>age4&age.bmi<=age5, paste0(age4+1, "-", age5),
#                                        ifelse(age.bmi>age5&age.bmi<=age6, paste0(age5+1, "-", age6),
#                                               ifelse(age.bmi>age6&age.bmi<=age7, paste0(age6+1, "-", age7),
#                                                      ifelse(age.bmi>age7&age.bmi<=95, paste0(age7+1, "-",95),
#                                                             ifelse(age.bmi>95&age.bmi<=maxage, paste0(95+1, "-",maxage),
#                                                                    "Missing")))))))))
# 
# a <- prop.table(table(data$agebmi.gr))*100
# a <- tibble(a, row.names(a))
# a %>% htmlTable()

# plot of recording age ------
data.bmi <- data.bmi %>% mutate(age.bmi=fun.age(dnaix, date.bmi)) 
test <- data.bmi %>% count(age.bmi, sexe) %>% na.omit() %>% 
  pivot_wider(names_from=sexe, values_from=n)%>%
  mutate(N_D=sum(D, na.rm = TRUE), N_H = sum(H, na.rm = TRUE)) %>% 
  mutate(percH=H/N_H*100, percD=D/N_D*100)

ggplot(test)+
  geom_line(aes(age.bmi, percH))+
  geom_line(aes(age.bmi, percD))

test2 <- data %>% count(age.bmi, sexe) %>% na.omit() %>% 
  pivot_wider(names_from=sexe, values_from=n)%>%
  mutate(N_D= nrow(data %>% filter(sexe=="D") %>% distinct(id)),
         N_H = nrow(data %>% filter(sexe=="H") %>% distinct(id))) %>% 
  mutate(percH=H/N_H*100, percD=D/N_D*100)

ggplot(test2)+
  geom_line(aes(age.bmi, percH))+
  geom_line(aes(age.bmi, percD))

test3 <- data %>% count(age.bmi, sexe) %>% na.omit() %>% 
  pivot_wider(names_from=sexe, values_from=n)%>%
  mutate(N_D=sum(D, na.rm = TRUE), N_H = sum(H, na.rm = TRUE)) %>% 
  mutate(percH=H/N_H*100, percD=D/N_D*100)

ggplot(test3)+
  geom_line(aes(age.bmi, percH))+
  geom_line(aes(age.bmi, percD))
