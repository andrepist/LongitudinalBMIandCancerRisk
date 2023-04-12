# analysis
library(lubridate)
library(tidyverse)
library(tableone)
library(kableExtra)

fun.age<-function(date1,date2){
  #date1 birth
  output<-year(date2)-year(date1)
  ifelse(month(date2)<month(date1)|(month(date2)==month(date1)&day(date2)<day(date1)),output-1,
         output)
}

# convert to survival data

pobl <- readRDS("//epofs/apistillo/EPOFS/SIDIAP/20180514 OBECAN - Talita/Bases de datos/poblacio.rds")
expl <- readRDS("//epofs/apistillo/EPOFS/SIDIAP/20180514 OBECAN - Talita/Bases de datos/exploracions.rds")

pop <- pobl
# index.date 1 year after entering
index.date <- ymd("2009-01-01")

# adults at index date
pop <- pop %>% mutate(age.index=fun.age(dnaix, index.date))


# cancer -----
cancer_sidiap <- readRDS("E:/Andrea/Trajectories OBECAN/cancer_sidiap.rds")
cancer_cmbd <- readRDS("E:/Andrea/Trajectories OBECAN/cancer_cmbd.rds")

cancer <- cancer_sidiap  %>% select(id, dat, group, general, cod) %>% mutate(sidiap=1) %>% 
  rbind(cancer_cmbd %>% select(id, dat=date_cancer_cmbd, group=group_cmbd, general=general_cmbd, cod) %>% mutate(sidiap=2))

cancer1 <- cancer %>% arrange(dat) %>% distinct(id, .keep_all=TRUE)

pop <- pop %>% left_join(cancer1)

# followup

pop <- pop %>% select(-general, -sidiap)
pop <- pop %>% mutate(dat.fake=if_else(is.na(dat), ymd("2018/12/31"), dat),
                      end = pmin(sortida, dat.fake, ymd("2018/12/31")))
pop <- pop%>%mutate(follow_up=as.numeric((end-index.date)/365.25))

# at least 1 year
n0 <- nrow(pop %>% distinct(id))
pop <- pop %>% filter(age.index>=40)
n1 <- nrow(pop %>% distinct(id))
pop <- pop %>% filter(entrada<=ymd("2008-01-01"))
n2 <- nrow(pop1 %>% distinct(id))
pop <- pop %>% filter(dat>=index.date|is.na(dat))
n3 <- nrow(pop %>% distinct(id))
pop <- pop %>% filter(follow_up>=1)
n4 <- nrow(pop %>% distinct(id))

pop1 <- pop %>% filter(follow_up<1)
n5 <- nrow(pop1 %>% distinct(id))
library(htmlTable)
count(pop, group)%>% htmlTable()

saveRDS(pop, "analysis.pop.rds")
nrow(analysis.pop %>% filter(follow_up>=2) %>% distinct(id))

id.ana <- analysis.pop$id
id.traj <- (data %>% distinct(id))$id
length(id.ana)
length(id.traj)

length(setdiff(id.ana, id.traj)) #ok
length(setdiff(id.traj, id.ana))

# flowchart -------
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
flow1 <- grViz("
digraph a_nice_graph
      {
      
      node [fontname = Helvetica, shape = box, color = black, penwidth = 1.5, style = lisrel]
      '@@1';'@@2';'@@3';'@@4';'@@5';'@@6';
       node [fontname = Helvetica, shape = box, color = '#FF8700', penwidth = 1.5, style = lisrel]
      '@@7';
      
      blank1[label = '', width = 0.01, height = 0.01]
      blank2[label = '', width = 0.01, height = 0.01]
      blank3[label = '', width = 0.01, height = 0.01]
      blank4[label = '', width = 0.01, height = 0.01]
      blank5[label = '', width = 0.01, height = 0.01]
     
    
      '@@1' -> blank1[ dir = none ];
      blank1 -> '@@2'[ minlen = 0 ,arrowhead = none, tailport = T];
      blank1 -> '@@3' ;  
      '@@3' -> blank2[dir = none];
      blank2 -> '@@4'[ minlen = 0 ,arrowhead = none, tailport = T];
      blank2 -> '@@5';
      '@@5' -> blank3[dir = none];
      blank3 -> '@@6'[ minlen = 0 ,arrowhead = none, tailport = T];
      blank3 -> '@@7'
     
      
      }
      [1]: paste0('SIDIAP population aged &#8805; 40 years on 01/01/2009','\\n','N = ', format(n1,big.mark = ',', decimal.mark='.'))
      [2]:  paste0('Individuals without 1 year of history in the SIDIAP', '\\n','N = ', format(n1-n2,big.mark = ',', decimal.mark='.'))
      [3]:  paste0('Eligible population', '\\n','N = ', format(n2,big.mark = ',', decimal.mark='.')) 
      [4]:  paste0('Individual with history of cancer', '\\n', 'N = ', format(n2-n3, big.mark = ',',decimal.mark='.'))
      [5]:  paste0('Eligible population', '\\n','N = ', format(n3,big.mark = ',', decimal.mark='.')) 
       [6]:  paste0('Individuals without 1 year of follow-up','\\n', 'N = ', format(n3-n4,big.mark = ',', decimal.mark='.')) 
       [7]:  paste0('Population for the time-to-event analysis', '\\n','N = ', format(n4,big.mark = ',', decimal.mark='.')) 
        ", width = 400,height = 1200)
flow1
flow1%>% export_svg %>% charToRaw %>% rsvg %>% png::writePNG("flowchart.png")
