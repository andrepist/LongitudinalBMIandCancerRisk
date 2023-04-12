##definition of DM


# Diagnostic codes

DM <- probl%>%filter(str_detect(cod,"E11"))%>%
  mutate(group="T2DM", general="Diabetes 2")

#For people with several DM diagnoses, we keep the first one (although this also keeps 2 dx on the same date)
DM <- DM %>% group_by(id) %>% filter(dat == min(dat)) %>%
  ungroup()
n1_DM <- nrow(DM %>% distinct(id))

#we only keep the first dx per person
DM <- DM %>%
  group_by(id) %>%
  mutate(n_diagnoses = 1:n())

DM <- DM %>% group_by(id) %>% filter(n_diagnoses<=1) %>%
  ungroup()

n2_DM <- nrow(DM %>% distinct(id))

#n1_DM and n2_DM now have the same N as the number of obs in the all_cancer dataframe

saveRDS(DM,"DM.rds")


## HTA definition

# Diagnostic codes

HTA <- probl%>%filter(str_detect(cod,"I10"))%>%
  mutate(group="HTA", general="Hypertension")

#For people with several HTA diagnoses, we keep the first one (although this also keeps 2 dx on the same date)
HTA <- HTA %>% group_by(id) %>% filter(dat == min(dat)) %>%
  ungroup()
n1_HTA <- nrow(HTA %>% distinct(id))

#we only keep the first dx per person
HTA <- HTA %>%
  group_by(id) %>%
  mutate(n_diagnoses = 1:n())

HTA <- HTA %>% group_by(id) %>% filter(n_diagnoses<=1) %>%
  ungroup()

n2_HTA <- nrow(HTA %>% distinct(id))

#n1_HTA and n2_HTA now have the same N as the number of obs in the HTA dataframe

saveRDS(HTA,"HTA.rds")
