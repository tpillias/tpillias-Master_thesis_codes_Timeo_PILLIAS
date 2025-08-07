library(dplyr)

#Reload data environment

# Find a special receiver in a specific individual dataset ####

LB7_50<- LB7 %>% filter(Receiver.x == "485850")

# After processing previous analysis (cf. XXXX) to count the number of individuals ####

unique_ids <- data_complet %>%
+   filter(Receiver.x == 485850) %>%
+   distinct(ID)

unique_ids

