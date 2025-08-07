library(ggplot2)
library(dplyr)
library(lubridate)
# Import Rdata environment
# Month filtering command ####

# PI spawning period
filtrer_october_2022<- function(df){
  if (!"Date.and.Time..UTC." %in% colnames(df)) {
  return(df)  # Si pas de colonne date, on laisse tel quel
}


df %>%
  mutate(Date.and.Time..UTC. = as.POSIXct(Date.and.Time..UTC., format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>% 
  filter(format(Date.and.Time..UTC., "%Y") %in% c( "2022")) %>%
  filter(format(Date.and.Time..UTC., "%m") %in% c("10"))
}

# LB spawning period 
filtrer_mois_automne <- function(df) {
  # Vérifier si la colonne "Date.and.Time..UTC." existe
  if (!"Date.and.Time..UTC." %in% colnames(df)) {
    return(df)  # Si pas de colonne date, on laisse tel quel
  }
  
  df %>%
    mutate(Date.and.Time..UTC. = as.POSIXct(Date.and.Time..UTC., format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>% 
    filter(format(Date.and.Time..UTC., "%Y") %in% c("2022")) %>%
    filter(
      (format(Date.and.Time..UTC., "%m") == "07" & as.integer(format(Date.and.Time..UTC., "%d")) >= 17) |
        (format(Date.and.Time..UTC., "%m") == "08" & as.integer(format(Date.and.Time..UTC., "%d")) <= 15)
    )}

  # Seasonal filtering ####

  filtrer_mois_saison <- function(df) {
  # Vérifier si la colonne "Date.and.Time..UTC." existe
  if (!"Date.and.Time..UTC." %in% colnames(df)) {
    return(df)  # Si pas de colonne date, on laisse tel quel
  }
  
  
  df %>%
    mutate(Date.and.Time..UTC. = as.POSIXct(Date.and.Time..UTC., format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>% 
    filter(format(Date.and.Time..UTC., "%m") %in% c("06", "07", "08"))
}

# Morph selection ####

objets <- ls(pattern = "^LB")
objets <- ls(pattern = "^PI")

# Dataset filtering #### 

for (obj in objets) {
  if (is.data.frame(get(obj))) {
    assign(obj, filtrer_octobre_2022(get(obj)))  # Met à jour l'objet filtré
  }
}

# Dataset binding and exporting ####

data_complet <- bind_rows(lapply(objets, function(obj) {
  df <- get(obj)  
  df$Individu <- obj  
  return(df)
}))

# Vérification
head(data_complet)
tail(data_complet)
dim(data_complet)

write.table(data_complet, file = "PI_spawning_month.csv", 
            sep = ",", dec = ".", 
            row.names = FALSE, col.names = TRUE, 
            fileEncoding = "UTF-8")
