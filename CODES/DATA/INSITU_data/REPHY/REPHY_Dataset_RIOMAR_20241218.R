#### Script REPHY : curating data tables for hydrological measurements
# V. POCHIC 2024/12/18

### Required packages ####

library(tidyverse)

work_dir <- "/home/terrats/Desktop/RIOMAR/DATA/REPHY/"

#### Import data ####

# Importing data tables with the **correct file encoding** allows to preserve accents
DataREPHY_MA <- read.csv2( work_dir %>% file.path('SEANOE_extraction', 'REPHY_Manche_Atlantique_1987-2022.csv'), fileEncoding = "ISO-8859-1")
DataREPHY_Med <- read.csv2( work_dir %>% file.path('SEANOE_extraction', 'REPHY_Med_1987-2022.csv'), fileEncoding = "ISO-8859-1")

# Creating a vector to import all columns as characters
classes_char <- rep('character', 56)
# CAREFUL! ALL 2023 DATA LACK TEMPERATURE! Problem seen with Maud Lemoine 
# on 2024/01/09. To be continued. These data are not yet published.
DataREPHY_2023 <- read.csv2( work_dir %>% file.path('Extraction SEANOE_REPHY_phyto-Hydro_Manche_Atl-Med validé 19122023.csv'), 
                            fileEncoding = "ISO-8859-1", colClasses = classes_char)

# Merge the first 2 datasets
DataREPHY <- bind_rows(DataREPHY_MA, DataREPHY_Med) %>%
  # This line allows to remove the problematic row that separates the hydro and phyto datasets
  filter(Passage...Mois != 'Passage : Mois')

# Binding the 2 datasets
DataREPHY_8723 <- bind_rows(DataREPHY, DataREPHY_2023)


### Load tables used to enrich the dataset
# Load table "Zones_marines"
ZM <- read.csv( work_dir %>% file.path('Zones_marines.csv') , sep = ';', header = TRUE)

#### Formatting the dataframe with better column names ####

# Extracting the numeric code for the ZM
Table1 <- DataREPHY_8723 %>%
  mutate(ZM_Quadrige_Numero = as.numeric(str_extract(Lieu.de.surveillance...Entité.de.classement...Libellé, '[:alnum:]+')))

# Change column names (to match ZM)
colnames(Table1)[which(names(Table1) == "Lieu.de.surveillance...Mnémonique")] <- "Code_point_Mnemonique"
colnames(Table1)[which(names(Table1) == "Lieu.de.surveillance...Libellé")] <- "Code_point_Libelle"
colnames(Table1)[which(names(Table1) == "Passage...Date")] <- "Date"
colnames(Table1)[which(names(Table1) == "Coordonnées.passage...Coordonnées.minx")] <- "lon"
colnames(Table1)[which(names(Table1) == "Coordonnées.passage...Coordonnées.miny")] <- "lat"
colnames(Table1)[which(names(Table1) == "Résultat...Libellé.unité.de.mesure.associé.au.quintuplet")] <- "Mesure_Unite"
colnames(Table1)[which(names(Table1) == "Résultat...Symbole.unité.de.mesure.associé.au.quintuplet")] <- "Mesure_Symbole"
colnames(Table1)[which(names(Table1) == "Résultat...Nom.du.taxon.référent")] <- "Taxon"
colnames(Table1)[which(names(Table1) == "Résultat...Libellé.du.groupe.de.taxon")] <- "Groupe_Taxon"
colnames(Table1)[which(names(Table1) == "Résultat...Valeur.de.la.mesure")] <- "Valeur_mesure"
colnames(Table1)[which(names(Table1) == "Prélèvement...Immersion")] <- "Profondeur.metre"
colnames(Table1)[which(names(Table1) == "Prélèvement...Niveau")] <- "Prelevement.niveau"
colnames(Table1)[which(names(Table1) == "Résultat...Code.paramètre")] <- "Code.parametre"
colnames(Table1)[which(names(Table1) == "Résultat...Libellé.paramètre")] <- "Parametre"
colnames(Table1)[which(names(Table1) == "Résultat...Libellé.méthode")] <- "Méthode.analyse"
colnames(Table1)[which(names(Table1) == "Résultat...Niveau.de.qualité")] <- "Qualite.resultat"
colnames(Table1)[which(names(Table1) == "Prélèvement...Niveau.de.qualité")] <- "Qualite.prelevement"
colnames(Table1)[which(names(Table1) == "Passage...Service.saisisseur...Libellé")] <- "Service.saisie"
colnames(Table1)[which(names(Table1) == "Passage...Heure")] <- "Heure"
colnames(Table1)[which(names(Table1) == "Résultat...Service.analyste...Code")] <- "Service.analyse"
colnames(Table1)[which(names(Table1) == "Prélèvement...Service.préleveur...Code")] <- "Service.prelevement"
colnames(Table1)[which(names(Table1) == "Prélèvement...Identifiant.interne")] <- "ID.interne.prelevement"
colnames(Table1)[which(names(Table1) == "Passage...Identifiant.interne")] <- "ID.interne.passage"

#### Curate table to keep only desired variables ####
Table1 <- Table1 %>%
  dplyr::select(c('ZM_Quadrige_Numero', 'Code_point_Mnemonique', 'Code_point_Libelle', 'Date', 
                  'Heure', 'lon', 'lat', 'Mesure_Unite', 'Mesure_Symbole', 'Taxon', 'Valeur_mesure', 
                  'Prelevement.niveau', 'Profondeur.metre', 'Code.parametre', 'Parametre', 'Méthode.analyse',
                  'Qualite.prelevement', 'Qualite.resultat', 'ID.interne.prelevement', 'ID.interne.passage'))

# Modifying date format so that it gives the year and month, and getting rid of rows with no Year value
Table1 <- Table1 %>%
  mutate(Date = dmy(Date)) %>%
  # modifies the date format
  mutate(Day = day(Date)) %>%
  mutate(Month = month(Date, label = F)) %>%
  mutate(Year = year(Date)) %>%
  filter(!is.na(Year))

# Transform the measured values into numbers
# Caution! The decimal separator is a comma, 
# we need first to transform it to a full stop to avoid fuckery
Table1$Valeur_mesure <-  str_replace_all(Table1$Valeur_mesure, ',', '.')
Table1 <- Table1 %>%
  mutate(Valeur_mesure = as.numeric(Valeur_mesure))

# Build a dataframe with the list of parameters for which we have data
Table_parametres <- as_tibble(unique(Table1$Parametre)) %>%
  add_column(unique(Table1$Code.parametre))

# Save said table
write.csv2(Table_parametres, 'Parametres_REPHY.csv', row.names = FALSE, fileEncoding = 'ISO-8859-1')

#### Tidying table structure ####

## Associate a region with each ZM code
Table1 <- left_join(Table1, ZM, by='ZM_Quadrige_Numero', suffix=c('',''))

# Basically, we want the hydrological measurements as columns and the phytoplankton taxa as rows

# Separate the table into 2 : 1 for hydrology and the other for phytoplankton
Table1_hydro <- Table1 %>%
  filter(Taxon == "")

Table1_phyto <- Table1 %>%
  filter(Taxon != "")

### We will focus on the hydro table ###
Table1_hydro_select <- Table1_hydro %>%
  select(c('ID.interne.passage', 'Qualite.resultat', 'Code.parametre', 'Valeur_mesure', 
           'Mesure_Unite', 'Mesure_Symbole', 'Prelevement.niveau','Profondeur.metre', 
           'Méthode.analyse', 'Code.Region', 'Region', 'Date', 'Day',
           'Month', 'Year', 'Heure', 'Code_point_Libelle', 'Code_point_Mnemonique',
           'lon', 'lat')) %>%
  #filter out Code.Region = 0 (only a few mysterious events in 2011-2013)
  filter(Code.Region != 0)
  # Filter out data with suspect quality
  #filter(Qualite.resultat == 'Bon') %>%
  # Filter out data taken not at the surface
  #filter(Prelevement.niveau == 'Surface (0-1m)' | Prelevement.niveau == '2 metres' |
  #          Prelevement.niveau == 'de 3 à 5 metres')

# write the table to csv file
# write.csv2(Table1_hydro_select, work_dir %>% file.path('Table1_REPHY_hydro_RIOMAR.csv'), row.names = FALSE, fileEncoding = "ISO-8859-1")
write.csv2(Table1_hydro_select, gzfile( work_dir %>% file.path("Table1_REPHY_hydro_RIOMAR.csv.gz") ), row.names = FALSE, fileEncoding = "ISO-8859-1")

#### Selecting regions or sites ####

# Let's restrict the table to only sites in regions 21, 22 and 23 (Atlantic)
Table1_hydro_Atlantic <- Table1_hydro_select %>% filter(Code.Region %in% c(21, 22, 23))

# Writing the table as csv file
write.csv2(Table1_hydro_Atlantic, gzfile( work_dir %>% file.path('Table_REPHY_hydro_Atlantic_RIOMAR.csv.gz') ), row.names = FALSE, fileEncoding = 'ISO-8859-1')

# In this region, we create a table with all the sites and their coordinates
Table_sites <- Table1_hydro_Atlantic %>%
  group_by(Code_point_Libelle) %>%
  # exclude rows where coordinates are not entered
  filter(is.na(lon) == FALSE & is.na(lat) == FALSE) %>%
  # convert lon and lat variables to numeric
  mutate(lon = as.numeric(lon)) %>%
  mutate(lat = as.numeric(lat)) %>%
  summarise(# Longitude
    Longitude.mean = mean(lon),
    # Latitude
    Latitude.mean = mean(lat),
    n_results = n(),
    .groups = 'keep')

# This table can be used to select the appropriate sites.
write.csv2(Table_sites, gzfile( work_dir %>% file.path('Table_REPHY_hydro_sites_Atlantic_RIOMAR.csv.gz') ), row.names = FALSE, fileEncoding = 'ISO-8859-1')
