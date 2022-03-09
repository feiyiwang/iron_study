#install.packages("remotes")
#remotes::install_github("PheWAS/PheWAS")
setwd("~/Documents/Materials/iron")

library(PheWAS)
phecode_map <- PheWAS::phecode_map
write.csv(phecode_map, "phecode_map.csv", row.names = FALSE)
# upload phecode_map.csv to green_bucket
# gsutil cp ~/Documents/Materials/iron/phecode_map.csv gs://fg-production-sandbox-6_green/FeiyiWang

# https://phewascatalog.org/files/phecode_definitions1.2.csv.zip
phecode_sex <- read_csv("phecode_definitions1.2.csv", 
                                    col_types = cols(.default = col_character()))