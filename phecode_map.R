##### In local machine #####
# install.packages("remotes")
# remotes::install_github("PheWAS/PheWAS")
setwd("~/Documents/Materials/iron")

library(PheWAS)
phecode_map <- PheWAS::phecode_map
write.csv(phecode_map, "phecode_map.csv", row.names = FALSE)
# upload phecode_map.csv to green_bucket
# gsutil cp ~/Documents/Materials/iron/phecode_map.csv gs://fg-production-sandbox-6_green/FeiyiWang

# https://phewascatalog.org/files/phecode_definitions1.2.csv.zip
phecode_dict <- read_csv("phecode_definitions1.2.csv", 
                                    col_types = cols(.default = col_character()))


##### In Finngen sandbox #####
# check if all the required packages have been installed
#'dplyr', 'tidyr', 'ggplot2', 'meta', 'ggrepel', 'DT', 'logistf', 'lmtest'
# need to install meta
# check if meta is in the list
ListOfPackages = library()$results[,1]
'meta' %in% ListOfPackages
# meta is not in the list
# upload meta tar file and phewas tar file

install.packages("/finngen/green/FeiyiWang/mathjaxr_1.6-0.tar.gz", '/home/ivm/R/x86_64-pc-linux-gnu-library/4.1', repos = NULL, type="source")
install.packages("/finngen/green/FeiyiWang/pbapply_1.5-0.tar.gz", '/home/ivm/R/x86_64-pc-linux-gnu-library/4.1', repos = NULL, type="source")
install.packages("/finngen/green/FeiyiWang/metafor_3.0-2.tar.gz", '/home/ivm/R/x86_64-pc-linux-gnu-library/4.1', repos = NULL, type="source")
install.packages("/finngen/green/FeiyiWang/CompQuadForm_1.4.3.tar.gz", '/home/ivm/R/x86_64-pc-linux-gnu-library/4.1', repos = NULL, type="source")
install.packages("/finngen/green/FeiyiWang/meta_5.2-0.tar.gz", '/home/ivm/R/x86_64-pc-linux-gnu-library/4.1', repos = NULL, type="source")
install.packages("/finngen/green/FeiyiWang/PheWAS_0.99.5-5.tar.gz", '/home/ivm/R/x86_64-pc-linux-gnu-library/4.1', repos = NULL, type = "source")
library(PheWAS)

sex_df <- read.csv('sex_df.csv')
events <- read.csv('events.csv')
events <- events %>% group_by(id,vocabulary_id,code) %>% mutate(count=cumsum(count1))
events$count1 <- NULL
events <- as.data.frame(events)
events_pheno <- createPhenotypes(events, id.sex = sex_df, min.code.count = 1, full.population.ids = sex_df$id)
write.csv(events_pheno, "events_pheno.csv", row.names = FALSE)



