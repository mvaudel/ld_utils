
##
#
# This script populates a database with LD information from TopLD.
#
##


# Libraries

library(DBI)
library(glue)
library(dplyr)
library(janitor)


# Files

ld_db <- "/mnt/archive/topld/db/ld_db"
download_folder <- "/mnt/archive/topld/downloads"

ld_suffix <- "LD.csv.gz"
info_suffix <- "info_annotation.csv.gz"


# Parameters

super_populations <- c("AFR", "EAS", "EUR", "SAS")


# Set up db

db_connection <- dbConnect(RSQLite::SQLite(), ld_db)


# Gather results from TopLD

ld_tables <- list()
annotation_tables <- list()

for (population in super_populations) {
  
  folder <- file.path(download_folder, population, "SNV")
  
  for (file in list.files(folder)) {
    
    if (endsWith(file, ld_suffix)) {
      
      ld_file <- file.path(folder, file)
      
      annotation_file_name <- paste0(substr(file, 1, nchar(file) - nchar(ld_suffix)), "info_annotation.csv.gz")
      annotation_file <- file.path(folder, annotation_file_name)
      
      if (!file.exists(annotation_file)) {
        
        stop(glue("Annotation file '{annotation_file}' for LD file '{file}' not found."))
        
      }
      
      ld_table <- read.table(
        file = ld_file,
        header = T,
        sep = ",",
        stringsAsFactors = F
      ) %>% 
        clean_names() %>% 
        mutate(
          pop = population
        )
      
      ld_tables[[length(ld_tables)]] <- ld_tables
      
      annotation_table <- read.table(
        file = annotation_file,
        header = T,
        sep = ",",
        stringsAsFactors = F
      ) %>% 
        clean_names() %>% 
        select(
          uniq_id,
          position,
          rs_id,
          maf,
          ref,
          alt
        )
      
      annotation_tables[[length(annotation_tables)]] <- annotation_table
      
    }
  }
}

ld_table <- do.call("rbind", ld_tables)
annotation_table <- do.call("rbind", annotation_tables) %>% distinct()


# Write to DB

dbWriteTable(db_connection, "ld_table", ld_table, overwrite = T)
dbWriteTable(db_connection, "annotation_table", annotation_table, overwrite = T)

dbDisconnect(db_connection)

