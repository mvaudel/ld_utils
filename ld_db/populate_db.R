
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

if (file.exists(ld_db)) {
  
  removed <- file.remove(ld_db)
  
  if (!removed) {
    
    stop(glue("Current LD database {ld_db} could not be removed."))
    
  }
}

db_connection <- dbConnect(RSQLite::SQLite(), ld_db)


# Gather results from TopLD

unique_ids <- c()

for (population in super_populations) {
  
  print(glue("{Sys.time()}    Loading LD values for {population}"))
  
  folder <- file.path(download_folder, population, "SNV")
  
  for (file in list.files(folder)) {
    
    if (endsWith(file, ld_suffix)) {
      
      chromosome <- substring(strsplit(x = file, split = "_")[[1]][2], 4)
      
      ld_file <- file.path(folder, file)
      
      annotation_file_name <- paste0(substr(file, 1, nchar(file) - nchar(ld_suffix)), "info_annotation.csv.gz")
      annotation_file <- file.path(folder, annotation_file_name)
      
      if (!file.exists(annotation_file)) {
        
        stop(glue("Annotation file '{annotation_file}' for LD file '{file}' not found."))
        
      }
      
      print(glue("{Sys.time()}    Loading LD values for {population} - importing {file}"))
      
      ld_table <- read.table(
        file = ld_file,
        header = T,
        sep = ",",
        stringsAsFactors = F
      ) %>% 
        clean_names() %>% 
        mutate(
          uniq_id_1 = paste(chromosome, uniq_id_1, sep = ":"),
          uniq_id_2 = paste(chromosome, uniq_id_2, sep = ":")
        )
      
      print(glue("{Sys.time()}    Loading LD values for {population} - importing {annotation_file}"))
      
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
      
      annotation_table <- annotation_table %>% 
        mutate(
          uniq_id = paste(chromosome, uniq_id, sep = ":")
        ) %>% 
        filter(
          !uniq_id %in% unique_ids
        )
      
      unique_ids <- c(unique_ids, annotation_table$uniq_id)
      
      print(glue("{Sys.time()}    Loading LD values for {population} - Saving {file} to DB"))
      
      annotation_table_name <- paste("annotation", chromosome, sep = "_")
      dbWriteTable(db_connection, annotation_table_name, annotation_table, overwrite = F, append = T)
      
      ld_table_name <- paste("ld", population, chromosome, sep = "_")
      dbWriteTable(db_connection, ld_table_name, ld_table, overwrite = F, append = T)
      
    }
  }
}


# Done

dbDisconnect(db_connection)

print(glue("{Sys.time()}    Done!"))
