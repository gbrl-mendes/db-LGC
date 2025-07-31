---
  title: "BLASTr  NT + 12sLGCdb"
author: "Hilario, HO; Mendes, GA"
date: "09/2024"
---
  
  ## Carregando bibliotecas ----
{
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggbreak)
  library(phyloseq)
  library(Biostrings)
  library(Matrix)
  library(ShortRead)
  library(dada2)
  library(DECIPHER)
  library(future)
  library(ggh4x)
  library(vegan)
  library(plotly)
  library(ggtext)
  library(BLASTr)
  library(writexl)
  library(readxl)
}

## Caminhos ----
{
  prjct_path <- "~/projetos/eDNA_Doce/"
  results_path <- paste0(prjct_path,"results/")
  figs_path <- paste0(results_path,"figures/")
  tbl_path <- paste0(prjct_path,"tables/")
  env_path <- paste0(prjct_path,"environments/")
  blast_path <- paste0(results_path,"blastr/")
}

## Obtencao dos dados ----
{
  # final_results <- read_excel(paste0(tbl_path,"/", "eDNA_Doce-Complete_analysis_results-2024-02-19.xlsx")) %>% tibble()
  final_results <- read.csv(paste0(tbl_path,"/","eDNA_Doce-Complete_analysis_results-2024-02-19_Resultados completos.csv"), sep = ",", check.names = FALSE) 
}

## Alteracoes pos-pipe ----

less_final_results <- final_results %>% 
  select(Unique_File_name,
         "ASV header",
         "ASV (Sequence)",
         "Primer expected length",
         "ASV absolute abundance",
         "Relative abundance on sample",
         "Relative abundance to all samples",
         "Sample total abundance") %>% 
  tibble()

View(less_final_results)

## BLASTn ----

# select ASVs for BLASTn search 
asvs_blast_all <- less_final_results %>%
  pull("ASV (Sequence)") %>% 
  as.character() %>% 
  unique()

asvs_blast_all
length(asvs_blast_all)

# Visualizando os tamanhos das ASVs
asvs_blast_all %>% nchar() %>% table() %>% plot()

# fazendo para apenas 1 ASV

asvs_blast_4434 <- less_final_results %>%
  filter("ASV header" %in% ">ASV_44334_196bp")

asvs_blast_less <- asvs_blast_4434 %>%
  unique() %>% 
  pull("ASV (Sequence)") %>% 
  as.character()

## Usando o BLASTn no R em paralelo com o BLASTr

# BLASTr NT + 12slGCdb
{
  tictoc::tic("Parallel - Furrr 2 threads")
  Sys.time()
  
  blast_res_1 <- BLASTr::parallel_blast(
    # db_path = '"/data/databases/nt_jun2023/nt /home/gabriel/projetos/peixes-eDNA/databases/12sLGCdb/LGC12Sdb_complete_noGaps_2024-0403.fasta"',
    # db_path = '"/data/databases/nt_jun2023/nt /home/gabriel/projetos/peixes-eDNA/databases/12sLGCdb/old_db/LGC12Sdb_complete_noGaps.fasta"',
    # db_path = '"/data/databases/nt_jun2023/nt /home/gabriel/projetos/db-LGC/database/blastr/LGC12Sdb_complete_noGaps-2024-09-20.fasta"',
    db_path = '"/home/gabriel/projetos/db-LGC/database/blastr/LGC12Sdb_complete_noGaps-2024-09-20_FH.fasta"',
    asvs = asvs_blast_all[1:5],
    out_file = paste0(blast_path, "blast_out_res_1.csv"),
    out_RDS = paste0(blast_path, "blast_out_res_1.RDS"),
    total_cores = 70,
    perc_id = 80,
    num_threads = 2,
    perc_qcov_hsp = 80,
    num_alignments = 3,
    blast_type = "blastn"
  )
  
  blast_res_backup <- blast_res_1
  
  tictoc::toc()
  Sys.time()
  }


View(blast_res_1)

# Saving environment 
base::save.image(paste0(env_path, "env_", Sys.Date(), "_postBLAST.RData"))