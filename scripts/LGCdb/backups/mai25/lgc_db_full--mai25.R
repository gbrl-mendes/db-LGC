---
title: "Construção do LGC12Sdb - v2024"
author: "Igor Nascimento; Gabriel Mendes"
date: "05/2024"
---


##################################################################################
#                                                                                #
# Orientações para a construção de fastas a serem integrados no LGC12Sdb - v2024 #
#                                                                                #
#     O identificador deve estar no formato: >1111-Genero_especie-BACIA          #
#    O ID numérico deve ser composto de 4 dígitos e a bacia de 2 caractéres      #
#                                                                                #
#                 ex: >8125-Trichogenes_claviger-DC                              #
#                                                                                #
# Caso os indivíduos sequenciados não estejam depositados no banco de tecidos,   #
#   utilize um código único ao final do ID numérico para que não ocorra          #
#                    sobreposição de IDs depositados                             #
#                                                                                #
#                 ex: >0152A-Copella_nigrofasciata-AZ                            #
#                                                                                #
##################################################################################


# 1- Libraries needed ----
{
  library(dada2)
  library(Rcpp)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(spider)
  library(ape)
  library(ips)
  library(phyloseq)
  library(Biostrings)
  library(ShortRead)
  library(pheatmap)
  library(DECIPHER)
  library(ggdendro)
  library(readr)
  library(RColorBrewer)
  library(future)
  library(BLASTr)
}

# 2- Define directory ----

DB_folder <- "/home/igorhan/projetos/LGCdb/data/sql/"
fastas_folder <- "/home/igorhan/projetos/LGCdb/data/fastas/"
results_folder <- "/home/igorhan/projetos/LGCdb/LGC12Sdb/"

# Create db
LGCdb_full <- dbConnect(SQLite(), paste0(DB_folder, "LGC12Sdb_full.sql"))

# 3- Read and add new fastas ----
{
  # São Francisco
  Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/12S_full_unlign.fas",
          type = "FASTA",
          dbFile = LGCdb_full,
          identifier = "SF")

  # Jequi
  Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/alinhamentojequi57seq_ed_unalig.fasta",
          type = "FASTA",
          dbFile = LGCdb_full,
          identifier = "JQ")

  # Doce
  Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/12SDocedb.fasta",
          type = "FASTA",
          dbFile = LGCdb_full,
          identifier = "DC")

  # Amazônia - DESATUALIZADO
  # Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/fasta_amz.fas",
  #         type = "FASTA",
  #         dbFile = LGCdb_full,
  #         identifier = "AM")

  # Outras espécies JQ e SF - GRUPOS OUT
  # Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/new_Jq_species.fasta",
  #         type = "FASTA",
  #         dbFile = LGCdb_full,
  #         identifier = "JqSf")

  # Novas espécies + atualização db - TCC Igor
  Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/12Sdb_Seq3a_03nov21-10sp_curado.fas",
          type = "FASTA",
          dbFile = LGCdb_full,
          identifier = "JqSf1")

  # Novas espécies + atualização db - TCC Igor
  Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/12Sdb_Seq4b_24dez21-13sp_curado.fas",
          type = "FASTA",
          dbFile = LGCdb_full,
          identifier = "JqSf2")

  # Novas espécies + atualização db - TCC Igor
  Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/12Sdb_Seq5c_17jan22-9sp_curado.fas",
          type = "FASTA",
          dbFile = LGCdb_full,
          identifier = "JqSf3")

  # Trichogenes | Doce | Juliana/Igor
  Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/Trichogenes_12S.fasta",
          type = "FASTA",
          dbFile = LGCdb_full,
          identifier = "DC1")

  # Brycon howesi e B. vonoi
  Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/Seq12S_B_vonoi_B_howesi.fasta",
          type = "FASTA",
          dbFile = LGCdb_full,
          identifier = "PD")

  # Engraulidae
  Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/3_engraulidae.fasta",
          type = "FASTA",
          dbFile = LGCdb_full,
          identifier = "Eng")

  # Referências
  # Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/non_BRfish.fasta",
  #         type = "FASTA",
  #         dbFile = LGCdb_full,
  #         identifier = "Ref")

  # Surubim do Doce
  Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/surubim_dc_04mar2024.fasta",
          type = "FASTA",
          dbFile = LGCdb_full,
          identifier = "Sur")

  # AMZ 2 DESATUALIZADO
  # Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/amz_seq2.fasta",
  #         type = "FASTA",
  #         dbFile = LGCdb_full,
  #         identifier = "Am2")

  # Lycengraulis sp. - PERD
  Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/lycengraulis_PERD_TCP-fev2024.fasta",
          type = "FASTA",
          dbFile = LGCdb_full,
          identifier = "Eng2")
  
  # Seq AMZ completa versão 24/09/2024
  Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/stmpgmfasta6.fasta",
          type = "FASTA",
          dbFile = LGCdb_full,
          identifier = "AmS")
}

# 4- Calculate seq lengths ----
seq_lengths <- IdLengths(dbFile = LGCdb_full)

# Add seq lengths to DB
Add2DB(myData = seq_lengths, dbFile = LGCdb_full, verbose = TRUE)

BrowseDB(LGCdb_full)

# 5- Retrieve all seqs from DB ----

{
  dna_sf <- SearchDB(dbFile = LGCdb_full,identifier = "SF",nameBy = "description")
  dna_jq <- SearchDB(dbFile = LGCdb_full, identifier = "JQ",nameBy = "description")
  dna_dc <- SearchDB(dbFile = LGCdb_full, identifier = "DC",nameBy = "description")
  # dna_am <- SearchDB(dbFile = LGCdb_full, identifier = "AM",nameBy = "description")
  # dna_jqsf <- SearchDB(dbFile = LGCdb_full, identifier = "JqSf",nameBy = "description")
  dna_jqsf1 <- SearchDB(dbFile = LGCdb_full, identifier = "JqSf1",nameBy = "description")
  dna_jqsf2 <- SearchDB(dbFile = LGCdb_full, identifier = "JqSf2",nameBy = "description")
  dna_jqsf3 <- SearchDB(dbFile = LGCdb_full, identifier = "JqSf3",nameBy = "description")
  dna_trich <- SearchDB(dbFile = LGCdb_full, identifier = "DC1",nameBy = "description")
  dna_bryc <- SearchDB(dbFile = LGCdb_full, identifier = "PD",nameBy = "description")
  dna_eng <- SearchDB(dbFile = LGCdb_full, identifier = "Eng",nameBy = "description")
  # dna_ref <- SearchDB(dbFile = LGCdb_full, identifier = "Ref",nameBy = "description")
  dna_sur <- SearchDB(dbFile = LGCdb_full, identifier = "Sur",nameBy = "description")
  # dna_am2 <- SearchDB(dbFile = LGCdb_full, identifier = "Am2",nameBy = "description")
  dna_eng2 <- SearchDB(dbFile = LGCdb_full, identifier = "Eng2",nameBy = "description")
  dna_ams <- SearchDB(dbFile = LGCdb_full, identifier = "AmS",nameBy = "description")
}

# Merge seqs set
dna_all <- c(dna_sf,
             dna_jq, 
             dna_dc,
             # dna_am,
             # dna_jqsf, 
             dna_jqsf1, 
             dna_jqsf2,
             dna_jqsf3,
             dna_trich,
             dna_bryc, 
             dna_eng, 
             # dna_ref,
             dna_sur,
             # dna_am2,
             dna_eng2,
             dna_ams
)

# 6- create and assign new names ----

DB_tbl <- tibble("Original names" = names(dna_all),
                 "Identifier" = as.character(""),
                 "Names" = as.character(""),
                 "Cluster" = as.numeric(""),
                 "Composed name"= as.character(""),
                 # "Final Identifier" = as.character(""),
                 "River basin"= as.character(""),
                 "genus"= as.character("")) # it will be used to search taxonomy on NCBI

# Identify basin from origin file
for (seq in 1:nrow(DB_tbl)) {
  DB_tbl$`River basin`[seq] <- str_split_fixed(string = DB_tbl$`Original names`[seq],pattern = "-",n = 3)[3]
}

# Correct species name
for (seq in 1:nrow(DB_tbl)) {
  
  DB_tbl$Identifier[seq] <-  str_split_fixed(string = DB_tbl$`Original names`[seq], pattern = "-", n = 2)[1]
  DB_tbl$Names[seq] <- str_split_fixed(string = DB_tbl$`Original names`[seq], pattern = "-", n = 2)[2] %>%
    str_replace_all(pattern = "_",replacement = " ") %>%
    # str_remove_all(pattern = "aff.") %>% # aff, cf e gr são informativos para curadoria, não remover!!
    str_remove_all(pattern = "-.*.$")
  DB_tbl$genus[seq] <- str_split(string = DB_tbl$Names[seq],pattern = " ")[[1]][1]
  
}

# Rewrite identifiers with 2024 model
# names(dna_all) %>% sort()
# which(names(dna_all) == "LGC_AZ0026A Moenkhausia celibela")
# 
# 
# 
# as.character(dna_all[401]) == as.character(dna_all[402])

for (seq in 1:nrow(DB_tbl)) {
  
  for (seq2 in 1:length(dna_all)) {
    
    if (names(dna_all)[seq2] == DB_tbl$`Original names`[seq]) {
      
      names(dna_all)[seq2] <- paste0("LGC_", DB_tbl$`River basin`[seq], DB_tbl$Identifier[seq], " ", DB_tbl$Names[seq])
      # %>% 
        # str_replace_all(pattern = "EC", replacement = "") %>% 
        # str_replace_all(pattern = "AMZ", replacement = "")
      DB_tbl$`Composed name`[seq] <- names(dna_all)[seq2] }
    
  }
}


DB_tbl %>% 
  select(`Composed name`) %>% unique() %>% View()

# 7- rewriting DB fasta with corrected names ----
writeXStringSet(x = dna_all,
                filepath = paste0(results_folder, "LGC12Sdb_new_names_unaligned_no_trim-",
                                  Sys.Date(), format = ".fasta"))

# 8- align DB seqs ----

# Must remove gaps and orient nucleotides, before proceed

dna_all <- RemoveGaps(dna_all)
dna_all <- OrientNucleotides(dna_all)

dna_aligned <- AlignSeqs(myXStringSet = dna_all, refinements = 700,iterations = 700, verbose = TRUE)

# 9- Cluster and Distance matrix, for dendrograms ----

## Distance matrix ----

# Select core alignment - corta o alinhamento
BrowseSeqs(dna_aligned, colorPatterns = FALSE)

dna_alg_sub <- subseq(x = dna_aligned,start = 21, end = 825)
BrowseSeqs(dna_alg_sub, colorPatterns = FALSE)

## Clusters ----

# Identify clusters for primer design
dna_alg_sub_dist <- DistanceMatrix(myXStringSet = dna_alg_sub,
                                   includeTerminalGaps = FALSE,
                                   correction = "Jukes-Cantor",
                                   processors = 20)
dim(dna_alg_sub_dist) # a symmetric matrix

pheatmap(mat = dna_alg_sub_dist,
         fontsize = 4,
         color = rev(brewer.pal(n = 30,name = "Greens")))

# Identify clusters based on matrix distance

dna_clust <- DECIPHER::TreeLine(
  myDistMatrix=dna_alg_sub_dist,
  method = "UPGMA",
  cutoff = 0.05, # use `cutoff = 0.03` for a 97% OTU
  type = "both",
  show = TRUE)

head(dna_clust) # cluster numbers

dev.off() #deleta os dendrogramas adicionados#

# Add clusters info to DB

Add2DB(data.frame(
  cluster=dna_clust[[1]]$cluster),
  dbFile = LGCdb_full,
  verbose = TRUE)

BrowseDB(LGCdb_full, orderBy = "cluster")

# Add clusters to DB_tbl

for (seq in 1:nrow(DB_tbl)) {
  for (seq2 in 1:nrow(dna_clust[[1]])) {
    
    if ((rownames(dna_clust[[1]])[seq2]) == (DB_tbl$`Composed name`[seq])) {
      DB_tbl$Cluster[seq] <- (dna_clust[[1]][seq2,])
    }
  }
}

# Order seqs by cluster
dna_order <- dna_alg_sub[order(DB_tbl$Cluster),]

BrowseSeqs(myXStringSet = dna_order,
           colorPatterns = FALSE)

# 10- FUNCTION to retrieve tax ranks using taxid ----

## Extraindo os generos

db_genus <- c(DB_tbl$genus)

db_genus

gen2search <- db_genus %>% unique() %>% na.omit() %>% as.character() %>% sort()

gen2search %>% sort()

#get taxID from name----

get_tax_by_name <- function(org_name){

  if (stringr::str_detect(org_name, "[^a-zA-Z ]")) {
    taxID <- ""
    message(glue::glue("\t\tSpecial character detected :\t\t{org_name}"))
    return(taxID)
  }

  taxID <- system(command = paste0("esearch -db taxonomy -query '",org_name,"' | efetch -format uid"), intern = TRUE)

  message(glue::glue("Done for:\t\t{org_name}"))

  if (length(taxID) >1) {
    taxID <- taxID[[1]]
  }

  if (identical(taxID,character(0))) {
    taxID <- ""
    message(glue::glue("\t\tProblem for:\t\t{org_name}"))
    }

  message(paste0(org_name, "  ",taxID))
  if (is.na(taxID)) {

    return("")
  }
  return(taxID)
}
###################________


gen2search[!gen2search %in% trad_tax_miss$Genus]


#so da primeira vez
            trad_tax_miss <- tibble("Genus" = gen2search) %>%
              dplyr::rowwise() %>%
              dplyr::mutate("query_taxID" = get_tax_by_name(`Genus`)) %>%
              dplyr::ungroup() 
            
            
            
            trad_tax_miss_2 <- tibble("Genus" = gen2search[!gen2search %in% trad_tax_miss$Genus]) %>%
              dplyr::rowwise() %>%
              dplyr::mutate("query_taxID" = get_tax_by_name(`Genus`)) %>%
              dplyr::ungroup()
            

#só para os que faltam
trad_tax_miss <- trad_tax_miss %>%
  dplyr::filter(!query_taxID == "") %>% 
  dplyr::bind_rows(trad_tax_miss_2)


trad_tax_miss$query_taxID %>% sort()


trad_taxIDs2search <- trad_tax_miss$query_taxID %>%
  BiocGenerics::unique() 

# Extract taxonomy

taxonomy_df <- parallel_get_tax(organisms_taxIDs = trad_taxIDs2search,
                                # total_cores = 20,
                 retry_times = 5,
                 verbose = F)


taxonomy_df_backup <- taxonomy_df

## Filtrar apenas clados informativos
# taxonomy_df <- taxonomy_df %>% 
#   filter(!Rank %in% c("no rank","clade")) %>% 
#   unique()

## Criar nova tabela de taxonomia wider
# taxonomy_tbl_df <- taxonomy_df %>% 
#   dplyr::select(-c("TaxId", "query_taxID")) %>%
#   unique() %>%
#   mutate("genus" = Sci_name ) %>% 
#   filter(Rank %in% c("superkingdom","kingdom","phylum","subphylum","class","subclass","order","suborder","family","subfamily","genus")) %>% 
#   tidyr::pivot_wider(
#     id_cols = c(genus,Sci_name),
#     names_from = Rank,
#     values_from = c(ScientificName),names_repair = "minimal") %>%
#   relocate("Sci_name","genus","superkingdom","kingdom","phylum","subphylum","class","subclass","order","suborder","family","subfamily","genus")

# View(taxonomy_tbl_df)

taxonomy_tbl_df <- taxonomy_df

taxonomy_tbl_df %>% colnames() %>% paste0(collapse = "\n") %>% cat()

## fill NA tax with combination of max_tax and rank
for (line in 1:nrow(taxonomy_tbl_df)) {
  # if (taxonomy_tbl_df$genus[line] %in% c("NA",NA,"")) {
  if (taxonomy_tbl_df$`Superkingdom (NCBI)`[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_df$`Superkingdom (NCBI)`[line] <- paste0("superkingdom of ", taxonomy_tbl_df$`Kingdom (NCBI)`[line]) }
  
  if (taxonomy_tbl_df$`Kingdom (NCBI)`[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_df$`Kingdom (NCBI)`[line] <- paste0("kingdom of ", taxonomy_tbl_df$`Superkingdom (NCBI)`[line]) }
  
  if (taxonomy_tbl_df$`Phylum (NCBI)`[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_df$`Phylum (NCBI)`[line] <- paste0("phylum of ", taxonomy_tbl_df$`Kingdom (NCBI)`[line]) }
  
  if (taxonomy_tbl_df$`Subphylum (NCBI)`[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_df$`Subphylum (NCBI)`[line] <- paste0("subphylum of ", taxonomy_tbl_df$`Phylum (NCBI)`[line]) }
  
  if (taxonomy_tbl_df$`Class (NCBI)`[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_df$`Class (NCBI)`[line] <- paste0("class of ", taxonomy_tbl_df$`Subphylum (NCBI)`[line]) }
  
  if (taxonomy_tbl_df$`Subclass (NCBI)`[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_df$`Subclass (NCBI)`[line] <- paste0("subclass of ", taxonomy_tbl_df$`Class (NCBI)`[line]) }
  
  if (taxonomy_tbl_df$`Order (NCBI)`[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_df$`Order (NCBI)`[line] <- paste0("order of ", taxonomy_tbl_df$`Subclass (NCBI)`[line]) }
  
  if (taxonomy_tbl_df$`Suborder (NCBI)`[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_df$`Suborder (NCBI)`[line] <- paste0("suborder of ", taxonomy_tbl_df$`Order (NCBI)`[line]) }
  
  if (taxonomy_tbl_df$`Family (NCBI)`[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_df$`Family (NCBI)`[line] <- paste0("family of ", taxonomy_tbl_df$`Suborder (NCBI)`[line]) }
  
  if (taxonomy_tbl_df$`Subfamily (NCBI)`[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_df$`Subfamily (NCBI)`[line] <- paste0("subfamily of ", taxonomy_tbl_df$`Family (NCBI)`[line]) }
  
  if (is.na(taxonomy_tbl_df$`Genus (NCBI)`[line])) {
    taxonomy_tbl_df$`Genus (NCBI)`[line] <- paste0("genus of ", taxonomy_tbl_df$`Subfamily (NCBI)`[line]) }
  # 
}

View(taxonomy_tbl_df)

taxonomy_tbl_df <- taxonomy_tbl_df %>% 
  dplyr::rename("genus" = "Sci_name")

# save file
write.csv(x = taxonomy_tbl_df,
          file = "/home/igorhan/projetos/LGCdb/data/db_tbl/db_taxonomy.csv",
          row.names = F)

# join taxonomy with db_tbl
DB_tbl_final <- left_join(x = DB_tbl,
                          y = taxonomy_tbl_df,
                          by = "genus")

View(DB_tbl_final)

# melhorando o output do DB_tbl_final
DB_tbl_final <- DB_tbl_final %>% 
  select(-c(`Genus (NCBI)`, `Original names`, `query_taxID`)) %>%
  dplyr::rename("subfamily" = "Subfamily (NCBI)",
                "family" = "Family (NCBI)", 
                "suborder" = "Suborder (NCBI)",
                "order" = "Order (NCBI)", 
                "subclass" = "Subclass (NCBI)", 
                "class" = "Class (NCBI)", 
                "subphylum" = "Subphylum (NCBI)", 
                "phylum" = "Phylum (NCBI)", 
                "kingdom" = "Kingdom (NCBI)", 
                "superkingdom" = "Superkingdom (NCBI)") %>% 
  relocate(`Composed name`, `River basin`, Identifier, Names, genus, subfamily, family, 
           suborder, order, subclass, class, subphylum, phylum, kingdom, superkingdom, Cluster)


# 11- Add sequence to db_tbl ----

DB_seqs_df <-
  dna_all %>%
  as.data.frame() %>%
  rownames_to_column("Composed name") %>% 
  dplyr::rename("Sequence" = "x") %>%
  relocate("Composed name", Sequence) %>% 
  tibble()

View(DB_seqs_df)

DB_tbl_final_seqs <- DB_tbl_final %>% 
  left_join(DB_seqs_df, by = "Composed name") %>% 
  relocate(`Composed name`, `Sequence`, `River basin`, Identifier, Names, genus, subfamily, family, 
           suborder, order, subclass, class, subphylum, phylum, kingdom, superkingdom, Cluster)

View(DB_tbl_final_seqs)

# save file
write.csv(x = DB_tbl_final,
          file = "/home/igorhan/projetos/LGCdb/data/db_tbl/DB_tbl_final.csv",
          row.names = F)

# 12- Write DB ordered by cluster to fasta file ----

# No Gaps
dna_order_noGaps <- RemoveGaps(dna_order)
writeXStringSet(x = dna_order_noGaps,
                filepath = paste0(results_folder,"LGC12Sdb_complete_noGaps-", Sys.Date(), format = ".fasta"))

# With Gaps
writeXStringSet(x = dna_order,
                filepath = paste0(results_folder,"LGC12Sdb_complete_withGaps-", Sys.Date(), format = ".fasta"))

# 13- Create nwk dendrograms ----

LGC12Sdb_tree <- ape::njs(dna_alg_sub_dist)

LGC12Sdb_tree$tip.label

ape::write.tree(phy = LGC12Sdb_tree,
                digits = 5,
                file = paste0(results_folder,"LGC12Sdb_NJtree-", Sys.Date(), format = ".nwk"))

# 14- Check for problems ----

# ######################## ATENÇÃO!!! ########################

# os conflitos foram removidos diretamente no arquivo .fasta, pq estava dando
# erro em outras análises pelo pipe, caso utilize os arquivos fasta localizados em qualquer
# lugar que nao seja em "/home/igorhan/projetos/LGCdb/data/fastas/", será necessário
# realizar alterações!!

# conferir se contém identificador duplicados e grupos out
names(dna_order) %>% as_tibble() %>% View() 


# 15- LGC12Sdb 2024 Final ----

View(DB_tbl_final_seqs)
  
DB_tbl_final_seqs <-
    DB_tbl_final_seqs %>% 
    mutate(`Composed name` = case_when(`Composed name` == "LGC_JQ231 Steindachnerina elegans" ~ "LGC_JQ231 Steindachnerina elegans ERRADO!!!!!!!!!!",
                                       `Composed name` == "LGC_JQ6586 Steindachnerina elegans" ~ "LGC_JQ6586 Steindachnerina elegans ERRADO!!!!!!!!!!",
                                       `Composed name` == "LGC_SF1284 Steindachnerina elegans" ~ "LGC_SF1284 Steindachnerina elegans ERRADO!!!!!!!!!!",
                                       TRUE ~ `Composed name`))

## 15.1- write tbl xlsx ----
writexl::write_xlsx(x = DB_tbl_final_seqs,
           path = paste0("~/projetos/LGCdb/LGC12Sdb/","LGC12Sdb_tbl","-", Sys.Date(), ".xlsx"),
           col_names = T,
           format_headers = T)

LGC12Sdb_tax <- read.csv(file = "/data/databases/LGC12Sdb/mai25/LGC12Sdb--db_taxonomy--mai25.csv",
                         sep = ",",
                         header = T,
                         check.names = F) %>% as_tibble()

LGC12Sdb_tbl <- read.csv(file = "/data/databases/LGC12Sdb/mai25/LGC12Sdb--DB_tbl_final--mai25.csv",
                         sep = ";",
                         header = T,
                         check.names = F) %>% as_tibble()


LGC12Sdb_all <- LGC12Sdb_tax %>% left_join(LGC12Sdb_tbl,
                                           by = "genus")

# tem que salvar uma tabela com o mesmo primeiro nome antes do espaço e com o query_taxID.
# esse arquivo tem que ter o nome exatamente igual ao do banco (sempre usamos o no gaps) e tem que terminar com .txt
# isso é necessário pra formatar o db do BLAST

IDS_and_tax <- LGC12Sdb_all %>%
  dplyr::select(`Composed name`,query_taxID) %>%
  dplyr::mutate(`Composed name` = stringr::str_remove_all(`Composed name`, pattern = " .*$")) %>%
  dplyr::mutate(`Composed name` = stringr::str_replace(`Composed name`, pattern = "^",replacement = "lcl|")) %>%
  BiocGenerics::unique() %>%
  tidyr::unite("Junto", `Composed name`,query_taxID,sep = " ") %>%
  dplyr::pull("Junto")

IDS_and_tax %>%
  write_lines(file = "~/projetos/LGCdb/LGC12Sdb/LGC12Sdb_complete_noGaps-2025-05-21.txt",
              sep = "\n")

#aqui explica como formatar o banco do blast
# https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.makeblastdb_application_opt/


# os arquivos corretos estãoa qui

# /data/databases/LGC12Sdb/mai25

# 16- plot dendrogram db ----
  
db_doce <- c(dna_dc, dna_trich,dna_sur) %>% 
  unique()

db_doce <- db_doce %>% 
  RemoveGaps() %>% 
  DECIPHER::OrientNucleotides() %>%
  sort(decreasing = T)


db_doce_alg <- msa::msa(inputSeqs = db_doce,
            method = "ClustalW",
            substitutionMatrix = "clustalw",
            type = "dna",
            order = "input")


## 16.1- calculate distance matrix ----
#convert object to phy
db_doce_alg_phy <- phangorn::as.phyDat(db_doce_alg,
                             type = "DNA",
                             names = names(db_doce_alg))
#calculate dist
ASVs_seqs_dist <- dist.ml(x = db_doce_alg_phy,
                          model = "JC69")

## 16.2- built tree ----
ASVs_tree <- phangorn::NJ(x = ASVs_seqs_dist) # Note, tip order != sequence order
ASVs_tree


## 16.3- sequence, names and tips ----
tips_labels <- db_doce_alg %>%
  as("DNAStringSet") %>% 
  names()

add.tips(tree = ASVs_tree, tips = tips_labels,where = 10, edge.length = NULL)


# fit = pml(ASVs_tree, data = ASVs_seqs_align_phy)

## 16.4- write/read tree ----
ape::write.tree(phy = ASVs_tree,
                file = "/home/igorhan/projetos/Doce12Sdb/DoceDB_ASV_tree.nwk")

ASVs_tree <- read.tree(file = "/home/igorhan/projetos/Doce12Sdb/DoceDB_ASV_tree.nwk")



## 16.5- add labels ----

# read metadata table
metadata_tbl <- DB_tbl_final_seqs


tips_metadata <- metadata_tbl %>%
  tidyr::unite(col = "tip_label",Identifier,Names, `River basin`,sep = "-",remove = F) %>%
  mutate("tip_label" = stringr::str_replace_all(tip_label,pattern = " ",replacement = "_")) %>%
  relocate("tip_label" )

# rownames(tips_metadata) <- tips_metadata$node

## 16.6- plot tree ----

options(ignore.negative.edge = TRUE)

#set anotations position relative to tip
tip_alignment <- FALSE


tree_plot <- ggtree(tr = ASVs_tree,
                    branch.length = 3,
                    ladderize = T)  %<+%
  tips_metadata +                           #adicionar metadados
  geom_tiplab(align = tip_alignment,
              linesize = 0.5)  +
  geom_treescale(width = 0.4)  +
  theme_tree2() +
  
  
  ## `Curated ID`              #####################
geom_tiplab(
  aes(label = `Names`, col = `Names`),
  offset = 0.027,
  # size = 3,
  linetype = "blank" ,
  geom = "text",
  align = tip_alignment) +
  scale_color_manual(values = viridis::turbo(n = 10)) +


  ## Order      #####################
geom_tiplab(
  aes(label = `order`, col = `order`),
  offset = 0.09,
  # size = 3,
  linetype = "blank" ,
  geom = "text",
  align = tip_alignment) +
  scale_color_manual(values = viridis::turbo(n = 10)) +

  ## Family      #####################
geom_tiplab(
  aes(label = `family`, col = `family`),
  offset = 0.12,
  # size = 3,
  linetype = "blank" ,
  geom = "text",
  align = tip_alignment) +
  scale_color_manual(values = viridis::turbo(n = 105)) +
  guides(col = "none")


dev.off()
tree_plot


## 16.7- save tree plot ----

ggsave(file = paste0("/home/igorhan/projetos/docedb_tree.pdf",
                     collapse = ""),
       plot = tree_plot,
       device = "pdf",
       width = 40,
       height = 100,
       limitsize = FALSE,
       units = "cm",
       dpi = 300)



# referencias
#https://yulab-smu.top/treedata-book/chapter7.html
#https://www.youtube.com/watch?v=3swFCSt2_x4
#https://f1000research.com/articles/5-1492/v1
#https://www.molecularecologist.com/2017/02/08/phylogenetic-trees-in-r-using-ggtree/
#https://epirhandbook.com/en/phylogenetic-trees-1.html
#https://boopsboops.blogspot.com/2010/10/negative-branch-lengths-in-neighbour.html
