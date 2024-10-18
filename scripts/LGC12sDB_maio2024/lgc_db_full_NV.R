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
}

# 2- Define directory ----

DB_folder <- "/home/igorhan/projetos/LGCdb/data/sql/"
fastas_folder <- "/home/igorhan/projetos/LGCdb/data/fastas/"
results_folder <- "/home/igorhan/projetos/LGCdb/LGC12Sdb/"

# Create db
LGCdb_full <- dbConnect(SQLite(), paste0(DB_folder, "LGC12Sdb_full.sql"))

# 3- Read and add new fastas ----
# {
#   # São Francisco
#   Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/12S_full_unlign.fas",
#           type = "FASTA",
#           dbFile = LGCdb_full,
#           identifier = "SF")
#   
#   # Jequi
#   Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/alinhamentojequi57seq_ed_unalig.fasta",
#           type = "FASTA",
#           dbFile = LGCdb_full,
#           identifier = "JQ")
#   
#   # Doce
#   Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/12SDocedb.fasta",
#           type = "FASTA",
#           dbFile = LGCdb_full,
#           identifier = "DC")
#   
#   # Amazônia
#   Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/fasta_amz.fas",
#           type = "FASTA",
#           dbFile = LGCdb_full,
#           identifier = "AM")
#   
#   # Outras espécies JQ e SF - GRUPOS OUT
#   # Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/new_Jq_species.fasta",
#   #         type = "FASTA",
#   #         dbFile = LGCdb_full,
#   #         identifier = "JqSf")
#   
#   # Novas espécies + atualização db - TCC Igor
#   Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/12Sdb_Seq3a_03nov21-10sp_curado.fas",
#           type = "FASTA",
#           dbFile = LGCdb_full,
#           identifier = "JqSf1")
#   
#   # Novas espécies + atualização db - TCC Igor
#   Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/12Sdb_Seq4b_24dez21-13sp_curado.fas",
#           type = "FASTA",
#           dbFile = LGCdb_full,
#           identifier = "JqSf2")
#   
#   # Novas espécies + atualização db - TCC Igor
#   Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/12Sdb_Seq5c_17jan22-9sp_curado.fas",
#           type = "FASTA",
#           dbFile = LGCdb_full,
#           identifier = "JqSf3")
#   
#   # Trichogenes | Doce | Juliana/Igor
#   Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/Trichogenes_12S.fasta",
#           type = "FASTA",
#           dbFile = LGCdb_full,
#           identifier = "DC1")
#   
#   # Brycon howesi e B. vonoi
#   Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/Seq12S_B_vonoi_B_howesi.fasta",
#           type = "FASTA",
#           dbFile = LGCdb_full,
#           identifier = "PD")
#   
#   # Engraulidae
#   Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/3_engraulidae.fasta",
#           type = "FASTA",
#           dbFile = LGCdb_full,
#           identifier = "Eng")
#   
#   # Referências
#   # Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/non_BRfish.fasta",
#   #         type = "FASTA",
#   #         dbFile = LGCdb_full,
#   #         identifier = "Ref")
#   
#   # Surubim do Doce
#   Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/surubim_dc_04mar2024.fasta",
#           type = "FASTA",
#           dbFile = LGCdb_full,
#           identifier = "Sur")
#   
#   # AMZ 2
#   Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/amz_seq2.fasta",
#           type = "FASTA",
#           dbFile = LGCdb_full,
#           identifier = "Am2")
#   
#   # Lycengraulis sp. - PERD
#   Seqs2DB(seqs = "/home/igorhan/projetos/LGCdb/data/fastas/lycengraulis_PERD_TCP-fev2024.fasta",
#           type = "FASTA",
#           dbFile = LGCdb_full,
#           identifier = "Eng2")
#   
# }

# 4- Calculate seq lengths ----
# seq_lengths <- IdLengths(dbFile = LGCdb_full)
# 
# # Add seq lengths to DB
# Add2DB(myData = seq_lengths, dbFile = LGCdb_full, verbose = TRUE)
# 
# BrowseDB(LGCdb_full)

# 5- Retrieve all seqs from DB ----

{
  dna_sf <- SearchDB(dbFile = LGCdb_full,identifier = "SF",nameBy = "description")
  dna_jq <- SearchDB(dbFile = LGCdb_full, identifier = "JQ",nameBy = "description")
  dna_dc <- SearchDB(dbFile = LGCdb_full, identifier = "DC",nameBy = "description")
  dna_am <- SearchDB(dbFile = LGCdb_full, identifier = "AM",nameBy = "description")
  # dna_jqsf <- SearchDB(dbFile = LGCdb_full, identifier = "JqSf",nameBy = "description")
  dna_jqsf1 <- SearchDB(dbFile = LGCdb_full, identifier = "JqSf1",nameBy = "description")
  dna_jqsf2 <- SearchDB(dbFile = LGCdb_full, identifier = "JqSf2",nameBy = "description")
  dna_jqsf3 <- SearchDB(dbFile = LGCdb_full, identifier = "JqSf3",nameBy = "description")
  dna_trich <- SearchDB(dbFile = LGCdb_full, identifier = "DC1",nameBy = "description")
  dna_bryc <- SearchDB(dbFile = LGCdb_full, identifier = "PD",nameBy = "description")
  dna_eng <- SearchDB(dbFile = LGCdb_full, identifier = "Eng",nameBy = "description")
  # dna_ref <- SearchDB(dbFile = LGCdb_full, identifier = "Ref",nameBy = "description")
  dna_sur <- SearchDB(dbFile = LGCdb_full, identifier = "Sur",nameBy = "description")
  dna_am2 <- SearchDB(dbFile = LGCdb_full, identifier = "Am2",nameBy = "description")
  dna_eng2 <- SearchDB(dbFile = LGCdb_full, identifier = "Eng2",nameBy = "description")
}

# Merge seqs set
dna_all <- c(dna_sf,
             dna_jq, 
             dna_dc,
             dna_am,
             # dna_jqsf, 
             dna_jqsf1, 
             dna_jqsf2,
             dna_jqsf3,
             dna_trich,
             dna_bryc, 
             dna_eng, 
             # dna_ref,
             dna_sur,
             dna_am2,
             dna_eng2
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

for (seq in 1:nrow(DB_tbl)) {
  
  for (seq2 in 1:length(dna_all)) {
    
    if (names(dna_all)[seq2] == DB_tbl$`Original names`[seq]) {
      
      names(dna_all)[seq2] <- paste0("LGC_", DB_tbl$`River basin`[seq], DB_tbl$Identifier[seq], " ", DB_tbl$Names[seq]) %>% 
        str_replace_all(pattern = "EC", replacement = "") %>% 
        str_replace_all(pattern = "AMZ", replacement = "")
      DB_tbl$`Composed name`[seq] <- names(dna_all)[seq2] }
    
  }
}

DB_tbl %>% 
  select(`Composed name`) %>% unique() %>% View()

# 7- rewriting DB fasta with corrected names ----
# writeXStringSet(x = dna_all,filepath = paste0(results_folder,"LGC12Sdb_new_names_unaligned_no_trim-", Sys.Date(), format = ".fasta"))

# 8- align DB seqs ----

# Must remove gaps and orient nucleotides, before proceed

dna_all <- RemoveGaps(dna_all)
dna_all <- OrientNucleotides(dna_all)

dna_aligned <- AlignSeqs(myXStringSet = dna_all, refinements = 700,iterations = 700, verbose = TRUE)

# BrowseSeqs(dna_aligned, colorPatterns = F)

# 9- Cluster and Distance matrix, for dendrograms ----

## Distance matrix ----

# Select core alignment - corta o alinhamento
BrowseSeqs(dna_aligned, colorPatterns = FALSE)

dna_alg_sub <- subseq(x = dna_aligned,start = 22, end = 825)
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

# 10- FUNCTION to retrieve tax ranks using genus ----
source("/home/igorhan/projetos/Doce12Sdb/src/extract_taxonomy_name-copy.R")

# example
extract_taxonomy_name("Bunocephalus")

# extract
db_genus <- c(DB_tbl$genus)

db_genus

taxIDs2search <- db_genus %>% unique() %>% na.omit() %>% as.character() %>% sort()

taxIDs2search

#buscando as classificações

taxonomy_df_db <- tibble()

for (genero in 1:length(taxIDs2search[1:4])) {
  
  taxonomy_df_gen_db <- extract_taxonomy_name(taxIDs2search[genero])
  
  taxonomy_df_db <- bind_rows(taxonomy_df_db,taxonomy_df_gen_db)
  
}

# lista quais generos n foram encontrados
taxIDs2search[!taxIDs2search %in% unique(taxonomy_df_db$max_tax)]


# repete a busca para os nao encontrados
for (genero in which(!taxIDs2search %in% unique(taxonomy_df_db$max_tax))) {
  
  taxonomy_df_gen_db <- extract_taxonomy_name(taxIDs2search[genero])
  
  taxonomy_df_db <- bind_rows(taxonomy_df_db,taxonomy_df_gen_db)
  
}

taxonomy_df_db

# taxonomy_df_bckp_db <- taxonomy_df_db
# taxonomy_df_db <- taxonomy_df_bckp_db


# filter tbl for informative clades only
taxonomy_df_db <- taxonomy_df_db %>% 
  filter(!Rank %in% c("no rank","clade")) %>% 
  unique()

taxonomy_tbl_db <- taxonomy_df_db %>% 
  dplyr::select(-c("TaxId")) %>%
  unique() %>%
  rename("genus" = "max_tax") %>% 
  # dplyr::filter(Rank %in% c("kingdom","phylum","class","order","family")) %>%
  filter(Rank %in% c("superkingdom","kingdom","phylum","subphylum","class","subclass","order","suborder","family","subfamily","genus")) %>% 
  tidyr::pivot_wider(
    id_cols = c(genus,Sci_name),
    names_from = Rank,
    values_from = c(ScientificName),names_repair = "minimal") %>%
  # tidyr::pivot_wider(id_cols = c(genus,Sci_name),
  #                    names_from = Rank,
  #                    values_from = c(ScientificName,TaxId),
  #                    names_repair = "minimal") %>% View()
  # dplyr::select(genus,dplyr::starts_with("Scie")) %>% 
  relocate("Sci_name","genus","superkingdom","kingdom","phylum","subphylum","class","subclass","order","suborder","family","subfamily","genus")


# fill NA tax with combination of max_tax and rank

for (line in 1:nrow(taxonomy_tbl_db)) {
  # if (taxonomy_tbl_db$genus[line] %in% c("NA",NA,"")) {
  if (taxonomy_tbl_db$superkingdom[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_db$superkingdom[line] <- paste0("superkingdom of ", taxonomy_tbl_db$kingdom[line]) }
  
  if (taxonomy_tbl_db$kingdom[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_db$kingdom[line] <- paste0("kingdom of ", taxonomy_tbl_db$superkingdom[line]) }
  
  if (taxonomy_tbl_db$phylum[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_db$phylum[line] <- paste0("phylum of ", taxonomy_tbl_db$kingdom[line]) }
  
  if (taxonomy_tbl_db$subphylum[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_db$subphylum[line] <- paste0("subphylum of ", taxonomy_tbl_db$phylum[line]) }
  
  if (taxonomy_tbl_db$class[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_db$class[line] <- paste0("class of ", taxonomy_tbl_db$subphylum[line]) }
  
  if (taxonomy_tbl_db$subclass[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_db$subclass[line] <- paste0("subclass of ", taxonomy_tbl_db$class[line]) }
  
  if (taxonomy_tbl_db$order[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_db$order[line] <- paste0("order of ", taxonomy_tbl_db$subclass[line]) }
  
  if (taxonomy_tbl_db$suborder[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_db$suborder[line] <- paste0("suborder of ", taxonomy_tbl_db$order[line]) }
  
  if (taxonomy_tbl_db$family[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_db$family[line] <- paste0("family of ", taxonomy_tbl_db$suborder[line]) }
  
  if (taxonomy_tbl_db$subfamily[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_db$subfamily[line] <- paste0("subfamily of ", taxonomy_tbl_db$family[line]) }
  
  if (is.na(taxonomy_tbl_db$genus[line])) {
    taxonomy_tbl_db$genus[line] <- paste0("genus of ", taxonomy_tbl_db$subfamily[line]) }
  # 
}


# save file
# write.csv(x = taxonomy_tbl_db,
#           file = "/home/igorhan/projetos/LGCdb/data/db_tbl/db_taxonomy.csv", row.names = F)

# join taxonomy with db_tbl
DB_tbl_final <- left_join(x = DB_tbl,
                          y = taxonomy_tbl_db,
                          by = "genus")

# melhorando o output do DB_tbl_final
DB_tbl_final <- DB_tbl_final %>% 
  select(-c(Sci_name, `Original names`)) %>%
  relocate(`Composed name`, `River basin`, Identifier, Names, genus, subfamily, family, 
           suborder, order, subclass, class, subphylum, phylum, kingdom, superkingdom, Cluster)

# 10a (Gabriel M.) Retrieve TaxIDs ----

## Função para recuperar o TaxID
retrieve_taxid <- function(organism_name) {
  Sys.sleep(0.5)
  result <- tryCatch({
    myTAI::taxonomy(
      organism = organism_name, 
      db = "ncbi", 
      output = "taxid"
    )
  }, error = function(e) {
    NULL
  })
  
  if (!is.null(result) && !inherits(result, "try-error")) {
    return(result$id[1])
  } else {
    return(NA)
  }
}

## Extraindo os generos
db_genus <- c(DB_tbl$genus)

db_genus

gen2search <- db_genus %>% unique() %>% na.omit() %>% as.character() %>% sort()

gen2search

## Loop para recuperar TaxIDs usando Genero

taxid_gen_db <- tibble() %>%
  mutate("taxid" = as.character(),
         "genero" = as.character())

for (genero in 1:length(gen2search)) {
  
  taxid_gen_db_tax <- tibble(taxid = ifelse(is.na(retrieve_taxid(gen2search[genero])), NA, as.character(retrieve_taxid(gen2search[genero]))),
                             genero = gen2search[genero])
  
  taxid_gen_db <- bind_rows(taxid_gen_db,taxid_gen_db_tax)
  
  # Verifica se o valor de taxid é NA
  if (is.na(taxid_gen_db_tax$taxid)) {
    # Se for NA, busca novamente o gênero
    genero_sem_taxid <- taxid_gen_db_tax$genero
    taxid_novo <- retrieve_taxid(genero_sem_taxid)
    
    # Atualiza a tabela com o novo valor de taxid
    taxid_gen_db <- taxid_gen_db %>%
      mutate(taxid = ifelse(genero == genero_sem_taxid, taxid_novo, taxid))
  }
}

# 10b (Gabriel M.) Retrieve tax ranks ----
library()
source("/home/heron/prjcts/ecomol/R/extract_taxonomy_taxID.R")
extract_taxonomy_taxID("42548")

## Loop para recuperar as taxonomias utilizando os TaxIds
plan(future::multisession(workers = 78))

## Função para extrair a taxonomia de um conjunto de taxIDs e repetir a busca até que todos os taxIDs sejam encontrados
search_tax <- function(taxIDs) {
  taxIDs <- unique(na.omit(as.character(taxIDs)))
  taxonomy_df <- furrr::future_map_dfr(taxIDs, extract_taxonomy_taxID, .options = furrr::furrr_options(seed = TRUE))
  missing_taxIDs <- taxIDs[!(taxIDs %in% unique(taxonomy_df$query_taxID))]
  if (length(missing_taxIDs) > 0) {
    return(extract_and_repeat_search(missing_taxIDs))
  } else {
    return(taxonomy_df)
  }
}

## Executar para extrair e repetir a busca até que todos os TaxIDs sejam encontrados
taxonomy_df <- search_tax(taxid_gen_db$taxid)

taxonomy_df_backup <- taxonomy_df

## Filtrar apenas clados informativos
taxonomy_df <- taxonomy_df %>% 
  filter(!Rank %in% c("no rank","clade")) %>% 
  unique()

## Criar nova tabela de taxonomia wider
taxonomy_tbl_df <- taxonomy_df %>% 
  dplyr::select(-c("TaxId", "query_taxID")) %>%
  unique() %>%
  mutate("genus" = Sci_name ) %>% 
  filter(Rank %in% c("superkingdom","kingdom","phylum","subphylum","class","subclass","order","suborder","family","subfamily","genus")) %>% 
  tidyr::pivot_wider(
    id_cols = c(genus,Sci_name),
    names_from = Rank,
    values_from = c(ScientificName),names_repair = "minimal") %>%
  relocate("Sci_name","genus","superkingdom","kingdom","phylum","subphylum","class","subclass","order","suborder","family","subfamily","genus")

View(taxonomy_tbl_df)

taxonomy_tbl_df %>% colnames() %>% paste0(collapse = "\n") %>% cat()

## fill NA tax with combination of max_tax and rank
for (line in 1:nrow(taxonomy_tbl_df)) {
  # if (taxonomy_tbl_df$genus[line] %in% c("NA",NA,"")) {
  if (taxonomy_tbl_df$superkingdom[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_df$superkingdom[line] <- paste0("superkingdom of ", taxonomy_tbl_df$kingdom[line]) }
  
  if (taxonomy_tbl_df$kingdom[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_df$kingdom[line] <- paste0("kingdom of ", taxonomy_tbl_df$superkingdom[line]) }
  
  if (taxonomy_tbl_df$phylum[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_df$phylum[line] <- paste0("phylum of ", taxonomy_tbl_df$kingdom[line]) }
  
  if (taxonomy_tbl_df$subphylum[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_df$subphylum[line] <- paste0("subphylum of ", taxonomy_tbl_df$phylum[line]) }
  
  if (taxonomy_tbl_df$class[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_df$class[line] <- paste0("class of ", taxonomy_tbl_df$subphylum[line]) }
  
  if (taxonomy_tbl_df$subclass[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_df$subclass[line] <- paste0("subclass of ", taxonomy_tbl_df$class[line]) }
  
  if (taxonomy_tbl_df$order[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_df$order[line] <- paste0("order of ", taxonomy_tbl_df$subclass[line]) }
  
  if (taxonomy_tbl_df$suborder[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_df$suborder[line] <- paste0("suborder of ", taxonomy_tbl_df$order[line]) }
  
  if (taxonomy_tbl_df$family[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_df$family[line] <- paste0("family of ", taxonomy_tbl_df$suborder[line]) }
  
  if (taxonomy_tbl_df$subfamily[line] %in% c("NA",NA,"")) {
    taxonomy_tbl_df$subfamily[line] <- paste0("subfamily of ", taxonomy_tbl_df$family[line]) }
  
  if (is.na(taxonomy_tbl_df$genus[line])) {
    taxonomy_tbl_df$genus[line] <- paste0("genus of ", taxonomy_tbl_df$subfamily[line]) }
  # 
}

View(taxonomy_tbl_df)

# save file
# write.csv(x = taxonomy_tbl_db,
#           file = "/home/igorhan/projetos/LGCdb/data/db_tbl/db_taxonomy.csv", row.names = F)

# join taxonomy with db_tbl
DB_tbl_final <- left_join(x = DB_tbl,
                          y = taxonomy_tbl_df,
                          by = "genus")

View(DB_tbl_final)

# melhorando o output do DB_tbl_final
DB_tbl_final <- DB_tbl_final %>% 
  select(-c(Sci_name, `Original names`)) %>%
  relocate(`Composed name`, `River basin`, Identifier, Names, genus, subfamily, family, 
           suborder, order, subclass, class, subphylum, phylum, kingdom, superkingdom, Cluster)


# 11- Add sequence to db_tbl ----

DB_seqs_df <- dna_all %>%
  as.data.frame() %>%
  rownames_to_column("Composed name") %>%
  rename("Sequence" = x) %>%
  relocate("Composed name", Sequence) %>% 
  tibble()

View(DB_seqs_df)

DB_tbl_final_seqs <- DB_tbl_final %>% 
  left_join(DB_seqs_df, by = "Composed name") %>% 
  relocate(`Composed name`, `Sequence`, `River basin`, Identifier, Names, genus, subfamily, family, 
           suborder, order, subclass, class, subphylum, phylum, kingdom, superkingdom, Cluster)

View(DB_tbl_final_seqs)

# save file
# write.csv(x = DB_tbl_final,
#           file = "/home/igorhan/projetos/LGCdb/data/db_tbl/DB_tbl_final.csv", row.names = F)

# 12- Write DB ordered by cluster to fasta file ----

# No Gaps
dna_order_noGaps <- RemoveGaps(dna_order)
writeXStringSet(x = dna_order_noGaps, filepath = paste0(results_folder,"LGC12Sdb_complete_noGaps-", Sys.Date(), format = ".fasta"))

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
# incluir o código abaixo!!!!

# conferir se contém identificador duplicados e grupos out
names(dna_order) %>% as_tibble() %>% View() 

# caso haja problemas, rodar o código abaixo:

# DB_tbl_final_semConflitos <-
#   DB_tbl_final %>% 
#   select(!c(`Original names`, Sci_name)) %>%
#   relocate(c("Identifier",
#              "Names",
#              "Cluster",
#              "Composed name",
#              "River basin",
#              "genus",
#              "subfamily",
#              "family",
#              "suborder",
#              "order",
#              "subclass",
#              "class",
#              "subphylum",
#              "phylum",
#              "kingdom",
#              "superkingdom")) %>% 
#   filter(!`Identifier` %in% c("0439",
#                               "1382",
#                               "1398",
#                               "1399",
#                               "0905",
#                               "0906",
#                               "JX899737.1",
#                               "KF209618.1",
#                               "LC104399.1")) %>%
#   filter(!`Composed name` %in% "LGC_SF1313 Odontostilbe sp")

# 15- LGC12Sdb 2024 Final ----

View(DB_tbl_final)
