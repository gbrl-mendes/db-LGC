#info ----

# author: Heron O. Hilário; Gabriel M. Brito
# purpose: integrate sequence files from different origins
#         to construct the LGC 12S sequence database and format it
#         for usage in DADA2.
#
# references:
# https://www.bioconductor.org/help/course-materials/2016/BioC2016/ConcurrentWorkshops2/Wright/BigBioSeqData.pdf
# https://www.bioconductor.org/packages/release/bioc/vignettes/DECIPHER/inst/doc/DECIPHERing.pdf
# http://www2.decipher.codes/RLessons/RLesson14.html
# https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12687
# http://www2.decipher.codes/OligoDesign.html
# http://www2.decipher.codes/Bioinformatics/BigBioSeqData/BigBioSeqData2.html
# https://cran.r-project.org/web/packages/ggdendro/vignettes/ggdendro.html

#0- load libraries ----
{
  library(dada2)
  library(Rcpp)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(spider)
  library(ape)
  library(RColorBrewer)
  library(ggplot2)
  library(ggseqlogo)
  library(ips)
  library(ggrepel)
  library(viridis)
  library(phyloseq)
  library(Biostrings)
  library(ShortRead)
  library(pheatmap)
  library(DECIPHER)
  # library(jpeg)
  # library(ggdendro)
  library(readr)
}

#1- Import seqs to db ----

## Creat/load SQLdb
db_fish12S <- dbConnect(SQLite(), "/home/gabriel/projetos/db-LGC/database/fish12S_all_seqs.sql")

## Loading seqs to db
{
#sequências do São Francisco
Seqs2DB(seqs  = "/home/heron/prjcts/fish_eDNA/DB/mai22/refs/12S_full_unlign.fas",
        type =  "FASTA",dbFile =  db_fish12S, identifier =  "São Francisco")
#sequências do Jequi
Seqs2DB(seqs = "/home/heron/prjcts/fish_eDNA/DB/mai22/refs/alinhamentojequi57seq_ed_unalig.fasta",
        type =  "FASTA",dbFile =  db_fish12S, identifier =  "Jequitinhonha")
#sequências de referência
Seqs2DB(seqs = "/home/heron/prjcts/fish_eDNA/DB/mai22/refs/non_BRfish.fasta",
        type =  "FASTA",dbFile =  db_fish12S, identifier =  "nonBRfish")
#sequências de outras espécies para completar o DB
Seqs2DB(seqs = "/home/heron/prjcts/fish_eDNA/DB/mai22/refs/new_Jq_species.fasta",
        type =  "FASTA",dbFile =  db_fish12S, identifier =  "SFJq2")

#sequências Igor HAN nov21
Seqs2DB(seqs = "/home/heron/prjcts/fish_eDNA/DB/mai22/refs/12Sdb_Seq3a_03nov21-10sp_curado.fas",
        type =  "FASTA",dbFile =  db_fish12S, identifier =  "SFJq3")
#sequências Igor HAN dez21
Seqs2DB(seqs = "/home/heron/prjcts/fish_eDNA/DB/mai22/refs/12Sdb_Seq4b_24dez21-13sp_curado.fas",
        type =  "FASTA",dbFile =  db_fish12S, identifier =  "SFJq4")
#sequências Igor HAN jan22
Seqs2DB(seqs = "/home/heron/prjcts/fish_eDNA/DB/mai22/refs/12Sdb_Seq5c_17jan22-9sp_curado.fas",
        type =  "FASTA",dbFile =  db_fish12S, identifier =  "SFJq5")
#sequências Júnio mai22 (3 engraulidae)
Seqs2DB(seqs = "/home/heron/prjcts/fish_eDNA/DB/mai22/refs/3_engraulidae.fasta",
        type =  "FASTA",dbFile =  db_fish12S, identifier =  "SFJq6")

}

BrowseDB(db_fish12S)

#2- Add seq info to DB ----

## Seqlengths to DB
l <- DECIPHER::IdLengths(dbFile = db_fish12S)

Add2DB(l, db_fish12S, verbose=TRUE)
BrowseDB(db_fish12S)

#3- Get all seqs from DB ----

{
dna_sf <- SearchDB(dbFile = db_fish12S,identifier = "São Francisco",nameBy = "description")
dna_jq <- SearchDB(dbFile = db_fish12S,identifier = "Jequitinhonha",nameBy = "description")
dna_non <- SearchDB(dbFile = db_fish12S,identifier = "nonBRfish",nameBy = "description")
dna_sfjq2 <- SearchDB(dbFile = db_fish12S,identifier = "SFJq2",nameBy = "description")
dna_sfjq3 <- SearchDB(dbFile = db_fish12S,identifier = "SFJq3",nameBy = "description")
dna_sfjq4 <- SearchDB(dbFile = db_fish12S,identifier = "SFJq4",nameBy = "description")
dna_sfjq5 <- SearchDB(dbFile = db_fish12S,identifier = "SFJq5",nameBy = "description")

## Merge seq sets
dna_all <- c(dna_sf,dna_jq, dna_non,dna_sfjq2,dna_sfjq3,dna_sfjq4,dna_sfjq5)
}

#4- create and assign new names ----

DB_tbl <- tibble("Original names" = names(dna_all),
                 "Identifier" = as.character(""),
                 "Names" = as.character(""),
                 "River basin" = as.character(""),
                 "Cluster" = as.numeric(""),
                 "Composed name"= as.character(""),
                 "genus"= as.character("")) # it will be used to search NCBI

for (seq in 1:nrow(DB_tbl)) {
  DB_tbl$`River basin`[seq] <- str_split_fixed(string = DB_tbl$`Original names`[seq],pattern = "-",n = 3)[3]
}

## Correct species name
for (seq in 1:nrow(DB_tbl)) {

  DB_tbl$Identifier[seq] <-  str_split_fixed(string = DB_tbl$`Original names`[seq], pattern = "-", n = 2)[1]
  DB_tbl$Names[seq] <- str_split_fixed(string = DB_tbl$`Original names`[seq], pattern = "-", n = 3)[2] %>%
    str_replace_all(pattern = "_",replacement = " ") #make replacement to remov aff, cf...
  DB_tbl$genus[seq] <- str_split(string = DB_tbl$Names[seq],pattern = " ")[[1]][1]

}

## Create new name and assign to original sequences

for (seq in 1:nrow(DB_tbl)) {
  for (seq2 in 1:length(dna_all)){

    if (names(dna_all)[seq2] == DB_tbl$`Original names`[seq]) {

      names(dna_all)[seq2] <- paste0(DB_tbl$`River basin`[seq]," | ",DB_tbl$Names[seq]," | ",DB_tbl$Identifier[seq])
      DB_tbl$`Composed name`[seq] <- names(dna_all)[seq2]
    }
  }
}

dna_all <- RemoveGaps(myXStringSet = dna_all,
                      processors = 10,
                      removeGaps = "all")

dna_all <- OrientNucleotides(dna_all)

#5- Rewriting DB with corrected names ----

writeXStringSet(x = dna_all,filepath = "/home/gabriel/projetos/db-LGC/data/mai22/LGC12Sdb_252seqs-mai22-pretty_names.fasta",format = "fasta")

#6- Align DB seqs ----

dna_all %>% unique()
dna_all[dna_all %>% duplicated()]
names(dna_all) %>% unique()
names(dna_all)

dna_aligned <- AlignSeqs(myXStringSet = dna_all,refinements = 700,iterations = 700,verbose = TRUE)

BrowseSeqs(dna_aligned)

#7 - Distance matrix ----
{
  dna_alg_sub<- subseq(x = dna_aligned,start = 304,end = 1118)   #including MiFish FWD and Teleo REV
dna_alg_sub<- subseq(x = dna_aligned,start = 326,end = 1089)   #excluding MiFish FWD and Teleo REV
}

## Identify clusters for primer design
dna_alg_sub_dist <- DistanceMatrix(myXStringSet = dna_alg_sub,
                           includeTerminalGaps = FALSE,
                           correction = "Jukes-Cantor",
                           processors = 20)
dim(dna_alg_sub_dist) # a symmetric matrix

#8 - Clusters ----

## Identificar clusters na sequencais em função de suas distâncias

dna_clust <- IdClusters(dna_alg_sub_dist,
                         method="UPGMA",
                         cutoff=0.05,
                         show=TRUE, #gera um cladograma, muito dificil de configurar/interpretar
                         type = "both")
head(dna_clust) # cluster numbers

## Identify sequences by cluster name in the database
Add2DB(data.frame(
  cluster=dna_clust[[1]]$cluster),
  dbFile = db_fish12S,
  verbose = TRUE)

BrowseDB(db_fish12S,orderBy = "cluster")


## Add clusters to DB_tbl

for (seq in 1:nrow(DB_tbl)) {
  for (seq2 in 1:nrow(dna_clust[[1]])) {

    if ((rownames(dna_clust[[1]])[seq2]) == (DB_tbl$`Composed name`[seq])) {
      DB_tbl$Cluster[seq] <- (dna_clust[[1]][seq2,])
    }
  }
}

dna_order <- dna_alg_sub[order(DB_tbl$Cluster),]

writeXStringSet(x = dna_order,
                filepath = "/home/gabriel/projetos/db-LGC/data/mai22/LGC12Sdb_252seqs-mai22-order_clust.fasta",
                format = "fasta")

## This file was used to generate a MLtree on MEGAX, to identify identic sequences


#9 - Get taxonomic information to construct DADA2 formated headers----

## load entrez function
extract_taxonomy <- function(organism_name) {
  `%>%` <- dplyr::`%>%`
  organism_xml <- sys::exec_internal("/usr/bin/bash",
                                     args = c("/home/heron/prjcts/fish_eDNA/src/entrez_taxonomy.sh", organism_name),
                                     timeout = 1000)
  if (length(sys::as_text(organism_xml$stdout)) ==  0) {
    message(paste0("deu errado para: ",organism_name))
    return(tibble())
  }else{
    organism_df <- sys::as_text(organism_xml$stdout) %>%
      jsonlite::fromJSON()
    message(paste0("deu certo para: ",organism_name))

    return(organism_df)
  }
}

## Create and fill taxonomy table for the function
taxonomy_df <- dplyr::tibble()

orgs2search <- unique(sort.default(DB_tbl$genus))

taxonomy_df <- purrr::map_df(orgs2search, extract_taxonomy)

## Tentando com um for

taxons_tbl <- tibble::tibble()

for (taxa in 1:length(orgs2search)) {

  tax_df <- taxonomy_df[1]$TaxaSet[[taxa]]$LineageEx$Taxon

  tax_df <- tax_df %>%
    mutate("genus" = taxonomy_df[1]$TaxaSet[[taxa]]$ScientificName)

  taxons_tbl <- dplyr::bind_rows(taxons_tbl,tax_df)

}

## Orgs2search

taxons_tbl_bckp <- taxons_tbl

## Create a and fill a dataframe for taxonomic ranks
orgs_tbl <- tibble::tibble()

for (org in 1:length(orgs2search)) {

  org_tbl <- taxons_tbl %>%
    filter(genus == orgs2search[org]) %>%
    dplyr::filter(Rank %in% c("kingdom","phylum","class","order","family")) %>%
    tidyr::pivot_wider(names_from = Rank,values_from = c(ScientificName,TaxId)) %>%
    dplyr::select(genus,dplyr::starts_with("Scie"))

  orgs_tbl <- dplyr::bind_rows(orgs_tbl,org_tbl)
}

#10 - Bind tax rank cols to DB_tbl ----

colnames(orgs_tbl) <- str_remove(string = colnames(orgs_tbl),
                                 pattern = "ScientificName_")
orgs_tbl <- unique(orgs_tbl)

DB_tbl <- left_join(x = DB_tbl, y = orgs_tbl,by = "genus")
DB_tbl <- DB_tbl %>% mutate("headers" = as.character(""))

#11- Creat fasta headers with ranks for DADA2 ----

for (seq in 1:nrow(DB_tbl)) {
  DB_tbl$headers[seq] <- paste0(DB_tbl$kingdom[seq],";",
                      DB_tbl$phylum[seq],";",
                      DB_tbl$class[seq],";",
                      DB_tbl$order[seq],";",
                      DB_tbl$family[seq],";",
                      DB_tbl$genus[seq],";",
                      DB_tbl$Names[seq],";",
                      DB_tbl$Identifier[seq],";",
                      DB_tbl$`River basin`[seq])

}

#12 - Generate final DB fasta ----

final_dna <- RemoveGaps(myXStringSet = dna_alg_sub,
                        processors = 10,
                        removeGaps = "all")

## Rename final dna header
for (seq in 1:nrow(DB_tbl)) {
  for (seq2 in 1:length(final_dna)){

    if (names(final_dna)[seq2] == DB_tbl$`Composed name`[seq]) {

      names(final_dna)[seq2] <- DB_tbl$headers[seq]
    }
  }
}


## Order seqs by cluster and write fasta
final_dna_order <- final_dna[order(DB_tbl$Cluster),]



## Identify problematic species to rule out

names(final_dna_order)

## The species to be removed from the final database are:

names(final_dna_order[252])
#193 - Trachelyopterus striatulus 6583

names(final_dna_order[41])
#193 - Orthospinus franciscensis 1071

final_dna_clean <- final_dna_order

## Write fasta
writeXStringSet(x = final_dna_clean,filepath = "/home/gabriel/projetos/db-LGC/data/mai22/LGC12Sdb-mai22-dada_tax_fullDB_order.fasta",format = "fasta")

## Create DB seqs for assigSpecies
for (seq in 1:nrow(DB_tbl)) {
  for (seq2 in 1:length(final_dna)){

    if (names(final_dna)[seq2] == DB_tbl$headers[seq]) {

      names(final_dna)[seq2] <- paste(DB_tbl$Identifier[seq],DB_tbl$Names[seq])
    }
  }
}

## Order seqs by cluster and write fasta
final_dna_order <- final_dna[order(DB_tbl$Cluster),]
final_dna_clean <- final_dna_order

## Write fasta
writeXStringSet(x = final_dna_clean,filepath = "/home/gabriel/projetos/db-LGC/data/mai22/LGC12Sdb-mai22-dada_SPs_fullDB_order.fasta",format = "fasta")

## Recover Database region for primers
BrowseSeqs(dna_aligned[order(DB_tbl$Cluster),]) #all seqs

#neo FWD: 741-760          CGCCGTCGCAAGCTTACCCT
{
  neo_FWD_seqs <- subseq(dna_aligned[order(DB_tbl$Cluster),],742,760)
  neo_FWD_seqs <- c(DNAStringSet(neo_FWD),neo_FWD_seqs)
  BrowseSeqs(neo_FWD_seqs) #neo FWD seqs
  BrowseSeqs(neo_FWD_seqs,colorPatterns = FALSE,highlight = 1) #neo FWD seqs
}

#neo REV: 973-992          GC-ACACACCGCCCGTCACT
{
  neo_REV_seqs <- subseq(dna_aligned[order(DB_tbl$Cluster),],972,991)
  neo_REV_seqs <- c(DNAStringSet("GC-ACACACCGCCCGTCACT"),neo_REV_seqs)
  BrowseSeqs(neo_REV_seqs)
  BrowseSeqs(neo_REV_seqs,colorPatterns = FALSE,highlight = 1) #neo REV seqs
}

#mif FWD: 304-325
{
  mif_FWD_seqs <- subseq(dna_aligned[order(DB_tbl$Cluster),],304,325)
  mif_FWD_seqs <- c(DNAStringSet("GTCGGTAAAACTCGT-GCCAGC"),mif_FWD_seqs)
  BrowseSeqs(mif_FWD_seqs) #mif FWD seqs
  BrowseSeqs(mif_FWD_seqs,colorPatterns = FALSE,highlight = 1) #mif FWD seqs
  BrowseSeqs(mif_FWD_seqs,colorPatterns = TRUE,highlight = 1) #mif FWD seqs
}

#mif REV: 533-559
{
  mif_REV_seqs <- subseq(dna_aligned[order(DB_tbl$Cluster),],532,559)
  mif_REV_seqs <- c(DNAStringSet("CAAACTGGGATTAGATACCCCACTATGT"),mif_REV_seqs)
  BrowseSeqs(mif_REV_seqs) #mif REV seqs
  BrowseSeqs(mif_REV_seqs,colorPatterns = FALSE,highlight = 1) #mif REV seqs
  BrowseSeqs(mif_REV_seqs,colorPatterns = TRUE,highlight = 1) #mif REV seqs
}

BrowseSeqs(neo_FWD_seqs,colorPatterns = FALSE,highlight = 1) #neo seqs
BrowseSeqs(subseq(dna_aligned[order(DB_tbl$Cluster),],761,972)) #neo seqs
BrowseSeqs(subseq(dna_aligned[order(DB_tbl$Cluster),],761,972)) #neo seqs
#mif FWD: 304-325
#mif REV: 533-560
BrowseSeqs(subseq(dna_aligned[order(DB_tbl$Cluster),],326,532)) #mif seqs
#tel FWD: 977-993
#tel REV: 1078-1098
BrowseSeqs(subseq(dna_aligned[order(DB_tbl$Cluster),],994,1077)) #tel seqs