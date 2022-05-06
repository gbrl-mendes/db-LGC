#info ----

# author: Heron O. Hilário
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
#1- import seqs to db ----

#creat/load SQLdb
# db_fish12S <- dbConnect(SQLite(), "/home/heron/prjcts/fish_eDNA/data/refs/db/sql/fish12S_all_seqs.sql")
db_fish12S <- dbConnect(SQLite(), "/home/gabriel/projetos/fish_eDNA/DB/fish12S_all_seqs.sql")

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



#
# #sequências de Prochilodus do NCBI
# Seqs2DB(seqs = "/home/heron/prjcts/proch/data/gb/",
#         type =  "FASTA",dbFile =  db_fish12S, identifier =  "nonBRfish")


BrowseDB(db_fish12S)

#2- add seq info to DB ----

#seqlengths to DB
l <- DECIPHER::IdLengths(dbFile = db_fish12S)

Add2DB(l, db_fish12S, verbose=TRUE)
BrowseDB(db_fish12S)

#3- get all seqs from DB ----

{
dna_sf <- SearchDB(dbFile = db_fish12S,identifier = "São Francisco",nameBy = "description")
dna_jq <- SearchDB(dbFile = db_fish12S,identifier = "Jequitinhonha",nameBy = "description")
dna_non <- SearchDB(dbFile = db_fish12S,identifier = "nonBRfish",nameBy = "description")
dna_sfjq2 <- SearchDB(dbFile = db_fish12S,identifier = "SFJq2",nameBy = "description")
dna_sfjq3 <- SearchDB(dbFile = db_fish12S,identifier = "SFJq3",nameBy = "description")
dna_sfjq4 <- SearchDB(dbFile = db_fish12S,identifier = "SFJq4",nameBy = "description")
dna_sfjq5 <- SearchDB(dbFile = db_fish12S,identifier = "SFJq5",nameBy = "description")

#merge seq sets
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


# identify basin from origin file
          # for (seq in 1:nrow(DB_tbl)) {
          #   if (DB_tbl$`Original names`[seq] %in% names(dna_sf)) {
          #     DB_tbl$`River basin`[seq] = "SF"}
          #
          #   if (DB_tbl$`Original names`[seq] %in% names(dna_jq)) {
          #     DB_tbl$`River basin`[seq] = "Jq"}
          #
          #   if (DB_tbl$`Original names`[seq] %in% names(dna_non)){
          #     DB_tbl$`River basin`[seq] = "Out"}
          #
          #   if (DB_tbl$`Original names`[seq] %in% names(dna_sfjq2)){
          #     DB_tbl$`River basin`[seq] = "SFJq2"}
          # }
for (seq in 1:nrow(DB_tbl)) {
  DB_tbl$`River basin`[seq] <- str_split_fixed(string = DB_tbl$`Original names`[seq],pattern = "-",n = 3)[3]
}




# correct species name
for (seq in 1:nrow(DB_tbl)) {

  DB_tbl$Identifier[seq] <-  str_split_fixed(string = DB_tbl$`Original names`[seq], pattern = "-", n = 2)[1]
  DB_tbl$Names[seq] <- str_split_fixed(string = DB_tbl$`Original names`[seq], pattern = "-", n = 3)[2] %>%
    str_replace_all(pattern = "_",replacement = " ") #make replacement to remov aff, cf...
  DB_tbl$genus[seq] <- str_split(string = DB_tbl$Names[seq],pattern = " ")[[1]][1]

}

# create new name and assign to original sequences

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









#5- rewriting DB with corrected names ----
writeXStringSet(x = dna_all,filepath = "/home/heron/prjcts/fish_eDNA/data/refs/db/LGC/fev22/DB/LGC12Sdb_252seqs-fev22-pretty_names.fasta",format = "fasta")

#6- align DB seqs ----
#remove non-fish sequences

            # dna_all <- dna_all[!(names(dna_all) %in% c("Out | Homo sapiens | JN034109.1", #removing
            #                                            "Out | Bos taurus | AJ885201.1"))] #removing

#orinentarsequencias antes do alinhamento caso tenha alguma invertida






dna_all %>% unique()
dna_all[dna_all %>% duplicated()]
names(dna_all) %>% unique()






names(dna_all)

dna_aligned <- AlignSeqs(myXStringSet = dna_all,refinements = 700,iterations = 700,verbose = TRUE)

# BrowseSeqs(dna_aligned,colorPatterns = c(300,600))
BrowseSeqs(dna_aligned)

#7 - distance matrix ----



# identify and remove
# V05F_898 (5′-AAACTCGTGCCAGCCACC-3′)
# AAACTCGT-GCCAGCC-ACC

#MiFish fwd (GTCGGTAAAACTCGTGCCAGC)

# and teleoR (5′-CTTCCGGTACACTTACCATG-3′)
  # rev comp               # CATGGTAAGTGTACCGGAAG


dna_alg_sub<- subseq(x = dna_aligned,start = 304,end = 1118)   #including MiFish FWD and Teleo REV
dna_alg_sub<- subseq(x = dna_aligned,start = 326,end = 1089)   #excluding MiFish FWD and Teleo REV




# Identify clusters for primer design
dna_alg_sub_dist <- DistanceMatrix(myXStringSet = dna_alg_sub,
                           includeTerminalGaps = FALSE,
                           correction = "Jukes-Cantor",
                           processors = 20)
dim(dna_alg_sub_dist) # a symmetric matrix

#8 - Clusters ----
#identificar clusters na sequencais em função de suas distâncias

dna_clust <- IdClusters(dna_alg_sub_dist,
                         method="UPGMA",
                         cutoff=0.05,
                         show=TRUE, #gera um cladograma, muito dificil de configurar/interpretar
                         type = "both")
head(dna_clust) # cluster numbers


#Adicionar a info de clusters no DB

# Identify sequences by cluster name in the database
Add2DB(data.frame(
  cluster=dna_clust[[1]]$cluster),
  dbFile = db_fish12S,
  verbose = TRUE)

BrowseDB(db_fish12S,orderBy = "cluster")


#add clusters to DB_tbl

for (seq in 1:nrow(DB_tbl)) {
  for (seq2 in 1:nrow(dna_clust[[1]])) {

    if ((rownames(dna_clust[[1]])[seq2]) == (DB_tbl$`Composed name`[seq])) {
      DB_tbl$Cluster[seq] <- (dna_clust[[1]][seq2,])
    }
  }
}

# DB_tbl %>%
#   arrange(by = Cluster)

#order seqs by cluster and write fasta
dna_order <- dna_alg_sub[order(DB_tbl$Cluster),]





writeXStringSet(x = dna_order,
                filepath = "/home/heron/prjcts/fish_eDNA/data/refs/db/LGC/fev22/DB/LGC12Sdb_252seqs-fev22-order_clust.fasta",
                format = "fasta")

# this file was used to generate a MLtree on MEGAX, to identify identic sequences


#9- get taxonomic information to construct DADA2 formated headers----

#load entrez function
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

#clear warnings
#assign("last.warning", NULL, envir = baseenv())


# create and fill taxonomy table for the function
taxonomy_df <- dplyr::tibble()

orgs2search <- unique(sort.default(DB_tbl$genus))

taxonomy_df <- purrr::map_df(orgs2search, extract_taxonomy)
# #taxonomy_df <- purrr::map_df(DB_tbl$genus, extract_taxonomy)


#tentando com um for

taxons_tbl <- tibble::tibble()

for (taxa in 1:length(orgs2search)) {

  # taxonomy_df <- purrr::map_df(orgs2search[taxa], extract_taxonomy) #use only if map_df doesnt work for all together

  tax_df <- taxonomy_df[1]$TaxaSet[[taxa]]$LineageEx$Taxon

  tax_df <- tax_df %>%
    mutate("genus" = taxonomy_df[1]$TaxaSet[[taxa]]$ScientificName)

  taxons_tbl <- dplyr::bind_rows(taxons_tbl,tax_df)

}



# orgs2search

# unique(taxons_tbl$genus)


taxons_tbl_bckp <- taxons_tbl


#o map_df ta dando erro, vou fazer com um for
#
# taxonomy_df <- vector(mode = "list", length = length(orgs2search))
#
# for (org in 1:length(orgs2search)) {
#
#   taxonomy_df[[org]]<- extract_taxonomy(organism_name = orgs2search[org])
#
# }
#





#
#
#
# # create a and fill a dataframe for taxonomic ranks
# orgs_tbl <- tibble::tibble()
#
# for (org in 1:nrow(taxonomy_df)) {
#
#   scin_name <- taxonomy_df$TaxaSet[[org]]$ScientificName
#
#   org_tbl <- tibble::as_tibble(taxonomy_df$TaxaSet[[org]]$LineageEx$Taxon) %>%
#     dplyr::mutate("genus" = scin_name)%>%
#     dplyr::filter(Rank %in% c("kingdom","phylum","class","order","family")) %>%
#     tidyr::pivot_wider(names_from = Rank,values_from = c(ScientificName,TaxId)) %>%
#     dplyr::select(genus,dplyr::starts_with("Scie"))
#
#   orgs_tbl <- dplyr::bind_rows(orgs_tbl,org_tbl)
# }


# create a and fill a dataframe for taxonomic ranks
orgs_tbl <- tibble::tibble()

for (org in 1:length(orgs2search)) {

  org_tbl <- taxons_tbl %>%
    filter(genus == orgs2search[org]) %>%
    dplyr::filter(Rank %in% c("kingdom","phylum","class","order","family")) %>%
    tidyr::pivot_wider(names_from = Rank,values_from = c(ScientificName,TaxId)) %>%
    dplyr::select(genus,dplyr::starts_with("Scie"))

  orgs_tbl <- dplyr::bind_rows(orgs_tbl,org_tbl)
}




#10- bind tax rank cols to DB_tbl ----

colnames(orgs_tbl) <- str_remove(string = colnames(orgs_tbl),
                                 pattern = "ScientificName_")
orgs_tbl <- unique(orgs_tbl)

# DB_tbl_bckp <- DB_tbl
# DB_tbl <- DB_tbl_bckp


DB_tbl <- left_join(x = DB_tbl, y = orgs_tbl,by = "genus")

DB_tbl <- DB_tbl %>% mutate("headers" = as.character(""))

#11- creat fasta headers with ranks for DADA2 ----

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

#12- generate final DB fasta ----



final_dna <- RemoveGaps(myXStringSet = dna_alg_sub,
                        processors = 10,
                        removeGaps = "all")


#rename final dna header
for (seq in 1:nrow(DB_tbl)) {
  for (seq2 in 1:length(final_dna)){

    if (names(final_dna)[seq2] == DB_tbl$`Composed name`[seq]) {

      names(final_dna)[seq2] <- DB_tbl$headers[seq]
    }
  }
}


#order seqs by cluster and write fasta
final_dna_order <- final_dna[order(DB_tbl$Cluster),]



#13 - Identify problematic species to rule out

names(final_dna_order)


#the species to be removed from the final database are:



names(final_dna_order[252])
#193 - Trachelyopterus striatulus 6583
names(final_dna_order[41])
#193 - Orthospinus franciscensis 1071

# final_dna_clean <- final_dna_order[-c(41,252)]
final_dna_clean <- final_dna_order



#write fasta
writeXStringSet(x = final_dna_clean,filepath = "/home/heron/prjcts/fish_eDNA/data/refs/db/LGC/fev22/DB/LGC12Sdb-fev22-dada_tax_fullDB_order.fasta",format = "fasta")
# writeXStringSet(x = final_dna_clean,filepath = "/home/heron/prjcts/fish_eDNA/data/refs/db/LGC/dada_tax_fullDB_order.fasta",format = "fasta")














#create DB seqs for assigSpecies


#rename final dna header
for (seq in 1:nrow(DB_tbl)) {
  for (seq2 in 1:length(final_dna)){

    if (names(final_dna)[seq2] == DB_tbl$headers[seq]) {

      names(final_dna)[seq2] <- paste(DB_tbl$Identifier[seq],DB_tbl$Names[seq])
    }
  }
}


#order seqs by cluster and write fasta
final_dna_order <- final_dna[order(DB_tbl$Cluster),]

# final_dna_clean <- final_dna_order[-c(41,252)]
final_dna_clean <- final_dna_order

#write fasta
# writeXStringSet(x = final_dna_clean,filepath = "/home/heron/prjcts/fish_eDNA/data/refs/db/LGC/dada_tax_fullDB_order_SPs.fasta",format = "fasta")
writeXStringSet(x = final_dna_clean,filepath = "/home/heron/prjcts/fish_eDNA/data/refs/db/LGC/fev22/DB/LGC12Sdb-fev22-dada_SPs_fullDB_order.fasta",format = "fasta")








# rename "aff" and "cf"

# sed -i '/>/ s/aff /aff_/'
# sed -i '/>/ s/cf /cf_/'








#recover Database region for primers

BrowseSeqs(dna_aligned[order(DB_tbl$Cluster),]) #all seqs

#neo FWD: 741-760          CGCCGTCGCAAGCTTACCCT
neo_FWD_seqs <- subseq(dna_aligned[order(DB_tbl$Cluster),],742,760)
neo_FWD_seqs <- c(DNAStringSet(neo_FWD),neo_FWD_seqs)
BrowseSeqs(neo_FWD_seqs) #neo FWD seqs
BrowseSeqs(neo_FWD_seqs,colorPatterns = FALSE,highlight = 1) #neo FWD seqs

#neo REV: 973-992          GC-ACACACCGCCCGTCACT
neo_REV_seqs <- subseq(dna_aligned[order(DB_tbl$Cluster),],972,991)
neo_REV_seqs <- c(DNAStringSet("GC-ACACACCGCCCGTCACT"),neo_REV_seqs)
BrowseSeqs(neo_REV_seqs)
BrowseSeqs(neo_REV_seqs,colorPatterns = FALSE,highlight = 1) #neo REV seqs


#mif FWD: 304-325
mif_FWD_seqs <- subseq(dna_aligned[order(DB_tbl$Cluster),],304,325)
mif_FWD_seqs <- c(DNAStringSet("GTCGGTAAAACTCGT-GCCAGC"),mif_FWD_seqs)
BrowseSeqs(mif_FWD_seqs) #mif FWD seqs
BrowseSeqs(mif_FWD_seqs,colorPatterns = FALSE,highlight = 1) #mif FWD seqs
BrowseSeqs(mif_FWD_seqs,colorPatterns = TRUE,highlight = 1) #mif FWD seqs

#mif REV: 533-559
mif_REV_seqs <- subseq(dna_aligned[order(DB_tbl$Cluster),],532,559)
mif_REV_seqs <- c(DNAStringSet("CAAACTGGGATTAGATACCCCACTATGT"),mif_REV_seqs)
BrowseSeqs(mif_REV_seqs) #mif REV seqs
BrowseSeqs(mif_REV_seqs,colorPatterns = FALSE,highlight = 1) #mif REV seqs
BrowseSeqs(mif_REV_seqs,colorPatterns = TRUE,highlight = 1) #mif REV seqs








BrowseSeqs(neo_FWD_seqs,colorPatterns = FALSE,highlight = 1) #neo seqs
BrowseSeqs(subseq(dna_aligned[order(DB_tbl$Cluster),],761,972)) #neo seqs
BrowseSeqs(subseq(dna_aligned[order(DB_tbl$Cluster),],761,972)) #neo seqs
#mif FWD: 304-325
#mif REV: 533-560
BrowseSeqs(subseq(dna_aligned[order(DB_tbl$Cluster),],326,532)) #mif seqs
#tel FWD: 977-993
#tel REV: 1078-1098
BrowseSeqs(subseq(dna_aligned[order(DB_tbl$Cluster),],994,1077)) #tel seqs





#
#
# how many times each species is in pool
# (jq_tbl$Species %>% str_replace_all(pattern = " ",replacement = "_")) %in%
#   (DB_tbl$`Original names` %>% str_remove(pattern = "^.*.[[:digit:]]_") )



# ADENDUM----
#13 - map primers on DB ----

# set primers sequences
neo_FWD <- "CGCCGTCGCAAGCTTACCCT" #Mini_bar3
neo_REV <- "AGTGACGGGCGGTGTGTGC"  #Mini_bar3

mif_FWD <- "GTCGGTAAAACTCGTGCCAGC"
mif_REV <-    "CATAGTGGGGTATCTAATCCCAGTTTG"
# mif_REV <- "ACATAGTGGGGTATCTAATCCCAGTTTG"

tel_FWD <- "ACACCGCCCGTCACTCT" #
tel_REV <-  "CTTCCGGTACACTTACCATG"
# tel_REV <- "ACTTCCGGTACACTTACCATG"

primers <- DNAStringSet(c(neo_FWD, neo_REV, mif_FWD, mif_REV, tel_FWD, tel_REV))
primers_rev <- reverseComplement(primers)

names(primers) <- c("neo_FWD","neo_REV","mif_FWD","mif_REV","tel_FWD","tel_REV")

##13a -generate primer orientations----

#function to get all possible primer orientations
allOrients <- function(primers) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primers)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}


#gerar todas orientações possíveis dos primers
neo_FWD.orients <- allOrients(neo_FWD)
names(neo_FWD.orients) <- paste0("neo_FWD-", names(neo_FWD.orients))
neo_REV.orients <- allOrients(neo_REV)
names(neo_REV.orients) <- paste0("neo_REV-", names(neo_REV.orients))
mif_FWD.orients <- allOrients(mif_FWD)
names(mif_FWD.orients) <- paste0("mif_FWD-", names(mif_FWD.orients))
mif_REV.orients <- allOrients(mif_REV)
names(mif_REV.orients) <- paste0("mif_REV-", names(mif_REV.orients))
tel_FWD.orients <- allOrients(tel_FWD)
names(tel_FWD.orients) <- paste0("tel_FWD-", names(tel_FWD.orients))
tel_REV.orients <- allOrients(tel_REV)
names(tel_REV.orients) <- paste0("tel_REV-", names(tel_REV.orients))

primers_all <- DNAStringSet(c(neo_FWD.orients, neo_REV.orients,
                              mif_FWD.orients, mif_REV.orients,
                              tel_FWD.orients, tel_REV.orients))

primers_all
as.list(primers_all)

#13b - find primers in aligned DB ----


neo_db_html <- BrowseSeqs(myXStringSet = dna_order,
                          patterns=primers_all, #neo primers
                          colorPatterns = TRUE,
                          colors = c("#00ff00","#99ff33","#00ff00","#99ff33","#00ff00","#99ff33","#00ff00","#99ff33",  # neo verde
                                     "#0000ff","#0066cc","#0000ff","#0066cc","#0000ff","#0066cc","#0000ff","#0066cc",  # mif azul
                                     "#ff0000","#ff3300","#ff0000","#ff3300","#ff0000","#ff3300","#ff0000","#ff3300")) # tel vermelho



#14- isolate neo amplicon ----
                 # changes every time is realigned

neoFWDstart <- 645
neoFWDend <- 664
neoREVstart <- 871
neoREVend <- 890




#identify mismatches in the primer sequence alignment




pool_sps[!(pool_sps %in% DB_tbl$Names)]


DB_tbl$`Composed name`[DB_tbl$Names %in% pool_sps]

BrowseSeqs(dna_order[DB_tbl$`Composed name`[DB_tbl$Names %in% pool_sps]],colorPatterns = FALSE)
















neo_ampli<- subseq(x = dna_aligned,685,896)

BrowseSeqs(neo_ampli)

#rename neo ampli seqs for mega

neo_ampli




neo_dna_order <- subseq(dna_order,367,577)


names(neo_dna_order)

BrowseSeqs(dna_order[126:136],highlight = 0,threshold=0.2,colorPatterns = FALSE)
BrowseSeqs(neo_dna_order[c(184,186)],highlight = 0)


writeXStringSet(x = neo_dna_order,filepath = "/home/heron/prjcts/fish_eDNA/data/refs/db/neo_DB_order_clust.fasta",format = "fasta")









names(neo_ampli)

BrowseSeqs(neo_ampli[c("Jq | Trachelyopterus cf galeatus | 5675","Jq | Trachelyopterus striatulus | 6583")],highlight = 1)














#pegar apenas a sequencia do primer e ordenar

neo_FWD_DNAFISH <- subseq(DNAFISH, neoFWDstart, neoFWDend,) %>% BiocGenerics::sort()
neo_REV_DNAFISH <- subseq(DNAFISH, neoREVstart, neoREVend,) %>% BiocGenerics::sort()


consensus_neo_FWD <- ConsensusSequence(neo_FWD_DNAFISH)
consensus_neo_REV <- ConsensusSequence(neo_REV_DNAFISH)


# NeoFIsh ggseqlogo ----





writeXStringSet(x = neo_ampli,filepath = "/home/heron/prjcts/fish_eDNA/data/refs/db/neo_ampli.fasta",format = "fasta")


names(dna_all)







#BOLD ----





# https://www.boldsystems.org/index.php/datarelease
#https://www.boldsystems.org/data/datarelease/NewPackages/iBOL_phase_6.50_COI.tsv.zip

BOLD_seqs <- read_tsv(file = "/home/heron/prjcts/fish_eDNA/data/refs/db/BOLD/iBOL_phase_6.50_COI.tsv",col_names = TRUE)

BOLD_seqs <- BOLD_seqs %>%
  mutate(Kingdom = "Metazoa",
         Specimen = "NA")

colnames(BOLD_seqs)
# taxLevel = c("Kingdom","Phylum","Class","Order","Family", "Genus", "Species","Specimen","Basin")
BOLD_seqs <- BOLD_seqs %>%  select(nucraw,Kingdom,phylum_reg ,class_reg, order_reg, family_reg, genus_reg, species_reg,Specimen,country_reg)

BOLD_dada <- tibble(Seq = BOLD_seqs$nucraw,
                    Header=paste0(">",BOLD_seqs$Kingdom,";",
                                  BOLD_seqs$phylum_reg, ";",
                                  BOLD_seqs$class_reg,";",
                                  BOLD_seqs$order_reg, ";",
                                  BOLD_seqs$family_reg, ";",
                                  BOLD_seqs$genus_reg, ";",
                                  BOLD_seqs$species_reg,";",
                                  BOLD_seqs$Specimen,";",
                                  BOLD_seqs$country_reg))




# BOLD_dada_fasta <- write.table(BOLD_dada)

BOLD_fasta <- paste0(BOLD_dada$Header,"\n",BOLD_dada$Seq)

write_lines(x = BOLD_fasta,path = "/home/heron/prjcts/fish_eDNA/data/refs/db/BOLD/BOLD_dada_tax.fasta")







BOLD_dada_sp <- tibble(Seq = BOLD_seqs$nucraw,
                    Header=paste0(">",BOLD_seqs$species_reg,";"))




# BOLD_dada_fasta <- write.table(BOLD_dada)

BOLD_fasta_sp <- paste0(BOLD_dada_sp$Header,"\n",BOLD_dada_sp$Seq)

write_lines(x = BOLD_fasta_sp,path = "/home/heron/prjcts/fish_eDNA/data/refs/db/BOLD/BOLD_dada_tax_Sp.fasta")





library(bold)

bold_tax_id(id = "ACG2675")
bold::bold_tax_name("Metazoa")




