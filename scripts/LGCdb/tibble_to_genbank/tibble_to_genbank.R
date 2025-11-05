---
title: "Formatting sequences from tibble to .fasta for Genbank"
author: "Gabriel Mendes"
date: "07/2025"
---
  
# 1- Libraries needed ----
{
  library(tidyverse)
  library(dplyr)
  library(stringr)
  library(readxl)
}

# 2- Define output path ----

output_folder <- "/home/gabriel/projetos/db-LGC/scripts/LGCdb/tibble_to_genbank/output"

# Does this directory exists?
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
  message("Directory created: ", output_folder)
} else {
  message("The directory already exists: ", output_folder)
}

# 3- Define input ----

input_tbl <- DB_tbl_final_seqs # objeto final do lgc_db_full--mai25.R, input deste trabalho
metadata_tbl <- read_excel("/home/igorhan/projetos/LGCdb/data/db_tbl/Banco de Tecidos LGC 2023.xlsx")

# 4- Workaround on metadata ----

# Filtrando as bacias e correcoes menores
metadata_tbl_less <- metadata_tbl %>% 
  # filter(Bacia %in% c("Rio Doce",
  #                     "Jequitinhonha",
  #                     "São Francisco",
  #                     "Pardo")) %>%
  filter(Gênero == "Trichogenes") %>% 
  mutate(Identifier = str_pad(`Número LGC`, 4, "left", 0)) %>% 
  mutate(País = case_when(País == "Brasil" ~ "Brazil",
                          TRUE ~ País)) %>%
  select(#`Status da amostra`,
         `Número LGC`,
         Identifier,
         `Espécie 2023`,
         País,
         Estado,
         Bacia,
         # `Longitude (em graus decimais)`,
         # `Latitude (em graus decimais)`,
         # `Coletor (es)`,
         `Data da coleta`)
  
View(metadata_tbl_less)

# Padronizando os multiplos formatos de datas
{
  # Vetor de meses ao início para reusar nas regras relacionadas
  meses <- c(
    jan="01", janeiro="01", fev="02", fevereiro="02", mar="03", marco="03", março="03",
    abr="04", abril="04", mai="05", maio="05", jun="06", junho="06", jul="07", julho="07",
    ago="08", agosto="08", set="09", setembro="09", out="10", outubro="10", nov="11", novembro="11",
    dez="12", dezembro="12"
  )
  # Pipeline para padronizar os formatos de data
  metadata_tbl_final_less <- metadata_tbl_less %>% 
    mutate(
      data_simplificada = case_when(
        is.na(`Data da coleta`) ~ NA_character_,
        str_detect(`Data da coleta`, "^\\d{4}$") ~ `Data da coleta`,
        str_detect(`Data da coleta`, "^\\d{4}-\\d{2}-00\\b") ~ paste0(str_sub(`Data da coleta`, 1, 4), "/", str_sub(`Data da coleta`, 6, 7)),
        str_detect(`Data da coleta`, "^\\d{5}$") ~ format(as.Date(as.numeric(`Data da coleta`), origin = "1899-12-30"), "%Y/%m"),
        str_detect(`Data da coleta`, "^\\d{2}/\\d{2}/\\d{4}$") ~ format(dmy(`Data da coleta`), "%Y/%m"),
        str_detect(`Data da coleta`, "^\\d{4}-\\d{2}-\\d{2}$") ~ format(ymd(`Data da coleta`), "%Y/%m"),
        
        # Texto por extenso: "7 de junho de 2015"
        str_detect(`Data da coleta`, "^\\d{1,2} de \\w+ de \\d{4}$") ~ {
          mes <- str_to_lower(str_extract(`Data da coleta`, "(?<= de )\\w+(?= de )"))
          ano <- str_extract(`Data da coleta`, "\\d{4}$")
          paste0(ano, "/", meses[mes])
        },
        
        # Formatos "25maio2009", "dez2008", "agosto2008"
        str_detect(`Data da coleta`, "^(\\d{1,2})?(jan|fev|mar|abr|mai|maio|jun|jul|ago|set|out|nov|dez|janeiro|fevereiro|marco|março|abril|junho|julho|agosto|setembro|outubro|novembro|dezembro)\\d{4}$") ~ {
          ano <- str_extract(`Data da coleta`, "\\d{4}$")
          mes <- str_extract(`Data da coleta`, "jan|fev|mar|abr|mai|maio|jun|jul|ago|set|out|nov|dez|janeiro|fevereiro|marco|março|abril|junho|julho|agosto|setembro|outubro|novembro|dezembro")
          paste0(ano, "/", meses[mes])
        },
        
        # "24/set/2009"
        str_detect(`Data da coleta`, "^\\d{2}/[a-z]{3}/\\d{4}$") ~ {
          ano <- str_extract(`Data da coleta`, "\\d{4}$")
          mes <- str_extract(`Data da coleta`, "(?<=/)[a-z]{3}(?=/|\\d)")
          paste0(ano, "/", meses[mes])
        },
        
        # Ano/mês invertido: "2016/01"
        str_detect(`Data da coleta`, "^\\d{4}/\\d{2}$") ~ `Data da coleta`,
        
        # Período de ano: "2005-2006"
        str_detect(`Data da coleta`, "^\\d{4}-\\d{4}$") ~ str_sub(`Data da coleta`, 6, 9),
        
        # Data com ponto: "24.05.2014"
        str_detect(`Data da coleta`, "^\\d{2}\\.\\d{2}\\.\\d{4}$") ~ format(dmy(str_replace_all(`Data da coleta`, "\\.", "/")), "%Y/%m"),
        
        # "1 semestre 2009"
        str_detect(`Data da coleta`, "semestre \\d{4}$") ~ str_extract(`Data da coleta`, "\\d{4}$"),
        
        # Finais DD/MM/AAAA ou DD-MM-AAAA
        str_detect(`Data da coleta`, "\\d{2}([-\\./])\\d{2}\\1\\d{4}$") ~ format(dmy(str_replace_all(str_extract(`Data da coleta`, "\\d{2}([-\\./])\\d{2}\\1\\d{4}$"), "[-\\.]","/")), "%Y/%m"),
        
        # Finais MM/AAAA ou MM-AAAA
        str_detect(`Data da coleta`, "\\d{2}([-\\./])\\d{4}$") ~ {
          final <- str_extract(`Data da coleta`, "\\d{2}([-\\./])\\d{4}$")
          parts <- unlist(str_split(str_replace_all(final, "[-\\.]", "/"), "/"))
          paste0(parts[2], "/", parts[1])
        },
        
        # "2012-07013" (ano-tracinho-cinco dígitos começando pelo mês)
        str_detect(`Data da coleta`, "^\\d{4}-\\d{5}$") ~ {
          y <- str_extract(`Data da coleta`, "^\\d{4}")
          m <- str_sub(str_extract(`Data da coleta`, "-\\d{5}$"), 2, 3)
          paste0(y, "/", m)
        },
        
        TRUE ~ NA_character_
      )
    )
}

View(metadata_tbl_final_less)

metadata_tbl_final_less %>% 
  select(data_simplificada, `Data da coleta`) %>% 
  filter(is.na(data_simplificada)) %>%  
  distinct()

# 5- Workaround on input ----

#Resizing input_tbl
View(input_tbl)

colnames(input_tbl)

cat(paste0('"', colnames(input_tbl), '"', collapse = ","))

input_tbl_less <- input_tbl %>% 
  filter(`River basin` %in% c("SF", 
                              "DC", 
                              "JQ", 
                              "PD")) %>% 
  select("Composed name",
         "Sequence",
         "River basin",
         "Identifier",
         "Names")

View(input_tbl_less)  

#Adding metadata
input_tbl_format <- 
  input_tbl_less %>%
  # input_tbl_less[1,] %>%
  mutate(SeqID = str_extract(`Composed name`, "^[^ ]+")) %>% 
  rename("organism" = Names) %>%
  mutate("isolate" = Identifier) %>% 
  mutate(isolation_source = case_when(
    `River basin` == "SF" ~ "Sao Francisco River Basin",
    `River basin` == "DC" ~ "Doce River Basin",
    `River basin` == "JQ" ~ "Jequitinhonha River Basin",
    TRUE ~ "Pardo River Basin"
  )) %>%
  left_join(metadata_tbl_final_less,
            by = "Identifier") %>% 
  mutate(header = paste0(
    ">", SeqID,
    " [", "organism", "=", organism, "]",
    " [", "isolate", "=", isolate, "]",
    " [", "isolation_source", "=", isolation_source, "]",
    " [", "country", "=", País, "]",
    " [", "collection_date", "=", data_simplificada, "]"
  # )) 
  )) %>%
  select(`Número LGC`,
         Identifier,
         organism,
         `Espécie 2023`,
         `Data da coleta`,
         data_simplificada,
         SeqID,
         isolate,
         isolation_source,
         header,
         Sequence)

View(input_tbl_format)

#Sequencias sem collection_date
no_date <- 
  input_tbl_format %>%
  select(`Número LGC`, `Espécie 2023`, Identifier, organism, data_simplificada, `Data da coleta`) %>% 
  filter(is.na(data_simplificada)) 

View(no_date)
View(input_tbl_less)
View(input_tbl)

# 6- Creating the LGC12sDB's .fasta ----

#Filtering seqs without date
input_tbl_format_filt <- input_tbl_format %>% 
  filter(!is.na(data_simplificada)) %>% 
  arrange(header)
  # head() #pegando apenas as 6 primeiras para o teste

View(input_tbl_format_filt)

#Creating test .fasta for the test submission
fasta_output <- c(rbind(input_tbl_format_filt$header, input_tbl_format_filt$Sequence)) 
# write(fasta_output, paste0(output_folder, "/LGC12sDB_genbank.fasta")) 
write(fasta_output, paste0(output_folder, "/LGC12sDB_genbank_test.fasta"))

#Reading the LGC12sDB's .fasta 
DB_fasta <- readDNAStringSet(filepath = paste0(output_folder, 
                                      # "/LGC12sDB_genbank.fasta"))
                                      "/LGC12sDB_genbank_test.fasta"))

BrowseSeqs(DB_fasta)
  
  
  
  
  
  




