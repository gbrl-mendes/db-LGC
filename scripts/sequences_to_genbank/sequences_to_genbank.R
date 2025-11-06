## LGC12sDB to Genbank function ##
## Gabriel A. Mendes
## Nov 2025

sequences_to_genbank <- function(
    input_tbl, 
    metadata_tbl, 
    output_folder = NULL, 
    output_filename = "output.fasta",
    species_filter = NULL, 
    basin_filter = NULL, 
    add_extra_fastas = list(),
    verbose = TRUE,
    which_species = FALSE,
    which_genus = FALSE,
    which_basin = FALSE
) {
  library(Biostrings)
  library(tidyverse)
  library(lubridate)
  library(readxl)
  
  read_table_if_path <- function(tbl) {
    if (is.character(tbl) && length(tbl) == 1) {
      ext <- tools::file_ext(tbl)
      if (tolower(ext) %in% c("csv")) {
        if (verbose) message(paste("Reading CSV file:", tbl))
        return(read_csv(tbl, show_col_types = FALSE))
      } else if (tolower(ext) %in% c("xls", "xlsx")) {
        if (verbose) message(paste("Reading Excel file (suppressing warnings):", tbl))
        # Suppress warnings here so they don’t overflow console
        return(suppressWarnings(read_excel(tbl)))
      } else {
        stop("Unsupported file extension for table input: ", ext)
      }
    }
    return(tbl)
  }
  
  # After reading metadata_tbl, convert all logical columns to character to avoid those 'expecting logical' warnings later:
  logical_cols <- sapply(metadata_tbl, is.logical)
  if (any(logical_cols)) {
    if (verbose) message("Converting logical columns to character in metadata to prevent warnings...")
    metadata_tbl[logical_cols] <- lapply(metadata_tbl[logical_cols], as.character)
  }
  
  input_tbl <- read_table_if_path(input_tbl)
  metadata_tbl <- read_table_if_path(metadata_tbl)
  
  # Convert key metadata columns to character to avoid type mismatch warnings
  metadata_tbl <- metadata_tbl %>%
    mutate(
      `Número LGC` = as.character(`Número LGC`),
      País = as.character(País),
      Bacia = as.character(Bacia)
    )
  
  available_species <- if("species" %in% colnames(input_tbl)) unique(input_tbl$species) else NULL
  available_genus <- if("genus" %in% colnames(input_tbl)) unique(input_tbl$genus) else NULL
  available_basins <- if("River basin" %in% colnames(input_tbl)) unique(input_tbl$`River basin`) else NULL
  
  if (which_species || which_genus || which_basin) {
    res <- list()
    if (which_species) res$species <- sort(available_species)
    if (which_genus) res$genus <- sort(available_genus)
    if (which_basin) res$basin <- sort(available_basins)
    if (verbose) message("Returning available filter options as requested.")
    return(res)
  }
  
  if (verbose) message("Checking input data...")
  if (missing(input_tbl) || is.null(input_tbl) || nrow(input_tbl) == 0) warning("Input sequence table is missing or empty.")
  if (missing(metadata_tbl) || is.null(metadata_tbl) || nrow(metadata_tbl) == 0) warning("Metadata table is missing or empty.")
  
  if (verbose) message("Checking or creating output folder...")
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
    if (verbose) message("Output folder created.")
  } else {
    if (verbose) message("Output folder already exists.")
  }
  
  if (verbose) message("Filtering input sequences...")
  filtered_tbl <- input_tbl
  if (!is.null(species_filter)) {
    filtered_tbl <- filtered_tbl %>% filter(genus %in% species_filter)
    if (verbose) message(paste("Filtered by species:", paste(species_filter, collapse = ", ")))
  }
  if (!is.null(basin_filter)) {
    filtered_tbl <- filtered_tbl %>% filter(`River basin` %in% basin_filter)
    if (verbose) message(paste("Filtered by basin:", paste(basin_filter, collapse = ", ")))
  }
  
  # Custom warning if filtered results are empty—this reliably triggers regardless of whether data was read from file or passed as object
  if (nrow(filtered_tbl) == 0) {
    species_text <- ifelse(is.null(species_filter), "any species/genera", paste0("species/genera filter: '", paste(species_filter, collapse =", "), "'"))
    basin_text <- ifelse(is.null(basin_filter), "any basin", paste0("basin filter: '", paste(basin_filter, collapse = ", "), "'"))
    warning(paste("Warning: No sequences found matching", species_text, "within", basin_text, ".",
                  "Check your filters or use which_species, which_genus, or which_basin options to see available filters."))
    return(invisible(NULL))
  }
  
  if (verbose) message("Adding extra FASTA sequences if provided...")
  if (length(add_extra_fastas) > 0) {
    for (nm in names(add_extra_fastas)) {
      fasta_obj <- add_extra_fastas[[nm]]
      if (length(fasta_obj) == 0) warning(paste0("FASTA sequence object for '", nm, "' is empty."))
      fasta_tbl <- tibble(header = names(fasta_obj),
                          sequence = as.character(fasta_obj)) %>%
        mutate(
          Sequence = str_remove_all(sequence, "-"),
          Identifier = sapply(str_split(header, "_", n = 2), `[`, 1),
          Names = nm,
          `Composed name` = paste0("LGC_DC", Identifier, " ", Names)
        ) %>%
        select(-c(header, sequence))
      filtered_tbl <- bind_rows(filtered_tbl, fasta_tbl)
      if (verbose) message(paste("Added extra fasta sequences for:", nm))
    }
  }
  
  meses <- c(
    jan="01", janeiro="01", fev="02", fevereiro="02", mar="03", marco="03", março="03",
    abr="04", abril="04", mai="05", maio="05", jun="06", junho="06", jul="07", julho="07",
    ago="08", agosto="08", set="09", setembro="09", out="10", outubro="10", nov="11", novembro="11",
    dez="12", dezembro="12"
  )
  
  if (verbose) message("Processing metadata and formatting dates...")
  md <- suppressWarnings({
    metadata_tbl %>%
      mutate(
        Identifier = str_pad(`Número LGC`, 4, "left", "0"),
        País = case_when(País == "Brasil" ~ "Brazil", TRUE ~ País)
      ) %>%
      mutate(
        data_simplificada = case_when(
          is.na(`Data da coleta`) ~ NA_character_,
          str_detect(`Data da coleta`, "^\\d{4}$") ~ `Data da coleta`,
          str_detect(`Data da coleta`, "^\\d{4}-\\d{2}-00\\b") ~ paste0(str_sub(`Data da coleta`, 1, 4), "/", str_sub(`Data da coleta`, 6, 7)),
          str_detect(`Data da coleta`, "^\\d{5}$") ~ format(as.Date(as.numeric(`Data da coleta`), origin = "1899-12-30"), "%Y/%m"),
          str_detect(`Data da coleta`, "^\\d{2}/\\d{2}/\\d{4}$") ~ format(lubridate::dmy(`Data da coleta`), "%Y/%m"),
          str_detect(`Data da coleta`, "^\\d{4}-\\d{2}-\\d{2}$") ~ format(lubridate::ymd(`Data da coleta`), "%Y/%m"),
          str_detect(`Data da coleta`, "^\\d{1,2} de \\w+ de \\d{4}$") ~ {
            mes <- str_to_lower(str_extract(`Data da coleta`, "(?<= de )\\w+(?= de )"))
            ano <- str_extract(`Data da coleta`, "\\d{4}$")
            paste0(ano, "/", meses[mes])
          },
          str_detect(`Data da coleta`, "^(\\d{1,2})?(jan|fev|mar|abr|mai|maio|jun|jul|ago|set|out|nov|dez|janeiro|fevereiro|marco|março|abril|junho|julho|agosto|setembro|outubro|novembro|dezembro)\\d{4}$") ~ {
            ano <- str_extract(`Data da coleta`, "\\d{4}$")
            mes <- str_extract(`Data da coleta`, "jan|fev|mar|abr|mai|maio|jun|jul|ago|set|out|nov|dez|janeiro|fevereiro|marco|março|abril|junho|julho|agosto|setembro|outubro|novembro|dezembro")
            paste0(ano, "/", meses[mes])
          },
          str_detect(`Data da coleta`, "^\\d{2}/[a-z]{3}/\\d{4}$") ~ {
            ano <- str_extract(`Data da coleta`, "\\d{4}$")
            mes <- str_extract(`Data da coleta`, "(?<=/)[a-z]{3}(?=/|\\d)")
            paste0(ano, "/", meses[mes])
          },
          str_detect(`Data da coleta`, "^\\d{4}/\\d{2}$") ~ `Data da coleta`,
          str_detect(`Data da coleta`, "^\\d{4}-\\d{4}$") ~ str_sub(`Data da coleta`, 6, 9),
          str_detect(`Data da coleta`, "^\\d{2}\\.\\d{2}\\.\\d{4}$") ~ format(lubridate::dmy(str_replace_all(`Data da coleta`, "\\.", "/")), "%Y/%m"),
          str_detect(`Data da coleta`, "semestre \\d{4}$") ~ str_extract(`Data da coleta`, "\\d{4}$"),
          str_detect(`Data da coleta`, "\\d{2}([-\\./])\\d{2}\\1\\d{4}$") ~ format(lubridate::dmy(str_replace_all(str_extract(`Data da coleta`, "\\d{2}([-\\./])\\d{2}\\1\\d{4}$"), "[-\\.]", "/")), "%Y/%m"),
          str_detect(`Data da coleta`, "\\d{2}([-\\./])\\d{4}$") ~ {
            final <- str_extract(`Data da coleta`, "\\d{2}([-\\./])\\d{4}$")
            parts <- unlist(str_split(str_replace_all(final, "[-\\.]", "/"), "/"))
            paste0(parts[2], "/", parts[1])
          },
          str_detect(`Data da coleta`, "^\\d{4}-\\d{5}$") ~ {
            y <- str_extract(`Data da coleta`, "^\\d{4}")
            m <- str_sub(str_extract(`Data da coleta`, "-\\d{5}$"), 2, 3)
            paste0(y, "/", m)
          },
          TRUE ~ NA_character_
        )
      )
  })
  
  if (verbose) message("Joining sequences with metadata and formatting FASTA headers...")
  final_tbl <- filtered_tbl %>%
    rename("organism" = Names) %>%
    left_join(md, by = "Identifier") %>%
    mutate(
      SeqID = str_extract(`Composed name`, "^[^ ]+"),
      isolate = Identifier,
      isolation_source = case_when(
        Bacia == "São Francisco" ~ "Sao Francisco River Basin",
        Bacia == "Rio Doce" ~ "Doce River Basin",
        Bacia == "Jequitinhonha" ~ "Jequitinhonha River Basin",
        Bacia == "Pardo" ~ "Pardo River Basin",
        TRUE ~ "BASIN MISSING"
      ),
      header = paste0(">", SeqID,
                      " [organism=", organism, "]",
                      " [isolate=", isolate, "]",
                      " [isolation_source=", isolation_source, "]",
                      " [country=", País, "]",
                      " [collection_date=", data_simplificada, "]")
    ) %>%
    filter(!is.na(data_simplificada)) %>%
    arrange(header)
  
  if (verbose) message("Writing sequences to FASTA file: ", output_filename)
  fasta_out <- c(rbind(final_tbl$header, final_tbl$Sequence))
  write(fasta_out, file = file.path(output_folder, output_filename))
  
  if (verbose) message("FASTA export completed successfully.")
  
  return(invisible(fasta_out))
}
