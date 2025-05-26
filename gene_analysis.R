rm(list=ls())

library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggpubr)
library(ggplot2)
library(writexl)
library(stringr)
library(dplyr)
library(ComplexHeatmap)

config <- list(
  dirs = list(
    base = "D:/CNVkit/tumor/PTJ_WES_IDT-30802789",
    out = "D:/CNVkit/tumor/tumor_imgs",
    res = "D:/CNVkit/tumor/SBALLATE_tumor_res"
  ), 
  dir_pattern = ".*-t.*$"
)

config <- list(
  dirs = list(
    base = "D:/CNVkit/model/WES_modelli",
    out = "D:/CNVkit/model/with_pooledref/model_imgs",
    res = "D:/CNVkit/model/with_pooledref/model_res"
  ),
  dir_pattern = ".*tum-001.*$"
)

# =====================================================================
# Function to import CNVkit files
# =====================================================================
import_cnvkit_files <- function(file_type, cols) {
  safely_read <- safely(read_tsv)
  
  samples <- list.files(config$dirs$res, pattern = config$dir_pattern) %>% 
    setdiff(c("multi_sample", "nfr", "nfrmulti_sample"))
  
  map_dfr(samples, ~{
    file_path <- file.path(config$dirs$res, .x, paste0(.x, "_", file_type))
    message("Importing sample: ", file_path)
    if(!file.exists(file_path)) {
      message("Missing file: ", file_path)
      return(NULL)
    }
    
    # Verify BAF column 
    #header <- read_lines(file_path, n_max = 1) %>% str_split("\t") %>% unlist()
    #if(! "baf" %in% tolower(header)) stop("'baf' column not found in ", file_path)
    
    res <- safely_read(file_path, col_types = cols)$result
    if(!is.null(res)) mutate(res, sample_id = .x) else NULL
  })
}

# =====================================================================
# Import data
# =====================================================================
genes <- import_cnvkit_files("genemetrics.txt", cols = cols(
  gene = "c", chromosome = "c", start = "d", end = "d", 
  log2 = "d", depth = "d", weight = "d", probes = "i", baf = "d"
)) %>%
  #filter(!is.na(baf), is.numeric(baf)) %>%
  mutate(
    chromosome = gsub("chr", "", chromosome),  
    chr = factor(
      chromosome,
      levels = c(1:22, "X", "Y"),  
      labels = c(1:22, "X", "Y")   
    ),
    #baf = if_else(baf > 0.5, 1 - baf, baf)
  ) %>%
  filter(!is.na(chr))  # Deletes non valid chromosomes 

# Initial check 
if(nrow(genes) == 0) stop("No data available after filtering!")

# =====================================================================
# Identify key genes 
# =====================================================================

extract_top_altered_genes <- function(genes, log2_threshold = 0.7) {
  top_altered_genes <- genes %>%
    # Filtra per soglia di log2
    dplyr::filter(abs(log2) > log2_threshold) %>%
    # Rimuovi suffissi dopo underscore, ma mantieni altre informazioni
    dplyr::mutate(gene = sub("_.*", "", gene)) %>%
    # Raggruppa per gene
    dplyr::group_by(gene) %>%
    # Aggregazione con più metriche
    dplyr::summarize(
      frequency = dplyr::n_distinct(sample_id),
      # Tipo di alterazione basato sul segno di log2
      type = ifelse(median(log2, na.rm = TRUE) > 0, "AMP", "DEL"),
      # Metriche di base
      log2_median = median(log2, na.rm = TRUE),
      cn_median = median(cn, na.rm = TRUE),
      baf_median = median(baf, na.rm = TRUE),
      # Aggiungi informazioni genomiche
      start = min(start, na.rm = TRUE),
      end = max(end, na.rm = TRUE)) %>%
    
    # Ordina per frequenza decrescente
    dplyr::arrange(dplyr::desc(frequency))
  return(top_altered_genes)
}

log2_threshold = 0.7
new_top_altered_genes <- extract_top_altered_genes(genes, log2_threshold)

  #' Pulisce e normalizza i simboli dei geni da CNVKit genemetrics
  #'
  #' @param data DataFrame con dati di CNVKit contenenti una colonna 'gene'
  #' @param validate_symbols Booleano, se TRUE verifica i simboli contro HGNC (richiede org.Hs.eg.db)
  #' @return DataFrame con simboli dei geni puliti e metriche aggregate
  #' @export
clean_gene_symbols <- function(data, validate_symbols = FALSE) {
    # Verifica la presenza della libreria se richiesta la validazione
    if (validate_symbols && !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      warning("Package org.Hs.eg.db non disponibile. Validazione simboli disattivata.")
      validate_symbols <- FALSE
    }
    
    # Estrai i simboli dei geni con un approccio semplificato
    cleaned_data <- data %>%
      # Estrai il gene symbol con un'unica regex
      dplyr::mutate(
        clean_gene_symbol = gsub("^Intron:", "", gene) %>%  # Rimuovi prefisso Intron:
          gsub(";.*$", "", .)              # Prendi solo la parte prima del primo ;
      ) %>%
      # Filtra via identificatori non-gene e righe vuote
      dplyr::filter(
        !grepl("^ENSG\\d+$|^ENST\\d+|^XM_\\d+|^NM_\\d+|^NR_\\d+$", clean_gene_symbol) &
          nchar(clean_gene_symbol) > 0
      )
    
    # Opzionale: Valida i simboli dei geni contro HGNC
    if (validate_symbols) {
      hgnc_symbols <- as.character(org.Hs.eg.db::org.Hs.egSYMBOL)
      cleaned_data <- cleaned_data %>%
        dplyr::filter(clean_gene_symbol %in% hgnc_symbols)
    }
    
    # Aggrega i dati per gene symbol
    aggregated_data <- cleaned_data %>%
      dplyr::group_by(clean_gene_symbol) %>%
      dplyr::summarize(
        frequency = sum(frequency),
        type = names(which.max(table(type))),
        log2_median = mean(log2_median, na.rm = TRUE),
        cn_median = mean(cn_median, na.rm = TRUE),
        baf_median = mean(baf_median, na.rm = TRUE),
        # Gestisci le colonne genomiche se presenti
        start = if("start" %in% colnames(cleaned_data)) min(start, na.rm = TRUE) else NA,
        end = if("end" %in% colnames(cleaned_data)) max(end, na.rm = TRUE) else NA
      ) %>%
      # Rinomina per coerenza
      dplyr::rename(gene_symbol = clean_gene_symbol) %>%
      # Mantieni l'ordinamento per frequenza
      dplyr::arrange(dplyr::desc(frequency))
    
    return(aggregated_data)
  }
  
clean_genes <- clean_gene_symbols(new_top_altered_genes)
  
output_file <- file.path(config$dirs$out, paste0("log2threshold_", as.character(log2_threshold), "_annotated_genes.xlsx"))
writexl::write_xlsx(annotated_genes, path = output_file)
# ----------------------------------------------------------------------------
# FUNZIONE 1: Pulizia e estrazione informazioni sui geni - CORRETTA DEFINITIVA
# ----------------------------------------------------------------------------
clean_gene_symbols <- function(data, validate_symbols = FALSE, exclude_sex_chr = TRUE) {
  
  if (validate_symbols && !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    warning("org.Hs.eg.db non disponibile. Procedendo senza validazione simboli.")
    validate_symbols <- FALSE
  }
  
  # Verifica che esista la colonna chromosome
  if (!"chromosome" %in% names(data)) {
    stop("La colonna 'chromosome' non è presente nei dati!")
  }
  
  cleaned <- data %>%
    mutate(
      # Estrazione gene symbol - pulizia standard
      gene_symbol = str_replace(gene, "^Intron:", ""),
      gene_symbol = str_replace(gene_symbol, "_.*$", ""),
      gene_symbol = str_replace(gene_symbol, ";.*$", ""),
      
      # USA LA COLONNA CHROMOSOME ESISTENTE - semplice!
      chr = as.character(chromosome)  # Converte in carattere per uniformità
    ) %>%
    # FILTRI di qualità
    filter(
      !is.na(gene_symbol),
      gene_symbol != "",
      gene_symbol != "Deprecated",  # Rimuovi esplicitamente "Deprecated"
      !str_detect(gene_symbol, "^[0-9]+$"),  # Solo numeri
      !str_detect(gene_symbol, "^ENSG|^ENST|^[NX]M_|^[NX]R_"),  # ID database
      !is.na(chr)  # Rimuovi righe senza info cromosoma
    )
  
  # Applica filtro cromosomi sessuali se richiesto
  if (exclude_sex_chr) {
    before_filter <- nrow(cleaned)
    sex_chr_genes <- cleaned %>% filter(chr %in% c("X", "Y"))
    cleaned <- cleaned %>% filter(!chr %in% c("X", "Y"))
    
    cat("Geni cromosomi sessuali rimossi:", before_filter - nrow(cleaned), "\n")
    if (nrow(sex_chr_genes) > 0) {
      cat("Esempi geni X/Y rimossi:", paste(head(sex_chr_genes$gene_symbol, 5), collapse = ", "), "\n")
    }
  }
  
  # Validazione simboli genici se richiesta
  if (validate_symbols) {
    valid_symbols <- AnnotationDbi::mappedkeys(org.Hs.eg.db::org.Hs.egSYMBOL)
    before_validation <- nrow(cleaned)
    cleaned <- filter(cleaned, gene_symbol %in% valid_symbols)
    cat("Geni rimossi per validazione simboli:", before_validation - nrow(cleaned), "\n")
  }
  
  # Report dettagliato
  cat("=== STATISTICHE PULIZIA ===\n")
  cat("Righe originali:", nrow(data), "\n")
  cat("Righe dopo pulizia:", nrow(cleaned), "\n")
  cat("Geni unici:", length(unique(cleaned$gene_symbol)), "\n")
  cat("Campioni:", length(unique(cleaned$sample_id)), "\n")
  
  # Distribuzione cromosomi
  chr_table <- table(cleaned$chr, useNA = "always")
  cat("Distribuzione cromosomi:\n")
  print(chr_table)
  
  # Esempi geni per cromosoma (primi 3 cromosomi)
  cat("\nEsempi geni per cromosoma:\n")
  for (chr_num in head(sort(unique(cleaned$chr)), 3)) {
    genes_chr <- cleaned %>% filter(chr == chr_num) %>% pull(gene_symbol) %>% unique()
    cat("Chr", chr_num, ":", paste(head(genes_chr, 3), collapse = ", "), "\n")
  }
  
  return(cleaned)
}

# ----------------------------------------------------------------------------
# FUNZIONE 2: Estrazione geni alterati - CORRETTA
# ----------------------------------------------------------------------------
extract_top_altered_genes <- function(cleaned_data, 
                                      log2_threshold = 0.7, 
                                      min_samples = 2,
                                      top_n = 50) {
  
  altered_genes <- cleaned_data %>%
    filter(abs(log2) > log2_threshold) %>%
    # Prima aggregazione: per gene/sample
    group_by(gene_symbol, sample_id) %>%
    summarise(
      log2 = mean(log2, na.rm = TRUE),
      cn = mean(cn, na.rm = TRUE),
      baf = mean(baf, na.rm = TRUE),
      genomic_length = sum(end - start + 1, na.rm = TRUE),
      chr = chr[1],  # Prende il primo valore non-NA
      .groups = 'drop'
    ) %>%
    # Seconda aggregazione: per gene
    group_by(gene_symbol) %>%
    summarise(
      chr = chr[1],
      total_samples = n(),
      amp_samples = sum(log2 > 0),
      del_samples = sum(log2 < 0),
      frac_amp = amp_samples / total_samples,
      frac_del = del_samples / total_samples,
      frac_altered = total_samples / length(unique(cleaned_data$sample_id)),
      
      log2_median = median(log2, na.rm = TRUE),
      log2_mean = mean(log2, na.rm = TRUE),
      cn_median = median(cn, na.rm = TRUE),
      baf_median = median(baf, na.rm = TRUE),
      
      genomic_length_total = sum(genomic_length, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      type = case_when(
        frac_amp > 0.7 ~ "AMP",
        frac_del > 0.7 ~ "DEL",
        amp_samples > del_samples ~ "AMP_MIXED",
        del_samples > amp_samples ~ "DEL_MIXED",
        TRUE ~ "BALANCED"
      ),
      alteration_score = frac_altered * abs(log2_median),
      intensity = case_when(
        abs(log2_median) > 2 ~ "HIGH",
        abs(log2_median) > 1 ~ "MODERATE", 
        TRUE ~ "LOW"
      )
    ) %>%
    filter(total_samples >= min_samples) %>%
    arrange(desc(alteration_score)) %>%
    slice_head(n = top_n)
  
  return(altered_genes)
}

# ----------------------------------------------------------------------------
# FUNZIONE 3: Classificazione oncogeni/tumor suppressors
# ----------------------------------------------------------------------------
classify_cancer_genes <- function(altered_genes) {
  
  oncogenes <- c("MYC", "MYCN", "CCND1", "MDM2", "CDK4", "EGFR", "ERBB2", "PIK3CA", "AKT1")
  tumor_suppressors <- c("TP53", "RB1", "CDKN2A", "PTEN", "APC", "BRCA1", "BRCA2", "ATM")
  
  classified <- altered_genes %>%
    mutate(
      gene_class = case_when(
        gene_symbol %in% oncogenes ~ "ONCOGENE",
        gene_symbol %in% tumor_suppressors ~ "TUMOR_SUPPRESSOR",
        TRUE ~ "OTHER"
      )
    )
  
  summary_table <- classified %>%
    count(gene_class, type) %>%
    pivot_wider(names_from = type, values_from = n, values_fill = 0)
  
  return(list(
    classified_genes = classified,
    summary = summary_table,
    top_amplified = filter(classified, type %in% c("AMP", "AMP_MIXED")) %>% slice_head(n = 10),
    top_deleted = filter(classified, type %in% c("DEL", "DEL_MIXED")) %>% slice_head(n = 10)
  ))
}

# ----------------------------------------------------------------------------
# FUNZIONE 4: Grafici - CORRETTI E TESTATI
# ----------------------------------------------------------------------------
create_cnv_plots <- function(altered_genes, genes_data) {
  
  # Plot 1: Volcano plot
  p1 <- altered_genes %>%
    ggplot(aes(x = log2_median, y = frac_altered)) +
    geom_point(aes(color = type, size = total_samples), alpha = 0.7) +
    geom_text(data = filter(altered_genes, frac_altered > 0.3 | abs(log2_median) > 2),
              aes(label = gene_symbol), vjust = -0.5, size = 3, check_overlap = TRUE) +
    scale_color_manual(values = c("AMP" = "red", "DEL" = "blue", 
                                  "AMP_MIXED" = "orange", "DEL_MIXED" = "purple",
                                  "BALANCED" = "gray")) +
    labs(title = "CNV Analysis - Volcano Plot",
         x = "Log2 Copy Number (median)",
         y = "Fraction of Samples Altered",
         color = "Type", size = "Samples") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Plot 2: Sample burden
  sample_stats <- genes_data %>%
    group_by(sample_id) %>%
    summarise(
      amplifications = sum(log2 > 0.7, na.rm = TRUE),
      deletions = sum(log2 < -0.7, na.rm = TRUE),
      total_alterations = amplifications + deletions,
      .groups = 'drop'
    )
  
  p2 <- sample_stats %>%
    pivot_longer(cols = c(amplifications, deletions), 
                 names_to = "type", values_to = "count") %>%
    ggplot(aes(x = reorder(sample_id, total_alterations), y = count, fill = type)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = c("amplifications" = "red", "deletions" = "blue")) +
    labs(title = "CNV Burden per Sample",
         x = "Sample ID", y = "Number of Alterations",
         fill = "Alteration Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Plot 3: Top genes barplot
  p3 <- altered_genes %>%
    slice_head(n = 20) %>%
    ggplot(aes(x = reorder(gene_symbol, alteration_score), y = alteration_score, fill = type)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("AMP" = "red", "DEL" = "blue", 
                                 "AMP_MIXED" = "orange", "DEL_MIXED" = "purple")) +
    labs(title = "Top 20 Altered Genes",
         x = "Gene", y = "Alteration Score",
         fill = "Type") +
    theme_minimal()
  
  return(list(
    volcano_plot = p1, 
    sample_burden = p2,
    top_genes = p3
  ))
}

# ----------------------------------------------------------------------------
# FUNZIONE 5: Report riassuntivo
# ----------------------------------------------------------------------------
generate_cnv_report <- function(altered_genes, gene_classification) {
  
  cat("=== CNV ANALYSIS REPORT ===\n\n")
  
  cat("OVERVIEW:\n")
  cat("- Total altered genes:", nrow(altered_genes), "\n")
  cat("- Amplifications:", sum(str_detect(altered_genes$type, "AMP")), "\n")
  cat("- Deletions:", sum(str_detect(altered_genes$type, "DEL")), "\n\n")
  
  cat("CANCER-RELATED GENES:\n")
  cancer_genes <- filter(gene_classification$classified_genes, gene_class != "OTHER")
  if (nrow(cancer_genes) > 0) {
    for (i in 1:nrow(cancer_genes)) {
      gene <- cancer_genes[i, ]
      cat(sprintf("- %s (%s): %s (%.1f%% samples, log2=%.2f)\n", 
                  gene$gene_symbol, gene$gene_class, gene$type, 
                  gene$frac_altered * 100, gene$log2_median))
    }
  } else {
    cat("- No known cancer genes found in top alterations\n")
  }
  
  cat("\nTOP 10 MOST ALTERED GENES:\n")
  top10 <- head(altered_genes, 10)
  for (i in 1:nrow(top10)) {
    gene <- top10[i, ]
    cat(sprintf("%2d. %s: %s (%.1f%% samples, score=%.3f)\n", 
                i, gene$gene_symbol, gene$type, 
                gene$frac_altered * 100, gene$alteration_score))
  }
}

# ----------------------------------------------------------------------------
# ESEMPIO COMPLETO DI UTILIZZO
# ----------------------------------------------------------------------------

# Esegui l'analisi completa:
run_complete_cnv_analysis <- function(genes_data) {
  
  cat("Inizio analisi CNV...\n")
  
  # Step 1: Pulizia dati
  clean_genes <- clean_gene_symbols(genes_data, validate_symbols = FALSE)
  
  # Step 2: Estrazione geni alterati
  top_altered <- extract_top_altered_genes(clean_genes, log2_threshold = 0.7)
  
  # Step 3: Classificazione geni oncologici
  gene_classification <- classify_cancer_genes(top_altered)
  
  # Step 4: Grafici
  plots <- create_cnv_plots(top_altered, clean_genes)
  
  # Step 5: Report
  generate_cnv_report(top_altered, gene_classification)
  
  return(list(
    clean_data = clean_genes,
    altered_genes = top_altered,
    gene_classification = gene_classification,
    plots = plots
  ))
}

# USO:
results <- run_complete_cnv_analysis(genes)
print(results$plots$volcano_plot)
View(results$altered_genes)
