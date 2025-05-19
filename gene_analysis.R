rm(list=ls())

library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggpubr)
library(ggplot2)
library(writexl)

config <- list(
  dirs = list(
    base = "D:/CNVkit/tumor/PTJ_WES_IDT-30802789/",
    out = "D:/CNVkit/tumor/tumor_imgs/",
    res = "D:/CNVkit/tumor/tumor_res/"
  ) 
)

config <- list(
  dirs = list(
    base = "D:/CNVkit/model/WES_modelli/",
    out = "D:/CNVkit/model/with_pooledref/model_imgs/",
    res = "D:/CNVkit/model/with_pooledref/model_res//"
  ) 
)

# =====================================================================
# Function to import CNVkit files
# =====================================================================
import_cnvkit_files <- function(file_type, cols) {
  safely_read <- safely(read_tsv)
  
  samples <- list.files(config$dirs$res, pattern = ".*2D-001.*$") %>% 
    setdiff(c("multi_sample", "nfr", "nfrmulti_sample"))
  
  map_dfr(samples, ~{
    file_path <- file.path(config$dirs$res, .x, paste0(.x, "_", file_type))
    if(!file.exists(file_path)) {
      message("File mancante: ", file_path)
      return(NULL)
    }
    
    # Verify BAF column 
    header <- read_lines(file_path, n_max = 1) %>% str_split("\t") %>% unlist()
    if(! "baf" %in% tolower(header)) stop("Colonna 'baf' non trovata in ", file_path)
    
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
    baf = if_else(baf > 0.5, 1 - baf, baf)
  ) %>%
  filter(!is.na(chr))  # Deletes non valid chromosomes 

# Initial check 
if(nrow(genes) == 0) stop("Nessun dato dopo il filtraggio BAF! Verifica i file di input.")

# =====================================================================
# Identify key genes 
# =====================================================================

extract_top_altered_genes <- function(genes, log2_threshold = 0.7) {
  top_altered_genes <- genes %>%
    dplyr::filter(abs(log2) > log2_threshold) %>%
    dplyr::mutate(gene = sub("_.*", "", gene)) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarize(
      frequency = dplyr::n_distinct(sample_id),
      # Label alteration type (amplification, deletion) based on log2 sign 
      type = ifelse(median(log2, na.rm = TRUE) > 0, "AMP", "DEL"),
      log2_median = median(log2, na.rm = TRUE),
      cn_median = median(cn, na.rm = TRUE),
      baf_median = median(baf, na.rm = TRUE)
    ) %>%
    dplyr::arrange(dplyr::desc(frequency))
  
  return(top_altered_genes)
}

log2_threshold = 0
top_altered_genes <- extract_top_altered_genes(genes, log2_threshold)

annotated_genes <- top_altered_genes %>%
  filter(gene %in% keys(org.Hs.eg.db, keytype="SYMBOL"))

output_file <- file.path(config$dirs$out, paste0("log2threshold_", as.character(log2_threshold), "_annotated_genes.xlsx"))
writexl::write_xlsx(annotated_genes, path = output_file)

# =====================================================================
# Domanda Biologica 1: Le alterazioni copiche sono associate a perdita di eterozigosi (LOH)?
# =====================================================================
plot_baf_vs_cnv <- genes %>%
  filter(abs(log2) > 0.7) %>%
  ggplot(aes(x = log2, y = baf, color = chr)) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  geom_vline(xintercept = c(-0.7, 0.7), linetype = "dashed") +
  geom_hline(yintercept = 0.1, color = "red", linetype = "dotted") + # Soglia LOH
  facet_wrap(~gene, scales = "free_x") +
  labs(title = "BAF vs Alterazioni Copiche",
       subtitle = "BAF < 0.1 suggerisce LOH nelle delezioni (log2 < -0.7)",
       x = "log2 ratio (CNVkit)",
       y = "BAF")

plot_baf_vs_cnv <- genes %>%
  filter(abs(log2) > 0.7) %>%
  ggplot(aes(x = log2, y = baf, color = chr)) +
  geom_jitter(width = 0.1, alpha = 0.6, na.rm = TRUE) + # Doppio controllo
  geom_vline(xintercept = c(-0.7, 0.7), linetype = "dashed") +
  geom_hline(yintercept = 0.1, color = "red", linetype = "dotted") +
  facet_wrap(~gene, scales = "free_x") +
  labs(title = "BAF vs Alterazioni Copiche (NA rimossi)")

plot_baf_vs_cnv <- genes %>%
  filter(abs(log2) > 0.7, gene %in% c("CFC1B", "CSAG3", "GPRIN2")) %>%  # Filtra geni
  ggplot(aes(x = log2, y = baf, color = chr)) +
  geom_point(alpha = 0.6, size = 1) +  # Più efficiente di jitter
  geom_vline(xintercept = c(-0.7, 0.7), linetype = "dashed") +
  geom_hline(yintercept = 0.1, color = "red", linetype = "dotted") +
  labs(title = "BAF vs CNV (Geni selezionati)") +
  theme_minimal()

print(plot_baf_vs_cnv)  # Forza il rendering esplicito

# =====================================================================
# Domanda Biologica 2: Esistono cromosomi con BAF anormale diffuso?
# =====================================================================
plot_chr_baf <- genes %>%
  group_by(sample_id, chr) %>%
  summarise(
    median_baf = median(baf, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = chr, y = median_baf, fill = chr)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.25, color = "blue", linetype = "dashed") +
  labs(title = "BAF Medio per Cromosoma",
       y = "BAF mediano",
       x = "Cromosoma") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_chr_baf <- genes %>%
  group_by(sample_id, chr) %>%
  summarise(
    median_baf = median(baf, na.rm = TRUE), 
    .groups = "drop"
  ) %>%
  ggplot(aes(x = chr, y = median_baf, fill = chr)) +
  geom_boxplot(na.rm = TRUE)

plot_chr_baf <- genes %>%
  group_by(sample_id, chr) %>%  # Usa "chr" corretto
  summarise(median_baf = median(baf, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = chr, y = median_baf, fill = chr)) +  # Asse x = chr
  geom_boxplot() +
  labs(x = "Cromosoma") +
  theme(legend.position = "none") 
# ====================================================================
# Visualizzazione e salvataggio file                                                                                       
# ====================================================================
print(plot_baf_vs_cnv)
print(plot_chr_baf)

ggsave(file.path(config$dirs$out, "baf_vs_cnv.pdf"), plot_baf_vs_cnv, width = 14, height = 8)
ggsave(file.path(config$dirs$out, "chr_baf.pdf"), plot_chr_baf, width = 14, height = 8)
