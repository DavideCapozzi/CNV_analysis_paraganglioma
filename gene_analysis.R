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
    out = "D:/CNVkit/tumor/tumor_imgs",
    res = "D:/CNVkit/tumor/tumor_res/"
  ) 
)

# =====================================================================
# Funzione modificata per importare genemetrics CON BAF
# Verifica struttura file genemetrics (CNVkit v0.9.11):
# gene, chromosome, start, end, log2, depth, weight, probes, baf (se presente)
# =====================================================================
import_cnvkit_files <- function(file_type, cols) {
  safely_read <- safely(read_tsv)
  
  samples <- list.files(config$dirs$res, pattern = ".*-t.*$") %>% 
    setdiff(c("multi_sample", "nfr", "nfrmulti_sample"))
  
  map_dfr(samples, ~{
    file_path <- file.path(config$dirs$res, .x, paste0(.x, "_", file_type))
    if(!file.exists(file_path)) {
      message("File mancante: ", file_path)
      return(NULL)
    }
    
    # Verifica presenza colonna BAF
    header <- read_lines(file_path, n_max = 1) %>% str_split("\t") %>% unlist()
    if(! "baf" %in% tolower(header)) stop("Colonna 'baf' non trovata in ", file_path)
    
    res <- safely_read(file_path, col_types = cols)$result
    if(!is.null(res)) mutate(res, sample_id = .x) else NULL
  })
}

# =====================================================================
# Importazione dati 
# =====================================================================
genes <- import_cnvkit_files("genemetrics.txt", cols = cols(
  gene = "c", chromosome = "c", start = "d", end = "d", 
  log2 = "d", depth = "d", weight = "d", probes = "i", baf = "d"
)) %>%
  #filter(!is.na(baf), is.numeric(baf)) %>%
  mutate(
    # Correzione critica: rimuovi "chr" e converti in factor
    chromosome = gsub("chr", "", chromosome),  # Rimuove "chr" dai valori
    chr = factor(
      chromosome,
      levels = c(1:22, "X", "Y"),  # Livelli corretti SENZA "chr"
      labels = c(1:22, "X", "Y")   # Etichette esplicite (opzionale)
    ),
    baf = if_else(baf > 0.5, 1 - baf, baf)
  ) %>%
  filter(!is.na(chr))  # Rimuove eventuali cromosomi non validi (es. "M")

# Verifica iniziale
if(nrow(genes) == 0) stop("Nessun dato dopo il filtraggio BAF! Verifica i file di input.")

# =====================================================================
# Identificazione geni chiave
# =====================================================================

top_altered_genes <- genes %>%
  filter(abs(log2) > 0.7) %>%
  mutate(gene = sub("_.*", "", gene)) %>%
  group_by(gene) %>%
  summarize(
    frequency = n_distinct(sample_id),
    # Determina il tipo di alterazione basandosi sulla mediana di log2
    type = ifelse(median(log2, na.rm = TRUE) > 0, "AMP", "DEL"),
    log2_median = median(log2, na.rm = TRUE),
    cn_median = median(cn, na.rm = TRUE),
    baf_median = median(baf, na.rm = TRUE)
  ) %>%
  arrange(desc(frequency))

geni_annotati <- top_altered_genes %>%
  filter(gene %in% keys(org.Hs.eg.db, keytype="SYMBOL"))

output_file <- file.path(config$dirs$out, "geni_annotati.xlsx")
writexl::write_xlsx(geni_annotati, path = output_file)

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
