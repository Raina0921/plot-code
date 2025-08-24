library(readr)
library(data.table)
library(ggplot2)

dt <- fread("D:\\Project\\s_k\\HCT116_Test_transcript_abundances.tsv", sep = "\t")

thresholds <- c(2, 5, 10)

res <- rbindlist(
  lapply(thresholds, function(th) {
    dt[kallisto_tpm >= th, 
       .N,            
       by = biotype
    ][, threshold := th]
  }))

setDT(res)
res[, label_pos := N * 1]       
max_lim <- max(res$label_pos) * 1.02 
res[, threshold_lbl := factor(paste0("TPM\u2265", threshold),
                              levels = paste0("TPM\u2265", thresholds))]
kallisto_transcript_counts <- ggplot(res, aes(x = reorder(biotype, N), y = N)) +
  geom_col() +
  geom_text(aes(y = label_pos, label = N),
            hjust = -0.1, size = 2.5) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    limits = c(0, max_lim),
    expand  = expansion(mult = c(0, 0))
  ) +
  facet_wrap(~threshold_lbl) +
  labs(
    title = "Transcript counts per biotype above TPM thresholds (kallisto)",
    x = "Biotype", y = "Number of transcripts"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text  = element_text(size = 8),
    axis.title = element_text(size = 9),
    plot.title = element_text(size = 11),
    plot.margin = margin(5.5, 32, 5.5, 5.5)
  )

print(kallisto_transcript_counts)
ggsave(
  filename = "D:/project/kallisto_transcript_counts.pdf",
  plot     = kallisto_transcript_counts,
  width    = 10, height = 5,
  device   = cairo_pdf    
)



res <- rbindlist(
  lapply(thresholds, function(th) {
    dt[stringtie_tpm >= th, 
       .N,            
       by = biotype
    ][, threshold := th]
  }))

setDT(res)
res[, label_pos := N * 1]       
max_lim <- max(res$label_pos) * 1.02 
res[, threshold_lbl := factor(paste0("TPM\u2265", threshold),
                              levels = paste0("TPM\u2265", thresholds))]
stringtie_transcript_counts <- ggplot(res, aes(x = reorder(biotype, N), y = N)) +
  geom_col() +
  geom_text(aes(y = label_pos, label = N),
            hjust = -0.1, size = 2.5) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    limits = c(0, max_lim),
    expand  = expansion(mult = c(0, 0))
  ) +
  facet_wrap(~threshold_lbl) +
  labs(
    title = "Transcript counts per biotype above TPM thresholds (stringtie)",
    x = "Biotype", y = "Number of transcripts"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text  = element_text(size = 8),
    axis.title = element_text(size = 9),
    plot.title = element_text(size = 11),
    plot.margin = margin(5.5, 32, 5.5, 5.5)
  )

print(stringtie_transcript_counts)
ggsave(
  filename = "D:/project/stringtie_transcript_counts.pdf",
  plot     = stringtie_transcript_counts,
  width    = 10, height = 5,
  device   = cairo_pdf    
)

pkgs <- c("data.table", "ggplot2", "ggpubr", "eulerr", "janitor")
new  <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(new)) install.packages(new)
lapply(pkgs, library, character.only = TRUE)

dt <- fread("D:\\Project\\s_k\\HCT116_Test_transcript_abundances.tsv", sep = "\t")
setDT(dt)                       
dt <- clean_names(dt)          


scatter_dt <- dt[!(kallisto_tpm == 0 & stringtie_tpm == 0)]

scatter_dt[, `:=`(log_kallisto = log10(kallisto_tpm + 1),
                  log_stringtie = log10(stringtie_tpm + 1))]

corr <- cor(scatter_dt$kallisto_tpm, scatter_dt$stringtie_tpm, method = "pearson")


tpm_cutoff <- 10 # ← 2、5、10
setA <- dt[kallisto_tpm   >= tpm_cutoff, transcript_id] 
setB <- dt[stringtie_tpm  >= tpm_cutoff, transcript_id] 

venn_fit <- eulerr::euler(list(Kallisto = setA,
                               Stringtie = setB))

p_venn <- plot(venn_fit,
     fills  = list(fill = c("#FC8D62", "#8DA0CB"), alpha = 0.6),
     labels = c("Kallisto", "StringTie"),
     legend = TRUE,
     quantities = TRUE,
     main = paste0("Transcripts with TPM ≥ ", tpm_cutoff))


print(p_venn)

ggsave(
  filename = "D:/project/k_vs_s_venn10.pdf",
  plot     = p_venn,
  width    = 5, height = 5,
  device   = cairo_pdf    
)


# 5 TPM threshold, protein coding and non-protein coding transcripts separately
         
library(tidyverse)
library(scales)

df <- read_tsv("D:\\Project\\s_k\\HCT116_Test_transcript_abundances.tsv", col_types = cols())

df <- df %>% 
  mutate(
    coding_status = case_when(
      str_detect(biotype, regex("protein_coding", ignore_case = TRUE)) ~ "Protein-coding",
      TRUE ~ "Non-coding"
    ),
    across(c(kallisto_tpm, stringtie_tpm), as.numeric)
  )


tpm_threshold <- 5
df_filt <- df %>% 
  filter(kallisto_tpm >= tpm_threshold | stringtie_tpm >= tpm_threshold)

make_scatter <- function(data, title_suffix) {
  
  r_val <- cor(data$kallisto_tpm, data$stringtie_tpm, use = "complete.obs")
  
  annotate("text",
           x = Inf, y = Inf,                      
           label = paste0("r = ", round(r_val, 3)),
           hjust = 1.1, vjust = 1.2, size = 4)  
  
  x_annot <- min(data$kallisto_tpm, na.rm = TRUE) * 1.05
  y_annot <- max(data$stringtie_tpm,  na.rm = TRUE) / 1.05
  
  ggplot(data, aes(kallisto_tpm, stringtie_tpm)) +

    geom_point(size = 0.4, alpha = 0.4, shape = 16) +
    
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    
  geom_smooth(
    method  = "lm",
    formula = y ~ x,          
    se      = FALSE,
    linewidth = 0.7,
    colour    = "#d73027"    
  ) +
    scale_x_log10(
      limits = c(tpm_threshold, NA),
      expand = expansion(mult = c(0.02, 0.05)),
      labels = comma_format(accuracy = 1),
      name   = "Kallisto TPM (log10)"
    ) +
    scale_y_log10(
      limits = c(tpm_threshold, NA),
      expand = expansion(mult = c(0.02, 0.05)),
      labels = comma_format(accuracy = 1),
      name   = "StringTie TPM (log10)"
    ) +
    coord_equal() +
    annotate("text",
             x = Inf, y = Inf,
             label = paste0("r = ", round(r_val, 3)),
             hjust = 1.05, vjust = 1.2, size = 4) +
    labs(
      title    = paste(title_suffix, "transcripts with TPM ≥ ", tpm_threshold),
    ) +
    
    theme_bw(base_size = 12) +
    theme(
      plot.title   = element_text(size = 12, face = "bold"),
      axis.title.x  = element_text(size = 10), 
      axis.title.y  = element_text(size = 10),
      plot.margin  = margin(15, 15, 15, 15),
      panel.grid.minor = element_blank()
    )
}

p_coding    <- make_scatter(df_filt %>% filter(coding_status == "Protein-coding"),
                            "Protein-coding")
p_noncoding <- make_scatter(df_filt %>% filter(coding_status == "Non-coding"),
                            "Non-coding")

library(patchwork)

p_both <- p_coding + p_noncoding +            
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a",
                  theme = theme(
                    plot.tag = element_text(size = 8, face = "bold")) 
                  )
print(p_both)
ggsave("D:/Project/s_k/tpm_scatter_both.pdf",
       p_both, width = 10, height = 5, device = cairo_pdf)
