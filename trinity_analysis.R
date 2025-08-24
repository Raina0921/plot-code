library(Biostrings) 
library(tidyverse) 

tr_fa   <- "D:/Project/HCT116_trinity_output/trinity.fasta"
gg_tr_fa<- "D:/Project/HCT116_trinity_output/GGtrinity.fasta"
orf_fa  <- "D:/Project/HCT116_trinity_output/trinity_orf.fasta"
gg_orf_fa<- "D:/Project/HCT116_trinity_output/GGtrinity_orf.fasta"

transcript   <- readRNAStringSet(tr_fa)
GGtranscript <- readRNAStringSet(gg_tr_fa)
ORF          <- readAAStringSet(orf_fa)
GGORF        <- readAAStringSet(gg_orf_fa)

head(names(transcript))
head(names(GGtranscript))
head(names(ORF))
head(names(GGORF))

print(
  tibble(
    object      = c("transcript", "GGtranscript", "ORF", "GGORF"),
    n_sequences = c(length(transcript),
                    length(GGtranscript),
                    length(ORF),
                    length(GGORF))
  )
)

counts <- tibble(
  Method = c("Trinity (de novo)", "Trinity (genome-guided)"),
  Transcript = c(length(transcript),  length(GGtranscript)),
  ORF    = c(length(ORF),         length(GGORF))
) %>% 
  pivot_longer(cols = c(Transcript, ORF),
               names_to = "Type",
               values_to = "Count")

library(stringr)
library(patchwork)

counts_tbl <- tibble(
  Method     = c("De novo", "Genome-guided"),
  Transcript = c(length(transcript),  length(GGtranscript)),
  ORF        = c(length(ORF),         length(GGORF))
) |>
  pivot_longer(c(Transcript, ORF),
               names_to = "Type", values_to = "Count") |>
  mutate(
    Method = stringr::str_replace_all(Method, "\u00A0", " "),
    Method = factor(Method, levels = c("De novo", "Genome-guided")),
    Type   = factor(Type,   levels = c("Transcript", "ORF"))
  )

pal <- c("De novo"       = "#8DA0CB",
         "Genome-guided" = "#A6D854")

p_tx <- counts_tbl |>
  dplyr::filter(Type == "Transcript") |>
  ggplot(aes(x = Method, y = Count, fill = Method)) +
  geom_col(width = 0.35, colour = "grey30", linewidth = 0.3) +
  geom_text(aes(label = Count), vjust = -0.25, size = 3) +
  scale_fill_manual(values = pal, name = NULL) +
  labs(title = "Transcript counts", x = "Trinity method", y = "Sequence number (n)") +
  theme_minimal(base_size = 10) +
  theme(plot.title  = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 0, vjust = 0.8),
        legend.position = "none")

p_orf <- counts_tbl |>
  dplyr::filter(Type == "ORF") |>
  ggplot(aes(x = Method, y = Count, fill = Method)) +
  geom_col(width = 0.35, colour = "grey30", linewidth = 0.3) +
  geom_text(aes(label = Count), vjust = -0.25, size = 3) +
  scale_fill_manual(values = pal, name = NULL) +
  labs(title = "ORF counts", x = "Trinity method", y = "Sequence number (n)") +
  theme_minimal(base_size = 10) +
  theme(plot.title  = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 0, vjust = 0.8),
        legend.position = "bottom")

y_lim <- ceiling(max(counts_tbl$Count, na.rm = TRUE) * 1.10)  

brks  <- scales::pretty_breaks(n = 6)(c(0, y_lim))

p_tx  <- p_tx  + scale_y_continuous(limits = c(0, y_lim),
                                    breaks = brks,
                                    expand = expansion(mult = c(0, 0)))
p_orf <- p_orf + scale_y_continuous(limits = c(0, y_lim),
                                    breaks = brks,
                                    expand = expansion(mult = c(0, 0)))

library(patchwork)

p_both <- (p_tx | p_orf) +
  plot_layout(ncol = 2) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(size = 9, face = "bold"),
      plot.tag.position = c(0.01, 0.99)
    )
  )

print(p_both)

ggsave("D:/Project/Plots/trinity_counts.pdf", p_both, width = 8, height = 5, device = cairo_pdf)



library(readxl)
library(janitor)
library(eulerr)
library(ggplot2)
library(dplyr)

xlsx <- "D:/project/overlap.xlsx"
df <- read_excel(xlsx, sheet = 1) |> clean_names()

if ("count" %in% names(df)) df <- select(df, -count)

vals <- df |>
  select(de_novo_only, intersection, genome_guided_only) |>
  slice(1) |>
  unlist(use.names = TRUE)


fit <- euler(c(
  de_novo      = vals[["de_novo_only"]],
   genome_guided         = vals[["genome_guided_only"]],
  "de_novo&genome_guided" = vals[["intersection"]]
))

p <- plot(
  fit,
  fills  = list(fill = c("#4E79A7", "#F28E2B"), alpha = 0.6),
  edges  = list(col = "grey30", lwd = 1),
  labels = FALSE,         # 圈旁标签（denovo / gg）
  quantities = list(cex = 0.8),     # 区域里的数字
  legend = list(cex = 0.8)         # 图例文字
)

print(p)
ggsave("D:/project/trinity_venn_from_excel.pdf", p, width = 5, height = 5, device = cairo_pdf)