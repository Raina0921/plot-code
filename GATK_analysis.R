if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("VariantAnnotation")
install.packages(c("dplyr","tidyr","ggplot2"))


vcf_file <- "D:/Project/HCT116_Variants.vcf"
keep_pass <- TRUE          
treat_dot_as_pass <- FALSE 
plot_title <- "Variant counts in HCT116"
x_label <- "Count"
y_label <- "Variant type"
out_counts_tsv <- "variant_type_counts.tsv"
out_plot_png   <- "variant_type_counts_like_example.png"

library(VariantAnnotation)
library(S4Vectors)
library(ggplot2)

vcf <- readVcf(vcf_file)

ALT_list   <- alt(vcf)
ALT_vec    <- toupper(as.character(unlist(ALT_list)))
REF_vec    <- toupper(as.character(ref(vcf)))
REF_rep    <- rep(REF_vec, elementNROWS(ALT_list))
FILTER_vec <- fixed(vcf)$FILTER
FILTER_rep <- rep(FILTER_vec, elementNROWS(ALT_list))

df <- data.frame(REF = REF_rep, ALT = ALT_vec, FILTER = FILTER_rep,
                 stringsAsFactors = FALSE)

if (keep_pass) {
  pass_set <- if (treat_dot_as_pass) c("PASS",".") else "PASS"
  df <- df[is.na(df$FILTER) | df$FILTER %in% pass_set, , drop = FALSE]
}


df <- df[!is.na(df$ALT) & df$ALT != "<NON_REF>", , drop = FALSE]


is_symbolic <- grepl("^<.*>$", df$ALT) | grepl("\\[|\\]", df$ALT)
sv_from_alt <- ifelse(grepl("^<.*>$", df$ALT),
                      sub("^<([^>]+)>$", "\\1", df$ALT),
                      ifelse(grepl("\\[|\\]", df$ALT), "BND", NA_character_))
svtype_info <- if ("SVTYPE" %in% colnames(info(vcf))) as.character(info(vcf)$SVTYPE) else NULL
svtype_rep  <- if (!is.null(svtype_info)) rep(svtype_info, elementNROWS(ALT_list)) else NA_character_


rlen <- nchar(df$REF); alen <- nchar(df$ALT)
TYPE <- rep(NA_character_, nrow(df))
TYPE[df$REF!="*" & df$ALT!="*" & !is_symbolic & rlen==1 & alen==1] <- "SNP"
TYPE[df$REF!="*" & df$ALT!="*" & !is_symbolic & alen>rlen]        <- "Insertion"
TYPE[df$REF!="*" & df$ALT!="*" & !is_symbolic & alen<rlen]        <- "Deletion"
TYPE[df$REF!="*" & df$ALT!="*" & !is_symbolic & rlen>1 & alen>1 & rlen==alen] <- "MNP"

fill_idx <- is.na(TYPE)
TYPE[fill_idx & !is.na(svtype_rep)]  <- svtype_rep[fill_idx & !is.na(svtype_rep)]
TYPE[fill_idx & !is.na(sv_from_alt)] <- sv_from_alt[fill_idx & !is.na(sv_from_alt)]
TYPE[is.na(TYPE)] <- "Unknown"
TYPE <- toupper(sub(":.*$", "", TYPE))  # "DUP:TANDEM" → "DUP"


tab <- table(TYPE)
counts <- data.frame(TYPE = names(tab), Count = as.integer(tab), row.names = NULL)

write.table(counts, out_counts_tsv, sep = "\t", quote = FALSE, row.names = FALSE)


counts$TYPE <- reorder(counts$TYPE, counts$Count)


lab_y <- function(x){
  x <- as.character(x)
  x[x == "INSERTION"] <- "Insertion"
  x[x == "DELETION"]  <- "Deletion"
  x
}

p <- ggplot(counts, aes(y = reorder(TYPE, Count), x = Count)) +
  geom_col(width = 0.5, fill = "grey40") +
  geom_text(aes(label = Count), hjust = -0.15, size = 4, color = "grey10") +
  labs(title = plot_title, x = x_label, y = y_label) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line.y.left = element_line(size = 1.2),
    axis.line.x = element_line(size = 0.6),
    axis.ticks.y = element_blank()
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  scale_y_discrete(labels = lab_y) +  
  coord_cartesian(clip = "off")

print(p)

out_plot_png   <- "D:/variant_type_counts.png"
ggsave(out_plot_png, p, width = 7.5, height = 4.5, dpi = 300)
