library(readr)
library(dplyr)
library(stringr)
library(ggplot2)

gff <- read_tsv("D:/Project/Spladder.gff3", 
  comment = "#",
  col_names = c("seqid","source","type","start","end","score","strand","phase","attributes")
)


events <- gff %>%
  filter(type == "gene") %>%
  transmute(
    event_id  = str_match(attributes, "(?<=ID=)[^;]+")[,1],  # e.g. alt_3prime.6
    event_src = source
  ) %>%
  distinct(event_id, .keep_all = TRUE)

events <- events %>%
  mutate(event_type = case_when(
    event_src %in% c("exon_skip","mult_exon_skip")          ~ "Exon skipping",
    event_src == "intron_retention"  ~ "Intron retention",
    event_src %in% c("mutually_exclusive","mutex_exons","MXE","mxe") ~ "Mutually exclusive exons",
    event_src == "alt_3prime"        ~ "Alt-3'splice site",
    event_src == "alt_5prime"        ~ "Alt-5'splice site",
  
    TRUE                             ~ event_src              # 其他保持原名
  ))


df_counts <- events %>%
  count(event_type, name = "count") %>%
  arrange(count)

print(df_counts)

p <- ggplot(df_counts, aes(x = reorder(event_type, count), y = count)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = count), hjust = -0.2, size = 3) +
  coord_flip(ylim = c(0, max(df_counts$count)*1.15)) +
  labs(x = "Event type", y = "Count", title = "Alt-splicing event counts in HCT116") +
  theme_classic(base_size = 11)

print(p)

ggsave("D:/Project/event_counts.pdf", p, width = 6, height = 4)



events_pos <- gff %>%
  filter(type == "gene") %>%
  transmute(
    event_id = str_match(attributes, "(?<=ID=)[^;]+")[,1],
    seqid,
    start = as.integer(start),
    end   = as.integer(end)
  ) %>%
  distinct(event_id, .keep_all = TRUE) %>%
  inner_join(events %>% select(event_id, event_type), by = "event_id")

library(readr)

df_chr <- events_pos %>% count(seqid, event_type, name = "count")

chr_levels <- df_chr %>% distinct(seqid) %>% arrange(parse_number(seqid)) %>% pull(seqid)
df_chr <- df_chr %>% mutate(seqid = factor(seqid, levels = chr_levels))

df_chr <- df_chr %>% mutate(seqid = factor(seqid, levels = chr_levels))

totals <- df_chr %>%
  group_by(seqid) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  mutate(seqid = factor(seqid, levels = levels(df_chr$seqid)))

library(ggplot2)

p_chr <- ggplot(df_chr, aes(x = seqid, y = count, fill = event_type)) +
  geom_col() +

  geom_text(data = totals,
            aes(x = seqid, y = total, label = total),
            vjust = -0.3, size = 2.2, inherit.aes = FALSE) +

  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_fill_brewer(palette = "Pastel1") +
  labs(x = "Chromosome", y = "Event count", fill = "Event type",
       title = "Alt-splicing events per chromosome in HCT116") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_chr)

ggsave("D:/Project/event_counts_by_chr_pastel.pdf",  p_chr, width = 8, height = 4)






