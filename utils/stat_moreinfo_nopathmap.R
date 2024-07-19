library(readr)
library(dplyr)
library(stringr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No input files provided as an argument.")
}

dir <- normalizePath(args[1])
gap <- normalizePath(args[2])
telom <- normalizePath(args[3])
telom_cutoff <- args[4]
telom_length <- normalizePath(args[7])

t_cut <- as.numeric(telom_cutoff)
print(t_cut)

setwd(dir)
print(getwd())

df <- read_delim("all_chr_assembly.tsv")
coor <- read_delim(gap, col_names = c("chr","gapstart","gapend"))
telomere <- read_delim(telom, col_names = c("chr","start","end","value"))

## for gaps
coor <- coor %>% mutate(gaplen = gapend - gapstart +1) %>%
  group_by(chr) %>%
  summarise(gapcount = n(), gapsize = sum(gaplen))

## for telomeres
tlm <- telomere %>%
  group_by(chr) %>%
  slice(c(1:2, (n() - 1):n())) %>%
  summarise(p = sum(value[1:2]),
            q = sum(value[(n() - 1):n()]))

## for telomere length with telomere_analysis.sh
tlm_len <- telomere_length %>%
  mutate(p = ifelse(start < 1, end - start, NA),
         q = ifelse(start > 0, end - start, NA),
         p_start = ifelse(start < 1, Start, NA),
         p_end = ifelse(start < 1, end, NA),
         q_start = ifelse(start > 0, start, NA),
         q_end = ifelse(start > 0, end, NA)) %>%
  group_by(Chromosome) %>%
  summarise(p_length = sum(p, na.rm = TRUE),
            q_length = sum(q, na.rm = TRUE),
            p_start = first(p_start[!is.na(p_start)]),
            p_end = first(p_end[!is.na(p_end)]),
            q_start = first(q_start[!is.na(q_start)]),
            q_end = first(q_end[!is.na(q_end)])) %>%
  ungroup()

df_new <- merge(df, coor, by = "chr", all.x = TRUE)
df_new <-  merge(df_new, tlm, by = "chr", all.x = TRUE)
df_new <- merge(df_new, tlm_len, by = "chr", all.x = TRUE)

# creating the table

df_final <- df_new %>%
  mutate(telomere = case_when(
    !is.na(p_len) & !is.na(q_len) ~ "pq",
    !is.na(p_len) ~ "p",
    !is.na(q_len) ~ "q",
    TRUE ~ "0"
    ),
    completion = case_when(
      telomere == "pq" & is.na(gapcount) ~ "T_2_T",
      telomere == "pq" & !is.na(gapcount) ~ "T_gap_T",
      telomere == "p" & is.na(gapcount) ~ "T_nogap_noT",
      telomere == "p" & !is.na(gapcount) ~ "T_gap_noT",
      telomere == "q" & is.na(gapcount) ~ "noT_nogap_T",
      telomere == "q" & !is.na(gapcount) ~ "noT_gap_T",
      telomere == "0" & is.na(gapcount) ~ "noT_nogap_noT",
      TRUE ~ "noT_gap_noT"
      ),
      score = case_when(
        telomere == "pq" ~ "2",
        telomere == "p" ~ "1",
        telomere == "q" ~ "1",
        TRUE ~ "0"
      )) %>%
  select(chr,
         completion,
         score,
         contig_name = query_name, 
         contig, 
         query_length, 
         ref_length, 
         contig_align_proportion, 
         queryref_proportion,
         missing_contig,
         gapcount,
         total_gap_bp = gapsize,
         telomere,
         p_count = p,
         q_count = q) %>%
  arrange(factor(chr, levels = c(1:29, "X", "Y")))

df_final[is.na(df_final)] <- 0

dirbase <- basename(dir)
write_tsv(df_final, file = paste0(dir,"/",dirbase,"_all_STATS.tsv"))
