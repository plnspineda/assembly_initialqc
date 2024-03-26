library(readr)
library(dplyr)
library(stringr)
library(tidyr)

## This is version 2.1 where I changed the categories to three instances e.g. T2T, TgapT, noT2T, and so on

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No input files provided as an argument.")
}

dir <- normalizePath(args[1])
gap <- normalizePath(args[4])
telom <- normalizePath(args[5])

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
df_new <- merge(df, coor, by = "chr", all.x = TRUE)
df_new <-  merge(df_new, tlm, by = "chr", all.x = TRUE)

# creating the table

df_final <- df_new %>%
  mutate(telomere = case_when(
    p > 50 & q > 50 ~ "pq",
    p > 50 ~ "p",
    q > 50 ~ "q",
    TRUE ~ "0"),
    t2t = case_when(
      telomere == "pq" & is.na(gapcount) ~ "T2T",
      telomere == "pq" & !is.na(gapcount) ~ "TgapT"
      telomere == "p" & is.na(gapcount) ~ "T2noT",
      telomere == "p" & !is.na(gapcount) ~ "TgapnoT",
      telomere == "q" & is.na(gapcount) ~ "noT2T",
      telomere == "q" & !is.na(gapcount) ~ "noTgapT",
      TRUE ~ "no")) %>%
  select(chr,
         completion = t2t,
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

write_tsv(df_final, file = paste0(dir,"/",dir,"_all_STATS.tsv"))
