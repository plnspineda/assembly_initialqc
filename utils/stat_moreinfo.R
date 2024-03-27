library(readr)
library(dplyr)
library(stringr)
library(tidyr)

## This is version 2.1 where I changed the categories to three instances e.g. T2T, TgapT, noT2T, and so on
## add to the script the contigs into individual files
## 

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No input files provided as an argument.")
}

dir <- normalizePath(args[1])
asm_path <- normalizePath(args[2])
asm_map <- normalizePath(args[3])
gap <- normalizePath(args[4])
telom <- normalizePath(args[5])

setwd(dir)
print(getwd())

# reading files
df <- read_delim("all_chr_assembly.tsv")
coor <- read_delim(gap, col_names = c("chr","gapstart","gapend"))
path <- read_delim(asm_path)
map <- read_csv(asm_map, col_names = c("col1"))
telomere <- read_delim(telom, col_names = c("chr","start","end","value"))

## for gaps
coor <- coor %>% mutate(gaplen = gapend - gapstart +1) %>%
  group_by(chr) %>%
  summarise(gapcount = n(), gapsize = sum(gaplen))

## for node paths in the graph
map <- map %>% filter(str_detect(col1,"^path")) %>% 
  separate(col1, c("remark", "query_name", "name"), sep = " ")

## for telomeres
tlm <- telomere %>%
  group_by(chr) %>%
  slice(c(1:2, (n() - 1):n())) %>%
  summarise(p = sum(value[1:2]),
            q = sum(value[(n() - 1):n()]))

key <- merge(path, map, by = "name")
df_new <- merge(df, coor, by = "chr", all.x = TRUE)
df_new <- merge(df_new, key, by = "query_name")
df_new <-  merge(df_new, tlm, by = "chr", all.x = TRUE)

# creating the table

df_final <- df_new %>%
  mutate(telomere = case_when(
    p > 50 & q > 50 ~ "pq",
    p > 50 ~ "p",
    q > 50 ~ "q",
    TRUE ~ "0"),
    completion = case_when(
      telomere == "pq" & is.na(gapcount) ~ "T_2_T",
      telomere == "pq" & !is.na(gapcount) ~ "T_gap_T",
      telomere == "p" & is.na(gapcount) ~ "T_nogap_noT",
      telomere == "p" & !is.na(gapcount) ~ "T_gap_noT",
      telomere == "q" & is.na(gapcount) ~ "noT_nogap_T",
      telomere == "q" & !is.na(gapcount) ~ "noT_gap_T",
      is.na(telomere) & !is.na(gapcount) ~ "noT_nogap_noT",
      TRUE ~ "noT_gap_noT"),
      score = case_when(
        telomere == "pq" ~ "1",
        telomere == "p" ~ "p0.5",
        telomere == "q" ~ "q0.5",
        TRUE ~ "0"
      )) %>%
  select(chr,
         completion,
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
         q_count = q,
         path_name = name,
         nodes = path) %>%
  arrange(factor(chr, levels = c(1:29, "X", "Y")))

df_final[is.na(df_final)] <- 0

write_tsv(df_final, file = paste0(dir,"/","all_STATS.tsv"))
