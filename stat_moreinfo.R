library(readr)
library(dplyr)
library(stringr)
library(tidyr)

## This is a version where I changed the criteria for the telomere. 2024.07.22

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No input files provided as an argument.")
}

dir <- normalizePath(args[1])
asm_path <- normalizePath(args[2])
asm_map <- normalizePath(args[3])
gap <- normalizePath(args[4])
telom <- normalizePath(args[5])
telom_cutoff <- args[6]
telom_length <- normalizePath(args[7])

t_cut <- as.numeric(telom_cutoff)
print(t_cut)

setwd(dir)
print(getwd())

# reading files
df <- read_delim("all_chr_assembly.tsv")
coor <- read_delim(gap, col_names = c("chr","gapstart","gapend"))
path <- read_delim(asm_path)
map <- read_csv(asm_map, col_names = c("col1"))
telomere <- read_delim(telom, col_names = TRUE)
telomere <- telomere %>%
  separate(chr, into = c("chr", "region"), sep = ":")
names(telomere) <- c("chr", "region", "q_count", "p_count")
telomere_length <- read_delim(telom_length, col_names = c("chr","start","end"))

## for gaps
coor <- coor %>% mutate(gaplen = gapend - gapstart +1) %>%
  group_by(chr) %>%
  summarise(gapcount = n(), gapsize = sum(gaplen))

## for node paths in the graph
map <- map %>% filter(str_detect(col1,"^path")) %>% 
  separate(col1, c("remark", "query_name", "name"), sep = " ")

# ## for telomeres
# tlm <- telomere %>%
#   group_by(chr) %>%
#   slice(c(1:2, (n() - 1):n())) %>%
#   summarise(p = sum(value[1:2]),
#             q = sum(value[(n() - 1):n()]))

## for telomere length with telomere_analysis.sh
tlm_len <- telomere_length %>%
  mutate(p = ifelse(start < 1, end - start, NA),
         q = ifelse(start > 0, end - start, NA),
         p_start = ifelse(start < 1, start, NA),
         p_end = ifelse(start < 1, end, NA),
         q_start = ifelse(start > 0, start, NA),
         q_end = ifelse(start > 0, end, NA)) %>%
  group_by(chr) %>%
  summarise(p_length = sum(p, na.rm = TRUE),
            q_length = sum(q, na.rm = TRUE),
            p_start = first(p_start[!is.na(p_start)]),
            p_end = first(p_end[!is.na(p_end)]),
            q_start = first(q_start[!is.na(q_start)]),
            q_end = first(q_end[!is.na(q_end)])) %>%
  ungroup()

## for telomere count
tlm <- telomere %>%
  group_by(chr) %>%
  summarize(
    p_count = sum(p_count, na.rm = TRUE),
    q_count = sum(q_count, na.rm = TRUE)
  )

key <- merge(path, map, by = "name")
df_new <- merge(df, coor, by = "chr", all.x = TRUE)
df_new <- merge(df_new, key, by = "query_name", all.x = TRUE)
df_new <-  merge(df_new, tlm, by = "chr", all.x = TRUE)
df_new <- merge(df_new, tlm_len, by = "chr", all.x = TRUE)

# creating the table

df_final <- df_new %>%
  mutate(telomere = case_when(
    !is.na(p_start) & !is.na(q_start) ~ "pq",
    !is.na(p_start) ~ "p",
    !is.na(q_start) ~ "q",
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
         p_count,
         p_length,
         q_count,
         q_length,
         p_start,
         p_end,
         q_start,
         q_end,
         path_name = name,
         nodes = path) %>%
  arrange(factor(chr, levels = c(1:29, "X", "Y")))

df_final[is.na(df_final)] <- 0

dirbase <- basename(dir)
write_tsv(df_final, file = paste0(dir,"/",dirbase,"_all_STATS.tsv"))