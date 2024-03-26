library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No directory path provided as an argument.")
}

dir <- args[1]

setwd(dir)
print(getwd())

contig_list <- "contig_list/"
archived <- "archived/"

# analyse paf file to get homologous chromosomes
print("Reading and preparing data to analyse and get homologous chromosome from contigs")
chrom_list = c(1:29, "X","Y")
minimap_paf <- read_delim("minimap.paf", "\t", col_names = FALSE)
names(minimap_paf) <- c("query_name","query_length","query_start","query_end","orientation",
                        "ref_name","ref_length","ref_start","ref_end","N_match","align_length","mapq","tp", "cm", "s1", "s2", "dv", "rl")

minimap_paf_mutate <- minimap_paf %>% mutate(query_align = query_end - query_start, proportion_ref = query_align/ref_length, proportion_query = query_align/query_length, proportion_match = N_match/align_length) %>%
  select(query_name,orientation,query_length,ref_name,ref_length,ref_start,ref_end,query_align,proportion_query,proportion_match,proportion_ref,query_start,query_end,align_length,N_match)

## for details of the contigs per chromosomes
print("Running analysis to get homologous chromosome from contigs: output in contig_list and archived")
group_func <- function(chrom_list) {
  df <- minimap_paf_mutate %>%
    filter(ref_name == {{chrom_list}}) %>% filter(proportion_match >= 0.1, query_align >= 5000, proportion_query >= 0.01) %>%
    group_by(query_name) %>% summarise(orientation = names(which.max(table(orientation))), min_ref_start = min(ref_start), max_ref_end = max(ref_end),
                                       sum_query_align = sum(query_align), sum_N_match = sum(N_match), sum_align_length = sum(align_length), min_query_start = min(query_start),
                                       max_query_end = max(query_end), query_length = mean(query_length), ref_length = mean(ref_length)) %>%
    mutate(new_proportion_ref = round((sum_query_align/ref_length*100), digits = 4), new_proportion_query = sum_query_align/query_length, new_proportion_match = sum_N_match/sum_align_length, queryref_proportion = round((query_length/ref_length*100), digits = 4)) %>%
    filter(new_proportion_match >= 0.1, sum_query_align >=1000000, new_proportion_query >= 0.4, new_proportion_ref >= 0.01) %>% arrange(min_ref_start)
  df2 <- df %>% select(query_name,orientation,new_proportion_ref,queryref_proportion)
  df_test <- df
  df_test$chr <- {{chrom_list}}
  assign(paste0("detailed", chrom_list), df_test, envir = .GlobalEnv) #saving data frame in the environment
  write_tsv(df, gsub(" ", "", paste(archived, "detailed", chrom_list, ".list")), col_names = TRUE) #saving as TSV in the path folder
  write_tsv(df2, paste0(contig_list, chrom_list, ".list"), col_names = FALSE)
}

invisible(lapply(chrom_list, group_func)) # loop the function

# all_chr <- rbindlist(mget(ls(pattern = "^detailed\\d+")))
all_chr <- rbindlist(mget(ls(pattern = "^detailed")))

print(all_chr$chr)

all_chr_summary <- all_chr %>% select(chr, query_name, query_length, ref_length, queryref_proportion) %>%
  group_by(chr) %>%
  summarise(
    contig = n(),
    query_name = paste(query_name, collapse = ","),
    query_length = sum(query_length),
    ref_length = max(ref_length),
    contig_align_proportion = paste(queryref_proportion, collapse = ","),
    queryref_proportion = sum(queryref_proportion),
    gaps = contig - 1,
    missing_contig = ifelse(queryref_proportion < 100, "likely", "no")
    )

write_tsv(all_chr_summary, file=("all_chr_assembly.tsv"))
