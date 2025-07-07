library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)

# Repeat masker out file as input needs to be cleaned first.
## tail -n +4 17.fasta.out | tr -s ' ' '\t' > tmp && sed 's/^[ \t]*//' tmp > 17.fasta.out.tsv

args = commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  stop("Usage: Rscript find_centromere.R [output_directory] [input_repeat_file] [output_prefix]")
}

dir <- normalizePath(args[1])
infile <- normalizePath(args[2])
basefile <- basename(infile)
inname <- gsub("\\.fasta\\.out$", "", basefile)

setwd(dir)

# Clean the repeat bed file
system(paste("tail -n +4", infile, "| tr -s ' ' '\t' > tmp"))
system(paste("sed 's/^[ \t]*//' tmp > tmp2"))
tsvfile <- paste(basefile, ".tsv", sep="")
new_dir <- file.path(dir, "ctr")
dir.create(new_dir)
setwd(new_dir)
system(paste0("mv ", dir, "/tmp2 ", new_dir, "/", tsvfile))
system(paste0("rm ", dir, "/tmp"))

df <- read_delim(paste0(basefile,".tsv"), 
                 col_names = c("score", "divergence", "deletion", "insertion", "query_sequence", "query_begin", 
                               "query_end", "query_left", "strand","repeat", "family", "repeat_begin", 
                               "repeat_end","repeat_left","ID", "remarks"))

df_filter <- df %>% filter(family == "Satellite/centr")

head(df_filter)
chr <- as.character(1:29)
chrX <- "X"
chrY <- "Y"
chr1_X <- c(chr,chrX,chrY)

df_filter$query_sequence[!df_filter$query_sequence %in% chr1_X] <- "Unplaced"

print("Step1")
#repeat_begin and repeat_left have silly ()
# Clean and convert repeat_begin and repeat_left columns
df_filter$new_repeat_begin <- as.integer(gsub("[()]", "", df_filter$repeat_begin))
df_filter$new_repeat_end <- as.integer(gsub("[()]", "", df_filter$repeat_end))
df_filter$new_repeat_left <- as.integer(gsub("[()]", "", df_filter$repeat_left))

print("Step2")
#loop thro and add query_seq_align_len and perc_id
infile_subset <- df_filter %>% 
  mutate(query_align_len = query_end - query_begin + 1, perc_id = 1 - (divergence/100)) %>%
  filter(perc_id > 0.6)
infile_subset$assembly <- rep(inname,nrow(infile_subset))

print("Step3")
#top align length of repeat family with highest presence
infile_subset_topRepeatLen <- infile_subset %>% 
  group_by(family) %>% summarise(total_len = sum(query_align_len)) %>% 
  arrange(desc(total_len))

print("Step4")
### use l-apply to all the autosomes
chrom_list <- c(1:29,"X","Y")
group_func <- function(chrom) {
  df <- infile_subset %>%
    filter(family == "Satellite/centr" & query_sequence == chrom) %>%
    mutate(query_align_len = query_end - query_begin + 1,
           perc_id = 1 - (divergence / 100),
           diff = query_end - lag(query_begin)) %>%
    filter(perc_id > 0.6) %>%
    select("query_sequence", "query_begin", "query_end", "query_left", "strand", "repeat", "family", "query_align_len", "diff")
  
  df$group <- cumsum(!is.na(df$diff) & df$diff > 1000000)
  
  df_summary <- df %>%
    group_by(group) %>%
    summarise(chromosome = chrom,
              start = first(query_begin),
              end = last(query_end),
              rep_size = sum(query_align_len),
              loc_size = last(query_end) - first(query_begin) + 1,
              asm = inname)
  
  assign(paste("HS_ctrloc_", chrom, sep = ""), df_summary, envir = .GlobalEnv)
}

invisible(lapply(chrom_list, group_func))

print("Step5")
infile_ctr <- do.call(rbind, lapply(chrom_list, function(chrom) get(paste("HS_ctrloc_", chrom, sep = ""))))
infile_ctr_final <- infile_ctr %>%
  group_by(chromosome) %>%
  top_n(1, rep_size)
infile_plot <- ggplot(infile_ctr_final, aes(x = factor(chromosome, levels = c(1:29, "X", "Y")), y = rep_size)) +
  geom_col(position = "identity", alpha = 0.3) +
  labs(title = "Size of centromeres per chromosomes", x = "Chromosome", y = "Size (bp)")
print(infile_plot)

print("Step6")
## plot the centromeres
infile_plot <- ggplot(infile_ctr_final, aes(x = factor(chromosome, levels = c(1:29, "X", "Y")))) +
  geom_col(aes(y = rep_size, fill = "Satellite repeats"), position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_col(aes(y = loc_size, fill = "Centromere length"), position = position_dodge(width = 0.8), alpha = 0.7) +
  labs(title = "Size of centromeres per chromosomes", x = "Chromosome", y = "Values") +
  scale_fill_manual(name = "Category", values = c("Satellite repeats" = "blue", "Centromere length" = "gray")) +
  theme(plot.title = element_text(hjust = 0.5))
png(paste0(inname,"_ctr_count.png"), width = 10, height = 8, units = "in", res = 300)
print(infile_plot)
dev.off()

ctr_type_plot <- infile_subset %>% 
  ggplot(aes(x = factor(query_sequence, levels = c(1:29, "X", "Y")), y = query_align_len, fill = `repeat`)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(title = "Centromere Length per Chromosomes", x = "Chromosomes", y = "Size", fill = "Repeat Families") +
  theme(plot.title = element_text(hjust = 0.5))
png(paste0(inname,"_ctr_reptype.png"), width = 10, height = 6, units = "in", res = 300)
print(ctr_type_plot)
dev.off()

print("Step7")
repeat_classes <- infile_subset %>% select(`repeat`) %>% group_by(`repeat`) %>% summarise(count=n())
write_tsv(repeat_classes, file = paste0(inname,"_ctr_repeatclasses.txt"))

infile_ctr_final <- infile_ctr_final %>% mutate(perc = rep_size/loc_size*100)
write_tsv(infile_ctr_final, file = paste0(inname,"_ctr_allchr.txt"))

write_tsv(infile_ctr, file = paste0(inname,"_candidate_ctr.txt"))

repeat_classes <- infile_subset %>% select(`repeat`) %>% group_by(`repeat`) %>% summarise(count=n())
write_tsv(repeat_classes, file = paste0(inname,"_ctr_repeatclasses.txt"))

infile_ctr_final <- infile_ctr_final %>% mutate(perc = rep_size/loc_size*100)
write_tsv(infile_ctr_final, file = paste0(inname,"_ctr_allchr.txt"))

write_tsv(infile_ctr, file = paste0(inname,"_candidate_ctr.txt"))
