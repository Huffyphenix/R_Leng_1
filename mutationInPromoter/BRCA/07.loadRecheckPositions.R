#!usr/bin/Rscript

###################
#Load library
##################
require(methods)
require(tidyverse)
require(stringr)

######################
# Load anno data
#####################
dataRootPath = "/home/cliu18/liucj/projects/1.Mutation_calling_in_non-condig_region_through_EXOME/3.calling/BRCA_reanalysis"

recheckDir = file.path(dataRootPath, "03.somatic/03.somaticForAnalysis.recheckPositions")

tumorDir = file.path(recheckDir, "tumor_tmp")
normalDir = file.path(recheckDir, "normal_tmp")
tumorList = list.files(path = tumorDir)
normalList = list.files(path = normalDir)

maniDir = file.path(dataRootPath, "03.somatic/03.somaticForAnalysis.saveFiles")
manifest <- read_rds(path = file.path(maniDir, "manifestData.info.rds")) %>%
  dplyr::select(barcode, type, bam, uuid)

forAnalysisDir <- file.path(dataRootPath, '03.somatic/03.somaticForAnalysis.saveFiles')
# Load somatic mutation with recur >= 3
somaticMutation <- 
  read_rds(path = file.path(forAnalysisDir, 
                            'somatic_mutation_recur3_for_analysis.rds')) %>%
  separate(mutation, into = c("mutation", "alt2"), sep = "\\/") 

##Mutation target
mutationTarget <- 
  read_rds(path = file.path(forAnalysisDir,
                            "nearest_target_coding_genes_3.rds")) %>% 
  as_tibble()

getPileupInfo <- function(x, path){
  print(x)
  read_tsv(file.path(path, x), 
           col_names = c("pileup","pos", "ref", "count", "base", "qual"),
           col_types = cols(ref = "c", count = 'c')) %>%
    separate(pileup, into = c("bam","chrom"), sep = ":") %>%
    mutate(bam = str_replace(bam, "/extraspace/TCGA/WXS_RAW/BRCA/regulatoryBam/",""),
           alt_pile = str_replace_all(base, pattern = "[\\,\\.\\$\\^\\-\\+\\:\\d\\@QP\"\\<\\]\\*\\!\\#\\/]","")) %>%
    mutate(bam = str_replace(bam, "tumor/", ""), alt_count = str_count(alt_pile)) %>%
    mutate(bam = str_replace(bam, "extracted\\.bam\\.mpileup","bam")) %>%
    unite(mutation, chrom, pos, sep=":") %>%
    unite(mutation, mutation, ref) %>%
    dplyr::select(-qual) %>%
    left_join(manifest, by = "bam")
}

tumorMutationBase <-
  tumorList %>%
  lapply(getPileupInfo, path = tumorDir) %>%
  bind_rows()

normalMutationBase <-
  normalList %>%
  lapply(getPileupInfo, path = normalDir) %>%
  bind_rows()

totalMutationBase <- 
  bind_rows(normalMutationBase, tumorMutationBase) 


totalMutationBase %>%
  write_rds(file.path(recheckDir, "recur3_totalMutationBase.rds"))

totalMutationBase %>%
  write_tsv(file.path(recheckDir, "recur3_totalMutationBase.tsv"))

# save the pileup info
write_tsv(totalMutationBase, path = file.path(dataRootPath, "03.somatic/03.somaticForAnalysis.saveFiles/recur3_totalMutationBase.tsv"))
write_rds(totalMutationBase, path = file.path(dataRootPath, "03.somatic/03.somaticForAnalysis.saveFiles/recur3_totalMutationBase.rds"))

##################################
# Calculate the refined recurrent#
##################################
# In fact, the mpileup is a way samtools call mutation.
# The mutation point import from mpileup output is a validation for the GATK output.

tumor_mutation_base_count_10 <- 
  tumorMutationBase %>% 
  filter(as.integer(count) >= 10) %>%
  group_by(mutation) 

# this is the false positive
# anti_join false positive mutation in normal
filter_out_mutation <- 
  totalMutationBase %>%
  filter(as.integer(count) >= 10, alt_count >= 3, type == "normal") %>%
  group_by(mutation) %>%
  count() %>%
  filter(n >=3)

mutation_feature <- 
  somaticMutation %>%
  anti_join(filter_out_mutation, by = "mutation") %>%
  dplyr::select(mutation, alt, feature, ensr) %>%
  distinct()

label_observation_count <- function(x){
  return(c(y = 180, label = length(x)))
}

totalMutationBase %>% 
  drop_na() %>% 
  mutate(
    count = as.integer(count)
  ) %>% 
  drop_na() ->
  totalMutationBase
# Figure 
depth_for_tumor_normal_barplot <- 
  totalMutationBase %>%
  anti_join(filter_out_mutation, by = "mutation") %>%
  left_join(mutation_feature, by = "mutation") %>%
  unite(mutation, mutation, alt, sep = "/") %>%
  filter(as.integer(count) >= 10) %>%
  ggplot(aes(type, as.integer(count))) +
  geom_boxplot(aes(color = type), outlier.shape  = NA) +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 200)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  stat_summary(fun.data = label_observation_count, geom = "text") +
  labs(title = "SAMtools mpileup check for candidate mutations",
       x = "Type", y = "Depth") +
  facet_wrap(~mutation)
# print(depth_for_tumor_normal_barplot)
ggsave(filename = file.path(recheckDir, 
                            "01.recur3_depth_for_tumor_normal_barplot.png"),
       device = "png", width = 20, height = 20)
depth_for_tumor_normal_barplot %>%
  write_rds(path = file.path(recheckDir, 
                             "01.recur3_depth_for_tumor_normal_barplot.rds"))

#######################
#filter false positive#
#######################
# Narrow down candidate mutation
# The filter criteria
# Mutation point has count >= 10 and alt_count >=3 in normal


# filter out normal false mtuation
# get depth >= 10 normal sample 
totalMutationBase %>% 
  anti_join(filter_out_mutation, by = "mutation") %>% 
  filter(type == "tumor", count >= 10) %>%
  tbl_df() %>%
  dplyr::select(mutation, depth = count, alt_depth = alt_count, barcode) %>%
  ggplot(aes(x = reorder(mutation, mutation, function(x) -length(x)))) +
  geom_bar() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Mutation", y = "# of Samples")
  

mpileup_mutation <-
  totalMutationBase %>% 
  anti_join(filter_out_mutation, by = "mutation") %>% 
  filter(type == "tumor", count >= 10) %>%
  tbl_df() %>%
  dplyr::select(mutation, depth = count,alt_depth = alt_count, barcode)

# Union the mpieup result
# For the same mutation and sample with different depth and alt_depth, I will chose GATK result.


merged_mutation <- 
  somaticMutation %>% 
  anti_join(filter_out_mutation, by = "mutation") %>%
  dplyr::select(mutation, depth, alt_depth = altdepth, barcode) %>%
  dplyr::union(mpileup_mutation) %>%
  mutate(type = ifelse(alt_depth >= 3, "MT", "WT")) %>%
  left_join(mutation_feature, by = "mutation") %>%
  unite(mutation, mutation, alt, sep = "/") %>%
  arrange(type) %>%
  distinct(mutation, barcode, .keep_all = T)

merged_mutation %>% 
  write_tsv(path=file.path(forAnalysisDir, "recur3_candidate_merged_mutation.tsv"))
merged_mutation %>%
  write_rds(path = file.path(forAnalysisDir, "recur3_candidate_merged_mutation.rds"))

mutated_samples_barplot <- 
  merged_mutation %>%
  group_by(mutation, type) %>%
  filter(type == "MT") %>%
  ggplot(aes(x = reorder(mutation, mutation, function(x) - length(x)))) +
  geom_bar(fill = "royalblue3") +
  theme_minimal() +
  geom_text(stat = 'count', aes(label = ..count..), 
            vjust = -0.2, hjust = 0.5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Mutation", y = "# of Samples", title = "Distribution of mutated samples") 
print(mutated_samples_barplot)
ggsave(filename = file.path(forAnalysisDir, 
                            "02.mutated_samples_barplot.png"),
       device = "png", width = 20, height = 20)
mutated_samples_barplot %>%
  write_rds(path = file.path(forAnalysisDir, 
                             "02.mutated_samples_barplot.rds"))

mutation_reserved <- 
  merged_mutation %>% distinct(mutation)

label_observation_count <- function(x){
  return(c(y = 180, label = length(x)))
}


candidate_mutation_coverage <- 
  totalMutationBase %>%
  left_join(mutation_feature, by = "mutation") %>%
  unite(mutation, mutation, alt, sep = "/") %>%
  semi_join(mutation_reserved,by = "mutation") %>% 
  filter(count >= 10) %>%
  ggplot(aes(type, count)) +
  geom_boxplot(aes(color = type), outlier.shape = NA) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  stat_summary(fun.data = label_observation_count, geom = "text") +
  labs(title = "SAMtools mpileup check for candidate mutations",
       x = "Type", y = "Depth") +
  scale_y_continuous(limits = c(0, 200)) +
  facet_wrap(~mutation)
print(candidate_mutation_coverage)

candidate_mutation_coverage %>%
  ggsave(filename = file.path(forAnalysisDir, 
                            "03.recur3_candidate_mutation_coverage.png"),
       device = "png", width = 20, height = 20)
candidate_mutation_coverage %>%
  write_rds(path = file.path(forAnalysisDir, 
                             "03.candidate_mutation_coverage.rds"))



# save workspace
save(list = ls(), file = file.path(recheckDir, "07.recur3_loadRecheckPositions.RData"))
load(file.path(recheckDir, "07.recur3_loadRecheckPositions.RData"))
