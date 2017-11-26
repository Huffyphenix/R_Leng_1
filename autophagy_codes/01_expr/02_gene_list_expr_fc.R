#load(file = file.path(out_path, "rda_00_gene_expr.rda"))
library(magrittr)
library(dplyr)
library(RDS)
library(gtools)

# pac_check<-function(x){
#   if(require(x)){
#   print(paste(x,"is loaded correctly."))} else {
#   print(paste("trying to install",x))
#   install.packages(x)
#   if(require(x)){
#     print(paste(x,"installed and loaded."))
#   } else {
#     stop(paste("could not install",x))
#   }
#   }
# }
#pac_check("dplyr")
# processed path
tcga_path = "S:/study/生存分析/免疫检查点project/liucj_tcga_process_data"
expr <- readr::read_rds(file.path(tcga_path, "pancan33_expr.rds.gz",fsep=.Platform$file.sep))

# Read gene list
# Gene list was compress as rds
gene_list_path <- "/home/huff/immune_checkpoint/checkpoint"
#gene_list <- readr::read_rds(file.path(gene_list_path, "all.entrez_id"))
gene_list <- read.table(file.path(gene_list_path, "all.entrez_id-gene_id"),header=T)
gene_list$symbol %>% as.character() ->gene_list$symbol
gene_type<-read.table(file.path(gene_list_path,"checkpoint.type"),header=T)
gene_list<-dplyr::left_join(gene_list,gene_type,by="symbol")
#output path
out_path <- "S:/study/生存分析/免疫检查点project/result/1.DE"

#######################
# filter out genes
#######################
filter_gene_list <- function(.x, gene_list) {
  gene_list %>%
    dplyr::select(symbol) %>%
    dplyr::left_join(.x, by = "symbol")
}
expr %>%
  dplyr::mutate(filter_expr = purrr::map(expr, filter_gene_list, gene_list = gene_list)) %>%
  dplyr::select(-expr) -> gene_list_expr

#################################
# Caculate mean expr of genes in each cancer.
##################################
calculate_mean_expr <- function(.x, .y) {
  .y %>%
    tibble::add_column(cancer_types = .x, .before = 1) -> df
  
  # get cancer types and get # of smaple >= 10
  samples <-
    tibble::tibble(barcode = colnames(df)[-c(1:3)]) %>%
    dplyr::mutate(
      sample = stringr::str_sub(
        string = barcode,
        start = 1,
        end = 12
      ),
      type = stringr::str_split(barcode, pattern = "-", simplify = T)[, 4] %>% stringr::str_sub(1, 2)
    ) %>%
    dplyr::filter(type %in% c("01", "11")) %>%
    dplyr::mutate(type = plyr::revalue(
      x = type,
      replace = c("01" = "Tumor", "11" = "Normal"),
      warn_missing = F
    )) %>%
    dplyr::group_by(sample) %>%
    dplyr::filter(n() >= 2, length(unique(type)) == 2) %>%
    dplyr::ungroup()
  sample_type_summary <- table(samples$type) %>% as.numeric()
  if (gtools::invalid(sample_type_summary) ||
      any(sample_type_summary < c(10, 10))) {
    return(NULL)
  }
  
  # filter out cancer normal pairs
  df_f <-
    df %>%
    dplyr::select(c(1, 2, 3), samples$barcode) %>%
    tidyr::gather(key = barcode, value = expr, -c(1, 2, 3)) %>%
    dplyr::left_join(samples, by = "barcode")
  #mean exp
  df_f %>%
    dplyr::group_by(cancer_types, symbol, type) %>%
    tidyr::drop_na(expr) %>%
    dplyr::summarise(mean=mean(expr))%>%
    dplyr::ungroup() ->df_mean_exp
  return(df_mean_exp)
}

purrr::map2(.x = gene_list_expr$cancer_types,
            .y = gene_list_expr$filter_expr,
            .f = calculate_mean_expr) -> gene_mean_exp
names(gene_mean_exp) <- gene_list_expr$cancer_types
gene_mean_exp %>% dplyr::bind_rows() -> gene_mean_exp_simplified
readr::write_rds(
  x = gene_mean_exp_simplified,
  path = file.path(out_path, "rds_01_gene_list_mean_exp_simplified.rds.gz"),
  compress = "gz"
)
readr::write_tsv(
  x = gene_mean_exp_simplified,
  path = file.path(out_path, "tsv_01_gene_list_mean_exp_simplified.tsv")
)
###############
#Draw pictures
################
gene_mean_exp_simplified %>%
  dplyr::mutate(type = as.factor(type)) %>% 
  dplyr::mutate(type = plyr::revalue(type, replace = c("Normal" = "N", "Tumor" = "T"))) %>%
  tidyr::unite(col = symbol, symbol, type) ->gene_mean_exp_simplified_plot
#heat map
ggplot(gene_mean_exp_simplified_plot,
       aes(x = symbol, y = cancer_types, fill = log2(mean))) +
  geom_tile(color = "black") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "white",
    breaks = seq(-5, 15, length.out = 5),
    #labels = c("<= -3", "-1.5", "0", "1.5", ">= 3"),
    name = "log2(MeanExp)"
  ) +
  #scale_y_discrete(limit = gene_rank$symbol) +
  scale_y_discrete(limit = cancer_types_rank$cancer_types, expand = c(0, 0)) +
  theme(
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
   # axis.text.y = element_text(color = gene_rank$color,face=gene_rank$size),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text.x = element_text(angle = 270, hjust = 0,vjust=0.5)
  ) -> p;p
ggsave(
  filename = "fig_01_mean_expr_pattern.pdf",
  plot = p,
  device = "pdf",
  width = 15,
  height = 6,
  path = out_path
)
readr::write_rds(
  p,
  path = file.path(out_path, "fig_01_mean_expr_pattern.pdf.rds.gz"),
  compress = "gz"
)
#################################
# Caculate p-value and fold-change.
##################################
calculate_fc_pvalue <- function(.x, .y) {
  .y %>%
    tibble::add_column(cancer_types = .x, .before = 1) -> df
  
  # get cancer types and get # of smaple >= 10
  samples <-
    tibble::tibble(barcode = colnames(df)[-c(1:3)]) %>%
    dplyr::mutate(
      sample = stringr::str_sub(
        string = barcode,
        start = 1,
        end = 12
      ),
      type = stringr::str_split(barcode, pattern = "-", simplify = T)[, 4] %>% stringr::str_sub(1, 2)
    ) %>%
    dplyr::filter(type %in% c("01", "11")) %>%
    dplyr::mutate(type = plyr::revalue(
      x = type,
      replace = c("01" = "Tumor", "11" = "Normal"),
      warn_missing = F
    )) %>%
    dplyr::group_by(sample) %>%
    dplyr::filter(n() >= 2, length(unique(type)) == 2) %>%
    dplyr::ungroup()
  sample_type_summary <- table(samples$type) %>% as.numeric()
  if (gtools::invalid(sample_type_summary) ||
      any(sample_type_summary < c(10, 10))) {
    return(NULL)
  }
  
  # filter out cancer normal pairs
  df_f <-
    df %>%
    dplyr::select(c(1, 2, 3), samples$barcode) %>%
    tidyr::gather(key = barcode, value = expr, -c(1, 2, 3)) %>%
    dplyr::left_join(samples, by = "barcode")
  
  # pvalue & fdr
  df_f %>%
    dplyr::group_by(cancer_types, symbol, entrez_id) %>%
    tidyr::drop_na(expr) %>%
    dplyr::do(broom::tidy(t.test(expr ~ type, data = .))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
    dplyr::select(cancer_types, symbol, entrez_id, p.value, fdr) -> df_pvalue
  
  # log2 fold change mean
  df_f %>%
    dplyr::group_by(cancer_types, symbol, entrez_id, type) %>%
    tidyr::drop_na(expr) %>%
    dplyr::summarise(mean = mean(expr)) %>%
    tidyr::spread(key = type, mean) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(fc = (Tumor + 0.1) / (Normal + 0.1)) -> df_fc
  
  df_fc %>%
    dplyr::inner_join(df_pvalue, by = c("cancer_types", "symbol", "entrez_id")) %>%
    dplyr::mutate(n_normal = sample_type_summary[1], n_tumor = sample_type_summary[2]) -> res
  return(res)
}

purrr::map2(.x = gene_list_expr$cancer_types,
            .y = gene_list_expr$filter_expr,
            .f = calculate_fc_pvalue) -> gene_list_fc_pvalue
names(gene_list_fc_pvalue) <- gene_list_expr$cancer_types

gene_list_fc_pvalue %>% dplyr::bind_rows() -> gene_list_fc_pvalue_simplified
readr::write_rds(
  x = gene_list_fc_pvalue_simplified,
  path = file.path(out_path, "rds_01_gene_list_fc_pvalue_simplified.rds.gz"),
  compress = "gz"
)
readr::write_tsv(
  x = gene_list_fc_pvalue_simplified,
  path = file.path(out_path, "tsv_01_gene_list_fc_pvalue_simplified.tsv")
)

# write sample pairs number
gene_list_fc_pvalue_simplified %>%
  dplyr::select(cancer_types, n_normal, n_tumor) %>%
  dplyr::distinct() -> pancan_samples_pairs
readr::write_rds(
  x = pancan_samples_pairs,
  path = file.path(out_path, "rds_02_pancan_samples_pairs.rds.gz"),
  compress = "gz"
)
readr::write_tsv(x = pancan_samples_pairs,
                 path = file.path(out_path, "tsv_02_pancan_samples_pairs.tsv"))

###############
#Draw pictures
################
gene_list_fc_pvalue_simplified %>%
  # filter |log2(fc)| >= log2(1.5), fdr <= 0.05
  dplyr::filter(abs(log2(fc)) >= log2(1.5), fdr <= 0.05) %>%
  dplyr::mutate(p.value = -log10(p.value)) -> gene_list_fc_pvalue_simplified_filter

# expression pattern
# significant high expression is 1
# significant low expression is -1
# not significant is 0
get_pattern <- function(fc, p.value) {
  if ((fc > 1.5) && (p.value < 0.05)) {
    return(1)
  } else if ((fc < 2 / 3) && (p.value < 0.05)) {
    return(-1)
  } else {
    return(0)
  }
}

gene_list_fc_pvalue_simplified %>%
  dplyr::mutate(expr_pattern = purrr::map2_dbl(fc, p.value, get_pattern)) %>%
  dplyr::select(cancer_types, symbol, expr_pattern) %>%
  tidyr::spread(key = cancer_types, value = expr_pattern) ->
  gene_expr_pattern

gene_expr_pattern %>%
  dplyr::rowwise() %>%
  dplyr::do(
    symbol = .$symbol,
    rank =  unlist(.[-1], use.names = F) %>% sum(),
    up = (unlist(.[-1], use.names = F) == 1) %>% sum(),
    down = (unlist(.[-1], use.names = F) == -1) %>% sum()
  ) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() %>%
  dplyr::arrange(rank) %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  dplyr::mutate(color = plyr::revalue(x = functionWithImmune, replace = c("Inhibit" = "red", "Activate" = "green", "TwoSide"="blue"))) %>% 
  dplyr::mutate(size = plyr::revalue(type,replace = c('Receptor'="bold.italic",'Ligand'="plain"))) %>%
  dplyr::select(symbol, rank, up, down, type, color,size) -> gene_rank


gene_expr_pattern %>%
  dplyr::summarise_if(.predicate = is.numeric, dplyr::funs(sum(.))) %>%
  tidyr::gather(key = cancer_types, value = rank) %>%
  dplyr::arrange(-rank) -> cancer_types_rank

library(ggplot2)
gene_rank$color %>% as.character() ->gene_rank$color
gene_rank$size %>% as.character() ->gene_rank$size
#####################
#heat map
####################
ggplot(gene_list_fc_pvalue_simplified_filter,
       aes(x = cancer_types, y = symbol, fill = log2(fc))) +
  geom_tile(color = "black") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "white",
    breaks = seq(-3, 3, length.out = 5),
    labels = c("<= -3", "-1.5", "0", "1.5", ">= 3"),
    name = "log2(FC)"
  ) +
  scale_y_discrete(limit = gene_rank$symbol) +
  scale_x_discrete(limit = cancer_types_rank$cancer_types, expand = c(0, 0)) +
  theme(
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(color = gene_rank$color,face=gene_rank$size),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text.x = element_text(angle = 315, hjust = 0,vjust=1)
  ) -> p;p
ggsave(
  filename = "fig_01_expr_pattern.pdf",
  plot = p,
  device = "pdf",
  width = 6,
  height = 6,
  path = out_path
)
readr::write_rds(
  p,
  path = file.path(out_path, "fig_01_expr_pattern.pdf.rds.gz"),
  compress = "gz"
)
################
#receptor and ligand heatmap
################
gene_list_fc_pvalue_simplified_filter %>%
  dplyr::filter(symbol %in% c("TNFRSF4","TNFSF4","TNFRSF9","TNFSF9","CTLA4","CD80","CD86")) %>%
  ggplot(aes(x = cancer_types, y = symbol, fill = log2(fc))) +
  geom_tile(color = "black") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "white",
    breaks = seq(-3, 3, length.out = 5),
    labels = c("<= -3", "-1.5", "0", "1.5", ">= 3"),
    name = "log2(FC)"
  ) +
  scale_y_discrete(limit = c("TNFRSF4","TNFSF4","TNFRSF9","TNFSF9","CTLA4","CD80","CD86")) +
  scale_x_discrete(limit = cancer_types_rank$cancer_types, expand = c(0, 0)) +
  theme(
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(color = c("green","green","green","green","red","blue","blue"),
                               face=c("bold.italic","plain","bold.italic","plain","bold.italic","plain","plain")),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text.x = element_text(angle = 315, hjust = 0,vjust=1)
  ) -> p;p

ggsave(
  filename = "fig_02_ligand-receptor_heatmap.pdf",
  plot = p,
  device = "pdf",
  width = 6,
  height = 3,
  path = out_path
)
readr::write_rds(
  p,
  path = file.path(out_path, "fig_02_expr_pattern_fc_pval.pdf.rds.gz"),
  compress = "gz"
)


###########################
#point and heat map
###########################
ggplot(
  dplyr::mutate(
    gene_list_fc_pvalue_simplified_filter,
    alt = ifelse(log2(fc) > 0,  "up", "down")),
  aes(x = symbol, fill = factor(alt))) +
  geom_bar(color = NA, width = 0.5) +
  scale_fill_manual(
    limit = c("down", "up"),
    values = c("royalblue1", "tomato1"),
    guide = FALSE
  ) +
  scale_y_continuous(
    limit = c(-0.1, 12.5),
    expand = c(0, 0),
    breaks = seq(0, 12, length.out = 5)
  ) +
  scale_x_discrete(limit = gene_rank$symbol, expand = c(0.01, 0.01)) +
  theme(
    panel.background = element_rect(
      colour = "black",
      fill = "white",
      size = 1
    ),
    panel.grid.major = element_line(linetype = "dashed", color = "lightgray"),
    axis.title = element_blank(),
    axis.text.y = element_text(color = gene_rank$color,face=gene_rank$size,size = 16),
    axis.text.x = element_text(size = 16),
    axis.ticks.x = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key = element_rect(fill = "white", colour = "black")
  ) +
  coord_flip() -> p;p
ggsave(
  filename = "fig_03_expr_pattern_cancer_counts.pdf",
  plot = p,
  device = "pdf",
  width = 15,
  height = 15,
  path = out_path
)
readr::write_rds(
  p,
  path = file.path(out_path, "fig_03_expr_pattern_cancer_counts.pdf.gz"),
  compress = "gz"
)

save.image(file = file.path(out_path, "rda_00_gene_expr.rda"))
rm(list=ls())

