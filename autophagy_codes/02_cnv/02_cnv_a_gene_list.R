library(magrittr)
tcga_path = "S:/study/生存分析/免疫检查点project/liucj_tcga_process_data"
expr_path <- "S:/study/生存分析/免疫检查点project/result"
expr_path_a <- file.path(expr_path, "all_expr")
cnv_path <- "S:/study/生存分析/免疫检查点project/result/5.cnv"

# load cnv and gene list
cnv <- readr::read_rds(file.path(tcga_path, "pancan34_cnv.rds.gz"))
<<<<<<< HEAD
gene_list_path <- "S:/study/生存分析/免疫检查点project/免疫检查点"
gene_list <- read.table(file.path(gene_list_path, "all.entrez_id-gene_id"),header=T)
gene_list$symbol %>% as.character() ->gene_list$symbol
gene_type<-read.table(file.path(gene_list_path,"checkpoint.type"),header=T)
gene_list<-dplyr::left_join(gene_list,gene_type,by="symbol")
=======
gene_list <- readr::read_rds(file.path(expr_path_a, "rds_03_a_atg_lys_gene_list.rds.gz"))
>>>>>>> c132521682f5114b521351d1257425ece5a06a56

filter_gene_list <- function(.x, gene_list) {
  gene_list %>%
    dplyr::select(symbol) %>%
    dplyr::left_join(.x, by = "symbol")
}

cnv %>%
  dplyr::mutate(filter_cnv = purrr::map(cnv, filter_gene_list, gene_list = gene_list)) %>%
  dplyr::select(-cnv) -> gene_list_cnv

readr::write_rds(x = gene_list_cnv, path = file.path(cnv_path, ".rds_02_cnv_a_gene_list.rds.gz"), compress = "gz")

fn_get_amplitue_threshold <- function(.x){
  ifelse(abs(.x) < log2(4 / 2) , 0, .x) -> .y
  tibble::tibble(a = sum(.y > 0) / length(.y), d = sum(.y < 0) / length(.y)) 
}
fn_get_ad <- function(.d){
  .d %>% 
    unlist(use.name = F) %>% 
    fn_get_amplitue_threshold()
}
fn_get_percent <- function(cancer_types, filter_cnv){
  filter_cnv %>%
    tidyr::nest(-symbol) %>% 
    dplyr::mutate(ad = purrr::map(data, .f = fn_get_ad)) %>% 
    dplyr::select(-data) %>% 
    tidyr::unnest(ad) %>% 
    tibble::add_column(cancer_types = cancer_types, .before = 1)
}

gene_list_cnv %>% dplyr::filter(cancer_types == "KIRC") %>% 
  dplyr::mutate(rs = purrr::map2(cancer_types, filter_cnv, fn_get_percent)) %>% 
  dplyr::select(-cancer_types, -filter_cnv) %>% 
  tidyr::unnest(rs)

cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
gene_list_cnv %>% 
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fn_get_amplitue_threshold", fn_get_amplitue_threshold)  %>%
  multidplyr::cluster_assign_value("fn_get_ad", fn_get_ad) %>% 
  multidplyr::cluster_assign_value("fn_get_percent", fn_get_percent) %>% 
  dplyr::mutate(rs = purrr::map2(cancer_types, filter_cnv, fn_get_percent)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) %>% 
  dplyr::select(-cancer_types, -filter_cnv) %>% 
  tidyr::unnest(rs) -> gene_list_cnv_per
on.exit(parallel::stopCluster(cluster))

library(ggplot2)
gene_list_cnv_per %>% 
  tidyr::drop_na() %>% 
  tidyr::gather(key = effect, value = per, a, d) %>% 
  dplyr::mutate(effect = plyr::revalue(effect, replace = c("a" = "Amplification", "d" = "Deletion"))) %>% 
  dplyr::mutate(per = ifelse(per > 0.6, 0.6, per)) -> plot_ready
plot_ready %>%
  dplyr::left_join(gene_list,by="symbol") %>%
  dplyr::mutate(color=plyr::revalue(functionWithImmune, replace = c('TwoSide' = "blue", "Inhibit" = "red", "Activate" = "green"))) %>% 
  dplyr::mutate(size = plyr::revalue(type,replace = c('Receptor'="bold.italic",'Ligand'="plain"))) %>%
  dplyr::select(symbol,color,size) %>%
  unique() ->gene_col_size
gene_col_size$size %>% as.character() ->gene_col_size$size
gene_col_size$color %>% as.character() ->gene_col_size$color

plot_ready %>% 
  ggplot(aes(y = symbol, x = cancer_types)) +
  geom_point(aes(size = per, color = effect)) +
  xlab("Cancer type") +
  ylab("Symbol") +
  scale_size_continuous(
    name = "CNV%",
    breaks = c(0.1, 0.2, 0.4, 0.6),
    limits = c(0.1, 0.6),
    labels = c("10", "20", "40", "60")
  ) +
  ggthemes::scale_color_gdocs(
    name = "SCNA Type"
  ) +
  theme(axis.text.y = element_text(color = gene_col_size$color,face = gene_col_size$size),
        axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
  facet_wrap(~ effect) +
  theme(strip.text.x=element_text(size = 15))-> p;p
ggsave(filename = "01_SCNV_all.pdf", plot = p, device = "pdf", path = cnv_path, width = 15, height = 10)




plot_ready %>% dplyr::filter(effect == "Amplification") -> a_ready
a_ready %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::summarise(s = sum(per)) %>% 
  dplyr::arrange(dplyr::desc(s)) -> cancer_rank
a_ready %>% 
  dplyr::filter(per > 0.1) %>% 
  dplyr::group_by(symbol) %>% 
  dplyr::summarise(s = n()) %>% 
  dplyr::filter(s >= 5) %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  # dplyr::filter(status %in% c("l")) %>% 
  dplyr::mutate(color=plyr::revalue(functionWithImmune, replace = c('TwoSide' = "blue", "Inhibit" = "red", "Activate" = "green"))) %>% 
  dplyr::mutate(size = plyr::revalue(type,replace = c('Receptor'="bold.italic",'Ligand'="plain"))) %>%
  # dplyr::filter(s > 1.8) %>% 
  dplyr::arrange(color,s) -> a_gene_rank
a_gene_rank$size %>% as.character() ->a_gene_rank$size
a_gene_rank$color %>% as.character() ->a_gene_rank$color

a_ready %>% 
  ggplot(aes(y = symbol, x = cancer_types)) +
  geom_point(aes(size = per, color = effect)) +
  xlab("Cancer type") +
  ylab("Symbol") +
  scale_size_continuous(
    name = "CNV percentage",
    breaks = c(0.1, 0.2, 0.4, 0.6),
    limits = c(0.1, 0.6),
    labels = c("10", "20", "40", "60")
  ) +
  scale_y_discrete(limit = a_gene_rank$symbol) +
  scale_x_discrete(limit = cancer_rank$cancer_types)+
  ggthemes::scale_color_gdocs(
    name = "SCNA Type"
  ) +
  ggthemes::theme_gdocs() +
  theme(axis.text.y = element_text(color = a_gene_rank$color,face = a_gene_rank$size),
        axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain")) -> p;p
ggsave(filename = "02_SCNV_amplification_seminar.pdf", plot = p, device = "pdf", path = cnv_path, width = 8, height = 6)


plot_ready %>% dplyr::filter(effect == "Deletion") -> d_ready
d_ready %>% 
  dplyr::filter(per >= 0.1) %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::summarise(s = sum(per)) %>% 
  dplyr::arrange(dplyr::desc(s)) -> cancer_rank
d_ready %>%   
  dplyr::filter(per >= 0.1) %>% 
  dplyr::group_by(symbol) %>% 
  dplyr::summarise(s = n()) %>% 
  dplyr::filter(s >= 5) %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  dplyr::mutate(color=plyr::revalue(functionWithImmune, replace = c('TwoSide' = "blue", "Inhibit" = "red", "Activate" = "green"))) %>% 
  dplyr::mutate(size = plyr::revalue(type,replace = c('Receptor'="bold.italic",'Ligand'="plain"))) %>%
  # dplyr::filter(s > 5) %>%
  dplyr::arrange(color, s) -> d_gene_rank
d_gene_rank$size %>% as.character() ->d_gene_rank$size
d_gene_rank$color %>% as.character() ->d_gene_rank$color

d_ready %>% 
  ggplot(aes(y = symbol, x = cancer_types)) +
  geom_point(aes(size = per, color = effect)) +
  xlab("Cancer type") +
  ylab("Symbol") +
  scale_size_continuous(
    name = "CNV percentage",
    breaks = c(0.1, 0.2, 0.4, 0.6),
    limits = c(0.1, 0.6),
    labels = c("10", "20", "40", "60")
  ) +
  scale_y_discrete(limit = d_gene_rank$symbol) +
  scale_x_discrete(limit = cancer_rank$cancer_types)+
  ggthemes::scale_color_gdocs(
    name = "SCNA Type"
  ) +
  ggthemes::theme_gdocs() +
  theme(axis.text.y = element_text(color = d_gene_rank$color,face = d_gene_rank$size),
        axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain")) -> p;p
ggsave(filename = "03_SCNV_deletion_seminar.pdf", plot = p, device = "pdf", path = cnv_path, width = 8, height = 6)


save.image(file = file.path(cnv_path, ".rda_02_cnv_a_gene_list.rda"))
load(file = file.path(cnv_path, ".rda_02_cnv_a_gene_list.rda"))

rm(list=ls())


