library(ggplot2)
`%>%` <- magrittr::`%>%`


# Path
tcga_path = "S:/study/生存分析/免疫检查点project/liucj_tcga_process_data"
expr_path <- "S:/study/生存分析/免疫检查点project/result"
expr_path_a <- file.path(expr_path, "all_expr")
methy_path <- "S:/study/生存分析/免疫检查点project/result/8.methylation"
methy_box <- file.path(methy_path, "boxplot")


# load methylation and gene list
methy <- readr::read_rds(file.path(tcga_path, "pancan33_meth.rds.gz"))
gene_list_path <- "S:/study/生存分析/免疫检查点project/免疫检查点"
gene_list <- read.table(file.path(gene_list_path, "all.entrez_id-gene_id"),header=T)
gene_list$symbol %>% as.character() ->gene_list$symbol
gene_type<-read.table(file.path(gene_list_path,"checkpoint.type"),header=T)
gene_list<-dplyr::left_join(gene_list,gene_type,by="symbol")

# functions
filter_gene_list <- function(.x, gene_list) {
  gene_list %>%
    dplyr::select(symbol) %>%
    dplyr::left_join(.x, by = "symbol")
}

methy %>% 
  # dplyr::slice(2:5) %>%  # tidyr::unnest()
  dplyr::mutate(filter_methy = purrr::map(methy, filter_gene_list, gene_list = gene_list)) %>% 
  dplyr::select(-methy) -> gene_list_methy
#------------------------------------------------------------
fun_barcode <- function(.b){
  stringr::str_sub(
    string = .b,
    start = 1,
    end = 12
  )
}
fun_tn_type <- function(.b){
  type <- .b %>% 
    stringr::str_split(pattern = "-", simplify = T) %>% 
    .[, 4] %>% 
    stringr::str_sub(1, 2)
}
fun_boxplot <- function(fig_name, data, path = methy_box){
  data %>% 
    ggpubr::ggboxplot(x = 'type', y = 'meth', color = 'type', pallete = 'jco' ) +
    ggpubr::stat_compare_means(
      method = "t.test", 
      label.y = 1, 
      label.x = 1.2) + 
    ggthemes::theme_gdocs() +
    scale_color_manual(values = c("#DC3912", "#3366CC")) +
    labs(y = "B-value", x = "", title = fig_name) -> p
  ggsave(filename = paste(fig_name, "pdf", sep = "."), p, device = "pdf", path = path, width = 4, height = 3)
}
fun_compare <- function(.x, .y ){
  .x %>% 
    dplyr::mutate(gene = as.character(gene )) %>% 
    tidyr::gather(key = barcode, value = meth, -symbol, -gene) %>% 
    dplyr::mutate(type = fun_tn_type(barcode)) %>% 
    dplyr::filter(type %in% c("01", "11")) %>% 
    dplyr::mutate(barcode = fun_barcode(barcode)) %>% 
    dplyr::select(symbol, gene, barcode, meth, type) %>% 
    dplyr::mutate(type = dplyr::case_when(
      type == "01" ~ "Tumor",
      type == "11" ~ "Normal"
    )) %>% 
    dplyr::filter(!is.na(gene)) -> .d
  if(nrow(.d) < 20 || length(unique(.d$type)) != 2){return(tibble::tibble())
    print(paste(.y,":no normal/tumor sample."))}
  # at least 10 samples
  .d %>% 
    dplyr::select(barcode, type) %>% 
    dplyr::distinct() %>%
    dplyr::group_by(type) %>% 
    dplyr::count() %>% 
    dplyr::pull(n) -> sample_num
  if(any(sample_num < 10)){return(tibble::tibble())
    print(paste(.y,":has less than 10 samples."))}
  
  .d %>% 
    dplyr::group_by(symbol, gene) %>% 
    dplyr::do(
      broom::tidy(
        t.test(meth ~ type, data = .)
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>% 
    dplyr::filter(fdr < 0.05) %>%
    dplyr::mutate(
      direction = dplyr::case_when(
          estimate > 0 ~ 0, # normal high
          estimate < 0 ~ 1 # tumor high
        )) %>% 
    dplyr::select(symbol, gene, direction, p_val = p.value, fdr,estimate) -> .d_out
  
  # draw every pic, I have draw it, so annotate it to get data again.
  # .d %>%
  #   dplyr::semi_join(.d_out, by = c("symbol", "gene")) %>%
  #   # dplyr::filter(symbol %in% c("ATP6V0D1", "ATP6V0A4")) %>%
  #   tidyr::nest(-symbol, -gene) %>%
  #   dplyr::mutate(fig_name = paste(.y, symbol, sep = "_")) %>% 
  #   dplyr::select(fig_name, data) %>% 
  #   purrr::pwalk(.f = fun_boxplot, path = methy_box)
  
  return(.d_out)
}


# gene_list_methy %>% 
#   dplyr::slice(1) %>%
#   dplyr::mutate(methy_comparison = purrr::map(.x = filter_methy, .y = cancer_types, .f = fun_compare)) %>% 
#   dplyr::select(-filter_methy) %>% 
#   tidyr::unnest() -> gene_list_methy_fdr

cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
gene_list_methy %>% 
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_library("ggplot2") %>%
  multidplyr::cluster_assign_value("fun_barcode", fun_barcode)  %>%
  multidplyr::cluster_assign_value("fun_tn_type", fun_tn_type)  %>%
  multidplyr::cluster_assign_value("fun_boxplot", fun_boxplot)  %>%
  multidplyr::cluster_assign_value("fun_compare", fun_compare)  %>%
  multidplyr::cluster_assign_value("methy_box", methy_box)  %>%
  dplyr::mutate(methy_comparison = purrr::map(.x = filter_methy, .y = cancer_types, .f = fun_compare)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) %>% 
  dplyr::select(-filter_methy) %>% 
  tidyr::unnest() -> gene_list_methy_fdr
on.exit(parallel::stopCluster(cluster))

readr::write_rds(gene_list_methy_fdr, path = file.path(methy_path, ".rds_02_gene_list_methy_fdr.rds.gz"), compress = "gz")
readr::write_tsv(gene_list_methy_fdr, path = file.path(methy_path, "tsv_02_gene_list_methy_fdr.tsv"))

gene_list_methy_fdr %>% 
  dplyr::mutate(fdr = -log10(fdr)) %>% 
  dplyr::mutate(fdr = ifelse(fdr > 50, 50, fdr)) %>%
  dplyr::mutate(Direction = ifelse(direction > 0, "Up", "Down"))-> plot_ready

plot_ready %>%
  dplyr::mutate(dir_a = ifelse(direction < 1, -1, direction)) %>%
  dplyr::group_by(symbol) %>%
  dplyr::summarise(rank = sum(dir_a)) %>% 
  dplyr::left_join(gene_list, by = "symbol") %>%
  dplyr::mutate(color = plyr::revalue(functionWithImmune, replace = c('TwoSide' = "blue", "Inhibit" = "red", "Activate" = "green"))) %>% 
  dplyr::mutate(size = plyr::revalue(type,replace = c('Receptor'="bold.italic",'Ligand'="plain"))) %>%
  dplyr::arrange(rank) ->gene_rank

plot_ready %>%
  dplyr::mutate(dir_a = ifelse(direction < 1, -1, direction)) %>%
  dplyr::group_by(symbol) %>%
  dplyr::summarise(rank = sum(dir_a)) %>% 
  dplyr::left_join(gene_list, by = "symbol") %>%
  dplyr::mutate(color = plyr::revalue(functionWithImmune, replace = c('TwoSide' = "blue", "Inhibit" = "red", "Activate" = "green"))) %>% 
  dplyr::mutate(size = plyr::revalue(type,replace = c('Receptor'="bold.italic",'Ligand'="plain"))) %>%
  dplyr::arrange(rank) ->gene_rank
gene_rank$size %>% as.character() ->gene_rank$size
gene_rank$color %>% as.character()-> gene_rank$color
plot_ready %>% 
  dplyr::mutate(dir_a = ifelse(direction < 1, -1, direction)) %>%
  dplyr::group_by(cancer_types) %>% 
  dplyr::summarise(s = sum(dir_a)) %>% 
  dplyr::arrange(dplyr::desc(s)) -> cancer_rank

plot_ready %>% 
  ggplot(aes(x = cancer_types, y = symbol)) +
  geom_point(aes(size = fdr, color = Direction)) +
  scale_x_discrete(limit = cancer_rank$cancer_types) +
  scale_y_discrete(limit = gene_rank$symbol) +
  scale_size_continuous(name = "FDR") +
  theme(
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "grey",
      linetype = "dashed",
      size = 0.2
    ),
    
    axis.title = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 45,vjust=1,hjust = 1),
    axis.text.y = element_text(color = gene_rank$color,face = gene_rank$size),
    
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
#    legend.direction = "horizontal",
    legend.key = element_rect(fill = "white", colour = "black") 
  ) +
  ggthemes::scale_color_gdocs(name = "M(Tumor)\n—————\nM(Normal)") ->p
ggsave(filename = 'methy_direction.pdf', plot = p, device = "pdf", 
       path = methy_path, width = 8, height = 8)


#
save.image(file = file.path(methy_path, ".rda_02_methy_a_gene_list.rda"))
rm(list=ls())
load(file = file.path(methy_path, ".rda_02_methy_a_gene_list.rda"))

