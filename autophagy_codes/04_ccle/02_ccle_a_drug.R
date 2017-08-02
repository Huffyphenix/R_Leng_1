library(magrittr)
tcga_path = "S:/study/生存分析/免疫检查点project/liucj_tcga_process_data"
expr_path <- "S:/study/生存分析/免疫检查点project/result"
expr_path_a <- file.path(expr_path, "all_expr")
ccle_path <- "S:/study/生存分析/免疫检查点project/result/7.ccle"

# load data
drug_gdsc <- readr::read_rds(file.path(tcga_path, "drug_gdsc_exp_spearman.rds.gz"))
drug_ctrp <- readr::read_rds(file.path(tcga_path, "drug_ctrp_exp_spearman.rds.gz"))
gene_list_path <- "S:/study/生存分析/免疫检查点project/免疫检查点"
gene_list <- read.table(file.path(gene_list_path, "all.entrez_id-gene_id"),header=T)
gene_list$symbol %>% as.character() ->gene_list$symbol
gene_type<-read.table(file.path(gene_list_path,"checkpoint.type"),header=T)
gene_list<-dplyr::left_join(gene_list,gene_type,by="symbol")

t_gdsc <- readr::read_rds(file.path(tcga_path, "drug_target_gdsc.rds.gz")) %>% 
  tidyr::unnest() %>% 
  dplyr::select(drug_name, target_pathway) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(target_pathway) %>% 
  dplyr::mutate(count = n()) %>% 
  dplyr::ungroup()
t_ctrp <- readr::read_rds(file.path(tcga_path, "drug_target_ctrp.rds.gz")) %>% 
  tidyr::unnest() %>% 
  dplyr::select(drug_name, target_pathway) %>% 
  dplyr::distinct()

# GDSC analysis
gene_list %>% 
  dplyr::inner_join(drug_gdsc, by = "symbol") -> gdsc_gene_list

fn_filter_drug <- function(.x, .cor = 0.25, .fdr = 0.05) {
  .x %>% dplyr::filter(abs(cor_sprm) > .cor, fdr < .fdr)
}

gdsc_gene_list %>% 
  # dplyr::filter(symbol == "ATG7") %>% 
  # .$drug %>% .[[1]] -> .x
  dplyr::mutate(cor_drug = purrr::map(.x = drug, .f = fn_filter_drug)) %>% 
  tidyr::unnest(cor_drug) -> gdsc_gene_list_sig_drug

# gene rank
library(ggplot2)
gdsc_gene_list_sig_drug %>% 
  dplyr::mutate(
    fdr = ifelse(-log10(fdr) > 40, 40, -log10(fdr)),
    cor_sprm = ifelse(cor_sprm > 0.4, 0.4, cor_sprm),
    cor_sprm = ifelse(cor_sprm < -0.4, -0.4, cor_sprm)
  ) %>% 
  dplyr::left_join(t_gdsc, by = "drug_name") -> gdsc_plot_ready

gdsc_plot_ready %>% 
  dplyr::group_by(symbol) %>% 
  dplyr::summarise(cor_sum = sum(cor_sprm)) %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  dplyr::mutate(color = plyr::revalue(functionWithImmune, replace = c('TwoSide' = "blue", "Inhibit" = "red", "Activate" = "green"))) %>% 
  dplyr::mutate(size = plyr::revalue(type,replace = c('Receptor'="bold.italic",'Ligand'="plain"))) %>%
  dplyr::arrange(cor_sum) -> gdsc_gene_rank
gdsc_gene_rank$color %>% as.character() ->gdsc_gene_rank$color
gdsc_gene_rank$size %>% as.character() ->gdsc_gene_rank$size

gdsc_plot_ready %>% 
  dplyr::distinct(drug_name, target_pathway, count) %>% 
  dplyr::group_by(target_pathway) %>% 
  dplyr::mutate(per = n() / count) %>% print(width = Inf) %>% 
  dplyr::arrange(per) %>% 
  dplyr::ungroup() %>% dplyr::select(drug_name, target_pathway, per) -> gdsc_drug_per

gdsc_plot_ready %>% 
  dplyr::left_join(gdsc_drug_per, by = c("drug_name", "target_pathway")) %>% 
  dplyr::group_by(drug_name) %>% 
  dplyr::mutate(drug_count = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(per, target_pathway, drug_count) %>% 
  dplyr::select(drug_name, target_pathway, drug_count, per) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(target_pathway = stringr::str_to_title(target_pathway)) -> gdsc_drug_rank_pre

pathway_color <- gdsc_drug_rank_pre %>% 
  dplyr::distinct(target_pathway, per) %>% 
  dplyr::arrange(per, target_pathway) %>% 
  dplyr::mutate(color = ggthemes::gdocs_pal()(17))

gdsc_drug_rank_pre %>% 
  dplyr::select(-per) %>% 
  dplyr::left_join(pathway_color, by = "target_pathway") -> drug_rank

pathway_color %>%
  ggplot(aes(x = 1, y = per, fill = color)) +
  geom_bar(stat = "identity", position = "stack", width = 0.001) +
  scale_fill_manual(values = pathway_color$color)
  
gdsc_plot_ready %>% 
  ggplot(aes(x = symbol, y = drug_name, color = cor_sprm)) +
  geom_point(aes(size = fdr)) +
  scale_x_discrete(limits = gdsc_gene_rank$symbol, expand = c(0.012,0.012)) +
  scale_y_discrete(limits = drug_rank$drug_name, expand = c(0.012,0.012), position = "right") +
  scale_color_gradient2(
    name = "Spearman Correlation",
    high = "red",
    mid = "white",
    low = "blue"
  ) +
  scale_size_continuous(
    name = "FDR"
  ) +
  # ggthemes::theme_gdocs() +
  theme(
    panel.background = element_rect(color = "black", fill = "white", size = 0.1),
    panel.grid=element_line(colour="grey",linetype="dashed"),
    panel.grid.major=element_line(colour="grey",linetype="dashed",size=0.2),
    
    axis.title = element_blank(),
    axis.text.x = element_text(size = 9, angle = 90, hjust = 1, vjust = 0.5,
                               color = gdsc_gene_rank$color,
                               face = gdsc_gene_rank$size),
    axis.text.y = element_text(size = 10, color = drug_rank$color,angle = 45),
    
    axis.ticks = element_line(color = "black"),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    # legend.key.width = unit(1,"cm"),
    # legend.key.heigh = unit(0.3,"cm"),
    legend.key = element_rect(fill="white",colour = "black")
  ) + guides(
    color = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = 0.5,
      barwidth = 10
    )
  )->p;p

ggsave(filename = 'gdsc_sig_drug.pdf', plot = p, device = "pdf", path = ccle_path, width = 8, height = 20)

# CTRP analysis
gene_list %>% 
  dplyr::inner_join(drug_ctrp, by = "symbol") -> ctrp_gene_list

fn_filter_drug_ctrp <- function(.x, .cor = 0.25, .p_val = 0.05) {
  .x %>% dplyr::filter(abs(cor_sprm) > .cor, p_val < .p_val)
}
ctrp_gene_list %>% 
  dplyr::mutate(cor_drug = purrr::map(.x = drug, .f = fn_filter_drug_ctrp)) %>% 
  tidyr::unnest(cor_drug) -> ctrp_gene_list_sig_drug

ctrp_gene_list_sig_drug %>% 
  dplyr::mutate(
    p_val = ifelse(-log10(p_val) > 50, 50, -log10(p_val)),
    cor_sprm = ifelse(cor_sprm > 0.5, 0.5, cor_sprm),
    cor_sprm = ifelse(cor_sprm < -0.5, -0.5, cor_sprm)
  ) -> ctrp_plot_ready

ctrp_plot_ready %>% 
  dplyr::group_by(symbol) %>% 
  dplyr::summarise(cor_sum = sum(cor_sprm)) %>% 
  dplyr::left_join(gene_list, by = "symbol") %>%
  dplyr::mutate(color = plyr::revalue(functionWithImmune, replace = c('TwoSide' = "blue", "Inhibit" = "red", "Activate" = "green"))) %>% 
  dplyr::mutate(size = plyr::revalue(type,replace = c('Receptor'="bold.italic",'Ligand'="plain"))) %>%
  dplyr::arrange(cor_sum) -> ctrp_gene_rank
ctrp_gene_rank$color %>% as.character() ->ctrp_gene_rank$color
ctrp_gene_rank$size %>% as.character() ->ctrp_gene_rank$size


ctrp_plot_ready %>% 
  ggplot(aes(x = symbol, y = drug_name, color = cor_sprm)) +
  geom_point(aes(size = p_val)) +
  scale_x_discrete(limits = ctrp_gene_rank$symbol, expand = c(0.012,0.012)) +
  scale_y_discrete(
    # limits = drug_rank$drug_name, 
    expand = c(0.012,0.012), 
    position = "right") +
  scale_color_gradient2(
    name = "Spearman Correlation",
    high = "red",
    mid = "white",
    low = "blue"
  ) +
  scale_size_continuous(
    name = "P-value",
    #range = c(1,5)
  ) +
  # ggthemes::theme_gdocs() +
  theme(
    panel.background = element_rect(color = "black", fill = "white", size = 0.1),
    panel.grid=element_line(colour="grey",linetype="dashed"),
    panel.grid.major=element_line(colour="grey",linetype="dashed",size=0.2),
    
    axis.title = element_blank(),
    axis.text.x = element_text(
      size = 9, 
      angle = 90,
      hjust = 1, 
      vjust = 0.5, 
      color = ctrp_gene_rank$color,
      face = gdsc_gene_rank$size),
    axis.text.y = element_text(
      color = drug_rank$color, 
      size = 10),
    
    axis.ticks = element_line(color = "black"),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    # legend.key.width = unit(1,"cm"),
    # legend.key.heigh = unit(0.3,"cm"),
    legend.key = element_rect(fill="white",colour = "black")
  ) + guides(
    color = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = 0.5,
      barwidth = 10
    )
  )->p;p

ggsave(filename = 'ctrp_sig_drug.pdf', plot = p, device = "pdf", path = ccle_path, width = 8, height = 20)

save.image(file = file.path(ccle_path, ".rda_02_ccle_a_gene_list.rda"))
load(file = file.path(ccle_path, ".rda_02_ccle_a_gene_list.rda"))
rm(list=ls())
