library(ggplot2)
`%>%` <- magrittr::`%>%`


# Path
tcga_path <- "S:/study/GSCALite/data"
expr_path <- "S:/study/GSCALite"
expr_path_a <- file.path(expr_path, "all_expr")
mirna_path <- file.path(expr_path, "mirna")

# load methylation and gene list
mirna_target <- readr::read_rds(file.path(tcga_path, "mirna_gene_target.rds.gz"))
gene_list_path <- "S:/study/生存分析/免疫检查点project/免疫检查点"
gene_list <- read.table(file.path(gene_list_path, "all.entrez_id-gene_id"),header=T)
gene_list$symbol %>% as.character() ->gene_list$symbol
gene_type<-read.table(file.path(gene_list_path,"checkpoint.type"),header=T)
gene_list<-dplyr::left_join(gene_list,gene_type,by="symbol")
mirna_expr <- readr::read_rds(file.path(tcga_path, "pancan33_mirna_expr.rds.gz"))

gene_list %>% 
#  dplyr::filter(pathway == "autophagesome formation-core") %>% 
  dplyr::select(symbol) %>% 
  dplyr::left_join(mirna_target, by = "symbol") -> gene_list_mirna
gene_list_mirna %>% readr::write_tsv(path = file.path(mirna_path, "02_a_gene_list_mirna.tsv"))

fn_filter_gene_list <- function(.x, gene_list) {
  gene_list_mirna %>%
    dplyr::rename(name = mirna) %>% 
    tidyr::drop_na() %>% 
    dplyr::inner_join(.x, by = "name")
}

mirna_expr %>% 
  dplyr::mutate(mirna = purrr::map(.x = mirna, .f = fn_filter_gene_list)) -> gene_list_mirna_expr 

gene_list_mirna_expr %>% 
  tidyr::unnest()



save.image(file = file.path(mirna_path, ".rda_02_mirna_a_gene_list_target.rda"))
rm(list=ls())
load(file = file.path(mirna_path, ".rda_02_mirna_a_gene_list_target.rda"))
