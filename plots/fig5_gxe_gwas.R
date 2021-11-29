#!/usr/bin/Rscript
##
##  fig5_gxe_gwas.R
##
##  EDF 7/9/21
##

setwd("~/projects/MANUSCRIPT_revisions//")

library(dplyr)
library(tidyr)
library(ggvenn)


rep_data = read.table("~/data/Findley2021/cASE_GxE_rep.tab",
                      header=TRUE,sep='\t')
sig_genes = read.table("data_tables/tf_eqtls/dual_evidence_tfeqtls.txt",
                       header=TRUE, sep='\t') %>%
  separate(phenotype_id, c('ensg'), sep='[.]',
           remove=FALSE) %>%
  #separate_rows(variants, datasets, sep=',') %>%
  unite(tiss_gene, datasets, phenotype_id, remove=FALSE)

all_eqtls = read.table('~/projects/examples/input_files/GTEx_Analysis_v8_eQTL.all_tissues.sig_egenes.txt.gz',
                       header=TRUE,sep='\t')
all_vars = read.table("~/projects/examples/input_files/caviar_var.95set.eqtls.MAF05.overlap.genesort.txt",
                      header=TRUE,sep='\t')
all_eqtls_tested = all_eqtls %>%
  filter(gene_id %in% as.character(all_vars$gene)) %>%
  unite(tiss_gene, tissue, gene_id, remove=FALSE) %>%
  separate(gene_id, c('ensg'), sep='[.]',
           remove=FALSE)


genes_overlap_list = list(TF.eQTL = filter(sig_genes, 
                                           phenotype_id %in% all_eqtls_tested$gene_id) %>% 
                            pull(ensg) %>% as.character(),
                          GxE = filter(rep_data, 
                                       ensg %in% all_eqtls_tested$ensg) %>% 
                            pull(ensg) %>% as.character())
ggvenn( genes_overlap_list,
        fill_color=c('darkorchid4', 'gray') )
ggsave('plots/fig5_GxE_gene.pdf',
       width=4, height=3)

fisher.test(unique(all_eqtls_tested$ensg) %in% sig_genes$ensg,
            unique(all_eqtls_tested$ensg) %in% rep_data$ensg)



moc_hits = read.table("data_tables/tf_eqtls/dual_evidence_tfeqtls.txt",
                      header=TRUE,sep='\t') %>%
  separate(phenotype_id, c('ensg'), sep='[.]',
           remove=FALSE) %>%
  separate_rows(variants, datasets, sep=',') %>%
  unite(tiss_gene, datasets, phenotype_id, remove=FALSE)
coloc_table = read.table("~/data/gtex/v8/GWAS_enloc/enloc_ENLOC_rcp_gt_0.5_with_gwas_pval.tsv",
                         header=TRUE,sep='\t') %>%
  unite(tiss_gene, tissue, gene_id, remove=FALSE)

genes_overlap_list = list(TF.eQTL = filter(moc_hits, 
                                           phenotype_id %in% all_eqtls_tested$gene_id) %>% 
                            pull(phenotype_id) %>% as.character(),
                          GWAS = filter(coloc_table, 
                                        gene_id%in% all_eqtls_tested$gene_id) %>% 
                            pull(gene_id) %>% as.character())
ggvenn( genes_overlap_list,
        fill_color=c('darkorchid4', 'gray') )
ggsave('plots/fig5_GWAScoloc_gene.pdf',
       width=4, height=3)

fisher.test(unique(all_eqtls_tested$gene_id) %in% moc_hits$phenotype_id,
            unique(all_eqtls_tested$gene_id) %in% coloc_table$gene_id)

fisher.test(unique(all_eqtls_tested$gene_id) %in% moc_hits$phenotype_id,
            unique(all_eqtls_tested$gene_id) %in% coloc_table$gene_id)$p.value


all_eqtls_tested_tiss = all_eqtls_tested %>%
  filter(tissue %in% as.character(moc_hits$datasets))

tissue_genes_overlap_list = list(TF.eQTL = filter(moc_hits, tiss_gene %in% all_eqtls_tested_tiss$tiss_gene) %>% 
                                   pull(tiss_gene) %>% as.character(),
                                 GWAS = filter(coloc_table, tiss_gene %in% all_eqtls_tested_tiss$tiss_gene) %>% 
                                   pull(tiss_gene) %>% as.character())
ggvenn( tissue_genes_overlap_list,
        fill_color=c('darkorchid4', 'gray') )
ggsave('plots/fig5_GWAScoloc_tissgene.pdf',
       width=4, height=3)

fisher.test(unique(all_eqtls_tested_tiss$tiss_gene) %in% moc_hits$tiss_gene,
            unique(all_eqtls_tested_tiss$tiss_gene) %in% coloc_table$tiss_gene)





