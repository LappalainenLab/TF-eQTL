#!/usr/bin/Rscript
##
##  sfig_caviar.R
##
##
##

library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/projects/MANUSCRIPT_revisions/")

cav_hits = read.table("~/projects/TFi-eQTL/variant_sets/caviar_var.95set.eqtls.MAF05.overlap.genesort.txt",
                      header=TRUE,sep='\t')

table(cav_hits$gene) %>%
  hist(main="Variants in fine-mapped set per gene", 
       xlab="Variants",
       ylab="Genes")
median(table(cav_hits$gene))

table(cav_hits$var) %>%
  hist(main="Genes associated with each variant", 
       xlab="Genes",
       ylab="Variants")
median(table(cav_hits$var))

cross_hits = read.table("~/projects/MANUSCRIPT_revisions/data_tables/tf_eqtls/cross_tiss_expr_tfeqtls.fdr05.txt.gz",
                        header=TRUE,sep='\t')
prot_hits = read.table("~/projects/MANUSCRIPT_revisions/data_tables/tf_eqtls/sig_prot_corrs.fdr05.cols.filtered.txt",
                       header=TRUE,sep='\t')
with_hits = read.table("~/projects/MANUSCRIPT_revisions/data_tables/tf_eqtls/all.16tiss.fdr05.sig_eqtls.txt",
                       header=TRUE,sep='\t')
moc_hits = read.table("~/projects/MANUSCRIPT_revisions/data_tables/tf_eqtls/dual_evidence_tfeqtls.txt",
                      header=TRUE,sep='\t')

cav_hits %>%
  group_by(gene) %>%
  summarize(n_var=n()) %>%
  mutate(cross = gene %in% as.character(cross_hits$gene),
         prot = gene %in% as.character(prot_hits$phenotype_id),
         with = gene %in% as.character(with_hits$gene),
         moc = gene %in% as.character(moc_hits$phenotype_id)) %>%
  pivot_longer(cols=c(cross,prot,with,moc),
               values_to='corr') %>%
  mutate(name=factor(name,
                     levels = c('cross','prot','with','moc'))) %>%
  ggplot(aes(name,n_var)) +
  geom_boxplot(aes(col=corr)) +
  theme_classic() +
  xlab("Corr type") +
  ylab("Number of vars per gene")
ggsave("plots/sfig12.pdf")  




