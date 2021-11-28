#!/usr/bin/Rscript
##
##  fig3_asb.R
##
##  EDF 5/7/21
##


library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)

setwd("~/projects/MANUSCRIPT_revisions//")

asb_annot_filt1 = read.table("data_tables/overlap/asb/ASB_susan_filtexpASB1_rev2.txt",
                             header=TRUE,sep='\t')

asb_annot_comb = read.table("data_tables/overlap/asb/ASB_susan_allcombined_rev2.txt",
                            header=TRUE,sep='\t')

asb_annot_filt1 %>%
 ggplot(aes(fct_reorder(tf, real_stat),real_stat,
           size=exp_ASB_corr)) +
  geom_col(width=.05) +
  geom_hline(yintercept=0) +
  geom_point(aes(col = -log10(p_perm_tf))) +
  scale_color_gradient(low='gray',high='darkorchid4',
                       limits=c(0,2)) +
  coord_flip() +
  theme_classic()
ggsave("plots/fig3_ASB_tffilt1.pdf",
       width=4, height=3)

asb_annot_filt2 = read.table("data_tables/overlap/asb/ASB_susan_filtexpASB2_rev2.txt",
                             header=TRUE,sep='\t')

asb_annot_filt2 %>%
  ggplot(aes(fct_reorder(tf, real_stat),real_stat,
             size=exp_ASB_corr)) +
  geom_col(width=.05) +
  geom_point(aes(col = -log10(p_perm_tf))) +
  scale_color_gradient(low='gray',high='darkorchid4',
                       limits=c(0,2)) +
  coord_flip() +
  theme_classic()

asb_annot_filt1 %>%
  rbind(asb_annot_comb) %>%
  mutate(tf_order = factor(tf, 
                           levels = c(levels(fct_reorder(tf, real_stat))[-16],'all'))) %>%
  ggplot(aes(tf_order,real_stat,
             size=exp_ASB_corr)) +
  geom_col(width=.05) +
  geom_hline(yintercept=0) +
  geom_point(aes(col = -log10(p_perm_tf))) +
  scale_color_gradient2(low='gray',mid='darkorchid4',high='darkorchid4',
                        limits=c(0,4),
                        midpoint=2) +
  coord_flip() +
  theme_classic() +
  xlab('Transcription Factor') +
  ylab('ASB Enrichment')
ggsave("plots/fig3_ASB_tffilt1_pluscomb.pdf",
       width=4, height=3)

asb_annot_filt2 %>%
  rbind(asb_annot_comb) %>%
  mutate(tf_order = factor(tf, 
                           levels = c(levels(fct_reorder(tf, real_stat))[-8],'all'))) %>%
  ggplot(aes(tf_order,real_stat,
             size=exp_ASB_corr)) +
  geom_col(width=.05) +
  geom_hline(yintercept=0) +
  geom_point(aes(col = -log10(p_perm_tf))) +
  scale_color_gradient2(low='gray',mid='darkorchid4',high='darkorchid4',
                        limits=c(0,4),
                        midpoint=2) +
  coord_flip() +
  theme_classic() +
  xlab('Transcription Factor') +
  ylab('ASB Enrichment')
ggsave("plots/fig3_ASB_tffilt2_pluscomb.pdf",
       width=4, height=3)
