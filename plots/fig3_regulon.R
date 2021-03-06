#!/usr/bin/Rscript
##
##  fig3_regulon.R
##
##  EDF 5/7/21
##


library(dplyr)
library(ggplot2)
library(tidyr)


setwd("~/projects/MANUSCRIPT_revisions/")

reg_fi_tsts_bycat = read.table("data_tables/overlap/regulon/regulon_enrich_bycat.txt",
                               header=TRUE,sep='\t')
reg_fi_tsts_overall = read.table("data_tables/overlap/regulon/regulon_enrich_overall.txt",
                                 header=TRUE,sep='\t')

set.seed(100)
reg_fi_tsts_bycat %>%
  mutate(reg_lab = factor(ifelse(reg=='cure','curated',as.character(reg)),
                          levels = c('any','curated','chip','motif','expr'))) %>%
  filter(reg_lab!='curated') %>%
  arrange(-fi_p) %>%
  ggplot(aes(reg_lab,fi_logOR_plot)) +
  geom_hline(yintercept=0) +
  geom_violin() +
  geom_point(aes(col=-log10(fi_p)),
             position=position_jitter(.2),
             size=.5) +
  geom_point(data=reg_fi_tsts_overall %>% filter(reg_lab!='curated'),
             aes(reg_lab,fi_log_OR),
             size=2) +
  scale_color_gradientn(colors=c('gray','darkorchid4','darkorchid4','darkorchid4'),
                        limits=c(0,10),
                        breaks=c(0,3,6,9)) +
  theme_classic() +
  scale_y_continuous(labels=c('-Inf',-2,0,2,4,6),
                     breaks=c(-3,-2,0,2,4,6)) +
  xlab("Regulon Source") +
  ylab("log2(Odds Ratio)")
ggsave("plots/fig3_regulon.pdf",
       height=2,width=4)

