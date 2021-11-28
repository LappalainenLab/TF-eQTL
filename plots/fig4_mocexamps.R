#!/usr/bin/Rscript
##
##  fig4_mocexamps.R
##
##  EDF 5/12/2021
##

library(ggplot2)
library(dplyr)
library(tidyr)

setwd("~/projects/MANUSCRIPT_revisions/")

fi_tst_moc = read.table("data_tables/IRF1_kd/Fi_tst_combined_byvar_timecomb.mocgenes2.txt",
                        header=TRUE,sep='\t') %>%
  mutate(description=factor(description,
                            levels=rev(levels(description)))) %>%
  group_by(description) %>%
  mutate(count=n(),
         offsetg=ifelse(count==1,0,c(-.2,+.2)),
         gene_i = as.numeric(description))

fi_tst_moc %>%
  pivot_longer(cols=c(altp_C,altp_P,tot_C,tot_P)) %>%
  separate(name,sep='_',into=c('stat','cond'),) %>%
  pivot_wider(names_from=stat) %>%
  mutate(cond=factor(cond,levels=c('P','C'))) %>%
  group_by(gene,description,gene_i,exon,variant,chr,pos,ref,alt) %>%
  mutate(p_lab = ifelse(fi_p < 0.01, "**",
                        ifelse(fi_p < 0.05, "*",
                               ifelse(fi_p < 0.1, "+",""))),
         p_height = .95) %>%
  group_by(gene,description,gene_i,exon,variant,chr,pos,ref,alt) %>%
  mutate(gene_p_lab = paste(description,"p =",
                            ifelse(fi_p<0.01,round(fi_p,3),round(fi_p,2)))) %>%
  ggplot(aes(1,altp)) +
  geom_point(aes(alpha=cond,group=cond,
                 size=tot),
           position=position_dodge(.75),
           color='forestgreen') +
  geom_text(aes(1,p_height,
                label=p_lab)) +
  scale_alpha_discrete(range=c(.5,1)) +
  facet_wrap(~gene_p_lab,
             scales = 'free_y',
             nrow = 2) +
  theme_classic() +
  ylim(0,1) +
  ylab('Read Count') +
  xlab('Time Point')




fi_tst_moc %>%
  pivot_longer(cols=c(altp_C,altp_P,tot_C,tot_P)) %>%
  separate(name,sep='_',into=c('stat','cond'),) %>%
  pivot_wider(names_from=stat) %>%
  mutate(cond=factor(cond,levels=c('P','C'))) %>%
  mutate(p_lab = ifelse(fi_p < 0.01, "**",
                        ifelse(fi_p < 0.05, "*",
                               ifelse(fi_p < 0.1, "+",""))),
         p_height = .95) %>%
  ggplot() +
  geom_segment(data=fi_tst_moc,aes(x=gene_i+offsetg,xend=gene_i+offsetg,
                                   y=altp_P,yend=altp_C)) +
  geom_point(aes(gene_i+offsetg,altp,
                 alpha=cond,group=cond,
                 size=tot),
             color='forestgreen') +
  geom_text(aes(gene_i+offsetg,p_height,
                label=p_lab)) +
  scale_alpha_discrete(range=c(.5,1)) +
  theme_classic() +
  ylim(0,1) +
  ylab('ALT AF') +
  xlab('Gene') +
  coord_flip() +
  scale_x_continuous(#limits=c(8,1),
                     breaks=1:8,
                     labels=as.character(levels(fi_tst_moc$description))) 
ggsave("plots/fig4_mocexamps.pdf",
       width=4,height=2)


