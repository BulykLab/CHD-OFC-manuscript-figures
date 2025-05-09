library(ggplot2)
library(forcats)
library(cowplot)
library(patchwork)
library(dplyr)
library(ggsignif)
library(purrr)
library(stringr)

### Figure 2

## Figure 2 - TADA result

#### CHD
tfs <- read.table("data/Lambert_TF_list.gene.txt", header=F, stringsAsFactors=F)$V1


tada_chd <- read.table("data/chd.denovo_cc.with_loeuf.112823.txt", header=T)


tada_chd_filtered <- tada_chd %>% filter((qval < 0.1) & (oe_lof_upper < 1))

tada_chd_filtered$BF_pos_all <- apply(tada_chd_filtered[c('BF.misA', 'BF.misB','BF.pLoF')], 1, function(x) prod(x[x > 1], na.rm = TRUE))

###
BF.misA <- tada_chd_filtered %>% select(gene, BF.misA, BF_pos_all)  
BF.misA$type <- "MisA"
colnames(BF.misA) <- c("gene", "BF", "BF_pos_all", "type")
BF.misB <- tada_chd_filtered %>% select(gene, BF.misB, BF_pos_all) 
BF.misB$type <- "MisB"
colnames(BF.misB) <- c("gene", "BF", "BF_pos_all","type")
BF.plof <- tada_chd_filtered %>% select(gene, BF.pLoF, BF_pos_all) 
BF.plof$type <- "pLoF"
colnames(BF.plof) <- c("gene", "BF", "BF_pos_all","type")

chd_BF_merged <- rbind(BF.misA, BF.misB, BF.plof)
chd_BF_merged$TF <-  ifelse(chd_BF_merged$gene %in% tfs, "TF", "")
###

### Counts
count.misA <- tada_chd_filtered %>% select(gene, dn.misA, BF_pos_all) 
count.misA$type <- "denovo_misA"
colnames(count.misA) <- c("gene", "count", "BF_pos_all", "type")
count.misB <- tada_chd_filtered %>% select(gene, dn.misB, BF_pos_all)
count.misB$type <- "denovo_misB"
colnames(count.misB) <- c("gene", "count", "BF_pos_all","type")
count.plof <- tada_chd_filtered %>% select(gene, dn.pLoF, BF_pos_all)
count.plof$type <- "denovo_pLoF"
colnames(count.plof) <- c("gene", "count", "BF_pos_all","type")
count.caseplof <- tada_chd_filtered %>% select(gene, case.pLoF, BF_pos_all)
count.caseplof$type <- "case_pLoF"
colnames(count.caseplof) <- c("gene", "count", "BF_pos_all","type")
count.ctrlplof <- tada_chd_filtered %>% select(gene, ctrl.pLoF, BF_pos_all)
count.ctrlplof$type <- "control_pLoF"
colnames(count.ctrlplof) <- c("gene", "count", "BF_pos_all","type")

chd_count_merged <- rbind(count.misA, count.misB, count.plof, count.caseplof, count.ctrlplof)

chd_count_merged_filtered <- chd_count_merged %>% filter((count > 0))


title <- "Congenital heart defect"
p_2A_top <- chd_BF_merged %>%  mutate(type = factor(type, levels=c("MisB", "MisA", "pLoF"))) %>%
  ggplot(aes(x=fct_reorder(gene, log10(BF_pos_all), .desc=T))) + geom_col(aes(y=log10(BF), fill=type), color='black',size=0.4, width=1) +
  geom_point(aes(y=log10(BF_pos_all)+0.65, color=TF), size = 2) +
  theme_classic() +
  labs(title="Congenital heart defect", size=12) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=12), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),
        legend.title = element_blank(), legend.text = element_text(size=12),
        aspect.ratio = 0.5,
        plot.title = element_text(
          hjust = 1e-2,
          margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1))),
        plot.background = element_rect(fill='transparent', color=NA)
  ) +
  labs(y=expression(log[10]*"(BF)")) +
  scale_x_discrete(expand = expansion(mult = c(0,0))) +
  scale_y_continuous(expand = expansion(mult = c(0,0)), limits=c(0,26)) +
  scale_fill_manual(values=c("#ffe7b7", "#F7AA58", "#B84E43")) +
  scale_color_manual(values=c("white", "#00B0F0"))



chd_BF_merged %>%  mutate(type = factor(type, levels=c("MisB", "MisA", "pLoF"))) %>%
ggplot(aes(x=fct_reorder(gene, log10(BF_pos_all), .desc=T))) + geom_col(aes(y=log10(BF), fill=type), color='black',size=0.4, width=1) +
  geom_point(aes(y=log10(BF_pos_all)+0.65, color=TF), size = 2) +
  theme_classic() +
  labs(title="Congenital heart defect", size=12) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=12), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),
        legend.title = element_blank(), legend.text = element_text(size=12),
        aspect.ratio = 0.5,
        plot.title = element_text(
          hjust = 1e-2,
          margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1))),
        plot.background = element_rect(fill='transparent', color=NA)
  ) +
  labs(y=expression(log[10]*"(BF)")) +
  scale_x_discrete(expand = expansion(mult = c(0,0))) +
  scale_y_continuous(expand = expansion(mult = c(0,0)), limits=c(0,26)) +
  scale_fill_manual(values=c("#ffe7b7", "#F7AA58", "#E76254")) +
  scale_color_manual(values=c("white", "#00B0F0"))

p_2A_bottom <- chd_count_merged %>% 
  mutate(type = factor(type, levels=c("denovo_pLoF","denovo_misA", "denovo_misB", "case_pLoF", "control_pLoF"))) %>%
  ggplot(aes(x=fct_reorder(gene, log10(BF_pos_all), .desc=T), y=type)) + 
  geom_tile(aes(fill=type, alpha=(count>0)), color='black') +
  geom_text(aes(label=count, alpha=(count>0)), size=3) +
  geom_rect(xmin=0.5,xmax=55.5, ymin=2.5,ymax=5.5, colour="black", size=0.8, alpha=0) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, linewidth=1),
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=0.5, face ='italic'),
        axis.text.y = element_text(size=10),
        legend.title = element_blank(), legend.text = element_text(size=12),
        aspect.ratio = 5/46,
        plot.title = element_text(
          hjust = 1e-2,
          margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1))),
        plot.background = element_rect(fill='transparent', color=NA)
  ) +
  scale_fill_manual(values=c("#E76254","#F7AA58", "#ffe7b7",'#F3B0A9', "#FCECEA")) +
  scale_x_discrete(expand = expansion(mult = c(0,0))) +
  scale_y_discrete(expand = expansion(mult = c(0,0)), limits=rev) +
  scale_alpha_manual(values=c(0,1), guide = "none")
  



p_2A <- p_2A_top + p_2A_bottom + plot_layout(nrow=2)

p_2A

#### OC
tada_oc <- read.table("data/oc.denovo_cc.with_loeuf.112823.txt", header=T)


tada_oc_cc_filtered <- tada_oc %>% filter((qval < 0.1) & (oe_lof_upper < 1))

tada_oc_cc_filtered$BF_pos_all <- apply(tada_oc_cc_filtered[c('BF.misA', 'BF.misB','BF.pLoF')], 1, function(x) prod(x[x > 1], na.rm = TRUE))



###
BF.misA <- tada_oc_cc_filtered %>% select(gene, BF.misA, BF_pos_all)  
BF.misA$type <- "MisA"
colnames(BF.misA) <- c("gene", "BF", "BF_pos_all", "type")
BF.misB <- tada_oc_cc_filtered %>% select(gene, BF.misB, BF_pos_all)
BF.misB$type <- "MisB"
colnames(BF.misB) <- c("gene", "BF", "BF_pos_all","type")
BF.plof <- tada_oc_cc_filtered %>% select(gene, BF.pLoF, BF_pos_all)
BF.plof$type <- "pLoF"
colnames(BF.plof) <- c("gene", "BF", "BF_pos_all","type")

tada_oc_cc_BF_merged <- rbind(BF.misA, BF.misB, BF.plof)
tada_oc_cc_BF_merged$TF <-  ifelse(tada_oc_cc_BF_merged$gene %in% tfs, "TF", "")
###

### Counts
count.misA <- tada_oc_cc_filtered %>% select(gene, dn.misA, BF_pos_all) 
count.misA$type <- "denovo_misA"
colnames(count.misA) <- c("gene", "count", "BF_pos_all", "type")
count.misB <- tada_oc_cc_filtered %>% select(gene, dn.misB, BF_pos_all) 
count.misB$type <- "denovo_misB"
colnames(count.misB) <- c("gene", "count", "BF_pos_all","type")
count.plof <- tada_oc_cc_filtered %>% select(gene, dn.pLoF, BF_pos_all)
count.plof$type <- "denovo_pLoF"
colnames(count.plof) <- c("gene", "count", "BF_pos_all","type")
count.caseplof <- tada_oc_cc_filtered %>% select(gene, case.pLoF, BF_pos_all)
count.caseplof$type <- "case_pLoF"
colnames(count.caseplof) <- c("gene", "count", "BF_pos_all","type")
count.ctrlplof <- tada_oc_cc_filtered %>% select(gene, ctrl.pLoF, BF_pos_all)
count.ctrlplof$type <- "control_pLoF"
colnames(count.ctrlplof) <- c("gene", "count", "BF_pos_all","type")

oc_count_merged <- rbind(count.misA, count.misB, count.plof, count.caseplof, count.ctrlplof)
oc_count_merged_filtered <- oc_count_merged %>% filter((count > 0))

title <- "Orofacial cleft"
p_2B_top <- tada_oc_cc_BF_merged %>% mutate(type = factor(type, levels=c("MisB", "MisA", "pLoF"))) %>%
ggplot(aes(x=fct_reorder(gene, log10(BF_pos_all), .desc=T))) + geom_col(aes( y=log10(BF), fill=type), color='black',size=0.4, width=1) +
  theme_classic() +
  geom_point(aes(y=log10(BF_pos_all)+1.67, color=TF), size = 2) +
  labs(title="Orofacial cleft", size=12) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=12), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),
        legend.position = 'none', aspect.ratio = 0.3,
        plot.title = element_text(
          hjust = 1e-2,
          margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1))),
        plot.background = element_rect(fill='transparent', color=NA)
  ) +
  labs(y=expression(log[10]*"(BF)")) +
  scale_x_discrete(expand = expansion(mult = c(0,0))) +
  scale_y_continuous(expand = expansion(mult = c(0,0)), limits=c(0,24)) +
  scale_fill_manual(values=c("#ffe7b7", "#F7AA58", "#B84E43"))+
  scale_color_manual(values=c("white", "#00B0F0"))



p_2B_bottom <- oc_count_merged %>% 
  mutate(type = factor(type, levels=c("denovo_pLoF", "denovo_misA", "denovo_misB", "case_pLoF", "control_pLoF"))) %>%
  ggplot(aes(x=fct_reorder(gene, log10(BF_pos_all), .desc=T), y=type)) + 
  geom_tile(aes(fill=type, alpha=(count>0)), color='black') +
  geom_text(aes(label=count, alpha=(count>0)), size=3) +
  geom_rect(xmin=0.5,xmax=55.5, ymin=2.5,ymax=5.5, colour="black", size=0.8, alpha=0) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, size=1),
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=0.5, face ='italic'),
        axis.text.y = element_text(size=10),
        legend.position = 'none',
        aspect.ratio = 5/22,
        plot.title = element_text(
          hjust = 1e-2,
          margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1))),
        plot.background = element_rect(fill='transparent', color=NA)
  ) +
  scale_fill_manual(values=c("#E76254","#F7AA58", "#ffe7b7",'#F3B0A9', "#FCECEA")) +
  scale_x_discrete(expand = expansion(mult = c(0,0))) +
  scale_y_discrete(expand = expansion(mult = c(0,0)), limits=rev) +
  scale_alpha_manual(values=c(0,1), guide = "none") 

p_2B <- p_2B_top + p_2B_bottom + plot_layout(nrow=2)
p_2B




### TF enrichment

## CHD
fisher.test(matrix(c(14,1625,58,17488-1639 - 58), nrow = 2), alternative = 'greater')

## OC
fisher.test(matrix(c(8,1631,28,17488-7 - 1632 - 28), nrow = 2), alternative = 'greater')


### Separately
ggsave('figures/Figure2A.pdf', 
       plot = p_2A, bg='transparent',
       width=250, height=160, units="mm")


ggsave('figures/Figure2B.pdf', 
       plot = p_2B, bg='transparent',
       width=120, height=83, units="mm")
