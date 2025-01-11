library(ggplot2)
library(forcats)
library(cowplot)
library(patchwork)
library(dplyr)
library(qvalue)
library(ggsignif)
library(purrr)

### Figure 3

## Figure 3 - TADA result

tada_chd_comparison <- read.table("~/Library/CloudStorage/GoogleDrive-rjeong@g.harvard.edu/My Drive/Manuscript/KidsFirst/draft/data/tada_result/chd.tada_comparison.080623.txt", header=T)


p_4_3 <-
tada_chd_comparison %>% filter(category != "both") %>% 
  ggplot(aes(x=category, y=cds_length)) + theme_classic() +
  geom_violin(aes(color=category, fill=category), alpha = 0.5) +
  geom_boxplot(color='black', width=0.2) +
  scale_y_continuous(limits = c(0,6250)) +
  labs(y="CDS length (aa)") +
  scale_x_discrete(labels=c(expression(italic("de novo")*" & case/control"), expression(italic("de novo")*" only"))) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=14), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.position = 'none', aspect.ratio = 1.1
  ) +
  geom_signif(
    y_position = 6200, xmin = c(1), xmax = c(2),
    annotation = rep("*", 1), tip_length = 0.01, textsize = 6, vjust = 0.4
  ) 


ggsave('~/Library/CloudStorage/GoogleDrive-rjeong@g.harvard.edu/My Drive/Manuscript/KidsFirst/draft/figures/Figure3.pdf', 
       plot = p_4_3, 
       width=120, height=120, units="mm")

## p = 0.0197
wilcox.test( x = as.numeric((tada_chd_comparison %>% filter(category == "denovo_cc_only"))$cds_length),
             y = as.numeric((tada_chd_comparison %>% filter(category == "denovo_only"))$cds_length), alternative = 'greater')


### OC

tada_oc_comparison <- read.table("~/Library/CloudStorage/GoogleDrive-rjeong@g.harvard.edu/My Drive/Manuscript/KidsFirst/draft/data/tada_result/oc.tada_comparison.120123.txt", header=T)


p_3_supp <-
  tada_oc_comparison %>% filter(category != "both") %>% 
  ggplot(aes(x=category, y=cds_length)) + theme_classic() +
  geom_violin(aes(color=category, fill=category), alpha = 0.5) +
  geom_boxplot(color='black', width=0.2) +
  scale_y_continuous(limits = c(0,16800)) +
  labs(y="CDS length (aa)") +
  scale_x_discrete(labels=c(expression(italic("de novo")*" & case/control"), expression(italic("de novo")*" only"))) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=14), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.position = 'none', aspect.ratio = 1.1
  ) +
  geom_signif(
    y_position = 16800, xmin = c(1), xmax = c(2),
    annotation = rep("*", 1), tip_length = 0.01, textsize = 6, vjust = 0.4
  ) 
p_3_supp

ggsave('~/Library/CloudStorage/GoogleDrive-rjeong@g.harvard.edu/My Drive/Manuscript/KidsFirst/draft/figures/Figure_supp3.pdf', 
       plot = p_3_supp, 
       width=120, height=120, units="mm")


## p = 0.01449
wilcox.test( x = as.numeric((tada_oc_comparison %>% filter(category == "denovo_cc_only"))$cds_length),
             y = as.numeric((tada_oc_comparison %>% filter(category == "denovo_only"))$cds_length), alternative = 'greater')



