library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(data.table)
library(dplyr)
library(stringr)
library(janitor)
library(scales)
# library(plyr)
library(tidyverse)
library(cowplot)
library(ggpubr)

N1=read.table("S:/HCS_Platform/Data/JavierJarazo/WB/summary all lines _ used for quantification_SSS.csv",header=TRUE, sep=',',stringsAsFactors = FALSE)

assign_significance = function(pvals){
  significance = rep(NA,length(pvals))
  significance[pvals >= 0.05] = 'ns'
  significance[pvals <0.05 & pvals >0.01] = '*' 
  significance[pvals <=0.01 & pvals >0.001] = '**'
  significance[pvals <=0.001 & pvals >0.0001] = '***'
  significance[pvals <=0.0001] = '****'
  return(significance)
}

p1 = wilcox.test(LC3.2.LC3.1_norm ~ condition, data=N1)$p.value
p2 = wilcox.test(P62.b.actin_norm ~ condition, data=N1)$p.value
c1 = "Control - Patient"
mytable = data.frame(Comparison = c(c1,c1), Feature=c("LC3.2.LC3.1_norm","P62.b.actin_norm"), pvalue = c(p1,p2))
mytable$significance = assign_significance(mytable$pvalue)
write.table(mytable, "MW_WB.csv", sep = ',', row.names = FALSE)

N1 %>% 
  group_by(condition) %>% 
  summarise(sd = sd(LC3.2.LC3.1_norm),mean= mean(LC3.2.LC3.1_norm)) %>% 
  ungroup() %>% 
  full_join(N1,by="condition")->toto

p=ggplot(toto, aes(x=condition,y=mean, fill=condition))+
  geom_col() +
  geom_errorbar(aes(ymin = mean-sd, ymax= mean+sd), width=0.2,position="identity") +
  scale_x_discrete(limits =c("WT","PINK1"))+
  scale_fill_manual(values=c("firebrick1","cornflowerblue")) +
  scale_y_continuous (limits = c(0,25), expand = c(0,0),breaks = pretty_breaks(n = 4)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none") +
  theme(axis.text.x = element_text(size = 5),axis.text.y = element_text(size = 5),axis.title.x=element_blank(),axis.title.y=element_blank()) +
  theme(plot.title = element_text(size =5, hjust = 0.5))

ggdraw() +
  draw_plot(p, x = 0, y = .45, width = .15, height = .15)