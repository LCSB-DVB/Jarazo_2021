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

N1=read.table("S:/HCS_Platform/Data/JavierJarazo/WB/Normalized_values.csv",header=TRUE, sep=',',stringsAsFactors = FALSE)
N1 %>%
  dplyr::mutate(Handle = paste0(AreaName,"_",Replicate)) %>% 
  pivot_wider(names_from = "Marker",values_from="Value")->Table_all
Corrected =c("825_GC","2122_GC")
Table_all$Category = ifelse(Table_all$AreaName %in% Corrected,"Corrected","Patient")  

assign_significance = function(pvals){
  significance = rep(NA,length(pvals))
  significance[pvals >= 0.05] = 'ns'
  significance[pvals <0.05 & pvals >0.01] = '*' 
  significance[pvals <=0.01 & pvals >0.001] = '**'
  significance[pvals <=0.001 & pvals >0.0001] = '***'
  significance[pvals <=0.0001] = '****'
  return(significance)
}

p1 = wilcox.test(TH ~ Category, data=Table_all)$p.value
p2 = wilcox.test(Tuj1 ~ Category, data=Table_all)$p.value
p3 = wilcox.test(GFAP ~ Category, data=Table_all)$p.value
c1 = "Corrected - Patient"
mytable = data.frame(Comparison = c(c1,c1,c1), Feature=c("TH","Tuj1","GFAP"), pvalue = c(p1,p2,p3))
mytable$significance = assign_significance(mytable$pvalue)
write.table(mytable, "KW_WB_GA.csv", sep = ',', row.names = FALSE)

Table_all %>% 
  group_by(Category) %>% 
  summarise(across(where(is.numeric),list(mean =~mean(.x, na.rm = TRUE),sd =~sd(.x, na.rm = TRUE)))) -> toto

graph =list()
graph[[1]]=ggplot(toto, aes(x=Category,y=TH_mean, fill=Category))+
  geom_bar(stat="identity",  position=position_dodge(0.8),alpha=0.5) +
  geom_errorbar(aes(ymin = TH_mean-TH_sd, ymax= TH_mean+TH_sd),width = 0.2, position=position_dodge(width=0.8),size=0.5) +
  scale_fill_manual(values=c("lightgoldenrod1","firebrick1")) +
  scale_y_continuous (limits = c(0,1.5), expand = c(0,0),breaks = pretty_breaks(n = 4)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none") +
  theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),axis.title.x=element_blank(),axis.title.y=element_blank()) +
  theme(plot.title = element_text(size =10, hjust = 0.5))+
  ggtitle("TH")

graph[[2]]=ggplot(toto, aes(x=Category,y=Tuj1_mean, fill=Category))+
  geom_bar(stat="identity",  position=position_dodge(0.8),alpha=0.5) +
  geom_errorbar(aes(ymin = Tuj1_mean-Tuj1_sd, ymax= Tuj1_mean+Tuj1_sd),width = 0.2, position=position_dodge(width=0.8),size=0.5) +
  scale_fill_manual(values=c("lightgoldenrod1","firebrick1")) +
  scale_y_continuous (limits = c(0,2), expand = c(0,0),breaks = pretty_breaks(n = 4)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none") +
  theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),axis.title.x=element_blank(),axis.title.y=element_blank()) +
  theme(plot.title = element_text(size =10, hjust = 0.5))+
  ggtitle("TUBB3")

graph[[3]]=ggplot(toto, aes(x=Category,y=GFAP_mean, fill=Category))+
  geom_bar(stat="identity",  position=position_dodge(0.8),alpha=0.5) +
  geom_errorbar(aes(ymin = GFAP_mean-GFAP_sd, ymax= GFAP_mean+GFAP_sd),width = 0.2, position=position_dodge(width=0.8),size=0.5) +
  scale_fill_manual(values=c("lightgoldenrod1","firebrick1")) +
  scale_y_continuous (limits = c(0,5), expand = c(0,0),breaks = pretty_breaks(n = 4)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none") +
  theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),axis.title.x=element_blank(),axis.title.y=element_blank()) +
  theme(plot.title = element_text(size =10, hjust = 0.5))+
  ggtitle("GFAP")

p=ggdraw() +
  draw_plot(graph[[1]], x = 0, y = .2, width = .15, height = .2) +
  draw_plot(graph[[2]], x = .15, y = .2, width = .15, height = .2) +
  draw_plot(graph[[3]], x = .3, y = .2, width = .15, height = .2)

ggsave("Isogenics_WB.pdf", plot=p, dpi=320, width = 11.69 , height = 8.27, units = "in",device=cairo_pdf)


