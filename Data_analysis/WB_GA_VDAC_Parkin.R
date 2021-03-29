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

N1=read.table("S:/HCS_Platform/Data/JavierJarazo/WB/Normalized_values_Parkin_VDAC.csv",header=TRUE, sep=',',stringsAsFactors = FALSE)
N1 %>%
  dplyr::mutate(Handle = paste0(AreaName,"_",Replicate)) %>% 
  pivot_wider(names_from = "Marker",values_from="Value")->Table_all
Corrected =c("825_GC","2122_GC")
Table_all$Category = ifelse(Table_all$AreaName %in% Corrected,"Corrected","Patient")  

Table_all %>% 
  group_by(Category,Treatment) %>% 
  summarise(across(where(is.numeric),list(mean =~mean(.x, na.rm = TRUE),sd =~sd(.x, na.rm = TRUE)))) -> toto
toto$Handle_cat_tr = paste0(toto$Category,"_",toto$Treatment)

category_levels=c("Corrected_DMSO","Corrected_CCCP","Patient_DMSO","Patient_CCCP")
mycolors = c("Corrected_DMSO" = "lightgoldenrod1","Corrected_CCCP" = "lightgoldenrod4", "Patient_DMSO" = "firebrick1","Patient_CCCP" ="indianred4" )
toto$Handle_cat_tr = factor(toto$Handle_cat_tr, levels= category_levels)


graph =list()
graph[[1]]=ggplot(toto, aes(x=Handle_cat_tr,y=Parkin_mean, fill=Handle_cat_tr))+
  geom_bar(aes(group=interaction(Category,Treatment)),stat="identity",  position=position_dodge(1),alpha=0.5) +
  geom_errorbar(aes(group=interaction(Category,Treatment),ymin = Parkin_mean-Parkin_sd, ymax= Parkin_mean+Parkin_sd),width = 0.2, position=position_dodge(width=1),size=0.5) +
  scale_fill_manual(values=mycolors) +
  scale_x_discrete(limits=category_levels) +
  scale_y_continuous (limits = c(0,4), expand = c(0,0),breaks = pretty_breaks(n = 4)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none") +
  theme(axis.text.x = element_text(size = 6),axis.text.y = element_text(size = 8),axis.title.x=element_blank(),axis.title.y=element_blank()) +
  theme(plot.title = element_text(size =10, hjust = 0.5))+
  ggtitle("Parkin")

graph[[2]]=ggplot(toto, aes(x=Handle_cat_tr,y=UbVDAC_mean, fill=Handle_cat_tr))+
  geom_bar(aes(group=interaction(Category,Treatment)),stat="identity",  position=position_dodge(1),alpha=0.5) +
  geom_errorbar(aes(group=interaction(Category,Treatment),ymin = UbVDAC_mean-UbVDAC_sd, ymax= UbVDAC_mean+UbVDAC_sd),width = 0.2, position=position_dodge(width=1),size=0.5) +
  scale_fill_manual(values=mycolors) +
  scale_x_discrete(limits=category_levels) +
  scale_y_continuous (limits = c(0,15), expand = c(0,0),breaks = pretty_breaks(n = 4)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none") +
  theme(axis.text.x = element_text(size = 6),axis.text.y = element_text(size = 8),axis.title.x=element_blank(),axis.title.y=element_blank()) +
  theme(plot.title = element_text(size =10, hjust = 0.5))+
  ggtitle("Ub-VDAC")

#Get the legend as a separate plot
legend_plot = ggplot(toto, aes(x=Handle_cat_tr,y=Parkin_mean, fill=Handle_cat_tr))+
  geom_bar(aes(group=interaction(Category,Treatment)),stat="identity",  position=position_dodge(1),alpha=0.5) +
  geom_errorbar(aes(group=interaction(Category,Treatment),ymin = Parkin_mean-Parkin_sd, ymax= Parkin_mean+Parkin_sd),width = 0.2, position=position_dodge(width=1),size=0.5) +
  scale_fill_manual(values=mycolors) +
  scale_x_discrete(limits=category_levels) +
  scale_y_continuous (limits = c(0,4), expand = c(0,0),breaks = pretty_breaks(n = 4)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank())

my_legend = get_legend(legend_plot)

p=ggdraw() +
  draw_plot(graph[[1]], x = 0, y = .2, width = .15, height = .2) +
  draw_plot(graph[[2]], x = .15, y = .2, width = .15, height = .2) +
  draw_plot(my_legend, x = .63, y = .2, width = .00000001, height = .00000001)

ggsave("Isogenics_WB_Parkin_UbVDAC.pdf", plot=p, dpi=320, width = 11.69 , height = 8.27, units = "in",device=cairo_pdf)


