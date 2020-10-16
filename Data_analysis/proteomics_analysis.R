library(tidyverse)
library(readxl)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
# library(seriation)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

# file "raw" has all proteins, other files (starting in "differential") have DE proteins for each comparison.
##############
# Literature #
##############
interesting_proteins = read.table("Proteine_interesanti.txt",header=1,sep=";")[,1]

######################
### Intensity data ###
######################
filename = "rawData_E-values_ProteinLevel_BO_v1.xlsx"

# read headers
headers <- read_excel(filename, col_names = FALSE, n_max = 8)

# merge all sample metadata in a label
column_labels <- headers %>% summarize_all(str_c, collapse = "--")

# read values
raw_data <- read_excel(filename, col_names = TRUE, skip = 9)
# assign metadata as column name
colnames(raw_data)[6:41] = column_labels[2:37]
# Now the file is a wide format file, with the multiple headers all in one header

# transform to long format and divide the metadata in columns
raw_data_long <- raw_data %>%
  pivot_longer(cols=`BO001--Human--Organoide--1--WT_UN_T_10_R_1--0.5-1.106 cells--WT--10`:`BO036--Human--Organoide--36--PAT_CD_T_30_R_3--0.5-1.106 cells--Patient--30`,
               names_to = "metadata", values_to = "values") %>%
  separate(metadata, headers %>% pull(1),sep="--",remove=TRUE) %>%
  mutate(timepoint = as.numeric(timepoint)) %>%
  mutate(treatment = if_else(grepl("CD",SampleID_Customer),"CD","UN")) %>%
  mutate(protein_id = str_c(`Sciomics ID`,if_else(is.na(Name),"",Name),sep="::"))


############
### DEPs ###
############
# add to the raw data the logFC of DE proteins in each comparison
extract_DEP <- function(filename){
  comparison = filename %>% str_remove("differential-proteins_") %>% str_remove(".csv")
  data = read_csv2(filename,skip=6, col_names=FALSE)
  DEPs = data[,c(2,5)]
  colnames(DEPs) = c("Sciomics ID", comparison)
  return(DEPs)
}

for (f in list.files(pattern = "differential-proteins")){
  this_DEPs = extract_DEP(f)
  raw_data_long = raw_data_long %>% left_join(this_DEPs,by="Sciomics ID")
}

# select only proteins that are DE in some comparison
filterout = function(df, ...) setdiff(df, filter(df, ...))
DE_data_long = raw_data_long %>% 
  filterout(across(PAT_Day10__Treated_vs_Untreated:Untreated_Day30__PAT_vs_WT, is.na))

# calculate average across replicates
DEP_means_long = DE_data_long %>% 
  group_by(protein_id,treatment,sample_group,timepoint) %>% 
  summarise_if(is.numeric, mean, na.rm = TRUE)

# lfc in each comparison
DEP_lfc = DEP_means_long %>% 
  ungroup() %>%
  dplyr::select(protein_id,PAT_Day10__Treated_vs_Untreated:Untreated_Day30__PAT_vs_WT) %>%
  group_by(protein_id) %>%
  summarize_all(unique)

# proteins DE in any patient treated vs control comparison
DE_treatment = DEP_lfc %>%
  rowid_to_column("row_name") %>%
  filter(!is.na(PAT_Day10__Treated_vs_Untreated) |
           !is.na(PAT_Day20__Treated_vs_Untreated) |
           !is.na(PAT_Day30__Treated_vs_Untreated)) %>% pull(row_name,protein_id)

# lfc in comparisons for heatmap
lfc_col_fun = colorRamp2(c(-3, 0, +3), c("blue", "white", "red"))
lfc_df = DEP_lfc%>%
  dplyr::select(protein_id,Untreated_Day10__PAT_vs_WT:Untreated_Day30__PAT_vs_WT,PAT_Day10__Treated_vs_Untreated:PAT_Day30__Treated_vs_Untreated) %>%
  column_to_rownames("protein_id") %>%
  as.data.frame()

##############
## HEATMAPS ##
##############
# ------------------------------------ #
# Scale and cluster without WT treated #
# ------------------------------------ #
DEP_means_wide_condition = DEP_means_long %>%
  ungroup() %>%
  mutate(condition = str_c(sample_group,treatment,sep="_")) %>%
  dplyr::select(protein_id,timepoint,condition, values) %>%
  pivot_wider( names_from = timepoint,values_from = values)

sc_sel_cond_mat_wide = DEP_means_wide_condition %>%
  filter(condition !="WT_CD") %>%
  dplyr::select(`10`:`30`) %>%
  apply(1, scale) %>% 
  t() 

colnames(sc_sel_cond_mat_wide) = c("10","20","30")
sc_sel_cond_mat_wide = sc_sel_cond_mat_wide %>% 
  as_tibble()
sc_sel_cond_mat_wide['protein_id'] = DEP_means_wide_condition %>% 
  filter(condition !="WT_CD") %>% 
  pull('protein_id')

sc_sel_cond_mat_wide['condition'] = DEP_means_wide_condition %>% 
  filter(condition !="WT_CD") %>% 
  pull('condition')

sc_sel_cond_mat_wider = sc_sel_cond_mat_wide %>%
  pivot_wider(protein_id, names_from = condition,values_from = `10`:`30`,names_glue = "{condition}_{.value}")

sc_sel_cond_mat = sc_sel_cond_mat_wider %>%
  column_to_rownames('protein_id') %>%
  as.matrix()

#-------------------#
# Intensity heatmap #
#-------------------#
set.seed(2021)

# annotation of columns in the main heatmap, showing scaled intensity
trend_ha = HeatmapAnnotation(
  foo = anno_block(gp = gpar(fill = 2:4),
                   labels = c("WT_UN", "Patient_UN", "Patient_CD"), 
                   labels_gp = gpar(col = "white", fontsize = 10)),
  type=c("patient","patient","wt","patient","patient","wt","patient","patient","wt"),
  treatment=c("treated","untreated","untreated","treated","untreated","untreated","treated","untreated","untreated"),
  timepoint=rep(c(10,20,30),each=3),
  col = list(type = c("wt"="cornflowerblue", "patient"="firebrick1"),
             treatment = c("untreated" = "white", "treated" = "#006E82"),
             timepoint  =colorRamp2(c(10, 20, 30), brewer.pal(n=3, name="Greens"))),
  annotation_name_side = "left")

# scaled intensity heatmap
trend_heatmap = Heatmap(sc_sel_cond_mat,
                 col = colorRamp2(c(1.5, 0, -1.5), brewer.pal(n=3, name="PuOr")),
                 column_order = sort(colnames(sc_sel_cond_mat)),
                 column_split = factor(str_sub(colnames(sc_sel_cond_mat),1,-4),levels=c("WT_UN","Patient_UN","Patient_CD")),
                 row_km = 8,
                 row_title = "Cluster %s",
                 row_title_gp = gpar(fontsize = 10),
                 row_title_rot = 0,
                 row_dend_reorder = TRUE,
                 show_row_names = FALSE,
                 show_column_names = FALSE,
                 column_title =NULL,
                 top_annotation =trend_ha,
                 name="Abundance (scaled)"
)

#-------------#
# lfc heatmap #
#-------------#
# annotation of the lfc heatmap
lfc_ha = HeatmapAnnotation(
  foo = anno_block(gp = gpar(fill = 5:6),
                   labels = c("Patient vs WT", "Treated vs. untreated"), 
                   labels_gp = gpar(col = "white", fontsize = 10)),
  timepoint=rep(c(10,20,30),2),
  col = list(timepoint=colorRamp2(c(10, 20, 30), brewer.pal(n=3, name="Greens"))),
  show_annotation_name  = FALSE)

lfc_mat = as.matrix(lfc_df[rownames(sc_sel_cond_mat),])
lfc_heatmap = Heatmap(lfc_mat,
                      col=colorRamp2(c(-3, 0, +3), c("blue", "grey", "red")),
                      na_col="white",
                      column_order = sort(colnames(lfc_mat)),
                      column_split = factor(str_split_fixed(colnames(lfc_mat),"_",2)[,1],levels=c("Untreated","PAT")),
                      show_row_names = FALSE,
                      show_column_names = FALSE,
                      column_title =NULL,
                      border = TRUE,
                      width=4,
                      name="LogFC",
                      top_annotation = lfc_ha)

# row annotation of interesting proteins
# select among protein id DE in treatment vs control, the ones that are interesting according to literature
interesting_protein_id = raw_data_long %>%
  filter(Name %in% interesting_proteins) %>%
  dplyr::select(protein_id,Name) %>%
  distinct() %>%
  filter(protein_id %in% names(DE_treatment))

ind = which(rownames(lfc_mat) %in% interesting_protein_id$protein_id)
ind_names = interesting_protein_id$Name[match(interesting_protein_id$protein_id, rownames(lfc_mat)[ind])]
row_an_interesting <- rowAnnotation(Proteins = anno_mark(at = ind, labels = ind_names,labels_gp=gpar(fontsize = 10)))

# combine all
final_heatmap = draw(trend_heatmap + lfc_heatmap + row_an_interesting)
final_heatmap = draw(final_heatmap)
pdf("proteomics_heatmap_lfc_interestingproteins_2.pdf")
draw(final_heatmap)
dev.off()

##############
## Clusters ##
##############
groups = row_order(final_heatmap)

row_group = cbind(unlist(groups),rep(names(groups),lengths(groups))) %>%
  as_tibble() %>%
  mutate(row_order=as.numeric(V1), row_group=paste0("cluster_",V2)) %>%
  arrange(row_order) %>% pull(row_group)

sc_sel_cond_mat_wider['row_group'] = row_group

sc_sel_cond_mat_wider %>% 
  dplyr::select(-protein_id) %>%
  pivot_longer(-row_group,names_to = c("condition","timepoint"), names_sep=-2, values_to="Scaled_abundance") %>%
  mutate(timepoint=as.numeric(timepoint)) %>%
  ggplot(aes(x=timepoint,y=Scaled_abundance,group=condition,color=condition)) + 
  geom_smooth(alpha=0.2,method="loess",se=FALSE) + 
  facet_wrap(~row_group) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "transparent", colour = NA), 
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  xlab("Time, days") +
  ylab("Protein abundance, scaled") +
  scale_x_continuous(breaks = c(10,20,30)) +
  scale_color_manual(name="Condition", 
                     values=c("WT_UN_"="cornflowerblue","Patient_UN_"="firebrick1","Patient_CD_"="#006E82"),
                     labels=c("Wt untreated","Patient untreated","Patient treated"),
                     breaks=c("WT_UN_","Patient_UN_","Patient_CD_"))

ggsave("cluster_behaviour_color_1.pdf", plot=last_plot(), dpi=320, width = 11.69 , height = 8.27, units = "in")

#########################
## Interesting proteins #
#########################
# boxplot/dotplot of all interesting proteins in each cluster (from literature): wt_un, p_un, p_cd
interesting_proteins_data = raw_data_long %>% 
  full_join(sc_sel_cond_mat_wider,on="protein_id") %>%
  filter(protein_id %in% interesting_protein_id$protein_id) %>%
  mutate(condition=paste0(sample_group,".",treatment),sample=paste0(condition,".",timepoint)) %>%
  filter(condition!="WT.CD") %>%
  dplyr::select(protein_id,condition,timepoint,values,row_group,sample) 
  
# "boxplot" (no box)
graph_list = list()
for (cluster in unique(interesting_proteins_data$row_group)){
  interesting_proteins_data_in_cluster = interesting_proteins_data %>%
    filter(row_group == cluster)
  
  if(nrow(interesting_proteins_data_in_cluster)>0){
    
    for (protein in unique(interesting_proteins_data_in_cluster$protein_id)){
      lines = interesting_proteins_data_in_cluster %>%
        filter(protein_id==protein) %>%
        group_by(condition,sample,timepoint) %>%
        summarize(q50=median(values),std=sd(values)) %>%
        ggplot() +
        # geom_boxplot(aes(x=sample,
        #                  ymin =q50-std, 
        #                  lower=q50,
        #                  middle=q50,
        #                  upper=q50,
        #                  ymax=q50+std,
        #                  group=sample), stat = "identity") +
        geom_errorbar(aes(x=sample,
                          ymin =q50-std, 
                          ymax=q50+std,
                          group=sample), width=0.5)
    
        points = interesting_proteins_data_in_cluster %>%
          filter(protein_id==protein) %>%
          geom_point(data=.,aes(x=sample,
                                y=values,
                                group=sample,
                                color=condition),
                     size=1)
      graph = lines + points + 
          ggtitle(interesting_protein_id %>% filter(protein_id == protein) %>% pull(Name))
      graph_list[[cluster]][[protein]] = graph
    }
  }
}

my_legend = get_legend(graph_list[[1]][[1]] +
                         scale_color_manual(name="Condition", 
                                            values=c("WT.UN"="cornflowerblue","Patient.UN"="firebrick1","Patient.CD"="#006E82"),
                                            limits=c("WT.UN","Patient.UN","Patient.CD")) +
                         scale_x_discrete(limits=c("WT.UN.10","WT.UN.20","WT.UN.30",
                                                   "Patient.UN.10","Patient.UN.20","Patient.UN.30",
                                                   "Patient.CD.10","Patient.CD.20","Patient.CD.30"),
                                          labels=rep(c(10,20,30),3)) +
                         theme_classic() +
                         theme(axis.title.x=element_blank(),
                               axis.title.y=element_blank(),
                               plot.title= element_text(hjust=0.5, size=8),
                               axis.text.x = element_text(size = 5),
                               axis.text.y = element_text(size = 5))
                         )

common_theme=
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title= element_text(hjust=0.5, size=8),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        legend.position = "none")
  

complete_graph = function(g){
  complete_g  = g +
    common_theme +
    scale_color_manual(name="Condition", 
                       values=c("WT.UN"="cornflowerblue","Patient.UN"="firebrick1","Patient.CD"="#006E82"),
                       limits=c("WT.UN","Patient.UN","Patient.CD")) +
    scale_x_discrete(limits=c("WT.UN.10","WT.UN.20","WT.UN.30",
                              "Patient.UN.10","Patient.UN.20","Patient.UN.30",
                              "Patient.CD.10","Patient.CD.20","Patient.CD.30"),
                     labels=rep(c(10,20,30),3))
  return(complete_g)
}

initial_graph_list = graph_list
graph_list = lapply(initial_graph_list,function(x) lapply(x,complete_graph))


p=ggdraw() +
  # 1st row: cluster_1
  draw_plot(graph_list[['cluster_1']][[1]], x = 0, y = .6, width = .15, height = .15) +
  draw_plot(graph_list[['cluster_1']][[2]], x = .15, y = .6, width = .15, height = .15) +
  draw_plot(graph_list[['cluster_1']][[3]], x = .3, y = .6, width = .15, height = .15) +
  draw_plot(graph_list[['cluster_1']][[4]], x = .45, y = .6, width = .15, height = .15) +
  draw_plot(graph_list[['cluster_1']][[5]], x = .6, y = .6, width = .15, height = .15) +
  draw_plot(graph_list[['cluster_1']][[6]], x = .75, y = .6, width = .15, height = .15) +
  # 2nd row: cluster_4 + cluster_8
  draw_plot(graph_list[['cluster_4']][[1]], x = 0, y = .45, width = .15, height = .15) +
  draw_plot(graph_list[['cluster_4']][[2]], x = .15, y = .45, width = .15, height = .15) +
  draw_plot(graph_list[['cluster_4']][[3]], x = .3, y = .45, width = .15, height = .15) +
  draw_plot(graph_list[['cluster_4']][[4]], x = .45, y = .45, width = .15, height = .15) +
  
  draw_plot(graph_list[['cluster_8']][[1]], x = .6, y = .45, width = .15, height = .15) +
  
  # 3rd row: cluster_3 + cluster_5 + cluster_7 + legend
  draw_plot(graph_list[['cluster_3']][[1]], x = 0, y = .3, width = .15, height = .15) +
  draw_plot(graph_list[['cluster_3']][[2]], x = .15, y = .3, width = .15, height = .15) +
  draw_plot(graph_list[['cluster_3']][[3]], x = .3, y = .3, width = .15, height = .15) +
  
  draw_plot(graph_list[['cluster_5']][[1]], x = .45, y = .3, width = .15, height = .15) +
  
  draw_plot(graph_list[['cluster_7']][[1]], x = .6, y = .3, width = .15, height = .15) +
  
  draw_plot(my_legend, x = .75, y = .3, width = .15, height = .15)

ggsave("Interesting_proteins_abundance_1.pdf", plot=p, dpi=320, width = 11.69 , height = 8.27, units = "in",device=cairo_pdf)

# these were missing before
p=ggdraw() +
  draw_plot(graph_list[['cluster_8']][[2]], x = 0, y = .45, width = .15, height = .15) +
  draw_plot(graph_list[['cluster_8']][[3]], x = .15, y = .45, width = .15, height = .15)

ggsave("Interesting_proteins_abundance_2.pdf", plot=p, dpi=320, width = 11.69 , height = 8.27, units = "in",device=cairo_pdf)



# other proteins
graph_list = list()

id = "ab0711" # RARB
id = "ab0056" # BIRC

this_id_data = raw_data_long %>% 
  full_join(sc_sel_cond_mat_wider,on="protein_id") %>%
  filter(`Sciomics ID` ==id) %>%
  mutate(condition=paste0(sample_group,".",treatment),sample=paste0(condition,".",timepoint)) %>%
  filter(condition!="WT.CD") %>%
  dplyr::select(protein_id,condition,timepoint,values,row_group,sample) 

for (protein in unique(this_id_data$protein_id)){
  lines = this_id_data %>%
    filter(protein_id==protein) %>%
    group_by(condition,sample,timepoint) %>%
    summarize(q50=median(values),std=sd(values)) %>%
    ggplot() +
    geom_errorbar(aes(x=sample,
                      ymin =q50-std, 
                      ymax=q50+std,
                      group=sample), width=0.5)
  
  points = this_id_data %>%
    filter(protein_id==protein) %>%
    geom_point(data=.,aes(x=sample,
                          y=values,
                          group=sample,
                          color=condition),
               size=1)
  graph = lines + points + 
    ggtitle(raw_data_long %>% filter(protein_id == protein) %>% pull(Name)%>% unique())
  graph_list[[protein]] = graph
}

glist =lapply(graph_list,complete_graph)

p=ggdraw() +
  draw_plot(glist[[1]], x = 0, y = .45, width = .15, height = .15)

ggsave("RARB_abundance.pdf", plot=p, dpi=320, width = 11.69 , height = 8.27, units = "in",device=cairo_pdf)



#####
# improve format of DEP analysis tables
####
simplify_DEP <- function(filename){
  comparison = filename %>% str_remove("differential-proteins_") %>% str_remove(".csv")
  data = read_csv2(filename,skip=6, col_names=FALSE)
  colnames(data) = strsplit("Protein;AntibodyID;Uniprot-Entry-Name;UniprotID;logFC;AveExp;P-Value;adj.P.Val;HGNC;HGNC.alternatives",";")[[1]]
  write.csv(data,sprintf("diff_expr_proteins_%s.csv",comparison),row.names=FALSE)
}

sapply(list.files(pattern = "differential-proteins"), simplify_DEP)
