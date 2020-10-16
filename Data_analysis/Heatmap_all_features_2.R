library(data.table)
library(tidyverse)
library(janitor)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

###################
### User input ####
###################
###################
##LOADING OF DATA #
###################
Table_diff_percent_old = read.table('P:/Scripts/R/GGPlot_scripts/Heatmaps_all_features/TablesAll_diff_percent_old_2.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
Table_diff_percent_new = read.table('P:/Scripts/R/GGPlot_scripts/Heatmaps_all_features/TablesAll_diff_percent_new_2.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
Table_astro_percent = read.table('P:/Scripts/R/GGPlot_scripts/Heatmaps_all_features/TablesAll_astro_percent_2.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
Table_neuromorph = read.table('P:/Scripts/R/GGPlot_scripts/Heatmaps_all_features/TablesAll_neuromorph_2.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
Table_SNCA = read.table('P:/Scripts/R/GGPlot_scripts/Heatmaps_all_features/TablesAll_SNCA_2.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
Table_prolif = read.table('P:/Scripts/R/GGPlot_scripts/Heatmaps_all_features/TablesAll_prolif_2.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
Table_apop = read.table('P:/Scripts/R/GGPlot_scripts/Heatmaps_all_features/TablesAll_apop_2.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)

###################################
# Making columns names consistent #
###################################
compareddf =compare_df_cols(Table_diff_percent_old,Table_diff_percent_new,Table_astro_percent,Table_neuromorph,Table_SNCA,Table_prolif,Table_apop)
Table_neuromorph %>% 
  rename(Column = COL,field =Field, Row =ROW, TimePoint = Timepoint ) -> Table_neuromorph

Table_diff_percent_new$Category[Table_diff_percent_new$Category=="Controls"]<- "Control"
Table_astro_percent$Category[Table_astro_percent$Category=="Controls"]<- "Control"
Table_SNCA$Category[Table_SNCA$Category=="Controls"]<- "Control"
Table_neuromorph$Category[Table_neuromorph$Category=="Controls"]<- "Control"

myIDnew = function(x){
  tmp = do.call(rbind,strsplit(x$Barcode,"_"))
  x$Replicate = substr(tmp[,9],2,2)
  x$Handle = paste0(x$Category,"_",x$Timepoint,"_",x$Replicate)
  toremove = c("Row","Column","Barcode","Timepoint",             
               "Category","field","Replicate")
  x %>% 
  dplyr::filter(Category!="Corrected") %>% 
  rename(Handle_ID = Handle) %>% 
  dplyr::select(-toremove)->toto
}
myIDold = function(x){
  tmp = do.call(rbind,strsplit(x$Barcode,"_"))
  x$Replicate = substr(tmp[,7],2,2)
  x$Handle = paste0(x$Category,"_",x$TimePoint,"_",x$Replicate)
  toremove = c("Row","Column","Barcode","TimePoint",             
               "Category","field","Replicate")
  x %>% 
    rename(Handle_ID = Handle) %>% 
    dplyr::select(-toremove)->toto
}

Table_diff_percent_old_ID = myIDold(Table_diff_percent_old)
Table_neuromorph_ID = myIDold(Table_neuromorph)
Table_prolif_ID = myIDold(Table_prolif)
Table_apop_ID = myIDold(Table_apop)

Table_diff_percent_new_ID = myIDnew(Table_diff_percent_new)
Table_astro_percent_ID = myIDnew(Table_astro_percent)
Table_SNCA_ID = myIDnew(Table_SNCA)

###################
# Select features #
###################
#that are interesting from each of them ->avoid that there are reapeated columns between tables
compareddf_ID =compare_df_cols(Table_diff_percent_old_ID,Table_diff_percent_new_ID,Table_astro_percent_ID,Table_neuromorph_ID,Table_SNCA_ID,Table_prolif_ID,Table_apop_ID)

tr_diff_new = c("GFAPbyGFAPVol","GFAPbyNucVol","GFAPVol","GFAPVolByNucVol","GFAPVolByTHVol","NucVol","THbyNucVol","THbyTHVol","THVol","THVolByNucVol","THVolByTuj1Vol","Tuj1byNucVol","Tuj1byTuj1Vol","Tuj1Vol","Tuj1VolByNucVol")
tr_snca = c("GFAPbyGFAPVol","GFAPbyNucVol","GFAPVol","GFAPVolByNucVol","GFAPVolByTHVol","NucVol","THbyNucVol","THbyTHVol","THVol","THVolByNucVol","SynbySynVol","SynbyNucVol")
tr_prol =c("NucArea","THArea","THAreaByNucArea","THbyNucArea","THbyTHArea","Ki67byNucArea") 
tr_apop =c("NucArea","THArea","THAreaByNucArea","THbyNucArea","THbyTHArea","PARPbyNucArea")
tr_astro_percent = c("NucVol","THbyNucVol","THbyTHVol","THVol","THVolByNucVol","GFAPbyGFAPVol","GFAPbyNucVol")
tr_neuromorph = c("NucVol","THArea","MinTuj1Area","MinTHArea","Tuj1Area")
tr_diff_old = c("THbyNucVol","THbyTHVol","Tuj1byTuj1Vol","Tuj1byNucVol")

myrm = function(x,tr){
  x %>% 
    dplyr::select(-tr)-> x
}
Table_diff_percent_new_ID = myrm(Table_diff_percent_new_ID,tr_diff_new)
Table_SNCA_ID = myrm(Table_SNCA_ID,tr_snca)
Table_prolif_ID = myrm(Table_prolif_ID,tr_prol)
Table_apop_ID = myrm(Table_apop_ID,tr_apop)
Table_astro_percent_ID = myrm(Table_astro_percent_ID,tr_astro_percent)
Table_neuromorph_ID = myrm(Table_neuromorph_ID,tr_neuromorph)
Table_diff_percent_old_ID = myrm(Table_diff_percent_old_ID,tr_diff_old)

#######################
## Preparing HM table #
#######################

list(Table_diff_percent_old_ID,Table_neuromorph_ID,Table_prolif_ID,Table_apop_ID,Table_diff_percent_new_ID,Table_astro_percent_ID,Table_SNCA_ID) %>% 
  purrr::reduce(left_join, by = "Handle_ID") ->Table_All_HM

Table_All_HM = na.omit(Table_All_HM)

tmp = do.call(rbind,strsplit(Table_All_HM$Handle_ID,"_"))
Table_All_HM$Category = tmp[,1]
Table_All_HM$Timepoint = tmp[,2]
Table_All_HM$Replicate = tmp[,3]

Table_All_HM_data = Table_All_HM[,c(1:86)]
Table_All_HM_metadata = Table_All_HM[,c(87:89)]

Table_All_HM_data = as.data.frame(Table_All_HM_data)
toto = data.frame(Table_All_HM_data[,-7],row.names = Table_All_HM_data[,7])

scaled = scale(toto)
scaled_t = t(scaled)

#Annotations
mydataframeClass = data.frame(c("TH_diff","Tuj1_diff","TH_diff","Tuj1_diff","TH_diff","TH_diff",
                           rep("Tuj1_morphology",23),rep("TH_morphology",24),"Proliferation","Proliferation",rep("Apoptosis",4),
                           "Astrocyte_propor","TH_diff","TH_diff","TH_diff",rep("Astrocyte_propor",12),"SNCA_TH",
                           "SNCA_TH","SNCA_astro","SNCA_TH","SNCA_TH","SNCA_astro","SNCA_astro","SNCA_TH",
                           "SNCA_TH","SNCA_astro"))

myscalecol = circlize::colorRamp2(c(7, 14, 21), brewer.pal(n=3, name="Greens"))
myannot = HeatmapAnnotation(
  Timepoint = as.numeric(Table_All_HM_metadata$Timepoint),
  Category =Table_All_HM_metadata$Category,
  col = list(Category = c("Control" = "cornflowerblue","Patient" = "firebrick1"),
             Timepoint = myscalecol))

#######
## HM #
#######
Heatmap(scaled_t,  col = colorRamp2(c(2, 0, -2), brewer.pal(n=3, name="PuOr")),
        name = "Value", #title of legend
        column_title = "Samples", row_title = "Variables",
        row_names_gp = gpar(fontsize = 4),
        column_names_gp = gpar(fontsize = 2),
        column_dend_reorder =TRUE,
        top_annotation = myannot,
        row_dend_reorder = TRUE,
        clustering_distance_columns =  "canberra",
        row_split = mydataframeClass,
        column_split = Table_All_HM_metadata[,2],
        heatmap_width = unit(8,"cm"),
        heatmap_height = unit(10,"cm"))

