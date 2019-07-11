
library(tidyverse)

GMDandPRIMe_raw = read_lines("GMDandPRIMe.msp")
GMDandPRIMe_names = GMDandPRIMe_raw[stringr::str_detect(GMDandPRIMe_raw, "Name:") & !stringr::str_detect(GMDandPRIMe_raw, "Synonym:")]
GMDandPRIMe_names = str_remove(GMDandPRIMe_names, "Name: ")
GMDandPRIMe_names = str_remove(GMDandPRIMe_names, "Name:")

write.csv(GMDandPRIMe_names, "GMDandPRIMe_names.csv", row.names =FALSE)

library(stringr)
#命名規則は以下の３つ
#１．M00なんたらかんたら
GMDandPRIMe_M = str_subset(GMDandPRIMe_names, "M[0-9]{6}_A[0-9]{6}-101-xxx_NA_[+-]?([0-9]*[.])?[0-9]+_[:alnum:]{4,5}_[:alnum:]{4,5}_[:alnum:]{3}_")
#2594
#２．GCTMSなんたらかんたら
GMDandPRIMe_G = str_subset(GMDandPRIMe_names, "GCTMS_.*")
#237
#３．数字(整数４桁、時々少数第一位１桁まで)
GMDandPRIMe_N = str_subset(GMDandPRIMe_names, "^[+-]?([0-9]*[.])?[0-9]+_.*")
#99

#いらんとこ消そう
GMDandPRIMe_M_names = str_remove(GMDandPRIMe_M , pattern = "M[0-9]{6}_A[0-9]{6}-101-xxx_NA_[+-]?([0-9]*[.])?[0-9]+_[:alnum:]{4,5}_[:alnum:]{4,5}_[:alnum:]{3}_")
GMDandPRIMe_G_names = str_sub(str_remove(GMDandPRIMe_G, pattern = "GCTMS_N12C_STD_[0-9]+_[+-]?([0-9]*[.])?[0-9]+_[+-]?([0-9]*[.])?[0-9]+_"), start = 2, end = -2)
GMDandPRIMe_N_names = str_remove(GMDandPRIMe_N, pattern = "^[+-]?([0-9]*[.])?[0-9]+_")

#まとめる
tmpnames = c(GMDandPRIMe_M_names, GMDandPRIMe_G_names, GMDandPRIMe_N_names)

#誘導体(TMSとかMeOXとか)を消す
MetNames = tmpnames %>%
str_remove(pattern ="[:punct:][:digit:]*TMS[:punct:]") %>%  
str_remove(pattern ="[:punct:].?MeOX[:punct:]") %>%  
str_remove(pattern ="[:punct:].?MEOX[:punct:]") %>%
str_remove(pattern ="[:punct:].?MEOX[:punct:]") %>%
str_remove(pattern ="BP") %>%
str_remove(pattern ="BP1") %>%
str_remove(pattern ="BP2") %>%
str_remove(pattern ="MP") %>%
# 消しといたほうがいいやつ
str_remove(pattern =fixed("(1TFA)")) %>%
  str_remove(pattern =fixed("[+CO2]")) %>%
str_remove(pattern =fixed("[-H2O]")) %>%
str_remove(pattern =fixed("[-H20]")) %>%
str_remove(pattern =fixed("(TMS)")) %>%
str_remove(pattern =fixed("[-NH3]")) %>%
str_remove(pattern =fixed("(7TMS)")) %>%
str_remove(pattern =fixed("[-O]")) %>%
str_remove(pattern =fixed("[+H2]")) %>%
str_remove(pattern ="_[:alnum:]*TMS>") %>%
str_remove(pattern ="-[:digit:]*HCl>") %>%
str_remove(pattern =fixed("_1MeOX") )%>%
str_remove(pattern =fixed("_1MeO") )%>%
str_remove(pattern =fixed("(1 MEOX)") )%>%
str_remove(pattern ="\"") %>%
str_remove(pattern =fixed("(+-)-")) %>% 
str_remove(pattern = fixed("(+/-)-")) %>%
str_remove(pattern =fixed("(+/-)-")) %>%
str_remove(pattern = fixed("(FAME MIX)")) %>%
str_remove(pattern = fixed("(Derivate not found)")) %>%
str_remove(pattern = "_[:digit:]?TMS") %>%
  
  
str_trim(side = "right")
#D-Glucose #<- Ok
#D(-)Fructose	#<-アカン
#D(+)-Mannose	#<-Ok
#D-(+)-Glucose アカン
#D-Glucose　Ok
#Glucose　D-GlucoseとしてOk
write.csv(MetNames, file = "MetNames.csv", row.names = FALSE)


NamesTable = data.frame(GMDandPRIMe_names, MetNames)
NamesTable$GMDandPRIMe_names = as.character(NamesTable$GMDandPRIMe_names)
NamesTable$MetNames = as.character(NamesTable$MetNames)
write.csv(NamesTable, "NamesTable.csv", row.names = FALSE)

#要注意
#Glycinamide_1 などの糖(アミノ)の後ろについた数字





#MetaboAnalystのインストール
#参考:https://github.com/xia-lab/MetaboAnalystR
metanr_packages <- function(){
  metr_pkgs <- c("Rserve", "ellipse", "scatterplot3d", "Cairo", "randomForest", "caTools", "e1071", "som", "impute", "pcaMethods", "RJSONIO", "ROCR", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "pheatmap", "SSPA", "sva", "Rcpp", "pROC", "data.table", "limma", "car", "fitdistrplus", "lars", "Hmisc", "magrittr", "methods", "xtable", "pls", "caret", "lattice", "igraph", "gplots", "KEGGgraph", "reshape", "RColorBrewer", "tibble", "siggenes", "plotly", "xcms", "CAMERA", "fgsea", "MSnbase", "BiocParallel", "metap", "reshape2", "scales")
  list_installed <- installed.packages()
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  if(length(new_pkgs)!=0){
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(new_pkgs, version = "3.9")
    print(c(new_pkgs, " packages added..."))
  }
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
}

metanr_packages()

install.packages("pacman")
library(pacman)
pacman::p_load(Rserve, ellipse, scatterplot3d, Cairo, randomForest, caTools, e1071, som, impute, pcaMethods, RJSONIO, ROCR, globaltest, GlobalAncova, Rgraphviz, preprocessCore, genefilter, pheatmap, SSPA, sva, Rcpp, pROC, data.table, limma, car, fitdistrplus, lars, Hmisc, magrittr, methods, xtable, pls, caret, lattice, igraph, gplots, KEGGgraph, reshape, RColorBrewer, tibble, siggenes, plotly, xcms, CAMERA, fgsea, MSnbase, BiocParallel, metap, reshape2, scales)


install.packages("devtools")
library(devtools)
devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual", "--no-build-vignettes"))
library(MetaboAnalystR)

mSet<-InitDataObjects("NA", "utils", FALSE)
cmpd.vec<- MetNames
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "name", T, T, T, T, T);
mSet<-CreateMappingResultTable(mSet)

MetaboAnalystResult = data.frame(mSet$dataSet$map.table[,!colnames(mSet$dataSet$map.table)=="Comment"])
MetaboAnalystResult$Query = as.character(MetaboAnalystResult$Query)
MetaboAnalystResult$Match = as.character(MetaboAnalystResult$Match)
MetaboAnalystResult$HMDB = as.character(MetaboAnalystResult$HMDB)
MetaboAnalystResult$PubChem = as.character(MetaboAnalystResult$PubChem)
MetaboAnalystResult$ChEBI = as.character(MetaboAnalystResult$ChEBI)
MetaboAnalystResult$KEGG = as.character(MetaboAnalystResult$KEGG)
MetaboAnalystResult$METLIN = as.character(MetaboAnalystResult$METLIN)
MetaboAnalystResult$SMILES = as.character(MetaboAnalystResult$SMILES)


write.csv(MetaboAnalystResult, "MetaboAnalystResult.csv", row.names=TRUE)


# merge namesTable and MetaboAnalystResult
FullNamesAndMetaboAnalystResult =  cbind(NamesTable, MetaboAnalystResult[-1])


#aono ga yatta msp kara totta name to KEGG no list kara yatta metaboanalyst no kekka to merge

IDlist_kegg <- read.csv("IDlist_kegg.csv", stringsAsFactors = FALSE)
IDlist_kegg = data.frame( IDlist_kegg[,!colnames(IDlist_kegg)=="ID" & !colnames(IDlist_kegg)=="Query"])
IDlist_kegg$PubChem = as.character(IDlist_kegg$PubChem)
IDlist_kegg$ChEBI = as.character(IDlist_kegg$ChEBI)
IDlist_kegg$METLIN = as.character(IDlist_kegg$METLIN)


FinalMatrix =  left_join(IDlist_kegg[,2:8] %>% gather(key, val1) %>% rowid_to_column(),
          FullNamesAndMetaboAnalystResult[,3:9] %>% gather(key, val2) %>% rowid_to_column(),
          by = c("rowid", "key")) %>% 
  mutate(val = val1,
         val = ifelse(is.na(val1), val2, val),
         rowid = (rowid - 1) %% 2930) %>% 
  select(-c(val1, val2)) %>% 
  spread(key, val) %>% 
  select(-rowid)

FinalMatrix = cbind(FullNamesAndMetaboAnalystResult[,1:2], FinalMatrix, IDlist_kegg[,c(1,9:14)])



write.csv(FinalMatrix, "FinalMatrix.csv", row.names = FALSE)