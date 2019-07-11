
library(tidyverse)

#ディレクトリ指定の部分
Directry= "MTBLS40vsMTBLS45"


# ID読み込みの部分
files = list.files(Directry)
filename_1 = files[stringr::str_detect(files, "GetPubChem")][1]
file_1 = read_csv(file = paste(Directry, filename_1, sep = "/"), col_names = FALSE)
IDs_1 = as.character(file_1$X2[!is.na(file_1$X2)])
filename_2 = files[stringr::str_detect(files, "GetPubChem")][2]
file_2 = read_csv(file = paste(Directry, filename_2, sep = "/"), col_names = FALSE)
IDs_2 = as.character(file_2$X2[!is.na(file_2$X2)])



#整形の部分
VennList = list(IDs_1, IDs_2)
##pathlist読み込み
pathlist = readxl::read_excel("pathlist.xlsx")
##filenameと合うdatasetnameを探す
dsname_1 = pathlist$datasetname[pathlist$GetPubChem == filename_1]
dsname_2 = pathlist$datasetname[pathlist$GetPubChem == filename_2]
names(VennList) = c(dsname_1, dsname_2)


#描く部分
library(VennDiagram)
VennDiagram::venn.diagram(VennList, filename = paste(Directry, "/", Directry, ".png", sep = ""))

##tableの読み込み
IDlist = read_csv("IDlist.csv", col_types = cols("c","c","c","c","c","c","c","c","c","c","c","c","c","c","c"))

#差を求める部分
IDs1_Only = data.frame(setdiff(IDs_1, IDs_2))
IDs1_Only$setdiff.IDs_1..IDs_2. = as.character(IDs1_Only$setdiff.IDs_1..IDs_2.)
IDs1OnlyAndNames = inner_join(IDs1_Only, IDlist, by = c("setdiff.IDs_1..IDs_2." = "PubChem"))
IDs1OnlyAndNames = IDs1OnlyAndNames[,c("setdiff.IDs_1..IDs_2.", "Match")]

conv.id1 <- IDs1OnlyAndNames$Match
names(conv.id1) <- IDs1OnlyAndNames$setdiff.IDs_1..IDs_2.

IDs1OnlyAndNamesUnique = conv.id1[IDs1_Only$setdiff.IDs_1..IDs_2.]

IDs2_Only = data.frame(setdiff(IDs_2, IDs_1))
IDs2_Only$setdiff.IDs_2..IDs_1. = as.character(IDs2_Only$setdiff.IDs_2..IDs_1.)
IDs2OnlyAndNames = inner_join(IDs2_Only, IDlist, by = c("setdiff.IDs_2..IDs_1." = "PubChem"))
IDs2OnlyAndNames = IDs2OnlyAndNames[,c("setdiff.IDs_2..IDs_1.", "Match")]

conv.id2 <- IDs2OnlyAndNames$Match
names(conv.id2) <- IDs2OnlyAndNames$setdiff.IDs_2..IDs_1.

IDs2OnlyAndNamesUnique = conv.id2[IDs2_Only$setdiff.IDs_2..IDs_1.]


#積を求める部分
IDs_common = data.frame(intersect(IDs_2, IDs_1))
IDs_common$intersect.IDs_2..IDs_1. = as.character(IDs_common$intersect.IDs_2..IDs_1.)
IDsCommonAndNames = inner_join(IDs_common, IDlist, by = c("intersect.IDs_2..IDs_1." = "PubChem"))
IDsCommonAndNames = IDsCommonAndNames[,c("intersect.IDs_2..IDs_1.", "Match")]

conv.idcommon = IDsCommonAndNames$Match
names(conv.idcommon) <- IDsCommonAndNames$intersect.IDs_2..IDs_1.

IDsCommonAndNamesUnique = conv.idcommon[IDsCommonAndNames$intersect.IDs_2..IDs_1.]



#IDの結果を書き込む部分
openxlsx::write.xlsx(list(data.frame(IDs1OnlyAndNamesUnique), data.frame(IDs2_Only),data.frame(IDs_common)),file  =paste(Directry, "/", Directry, ".xlsx", sep = ""),row.names=TRUE  )


