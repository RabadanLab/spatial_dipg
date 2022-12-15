rm(list=ls())


############################################################################################################
#curated_data# 2022-10-27###################################################################################
############################################################################################################

root <- "/Users/cswc/Documents/spatial_dipg/"
# [1] annot = "column_annotation.txt" 
#   tma1 = "TMZ273_CellbyCell.xlsx" 
#   rownames(annot) = colnames(tma1) 
#   column 'feature' : feaure names
#   column 'feature_sub' : antibody-specific parameters e.g. min/ max/ mean/sd (for "Nucleus", range/ sum)

# [2] summary = "summary.txt"
#  summary statististic for tma1
#  column 'core' = tma1$TMA.core  #'TMA.core' has two distinct subset 'PathCellObject' & 'VEGFR2' in 'Name' columns (meaning?)
#  column 'count'= number of tma1 rows, in other words, number of cells within the core

#library(dplyr); library(reshape2)
#sum_df <- data.frame(dcast(as.data.frame(summary %>% group_by(sector, sector_sub) %>% count()), 
#                           sector ~ sector_sub))
#summary_df <- left_join(data.frame(sum_df, sum=rowSums(sum_df[,-1], na.rm=T))[,c("sector", "sum")],
#                        as.data.frame(summary %>% group_by(sector) %>% summarize(sum=sum(count))),
#                        by=c("sector"="sector"))
#colnames(summary_df)<- c("sector", "no_sub_sector", "no_cells")
#summary_df
#summry_df >> column 'sector'; column 'no_sub_sector'; column 'no_cells

# [3] DIR> 'data_transformed'
# 'tma1' original dataset has been splitted into 17 independent files according to TMA.core ('sector')

#NEXT STEP#
#1. check the meaning of column 'Name' ('PathCellObject' & 'VEGFR2')
#2. pathologic annotation for each TMA.core
#3. check which feature should be used for thresholding the antibody staining status 
#   -->> select the splitted dataset with subsetted feature sets

############################################################################################################
root<- "/Users/cswc/Documents/spatial_dipg/"

library(openxlsx)
tma1 <- read.xlsx(paste0(root, "TMA273_CellbyCell.xlsx"), sheet=1)
############################################################################################################
#annotation_file for column names
############################################################################################################
columns <- colnames(tma1)
#1_Image/ROI
#790_NAME
#2_Class
#404_TMA.core/ Parent
annot <- data.frame(pri=sapply(columns, function(x){ strsplit(x, ":.", fixed=T)[[1]][1]}),
           sec=sapply(columns, function(x){ strsplit(x, ":.", fixed=T)[[1]][2]}))

for(i in 1:nrow(annot)){
  x <- strsplit(as.character(annot[i,2]), ".", fixed=T)[[1]]
  annot[i,3] <- if(length(x)==2 & annot[i,2] %in% c("Max.caliper", "Min.caliper") ){ paste0(x[1], "_", x[2])
                  }else{ x[1]}
  annot[i,4] <- if(length(x) ==2){x[2]}else if(length(x) >=3){ x[3] }else{NA  }
  }

for(i in 1:nrow(annot)){
  list <- c("Area", "Perimeter", "Circularity", "Max_caliper", "Min_caliper", "Eccentricity", NA)
  annot[i,"feature"] <- if(annot[i,3] %in% list){annot[i,3]
                            }else if(annot[i,3] %!in% list & annot[i,4] %in% c("std", "dev")){paste0(annot[i,3], "_sd")}else{ paste0(annot[i,3], "_", annot[i,4])}
}
annot <- data.frame(annot, feature_sub= do.call(rbind, lapply(annot$feature, function(x){ strsplit(as.character(x), "_", fixed=T)[[1]][2]})))
annot[which(annot$feature_sub == "caliper"), "feature_sub"] <- NA

#feature_counts
#[1]Cell_90; [2]Nucleus_132; [3]Cytoplassm_84; [4]Nucleus/Cell.area.ratio
morpho <- c("Area", "Perimeter", "Circularity", "Max_caliper", "Min_caliper", "Eccentricity")

data.frame(Cell = data.frame(annot[which(annot$pri == "Cell" & annot$feature %!in% morpho), ] %>% group_by(V3) %>% count()),
           Cytoplasm = data.frame(annot[which(annot$pri == "Cytoplasm" & annot$feature %!in% morpho), ] %>% group_by(V3) %>% count())[,2],
           Nucleus = data.frame(annot[which(annot$pri == "Nucleus" & annot$feature %!in% morpho), ] %>% group_by(V3) %>% count())[,2])
annot[which(annot$pri == "Nucleus" & annot$feature %!in% morpho), "feature_sub"] %>% unique(.) #mean; sd; mas; min // range; sum (in Nucleus)

write.table(annot, paste0(root, "column_annotation.txt"), sep="\t", col.names=T, row.names=T, quote=F)


############################################################################################################
#summary_statistics for tma1
############################################################################################################
#790_NAME
#2_Class
#404_TMA.core/ Parent

Name <- unique(tma1$Name)
unique(do.call(rbind, lapply(Name, function(x) { strsplit(as.character(x), " - ", fixed=T)[[1]][2]}))) #PathCellObject; VEGFR2 ?????
all(unique(do.call(rbind, lapply(Name, function(x) { strsplit(as.character(x), " - ", fixed=T)[[1]][1]}))) == unique(tma1$TMA.core))

summary <- as.data.frame(tma1 %>% group_by(TMA.core) %>% count())
colnames(summary) <- c("core", "count")
summary <- data.frame(summary, sector=do.call(rbind, lapply(summary$core, function(x){ strsplit(as.character(x), "-", fixed=T)[[1]][1]})))
summary <- data.frame(summary, sector_sub=do.call(rbind, lapply(summary$core, function(x){ strsplit(as.character(x), "-", fixed=T)[[1]][2]})))


library(dplyr); library(reshape2)
sum_df <- data.frame(dcast(as.data.frame(summary %>% group_by(sector, sector_sub) %>% count()), 
      sector ~ sector_sub))
summary_df <- left_join(data.frame(sum_df, sum=rowSums(sum_df[,-1], na.rm=T))[,c("sector", "sum")],
          as.data.frame(summary %>% group_by(sector) %>% summarize(sum=sum(count))),
          by=c("sector"="sector"))
colnames(summary_df)<- c("sector", "no_sub_sector", "no_cells")
summary_df
write.table(summary, paste0(root, "summary.txt"), sep="\t", col.names=T, row.names=F, quote=F)

summary <- read.delim(file=paste0(root, "summary.txt"), sep="\t", header=T)


############################################################################################################
#data splitted by TMA.core
############################################################################################################
sector <- LETTERS[1:17]
sub_root <- paste0(root, "data_transformed/")

for(i in 1:length(sector)){
 sub_tma <- tma1[which( tma1$TMA.core %in% summary[which(summary$sector == sector[i]), "core"]), ]
 write.table(sub_tma, paste0(sub_root, "tma_", sector[i], ".txt"), sep="\t", col.names=T, row.names=F, quote=F)
}


for(i in 1:length(sector)){
sim <- read.delim(paste0(sub_root, "tma_", sector[i], ".txt"), sep="\t", header=T)
print(dim(sim)[1] == summary_df[which(summary_df$sector == sector[i]), "no_cells"])
} #SHOULD BE ALL TRUE

############################################################################################################
#check the annotation file#
############################################################################################################
root<- "/Users/cswc/Documents/spatial_dipg/"

list.files(root)

library(openxlsx)
clinical <- read.xlsx(paste0(root, "TMA_Analysis_all_Robyn_JN_2022-Nov.xlsx"), sheet=1)

#EDA for clinical annotation
'eda' <- function(x){ df <- data.frame(class = class(x),
                                       length = length(unique(x)),
                                       example = paste(head(unique(x)), collapse=";"))
             
                      return(df)}


sum_col <- grep("Sum", colnames(clinical))
eda_annot <- do.call(rbind, lapply(clinical[, -sum_col], function(x){ eda(x)}))

write.table(eda_annot, file=paste0(root, "eda_annot.txt"), sep="\t", col.names=T, row.names=T, quote=F)


####structure
sim <- clinical[,1:13]; colnames(sim)[c(1,5:7)] <- c("spot", "path", "H3K27M_score", "path_re")
dup_spot <- unique(as.data.frame(sim %>% count(spot) %>% subset(., n>1))$spot)
length(unique(dup_spot))

for(i in 1:10){
  x <- dup_spot[i]
  print(sim[which(sim$spot == x),c(1:7,11:13)])
}

colnames(clinical)[c(1:5)]
#####detailed inspection
###column_1 'Spot.#'
unique(do.call(rbind, lapply(clinical$`Spot.#`, function(x){ strsplit(x, "-")[[1]][2]})))

###column_2 'PID'
unique(clinical$PID[which(nchar(clinical$PID) >4)]) #6 PID annotated without "number" --> which patients? (spot )
unique(clinical$PID[which(nchar(clinical$PID) <=4)]) #74 patients annotated with "number"

clinical[which(clinical$PID %in% unique(clinical$PID[which(nchar(clinical$PID) >4)]) ),]
#UMPED67/ UMPED69/ PIDz101 are okay --> might be PATIENTS?
#T20-142 ? Frontal ? 

###column_3 'Diagnosis'
unique(clinical$Diagnosis)
clinical[which(is.na(clinical$Diagnosis) ==TRUE ),] #PID == 1418 lacking diagnosis but HAS DATA with pathologic annotation
clinical[which(colnames(clinical) %!in% c("Sum"))] %>% count(PID, Diagnosis) #%>% summarise(sum=sum(n))#single diagnosis per patient 

###column_4 'Location'
unique(clinical$Location)
#'Tumor'; "Normal" --> not 'location'
#'Brainstem'; 'Cortex' --> too general
#'CereTumor'; 'Cere Tumor'; 've' --> wrong annotation?

key <- c("Cere Tumor", "ve", "CereTumor")
clinical[which(clinical$Location %in% key),]
PIDs <- unique(clinical[which(clinical$Location %in% key), "PID"])
clinical[which(clinical$PID %in% PIDs),]
#options(max.print= 10000)

###column_5
unique(clinical$`Path.T:N`) #T; N; ?; numbers..
clinical[which(clinical$`Path.T:N` %in% c("?")),] #clear --> pathologic score was defined for this case
clinical[which(clinical$`Path.T:N` %!in% c("T", "N")),]
#for these cases, column_5 %!in% c("T", "N") & is.na(column_7) --> might be problematic! for classification!

###column_6
unique(clinical$Path.H3K27M.score) #NA, N, P, 0
clinical[which(clinical$Path.H3K27M.score == "0"),1:13] %>% count(`Path.T:N`, Path.virtual.H7E.Score)
clinical[which(is.na(clinical$Path.H3K27M.score) == TRUE),1:13] %>% count(`Path.T:N`, Path.virtual.H7E.Score)

###column_7
unique(clinical$Path.virtual.H7E.Score) #"M"; "B"; "?"; "CB"; "X"
key <- "?"
clinical[which(clinical$Path.virtual.H7E.Score == key),1:13] %>% count(`Path.T:N`, Path.H3K27M.score)

####sim
sim <- clinical[,1:13]
for(i in seq(5,7,1)){print(unique(sim[,i]))}

sim$cri_1 <- sim$`Path.T:N`; sim$cri_2 <- sim$Path.H3K27M.score; sim$cri_3 <- sim$Path.virtual.H7E.Score

sim$cri_1 <- ifelse(sim$cri_1 %in% "?", "?",
                    ifelse(sim$cri_1 == "T", "T",
                           ifelse(sim$cri_1 == "N", "N", "others")))
#xtabs(~cri_1, sim); xtabs(~`Path.T:N`, sim)

for(i in 1:nrow(sim)){
  x <- sim[i, "cri_2"]
  sim[i, "cri_2"] <- ifelse(is.na(x) == TRUE, NA, ifelse(x %in% c("N", "P"), "others", x))
}
#xtabs(~Path.H3K27M.score, sim); xtabs(~cri_2, sim)

for(i in 1:nrow(sim)){
  x <- sim[i, "cri_3"]
  sim[i, "cri_3"] <- ifelse( x %in% c("?"), "Y", ifelse(x %!in% c("T", "N", "N+"), "others", x))
}
xtabs(~Path.virtual.H7E.Score, sim); xtabs(~cri_3, sim)

sim %>% count(cri_1, cri_2, cri_3)




#
clinical <- clinical[,1:7]
colnames(clinical) <- clinical[2,]
clinical <- clinical[-c(1,2),]

colnames(clinical) <- c("Image", "TMA_core", "patients", "disease", "noidea", "x_coord", "y_coord")
data.frame(patient=clinical$patients, 
            pri=do.call(rbind, lapply(clinical$patients, function(x){ strsplit(as.character(x), " ", fixed=T)[[1]][1]})),
            sec=do.call(rbind, lapply(clinical$patients, function(x){ strsplit(as.character(x), " ", fixed=T)[[1]][2]})))

library(rjson)
files <- list.files(paste0(root, "threshold/"))
thr_path <- paste0(root, "threshold/")
thr_df <- data.frame()
for(i in 1:length(files)){
  result <- fromJSON(file = paste0(thr_path, files[[i]]))
  thr_df[i,1] <- result[[2]][2]
  thr_df[i,2] <- result[[2]][5]
}


annot <- read.delim(paste0(root, "column_annotation.txt"), sep="\t", header=T)
unique(annot$V3)[which(unique(annot$V3) %!in% do.call(rbind, lapply(thr_df$measurement, function(x){ strsplit(as.character(x), " ", fixed=T)[[1]][2]})))]
