##############################################################################################################################################
#library
packages <- c("openxlsx","dplyr", "tidyverse", "reshape2" , "ggplot2", "openxlsx", "ComplexHeatmap", "gridExtra", "ggpubr")
lapply(packages, function(x) require(x, character.only=T, quietly=T))

#function_basic
'%!in%' <- function(x,y) !'%in%'(x,y)
theme_sw <- theme_bw()+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
'custom_split' <- function(x, sep, n){ strsplit(as.character(x), sep)[[1]][n]}
'summary_stat' <- function(x){
  mean <- round(mean(x, na.rm=T),1)
  median <- round(median(x, na.rm=T),1)
  QR1 <- round(as.numeric(quantile(x, na.rm=T)[2]),1)
  QR4 <- round(as.numeric(quantile(x, na.rm=T)[4]),1)
  result <- list(mean=mean, median=median, QR1=QR1, QR4=QR4)
  return(unlist(result))
}
'name_celltype' <- function(x){
  x2 <- gregexpr("1", x)[[1]]
  x3 <- paste(abs[x2], collapse=";")
  return(x3)
}

#dir
root <- "/Users/cswc/Documents/spatial_dipg/"

##############################################################################################################################################
#threshold####################################################################################################################################
##############################################################################################################################################
library(rjson)
thres <- paste0(root ,"threshold/") #where all threshold .json files are stored

thres_files <- list.files(thres)
thres_df <- data.frame()
for(i in 1:length(thres_files)){
  
  df <-fromJSON(file=paste0(thres, thres_files[i]))
  measurement <- df[[2]]$measurement
  threshold <- df[[2]]$threshold
  
  thres_df[i, "cell_str"] <- strsplit(measurement, ": ")[[1]][1]
  thres_df[i, "Ab"] <-strsplit(strsplit(measurement, ": ")[[1]][2], " ")[[1]][1]
  thres_df[i, "stat"] <-strsplit(strsplit(measurement, ": ")[[1]][2], " ")[[1]][2]
  thres_df[i, "thres"] <- threshold
}
thres_df$col <- paste0(thres_df$cell_str, "..", thres_df$Ab, ".", thres_df$stat)

#write.table(thres_df, paste0(root, "threshold.txt"), sep="\t", col.names=T, row.names=F, quote=F)

##############################################################################################################################################
###data transformation########################################################################################################################
##############################################################################################################################################
trans <- paste0(root, "data_transformed/") #where are splitted TMA datas are stored
trans_files <- list.files(trans)
#length(trans_files) #17 files for TMA_272

for(i in 1:length(trans_files)){
  df <- read.delim(paste0(trans, trans_files[i]), sep="\t")
  df_2 <- df[,thres_df$col]
  
  df_3 <- df_2
  
  for(m in 1:ncol(df_2)){
    thres <- thres_df[which(thres_df$col == colnames(df_2)[m]), "thres"]
    for(n in 1:nrow(df_2)){
      df_3[n,m] <- ifelse(df_2[n,m] >= thres, 1, 0)
    }
  }
  colnames(df_3) <- thres_df$Ab
  df_4 <- data.frame(name=df$Name, x=df$Centroid.X..µm, y=df$Centroid.Y..µm, df_3)
  
  write.table(df_4, 
              paste0(root, "concat_transformed.txt"), 
              sep="\t", col.names=F, row.names=F, quote=F, append=T)
}


write.table(concat, 
            paste0(root, "concat_transformed.txt"), 
            sep="\t", col.names=T, row.names=F, quote=F)

##############################################################################################################################################
###explore the cell types#####################################################################################################################
##############################################################################################################################################

### below steps are necessary to assign cell identifier code per each row (cell)
#data <- read.table(paste0(root, "concat_transformed.txt"), sep="\t", header=T)

#

#data$tma <- do.call(rbind, lapply(data$name, function(x){ custom_split(x, "-", 1)}))
#data$sub_tma <- do.call(rbind, lapply(data$name, function(x){ custom_split(x, "-", 2)}))
#data$spot <- paste0(data$tma, "_", data$sub_tma)

#assign unique 'code' for each cell type
#for(i in 1:nrow(data)){
#  data[i, "code"]<- as.character(paste0(data[i, 4:19], collapse=""))
#}

#write.table(data, 
#            paste0(root, "concat_transformed_updated.txt"), 
#            sep="\t", col.names=T, row.names=F, quote=F)


###we will use below 'data' with cell identifier code
data <- read.table(paste0(root, "concat_transformed_updated.txt"), sep="\t", header=T)

#number of cells per each tma_spot
data %>% count(spot) %>%
  ggplot(aes(x=n))+geom_histogram()+ 
  ggtitle("total number of cells / TMA") +theme_sw

data %>% count(spot) %>% 
  #summarize(median=median(n), mean=mean(n), QR1=quantile(n)[2], QR3=quantile(n)[4]) %>%
  ggplot(aes(x=n))+geom_histogram()+ 
  ggtitle("total number of cells / TMA") +
  theme_sw

length(unique(data$code)) #2192 distinct cell types based on Ab staining status
data %>% count(code) %>% arrange(desc(n)) #%>%
ggplot(aes(x=n)) + geom_histogram(bins=50)

#unique cell types per each tma spot
df <- data %>% group_by(spot, code) %>% count() %>% dcast(spot ~ code) %>% as.data.frame()
#df >> replace 'na' with '0' (zero)
df2 <- df %>% replace(is.na(.), 0) %>% as.data.frame()

#summary table for tma-spot-wise cell counts
df_2 <- data.frame(spot=df$spot, total = rowSums(df2[,-1]))
for(i in 1:nrow(df)){
  tot <- df_2[i, 2]
  df_2[i, 3] <- length(df[i, 2:ncol(df)][which(is.na(df[i, 2:ncol(df)])== FALSE)])
  df_2[i, 4] <-length(df2[i, 2:ncol(df2)][which(df2[i, 2:ncol(df2)] >= tot * 0.1)])
  df_2[i, 5] <-length(df2[i, 2:ncol(df2)][which(df2[i, 2:ncol(df2)] >= tot * 0.05)])
  df_2[i, 6] <-length(df2[i, 2:ncol(df2)][which(df2[i, 2:ncol(df2)] >= tot * 0.03)])
  df_2[i, 7] <-length(df2[i, 2:ncol(df2)][which(df2[i, 2:ncol(df2)] >= tot * 0.01)])
  df_2[i, 8] <- paste(colnames(df2[i, 2:ncol(df2)][which(df2[i, 2:ncol(df2)] >= tot * 0.1)]), collapse=";")
  df_2[i, 9] <- paste(colnames(df2[i, 2:ncol(df2)][which(df2[i, 2:ncol(df2)] >= tot * 0.05)]), collapse=";")
  df_2[i, 10] <- paste(colnames(df2[i, 2:ncol(df2)][which(df2[i, 2:ncol(df2)] >= tot * 0.03)]), collapse=";")
  df_2[i, 11] <- paste(colnames(df2[i, 2:ncol(df2)][which(df2[i, 2:ncol(df2)] >= tot * 0.01)]), collapse=";")
  
}

colnames(df_2) <- c("spot", "no_total_cell", "no_unique_celltype", 
                    "no_celltype_10percent", "no_celltype_5percent", "no_celltype_3percent","no_celltype_1percent",
                    "celltype_10percent", "celltype_5percent", "celltype_3percent", "celltype_1percent")
#visualization
ggplot(df_2, aes(x=no_unique_celltype))+
  geom_histogram(bins=100)+
  theme_sw

#summary statistics
summary_stat(df_2$no_unique_celltype)
do.call(rbind,lapply(df_2[,c(2:7)], function(x){ summary_stat(x)}))

ggplot(df_2, aes(x=no_celltype_10percent))+geom_histogram(bins=30)

#cell types with pre-defined frequency 
celltype_10 <- unique(unlist(lapply(df_2$celltype_10percent, function(x){ strsplit(x, ";")[[1]]})))
celltype_5 <- unique(unlist(lapply(df_2$celltype_5percent, function(x){ strsplit(x, ";")[[1]]})))
celltype_3 <- unique(unlist(lapply(df_2$celltype_3percent, function(x){ strsplit(x, ";")[[1]]})))
celltype_1 <- unique(unlist(lapply(df_2$celltype_1percent, function(x){ strsplit(x, ";")[[1]]})))

#assign cell-type specific names per each code
abs <- colnames(data)[4:19] #16 antibodies used in TMA_272 datasett

celltypes <- rbind(
  data.frame( celltype=do.call(rbind,lapply(celltype_10, function(x){name_celltype(x) })),
              cat="celltype_10percent"),
  data.frame( celltype=do.call(rbind,lapply(celltype_5, function(x){name_celltype(x) })),
              cat="celltype_5percent"),
  data.frame( celltype=do.call(rbind,lapply(celltype_3, function(x){name_celltype(x) })),
              cat="celltype_3percent"),
  data.frame( celltype=do.call(rbind,lapply(celltype_1, function(x){name_celltype(x) })),
              cat="celltype_1percent")
)



