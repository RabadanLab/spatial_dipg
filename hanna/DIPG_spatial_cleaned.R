####### SPATIAL DIPG DATA PROCESSING - LAST UPDATED: 01-11-2023 #######################################

packages <- c("openxlsx","dplyr", "tidyverse", "ggplot2", "ComplexHeatmap", "tidyr")
lapply(packages, function(x) require(x, character.only=T, quietly=T))

#Read in TMA files
TMA_273_Dec2022 <- read.xlsx(paste0(path, "TMA273_CellbyCell_12_2022.xlsx"), sheet=1)
TMA_272_Dec2022 <- read.xlsx(paste0(path, "TMA272_CellbyCell_12_2022.xlsx"), sheet=1)

######## BINARIZATION ##################################################################################

#First, split up each TMA (272 and 273) into A-Q datasets. For example, "A_272"
#Second, I filtered these subsetted TMA files to include only the antibody-read columns that were used in assigning the thresholds, see below. I called these for example, "A_272_filter"

Threshold_list = list('Nucleus:.CD3.mean' = 845.8133503,
                      'Nucleus:.CD4.mean' = 200,
                      'Nucleus:.CD8.mean' = 528.5746,
                      'Nucleus:.CD20.mean' = 1000,
                      'Cell:.CD31.mean' = 300,
                      'Nucleus:.CD34.mean' = 1867.849396,
                      'Nucleus:.CD45.mean' = 364.2501998, 
                      'Nucleus:.CD68.mean' = 265.3169909,
                      'Nucleus:.CD163.mean' = 350,
                      'Cell:.GFAP.mean' = 5000,
                      'Nucleus:.H3K27M.mean' = 4760.845716,
                      'Nucleus:.KI67.mean' = 585.0132717,
                      'Nucleus:.PD1.mean' = 420.2534362,
                      'Nucleus:.PERK.mean' = 517.2183483,
                      'Cell:.VEGFR2.mean' = 1500,
                      'Nucleus:.VIMENTIN.mean' = 633.79527672,
                      'Nucleus:.DAPI.S001.mean' = 0)

df_list_272filt = list(A_272_filter, B_272_filter, C_272_filter, D_272_filter, E_272_filter, F_272_filter,
                   G_272_filter, H_272_filter, I_272_filter, J_272_filter, K_272_filter, L_272_filter, M_272_filter, 
                   N_272_filter, O_272_filter, P_272_filter, Q_272_filter)
df_list_273filt = list(A_273_filter, B_273_filter, C_273_filter, D_273_filter, E_273_filter, F_273_filter,
                   G_273_filter, H_273_filter, I_273_filter, J_273_filter, K_273_filter, L_273_filter, M_273_filter, 
                   N_273_filter, O_273_filter, P_273_filter, Q_273_filter)

#Assign binary numerals for 272 dfs, do the same for 273 list :: 1 if >= threshold, 0 if < threshold
for(m in 1:length(df_list_272filt)){
  df = as.data.frame(df_list_272filt[m])
  binary = df[1:14]
  for(i in 1:17){
    for(j in 1:nrow(df)){
      if(df[j, i+14] >= Threshold_list[[i]]){
        binary[j,i+14] <- 1
      }else{
        binary[j,i+14] <- 0
      }
    }
  }
  colnames(binary) = c('Image', 'Name', 'Class', 'TMA.core', 'Parent', 'ROI', 'X', 'Y', 'Nucleus_Area', 'Nucleus_Perimeter', 'Nucleus_Circularity', 
                                'Nucleus_Max_caliper', 'Nucleus_Min_caliper', 'Nucleus_Eccentricity', 'CD3', 'CD4', 'CD8','CD20','CD31',
                                'CD34','CD45', 'CD68','CD163','GFAP','H3K27M','KI67','PD1','PERK',
                                'VEGFR2','VIMENTIN','DAPI')
#  write.csv(binary, paste0(path6, alpha[m], "_272_binary.csv"), quote = F)
  assign(paste0(alpha[m], "_272_binary"), binary)
}

#Merge binary files by TMA and then 272+273 combined
df_list_272binary = list(A_272_binary, B_272_binary, C_272_binary, D_272_binary, E_272_binary, F_272_binary,
G_272_binary, H_272_binary, I_272_binary, J_272_binary, K_272_binary, L_272_binary, M_272_binary, 
N_272_binary, O_272_binary, P_272_binary, Q_272_binary)

BINARY_272_MERGED = bind_rows(df_list_272binary, .id = NULL)
nrow(BINARY_272_MERGED) == nrow(TMA_272_Dec2022)

# CREATE UNIQUE ID
BINARY_272_MERGED$UniqueID = paste0("272_",BINARY_272_MERGED$TMA.core)
#write.csv(BINARY_272_MERGED, paste0(path8,"MERGED_272_binary.csv"), quote = F)

#Merge both 272 and 273 into complete binary file
DIPG_MERGED = bind_rows(BINARY_272_MERGED, BINARY_273_MERGED, .id = NULL)

######################## CONCATINATION AND SUMMATION OF ANTIBODY POSITIVITY USING FILTERED DIPG DATASET #####################################################

#Read in filtered metadata for:
  # DIPG diagnosis only
  # With confirmed either normal/tumor tissue classification (based on path columns - see metadata sheet)
  # Confirmed all columns of TX info
  # Result = 37 patients
metadata_DIPG <- read.xlsx(paste0(path, "DIPG_metadata_12-20-2022.xlsx"), sheet=1)
metadata_DIPG_TUMOR <- metadata_DIPG[which(metadata_DIPG$TUMOR == 1), ]
metadata_DIPG_NORMAL <- metadata_DIPG[which(metadata_DIPG$NORMAL == 1), ]

DIPG_MERGED_TUMOR <- DIPG_MERGED[which(DIPG_MERGED$UniqueID %in% metadata_DIPG_TUMOR$UniqueID),]
# 207 spots
DIPG_MERGED_NORMAL <- DIPG_MERGED[which(DIPG_MERGED$UniqueID %in% metadata_DIPG_NORMAL$UniqueID),]
# 59 spots

#Exclude antibodies that aren't ready for analyzing. Keep: CD3, CD8, CD31, CD34, CD163, H3K27M, KI67, VIMENTIN
DIPG_MERGED_TUMOR = DIPG_MERGED_TUMOR[, -c(16,18,21,22,24,27,28,29)]
DIPG_MERGED_NORMAL = DIPG_MERGED_NORMAL[, -c(16,18,21,22,24,27,28,29)]

## Writing for processing just the tumor samples below, but can do with normal as well 

#Concatenate cell codes
DIPG_MERGED_TUMOR <- DIPG_MERGED_TUMOR %>% tidyr::unite(cellcode, CD3:VIMENTIN, sep = "_", remove = F)
length(unique(DIPG_MERGED_TUMOR$cellcode))
#With exclusion of Abs = 138

#Assess the number of each unique cell code
df_tumor <- DIPG_MERGED_TUMOR %>% 
  count(cellcode)
df_tumor <- df_tumor[order(df_tumor$n, decreasing =  TRUE),]

#In excel, I deconvoluted the cell code into antibody stain based on column location
cellcode_concat_tumor <- read.xlsx(paste0("DIPG_TUMOR_CELLCODE_CONCAT_01-04-2023.xlsx"), sheet=1)
#e.g. 
# cellcode	        CONCAT
# 0_0_0_0_0_0_0_0	NEGATIVE
# 0_0_0_0_0_1_0_0	H3K27M+

#Integrate into df
DIPG_MERGED_TUMOR_ANNOT = merge(DIPG_MERGED_TUMOR, cellcode_concat_tumor, by = "cellcode")
#confirm merge
length(grep("NEGATIVE", DIPG_MERGED_TUMOR_ANNOT$CONCAT))
#114308

#Sum antibody combinations for each sample (Unique ID)
#e.g. spot B-25 for TMA 272 has 18 cells that are "VIMENTIN+"

for(m in 1:length(unique(DIPG_MERGED_TUMOR_ANNOT$UniqueID))){
  df = DIPG_MERGED_TUMOR_ANNOT[which(DIPG_MERGED_TUMOR_ANNOT$UniqueID == unique(DIPG_MERGED_TUMOR_ANNOT$UniqueID)[m]),]
  sumdf = data.frame(sum = 0)
  n = length(which(df$CONCAT == unique(DIPG_MERGED_TUMOR_ANNOT$CONCAT)[1]))
  sumdf[1,] = n
  for(i in 2:length(unique(DIPG_MERGED_TUMOR_ANNOT$CONCAT))){
    t = length(which(df$CONCAT == unique(DIPG_MERGED_TUMOR_ANNOT$CONCAT)[i]))
    sumdf[nrow(sumdf) + 1, ] = t
  }
  row.names(sumdf) = unique(DIPG_MERGED_TUMOR_ANNOT$CONCAT)
  colnames(sumdf) = paste0(unique(DIPG_MERGED_TUMOR_ANNOT$UniqueID)[m])
  assign(paste0("sumdf_", m), sumdf)
#  write.csv(sumdf, paste0(path, unique(DIPG_MERGED_TUMOR_ANNOT$UniqueID)[m], "_df.csv"), quote = F)
}

#Create master df will all the antibody sums for each Unique ID
masterdf_tumor = data.frame(default = 1:61)
row.names(masterdf_tumor) = unique(DIPG_MERGED_TUMOR_ANNOT$CONCAT)

#Tranpose and add PT_LOCATION identifer = (PID_Location) from metadata
#E.g. df ::
#UniqueID	PT_LOCATION	    TOTAL	CD163+	CD163+H3K27M+, etc....
#272_B-9	1234_Frontal 	631	    17	    0            , etc....

MERGED_trans_tumor <- read.xlsx(paste0(path, "mergedsums_tumor_transposed_01-05-2023.xlsx"), sheet=1)

#Create averages for replicates across each PT_LOCATION
#61 tumor, 21 normal
ptdf = data.frame(default = 1:61)
row.names(ptdf) = unique(MERGED_trans_tumor$PT_LOCATION)
for(k in 1:nrow(ptdf)){
  for(l in 3:ncol(MERGED_trans_tumor)){
    t = mean(MERGED_trans_tumor[which(MERGED_trans_tumor$PT_LOCATION == row.names(ptdf)[k]), l])
    ptdf[k, l-1] = t
  }
}
ptdf = ptdf[, -1]
colnames(ptdf) = colnames(MERGED_trans_tumor[3:64])


