---
title: "Spatial_DIPG"
output: html_document
date: "2023-04-17"
---

```{r setup, include=FALSE}
library(spatstat)

vignette('getstart')
library(readxl)
library(ggplot2)
library(tidyverse)
#demo(data)

```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:


# overall process for algorithm
```{r cars}
# read in data
# QC
# annotate cells
# create clusters
  #run PCA

# compare between normal and tumor
```

# read in external files
```{r cars}
# read in data
full_data <- read_excel("/Users/agalianese/Downloads/1_TMA_Analysis_all_Adjusted.xlsx")
threshold <- read_tsv("/Users/agalianese/Downloads/threshold.txt")
antibodies <- read_excel("/Users/agalianese/Downloads/Antibodies.xlsx")

#replace NA w 0
antibodies[is.na(antibodies)] = 0
antibodies
          #Parra, E.R., Ferrufino-Schmidt, M.C., Tamegnon, A. et al. Immuno-profiling and cellular spatial analysis using five immune oncology multiplex immunofluorescence panels for paraffin tumor tissue. Sci Rep 11, 8511 (2021). https://doi.org/10.1038/s41598-021-88156-0
```

# remove NAs from data
```{r cars}
#CK:Cytokeratin
#DAPI:4′,6-diamidino-2-phenylindole

# remove columns full of NA
full_data <- full_data[,colSums(is.na(full_data))<nrow(full_data)]
#remove blank rows, if X and Y are empty, remove the row
full_data <- full_data %>%
  filter(!is.na(`Centroid X µm` & `Centroid Y µm`))
```

# Full code for threshold calculations, commented out for now
```{r cars}
# create gene+ col where 1 is present, and 0 is absent
# all values come from threshold
  # lines not mentioned by antibodies will be commented out for now 
# test <- df %>%
#   #mutate("CD163+" = as.integer(df$`Num CD163` > 350)) %>% #monocyte/macrophage lineage
#   #mutate("CD20+" = as.integer(df$`Num CD20` > 1000)) %>% #B-lymphocyte surface molecule, expressed on all B cells except for plasma cells / plasma blasts
#   mutate("CD3+" = as.integer(df$`Num CD3` > 845.813350340136)) %>% # T lymphocyte
#   #mutate("CD31+" = as.integer(df$`Num CD31` > 300)) %>% #endothelial cells, platelets, macrophages and Kupffer cells, granulocytes, lymphocytes (T cells, B cells, and NK cells), megakaryocytes, and osteoclasts.
#   #mutate("CD34+" = as.integer(df$`Num CD34` > 200)) %>% # stem cells
#   mutate("CD45+" = as.integer(df$`Num CD45` > 364.25019979137)) %>% # immune cells
#   mutate("CD68+" = as.integer(df$`Num CD68` > 265.316990858843)) %>% # monocyte lineage, circulating macrophages, and tissue macrophages
#   mutate("CD8+" = as.integer(df$`Num CD8` > 528.5746)) %>% # cytotoxic t cells
#   #mutate("GFAP+" = as.integer(df$`Num GFAP` > 5000)) %>% #expressed in CNS cells
#   #mutate("H3K27M+" = as.integer(df$`Num H3K27M` > 4760.84571641157)) %>% # histone binding??
#   #mutate("KI67+" = as.integer(df$`Num KI67` > 585.013271684108)) %>% # used to tell when cells are growing, low value indicates resting / quiescent (G0) cells
#   mutate("PD1+" = as.integer(df$`Num PD1` > 420.253436188766)) %>% #T and B cells, involved in apoptosis
#   #mutate("PERK+" = as.integer(df$`Num PERK` > 517.218348312131)) %>% # indicates ER stress?
#   #mutate("VEGFR2+" = as.integer(df$`Num VEGFR2` > 1500)) %>% #
#   #mutate("VIMENTIN+" = as.integer(df$`Num VIMENTIN` > 633.795276722027)) #Mesenchymal cells

```

# threshold calculations, only using antibody combos found in lit
```{r}

summary(full_data$`Num CD3`)
summary(full_data$`Num CD8`)
boxplot(full_data$`Num CD8`)

sum(df$`CD8+`)

df <- full_data %>%
  mutate("CD3+" = as.integer(full_data$`Num CD3` >= 845.813350340136)) %>% # T lymphocyte
  mutate("CD45+" = as.integer(full_data$`Num CD45` >= 364.25019979137)) %>% # immune cells
  mutate("CD68+" = as.integer(full_data$`Num CD68` >= 265.316990858843)) %>% # monocyte lineage, circulating macrophages, and tissue macrophages
  mutate("CD8+" = as.integer(full_data$`Num CD8` >= 528.5746)) %>% # cytotoxic t cells
  mutate("PD1+" = as.integer(full_data$`Num PD1` >= 420.253436188766))  #T and B cells, involved in apoptosis

# sum significant antibody status
df <- df %>%
  mutate(antibody_staining_status = rowSums(df[,31:35]))

# remove cells where nothing significant was stained
df <- df %>%
  filter(antibody_staining_status > 0)

df
#merge based on antibody staining status, then remove those columns as theyre redundant
data <- merge(df, antibodies, by.x = c("CD3+", "CD45+", "CD68+", "CD8+", "PD1+"), by.y = c("CD3+", "CD45+", "CD68+", "CD8+", "PD1+"))
table(data$`Marker co-expression`)
```

# following spatstat vignette('getstart')
```{r}

#range of x and y coords
summary(df$`Centroid X µm`)
summary(df$`Centroid Y µm`)


# general form is ppp(x.coord, y.coord, x.range, y.range)
  # i rounded y coords min and max
mypattern <- ppp(df$`Centroid X µm`, df$`Centroid Y µm`, c(2577, 28879), c(525, 17782))
plot(mypattern)
summary(mypattern)

# Ripley's K function
plot(Kest(mypattern))

#envelopes of K function
plot(envelope(mypattern,Kest))

#kernel smoother of density
plot(density(mypattern))

# additional data can be added as marks
marks(mypattern) <- df[,c(4,32)]

mypattern

plot(Smooth(mypattern))



```


```{r cars}


results <- prcomp(full_data[,14:30], scale=T)

test <- 
df
# QC - remove x% of 'low-quality' data - whats low quality, how to remove


head(full_data[,14:30])

df


# annotate cells
dput(head(returns, 10))

# create clusters of neighborhoods



# compare between normal and tumor using berman
# Berman's Test - compare the observed distribution of the values of a spatial covariate at the data points, and the predicted distribution of the same covariate under the model.
# https://rdrr.io/cran/spatstat.model/man/berman.test.html#:~:text=The%20test%20is%20performed%20by%20comparing%20the%20observed,distribution%20of%20the%20same%20covariate%20under%20the%20model.



```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
while DIPGs invariably exhibit tumour spread outside of the pons and brainstem 13 , there are to date no comprehensive molecular studies of these invading cells.
```
