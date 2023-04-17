##revising the metadata## \
\
[1] column 8,9,10
````sh
for(i in 9:11){ annot[,i] <- ifelse(annot[,i] %!in% c("N", "Y"), NA, annot[,i])}
````
````sh
colnames(annot)[c(2, 6:8)] <- c("spot", "path_TN", "path_H3K27M", "path_H7E")
````
\
[2] column 4 - histology
````sh
annot$path <- ifelse(annot[,4] %in% c("DMG spinal", "Thalamic DMG", "DMG Medulla"), "DMG",
                     ifelse(annot[,4] %in% c("Ependymoma","Ependymoma " ), "ependymoma",
                            ifelse(annot[,4] %in%  ("Anaplastic Astrocytoma" ), "AA", 
                                   ifelse(annot[,4] %in% c("HG" ,"HGG" ), "HGG", 
                                          ifelse(annot[,4] %in% c("DIPG (Use this data intarnally only!)"), "DIPG",annot[,4])))))
````

[3] column 5 - location
````sh
annot$location <- gsub(" ", "", tolower(annot$Location)) 
````

[4] column 8 - path_H7E
````sh
annot$path_H7E <- ifelse(annot$path_H7E %in% c("CB", "?", "B") , NA, ifelse(annot$path_H7E %in% c("M"), "N+", annot$path_H7E))
````
**revised metadata is uploaded as 'metadata_Apr-2023.txt'**

##histologyg classification -Tumor vs. Normal##


