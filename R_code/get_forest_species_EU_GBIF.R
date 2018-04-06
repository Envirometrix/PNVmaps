## Download all forest/tree species using GBIF DB
## tom.hengl@envirometrix.net

setwd("/mnt/nas/GBIF/forest_species")
library(rgdal)
library(rgbif)
library(plyr)
library(data.table)

## Targeted EU forest species:
eu.sp = read.csv("European_forest_species.csv")

## Download species occurrences using GBIF
getOption("gbif_user")
occ_count(taxonKey=eu.sp$GBIF_ID[1], "georeferenced = TRUE")
res <- occ_download(paste0("taxonKey = ", paste(eu.sp$GBIF_ID, collapse=",")), "hasCoordinate = TRUE")
occ_download_get(res, overwrite=TRUE)
## Occurrences 8,559,574
## Citation GBIF Occurrence Download doi:10.15468/dl.fhucwx accessed via GBIF.org on 24 Jan 2018
occ_download_cancel_staged()

x = fread("/mnt/nas/GBIF/forest_species/0012459-171219132708484/occurrence.txt", nrows=0)
names(x)
## TAKES >20 mins
occ.sp = fread("/mnt/nas/GBIF/forest_species/0012459-171219132708484/occurrence.txt", fill=TRUE, blank.lines.skip=TRUE, select=c("year","countryCode","decimalLongitude","decimalLatitude","taxonKey","hasGeospatialIssues","coordinateUncertaintyInMeters"))
str(occ.sp)
## 3908804 obs. of  7 variables
summary(as.factor(occ.sp$taxonKey))
summary(as.factor(occ.sp$hasGeospatialIssues))
summary(occ.sp$coordinateUncertaintyInMeters)
occ.xy = occ.sp[which(occ.sp$coordinateUncertaintyInMeters<2000 & occ.sp$hasGeospatialIssues=="false"),c("taxonKey","decimalLongitude","decimalLatitude","year")]
str(occ.xy)
occ.xy$w = 1
## 1417171 obs. of  4 variables

## doi:10.1038/sdata.2016.123
## NFI for National Forest Inventory dataset, FF for Forest Focus and BS for Biosoil
euf = fread("/mnt/nas/EU_forest/occurrences/EUForestspecies.csv")
coordinates(euf) = ~ X+Y
proj4string(euf) = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")
euf.ll = as.data.frame(spTransform(euf[c("SPECIES NAME")], CRS("+proj=longlat +datum=WGS84")))
names(euf.ll) = c("NameS","decimalLongitude","decimalLatitude")
## 588983 obs. of  3 variables
summary(as.factor(euf.ll$NameS))
write.csv(data.frame(NamesS=levels(as.factor(euf.ll$NameS))), "/mnt/nas/EU_forest/occurrences/EUForestspecies_names_list.csv")
euf.ll$taxonKey = plyr::join(euf.ll, eu.sp, match="first")$GBIF_ID
euf.ll$w = 4
str(euf.ll)

## Merge two files:
eu.trees = rbind.fill(occ.xy, euf.ll[,c("taxonKey","decimalLongitude","decimalLatitude","w")])
eu.trees = eu.trees[!is.na(eu.trees$taxonKey),]
str(eu.trees)
## 1,902,428 obs. of  4 variables
## Export for spatial analysis:
saveRDS(eu.trees, "/data/PNV/Data/EU_forest/EU_tree_sp_points.rds")
gc()
save.image()

library(maps)
library(maptools)
library(scales)
country.m <- map('world', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(country.m$names, ":"), function(x) x[1])
require(maptools)
country <- as(map2SpatialPolygons(country.m, IDs=IDs), "SpatialLines")
png(file="/data/PNV/img/Fig_EU_forest_species_worldmap.png", res=150, width=2000, height=780)
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(country, col="darkgrey", ylim=c(-60, 81))
points(x=eu.trees$decimalLongitude, y=eu.trees$decimalLatitude, pch=21, bg=alpha("yellow", 0.5), cex=.6, col="black")
dev.off()
