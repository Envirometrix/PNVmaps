## Mapping potential natural vegetation (PNV)
## Hengl, T., Sanderman, J., Walsh, M.G, Wheeler, I., Harrison, S.P., Colin Prentice, I., 2018. Global Maps of Potential Natural Vegetation: An Assessment of Machine Learning Algorithms for Operational Mapping. PeerJ, in review.
## tom.hengl@gmail.com

list.of.packages <- c("plyr", "parallel", "plotKML", "GSIF", "ranger", "raster", "rgdal", "htmlwidgets", "leaflet", "gbm", "nnet", "glmnet", "doParallel", "dismo", "caret", "devtools", "ggplot2", "Hmisc", "compositions", "factoextra")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE)

setwd("/data/PNV/R_code")
load(".RData")
## RandomForest with quantile regression:
#library(devtools)
#devtools::install_github("imbs-hl/ranger")
#packageDescription("ranger")
library(ranger)
library(GSIF)
library(rgdal)
library(raster)
library(caret)
library(plotKML)
library(scales)
library(parallel)
library(plyr)
source("PNV_mapping_functions.R")
source("/data/LandGIS/R/LandGIS_functions.R")
source("/data/LandGIS/R/saveRDS_functions.R")

## *** 
## PNV Biomes
## *** 

## BIOME 6000 ----
## available via http://www.bridge.bris.ac.uk/resources/Databases/BIOMES_data
## data set prepared by Sandy Sandy P. Harrison and Iain C. Prentice
biome = read.csv("/data/PNV/Data/Biomes/BIOME_6000_EMBSECBIO_DB_classified_plotfile_v2.csv", na.strings=c("","NA"))
str(biome)
## add training points for Brasil:
radam = readOGR("/data/PNV/Data/Biomes/Radam_Vegetacao_SIRGAS.shp")
radam.leg = read.csv("/data/PNV/Data/Biomes/BR_vegetation_BIOME.csv", na.strings=c("","NA"))
br.biome = spsample(radam, n = 550, type = "random")
br.biome.df = over(br.biome, radam)
br.biome.df$New.global.consolidated.biome.scheme = plyr::join(br.biome.df, radam.leg, by="NM_UVEG")$New.global.consolidated.biome.scheme
summary(as.factor(br.biome.df$New.global.consolidated.biome.scheme))
br.biome.df$Site.Name = paste0("SIM_radam_", 1:nrow(br.biome.df))
br.biome.df$Latitude = br.biome@coords[,2]
br.biome.df$Longitude = br.biome@coords[,1]
br.biome.df$Target.age..ka. = 0
br.biome.df$Remove_points = 0
rm(radam)
## Bind all training dfs
biome_b = plyr::rbind.fill(list(biome[biome$Remove_points==0,], br.biome.df[,c("Site.Name","Latitude","Longitude","Target.age..ka.","Remove_points","New.global.consolidated.biome.scheme")]))
## 11040 rows
## filter out duplicates
## location ID:
biome_b$LOC_ID = paste("ID", round(biome_b$Latitude,4), round(biome_b$Longitude,4), sep="_")
summary(duplicated(biome_b$LOC_ID))
#  Mode   FALSE    TRUE 
# logical    8057    2983
fmode <- function(x) names(table(x))[which.max(table(x))]
## take the dominant class:
dup.biome <- ddply(biome_b[!is.na(biome_b$New.global.consolidated.biome.scheme),], .(LOC_ID), .fun = function(x){fmode(x$New.global.consolidated.biome.scheme)}) ##, .parallel=TRUE)
str(dup.biome)
## if there are only 2 classes then take both
ndup.biome <- ddply(biome_b, .(LOC_ID), summarize, N=length(New.global.consolidated.biome.scheme))
str(ndup.biome)
summary(ndup.biome$N)
summary(ndup.biome$N==2)
## 902 points with 2 observations
## second dominant class:
sel.dup2 = biome_b$LOC_ID %in% ndup.biome$LOC_ID[which(ndup.biome$N==2)]
dup2.biome <- data.frame(LOC_ID=biome_b[sel.dup2,"LOC_ID"], V1=biome_b$New.global.consolidated.biome.scheme[sel.dup2])
str(dup2.biome)
dup.biome[dup.biome$LOC_ID=="ID_-0.25_35.33",]; dup2.biome[dup2.biome$LOC_ID=="ID_-0.25_35.33",]
biome_b[biome_b$LOC_ID=="ID_-0.25_35.33",]

## Bind everything together:
biome_f = plyr::join(do.call(rbind, list(dup.biome[dup.biome$LOC_ID %in% ndup.biome$LOC_ID[which(!ndup.biome$N==2)],], dup2.biome)), biome_b[,c("LOC_ID","Site.Name","Longitude","Latitude","MegaBiomes..Scheme.2.")], match="first")
biome_f = biome_f[!is.na(biome_f$V1),]
biome_f$Biome00k_c = make.names(biome_f$V1)
summary(as.factor(biome_f$Biome00k_c)) 
## cool.mixed.forest = 1560
## temperate.deciduous.broadleaf.forest = 989
## warm.temperate.evergreen.and.mixed.forest = 1012
## steppe = 905
## ...
## 20 generalized levels
## 8797 points / 8057 unique locations

plot(biome_f$Longitude, biome_f$Latitude)
coordinates(biome_f) = ~ Longitude + Latitude
proj4string(biome_f) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
unlink("/data/PNV/Data/Biomes/biome00k_points.gpkg")
writeOGR(biome_f, "/data/PNV/Data/Biomes/biome00k_points.gpkg", "biome00k_points", "GPKG")
plotKML(biome_f["Biome00k_c"], folder.name="Biome.6000.Consolidated.Name", file.name="Biome.6000.Consolidated.Name.kml", kmz=TRUE)
saveRDS(biome_f, "/data/PNV/Data/Biomes/biome_f.rds")

## World plot BIOMES ----
library(scales)
library(maps)
library(maptools)
country.m <- map('world', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(country.m$names, ":"), function(x) x[1])
require(maptools)
country <- as(map2SpatialPolygons(country.m, IDs=IDs), "SpatialLines")
png(file="/data/PNV/img/Fig_biome_points_worldmap.png", res=150, width=2000, height=780)
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(country, col="darkgrey", ylim=c(-60, 81))
points(biome_f, pch=21, bg=alpha("blue", 0.6), cex=.6, col="black")
dev.off()
save.image()

## Global covariates 1km ----
unlink(list.files("/data/LandGIS/layers1km", pattern="tif.aux.xml", full.names = TRUE))
unlink(paste0("/data/LandGIS/layers1km/clm_bioclim.var_chelsa.", c(8,9,15,18,19), "_m_1km_s0..0cm_1979..2013_v1.0.tif")) ## missing pixels in Sahara or strange artifacts
covs1km = unlist(sapply(c("merit.dem", "usgs.ecotapestry", "water.table.depth", "inundation.extent", "earthenv.modis", "snow.prob", "monthly.temp", "water.vapor", "precipitation", "admin0", "earthquakes", "esacci.lc.l4", "surface.water.occ", "bioclim.var"), function(x){ c(list.files("/data/LandGIS/layers1km", pattern=x, full.names = TRUE), list.files("/data/LandGIS/upscaled1km", pattern=x, full.names = TRUE)) })) ## "lst_mod11a2"
covs1km = covs1km[!duplicated(covs1km)]
## 157
#s = raster::stack(covs1km)
## Overlay Biome points and covs1km ----
## takes 10 mins
library(snowfall)
snowfall::sfInit(parallel=TRUE, cpus=parallel::detectCores())
snowfall::sfLibrary(rgdal)
snowfall::sfLibrary(raster)
snowfall::sfExport("covs1km", "biome_f")
ov.biome <- snowfall::sfClusterApplyLB(covs1km, function(i){ raster::extract(raster::raster(i), biome_f) })
snowfall::sfStop()
names(ov.biome) = basename(covs1km)
rm.biome0 <- cbind(as.data.frame(biome_f[,c("Site.Name","Biome00k_c")]), as.data.frame(ov.biome))
rm.biome = fix_biome(rm.biome0)
write.csv(rm.biome, "/data/PNV/Data/Biomes/Biome00k_regression_matrix.csv")
zip("/data/PNV/Data/Biomes/Biome00k_regression_matrix.zip", files="/data/PNV/Data/Biomes/Biome00k_regression_matrix.csv")

## Biomes RF model ----
#summary(rm.biome)
biome00k.fm <- as.formula(paste("Biome00k_c ~ ", paste0(names(ov.biome)[-unlist(sapply(c("admin0","esacci.lc.l4","wdpa","bioclimatic.zones"), function(i){grep(i, names(ov.biome))}))], collapse = "+")))
rm.biome.s = rm.biome[complete.cases(rm.biome[,all.vars(biome00k.fm)]),]
## some 300 points have missing values i.e. fall in the sea
m.biome00k <- ranger::ranger(biome00k.fm, rm.biome.s, importance="impurity", probability=TRUE, num.trees=151, mtry=19)
m.biome00k
# Type:                             Probability estimation 
# Number of trees:                  151 
# Sample size:                      8628 
# Number of independent variables:  139 
# Mtry:                             19 
# Target node size:                 10 
# Variable importance mode:         impurity 
# OOB prediction error:             0.3027906
## 70% accuracy
xl1.P <- as.list(ranger::importance(m.biome00k))
print(t(data.frame(xl1.P[order(unlist(xl1.P), decreasing=TRUE)[1:20]])))
#clm_precipitation_imerge.annual_m_1km_s0..0cm_1980..2017_v1.0.tif          138.23305
#clm_monthly.temp_worldclim.chelsa.jan_u.99_1km_s0..0cm_1979..2013_v1.0.tif 124.55348
#clm_precipitation_imerge.may_m_1km_s0..0cm_1980..2017_v1.0.tif             122.79032
#clm_monthly.temp_worldclim.chelsa.feb_u.99_1km_s0..0cm_1979..2013_v1.0.tif 121.91930
# clm_water.vapor_nasa.eo.oct_m_1km_s0..0cm_2000..2017_v1.0.tif              118.16639
# clm_monthly.temp_worldclim.chelsa.dec_u.99_1km_s0..0cm_1979..2013_v1.0.tif 115.98250
# clm_bioclim.var_chelsa.2_m_1km_s0..0cm_1979..2013_v1.0.tif                 113.96553
# clm_bioclim.var_chelsa.3_m_1km_s0..0cm_1979..2013_v1.0.tif                 113.30038
# clm_monthly.temp_worldclim.chelsa.feb_m_1km_s0..0cm_1979..2013_v1.0.tif    107.73351
# clm_bioclim.var_chelsa.12_m_1km_s0..0cm_1979..2013_v1.0.tif                105.62608
# clm_water.vapor_nasa.eo.sep_m_1km_s0..0cm_2000..2017_v1.0.tif              103.60185
# clm_monthly.temp_worldclim.chelsa.mar_u.99_1km_s0..0cm_1979..2013_v1.0.tif  97.70371
# clm_water.vapor_nasa.eo.nov_m_1km_s0..0cm_2000..2017_v1.0.tif               97.03918
# dtm_elevation_merit.dem_m_1km_s0..0cm_2017_v1.0.tif                         88.59517
# clm_monthly.temp_worldclim.chelsa.mar_m_1km_s0..0cm_1979..2013_v1.0.tif     86.63320
# clm_bioclim.var_chelsa.4_m_1km_s0..0cm_1979..2013_v1.0.tif                  83.63234
# clm_monthly.temp_worldclim.chelsa.oct_u.99_1km_s0..0cm_1979..2013_v1.0.tif  82.31959
# clm_bioclim.var_chelsa.16_m_1km_s0..0cm_1979..2013_v1.0.tif                 81.38991
# clm_precipitation_imerge.jun_m_1km_s0..0cm_1980..2017_v1.0.tif              81.28824
# clm_precipitation_imerge.oct_m_1km_s0..0cm_1980..2017_v1.0.tif              78.14946
saveRDS(m.biome00k, "m.biome00k.rds")
save.image()

## Prepare prediction data ----
obj <- GDALinfo("/data/LandGIS/layers1km/lcv_landmask_esacci.lc.l4_c_1km_s0..0cm_2000..2015_v1.0.tif")
tile.lst <- GSIF::getSpatialTiles(obj, block.x=2, return.SpatialPolygons=TRUE)
tile.tbl <- GSIF::getSpatialTiles(obj, block.x=2, return.SpatialPolygons=FALSE)
tile.tbl$ID = as.character(1:nrow(tile.tbl))
str(tile.tbl)
tile.pol = SpatialPolygonsDataFrame(tile.lst, tile.tbl)
writeOGR(tile.pol, "/data/PNV/Data/global/tiles_ll_200km.shp", "tiles_ll_200km", "ESRI Shapefile")
#tile.pol = readOGR("/data/PNV/Data/global/tiles_ll_200km.shp")
system(paste('gdal_translate /data/LandGIS/layers1km/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7_1km_ll.tif /data/tmp/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7_1km_ll.sdat -of \"SAGA\" -ot \"Byte\"'))
system(paste0('saga_cmd -c=24 shapes_grid 2 -GRIDS=\"/data/tmp/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7_1km_ll.sgrd\" -POLYGONS=\"/data/PNV/Data/global/tiles_ll_200km.shp\" -PARALLELIZED=1 -RESULT=\"/data/PNV/Data/global/ov_ADMIN_tiles.shp\"'))
ov_ADMIN = readOGR("/data/PNV/Data/global/ov_ADMIN_tiles.shp", "ov_ADMIN_tiles")
summary(sel.t <- !ov_ADMIN@data[,"ESACCI.LC.L.5"]==210)
ov_ADMIN = ov_ADMIN[sel.t,]
unlink("/data/PNV/Data/global/land_tiles_200km.shp")
writeOGR(ov_ADMIN, "/data/PNV/Data/global/land_tiles_200km.shp", "land_tiles_200km", "ESRI Shapefile")
str(ov_ADMIN@data)
## 5369 dirs
new.dirs <- paste0("/data/PNV/tiled/T", ov_ADMIN$ID)
x <- lapply(new.dirs, dir.create, recursive=TRUE, showWarnings=FALSE)
pr.dirs <- paste0("T", ov_ADMIN$ID)

## Biome legend file ----
col.legend = read.csv("/data/PNV/Data/Biomes/Biome_legend.csv")
col.legend = col.legend[!duplicated(col.legend$New.global.consolidated.biome.scheme),]
col.legend$Group = make.names(col.legend$New.global.consolidated.biome.scheme)
m.biome00k$forest$levels[which(!m.biome00k$forest$levels %in% col.legend$Group)]

## clean-up:
#rds.lst = list.files("/data/PNV/tiled", pattern=".rds", full.names = TRUE, recursive = TRUE)
#x = file.size(rds.lst)
#unlink(rds.lst[x==0])
#x = as.numeric(sapply(sapply(basename(rds.lst), function(i){strsplit(i, "\\.")[[1]][1]}), function(i){strsplit(i, "T")[[1]][2]}))
#unlink(rds.lst[x>12601])
#unlink(rds.lst)
## remove only northern latitudes:
#rds.lst = list.files("/data/PNV/tiled", pattern=".rds", full.names = TRUE, recursive = TRUE)
#rds.lst.t = sapply(basename(rds.lst), function(i){ as.numeric(tools::file_path_sans_ext(strsplit(i, "T")[[1]][2])) })
#summary(rds.lst.t>12268)
#unlink(rds.lst[rds.lst.t>12268])
#unlink(paste0("/data/PNV/tiled/T",1:710,"/T",1:710,".rds"))

## Test it:
#new_data(i="T8547", tif.sel=covs1km[-unlist(sapply(c("esacci.lc.l4","wdpa","bioclimatic.zones"), function(i){grep(i, covs1km)}))], tile.tbl=tile.tbl)
#new_data(i="T12664", tif.sel=covs1km[-unlist(sapply(c("esacci.lc.l4","wdpa","bioclimatic.zones"), function(i){grep(i, covs1km)}))], tile.tbl=tile.tbl)
#new_data(i="T7213", tif.sel=covs1km[-unlist(sapply(c("esacci.lc.l4","wdpa","bioclimatic.zones"), function(i){grep(i, covs1km)}))], tile.tbl=tile.tbl)
#new_data(i="T12378", tif.sel=covs1km[-unlist(sapply(c("esacci.lc.l4","wdpa","bioclimatic.zones"), function(i){grep(i, covs1km)}))], tile.tbl=tile.tbl)
#new_data(i="T596", tif.sel=covs1km[-unlist(sapply(c("esacci.lc.l4","wdpa","bioclimatic.zones"), function(i){grep(i, covs1km)}))], tile.tbl=tile.tbl)

## takes ca 2 hrs...
library(snowfall)
snowfall::sfInit(parallel=TRUE, cpus= parallel::detectCores())
snowfall::sfLibrary(rgdal)
snowfall::sfLibrary(raster)
snowfall::sfExport("new_data", "tile.tbl", "covs1km", "pr.dirs")
out <- snowfall::sfClusterApplyLB(pr.dirs, function(i){ new_data(i, tif.sel=covs1km[-unlist(sapply(c("esacci.lc.l4","wdpa","bioclimatic.zones"), function(i){grep(i, covs1km)}))], tile.tbl=tile.tbl) })
snowfall::sfStop()

## Clean up predictions
#tif.lst = list.files("/data/PNV/tiled", pattern=glob2rx("Biome00k_*_*.tif$"), full.names = TRUE, recursive = TRUE)
#unlink(tif.lst)
## remove northern latitudes
#tif.lst = list.files("/data/PNV/tiled", pattern=glob2rx("Biome00k_C_*.tif$"), full.names = TRUE, recursive = TRUE)
#tif.lst.t = sapply(basename(tif.lst), function(i){ as.numeric(tools::file_path_sans_ext(strsplit(strsplit(i, "_")[[1]][3], "T")[[1]][2])) })
#unlink(tif.lst[tif.lst.t>12268])
#unlink(paste0("/data/PNV/tiled/T",1:710,"/Biome00k_C_T",1:710,".tif"))
#unlink(paste0("/data/PNV/tiled/T",1:710,"/SSI_Biome00k_M_",1:710,".tif"))

## Predict PNV biomes ----
#pred_probs(i="T8547", gm=m.biome00k, tile.tbl=tile.tbl, col.legend=col.legend, varn="Biome00k")
#pred_probs(i="T12664", gm=m.biome00k, tile.tbl=tile.tbl, col.legend=col.legend, varn="Biome00k")
#pred_probs(i="T7213", gm=m.biome00k, tile.tbl=tile.tbl, col.legend=col.legend, varn="Biome00k")

## Takes >30 hours with Jacknifing >128GB RAM:
#cl <- parallel::makeCluster(8, type="FORK")
#x = parallel::clusterApplyLB(cl, pr.dirs, function(i){ try( pred_probs(i, gm=m.biome00k, tile.tbl=tile.tbl, col.legend=col.legend, varn="Biome00k") ) } )
#parallel::stopCluster(cl)
## Results in problems / missing tiles and is very computational
## Takes 30 mins without errors
snowfall::sfInit(parallel=TRUE, cpus= parallel::detectCores())
snowfall::sfLibrary(rgdal)
snowfall::sfLibrary(ranger)
snowfall::sfLibrary(stats)
snowfall::sfLibrary(plyr)
snowfall::sfExport("pred_probs", "tile.tbl", "m.biome00k", "pr.dirs", "col.legend")
out <- snowfall::sfClusterApplyLB(pr.dirs, function(i){ try( pred_probs(i, gm=m.biome00k, tile.tbl=tile.tbl, col.legend=col.legend, varn="Biome00k", with.se=FALSE) ) })
snowfall::sfStop()
## 235 nodes produced errors; first error: Error in data.frame(..., check.names = FALSE) :
## arguments imply differing number of rows: 1, 0

mos.lst = list.files("/data/LandGIS/predicted1km", pattern=glob2rx("pnv_biome.type_biome00k*.tif$"), full.names = TRUE)
unlink(mos.lst)

## make mosaics:
r = raster("/data/GEOG/TAXOUSDA_250m_ll.tif")
te = as.vector(extent(r))[c(1,3,2,4)]
d.lst = m.biome00k$forest$levels
filename = paste0("pnv_biome.type_biome00k.", d.lst, "_p_1km_s0..0cm_2000..2017_v0.1.tif")
library(snowfall)
sfInit(parallel=TRUE, cpus=length(filename))
sfExport("d.lst", "mosaick_ll", "filename", "te")
out <- sfClusterApplyLB(1:length(d.lst), function(x){ try( mosaick_ll(varn="Biome00k_M", i=d.lst[x], out.tif=filename[x], in.path="/data/PNV/tiled", out.path="/data/LandGIS/predicted1km", tr=1/120, te=paste(te, collapse = " "), ot="Byte", dstnodata=255, aggregate=FALSE) )})
sfStop()

mosaick_ll(varn="Biome00k", i="C", out.tif="pnv_biome.type_biome00k_c_1km_s0..0cm_2000..2017_v0.1.tif", dominant = TRUE, in.path="/data/PNV/tiled", out.path="/data/LandGIS/predicted1km", tr=1/120, te=paste(te, collapse = " "), ot="Byte", dstnodata=255, aggregate=FALSE)
save.image()
write.csv(col.legend, "/data/LandGIS/predicted1km/pnv_biome.type_biome00k_c_1km_s0..0cm_2000..2017_v0.1.tif.csv")

## Biome Shannon Entropy Index ----
tmp.lst <- list.files(path="/data/PNV/tiled", pattern=glob2rx(paste0("SSI_Biome00k_*.tif$")), full.names=TRUE, recursive=TRUE)
unlink(tmp.lst)

## takes 30 mins
library(entropy)
library(snowfall)
sfInit(parallel=TRUE, cpus=parallel::detectCores())
sfExport("entropy_tile", "d.lst", "pr.dirs")
sfLibrary(raster)
sfLibrary(rgdal)
sfLibrary(plyr)
sfLibrary(entropy)
out <- sfClusterApplyLB(pr.dirs, function(i){try( entropy_tile(i, in.path="/data/PNV/tiled", varn = "Biome00k_M", levs=d.lst) )})
## 178 nodes produced errors
sfStop()

tmp.lst <- list.files(path="/data/PNV/tiled", pattern=glob2rx(paste0("SSI_Biome00k_*.tif$")), full.names=TRUE, recursive=TRUE)
out.tmp <- tempfile(fileext = ".txt")
vrt.tmp <- tempfile(fileext = ".vrt")
cat(tmp.lst, sep="\n", file=out.tmp)
system(paste0('gdalbuildvrt -input_file_list ', out.tmp, ' ', vrt.tmp))
system(paste0('gdalwarp ', vrt.tmp, ' /data/LandGIS/predicted1km/pnv_biome.type_biome00k_sse_1km_s0..0cm_2000..2017_v0.1.tif -ot \"Byte\" -dstnodata \"255\" -co \"BIGTIFF=YES\" -multi -wm 2000 -co \"COMPRESS=DEFLATE\" -r \"near\"'))

## Biomes model CV ----
library(caret)
library(doParallel)
cl <- makeCluster(62)
registerDoParallel(cl)
tc <- trainControl(method="repeatedcv", number=5, repeats=2, classProbs=TRUE, allowParallel=TRUE, verboseIter=TRUE)
## model training -- takes ca 30 mins!
getModelInfo("ranger")
tg <- expand.grid(mtry=seq(10, 160, by=8), splitrule=c("gini"), min.node.size=10)
#getModelInfo("mxnet")

CP.rf <- train(form=biome00k.fm, data=rm.biome[complete.cases(rm.biome[,all.vars(biome00k.fm)]),], method="ranger", trControl=tc, tuneGrid=tg, na.action=na.omit) ## metric="ROC"
CP.rf
## 0.57
tg2 <- expand.grid(interaction.depth=c(1,2), # Depth of variable interactions
                   n.trees=c(10,20),	        # Num trees to fit
                   shrinkage=c(0.01,0.1),		# Try 2 values for learning rate 
                   n.minobsinnode = 20)
CP.gb <- train(form=biome00k.fm, data=rm.biome[complete.cases(rm.biome[,all.vars(biome00k.fm)]),], method="gbm", preProc=c("center","scale"), tuneGrid=tg2, trControl=tc, na.action=na.omit)
CP.gb
## 0.52
CP.nn <- train(form=biome00k.fm, data=rm.biome[complete.cases(rm.biome[,all.vars(biome00k.fm)]),], method="nnet", preProc=c("center","scale"), trControl=tc, na.action=na.omit)
CP.nn
## 0.43
CP.kn <- train(form=biome00k.fm, data=rm.biome[complete.cases(rm.biome[,all.vars(biome00k.fm)]),], method="kknn", preProc=c("center","scale"), trControl=tc, na.action=na.omit)
CP.kn
## 0.51
stopCluster(cl)
closeAllConnections()
#CP.mx <- train(form=biome00k.fm, data=rm.biome[complete.cases(rm.biome[,all.vars(biome00k.fm)]),], method="mxnet", trControl=tc, na.action=na.omit)
#CP.mx
## Error: Model mxnet is not in caret's built-in library
## Conclusions: ranger achieves best accuracy; second best is "kknn" / "gbm"
biome.results <- resamples(list(ranger=CP.rf, kknn=CP.kn, gbm=CP.gb, nnet=CP.nn))
pdf(file="/data/PNV/img/Fig_boxplot_biomes_accuracy.pdf", width=7, height=5)
bwplot(biome.results, fill="grey")
dev.off()


## Detailed CV for ranger
## drop all levels with <10 observations
keep.b <- levels(as.factor(rm.biome$Biome00k_c))[table(as.factor(rm.biome$Biome00k_c)) > 10]
rm.biome.t = rm.biome[rm.biome$Biome00k_c %in% keep.b,]
rm.biome.t = rm.biome.t[complete.cases(rm.biome.t[,all.vars(biome00k.fm)]),]
cv.biome00k = cv_factor(biome00k.fm, rm.biome.t, nfold=5, idcol="Site.Name", cpus=1)
cv.biome00k[["Cohen.Kappa"]]
cv.biome00k[["Classes"]]
cv.biome00k$Purity
## accuracy = 68%
saveRDS(cv.biome00k, file="cv.biome00k.rds")

## *** 
## FAPAR points ----
## ***

## Simulated desert points
load("/data/profs/TAXOUSDA/deserts.pnt.rda")
deserts.pnt = spTransform(deserts.pnt, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
deserts.pnt = as.data.frame(deserts.pnt)
deserts.pnt$Lon = deserts.pnt$x
deserts.pnt$Lat = deserts.pnt$y

Ni = 50000
mask.lst = c("ifl_2013.tif", "ifl_2000.tif", "WDPA_Dec2017-shapefile-polygons.tif")
## sample points randmoly using masks above (takes >10 mins)
rpnt.lst = parallel::mclapply(mask.lst, function(i){ raster::sampleRandom(raster(paste0("/data/protectedplanet/", i)), size=Ni, sp=TRUE) }, mc.cores=3)
rpnt.df.lst = do.call(rbind, lapply(rpnt.lst, function(i) as.data.frame(i)[,c("x","y")]))
names(rpnt.df.lst) = c("Lon","Lat")
rpnt.xy = do.call(rbind, list(rpnt.df.lst, deserts.pnt[,c("Lon","Lat")]))
## 30301 pnts.
rpnt.xy$ID = paste0(c(make.unique(rep("IFL_2013", nrow(rpnt.lst[[1]]))), make.unique(rep("IFL_2000", nrow(rpnt.lst[[2]]))), make.unique(rep("WDPA_Dec2017", nrow(rpnt.lst[[3]]))), make.unique(rep("Deserts", nrow(deserts.pnt)))))
str(rpnt.xy)
rpnts = rpnt.xy
coordinates(rpnts) = ~ Lon + Lat
proj4string(rpnts) = biome_f@proj4string
## remove points falling in water bodies:
ov.rpnts = raster::extract(raster("/data/LandGIS/layers250m/lcv_landmask_esacci.lc.l4_c_250m_s0..0cm_2000..2015_v1.0.tif"), rpnts)
rpnts = rpnts[!ov.rpnts==2,]
#plot(rpnts)
## Final number of points: 21020
unlink("/data/PNV/Data/FAPAR/intact_rpoints.gpkg")
writeOGR(rpnts, "/data/PNV/Data/FAPAR/intact_rpoints.gpkg", "intact_rpoints", "GPKG")
#rpnts = readOGR("/data/PNV/Data/FAPAR/intact_rpoints.gpkg")

library(rgdal)
GDALinfo("/data/protectedplanet/ifl_2000.tif")
system(paste0('gdalwarp /data/protectedplanet/ifl_2000.tif /data/protectedplanet/ifl_2000_5km.tif -tr 0.05 0.05 -r \"near\" -co \"COMPRESS=DEFLATE\" -overwrite'))
system(paste0('gdalwarp /data/protectedplanet/WDPA_Dec2017-shapefile-polygons.tif /data/protectedplanet/WDPA_Dec2017-shapefile-polygons_5km.tif -tr 0.05 0.05 -r \"near\" -co \"COMPRESS=DEFLATE\" -overwrite'))
g5km = stack(c("/data/protectedplanet/ifl_2000_5km.tif","/data/protectedplanet/WDPA_Dec2017-shapefile-polygons_5km.tif"))
g5km = as(g5km, "SpatialGridDataFrame")
g5km$c = ifelse(!is.na(g5km$ifl_2000_5km), "ifl", ifelse(!is.na(g5km$WDPA_Dec2017.shapefile.polygons_5km), "WDPA", NA))
png(file="/data/PNV/img/Fig_intactareas_worldmap.png", res=150, width=2000, height=780)
par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
image(raster(g5km["c"]), col=c("green","darkgrey"), ylim=c(-60, 81))
lines(country, col="black")
dev.off()

## Extract FAPAR values (takes 10 mins):
fapar.lst = list.files("/data/LandGIS/layers250m", pattern=glob2rx("veg_fapar_proba.v.*_d_250m_s0..0cm_2014..2017_v1.0.tif$"), full.names = TRUE)
fapar.lst = fapar.lst[-grep("annual", fapar.lst)] 
## apr, aug, dec, ...
ov.fapar <- parallel::mclapply(fapar.lst, function(i){ raster::extract(raster(i), rpnts) }, mc.cores = length(fapar.lst))
ov.fapar = as.data.frame(ov.fapar, col.names=basename(fapar.lst))
str(ov.fapar)
#'data.frame':	21020 obs. of  12 variables:
ov.fapar.m = data.frame(FAPAR=as.vector(unlist(ov.fapar[,1:12])), Month=as.vector(sapply(c(4,8,12,2,1,7,6,3,5,11,10,9), function(i){rep(i, nrow(ov.fapar))})), ID=rep(rpnts$ID, 12)) 
str(ov.fapar.m)
## 'data.frame':	252240 obs. of  4 variables
## Cosine function to represent seasonality
ov.fapar.m$cMonth = cos((ov.fapar.m$Month/12)*(2*pi))
## Extract environmental covs from RDS files / tiles (15 minutes):
ov.rds = extract.tiled(x=rpnts, tile.pol=tile.pol, path="/data/PNV/tiled", ID="ID")
str(ov.rds)
## some covariates also refer to months:
# clm_water.vapor_nasa.eo.*_m_1km_s0..0cm_2000..2017_v1.0.tif
# clm_precipitation_imerge.*_m_1km_s0..0cm_1980..2017_v1.0.tif
# clm_monthly.temp_worldclim.chelsa.*_u.99_1km_s0..0cm_1979..2013_v1.0.tif
# clm_monthly.temp_worldclim.chelsa.*_m_1km_s0..0cm_1979..2013_v1.0.tif
# clm_monthly.temp_worldclim.chelsa.*_l.01_1km_s0..0cm_1979..2013_v1.0.tif
# clm_cloud.fraction_earthenv.modis.*_m_1km_s0..0cm_2000..2015_v1.0.tif
# clm_snow.prob_esacci.jan_p_1km_s0..0cm_2000..2016_v1.0.tif
vars.t = c("clm_water.vapor_nasa.eo.","clm_precipitation_imerge.","u.99_1km_s0..0cm_1979..2013", "m_1km_s0..0cm_1979..2013","l.01_1km_s0..0cm_1979..2013","clm_cloud.fraction_earthenv.modis.","clm_snow.prob_esacci.")
ttn.lst = c("water.vapor","precipitation","temp_max","temp_mean","temp_min","could.fraction","snow.prob")
tt.lst = list(NULL)
tt.lstP = list(NULL)
for(p in 1:length(vars.t)){
  if(vars.t[p]=="clm_cloud.fraction_earthenv.modis."){ 
    mnts = c(1,10:12,2:9)
    sel.tt = grep(vars.t[p], names(ov.rds), fixed = TRUE)
    sel.tt = sel.tt[1:12]
  } else {
    mnts = c(4,8,12,2,1,7,6,3,5,11,10,9)
    sel.tt = grep(vars.t[p], names(ov.rds), fixed=TRUE)
    if(vars.t[p]=="m_1km_s0..0cm_1979..2013_v1.0.tif"){
      sel.tt = sel.tt[1:12]
    }
  }
  if(length(sel.tt)==13){ sel.tt = sel.tt[-1] }
  tt.lst[[p]] = data.frame(v1=as.vector(unlist(ov.rds[,sel.tt])), Month=as.vector(sapply(mnts, function(i){rep(i, nrow(ov.rds))})), ID=rep(ov.rds$ID, 12))
  names(tt.lst[[p]])[1] = ttn.lst[p]
  ov.fapar.m[,ttn.lst[p]] = tt.lst[[p]][,1]
  tt.lstP[[p]] = names(ov.rds)[sel.tt]
}
#unlist(tt.lstP)
## bind everything together:
rm.fapar = plyr::join_all(list(ov.fapar.m, ov.rds[,-which(names(ov.rds) %in% unlist(tt.lstP))]))
hist(rm.fapar$FAPAR, breaks=45, col="grey")
#str(rm.fapar)
saveRDS(rm.fapar, "rm.fapar.rds")
#rm.fapar = readRDS("rm.fapar.rds")

## RF model for potential FAPAR ----
#str(rm.fapar)
fm.FAPAR = as.formula(paste("FAPAR ~ cMonth + ", paste0(ttn.lst, collapse="+"), "+", paste0(names(ov.rds)[-which(names(ov.rds) %in% c("FAPAR.m", "X", "Y", "row.index", "ID", "LandCover", "lcv_admin0_fao.gaul_c_1km_s0..0cm_2015_v1.0.tif", unlist(tt.lstP)))], collapse = "+")))
m.FAPAR <- ranger::ranger(fm.FAPAR, rm.fapar[complete.cases(rm.fapar[,all.vars(fm.FAPAR)]),], importance="impurity", mtry = 29, num.trees=151, quantreg = TRUE)
m.FAPAR
# Type:                             Regression 
# Number of trees:                  151 
# Sample size:                      180990 
# Number of independent variables:  63 
# Mtry:                             29 
# Target node size:                 5 
# Variable importance mode:         impurity 
# OOB prediction error (MSE):       244.8929 
# R squared (OOB):                  0.9597724
xl2.P <- as.list(ranger::importance(m.FAPAR))
print(t(data.frame(xl2.P[order(unlist(xl2.P), decreasing=TRUE)[1:20]])))
# clm_precipitation_imerge.annual_m_1km_s0..0cm_1980..2017_v1.0.tif          305822747
# clm_bioclim.var_chelsa.12_m_1km_s0..0cm_1979..2013_v1.0.tif                142821570
# clm_bioclim.var_chelsa.2_m_1km_s0..0cm_1979..2013_v1.0.tif                  94543419
# clm_bioclim.var_chelsa.7_m_1km_s0..0cm_1979..2013_v1.0.tif                  73010931
# clm_bioclim.var_chelsa.4_m_1km_s0..0cm_1979..2013_v1.0.tif                  62760508
# clm_cloud.fraction_earthenv.modis.annual_m_1km_s0..0cm_2000..2015_v1.0.tif  55349672
# precipitation                                                               50751543
# temp_min                                                                    48065402
# temp_mean                                                                   28855837
# clm_bioclim.var_chelsa.3_m_1km_s0..0cm_1979..2013_v1.0.tif                  27168582
# clm_bioclim.var_chelsa.10_m_1km_s0..0cm_1979..2013_v1.0.tif                 18052554
# clm_bioclim.var_chelsa.16_m_1km_s0..0cm_1979..2013_v1.0.tif                 16176921
# clm_bioclim.var_chelsa.5_m_1km_s0..0cm_1979..2013_v1.0.tif                  15409033
# clm_bioclim.var_chelsa.13_m_1km_s0..0cm_1979..2013_v1.0.tif                 14787175
# temp_max                                                                    13096570
# clm_bioclim.var_chelsa.6_m_1km_s0..0cm_1979..2013_v1.0.tif                  12516379
# clm_bioclim.var_chelsa.14_m_1km_s0..0cm_1979..2013_v1.0.tif                 10651143
# water.vapor                                                                  9354118
# clm_bioclim.var_chelsa.1_m_1km_s0..0cm_1979..2013_v1.0.tif                   8742164
# clm_bioclim.var_chelsa.17_m_1km_s0..0cm_1979..2013_v1.0.tif                  8620853
saveRDS.gz(m.FAPAR, "m.FAPAR.rds")
save.image()

## Upper FAPAR ----
faparU.lst = list.files("/data/LandGIS/layers250m", pattern=glob2rx("veg_fapar_proba.v.*_u.975_250m_s0..0cm_2014..2017_v1.0.tif$"), full.names = TRUE)
ovU.fapar <- parallel::mclapply(faparU.lst, function(i){ raster::extract(raster(i), rpnts) }, mc.cores = length(faparU.lst))
ovU.fapar = as.data.frame(ovU.fapar, col.names=basename(faparU.lst))
ovU.fapar.m = data.frame(FAPAR=as.vector(unlist(ovU.fapar[,1:12])), Month=as.vector(sapply(c(4,8,12,2,1,7,6,3,5,11,10,9), function(i){rep(i, nrow(ovU.fapar))})), ID=rep(rpnts$ID, 12)) 
ovU.fapar.m$cMonth = cos((ovU.fapar.m$Month/12)*(2*pi))
for(p in 1:length(vars.t)){
  if(vars.t[p]=="clm_cloud.fraction_earthenv.modis."){ 
    mnts = c(1,10:12,2:9)
    sel.tt = grep(vars.t[p], names(ov.rds), fixed = TRUE)
    sel.tt = sel.tt[1:12]
  } else {
    mnts = c(4,8,12,2,1,7,6,3,5,11,10,9)
    sel.tt = grep(vars.t[p], names(ov.rds), fixed=TRUE)
    if(vars.t[p]=="m_1km_s0..0cm_1979..2013_v1.0.tif"){
      sel.tt = sel.tt[1:12]
    }
  }
  if(length(sel.tt)==13){ sel.tt = sel.tt[-1] }
  tt.lst[[p]] = data.frame(v1=as.vector(unlist(ov.rds[,sel.tt])), Month=as.vector(sapply(mnts, function(i){rep(i, nrow(ov.rds))})), ID=rep(ov.rds$ID, 12))
  names(tt.lst[[p]])[1] = ttn.lst[p]
  ovU.fapar.m[,ttn.lst[p]] = tt.lst[[p]][,1]
  tt.lstP[[p]] = names(ov.rds)[sel.tt]
}
rmU.fapar = plyr::join_all(list(ovU.fapar.m, ov.rds[,-which(names(ov.rds) %in% unlist(tt.lstP))]))
saveRDS.gz(rmU.fapar, "rmU.fapar.rds")
mU.FAPAR <- ranger::ranger(fm.FAPAR, rmU.fapar[complete.cases(rmU.fapar[,all.vars(fm.FAPAR)]),], importance="impurity", mtry=29, num.trees=151, quantreg = TRUE)
mU.FAPAR
# Type:                             Regression 
# Number of trees:                  151 
# Sample size:                      180990 
# Number of independent variables:  63 
# Mtry:                             29 
# Target node size:                 5 
# Variable importance mode:         impurity 
# OOB prediction error (MSE):       299.2672 
# R squared (OOB):                  0.9489668
saveRDS.gz(mU.FAPAR, "mU.FAPAR.rds")

## FAPAR predictions ----

#tif.lst = list.files("/data/PNV/tiled", pattern=glob2rx("FAPAR_*.tif"), full.names = TRUE, recursive = TRUE)
#unlink(tif.lst)
#pred_FAPAR(i="T9638", gm=m.FAPAR, tile.tbl=tile.tbl)
#pred_FAPAR(i="T5454", gm=m.FAPAR, tile.tbl=tile.tbl)
## takes ca 5 hrs
library(snowfall)
snowfall::sfInit(parallel=TRUE, cpus= parallel::detectCores())
snowfall::sfLibrary(rgdal)
snowfall::sfLibrary(ranger)
snowfall::sfExport("pred_FAPAR", "tile.tbl", "m.FAPAR", "pr.dirs")
out <- snowfall::sfClusterApplyLB(pr.dirs, function(i){ try( pred_FAPAR(i, gm=m.FAPAR, tile.tbl=tile.tbl, with.se=TRUE) ) })
#234 nodes produced errors; first error: Error in data.frame(..., check.names = FALSE) : 
#arguments imply differing number of rows: 1, 0
snowfall::sfStop()

## Upper FAPAR
snowfall::sfInit(parallel=TRUE, cpus= parallel::detectCores())
snowfall::sfLibrary(rgdal)
snowfall::sfLibrary(ranger)
snowfall::sfExport("pred_FAPAR", "tile.tbl", "mU.FAPAR", "pr.dirs")
out <- snowfall::sfClusterApplyLB(pr.dirs, function(i){ try( pred_FAPAR(i, gm=mU.FAPAR, level="U", tile.tbl=tile.tbl, with.se=TRUE) ) })
snowfall::sfStop()

m.lst = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
d2.lst = expand.grid(m.lst, c("M", "se"))
filename2 = paste0("pnv_fapar_proba.v.", tolower(d2.lst$Var1), "_", ifelse(d2.lst$Var2=="M","d","sd"),"_1km_s0..0cm_2014..2017_v0.1.tif")
## pnv_fapar_proba.v.sep_d_1km_s0..0cm_2014..2017_v0.1.tif
library(snowfall)
sfInit(parallel=TRUE, cpus=nrow(d2.lst)) # parallel::detectCores())
sfExport("d2.lst", "mosaick_ll", "filename2", "te")
out <- sfClusterApplyLB(1:nrow(d2.lst), function(x){ try( mosaick_ll(varn="FAPAR", i=paste(d2.lst$Var1, d2.lst$Var2, sep="_")[x], out.tif=filename2[x], in.path="/data/PNV/tiled", out.path="/data/LandGIS/predicted1km", tr=1/120, te=paste(te, collapse = " "), ot="Byte", dstnodata=255, aggregate=FALSE) )})
sfStop()
#save.image()

d3.lst = expand.grid(m.lst, c("U", "seU"))
filename3 = paste0("pnv_fapar_proba.v.", tolower(d3.lst$Var1), "_", ifelse(d3.lst$Var2=="U","u.975","sd.u.975"),"_1km_s0..0cm_2014..2017_v0.1.tif")
library(snowfall)
sfInit(parallel=TRUE, cpus= parallel::detectCores())
sfExport("d3.lst", "mosaick_ll", "filename3", "te")
out <- sfClusterApplyLB(1:nrow(d3.lst), function(x){ try( mosaick_ll(varn="FAPAR", i=paste(d3.lst$Var1, d3.lst$Var2, sep="_")[x], out.tif=filename3[x], in.path="/data/PNV/tiled", out.path="/data/LandGIS/predicted1km", tr=1/120, te=paste(te, collapse = " "), ot="Byte", dstnodata=255, aggregate=FALSE) )})
sfStop()
save.image()

## FAPAR mean annual ----
meanf <- function(x){calc(x, mean, na.rm=TRUE)}
sdf <- function(x){calc(x, sd, na.rm=TRUE)*10}
o.l = raster::stack(list.files("/data/LandGIS/predicted1km", glob2rx("pnv_fapar_proba.v.*_d_1km_s0..0cm_2014..2017_v0.1.tif"), full.names = TRUE))
## run in parallel:
beginCluster()
r1 <- clusterR(o.l, fun=meanf, filename="/data/LandGIS/predicted1km/pnv_fapar_proba.v.annual_d_1km_s0..0cm_2014..2017_v0.1.tif", datatype="INT2S", options=c("COMPRESS=DEFLATE"))
## takes ca 15 mins
r2 <- clusterR(o.l, fun=sdf, filename="/data/LandGIS/predicted1km/pnv_fapar_proba.v.annual_sd_1km_s0..0cm_2014..2017_v0.1.tif", datatype="INT2S", options=c("COMPRESS=DEFLATE"))
endCluster()
## vs actual FAPAR:
#o2.l = raster::stack(list.files("/data/LandGIS/upscaled1km", glob2rx("veg_fapar_proba.v.*_d_1km_s0..0cm_2014..2017_v0.1.tif"), full.names = TRUE))
#beginCluster()
#r3 <- clusterR(o2.l, fun=meanf, filename="/data/LandGIS/upscaled1km/veg_fapar_proba.v.annual_d_1km_s0..0cm_2014..2017_v0.1.tif", datatype="INT2S", options=c("COMPRESS=DEFLATE"))
## takes ca 15 mins
#r4 <- clusterR(o2.l, fun=sdf, filename="/data/LandGIS/upscaled1km/veg_fapar_proba.v.annual_sd_1km_s0..0cm_2014..2017_v0.1.tif", datatype="INT2S", options=c("COMPRESS=DEFLATE"))
#endCluster()

## FAPAR model CV ----
## takes >1 hrs
library(caret)
library(doParallel)
tc <- trainControl(method="repeatedcv", number=5, repeats=1, allowParallel=TRUE, verboseIter=TRUE)
tg <- expand.grid(mtry=seq(4, 30, by=5), splitrule=c("gini"), min.node.size=10)
gb.tuneGrid <- expand.grid(eta = c(0.3,0.4,0.5), nrounds = c(50,100,150), max_depth = 2:4, gamma = 0, colsample_bytree = 0.8, min_child_weight = 1, subsample=1)
tg3 = expand.grid(committees=1:2, neighbors=c(2:6))

cl <- makeCluster(62)
registerDoParallel(cl)
FAPAR.rf <- train(fm.FAPAR, rm.fapar[complete.cases(rm.fapar[,all.vars(fm.FAPAR)]),], method="ranger", trControl=tc, na.action=na.omit)
FAPAR.rf
FAPAR.gb <- train(fm.FAPAR, rm.fapar[complete.cases(rm.fapar[,all.vars(fm.FAPAR)]),], method="xgbTree", trControl=tc, tuneGrid=gb.tuneGrid)
FAPAR.gb
FAPAR.mgb <- train(fm.FAPAR, rm.fapar[complete.cases(rm.fapar[,all.vars(fm.FAPAR)]),], method="gbm", preProc=c("center","scale"), tuneGrid=tg2, trControl=tc, na.action=na.omit)
FAPAR.mgb
## very computational...
FAPAR.cub <- train(fm.FAPAR, rm.fapar[complete.cases(rm.fapar[,all.vars(fm.FAPAR)]),], method="cubist", preProc=c("center","scale"), tuneGrid=tg3, trControl=tc, na.action=na.omit)
FAPAR.cub
stopCluster(cl)
closeAllConnections()
save.image()

FAPAR.results <- resamples(list(cubist=FAPAR.cub, gbm=FAPAR.mgb, xgboost=FAPAR.gb, ranger=FAPAR.rf))
pdf(file="/data/PNV/img/Fig_boxplot_FAPAR_accuracy.pdf", width=7, height=5)
bwplot(FAPAR.results, metric=c("RMSE","Rsquared"), fill="grey", scales = list(relation = "free"), xlim = list(c(0, 40), c(0.8, 1.0)))
dev.off()


## LOO cross-validation (not needed in this case but good to check!)
cv.FAPAR = cv_ranger(fm=fm.FAPAR, rmatrix=rm.fapar, idcol="ID")
cv.FAPAR$Summary
## R-square 0.90
## RMSE = 25

range.FAPAR = range(rm.fapar$FAPAR, na.rm=TRUE)
library(hexbin)
library(lattice)
library(scales)
library(gridExtra)
pdf(file = "/data/PNV/img/Fig_correlation_plot_FAPAR_global.pdf", width=6, height=6)
par(oma=c(0,0,0,1), mar=c(0,0,0,2))
plt.FAPAR = hexbinplot(cv.FAPAR$CV_residuals$Observed~cv.FAPAR$CV_residuals$Predicted, colramp=colorRampPalette(plotKML::SAGA_pal[[1]][8:20]), main=paste0("Global prediction model for FAPAR (N = ", prettyNum(length(cv.FAPAR$CV_residuals$Observed), big.mark=","), ")"), ylab="observed", xlab="predicted (machine learning)", type="g", lwd=1, lcex=8, inner=.4, cex.labels=1, xbins=30, asp=1, xlim=c(0,250), ylim=c(0,250), colorcut=c(0,0.005,0.01,0.03,0.07,0.15,0.25,1), panel=pfun)
plt.FAPAR
dev.off()
save.image()

## *** 
## EU Forest species
## *** 

## EU forest points ----
eu.trees = readRDS("/data/PNV/Data/EU_forest/EU_tree_sp_points.rds")
eu.sp = read.csv("/data/PNV/Data/EU_forest/European_forest_species.csv")
names(eu.sp)[which(names(eu.sp)=="GBIF_ID")] = "taxonKey"
str(eu.sp)
eu.sp$sp = make.names(trimws(eu.sp$ScienceName, "right"))
eu.trees$sp = plyr::join(eu.trees, eu.sp, match="first")$sp
summary(as.factor(eu.trees$sp))
## 350787 NA values
eu.trees = eu.trees[!is.na(eu.trees$sp),]
str(eu.trees)
## 1,551,641 obs. of  5 variables
coordinates(eu.trees) = ~ decimalLongitude + decimalLatitude
proj4string(eu.trees) = CRS("+proj=longlat +datum=WGS84")
writeOGR(eu.trees, "EU_tree_sp_points.shp", "EU_tree_sp_points", "ESRI Shapefile")

## overlay takes ca 30 mins
library(snowfall)
snowfall::sfInit(parallel=TRUE, cpus=parallel::detectCores())
snowfall::sfLibrary(rgdal)
snowfall::sfLibrary(raster)
snowfall::sfExport("covs1km", "eu.trees")
ov.eu.trees <- snowfall::sfClusterApplyLB(covs1km, function(i){ raster::extract(raster::raster(i), eu.trees) })
snowfall::sfStop()
names(ov.eu.trees) = basename(covs1km)
rm.eu.trees0 <- as.data.frame(cbind(as.data.frame(eu.trees[,c("sp","w")]), ov.eu.trees))
rm.eu.trees = fix_biome(rm.eu.trees0, Lat="decimalLatitude")
## too large for a csv
#write.csv(rm.eu.trees, "/data/PNV/Data/EU_forest/Species_regression_matrix.csv")
saveRDS.gz(rm.eu.trees, "/data/PNV/Data/EU_forest/Species_regression_matrix.rds")
rm(rm.eu.trees0); gc()
save.image()
#rm.eu.trees = readRDS.gz("/data/PNV/Data/EU_forest/Species_regression_matrix.rds")

## EU forest RF ----
## Assign higher weights to the EU-Forest points
rm.eu.trees$sp.f = rm.eu.trees$sp
rm.eu.trees$sp.f[which(rm.eu.trees$sp.f %in% c("Chamaecyparis.lawsoniana", "Eucalyptus.globulus", "Pseudotsuga.menziesii"))] = NA
sp.fm <- as.formula(paste("sp.f ~ ", paste0(names(rm.eu.trees)[-c(which(names(rm.eu.trees) %in% c("sp","w","sp.f","decimalLongitude","decimalLatitude")), unlist(sapply(c("admin0","esacci.lc.l4","wdpa","bioclimatic.zones"), function(i){grep(i, names(rm.eu.trees))})))], collapse = "+")))
rm.eu.trees.s = rm.eu.trees[complete.cases(rm.eu.trees[,all.vars(sp.fm)]),]
m.eu.trees <- ranger::ranger(sp.fm, rm.eu.trees.s, importance="impurity", probability=TRUE, keep.inbag=TRUE, num.trees=85, mtry=19, case.weights=rm.eu.trees.s$w)
## TAKES >1hr
m.eu.trees
# Type:                             Probability estimation 
# Number of trees:                  85 
# Sample size:                      1532243 
# Number of independent variables:  139 
# Mtry:                             19 
# Target node size:                 10 
# Variable importance mode:         impurity 
# OOB prediction error:             0.772097
saveRDS.gz(m.eu.trees, "m.eu.trees.rds")
xlF.P <- as.list(ranger::importance(m.eu.trees))
print(t(data.frame(xlF.P[order(unlist(xlF.P), decreasing=TRUE)[1:20]])))
# clm_cloud.fraction_earthenv.modis.annual_m_1km_s0..0cm_2000..2015_v1.0.tif 12378.171
# clm_monthly.temp_worldclim.chelsa.feb_u.99_1km_s0..0cm_1979..2013_v1.0.tif 11710.816
# clm_monthly.temp_worldclim.chelsa.oct_u.99_1km_s0..0cm_1979..2013_v1.0.tif 11558.127
# clm_bioclim.var_chelsa.4_m_1km_s0..0cm_1979..2013_v1.0.tif                 10696.869
# dtm_elevation_merit.dem_m_1km_s0..0cm_2017_v1.0.tif                        10223.812
# clm_cloud.fraction_earthenv.modis.7_m_1km_s0..0cm_2000..2015_v1.0.tif       9851.090
# clm_cloud.fraction_earthenv.modis.8_m_1km_s0..0cm_2000..2015_v1.0.tif       8989.338
# clm_bioclim.var_chelsa.3_m_1km_s0..0cm_1979..2013_v1.0.tif                  8858.412
# clm_cloud.fraction_earthenv.modis.11_m_1km_s0..0cm_2000..2015_v1.0.tif      8693.738
# clm_cloud.fraction_earthenv.modis.1_m_1km_s0..0cm_2000..2015_v1.0.tif       8482.923
# clm_cloud.fraction_earthenv.modis.12_m_1km_s0..0cm_2000..2015_v1.0.tif      8481.782
# clm_cloud.fraction_earthenv.modis.6_m_1km_s0..0cm_2000..2015_v1.0.tif       8270.046
# clm_cloud.fraction_earthenv.modis.9_m_1km_s0..0cm_2000..2015_v1.0.tif       8222.094
# clm_cloud.fraction_earthenv.modis.2_m_1km_s0..0cm_2000..2015_v1.0.tif       8142.295
# clm_cloud.fraction_earthenv.modis.3_m_1km_s0..0cm_2000..2015_v1.0.tif       8140.595
# clm_cloud.fraction_earthenv.modis.4_m_1km_s0..0cm_2000..2015_v1.0.tif       8107.808
# dtm_twi_merit.dem_m_1km_s0..0cm_2017_v1.0.tif                               8068.634
# clm_cloud.fraction_earthenv.modis.5_m_1km_s0..0cm_2000..2015_v1.0.tif       7998.804
# clm_cloud.fraction_earthenv.modis.10_m_1km_s0..0cm_2000..2015_v1.0.tif      7981.267
# clm_precipitation_imerge.annual_m_1km_s0..0cm_1980..2017_v1.0.tif           7939.382
save.image()

## EU forest CV ----
library(caret)
tg.rf <- expand.grid(mtry=seq(10, 75, by=8), splitrule=c("gini"), min.node.size=10) 
tc <- trainControl(method="repeatedcv", number=3, repeats=1, classProbs=TRUE, allowParallel=TRUE, verboseIter=TRUE)
tg2 <- expand.grid(interaction.depth=c(1,2), # Depth of variable interactions
                   n.trees=c(10,20),	        # Num trees to fit
                   shrinkage=c(0.01,0.1),		# Try 2 values for learning rate 
                   n.minobsinnode = 20)
library(doParallel)
cl <- makeCluster(parallel::detectCores()) ## , type='PSOCK'
registerDoParallel(cl)
## Too computational hence - subset to ca 2-3% of data
## Too many classes with <5 observations which leads to CV breaking
round((6e4/nrow(rm.eu.trees))*100,2)
rm.eu.t = rm.eu.trees.s[sample.int(nrow(rm.eu.trees.s), size=6e4),]
## drop all levels with <10 observations
keep <- levels(as.factor(rm.eu.t$sp))[table(as.factor(rm.eu.t$sp)) > 10]
rm.eu.t <- rm.eu.t[rm.eu.t$sp %in% keep,]
rm.eu.t$ID = paste("ID", row.names(rm.eu.t), sep="_")
## takes 60 mins
EF.rf <- train(form=sp.fm, data=rm.eu.t, method="ranger", tuneGrid=tg.rf, trControl=tc, na.action=na.omit) ## metric="ROC"
EF.rf
## 0.25
EF.gb <- train(form=sp.fm, data=rm.eu.t, method="gbm", preProc=c("center","scale"), tuneGrid=tg2, trControl=tc, na.action=na.omit)
EF.gb
## 0.21
EF.nn <- train(form=sp.fm, data=rm.eu.t, method="nnet", preProc=c("center","scale"), trControl=tc, na.action=na.omit)
EF.nn
## 0.17
EF.kn <- train(form=sp.fm, data=rm.eu.t, method="kknn", preProc=c("center","scale"), trControl=tc, na.action=na.omit)
EF.kn
## 0.17
stopCluster(cl)

EF.results <- resamples(list(ranger=EF.rf, kknn=EF.kn, gbm=EF.gb, nnet=EF.nn))
pdf(file="/data/PNV/img/Fig_boxplot_EU_forest_accuracy.pdf", width=7, height=5)
bwplot(EF.results, fill="grey")
dev.off()

## Detailed CV for ranger
cv.eu.trees = cv_factor(sp.fm, rm.eu.t, nfold=5, idcol="ID", cpus=1)
cv.eu.trees[["Cohen.Kappa"]]
cv.eu.trees[["Classes"]]
saveRDS(cv.eu.trees, file="cv.eu.trees.rds")
save.image()

## EU forest predict ----

## EU selection of tiles:
system('gdal_translate /data/PNV/Data/EU_forest/Forest_TYPE_MAP_2006_250m.tif /data/PNV/Data/EU_forest/Forest_TYPE_MAP_2006_250m.sdat -of \"SAGA\" -ot \"Byte\"')
tile.pol.eu = spTransform(tile.pol, CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))
writeOGR(tile.pol.eu, "/data/PNV/Data/EU_forest/tiles_eu_200km.shp", "tiles_eu_200km", "ESRI Shapefile")
system(paste0('saga_cmd -c=64 shapes_grid 2 -GRIDS=\"/data/PNV/Data/EU_forest/Forest_TYPE_MAP_2006_250m.sgrd\" -POLYGONS=\"/data/PNV/Data/EU_forest/tiles_eu_200km.shp\" -PARALLELIZED=1 -RESULT=\"/data/PNV/Data/EU_forest/ov_tiles_eu.shp\"'))
ov_tiles_eu = readOGR("/data/PNV/Data/EU_forest/ov_tiles_eu.shp", "ov_tiles_eu")
str(ov_tiles_eu@data)
summary(selS.eu <- (!ov_tiles_eu$Forest_TYPE.5==0)&(!is.na(ov_tiles_eu$Forest_TYPE.5)))
pr.dirs.eu = paste0("T", ov_tiles_eu$ID[which(selS.eu)])
## 995 tiles only
#eu.lst = list.files("/data/PNV/tiled", pattern=glob2rx("EUforest_*.tif$"), full.names = TRUE, recursive = TRUE)
#unlink(eu.lst)

## Takes 6-7hrs mins without errors
m.eu.trees = readRDS.gz("m.eu.trees.rds")
col.legend.for = data.frame(Group=eu.sp$sp, Number=1:nrow(eu.sp))
#x = lapply(pr.dirs.eu, function(i){ try( pred_probs(i, gm=m.eu.trees, tile.tbl=tile.tbl, col.legend=col.legend.for, varn="EUforest", with.se=FALSE) ) } )
cl <- parallel::makeCluster(5, type="FORK")
## max 5 cores since otherwise not enough RAM
x = parallel::parLapply(cl, pr.dirs.eu, function(i){ try( pred_probs(i, gm=m.eu.trees, tile.tbl=tile.tbl, col.legend=col.legend.for, varn="EUforest", with.se=FALSE) ) } )
parallel::stopCluster(cl)

r.eu = raster("/data/PNV/Data/EU_forest/Forest_TYPE_MAP_2006_250m.tif")
# te.eu = as.vector(extent(r.eu))[c(1,3,2,4)]
eu.tiles = ov_tiles_eu[selS.eu,]
eu.tiles.ll = spTransform(eu.tiles, CRS("+proj=longlat +datum=WGS84"))
plot(eu.tiles.ll)
te.eu = as.vector(eu.tiles.ll@bbox)

f.lst = m.eu.trees$forest$levels
filename2 = tolower(paste0("pnv_forest.tree.sp_eu.forest00k.", f.lst, "_p_1km_s0..0cm_2000..2017_v0.1.tif"))
library(snowfall)
sfInit(parallel=TRUE, cpus= parallel::detectCores())
sfExport("f.lst", "mosaick_ll", "filename2", "te.eu")
out <- sfClusterApplyLB(1:length(f.lst), function(x){ try( mosaick_ll(varn="EUforest_M", i=f.lst[x], out.tif=filename2[x], in.path="/data/PNV/tiled", out.path="/data/PNV/predicted1km", tr=1/120, te=paste(te.eu, collapse = " "), ot="Byte", dstnodata=255, aggregate=FALSE) )})
sfStop()

## mask out ADMIN units without predictions
system(paste0('gdalwarp /data/countries/GAUL_ADMIN1_landmask_250m.tif EUforest_ADMIN1.tif -r \"near\" -co \"COMPRESS=DEFLATE\" -tr ', 1/120,' ', 1/120,' -te ', paste(te.eu, collapse = " ")))
eu1km = readGDAL("EUforest_ADMIN1.tif")
eu1km$band1 = as.factor(eu1km$band1)
eu1km$p = readGDAL("/data/PNV/predicted1km/pnv_forest.tree.sp_eu.forest00k.quercus.robur.and.quercus.petraea_p_1km_s0..0cm_2000..2017_v0.1.tif")$band1
str(eu1km@data)
library(plyr)
library(parallel)
doParallel::registerDoParallel(cores = parallel::detectCores())
p.s = ddply(eu1km@data, .(band1), summarize, mean.p=mean(p, na.rm=TRUE), tot.p=sum(is.na(p)), n.p=length(p), .parallel=TRUE)
str(p.s)
#'data.frame':	1353 obs. of  4 variables
sel.admin = paste(p.s$band1[which((p.s$tot.p / p.s$n.p)>0.1)])
str(sel.admin)
extra.p = c("3381","3382","1786","1791","1792","1787","1795","2452","2454","2995","1598","1590","1591","1599","846","842","2345")
sel.admin = c(sel.admin, extra.p)
## some parts of Iceland get filtered out
sel.admin = sel.admin[-which(sel.admin %in% c(1159,1173,1182,1162,1165,1174,1178,1167,1163,1177,1171,1159,2386))]
## 313
eu1km$f = ifelse(!eu1km$band1 %in% sel.admin, eu1km$p, NA)
plot(raster(eu1km["f"]))

GDALinfo("/data/PNV/Data/EU_forest/Forest_TYPE_MAP_2006_250m.tif")
system(paste0('gdalwarp /data/PNV/Data/EU_forest/Forest_TYPE_MAP_2006_250m.tif Forest_TYPE_MAP_2006_1km.tif -t_srs \"', proj4string(eu.tiles.ll), '\" -r \"near\" -co \"COMPRESS=DEFLATE\" -tr ', 1/120,' ', 1/120,' -te ', paste(te.eu, collapse = " ")))
#eu1km = readGDAL("Forest_TYPE_MAP_2006_1km.tif")
#plot(raster(eu1km))
#eu1km$p = readGDAL("/data/LandGIS/predicted1km/EUforest_M_Quercus.ilex_1km_ll.tif")$band1
#rm(eu1km)
gc()
save.image()

#tif.eu.lst = list.files("/data/PNV/predicted1km", pattern=glob2rx("pnv_forest.tree.sp_eu.forest00k.*_p_1km_s0..0cm_2000..2017_v0.1.tif$"), full.names = TRUE)
#snowfall::sfInit(parallel=TRUE, cpus=parallel::detectCores())
#snowfall::sfLibrary(rgdal)
#snowfall::sfExport("mask_admin", "tif.eu.lst", "sel.admin")
#out <- snowfall::sfClusterApplyLB(tif.eu.lst, function(i){ mask_admin(i, sel.admin) })
#snowfall::sfStop()

## EU forest sp correlation ----
eu.rpnt = raster::sampleRandom(raster("/data/PNV/predicted1km/pnv_forest.tree.sp_eu.forest00k.pinus.sylvestris_p_1km_s0..0cm_2000..2017_v0.1.tif"), size=2e4, sp=TRUE)
eu.sp.lst = list.files("/data/PNV/predicted1km", pattern=glob2rx("pnv_forest.tree.sp_eu.forest00k.*_p_1km_s0..0cm_2000..2017_v0.1.tif$"), full.names = TRUE)
eu.sp.rpnt <- parallel::mclapply(eu.sp.lst, function(i){ raster::extract(raster(i), eu.rpnt) }, mc.cores = 16)
eu.sp.rpnt = as.data.frame(eu.sp.rpnt, col.names=sapply(basename(eu.sp.lst), function(i){strsplit(strsplit(i, "_")[[1]][3], "eu.forest00k.")[[1]][2]}))
str(eu.sp.rpnt)
gc()
library("Hmisc")
library(compositions); library(factoextra)
eu.cor <- rcorr(as.matrix(eu.sp.rpnt))
write.csv(eu.cor$r, "/data/PNV/Data/EU_forest/rcorr_matrix_species.csv")
mx = eu.cor$r
mx[mx==1] = NA
mx = apply(mx, 1, which.max)
View(data.frame(x1=names(eu.sp.rpnt), x2=names(eu.sp.rpnt)[mx]))
library(corrplot)
pdf(file="/data/PNV/img/Fig_corrplot_EU_species.pdf", width=10, height=9)
corrplot(stats::cor(acomp(eu.sp.rpnt, total=100)), type="upper", order="hclust", tl.cex=0.6, col=c("black", "white"), bg="lightblue")
dev.off()
#dfc0 = acomp(eu.sp.rpnt, total=100)
#sum(is.na(dfc0))
#sum(is.infinite(dfc0))
#pc.eu <- princomp(dfc0)
#fviz_pca_var(pc.eu, title=NULL, labelsize=5)

## Single point test ----
## http://www.commonland.com/en/projects/187/altiplano-estepario
p1 = data.frame(latitude=37.957332, longitude=-2.163181)
coordinates(p1) = ~ longitude + latitude
proj4string(p1) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
t.lst = c(list.files("/data/LandGIS/predicted1km", pattern=glob2rx("*.tif$"), full.names = TRUE), list.files("/data/LandGIS/layers250m", pattern=glob2rx("veg_fapar_proba.v.*.tif$"), full.names = TRUE), list.files("/data/PNV/predicted1km", pattern=glob2rx("*.tif$"), full.names = TRUE), "/data/LandGIS/layers250m/lcv_land.cover_esacci.lc.l4_c_250m_s0..0cm_2015_v1.0.tif") 
ov.t = parallel::mclapply(t.lst, function(i){ raster::extract(raster(i), p1) }, mc.cores = parallel::detectCores())
ov.t = as.data.frame(ov.t, col.names=basename(t.lst))
write.csv(ov.t, "/data/PNV/R_code/example_altiplano_estepario.csv")
## manual edits:
ov.tf = read.csv("/data/PNV/R_code/example_altiplano_estepario_f.csv")
ov.tf = ov.tf[!is.na(ov.tf$Type_v),]
ov.tf = ov.tf[ov.tf$Type_v %in% c("Median","Median (actual)", "Upper", "Upper (actual)"),]
str(format(Sys.time(), format="%B"))
library(ggplot2)
ov.tf$com = as.character(ov.tf$Type_v) 
ov.tf$com[ov.tf$com =="Upper"] <- "Upper (PNV)" 
ov.tf$com[ov.tf$com =="Median"] <- "Median (PNV)" 
ov.tf$com = as.factor(ov.tf$com) 
ggplot(ov.tf, aes(x=as.POSIXct(paste0("01-", Month, "-20018"), format="%d-%B-%Y"), y=Percent, group=com, color = com, linetype = com)) + geom_line(size=1) + scale_color_manual(values = c("darkgoldenrod2", "darkgoldenrod2", "forestgreen","forestgreen")) + scale_linetype_manual(values = c(1,2,1,2))+ scale_x_datetime(labels = date_format("%b")) + labs(x = "Month", y = "Fraction of FAPAR")+ theme(legend.title=element_blank(), legend.position="bottom", legend.direction="horizontal")
#ov.tf$Value = ifelse(ov.tf$Type_v=="Upper"|ov.tf$Type_v=="Upper (actual)", "Upper", "Median")
#ov.tf$Type = ifelse(ov.tf$Type_v=="Median (actual)"|ov.tf$Type_v=="Upper (actual)", "Actual", "PNV")
#ggplot(ov.tf, aes(x=as.POSIXct(paste0("01-", Month, "-20018"), format="%d-%B-%Y"), y=Percent, group=Type_v)) + geom_line(aes(linetype=Type, col=Value)) +  geom_point() + scale_x_datetime(labels = date_format("%b")) + labs(x = "Month", y = "Fraction of FAPAR")
