## Mapping potential natural vegetation mapping at 250 m
## tom.hengl@envirometrix.net
## Data sources: Geo-wiki points (https://zenodo.org/record/3356758), LandPKS points (https://landpotential.org/data-portal/) and simulated points

library(rgdal)
library(raster)
library(plyr)
library(fastSave)
library(googledrive)
library(openxlsx)
library(landmap)
library(mlr)
load(".RData")
source("/data/LandGIS/R/saveRDS_functions.R")
source("/data/LandGIS/R/LandGIS_functions.R")
lc100m = raster("/mnt/DATA/Copernicus_vito/GLC100m/discrete-classification.vrt")

## GLC training points ----
gs.url = as_id("https://docs.google.com/spreadsheets/d/10Y8PryG5TiPZTiArteRXE0iFdKUhNYhj38qb9bRnpVU/edit")
xlsxFile = "Global_forest_management_and_LC_mapping_tables.xlsx"
drive_download(gs.url, xlsxFile, overwrite = TRUE)
esa.tbl = openxlsx::read.xlsx(xlsxFile, sheet = 1)
head(esa.tbl)

## tiling system
tile.tbl = readRDS("/data/LandGIS/models/stacked250m_tiles.rds")
pr.dirs = basename(dirname(list.files(path=path, pattern=glob2rx("*.rds$"), recursive=TRUE)))
tile.pol = readOGR("/data/LandGIS/models/tiles_ll_100km.shp")
tile.pol = tile.pol[paste0("T", tile.pol$ID) %in% pr.dirs,]

## GLC overlay ---- 
tot_pnts = readRDS("/data/LandGIS/training_points/vegetation/glc100m_natural_vegetation.pnts.rds")
## (takes 30 mins):
ov.lc <- extract.tiled(obj=tot_pnts, tile.pol=tile.pol, path="/data/tt/LandGIS/grid250m", ID="ID", cpus=64)
ov.lc$map_code_pnv.f = as.factor(make.names(ov.lc$map_code_pnv))
summary(ov.lc$map_code_pnv.f)
#head(ov.lc)
pr.vars = make.names(unlist(sapply(c("sm2rain","monthly.temp_worldclim.chelsa","bioclim.var_chelsa","irradiation_solar.atlas", "usgs.ecotapestry", "floodmap.500y", "water.table.depth_deltares", "snow.prob_esacci", "water.vapor_nasa.eo", "wind.speed_terraclimate", "merit.dem_m", "merit.hydro_m", "cloud.fraction_earthenv", "water.occurance_jrc", "wetlands.cw_upmc"), function(i){names(ov.lc)[grep(i, names(ov.lc))]})))
str(pr.vars)
## 203
saveRDS.gz(ov.lc, "./training_points/ov_glc100m_natural_vegetation.pnts.rds")
#ov.lc = readRDS.gz("./training_points/ov_glc100m_natural_vegetation.pnts.rds")

## Subset data ----
rm.lc = ov.lc[complete.cases(ov.lc[,pr.vars]),]
dim(rm.lc)
## 51675   376
formulaString.GLC = as.formula(paste('map_code_pnv.f ~ ', paste(pr.vars, collapse="+")))
saveRDS.gz(rm.lc, "./training_points/regression.matrix_tot_pnts.rds")
esa.tbl$response = make.names(esa.tbl$map_code)
esa.tbl$Number = esa.tbl$map_code
save.image.pigz(n.cores = 64)

## GLC testing RF ----
## Geographically distributed sample:
pnts.s <- GSIF::sample.grid(tot_pnts, cell.size=c(1,1), n=2)
length(pnts.s$subset)
## 13881
## determine Mtry / optimal subset of covs:
df = rm.lc[rm.lc$point_id %in% pnts.s$subset$point_id,]
## Remove smaller classes as they leads to errors in the train function
xg = summary(as.factor(df$map_code_pnv.f))
selg.levs = attr(xg, "names")[xg > 5]
df$map_code_pnv[which(!df$map_code_pnv.f %in% selg.levs)] <- NA
df$map_code_pnv <- droplevels(df$map_code_pnv.f)
df <- df[complete.cases(df[,all.vars(formulaString.GLC)]),all.vars(formulaString.GLC)]
## test run:
library(caret)
t.mrfX0 <- caret::train(formulaString.GLC, data=df, method="ranger", 
                   trControl = trainControl(method="repeatedcv", classProbs=TRUE, number=3, repeats=1),
                   na.action = na.omit, num.trees=85, importance="impurity",
                   tuneGrid=expand.grid(mtry = c(15,40,80,110,130,150), splitrule="gini", min.node.size=10))
t.mrfX0
## Very small difference - use default value
# mtry  Accuracy   Kappa    
# 15   0.7389448  0.7012920
# 40   0.7409291  0.7036902

## GLC ranger model ----
summary(as.factor(rm.lc$source_db))
case.weights = ifelse(rm.lc$source_db=="MangrovesDB" | rm.lc$source_db=="Biome6000" | rm.lc$source_db=="LandPKS_land_cover", 10, ifelse(rm.lc$source_db=="Geo-wiki_HIF" | rm.lc$source_db=="GeoFeedback", 3, 1))
summary(as.factor(case.weights))
mrfX_lc <- ranger::ranger(formulaString.GLC, rm.lc, importance="impurity", mtry=t.mrfX0$bestTune$mtry, probability=TRUE, num.trees=85, case.weights=case.weights) 
mrfX_lc
# Type:                             Probability estimation 
# Number of trees:                  85 
# Sample size:                      51675 
# Number of independent variables:  203 
# Mtry:                             40 
# Target node size:                 10 
# Variable importance mode:         impurity 
# Splitrule:                        gini 
# OOB prediction error (Brier s.):  0.2680825
saveRDS.gz(mrfX_lc, "./models/mrfX_lc.rds")
#mrfX_lc = readRDS.gz("./models/mrfX_lc.rds")
xl <- as.list(ranger::importance(mrfX_lc))
print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:25]])))
# dtm_elevation_merit.dem_m_250m_s0..0cm_2017_v1.0.tif                       1238.8743
# clm_water.vapor_nasa.eo.oct_m_1km_s0..0cm_2000..2017_v1.0.tif               966.2789
# clm_precipitation_sm2rain.annual_m_1km_s0..0cm_2007..2018_v0.2.tif          958.0307
# clm_water.vapor_nasa.eo.nov_m_1km_s0..0cm_2000..2017_v1.0.tif               885.0004
# clm_water.vapor_nasa.eo.may_m_1km_s0..0cm_2000..2017_v1.0.tif               807.6233
# clm_bioclim.var_chelsa.7_m_1km_s0..0cm_1979..2013_v1.0.tif                  778.6673
# clm_bioclim.var_chelsa.2_m_1km_s0..0cm_1979..2013_v1.0.tif                  673.8871
# clm_water.vapor_nasa.eo.apr_m_1km_s0..0cm_2000..2017_v1.0.tif               644.2846
# clm_direct.irradiation_solar.atlas.kwhm2.10_m_1km_s0..0cm_2016_v1.tif       612.3663
# clm_precipitation_sm2rain.sep_m_1km_s0..0cm_2007..2018_v0.2.tif             585.1126
# clm_bioclim.var_chelsa.4_m_1km_s0..0cm_1979..2013_v1.0.tif                  524.8560
# clm_bioclim.var_chelsa.12_m_1km_s0..0cm_1979..2013_v1.0.tif                 513.6335
# clm_water.vapor_nasa.eo.dec_m_1km_s0..0cm_2000..2017_v1.0.tif               510.7857

library(mlr)
tsk.C <- mlr::makeClassifTask(data = rm.lc[,all.vars(formulaString.GLC)], target = all.vars(formulaString.GLC)[1])
discrete_ps = makeParamSet( makeDiscreteParam("mtry", values = seq(5,80,by=5)) )
ctrl = makeTuneControlGrid()
rdesc = makeResampleDesc("CV", iters = 3L)
parallelMap::parallelStartSocket(parallel::detectCores())
#resC = tuneParams("classif.ranger", task = tsk.C, resampling = rdesc, par.set = discrete_ps, control = ctrl)
resC = tuneParams(mlr::makeLearner("classif.ranger", num.threads = round(parallel::detectCores()/16), num.trees=85), task = tsk.C, resampling = rdesc, par.set = discrete_ps, control = ctrl)
parallelMap::parallelStop()
resC
#p. pars: mtry=50
## GLC Feature selection ----
## TAKES >15mins
lrn.rf = mlr::makeLearner("classif.ranger", num.threads = parallel::detectCores(), mtry=resC$x$mtry, num.trees=85, predict.type = "prob")
lrn.rf
outer = makeResampleDesc("CV", iters = 3L)
inner = makeResampleDesc("Holdout")
ctrl = makeFeatSelControlRandom(maxit = 20)
lrn1 = makeFeatSelWrapper(lrn.rf, resampling = inner, control = ctrl, show.info=TRUE)
parallelMap::parallelStartSocket(parallel::detectCores())
glc.mod1 = train(lrn1, task = tsk.C, weights = case.weights)
parallelMap::parallelStop()
glc.sfeats1 = getFeatSelResult(glc.mod1)
str(glc.sfeats1$x)
## 109
lrn.xg = mlr::makeLearner("classif.xgboost")
lrn2 = makeFeatSelWrapper(lrn.xg, resampling = inner, control = ctrl, show.info=TRUE)
parallelMap::parallelStartSocket(parallel::detectCores())
glc.mod2 = train(lrn2, task = tsk.C, weights = case.weights)
parallelMap::parallelStop()
glc.sfeats2 = getFeatSelResult(glc.mod2)
str(glc.sfeats2$x)
## new shorter formula
formulaString.GLC0 = as.formula(paste('map_code_pnv.f ~ ', paste(unique(c(glc.sfeats1$x, glc.sfeats2$x)), collapse="+")))
length(all.vars(formulaString.GLC0))
## 163

## GLC final model ----
SL.library <- c("classif.ranger", "classif.xgboost", "classif.nnTrain")
## "classif.ksvm", "classif.kknn", "classif.nnet"
## "classif.multinom" --> 00005: Error in nnet.default(X, Y, w, mask = mask, size = 0, skip = TRUE, softmax = TRUE,  : too many (3075) weights
parallelMap::parallelStartSocket(parallel::detectCores())
tsk.C0 <- mlr::makeClassifTask(data = rm.lc[,all.vars(formulaString.GLC0)], target = all.vars(formulaString.GLC0)[1]) ## weights = case.weights 
#Error in checkLearnerBeforeTrain(task, learner, weights) : 
#  Weights vector passed to train, but learner 'stack' does not support that!
lrns <- list(lrn.rf, mlr::makeLearner(SL.library[2], verbose=1), mlr::makeLearner(SL.library[3]))
lrns <- lapply(lrns, setPredictType, "prob")
init.m <- mlr::makeStackedLearner(base.learners = lrns, predict.type = "prob", method = "stack.cv", super.learner = "classif.glmnet")
## takes 10+ minutes
system.time( m.C <- mlr::train(init.m, tsk.C0) ) ## weights = case.weights
#user  system elapsed 
#715.030   5.831 374.824 
parallelMap::parallelStop()
saveRDS.gz(m.C, "./models/eml_lc.rds")
save.image.pigz(n.cores = 64)

## GLC predictions ----
#x = list.files(path="/data/tt/LandGIS/grid250m", pattern=glob2rx("potential.landcover_C_*.tif$"), recursive=TRUE, full.names = TRUE)
#x = list.files(path="/data/tt/LandGIS/grid250m", pattern=glob2rx("^potential.landcover_*.tif$"), recursive=TRUE, full.names = TRUE)
#unlink(x)
#pred_probs.mlr(i="T39797", m.C, tile.tbl, col.legend=esa.tbl[,c("response","Number")], varn="potential.landcover", out.dir="/data/tt/LandGIS/grid250m")
#pred_probs.mlr(i="T39798", m.C, tile.tbl, col.legend=esa.tbl[,c("response","Number")], varn="potential.landcover", out.dir="/data/tt/LandGIS/grid250m")
#pred_probs.mlr(i="T38707", m.C, tile.tbl, col.legend=esa.tbl[,c("response","Number")], varn="potential.landcover", out.dir="/data/tt/LandGIS/grid250m")
system.time( pred_probs.mlr(i="T38715", m.C, tile.tbl, col.legend=esa.tbl[,c("response","Number")], varn="potential.landcover", out.dir="/data/tt/LandGIS/grid250m") )
## 100 secs
## total hours:
19000*77/60/60/62

## TAKES >10hrs
#cpus = unclass(round((400-25)/(3*(object.size(m)/1e9))))
library(parallel)
library(snowfall)
sfInit(parallel=TRUE, cpus=60)
sfExport("m.C", "pred_probs.mlr", "tile.tbl", "esa.tbl", "pr.dirs")
sfLibrary(plyr)
sfLibrary(ranger)
sfLibrary(xgboost)
sfLibrary(deepnet)
sfLibrary(mlr)
sfLibrary(rgdal)
sfLibrary(stats)
out <- sfClusterApplyLB(pr.dirs, function(i){ try( pred_probs.mlr(i, m.C, tile.tbl, col.legend=esa.tbl[,c("response","Number")], varn="potential.landcover", out.dir="/data/tt/LandGIS/grid250m") ) })
sfStop()

## GLC 250m mosaics ----
r = raster("/data/LandGIS/layers250m/lcv_admin0_fao.gaul_c_250m_s0..0cm_2015_v1.0.tif")
te = as.vector(extent(r))[c(1,3,2,4)]
cellsize = res(r)[1]

#x = list.files("./predicted250m", pattern=glob2rx("pnv_*_probav.lc100_*.tif$"), full.names = TRUE)
#unlink(x)
df.n = expand.grid(m.C$task.desc$class.levels, c("_p_", "_sd_"))
names(df.n) = c("response", "prob")
df.n$varn = paste0("potential.landcover", ifelse(df.n$prob=="_p_", "_M", "_sd"))
df.n$land_cover_class = tolower(make.names(join(df.n, esa.tbl)$land_cover_class))
filename = paste0("./predicted250m/pnv_potential.landcover_probav.lc100.", df.n$land_cover_class, df.n$prob, "250m_s0..0cm_2017_v0.1.tif")
#View(df.n)
save.image.pigz(n.cores = 64)

library(snowfall)
sfInit(parallel=TRUE, cpus=20) ## length(filename)
sfExport("df.n", "mosaick_ll", "filename", "te", "cellsize")
out <- sfClusterApplyLB(1:nrow(df.n), function(x){ try( mosaick_ll(varn=df.n$varn[x], i=paste0("prob.", df.n$response[x]), out.tif=filename[x], in.path="/data/tt/LandGIS/grid250m", out.path="/data/Geo-wiki/predicted250m", tr=cellsize, te=paste(te, collapse = " "), ot="Byte", dstnodata=255, aggregate=FALSE) )})
sfStop()
## most probable class
mosaick_ll(varn="potential.landcover", i="C", out.tif="/data/Geo-wiki/predicted250m/pnv_potential.landcover_probav.lc100_c_250m_s0..0cm_2017_v0.1.tif", dominant = TRUE, in.path="/data/tt/LandGIS/grid250m", tr=cellsize, te=paste(te, collapse = " "), ot="Byte", dstnodata=255, aggregate=FALSE)
write.csv(esa.tbl[,c("map_code","land_cover_class","pv_mapping","land_cover_class_agg")], "/data/Geo-wiki/predicted250m/pnv_potential.landcover_probav.lc100_c_250m_s0..0cm_2017_v0.1.tif.csv")

## Biomes ----
biome_f = readRDS("/mnt/DATA/PNV/Data/Biomes/biome_f.rds")
## 8797 points
ov.biome <- extract.tiled(obj=biome_f, tile.pol=tile.pol, path="/data/tt/LandGIS/grid250m", ID="ID", cpus=64)
summary(ov.biome$Biome00k_c)
saveRDS.gz(ov.biome, "./training_points/ov_biome_natural_vegetation.pnts.rds")
#ov.biome = readRDS.gz("./training_points/ov_biome_natural_vegetation.pnts.rds")
rm.biome = ov.biome[complete.cases(ov.biome[,pr.vars]),]
dim(rm.biome)
## 7633   366
formulaString.biome = as.formula(paste('Biome00k_c ~ ', paste(pr.vars, collapse="+")))
saveRDS.gz(rm.biome, "./training_points/regression.matrix_biome.rds")

## Biomes parameter tuning ----
library(mlr)
tsk.B <- mlr::makeClassifTask(data = rm.biome[,all.vars(formulaString.biome)], target = all.vars(formulaString.biome)[1])
discrete_ps = makeParamSet( makeDiscreteParam("mtry", values = seq(5,80,by=5)) )
ctrl = makeTuneControlGrid()
rdesc = makeResampleDesc("CV", iters = 3L)
parallelMap::parallelStartSocket(parallel::detectCores())
res = tuneParams("classif.ranger", task = tsk.B, resampling = rdesc, par.set = discrete_ps, control = ctrl)
parallelMap::parallelStop()
res
## Biomes Feature selection ----
outer = makeResampleDesc("CV", iters = 3L)
inner = makeResampleDesc("Holdout")
ctrl = makeFeatSelControlRandom(maxit = 20)
lrn.rf = mlr::makeLearner("classif.ranger", num.threads = parallel::detectCores(), mtry=res$x$mtry, num.trees=85)
lrn1 = makeFeatSelWrapper(lrn.rf, resampling = inner, control = ctrl, show.info=TRUE)
parallelMap::parallelStartSocket(parallel::detectCores())
biome.mod1 = train(lrn1, task = tsk.B)
parallelMap::parallelStop()
biome.sfeats1 = getFeatSelResult(biome.mod1)
str(biome.sfeats1$x)
## 103
lrn.mn = mlr::makeLearner("classif.glmnet")
lrn2 = makeFeatSelWrapper(lrn.mn, resampling = inner, control = ctrl, show.info=TRUE)
parallelMap::parallelStartSocket(parallel::detectCores())
biome.mod2 = train(lrn2, task = tsk.B)
parallelMap::parallelStop()
biome.sfeats2 = getFeatSelResult(biome.mod2)
str(biome.sfeats2$x)
## new shorter formula
formulaString.biome0 = as.formula(paste('Biome00k_c ~ ', paste(unique(c(biome.sfeats1$x, biome.sfeats2$x)), collapse="+")))
length(all.vars(formulaString.biome0))
## 157
tsk.B0 <- mlr::makeClassifTask(data = rm.biome[,all.vars(formulaString.biome0)], target = all.vars(formulaString.biome0)[1])

## Biomes mlr ----
SL2.library <- c("classif.ranger", "classif.glmnet", "classif.xgboost", "classif.nnTrain") 
## "classif.kknn" <- very slow for generating predictions
lrns2 <- list(lrn.rf, mlr::makeLearner(SL2.library[2]), mlr::makeLearner(SL2.library[3], verbose=1), mlr::makeLearner(SL2.library[4]))
lrns2 <- lapply(lrns2, setPredictType, "prob")
parallelMap::parallelStartSocket(parallel::detectCores())
init2.m <- mlr::makeStackedLearner(base.learners = lrns2, predict.type = "prob", method = "stack.cv", super.learner = "classif.glmnet") #, method="compress")
## takes 10+ minutes
system.time( m.B <- mlr::train(init2.m, tsk.B0) )
#[1] train-merror:0.310394
#user  system elapsed 
#316.962   1.502 544.755
## one multinomial or binomial class has fewer than 8  observations; dangerous ground
## mapping accuracy
test.set = seq(5, nrow(rm.biome), by = 5)
training.set = which(!1:nrow(rm.biome) %in% test.set)
rm.biome.t = rm.biome[training.set, all.vars(formulaString.biome0)]
tsk.B0t = mlr::makeClassifTask(data = rm.biome.t, target = all.vars(formulaString.biome0)[1])
m.Bt <- mlr::train(init2.m, tsk.B0t)
x = predict(m.Bt, newdata=rm.biome[test.set, all.vars(formulaString.biome0)[-1]])
x$data$truth = rm.biome$Biome00k_c[test.set]
xp = mlr::performance(x, measures=list(acc, logloss))
xp
##       acc   logloss 
## 0.6610059 1.1297139
## 0.6583932 1.1320812
parallelMap::parallelStop()
saveRDS.gz(m.B, "./models/eml_biome.rds")
save.image.pigz(n.cores = 64)

mrfX_bm <- ranger::ranger(formulaString.biome0, rm.biome, importance="impurity", probability=TRUE, num.trees=85, mtry=res$x$mtry) 
mrfX_bm
## OOB prediction error (Brier s.):  0.3096173
saveRDS.gz(mrfX_bm, "./models/mrfX_bm.rds")
#mrfX_bm = readRDS.gz("./models/mrfX_bm.rds")
xb <- as.list(ranger::importance(mrfX_bm))
print(t(data.frame(xb[order(unlist(xb), decreasing=TRUE)[1:25]])))
# clm_precipitation_sm2rain.may_m_1km_s0..0cm_2007..2018_v0.2.tif            267.78061
# clm_monthly.temp_worldclim.chelsa.feb_u.99_1km_s0..0cm_1979..2013_v1.0.tif 199.06684
# clm_monthly.temp_worldclim.chelsa.dec_u.99_1km_s0..0cm_1979..2013_v1.0.tif 176.32663
# clm_precipitation_sm2rain.annual_m_1km_s0..0cm_2007..2018_v0.2.tif         159.17973
# clm_diffuse.irradiation_solar.atlas.kwhm2.100_m_1km_s0..0cm_2016_v1.tif    150.86225
# clm_monthly.temp_worldclim.chelsa.jan_u.99_1km_s0..0cm_1979..2013_v1.0.tif 136.70924
# clm_bioclim.var_chelsa.2_m_1km_s0..0cm_1979..2013_v1.0.tif                 136.41685
# clm_water.vapor_nasa.eo.nov_m_1km_s0..0cm_2000..2017_v1.0.tif              112.57023
# clm_bioclim.var_chelsa.3_m_1km_s0..0cm_1979..2013_v1.0.tif                 101.02893

biome.tbl = read.csv("/mnt/DATA/PNV/Data/Biomes/Biome_legend.csv")
biome.tbl = biome.tbl[!duplicated(biome.tbl$New.global.consolidated.biome.scheme)&!biome.tbl$New.global.consolidated.biome.scheme=="",]
biome.tbl$response = make.names(biome.tbl$New.global.consolidated.biome.scheme)
## 20 classes
biome.tax <- m.B$task.desc$class.levels
biome.tbl$response[which(!biome.tbl$response %in% biome.tax)]
#x = list.files(path="/data/tt/LandGIS/grid250m", pattern=glob2rx("biome.type_C_*.tif$"), recursive=TRUE, full.names = TRUE)
#x = list.files(path="/data/tt/LandGIS/grid250m", pattern=glob2rx("biome.type_*.tif$"), recursive=TRUE, full.names = TRUE)
#unlink(x)
## test:
system.time( pred_probs.mlr(i="T38715", m.B, tile.tbl, col.legend=biome.tbl[,c("response","Number")], varn="biome.type", out.dir="/data/tt/LandGIS/grid250m") )
system.time( pred_probs.mlr(i="T38716", m.B, tile.tbl, col.legend=biome.tbl[,c("response","Number")], varn="biome.type", out.dir="/data/tt/LandGIS/grid250m") )
## 100 secs
## total hours:
19000*100/60/60/62
pred_probs.mlr(i="T45733", m.B, tile.tbl, col.legend=biome.tbl[,c("response","Number")], varn="biome.type", out.dir="/data/tt/LandGIS/grid250m")

## TAKES >12hrs
library(parallel)
library(snowfall)
sfInit(parallel=TRUE, cpus=62)
sfExport("m.B", "pred_probs.mlr", "tile.tbl", "biome.tbl", "pr.dirs")
sfLibrary(plyr)
sfLibrary(ranger)
sfLibrary(xgboost)
sfLibrary(kknn)
sfLibrary(glmnet)
sfLibrary(deepnet)
sfLibrary(mlr)
sfLibrary(rgdal)
sfLibrary(stats)
out <- sfClusterApplyLB(pr.dirs, function(i){ try( pred_probs.mlr(i, m.B, tile.tbl, col.legend=biome.tbl[,c("response","Number")], varn="biome.type", out.dir="/data/tt/LandGIS/grid250m") ) })
sfStop()

#x = list.files("./predicted250m", pattern=glob2rx("pnv_biome.*.tif$"), full.names = TRUE)
#unlink(x)
df2.n = expand.grid(m.B$task.desc$class.levels, c("_p_", "_sd_"))
names(df2.n) = c("response", "prob")
df2.n$varn = paste0("biome.type", ifelse(df2.n$prob=="_p_", "_M", "_sd"))
filename2 = paste0("./predicted250m/pnv_biome.type_biome00k.", df2.n$response, df2.n$prob, "250m_s0..0cm_2000..2017_v0.2.tif")
## pnv_biome.type_biome00k.graminoid.and.forb.tundra_p_250m_s0..0cm_2000..2017_v0.2.tif
save.image.pigz(n.cores = 64)

library(snowfall)
sfInit(parallel=TRUE, cpus=20) ## length(filename)
sfExport("df2.n", "mosaick_ll", "filename2", "te", "cellsize")
out <- sfClusterApplyLB(1:nrow(df2.n), function(x){ try( mosaick_ll(varn=df2.n$varn[x], i=paste0("prob.", df2.n$response[x]), out.tif=filename2[x], in.path="/data/tt/LandGIS/grid250m", out.path="/data/Geo-wiki/predicted250m", tr=cellsize, te=paste(te, collapse = " "), ot="Byte", dstnodata=255, aggregate=FALSE) )})
sfStop()
## most probable class
mosaick_ll(varn="biome.type", i="C", out.tif="/data/Geo-wiki/predicted250m/pnv_biome.type_biome00k_c_250m_s0..0cm_2000..2017_v0.2.tif", dominant = TRUE, in.path="/data/tt/LandGIS/grid250m", tr=cellsize, te=paste(te, collapse = " "), ot="Byte", dstnodata=255, aggregate=FALSE)
#file.copy("/mnt/DATA/LandGIS/predicted1km/pnv_biome.type_biome00k_c_1km_s0..0cm_2000..2017_v0.1.tif.csv", "/data/Geo-wiki/predicted250m/pnv_biome.type_biome00k_c_250m_s0..0cm_2000..2017_v0.2.tif.csv")
