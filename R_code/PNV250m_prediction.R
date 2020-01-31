## Mapping potential natural vegetation mapping at 250 m
## Code by: tom.hengl@envirometrix.net & jung@iiasa.ac.at 
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
lc100m = raster("/mnt/GIS/Copernicus_vito/GLC100m/discrete-classification.vrt")

## GLC training points ----
gs.url = as_id("https://docs.google.com/spreadsheets/d/10Y8PryG5TiPZTiArteRXE0iFdKUhNYhj38qb9bRnpVU/edit")
xlsxFile = "Global_forest_management_and_LC_mapping_tables.xlsx"
drive_download(gs.url, xlsxFile, overwrite = TRUE)
esa.tbl = openxlsx::read.xlsx(xlsxFile, sheet = 1)
head(esa.tbl)

## tiling system
tile.tbl = readRDS("/data/LandGIS/models/stacked250m_tiles.rds")
pr.dirs = basename(dirname(list.files(path="/data/tt/LandGIS/grid250m", pattern=glob2rx("*.rds$"), recursive=TRUE)))
tile.pol = readOGR("/data/LandGIS/models/tiles_ll_100km.shp")
tile.pol = tile.pol[paste0("T", tile.pol$ID) %in% pr.dirs,]

## GLC overlay ---- 
tot_pnts = readRDS("glc100m_natural_vegetation.pnts.rds")
## (takes 30 mins):
ov.lc <- extract.tiled(obj=tot_pnts, tile.pol=tile.pol, path="/data/tt/LandGIS/grid250m", ID="ID", cpus=64)
summary(!is.na(ov.lc$map_code_pnv))
ov.lc = ov.lc[!is.na(ov.lc$map_code_pnv),]
ov.lc$map_code_pnv.f = as.factor(make.names(ov.lc$map_code_pnv))
summary(ov.lc$map_code_pnv.f)
#head(ov.lc)
pr.vars = make.names(unlist(sapply(c("sm2rain","monthly.temp_worldclim.chelsa","bioclim.var_chelsa","irradiation_solar.atlas", "usgs.ecotapestry", "floodmap.500y", "water.table.depth_deltares", "snow.prob_esacci", "water.vapor_nasa.eo", "wind.speed_terraclimate", "merit.dem_m", "merit.hydro_m", "cloud.fraction_earthenv", "water.occurance_jrc", "wetlands.cw_upmc", "pb2002"), function(i){names(ov.lc)[grep(i, names(ov.lc))]})))
str(pr.vars)
## 230
saveRDS.gz(ov.lc, "./training_points/ov_glc100m_natural_vegetation.pnts.rds")
#ov.lc = readRDS.gz("./training_points/ov_glc100m_natural_vegetation.pnts.rds")

## Subset data ----
rm.lc = ov.lc[complete.cases(ov.lc[,pr.vars]),]
dim(rm.lc)
## 64293   403
formulaString.GLC = as.formula(paste('map_code_pnv.f ~ ', paste(pr.vars, collapse="+")))
saveRDS.gz(rm.lc, "./training_points/regression.matrix_tot_pnts.rds")
esa.tbl$response = make.names(esa.tbl$map_code)
esa.tbl$Number = esa.tbl$map_code
save.image.pigz(n.cores = 64)

## GLC testing RF ----
## Geographically distributed sample:
pnts.s <- GSIF::sample.grid(tot_pnts, cell.size=c(1,1), n=1)
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
# 21263 samples
# 230 predictor
# mtry  Accuracy   Kappa    
# 15   0.6432283  0.5933275
# 40   0.6448278  0.5954506
# 80   0.6458154  0.5967604
# 110   0.6468968  0.5979766
# 130   0.6438399  0.5946539
# 150   0.6436049  0.5944365

## GLC ranger model ----
summary(as.factor(rm.lc$source_db))
case.weights = ifelse(rm.lc$source_db=="MangrovesDB" | rm.lc$source_db=="Biome6000" | rm.lc$source_db=="LandPKS_land_cover" | rm.lc$source_db=="Phenocam" | rm.lc$source_db=="Cheng et al." | rm.lc$source_db=="Marinova et al." | rm.lc$source_db=="Lezine et al.", 10, ifelse(rm.lc$source_db=="Geo-wiki_HIF" | rm.lc$source_db=="GeoFeedback" | rm.lc$source_db=="Land_Cover_maps", 3, 1))
summary(as.factor(case.weights))
mrfX_lc <- ranger::ranger(formulaString.GLC, rm.lc, importance="impurity", mtry=t.mrfX0$bestTune$mtry, probability=TRUE, num.trees=85, case.weights=case.weights) 
mrfX_lc
# Type:                             Probability estimation 
# Number of trees:                  85 
# Sample size:                      69259 
# Number of independent variables:  230 
# Mtry:                             110 
# Target node size:                 10 
# Variable importance mode:         impurity 
# Splitrule:                        gini 
# OOB prediction error (Brier s.):  0.3039445
saveRDS.gz(mrfX_lc, "./models/mrfX_lc.rds")
#mrfX_lc = readRDS.gz("./models/mrfX_lc.rds")
xl <- as.list(ranger::importance(mrfX_lc))
print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:25]])))
# clm_water.vapor_nasa.eo.may_m_1km_s0..0cm_2000..2017_v1.0.tif              2109.1795
# dtm_elevation_merit.dem_m_250m_s0..0cm_2017_v1.0.tif                       1982.5696
# clm_bioclim.var_chelsa.2_m_1km_s0..0cm_1979..2013_v1.0.tif                 1539.9564
# clm_water.vapor_nasa.eo.oct_m_1km_s0..0cm_2000..2017_v1.0.tif              1393.8601
# clm_precipitation_sm2rain.annual_m_1km_s0..0cm_2007..2018_v0.2.tif         1059.3181
# clm_monthly.temp_worldclim.chelsa.oct_u.99_1km_s0..0cm_1979..2013_v1.0.tif  983.3433
# clm_precipitation_sm2rain.may_m_1km_s0..0cm_2007..2018_v0.2.tif             827.1182
# clm_diffuse.irradiation_solar.atlas.kwhm2.100_m_1km_s0..0cm_2016_v1.tif     681.0933
# clm_monthly.temp_worldclim.chelsa.apr_m_1km_s0..0cm_1979..2013_v1.0.tif     660.9366
# clm_bioclim.var_chelsa.7_m_1km_s0..0cm_1979..2013_v1.0.tif                  619.6422
# clm_bioclim.var_chelsa.3_m_1km_s0..0cm_1979..2013_v1.0.tif                  612.9946
# clm_bioclim.var_chelsa.12_m_1km_s0..0cm_1979..2013_v1.0.tif                 605.2164
# clm_water.vapor_nasa.eo.dec_m_1km_s0..0cm_2000..2017_v1.0.tif               602.4202
# clm_monthly.temp_worldclim.chelsa.feb_u.99_1km_s0..0cm_1979..2013_v1.0.tif  600.7344
# clm_bioclim.var_chelsa.4_m_1km_s0..0cm_1979..2013_v1.0.tif                  591.6802
# clm_precipitation_sm2rain.jun_m_1km_s0..0cm_2007..2018_v0.2.tif             576.1666

library(mlr)
## generic parameters
discrete_ps = makeParamSet( makeDiscreteParam("mtry", values = seq(5,80,by=5)) )
xg.model_Params <- makeParamSet(
  makeDiscreteParam("nrounds", value=c(10,50)),
  makeDiscreteParam("max_depth", value=c(1,2,3,4)),
  makeDiscreteParam("eta", value=c(0.3,0.4,0.5)),
  makeDiscreteParam("subsample", value=c(1)),
  makeDiscreteParam("min_child_weight", value=c(1,2,3)),
  makeDiscreteParam("colsample_bytree", value=c(0.8))
)
rdesc = makeResampleDesc("CV", iters = 3L)

tsk.C <- mlr::makeClassifTask(data = rm.lc[,all.vars(formulaString.GLC)], target = all.vars(formulaString.GLC)[1])
tsk.Cs <- mlr::makeClassifTask(data = df[,all.vars(formulaString.GLC)], target = all.vars(formulaString.GLC)[1])
parallelMap::parallelStartSocket(parallel::detectCores())
resC = tuneParams(mlr::makeLearner("classif.ranger", num.threads = round(parallel::detectCores()/16), num.trees=85), task = tsk.C, resampling = rdesc, par.set = discrete_ps, control = makeTuneControlGrid())
parallelMap::parallelStop()
resC
# Tune result:
#   Op. pars: mtry=55
# mmce.test.mean=0.3169119
## GLC Feature selection ----
## TAKES >15mins
lrn.rf = mlr::makeLearner("classif.ranger", num.threads = parallel::detectCores(), mtry=resC$x$mtry, num.trees=85, predict.type="prob")
#lrn.rf
outer = makeResampleDesc("CV", iters = 3L)
inner = makeResampleDesc("Holdout")
ctrl = makeFeatSelControlRandom(maxit = 20)
lrn1 = makeFeatSelWrapper(lrn.rf, resampling = inner, control = ctrl, show.info=TRUE)
parallelMap::parallelStartSocket(parallel::detectCores())
glc.mod1 = mlr::train(lrn1, task = tsk.C, weights = case.weights)
parallelMap::parallelStop()
glc.sfeats1 = getFeatSelResult(glc.mod1)
str(glc.sfeats1$x)
## 111
## https://www.kaggle.com/xanderhorn/train-r-ml-models-efficiently-with-mlr
parallelMap::parallelStartSocket(parallel::detectCores())
resC2 = tuneParams(mlr::makeLearner("classif.xgboost"), task = tsk.C, resampling = rdesc, par.set = xg.model_Params, control = makeTuneControlRandom(maxit = 2L)) # control = makeTuneControlGrid())
## [Tune] Result: nrounds=10; max_depth=3; eta=0.5; subsample=1; min_child_weight=1; colsample_bytree=0.8 : mmce.test.mean=0.3916169
lrn.xg = mlr::makeLearner("classif.xgboost")
lrn.xg = setHyperPars(lrn.xg, par.vals = resC2$x)
lrn2 = makeFeatSelWrapper(lrn.xg, resampling = inner, control = ctrl, show.info=TRUE)
glc.mod2 = mlr::train(lrn2, task = tsk.C, weights = case.weights)
parallelMap::parallelStop()
glc.sfeats2 = getFeatSelResult(glc.mod2)
str(glc.sfeats2$x)
## 106
## new shorter formula
formulaString.GLC0 = as.formula(paste('map_code_pnv.f ~ ', paste(unique(c(glc.sfeats1$x, glc.sfeats2$x)), collapse="+")))
length(all.vars(formulaString.GLC0))
## 178

## GLC final model ----
SL.library <- c("classif.ranger", "classif.xgboost", "classif.nnTrain")
## "classif.ksvm", "classif.kknn", "classif.nnet"
## "classif.multinom" --> 00005: Error in nnet.default(X, Y, w, mask = mask, size = 0, skip = TRUE, softmax = TRUE,  : too many (3075) weights
parallelMap::parallelStartSocket(parallel::detectCores())
tsk.C0 <- mlr::makeClassifTask(data = rm.lc[,all.vars(formulaString.GLC0)], target = all.vars(formulaString.GLC0)[1]) ## weights = case.weights 
#Error in checkLearnerBeforeTrain(task, learner, weights) : 
#  Weights vector passed to train, but learner 'stack' does not support that!
lrns <- list(lrn.rf, lrn.xg, mlr::makeLearner(SL.library[3]))
lrns <- lapply(lrns, setPredictType, "prob")
init.m <- mlr::makeStackedLearner(base.learners = lrns, predict.type = "prob", method = "stack.cv", super.learner = "classif.glmnet")
## takes 10+ minutes
system.time( m.C <- mlr::train(init.m, tsk.C0) ) ## weights = case.weights
#user  system elapsed 
#715.030   5.831 374.824 
parallelMap::parallelStop()
saveRDS.gz(m.C, "./models/eml_lc.rds")
save.image.pigz(n.cores = 64)
m.C = readRDS.gz("./models/eml_lc.rds")

## GLC predictions ----
#x = list.files(path="/data/tt/LandGIS/grid250m", pattern=glob2rx("potential.landcover_C_*.tif$"), recursive=TRUE, full.names = TRUE)
#library(pbmcapply)
#xi = pbmclapply(x, FUN=function(i){file.info(i)$ctime}, mc.cores = parallel::detectCores())
#hist(unlist(xi))
#unlink(x[unlist(xi)<1575960000])
#x = list.files(path="/data/tt/LandGIS/grid250m", pattern=glob2rx("^potential.landcover_*.tif$"), recursive=TRUE, full.names = TRUE)
#unlink(x)
#pred_probs.mlr(i="T39797", m.C, tile.tbl, col.legend=esa.tbl[,c("response","Number")], varn="potential.landcover", out.dir="/data/tt/LandGIS/grid250m")
#pred_probs.mlr(i="T39798", m.C, tile.tbl, col.legend=esa.tbl[,c("response","Number")], varn="potential.landcover", out.dir="/data/tt/LandGIS/grid250m")
#pred_probs.mlr(i="T38707", m.C, tile.tbl, col.legend=esa.tbl[,c("response","Number")], varn="potential.landcover", out.dir="/data/tt/LandGIS/grid250m")
system.time( pred_probs.mlr(i="T38715", m.C, tile.tbl, col.legend=esa.tbl[,c("response","Number")], varn="potential.landcover", out.dir="/data/tt/LandGIS/grid250m") )
## 80 secs
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

#x = list.files("./predicted250m", pattern=glob2rx("pnv_*_probav.lc100*.tif$"), full.names = TRUE)
#unlink(x)
df.n = expand.grid(m.C$task.desc$class.levels, c("_p_", "_sd_"))
names(df.n) = c("response", "prob")
df.n$varn = paste0("potential.landcover", ifelse(df.n$prob=="_p_", "_M", "_sd"))
df.n$land_cover_class = tolower(make.names(plyr::join(df.n, esa.tbl)$land_cover_class))
filename = paste0("./predicted250m/pnv_potential.landcover_probav.lc100.", df.n$land_cover_class, df.n$prob, "250m_s0..0cm_2017_v0.1.tif")
View(df.n)
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
biome.mod1 = mlr::train(lrn1, task = tsk.B)
parallelMap::parallelStop()
biome.sfeats1 = getFeatSelResult(biome.mod1)
str(biome.sfeats1$x)
## 103
lrn.mn = mlr::makeLearner("classif.glmnet")
lrn2 = makeFeatSelWrapper(lrn.mn, resampling = inner, control = ctrl, show.info=TRUE)
parallelMap::parallelStartSocket(parallel::detectCores())
biome.mod2 = mlr::train(lrn2, task = tsk.B)
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

## Canopy height ----
## ICESat-2 points
canopy_h = readRDS.gz("/mnt/DATA/NSIDC/icesat2_heights.rds")
## 260M points
summary(canopy_h$h_canopy)
## 1 to 150 m -- missing values are 0 canopy?
summary(canopy_h$h_canopy_uncertainty)
## 0 to 1435 with median = 4 / mean = 11
#sel.can = (canopy_h$h_canopy_uncertainty < 10 & !is.na(canopy_h$h_canopy))
sel.can = !is.na(canopy_h$h_canopy)
#summary(sel.can)
canopy_h = canopy_h[sel.can,]
dim(canopy_h)
# 32M points
saveRDS.gz(canopy_h, "./training_points/icesat2_h_canopy.rds")
## subset to native vegetation:
mask.lst = c(paste0("/mnt/DATA/protectedplanet/", c("ifl_2013.tif", "ifl_2000.tif", "WDPA_Dec2017-shapefile-polygons.tif")), "/mnt/GIS/Copernicus_vito/GLC100m/discrete-classification.vrt")
## takes >2hrs
can.test = sampleInt(nrow(canopy_h), size = 1.5e6)
## with 20M points - takes over 24hrs
system.time( ov.can <- parallel::mclapply(mask.lst, function(i){raster::extract(raster::raster(i), canopy_h[can.test,c("lon","lat")])}, mc.cores = length(mask.lst)) )
ov.can = data.frame(ov.can)
names(ov.can) = basename(mask.lst)
summary(as.factor(ov.can$`discrete-classification.vrt`))
## NA?
sel.lc = is.na(ov.can$`discrete-classification.vrt`) | ov.can$`discrete-classification.vrt` %in% esa.tbl$map_code_pnv[!is.na(esa.tbl$map_code_pnv)] | !is.na(ov.can$ifl_2013.tif) | !is.na(ov.can$ifl_2000.tif) | !is.na(ov.can$`WDPA_Dec2017-shapefile-polygons.tif`)
#sel.lc = ov.can$`discrete-classification.vrt` %in% c(esa.tbl$map_code_pnv[!is.na(esa.tbl$map_code_pnv)])
summary(sel.lc)
## 1,160,574 points
canopy_h.s = canopy_h[can.test,][sel.lc,]
canopy_h.s$glc_class = ov.can$`discrete-classification.vrt`[sel.lc]
canopy_h.s[30:35,]
## bare land = 0 canopy height
bareland = as.data.frame(tot_pnts[which(tot_pnts$map_code_pnv==60),])
## 4502 points
canopy_h.sb = dplyr::bind_rows(list(canopy_h.s, data.frame(h_canopy=0, lat=bareland$latitude_decimal_degrees, lon=bareland$longitude_decimal_degrees, h_canopy_uncertainty=2)))
coordinates(canopy_h.sb) = ~ lon + lat
proj4string(canopy_h.sb) = tot_pnts@proj4string
unlink("h_canopy_natural_vegetation.pnts.gpkg")
writeOGR(canopy_h.sb, "h_canopy_natural_vegetation.pnts.gpkg", "h_canopy_natural_vegetation.pnts", "GPKG")
dim(canopy_h.sb)

ov.canopy <- extract.tiled(obj=canopy_h.sb, tile.pol=tile.pol, path="/data/tt/LandGIS/grid250m", ID="ID", cpus=64)
summary(ov.canopy$h_canopy)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   4.734  10.323  13.721  17.435 149.771
ov.canopy$log.h_canopy = log1p(ov.canopy$h_canopy)
#summary(ov.canopy$log.h_canopy)
saveRDS.gz(ov.canopy, "./training_points/ov_canopy_natural_vegetation.pnts.rds")
#ov.canopy = readRDS.gz("./training_points/ov_canopy_natural_vegetation.pnts.rds")
rm.canopy = ov.canopy[complete.cases(ov.canopy[,pr.vars]),]
## TH: remove suspicious points e.g. high canopy height for bare land etc
## far nother latitudes show values >100 m --> artifacts?!
hist(rm.canopy$veg_fapar_proba.v.annnual_d_250m_s0..0cm_2014..2019_v1.0.tif, breaks=35)
hist(rm.canopy$lcv_bareground_usgs.landsat_p_250m_s0..0cm_2010_v1.0.tif, breaks=35)
sel.rm.a = (rm.canopy$veg_fapar_proba.v.annnual_d_250m_s0..0cm_2014..2019_v1.0.tif<25 | rm.canopy$Y > 75 | rm.canopy$lcv_bareground_usgs.landsat_p_250m_s0..0cm_2010_v1.0.tif > 70) & rm.canopy$h_canopy>2
summary(sel.rm.a)
rm.canopy = rm.canopy[-which(sel.rm.a),]
dim(rm.canopy)
## 1125395     408

formulaString.canopy = as.formula(paste('log.h_canopy ~ ', paste(pr.vars, collapse="+")))
saveRDS.gz(rm.canopy, "./training_points/regression.matrix_canopy.rds")
rm(ov.canopy); rm(canopy_h)
save.image.pigz(n.cores = 64)

## Canopy parameter tuning ----
library(mlr)
discrete_ps = makeParamSet(makeDiscreteParam("mtry", values = seq(5,80,by=5)))
ctrl = makeTuneControlGrid()
rdesc = makeResampleDesc("CV", iters = 2L)
outer = makeResampleDesc("CV", iters = 2L)
inner = makeResampleDesc("Holdout")
ctrlF = makeFeatSelControlRandom(maxit = 20)
xg.model_Params <- makeParamSet(
  makeDiscreteParam("nrounds", value=c(10,50)),
  makeDiscreteParam("max_depth", value=c(1,2,3,4)),
  makeDiscreteParam("eta", value=c(0.3,0.4,0.5)),
  makeDiscreteParam("subsample", value=c(1)),
  makeDiscreteParam("min_child_weight", value=c(1,2,3)),
  makeDiscreteParam("colsample_bytree", value=c(0.8))
)
## final RM
rC.sel = complete.cases(rm.canopy[,c(all.vars(formulaString.canopy), "X", "Y")])
rm.canopym = rm.canopy[rC.sel,all.vars(formulaString.canopy)]
## use only 40,000 training points for XGBoost --> should be enough to estimate the model pars
## otherwise R crashes / RAM problems 
rm.test = sampleInt(nrow(rm.canopym), size = 4e4)
rm.canopym.s = rm.canopym[rm.test,]
tsk0.s <- mlr::makeRegrTask(data = rm.canopym.s, target = "log.h_canopy", coordinates = rm.canopy[rC.sel,][rm.test,c("X","Y")])
tsk0 <- mlr::makeRegrTask(data = rm.canopym, target = "log.h_canopy", coordinates = rm.canopy[rC.sel,c("X","Y")])
parallelMap::parallelStartSocket(parallel::detectCores())
resR.rf = tuneParams(mlr::makeLearner("regr.ranger", num.threads = round(parallel::detectCores()/length(discrete_ps$pars$mtry$values)), num.trees=85), task = tsk0, resampling = rdesc, par.set = discrete_ps, control = ctrl)
#Mapping in parallel: mode = socket; level = mlr.tuneParams; cpus = 64; elements = 16.
#[Tune] Result: mtry=25 : mse.test.mean=0.2789476
## feature selection
lrn.rf = mlr::makeLearner("regr.ranger", num.threads = parallel::detectCores(), mtry=resR.rf$x$mtry, num.trees=85, importance="impurity")
lrn1 = makeFeatSelWrapper(lrn.rf, resampling = inner, control = ctrlF, show.info=TRUE)
var.mod1 = mlr::train(lrn1, task = tsk0)
saveRDS.gz(var.mod1, "./models/t.RF.h_canopy_mean.rds")
#var.mod1 = readRDS.gz("./models/t.RF.h_canopy_mean.rds")
var.sfeats1 = getFeatSelResult(var.mod1)
str(var.sfeats1$x.bit.names)
## 108
## fine-tune xgboost
resX.xg = tuneParams(mlr::makeLearner("regr.xgboost"), task = tsk0.s, resampling = rdesc, par.set = xg.model_Params, control = ctrl)
#[Tune] Result: nrounds=50; max_depth=4; eta=0.5; subsample=1; min_child_weight=3; colsample_bytree=0.8 : mse.test.mean=0.2982494
lrn.xg = mlr::makeLearner("regr.xgboost")
lrn.xg = setHyperPars(lrn.xg, par.vals = resX.xg$x)
lrn2 = makeFeatSelWrapper(lrn.xg, resampling = inner, control = ctrlF, show.info=TRUE)
var.mod2 = mlr::train(lrn2, task = tsk0)
saveRDS.gz(var.mod2, "./models/t.XG.h_canopy_mean.rds")
var.sfeats2 = getFeatSelResult(var.mod2)
str(var.sfeats2$x.bit.names)
## new shorter formula
formulaString.canopy0 = as.formula(paste('log.h_canopy ~ ', paste(unique(c(var.sfeats1$x, var.sfeats2$x)), collapse="+")))
str(all.vars(formulaString.canopy0))
## 171
parallelMap::parallelStop()
save.image.pigz(n.cores = 64)

## Canopy EML ----
SL.library <- c("regr.ranger", "regr.xgboost", "regr.nnet")
case.weights = 1/(rm.canopy[rC.sel, "h_canopy_uncertainty"]^2)
cw.q = quantile(case.weights, c(0.2), na.rm=TRUE)
summary(is.na(case.weights))
case.weights = ifelse(is.na(case.weights), cw.q, case.weights)
parallelMap::parallelStartSocket(parallel::detectCores())
tsk <- mlr::makeRegrTask(data = rm.canopym[,all.vars(formulaString.canopy0)], target = "log.h_canopy", weights = case.weights)
lrn.rf = setHyperPars(mlr::makeLearner(SL.library[1]), par.vals = getHyperPars(var.mod1$learner))
lrn.xg = setHyperPars(mlr::makeLearner(SL.library[2]), par.vals = getHyperPars(var.mod2$learner))
lrns <- list(lrn.rf, lrn.xg, mlr::makeLearner(SL.library[3]))
init.m <- mlr::makeStackedLearner(base.learners = lrns, predict.type = "response", method = "stack.cv", super.learner = "regr.glm")
m.canopy <- mlr::train(init.m, tsk)
parallelMap::parallelStop()
saveRDS.gz(m.canopy, "./models/eml_h_canopy_mean.rds")
m.canopy$learner.model$super.model$learner.model
# Coefficients:
#   (Intercept)   regr.ranger  regr.xgboost     regr.nnet  
#     0.28170       1.06611      -0.06295      -0.06500
# 
# Degrees of Freedom: 1125394 Total (i.e. Null);  1125391 Residual
# Null Deviance:	    719400 
# Residual Deviance: 367900 	AIC: 1935000
#m.canopy = readRDS.gz("./models/eml_h_canopy_mean.rds")

## model accuracy
r.file = paste0("./models/resultsFit_h_canopy_mean.txt")
cat("Results of model fitting 'ranger', 'xgboost', 'deepnet':\n", file=r.file)
cat("\n", file=r.file, append=TRUE)
cat(paste("Variable: h_canopy\n"), file=r.file, append=TRUE)
cat(paste("R-square:", round(1-m.canopy$learner.model$super.model$learner.model$deviance / m.canopy$learner.model$super.model$learner.model$null.deviance , 3), "\n"), file=r.file, append=TRUE)
sink(file=r.file, append=TRUE, type="output")
cat("Meta learner model:", file=r.file, append=TRUE)
print(m.canopy$learner.model$super.model$learner.model)
cat("\n", file=r.file, append=TRUE)
cat("Variable importance:\n", file=r.file, append=TRUE)
xl <- mlr::getFeatureImportance(var.mod1)$res
print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:25]])))
sink()
str(m.canopy$learner.model$pred.train)
summary(m.canopy$learner.model$pred.train$regr.ranger)
summary(m.canopy$learner.model$pred.train$regr.xgboost)
summary(rm.canopym$log.h_canopy)

## Canopy predict ----
x = list.files(path="/data/tt/LandGIS/grid250m", pattern=glob2rx("^canopy.height_*.tif$"), recursive=TRUE, full.names = TRUE)
unlink(x)
pred_response.mlr(i="T38715", m=m.canopy, tile.tbl, varn="canopy.height", sd = 0, zmin = 0, zmax = 5.2, multiplier=10)
pred_response.mlr(i="T38716", m=m.canopy, tile.tbl, varn="canopy.height", sd = 0, zmin = 0, zmax = 5.2, multiplier=10)
pred_response.mlr(i="T21098", m=m.canopy, tile.tbl, varn="canopy.height", sd = 0, zmin = 0, zmax = 5.2, multiplier=10)
#expm1(c(0,5,10,15,20,25,30,35,40)/10)
expm1(21/10)

library(snowfall)
cpus = 20
sfInit(parallel=TRUE, cpus=cpus)
sfExport("m.canopy", "pred_response.mlr", "tile.tbl", "pr.dirs")
sfLibrary(ranger)
sfLibrary(xgboost)
sfLibrary(deepnet)
sfLibrary(mlr)
sfLibrary(rgdal)
sfLibrary(stats)
out <- sfClusterApplyLB(pr.dirs, function(i){ try( pred_response.mlr(i, m=m.canopy, tile.tbl, varn="canopy.height", sd = 0, zmin = 0, zmax = 5.2, multiplier=10) ) })
sfStop()
gc()

## most probable class
mosaick_ll(varn="canopy.height", i="M_sl1", out.tif="/data/Geo-wiki/predicted250m/pnv_canopy.height_icesat2.atl08.r98_m_250m_s0..0cm_2018_v0.1.tif", in.path="/data/tt/LandGIS/grid250m", tr=cellsize, te=paste(te, collapse = " "), ot="Byte", dstnodata=255, aggregate=FALSE)
mosaick_ll(varn="canopy.height", i="sd_sl1", out.tif="/data/Geo-wiki/predicted250m/pnv_canopy.height_icesat2.atl08.r98_sd_250m_s0..0cm_2018_v0.1.tif", in.path="/data/tt/LandGIS/grid250m", tr=cellsize, te=paste(te, collapse = " "), ot="Byte", dstnodata=255, aggregate=FALSE)
