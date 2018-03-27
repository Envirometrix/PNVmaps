## PNV mapping functions
## tom.hengl@gmail.com

## prediction locations / data ----
new_data <- function(i, tif.sel, tile.tbl, tif.lc="/data/LandGIS/layers1km/lcv_landmask_esacci.lc.l4_c_1km_s0..0cm_2000..2015_v1.0.tif", tif.fapar="/data/LandGIS/upscaled1km/veg_fapar_proba.v.annual_d_1km_s0..0cm_2014..2017_v1.0.tif", out.dir="/data/PNV/tiled"){
  i.n = which(tile.tbl$ID == strsplit(i, "T")[[1]][2])
  out.rds <- paste0(out.dir, "/T", tile.tbl[i.n,"ID"], "/T", tile.tbl[i.n,"ID"], ".rds")
  if(!file.exists(out.rds)){
    m = readGDAL(fname=tif.lc, offset=unlist(tile.tbl[i.n,c("offset.y","offset.x")]), region.dim=unlist(tile.tbl[i.n,c("region.dim.y","region.dim.x")]), output.dim=unlist(tile.tbl[i.n,c("region.dim.y","region.dim.x")]), silent = TRUE)
    m$band2 = readGDAL(fname=tif.fapar, offset=unlist(tile.tbl[i.n,c("offset.y","offset.x")]), region.dim=unlist(tile.tbl[i.n,c("region.dim.y","region.dim.x")]), output.dim=unlist(tile.tbl[i.n,c("region.dim.y","region.dim.x")]), silent = TRUE)$band1
    m = as(m, "SpatialPixelsDataFrame")
    ## Filter out water bodies, permanent ice and shifting sands:
    sel.p = (m$band1==1|m$band1==3) & (!m$band2==0 | is.na(m$band2)) ## FAPAR annual not available for north lats
    if(sum(sel.p)>0){
      m = m[sel.p,]
      for(j in 1:length(tif.sel)){
        m@data[,j+2] = readGDAL(fname=tif.sel[j], offset=unlist(tile.tbl[i.n,c("offset.y","offset.x")]), region.dim=unlist(tile.tbl[i.n,c("region.dim.y","region.dim.x")]), output.dim=unlist(tile.tbl[i.n,c("region.dim.y","region.dim.x")]), silent=TRUE)$band1[m@grid.index]
      }
      names(m) = c("LandCover", "FAPAR.m", basename(tif.sel))
      ## Systematic fixes
      sn.sel = grep(pattern="snow.prob", names(m))
      for(p in sn.sel){ m@data[,p] = ifelse(m@data[,p]>100, NA, m@data[,p]) }
      ## Missing values in latitudes >61
      if(any(m@coords[,2]>61.4)){
        m$clm_snow.prob_esacci.jan_p_1km_s0..0cm_2000..2016_v1.0.tif = ifelse(m@coords[,2]>61.4 & is.na(m$clm_snow.prob_esacci.jan_p_1km_s0..0cm_2000..2016_v1.0.tif), 100, m$clm_snow.prob_esacci.jan_p_1km_s0..0cm_2000..2016_v1.0.tif)
        m$clm_snow.prob_esacci.dec_p_1km_s0..0cm_2000..2016_v1.0.tif = ifelse(m@coords[,2]>61.4 & is.na(m$clm_snow.prob_esacci.dec_p_1km_s0..0cm_2000..2016_v1.0.tif), 100, m$clm_snow.prob_esacci.dec_p_1km_s0..0cm_2000..2016_v1.0.tif)
        m$clm_snow.prob_esacci.nov_p_1km_s0..0cm_2000..2016_v1.0.tif = ifelse(m@coords[,2]>61.4 & is.na(m$clm_snow.prob_esacci.nov_p_1km_s0..0cm_2000..2016_v1.0.tif), 100, m$clm_snow.prob_esacci.nov_p_1km_s0..0cm_2000..2016_v1.0.tif)
      }
      ## Filter northern latitudes MODFC:
      m$clm_cloud.fraction_earthenv.modis.11_m_1km_s0..0cm_2000..2015_v1.0.tif = ifelse(is.na(m$clm_cloud.fraction_earthenv.modis.11_m_1km_s0..0cm_2000..2015_v1.0.tif), rowMeans(m@data[,c("clm_cloud.fraction_earthenv.modis.9_m_1km_s0..0cm_2000..2015_v1.0.tif", "clm_cloud.fraction_earthenv.modis.10_m_1km_s0..0cm_2000..2015_v1.0.tif", "clm_cloud.fraction_earthenv.modis.12_m_1km_s0..0cm_2000..2015_v1.0.tif", "clm_cloud.fraction_earthenv.modis.1_m_1km_s0..0cm_2000..2015_v1.0.tif")], na.rm=TRUE), m$clm_cloud.fraction_earthenv.modis.11_m_1km_s0..0cm_2000..2015_v1.0.tif)
      m$clm_cloud.fraction_earthenv.modis.12_m_1km_s0..0cm_2000..2015_v1.0.tif = ifelse(is.na(m$clm_cloud.fraction_earthenv.modis.12_m_1km_s0..0cm_2000..2015_v1.0.tif), rowMeans(m@data[,c("clm_cloud.fraction_earthenv.modis.10_m_1km_s0..0cm_2000..2015_v1.0.tif", "clm_cloud.fraction_earthenv.modis.11_m_1km_s0..0cm_2000..2015_v1.0.tif", "clm_cloud.fraction_earthenv.modis.1_m_1km_s0..0cm_2000..2015_v1.0.tif", "clm_cloud.fraction_earthenv.modis.2_m_1km_s0..0cm_2000..2015_v1.0.tif")], na.rm=TRUE), m$clm_cloud.fraction_earthenv.modis.12_m_1km_s0..0cm_2000..2015_v1.0.tif)
      m$clm_cloud.fraction_earthenv.modis.1_m_1km_s0..0cm_2000..2015_v1.0.tif = ifelse(is.na(m$clm_cloud.fraction_earthenv.modis.1_m_1km_s0..0cm_2000..2015_v1.0.tif), rowMeans(m@data[,c("clm_cloud.fraction_earthenv.modis.11_m_1km_s0..0cm_2000..2015_v1.0.tif", "clm_cloud.fraction_earthenv.modis.12_m_1km_s0..0cm_2000..2015_v1.0.tif", "clm_cloud.fraction_earthenv.modis.2_m_1km_s0..0cm_2000..2015_v1.0.tif", "clm_cloud.fraction_earthenv.modis.3_m_1km_s0..0cm_2000..2015_v1.0.tif")], na.rm=TRUE), m$clm_cloud.fraction_earthenv.modis.1_m_1km_s0..0cm_2000..2015_v1.0.tif)
      ## Filter Greenland
      m$dtm_lithology_usgs.ecotapestry.07_p_1km_s0..0cm_2014_v1.0.tif = ifelse(is.na(m$dtm_lithology_usgs.ecotapestry.07_p_1km_s0..0cm_2014_v1.0.tif) & m$lcv_admin0_fao.gaul_c_1km_s0..0cm_2015_v1.0.tif==99, 100, m$dtm_lithology_usgs.ecotapestry.07_p_1km_s0..0cm_2014_v1.0.tif)
      ## Global Surface Water
      m$lcv_surface.water.occ_gsw.jrc_p_1km_b0..200cm_1984..2016_v1.0.tif = ifelse(m$lcv_surface.water.occ_gsw.jrc_p_1km_b0..200cm_1984..2016_v1.0.tif>100 | is.na(m$lcv_surface.water.occ_gsw.jrc_p_1km_b0..200cm_1984..2016_v1.0.tif), 0, m$lcv_surface.water.occ_gsw.jrc_p_1km_b0..200cm_1984..2016_v1.0.tif)
      ## GIEMS
      m$dtm_inundation.extent_giems.d15_m_1km_s0..0cm_2015_v1.0.tif = ifelse(is.na(m$dtm_inundation.extent_giems.d15_m_1km_s0..0cm_2015_v1.0.tif), 0, m$dtm_inundation.extent_giems.d15_m_1km_s0..0cm_2015_v1.0.tif)
      ## Indicators
      us.sel = grep(pattern="usgs.ecotapestry", names(m))
      for(q in us.sel){ m@data[,q] = ifelse(is.na(m@data[,q]), 0, m@data[,q]) }
      ## Fill-in the remaining missing values (can be very tricky)
      sel.mis = sapply(m@data[,-unlist(sapply(c("usgs.ecotapestry", "admin0"), function(i){grep(i,names(m))}))], function(x){sum(is.na(x))>0})
      if(sum(sel.mis)>0){
        x = which(sel.mis)
        for(k in 1:length(x)){
          if(!is.factor(m@data[,attr(x, "names")[k]])){
            if(length(grep(pattern="occurrence", attr(x, "names")[k]))>0 | length(grep(pattern="snow.prob", attr(x, "names")[k]))>0 | length(grep(pattern="usgs.ecotapestry", attr(x, "names")[k]))>0 ){ 
              repn = rep(0, nrow(m)) 
            } else {
              r = raster::raster(m[attr(x, "names")[k]])
              ## 1 using proximity filter:
              rf = raster::focal(r, w=matrix(1,15,15), fun=mean, pad=TRUE, na.rm=TRUE, NAonly=TRUE)
              repn = as(rf, "SpatialGridDataFrame")@data[m@grid.index,1]
              ## 2 using dominant value:
              repn = ifelse(is.na(repn), quantile(repn, probs=.5, na.rm=TRUE), repn)
            }
            m@data[,attr(x, "names")[k]] = ifelse(is.na(m@data[,attr(x, "names")[k]]), repn, m@data[,attr(x, "names")[k]])
          }
        }
      }
      saveRDS(m, out.rds)
    }
  }  
}


pred_probs = function(i, gm, tile.tbl, col.legend, varn, out.dir="/data/PNV/tiled", with.se=FALSE){
  i.n = which(tile.tbl$ID == strsplit(i, "T")[[1]][2])
  out.rds <- paste0(out.dir, "/T", tile.tbl[i.n,"ID"], "/T", tile.tbl[i.n,"ID"], ".rds")
  out.c <- paste0(out.dir, "/", i, "/", varn, "_C_", i, ".tif")
  if(file.exists(out.rds) & !file.exists(out.c)){
    m = readRDS(out.rds)
    m = m[complete.cases(m@data[,gm$forest$independent.variable.names]),]
    if(with.se==TRUE){
      pred = predict(gm, m@data, type="se")
      tax = gm$forest$levels[gm$forest$class.values]
    } else {
      pred = predict(gm, m@data)
      tax = attr(pred$predictions, "dimnames")[[2]]
    }
    rs <- rowSums(pred$predictions, na.rm=TRUE)
    ## Write GeoTiffs:
    if(sum(rs,na.rm=TRUE)>0&length(rs)>0){
      ## predictions
      m@data <- data.frame(pred$predictions)
      x = m[1]
      for(j in 1:ncol(m)){
        out <- paste0(out.dir, "/", i, "/", varn, "_M_", tax[j], "_", i, ".tif")
        x@data[,1] <- round(m@data[,j]*100)
        writeGDAL(x[1], out, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
      }
      ## most probable class:
      col.tbl <- plyr::join(data.frame(Group=tax, int=1:length(tax)), col.legend, type="left")
      ## match most probable class
      m$cl <- col.tbl[match(apply(m@data,1,which.max), col.tbl$int),"Number"]  
      writeGDAL(m["cl"], out.c, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
      ## error maps
      if(with.se==TRUE){
        m@data <- data.frame(pred$se)
        x = m[1]
        for(j in 1:ncol(m)){
          out.se <- paste0(out.dir, "/", i, "/", varn, "_se_", tax[j], "_", i, ".tif")
          x@data[,1] <- round(m@data[,j]*100)
          writeGDAL(x[1], out.se, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
        }
      }
    }
    #gc()
  }
}

pred_FAPAR = function(i, gm, tile.tbl, varn="FAPAR", level="M", months=1:12, m.lst=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"), out.dir="/data/PNV/tiled", quantiles = c((1-.682)/2, 1-(1-.682)/2), with.se=TRUE, vars.t = c("clm_water.vapor_nasa.eo.","clm_precipitation_imerge.","u.99_1km_s0..0cm_1979..2013", "m_1km_s0..0cm_1979..2013","l.01_1km_s0..0cm_1979..2013","clm_cloud.fraction_earthenv.modis.","clm_snow.prob_esacci."), ttn.lst = c("water.vapor","precipitation","temp_max","temp_mean","temp_min","could.fraction","snow.prob")){
  i.n = which(tile.tbl$ID == strsplit(i, "T")[[1]][2])
  out.rds <- paste0(out.dir, "/T", tile.tbl[i.n,"ID"], "/T", tile.tbl[i.n,"ID"], ".rds")
  out.m <- paste0(out.dir, "/", i, "/", varn, "_", m.lst, "_", level, "_", i, ".tif")
  if(file.exists(out.rds) & any(!file.exists(out.m))){
    m = readRDS(out.rds)
    if(nrow(m@data)>0){
      for(j in months){
        m@data[,"cMonth"] = cos((j/12)*(2*pi))
        ## rename some columns:
        for(p in 1:length(vars.t)){
          if(vars.t[p]=="clm_cloud.fraction_earthenv.modis."){
            fname = grep(paste0(vars.t[p], j, "_"), names(m), fixed = TRUE)
          } else {
            if(length(grep("temp_", ttn.lst[p]))>0){
              fname = grep(paste0(tolower(m.lst[j]), "_", vars.t[p], "_"), names(m), fixed = TRUE)
            } else {
              fname = grep(paste0(vars.t[p], tolower(m.lst[j]), "_"), names(m), fixed = TRUE)
            }
          }
          m@data[,ttn.lst[p]] = m@data[,fname]
        }
        m = m[complete.cases(m@data[,gm$forest$independent.variable.names[-1]]),]
        if(with.se==TRUE){
          ## separately predictions of the mean value and lower/upper range
          m$v <- predict(gm, m@data)$predictions
          pred = predict(gm, m@data, type="quantiles", quantiles=quantiles)
          m$se = (pred$predictions[,2]-pred$predictions[,1])/2
          if(level=="M"){ 
            se.t = paste0("se")  
          } else {
            se.t = paste0("se", level)
          }
          writeGDAL(m["se"], paste0(out.dir, "/", i, "/", varn, "_", m.lst[j], "_", se.t, "_", i, ".tif"), type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
        } else {
          pred = predict(gm, m@data)
          m$v <- pred$predictions
        }
        writeGDAL(m["v"], out.m[j], type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
      }
    }
  }
}

## Fix biomes regression matrix:
fix_biome <- function(rm.biome, Lat="Latitude"){
  sn.sel = grep(pattern="snow.prob", names(rm.biome))
  ## problem with snow probability maps -- values above 100%
  for(p in sn.sel){ rm.biome[,p] = ifelse(rm.biome[,p]>100, NA, rm.biome[,p]) }
  ## Missing values in snow.prob latitudes >61
  rm.biome$clm_snow.prob_esacci.jan_p_1km_s0..0cm_2000..2016_v1.0.tif = ifelse(rm.biome[,Lat]>61.4 & is.na(rm.biome$clm_snow.prob_esacci.jan_p_1km_s0..0cm_2000..2016_v1.0.tif), 100, rm.biome$clm_snow.prob_esacci.jan_p_1km_s0..0cm_2000..2016_v1.0.tif)
  rm.biome$clm_snow.prob_esacci.dec_p_1km_s0..0cm_2000..2016_v1.0.tif = ifelse(rm.biome[,Lat]>61.4 & is.na(rm.biome$clm_snow.prob_esacci.dec_p_1km_s0..0cm_2000..2016_v1.0.tif), 100, rm.biome$clm_snow.prob_esacci.dec_p_1km_s0..0cm_2000..2016_v1.0.tif)
  rm.biome$clm_snow.prob_esacci.nov_p_1km_s0..0cm_2000..2016_v1.0.tif = ifelse(rm.biome[,Lat]>61.4 & is.na(rm.biome$clm_snow.prob_esacci.nov_p_1km_s0..0cm_2000..2016_v1.0.tif), 100, rm.biome$clm_snow.prob_esacci.nov_p_1km_s0..0cm_2000..2016_v1.0.tif)
  ## Missing values for norther latitudes for MODFC:
  rm.biome$clm_cloud.fraction_earthenv.modis.11_m_1km_s0..0cm_2000..2015_v1.0.tif = ifelse(is.na(rm.biome$clm_cloud.fraction_earthenv.modis.11_m_1km_s0..0cm_2000..2015_v1.0.tif), rowMeans(rm.biome[,c("clm_cloud.fraction_earthenv.modis.9_m_1km_s0..0cm_2000..2015_v1.0.tif", "clm_cloud.fraction_earthenv.modis.10_m_1km_s0..0cm_2000..2015_v1.0.tif", "clm_cloud.fraction_earthenv.modis.12_m_1km_s0..0cm_2000..2015_v1.0.tif", "clm_cloud.fraction_earthenv.modis.1_m_1km_s0..0cm_2000..2015_v1.0.tif")], na.rm=TRUE), rm.biome$clm_cloud.fraction_earthenv.modis.11_m_1km_s0..0cm_2000..2015_v1.0.tif)
  rm.biome$clm_cloud.fraction_earthenv.modis.12_m_1km_s0..0cm_2000..2015_v1.0.tif = ifelse(is.na(rm.biome$clm_cloud.fraction_earthenv.modis.12_m_1km_s0..0cm_2000..2015_v1.0.tif), rowMeans(rm.biome[,c("clm_cloud.fraction_earthenv.modis.10_m_1km_s0..0cm_2000..2015_v1.0.tif", "clm_cloud.fraction_earthenv.modis.11_m_1km_s0..0cm_2000..2015_v1.0.tif", "clm_cloud.fraction_earthenv.modis.1_m_1km_s0..0cm_2000..2015_v1.0.tif", "clm_cloud.fraction_earthenv.modis.2_m_1km_s0..0cm_2000..2015_v1.0.tif")], na.rm=TRUE), rm.biome$clm_cloud.fraction_earthenv.modis.12_m_1km_s0..0cm_2000..2015_v1.0.tif)
  rm.biome$clm_cloud.fraction_earthenv.modis.1_m_1km_s0..0cm_2000..2015_v1.0.tif = ifelse(is.na(rm.biome$clm_cloud.fraction_earthenv.modis.1_m_1km_s0..0cm_2000..2015_v1.0.tif), rowMeans(rm.biome[,c("clm_cloud.fraction_earthenv.modis.11_m_1km_s0..0cm_2000..2015_v1.0.tif", "clm_cloud.fraction_earthenv.modis.12_m_1km_s0..0cm_2000..2015_v1.0.tif", "clm_cloud.fraction_earthenv.modis.2_m_1km_s0..0cm_2000..2015_v1.0.tif", "clm_cloud.fraction_earthenv.modis.3_m_1km_s0..0cm_2000..2015_v1.0.tif")], na.rm=TRUE), rm.biome$clm_cloud.fraction_earthenv.modis.1_m_1km_s0..0cm_2000..2015_v1.0.tif)
  ## Missing values for Greenland:
  rm.biome$dtm_lithology_usgs.ecotapestry.07_p_1km_s0..0cm_2014_v1.0.tif = ifelse(is.na(rm.biome$dtm_lithology_usgs.ecotapestry.07_p_1km_s0..0cm_2014_v1.0.tif) & rm.biome$lcv_admin0_fao.gaul_c_1km_s0..0cm_2015_v1.0.tif==99, 100, rm.biome$dtm_lithology_usgs.ecotapestry.07_p_1km_s0..0cm_2014_v1.0.tif)
  ## GIEMS:
  rm.biome$dtm_inundation.extent_giems.d15_m_1km_s0..0cm_2015_v1.0.tif = ifelse(is.na(rm.biome$dtm_inundation.extent_giems.d15_m_1km_s0..0cm_2015_v1.0.tif), 0, rm.biome$dtm_inundation.extent_giems.d15_m_1km_s0..0cm_2015_v1.0.tif)
  ## All other snowfall should be 0
  sn.sel = grep(pattern="snow.prob", names(rm.biome))
  for(p in sn.sel){ rm.biome[,p] = ifelse(is.na(rm.biome[,p]), 0, rm.biome[,p]) }
  us.sel = grep(pattern="usgs.ecotapestry", names(rm.biome))
  for(q in us.sel){ rm.biome[,q] = ifelse(is.na(rm.biome[,q]), 0, rm.biome[,q]) }
  rm.biome$lcv_surface.water.occ_gsw.jrc_p_1km_b0..200cm_1984..2016_v1.0.tif = ifelse(rm.biome$lcv_surface.water.occ_gsw.jrc_p_1km_b0..200cm_1984..2016_v1.0.tif>100 | is.na(rm.biome$lcv_surface.water.occ_gsw.jrc_p_1km_b0..200cm_1984..2016_v1.0.tif), 0, rm.biome$lcv_surface.water.occ_gsw.jrc_p_1km_b0..200cm_1984..2016_v1.0.tif)
  s.rm = sapply(rm.biome, function(i){sum(is.na(i))})
  if(sum(s.rm>200)>0){
    message("Missing values detected:")
    print(s.rm[which(s.rm>200)])
  }
  return(rm.biome)
}


## Mask selected admin classes
mask_admin <- function(tif, sel.admin, admin.tif="/data/PNV/R_code/EUforest_ADMIN1.tif"){
  m = readGDAL(admin.tif)
  m$band1 = as.factor(m$band1)
  m$p = readGDAL(tif)$band1
  m$f = ifelse(!m$band1 %in% sel.admin, m$p, NA)
  writeGDAL(m["f"], tif, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
  gc()
}
