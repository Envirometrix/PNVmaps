# Global Maps of Potential Natural Vegetation at 1 km resolution

Update of the predictions at 250 m spatial resolution is available for **[download](https://doi.org/10.5281/zenodo.3526619)**.

![Biomes map at 250m](img/001_pnv_biome.type_biome00k_c_250m_s0..0cm_2000..2017_v0.2.png "Potential distribution of biomes (Potential Natural Vegetation) at 250 m spatial resolution.")

Improvements in the v0.2 of the PNV biomes map at 250 m:

- 40% of new covariates have been added (see [the variable importance list](R_code/Biome_randomForest_v02.txt))
- Predictions are based on the Ensemble Machine Learning (`"classif.ranger", "classif.glmnet", "classif.xgboost", "classif.nnTrain"`) as implemented in the [mlr package](https://mlr.mlr-org.com/),
- Model errors are provided per class (derived as the weighted standard deviation between multiple models),
- Model fine-tuning, feature selection and accuracy assessment is based on repeated cross-validation,

Summary: This repository contains R code and some outputs of spatial predictions related with the production of [Global Maps of Potential Natural Vegetation](https://www.arcgis.com/apps/MapJournal/index.html?appid=1856322400844a7cab348bccfa4bee76). Three case studies were considered: (1) global distribution of biomes based on the BIOME 6000 data set (8057 modern pollen-based site reconstructions), (2) distribution of forest tree species in Europe based on detailed occurrence records (1,546,435 ground observations), and (3) global monthly Fraction of Absorbed Photosynthetically Active Radiation (FAPAR) values (30,301 randomly-sampled points).

![alt text](https://github.com/envirometrix/PNVmaps/blob/master/img/Fig_global_biomes_map.png "Output predictions for global biomes.")

# Step-by-step tutorial

[This tutorial](https://github.com/Envirometrix/PNVmaps/tree/master/tutorial) explains how to fit models and produce predictions for smaller area in Europe. To run this tutorial you might need to install and customize some [R / OS GIS software](https://envirometrix.github.io/PredictiveSoilMapping/software.html).

*Please cite as:*

* Hengl T, Walsh MG, Sanderman J, Wheeler I, Harrison SP, Prentice IC. 2018. **[Global mapping of potential natural vegetation: an assessment of machine learning algorithms for estimating land potential](https://doi.org/10.7717/peerj.5457)**. PeerJ 6:e5457 https://doi.org/10.7717/peerj.5457

# Download Maps

The 250m resolution predictions of biomes are available for download from https://doi.org/10.5281/zenodo.3526619

The 1km resolution maps are available for download under the [Open Database License (ODbl) v1.0](https://opendatacommons.org/licenses/odbl/) and can be downloaded from http://dx.doi.org/10.7910/DVN/QQHCIK without restrictions.

# Disclaimer

These are premilimary maps of the Global Potential Natural Vegetation. Errors and artifacts are still possible. Training data sets BIOME 6000 and EU Forest are constantly being updated and could still contain erroneously geolocated points. Predictions of FAPAR are based on randomly simulated points and not on ground observations of FAPAR and classification of sites. Predictions of EU forest tree species are presented for experimental purposes only. To report an issue or artifact in maps, please use https://github.com/envirometrix/PNVmaps/issues.

