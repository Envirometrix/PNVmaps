# Future projections of Potential Natural Vegetation across different climatic scenarios based on Machine Learning

PNV predictions of the general IUCN classes and BIOME 6000 classes at 1 km spatial resolution are available for **[download](https://doi.org/10.5281/zenodo.7520813)**.

![Biomes map at 1km for RCP 2.6 epoch 2040-2060](img/001_iucn_biomes.png "Potential distribution of terrestrial biomes (Potential Natural Vegetation) at 1 km spatial resolution.")

Improvements in the future projections of the PNV biomes maps at 1km:

- Only 72 covariates were used (temperature, precipitation and topographical covariates), but model accuracy was doubled,
- Predictions are based on the Ensemble Machine Learning (`"classif.ranger", "classif.glmnet", "classif.xgboost"`) as implemented in the [mlr package](https://mlr.mlr-org.com/),
- Predictions are provided in the original BIOME 6000 classification system (20 classes) and the **[IUCN Global Ecosystem Typology](https://global-ecosystems.org/page/typology)** classification system
- Model errors are provided per class (derived as the weighted standard deviation between multiple models) for single class maps,
- Model errors are provided using the _margin of victory_ [(Calderón-Loor et al., 2021)](https://doi.org/10.1016/j.rse.2020.112148) for hard classes maps,
- Model fine-tuning and accuracy assessment is based on repeated 5-fold spatial cross-validation,
- Predictions are provided for current and future epochs (2040-2060 and 2060-2080) under three climatic scenarios (RCP 2.6, RCP 4.5 and RCP 8.5).


# Global Maps of Potential Natural Vegetation based on Machine Learning

PNV predictions of the general land cover classes at 250 m spatial resolution is available for **[download](https://doi.org/10.5281/zenodo.3631253)**.

![GLC map at 250m](img/001_pnv_predictions_glc100.png "Potential distribution of land cover classes (Potential Natural Vegetation) at 250 m spatial resolution.")


Update of the predictions at 250 m spatial resolution is available for **[download](https://doi.org/10.5281/zenodo.3526619)**.

![Biomes map at 250m](img/001_pnv_biome.type_biome00k_c_250m_s0..0cm_2000..2017_v0.2.png "Potential distribution of biomes (Potential Natural Vegetation) at 250 m spatial resolution.")

Improvements in the v0.2 of the PNV biomes map at 250 m:

- 40% of new covariates have been added (see [the variable importance list](R_code/Biome_randomForest_v02.txt))
- Predictions are based on the Ensemble Machine Learning (`"classif.ranger", "classif.glmnet", "classif.xgboost", "classif.nnTrain"`) as implemented in the [mlr package](https://mlr.mlr-org.com/),
- Model errors are provided per class (derived as the weighted standard deviation between multiple models),
- Model fine-tuning, feature selection and accuracy assessment is based on repeated cross-validation,

*Summary*: This repository contains R code and some outputs of spatial predictions related with the production of [Global Maps of Potential Natural Vegetation](https://www.arcgis.com/apps/MapJournal/index.html?appid=1856322400844a7cab348bccfa4bee76). Three case studies were considered: (1) global distribution of biomes based on the BIOME 6000 data set (8057 modern pollen-based site reconstructions), (2) distribution of forest tree species in Europe based on detailed occurrence records (1,546,435 ground observations), and (3) global monthly Fraction of Absorbed Photosynthetically Active Radiation (FAPAR) values (30,301 randomly-sampled points).

![alt text](https://github.com/envirometrix/PNVmaps/blob/master/img/Fig_global_biomes_map.png "Output predictions for global biomes.")

# Step-by-step tutorial

[This tutorial](https://github.com/Envirometrix/PNVmaps/tree/master/tutorial) explains how to fit models and produce predictions for smaller area in Europe. To run this tutorial you might need to install and customize some [R / OS GIS software](https://envirometrix.github.io/PredictiveSoilMapping/software.html).

*Please cite as:*

Future projections:
* Bonannella C, Hengl T, Parente L, de Bruin S. 2023. **[Biomes of the world under climate change scenarios: increasing aridity and higher temperatures lead to significant shifts in natural vegetation](https://doi.org/10.7717/peerj.15593)**. PeerJ 11:e15593 https://doi.org/10.7717/peerj.15593

Original PNV maps:
* Hengl T, Walsh MG, Sanderman J, Wheeler I, Harrison SP, Prentice IC. 2018. **[Global mapping of potential natural vegetation: an assessment of machine learning algorithms for estimating land potential](https://doi.org/10.7717/peerj.5457)**. PeerJ 6:e5457 https://doi.org/10.7717/peerj.5457

# Download Maps

The 250m resolution predictions of biomes are available for download from https://doi.org/10.5281/zenodo.3526619

The 1km resolution maps are available for download under the [Open Database License (ODbl) v1.0](https://opendatacommons.org/licenses/odbl/) and can be downloaded from http://dx.doi.org/10.7910/DVN/QQHCIK without restrictions.

# Disclaimer

These are premilimary maps of the Global Potential Natural Vegetation. Errors and artifacts are still possible. Training data sets BIOME 6000 and EU Forest are constantly being updated and could still contain erroneously geolocated points. Predictions of FAPAR are based on randomly simulated points and not on ground observations of FAPAR and classification of sites. Predictions of EU forest tree species are presented for experimental purposes only. To report an issue or artifact in maps, please use https://github.com/envirometrix/PNVmaps/issues.

