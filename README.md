# Predicting fluvial sedimentary geochemistry

This repository contains supporting code and data for the publication "River sediment geochemistry as a conservative mixture of source regions: Observations and predictions from the Cairngorms, UK" (under review)

This repository contains 4 files.

## Requisites 

Code was written in `python v3.7.3` using the module [`LandLab`](https://pypi.org/project/landlab/) `v1.9.0`. 

### Code 

`predict_sed_geochem.py` is a minimal working example for predicting the geochemical composition of higher-order river sediments from topographic data and a geochemical map of the sediment source region. This example models the magnesium concentration of sediment in rivers draining the Cairngorms, UK. Input files required to run this example are contained in `input_dir/`.`topo_CG.nc` is a DEM of the studied regions in the Cairngorms UK, created by gridding the `SRTM1s` DEM to 200 x 200 m squares. `CG_log_Mg.nc` is a map of the (log transformed) Magnesium concentration in first order stream sediments created by interpolating the [G-BASE](https://www.bgs.ac.uk/gbase/home.html) geochemical survey.   

### Data

`higher_order_geochemistry.csv` contains the geochemical composition of higher-order river sediments sampled at 67 points across the Cairngorms and the surrounding region.  

