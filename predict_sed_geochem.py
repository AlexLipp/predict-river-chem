#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import netCDF4
import datetime
import numpy as np
from landlab import RasterModelGrid
from landlab import imshow_grid
from landlab.components import FlowAccumulator, FastscapeEroder, SinkFillerBarnes
from landlab import FIXED_VALUE_BOUNDARY
from landlab.components.flow_accum import find_drainage_area_and_discharge


"""
Created on Wed May 20 13:07:05 2020


This script is a minimum working example for predicting the composition of higher order
fluvial sediments from topographic data and geochemical maps. 

Inputs: 
    1. DEM of Cairngorms as a netCDF file (input_dir/topo_CG.nc)
    2. Geochemical map of source region Mg concentration after log10 transform. 
    This was created by interpolating G-BASE geochemical survey (input_dir/CG_log_Mg.nc)

Outputs: 
    This script produces no outputs, but visualises the results as maps.

@author: Alex Lipp; a.lipp18@imperial.ac.uk
"""

print("#### Modelling landscape evolution and composition of sediment flux ####")
print("Starting")
print(datetime.datetime.now().time())


#### Set up model grid ####

print("Loading in topo")
print(datetime.datetime.now().time())

zr_nc=netCDF4.Dataset('input_dir/topo_CG.nc')
zr_ma = zr_nc['z'][:,:]
mg = RasterModelGrid(zr_ma.shape)
zr = mg.add_zeros('node', 'topographic__elevation')
zr += np.reshape(zr_ma.data.astype(float),zr.shape)
dx = 200
mg._dx,mg._dy = (dx,dx)

# Set up the boundary conditions on the square grid
mg.status_at_node[mg.nodes_at_left_edge] = FIXED_VALUE_BOUNDARY
mg.status_at_node[mg.nodes_at_top_edge] = FIXED_VALUE_BOUNDARY
mg.status_at_node[mg.nodes_at_bottom_edge] = FIXED_VALUE_BOUNDARY
mg.status_at_node[mg.nodes_at_right_edge] = FIXED_VALUE_BOUNDARY


# This natural data has some pits/depressions in which we want to fill
# using the sink filler algorithm

print("Filling pits")
print(datetime.datetime.now().time())
sfb = SinkFillerBarnes(mg, method='D8',fill_flat=False)
sfb.run_one_step()

# instantiate the flow routing:
frr = FlowAccumulator(
    mg,
    'topographic__elevation',
    flow_director = 'FlowDirectorD8')

# instantiate the stream power law
k=3.62
m=0.35
n=1
spr = FastscapeEroder(mg, K_sp=k, m_sp=m, n_sp=n)
# NB k is only needed for calculating incision, absolute comp values are 
# independent.

# Set up the incision field
incise = mg.add_zeros('node', 'incision')

# set model variables. dt is set by courant condition.

dt=0.1# Timestep in Ma
cell_area = mg.dx * mg.dy

######################################
### run erosion model for one step ###
######################################

print("Running flow router")
print(datetime.datetime.now().time())
frr.run_one_step()  # flow routing
zrold = zr[mg.nodes] # pre-incision elevations

print("Running fastscaper")
print(datetime.datetime.now().time())
spr.run_one_step(dt)    # erosion: stream power
zrnew = zr[mg.nodes]    # post-incision elevations
incise = zrold - zrnew  # incision per cell
qs = incise*cell_area # Volume sediment produced per cell
qsflat = qs.ravel()     # flatten qs for flow routing calculation

# extract cumulative flux (q) as function of flow length.
a,q = find_drainage_area_and_discharge(mg.at_node['flow__upstream_node_order'], 
                                       mg.at_node['flow__receiver_node'],
                                       runoff = qsflat) # a is number of nodes

area = mg.at_node['drainage_area']
mg.add_field('node','flux',q,noclobber=False)
area_threshold = 25 #float(sys.argv[1]) #km2
is_drainage = area>(area_threshold*1000000) #km2 to m2


q_rate = q/dt
incise_rate = (incise/dt)*0.001 #mm/yr
# Update fields
print("Updating fields")
print(datetime.datetime.now().time())
mg.add_field('node', 'elevation', zrold, noclobber=False)
# Outputs present day elevation
mg.add_field('node', 'incision', incise_rate, noclobber=False)
# Channels
mg.add_field('node','channels',is_drainage,noclobber=False)
# Incision rate
mg.add_field('node', 'sedflux', q_rate, noclobber=False)
# Sed flux rate

# Visualise results 
print("Topography")
imshow_grid(mg,'elevation',output= True, plot_name = "Topography" ,
            cmap='gist_earth',limits=(0,1000),var_name="Elevation, m ")
print("Drainage")
imshow_grid(mg,'channels',output= True ,plot_name = "Drainage", 
            cmap='gist_earth',limits=(0,1))
print("Incision Rate")
imshow_grid(mg,'incision',output= True,plot_name = "Incision Rate",
            cmap='Reds', var_name = "Incision Rate, mm/yr")
print("Cumulative Sediment Flux")
imshow_grid(mg,'sedflux',output= True,plot_name = "Cumulative sediment flux",
            cmap='Blues',limits=(0,1e11), var_name = "m3 / Ma ")

#############################################
#### Predicting geochemistry - Magnesium #### 
#############################################

# Now we use results above to predict Mg concentration in river sediments 

print("Loading in input Mg composition")
print(datetime.datetime.now().time())

## For the log10 interpolated values
comp_nc=netCDF4.Dataset('input_dir/CG_log_Mg.nc')
comp_ma = comp_nc['z'][:,:]
comp_log = mg.add_zeros('node','log_bdrck__Mg')
comp_log += np.reshape((comp_ma.data.astype(float)),comp_log.shape)
comp = mg.add_zeros('node','bdrck__Mg')
comp += np.reshape(((10**comp_ma).data.astype(float)),comp.shape) # Convert to raw values from log10


print('Calculating sediment Mg composition')
print(datetime.datetime.now().time())
a, q_component = find_drainage_area_and_discharge(mg.at_node['flow__upstream_node_order'], 
                                                  mg.at_node['flow__receiver_node'],
                                                  runoff = qsflat * comp)  
# a is number of nodes
# q_component is the total flux of the compositional component (e.g. kg Mg eroded)
    
q_comp_norm = q_component/q
q_comp_norm[q==0] = comp[q==0]
# To better visualise composition in channels we now multiply it by the grid 
# `channels'. 
q_comp_norm_channel = np.log10(q_comp_norm) * is_drainage

# q_comp_norm is the composition of the sedflux at each point component (e.g. wt% Mg)
# Where cumulative sedflux is 0, q_comp is set to 'bedrock' composition

# Attach omposition of sediment as grid layer
mg.add_field('node','sed__Mg',q_comp_norm,noclobber=False)
mg.add_field('node','log_sed__Mg',np.log10(q_comp_norm),noclobber=False)
mg.add_field('node','log_sed__Mg_channel',q_comp_norm_channel,noclobber=False)


print("Sediment Mg concentration (stream power law)")
imshow_grid(mg,'log_sed__Mg_channel',output= True,plot_name ="Sediment Mg concentration (SPL)",
            cmap='viridis',limits=(3.2,4.4),var_name="Mg conc. / log10")

#### Calculating sediment composition for homogeneous incision ####

# q_component_homo_incis is the total flux of the compositional component (e.g. kg of Mg eroded) if incision is homogeneous
a, q_component_homo_incis = find_drainage_area_and_discharge(mg.at_node['flow__upstream_node_order'], 
                                                             mg.at_node['flow__receiver_node'],
                                                             runoff = comp)  # composition but homogeneous erosion
a, q_homo_incis = find_drainage_area_and_discharge(mg.at_node['flow__upstream_node_order'], 
                                                   mg.at_node['flow__receiver_node']) # a is number of nodes

q_comp_homo_incis_norm = q_component_homo_incis/q_homo_incis
q_comp_homo_incis_norm[q_homo_incis==0] = comp[q_homo_incis==0]
# To better visualise composition in channels we now multiply it by the grid 
# `channels'. 
q_comp_homo_incis_norm_channel = np.log10(q_comp_homo_incis_norm) * is_drainage

# Where cumulative sedflux is 0, q_comp_homo_incis is set to bdrock comp
mg.add_field('node','homo_incis_sed__Mg',q_comp_homo_incis_norm,noclobber=False)
mg.add_field('node','homo_incis_log_sed__Mg',np.log10(q_comp_homo_incis_norm),noclobber=False)
mg.add_field('node','homo_incis_log_sed__Mg_channel',q_comp_homo_incis_norm_channel,noclobber=False)

# Composition of sediment with homogeneous incision assumption
print("Sediment Mg concentration (homogeneous incision)")
imshow_grid(mg,'homo_incis_log_sed__Mg_channel',output= True,plot_name ="Sediment Mg concentration (Homo. Incis.)",
            cmap='viridis',limits=(3.2,4.4),var_name="Mg conc. / log10")

print("Finished")
print(datetime.datetime.now().time())
