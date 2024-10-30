#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 23:39:24 2024

@author: gruma-r
"""

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import os
import cartopy.crs as ccrs
from cartopy.feature import ShapelyFeature
import geopandas as gpd
from datetime import datetime
# import glob

data_directory = '/home/gruma-r/Área de Trabalho/nc_era5_data/jato'
files = os.listdir(data_directory)

netcdf_files = [f for f in files if os.path.isfile(os.path.join(data_directory, f))]
netcdf_files.sort()

u_wind_var = 'u'
v_wind_var = 'v'
pressure_level = 200 # 200 hPa

#Shapefile path
shapefile_path = '/home/gruma-r/Área de Trabalho/Lim_america_do_sul_2021.shp'
# Load Brazilian states shapefile
brazil_states = gpd.read_file(shapefile_path)
brazil_feature = ShapelyFeature(brazil_states['geometry'], ccrs.PlateCarree(), edgecolor='dimgrey', facecolor='none')

# Folder where the figures will be saved
output_folder = '/home/gruma-r/Imagens'  # Change this to your desired output folder
os.makedirs(output_folder, exist_ok=True)

start_date_range = datetime(2024, 1, 16)
end_date_range = datetime(2024, 1, 17)

# Define a fixed range for wind speed
# vmin = 0  # Minimum wind speed
vmin = 30
vmax = 90 # tava 70 # Maximum wind speed (adjust this value based on your data)

# Define levels and ticks for colorbar
# ticks = np.arange(vmin, vmax + 5, 5)
ticks = np.arange(vmin, vmax + 10, 10)


# levels = np.arange(vmin, vmax + 1, 1)

# levels = np.arange(0, 75, 5)
levels = np.arange(30, 90 + 10, 10)

print(start_date_range)
print(end_date_range)

for file_name in netcdf_files:
    file_path = os.path.join(data_directory, file_name)
    
    with nc.Dataset(file_path, 'r') as netcdf_file:
        # Get the data of latitude, longitude, and the variable
        latitude = netcdf_file.variables['latitude'][:]
        longitude = netcdf_file.variables['longitude'][:]
                
        u_wind_data = netcdf_file.variables[u_wind_var][:, :, :]
        v_wind_data = netcdf_file.variables[v_wind_var][:, :, :]

        time_variable = netcdf_file.variables['time']
        # Convert the unit to datetime
        time_values = time_variable[:]
        time_units = time_variable.units
        calendar = 'gregorian'  # Ensure 'gregorian' calendar is used
        custom_time = nc.num2date(time_values, units=time_units, calendar=calendar)
            
        start_index = np.argmin(np.abs(custom_time - start_date_range))
        end_index = np.argmin(np.abs(custom_time - end_date_range))
        
        for time_step in range(start_index, end_index + 1):
            u_wind = np.array(u_wind_data[time_step, :, :])
            v_wind = np.array(v_wind_data[time_step, :, :])
            u_wind = np.round(u_wind).astype(int)
            v_wind = np.round(v_wind).astype(int)
            lon, lat = np.meshgrid(longitude, latitude)
            
            # Calcula a wind speed e direcão
            wind_speed = np.sqrt(u_wind**2 + v_wind**2)
            
            # Plot the data using Cartopy
            plt.figure(figsize=(10, 8))
            ax = plt.axes(projection=ccrs.PlateCarree())
            ax.set_extent([-30, -90, -20, -55])
            
            #Gridlines (lat e lon)
            gl = ax.gridlines(draw_labels=True, linewidth=0)
            gl.top_labels = False
            gl.right_labels = False
            gl.xlocator = plt.MultipleLocator(15)
            gl.ylocator = plt.MultipleLocator(15)


            # Draw coastlines and country/state borders
            ax.coastlines(linewidth=0.5, color='dimgrey', zorder=2)
            ax.add_feature(brazil_feature, zorder=2)
            
            # Plot wind speed
            wind_plot = ax.contourf(lon, lat, wind_speed, levels=levels, cmap='hot_r', vmin=vmin, vmax=vmax,  zorder=1)
            
            # Define colorbar ticks
            # ticks = np.arange(vmin, vmax + 1, 5)
            
            # Create colorbar
            cbar = plt.colorbar(wind_plot, ax=ax, orientation='horizontal', aspect=30, shrink=0.76, pad = 0.08)
            cbar.set_label('Wind Speed (m/s)')
            cbar.set_ticks(levels)
            
            # Ensure streamplot arguments are numpy arrays and not masked arrays
            u_wind = np.ma.filled(u_wind, np.nan)
            v_wind = np.ma.filled(v_wind, np.nan)
            
            # Plot wind direction arrows
            # ax.quiver(lon[::10, ::10], lat[::10, ::10], u_wind[::10, ::10], v_wind[::10, ::10], zorder=3, color='black')
            
            # Plot wind direction streamlines
            ax.streamplot(lon, lat, u_wind, v_wind, color='black', linewidth=0.5, density=1, zorder=3)    

            # Format the datetime information for the filename
            formatted_time = custom_time[time_step].strftime("%d-%m-%Y_%H:%M:%S")
            plt.title(f'Velocidade do vento em 200mb às {formatted_time}')

            # Construct the output filename
            output_filename = os.path.join(output_folder, f'wind_200hPa_{formatted_time}.png')

            # Save the figure
            plt.show()
            plt.savefig(output_filename, format='png', bbox_inches='tight')
            plt.close()
        