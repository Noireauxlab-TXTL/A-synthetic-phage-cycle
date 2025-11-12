# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 05:32:42 2024

@author: Ziane
"""

import os
import pandas as pd
import numpy as np
np.seterr(divide='ignore', invalid='ignore')
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.image import imread
import math
import copy

plt.close('all')
plt.ioff()

scale = 0.325 # pixel to um
background_0 = 500
H = 59/0.29*scale

matplotlib.rcParams.update({'font.size': 15})

adjust_flag = False # plot the details of the search of the fingerprint of the radial fluorescence intensity
process_flag = False # plot the details of the curve fitting process

root_directory = (((os.getcwd()).replace("\\", "/")).split('Script'))[0]

directory_data = f"{root_directory}Output/"
directory_data_details = f"{directory_data}details/"
if adjust_flag or process_flag:
    if not os.path.isdir(directory_data_details):
        os.mkdir(directory_data_details)

def height(r, R):
    return R*math.sqrt(max((1 - (r/R)**2), 0))

def fluo_sphere(r, u):
    I_o, I_i, radius = u
    return I_o + (I_i - I_o)*height(r, radius)/radius

def increments(a, b):
    N = 10
    return np.linspace(a, b, N)

def custom_norm(u):
    def weight_function(v):
        weight_of_pos = 3
        weight_of_neg = 1
        sign = math.copysign(1, v)
        return sign*(weight_of_pos - weight_of_neg)/2 + (weight_of_pos + weight_of_neg)/2
    return np.nanmean([abs(weight_function(i)*i) for i in u])

r_int_diff_screening_range = 5
excluded_center_range_pixels = 5
# excluded_outer_range_fraction = 1/7
excluded_outer_range_fraction = 0.05
int_diff_peak_cutoff = 0.8
int_diff_peak_edge_cutoff = 0.15
background_radius_cutoff = 2
localization_edge_cutoff  = 0.8
center_cutoff = 0.2
# localization_edge_cutoff  = 0.6
    
list_conditions = ['3000', '300', '30', '3']
column_names_full = ['Radius [um]'] + [f"mC concentration  - MOI {u} [a.u.]" for u in list_conditions] + ['Replicate', 'Object']
target_file = f"{directory_data}data.xlsx"
writer = pd.ExcelWriter(target_file, engine = 'xlsxwriter')
data_formated = {key: [] for key in column_names_full}
for condition in list_conditions:
    directory_condition = f"{root_directory}Source/{condition}/"
    
    column_names = ['Image', 'Object','GFP_Radius', 'GFP_Imax', 'GFP_Imin', 'GFP_Imax-Imin', 'GFP_fit_Imax', 'GFP_fit_BG', 'GFP_a_conc', 'TexRed_Imax', 'TexRed_Imin', 'TexRed_Imax-Imin', 'TexRed_Imax-Imin/radius']
    
    plt.figure(100)
    
    data = np.array([])
    img_labels = np.array([])
    
    radius_list = []
    mC_conc_list = []
    replicate_list = []
    object_list = []
    
    for fold_1 in [name for name in os.listdir(directory_condition) if os.path.isdir(os.path.join(directory_condition, name)) and 's' in name]: # Processing of each condition
        dir_1 = f"{directory_condition}{fold_1}/" 
        dir_rad_dist = f"{dir_1}Radial GFP intensities/"
        dir_crop_images = f"{dir_1}Cropped images in mC channel/"
        
        list_images = [u for u in os.listdir(dir_crop_images) if 'mC_crop_s' in u]
        list_rad_dist = [u for u in os.listdir(dir_rad_dist) if 'Values_s' in u]
        
        for idx in range(1, len(list_rad_dist)+1):
            if f"{fold_1} - {idx}" in ['s8 - 9', 's15 - 1']:
                continue
                    
            data_simple = np.array([0.0]*len(column_names))
            
            data_simple[1]=idx
            
            print(f"{fold_1} - {idx}")
            
            object_list += [idx]
            replicate_list += [fold_1]
            
            df_rad_dist = pd.read_csv(f"{dir_rad_dist}Values_{fold_1}_object_{idx}.csv")
            df_rad_dist.columns = ['Radius_[pixels]', 'Normalized_Integrated_Intensity']    
            intensity = df_rad_dist['Normalized_Integrated_Intensity']
            radius = df_rad_dist['Radius_[pixels]']
            
            intensity_max = np.nanmax(intensity)
            intensity_min = np.nanmin(intensity)
            
            data_simple[3] = intensity_max
            data_simple[4] = intensity_min
            data_simple[5] = intensity_max - intensity_min
            
            intensity_diff = np.diff(intensity)/np.diff(radius)/(np.nanmax(intensity) - np.nanmin(intensity))
            
            ## Find the main peak, where the fluorescence intensity shows the sharpest decrease
            radius_screening_range = [u for u in range(excluded_center_range_pixels,int(len(intensity_diff)*(1 - excluded_outer_range_fraction)))]
            min_intensity_diff = np.nanmin(intensity_diff[radius_screening_range])
            min_intensity_diff_idx = [u for u in range(len(intensity_diff)) if intensity_diff[u] == min_intensity_diff][0]
            
            ## Find all the possible peaks, where the fluorescence intensity's decrease is with (1-int_diff_peak_cutoff)% away from the sharpest decrease 
            rad_from_diff_idx = [u for u in range(len(intensity_diff)-1) if intensity_diff[u] <= int_diff_peak_cutoff*min_intensity_diff if u in radius_screening_range if np.diff(intensity_diff)[u] > 0 if np.diff(intensity_diff)[u-1] < 0 if u > min_intensity_diff_idx*0.5]
            
            ## Find the lower boundary of the peak, defined as the location where the decrease is the smallest within r_int_diff_screening_range away from the considered peak
            #### Intensities of those lower boundaries
            intensity_max_in_range_low = [max(intensity_diff[u - r_int_diff_screening_range: u]) for u in rad_from_diff_idx]
            
            #### Indices of those lower boundaries
            rad_cutoff_low_idx = [[u for u in range(rad_from_diff_idx[v] - r_int_diff_screening_range, rad_from_diff_idx[v]) if intensity_diff[u] == intensity_max_in_range_low[v]][0] for v in range(len(rad_from_diff_idx))]
                
            ## Find the higher boundary of the peak, defined as the location where the decrease is the smallest within r_int_diff_screening_range away from the considered peak
            #### Find indices of all the points where the fluorescence intensity decreases at a rate at least int_diff_peak_edge_cutoff% of its largest decrease, around each of the peaks foud
            peak_edge_high_onset_idx = [max([u for u in range(rad_from_diff_idx[v], rad_from_diff_idx[v] + r_int_diff_screening_range) if intensity_diff[u] <= int_diff_peak_edge_cutoff*min_intensity_diff]) for v in range(len(rad_from_diff_idx))]
            #### Find the first maximum of the derivative of the fluorescence intensity after the boundary found in the code line above
            rad_cutoff_high_idx = [0]*len(rad_from_diff_idx)
            for uu in range(len(rad_from_diff_idx)):
                aa = ([u for u in range(peak_edge_high_onset_idx[uu], min(peak_edge_high_onset_idx[uu] + r_int_diff_screening_range, len(intensity_diff)-2)) if np.diff(intensity_diff)[u]>0 if np.diff(intensity_diff)[u+1]<0])
                aa.append(min(peak_edge_high_onset_idx[uu] + r_int_diff_screening_range, len(intensity_diff)-1))
                rad_cutoff_high_idx[uu] = min(aa) + 1 #min(([u for u in range(peak_edge_high_onset_idx[uu], min(peak_edge_high_onset_idx[uu] + r_int_diff_screening_range, len(intensity_diff)-1)) if np.diff(intensity_diff)[u]>0 if np.diff(intensity_diff)[u+1]<0]).append(min(peak_edge_high_onset_idx[uu] + r_int_diff_screening_range, len(intensity_diff)-1)))+1
            
            if adjust_flag: # Display the results of the analysis of the derivative of the fluorecence intensity
                plt.figure()
                plt.plot(radius[0:-1], intensity_diff, color = 'g')
                
                for u in range(len(rad_from_diff_idx)):
                    plt.axvline(x = radius[rad_from_diff_idx[u]], color = 'r', linestyle = ':')
                    
                    plt.axvline(x = radius[rad_cutoff_low_idx[u]], color = 'k', linestyle = ':')
                    plt.axvline(x = radius[rad_cutoff_high_idx[u]], color = 'k', linestyle = ':')
                    
                plt.title(f"{condition} - {fold_1} - {idx}")
                plt.xlabel('Radius [pixel]')
                
                plt.savefig(f"{directory_data_details}{condition}_{fold_1}_{idx}_rad_search_diff.png", bbox_inches='tight', pad_inches = 0)
                plt.close()
                
                ## Plot and save the normalized derivative of the fluorescence intensity, with horizontal lines showing various levels, and vertical lines showing peaks and search boundaries
                plt.figure()
                plt.plot(radius[0:-1], intensity_diff/abs(min_intensity_diff), color = 'g')
                for uu in range(len(rad_from_diff_idx)):
                    plt.axvline(x = radius[rad_from_diff_idx[uu]], color = 'r', linestyle = ':') # peaks
                    plt.axvline(x = radius[rad_cutoff_low_idx[uu]], color = 'k', linestyle = ':') # lower search boundaries
                    plt.axvline(x = radius[rad_cutoff_high_idx[uu]], color = 'k', linestyle = ':') # higher search boundaries
                    
                plt.axhline(y = -int_diff_peak_edge_cutoff, color = 'b', linestyle = ':')
                
                plt.axhline(y = -int_diff_peak_cutoff, color = 'r', linestyle = ':')
                
                plt.title(f"{condition} - {fold_1} - {idx}")
                plt.xlabel('Radius [pixel]')
                plt.ylim((-1,0))
                
                plt.savefig(f"{directory_data_details}{condition}_{fold_1}_{idx}_rad_search_diff_norm.png", bbox_inches='tight', pad_inches = 0)
                plt.close()
                
                ## Plot and save the fluorescence intensity, with vertical lines showing inflection points and search boundaries    
                plt.figure()
                plt.plot(radius, intensity, color = 'g')
                for uu in range(len(rad_from_diff_idx)):
                    plt.axvline(x = radius[rad_from_diff_idx[uu] + 1], color = 'r', linestyle = ':')
                    plt.axvline(x = radius[rad_cutoff_low_idx[uu]], color = 'k', linestyle = ':')
                    plt.axvline(x = radius[rad_cutoff_high_idx[uu]], color = 'k', linestyle = ':')
                    
                plt.title(f"{condition} - {fold_1} - {idx}")
                plt.xlabel('Radius [pixel]')
                plt.savefig(f"{directory_data_details}{condition}_{fold_1}_{idx}_rad_search_intensity.png", bbox_inches='tight', pad_inches = 0)
                plt.close()

            ## Screening range for the radius as the collection of the ranges around all the inflecition points found
            r_range = []
            for uu in range(len(rad_from_diff_idx)):
                for vv in range(rad_from_diff_idx[uu] + 1, rad_cutoff_high_idx[uu] + 1):
                    if vv not in r_range:
                        r_range.append(radius[vv])
           
            ## Screening range for the background as the range of intensity with background_radius_cutoff% radius away from the higher boundary of the peaks found
            I_o_min = min(intensity[max(rad_cutoff_high_idx) : min(int(round(max(rad_cutoff_high_idx)*background_radius_cutoff)), len(df_rad_dist))])
            I_o_max = intensity[max(rad_cutoff_high_idx)]
            if math.isnan(I_o_min) or math.isnan(I_o_max):
                I_o_min = np.nanmin(intensity)
                I_o_max = I_o_min + 0.2*(np.nanmax(intensity)-I_o_min)
                
            I_o_range = increments(I_o_min, I_o_max)
            
            ## Screening range for the liposome intensity as the range of intensity between the minimum intensity over half the range of the highest boundary and the maximum intensity over half the range of the lowest boundary
            I_i_min = np.min(intensity[0:int(np.ceil(max(rad_cutoff_high_idx)/2))])
            I_i_max = np.max(intensity[0:int(np.ceil(min(rad_cutoff_low_idx)/2))])
            if np.isnan(I_i_max):
                I_i_max = I_i_min*1.1
            I_i_range = increments(I_i_min, I_i_max)
            
            distance_max = float('inf')
            i = 0
            for r_ in r_range:
                i += 1
                j = 0
                for I_o in I_o_range:
                    j += 1
                    k = 0
                    for I_i in I_i_range:
                        k += 1
                        tested_parameters = [I_o, I_i, r_]
                        fluo_ref = [fluo_sphere(u, tested_parameters) for u in radius]
                        
                        vector = [(fluo_ref[v]-intensity[v]) for v in range(max(rad_cutoff_high_idx))]
                        distance = custom_norm(vector)
                        
                        if distance < distance_max:
                            fitting_parameters = [I_o, I_i, r_]
                            distance_max = min(distance_max, distance)
            
            if process_flag: # Display the results of the object analysis
                # Plot fluorescence intensity vs radial position for each object, with the fitted uniform distribution and the fitting ranges used
                plt.figure()
                plt.ylabel('Fluorescence intensity [a.u.]')
                plt.xlabel('Radius [$\mu$m]')
                plt.title(f"{condition} - {fold_1} - {idx}")
                plt.plot(radius*scale, intensity, color='g', label = 'Experiment')
                plt.plot(radius*scale, [fluo_sphere(u, fitting_parameters) for u in radius], color='r', label =
                          'Best fit\n'+
                          'Background: ' + str(round(fitting_parameters[0],1)) + ' [a.u.]\n'+
                          'Center: ' + str(round(fitting_parameters[1],1)) + ' [a.u.]\n'+
                          'Radius: ' + str(round(fitting_parameters[2]*scale,1)) + ' $\mu$m\n'+
                          'Concentration difference: ' + str(round((fitting_parameters[1] - fitting_parameters[0])/fitting_parameters[2],3)))
                plt.axvline(x = min(r_range)*scale, color = 'b', linestyle = ':')
                plt.axvline(x = max(r_range)*scale, color = 'b', linestyle = ':')
                plt.axhline(y = min(I_o_range), color = 'k', linestyle = '--')
                plt.axhline(y = max(I_o_range), color = 'k', linestyle = '--')
                plt.axhline(y = min(I_i_range), color = 'k', linestyle = ':')
                plt.axhline(y = max(I_i_range), color = 'k', linestyle = ':')
                plt.legend(prop={'size': 10})
                plt.savefig(f"{directory_data_details}{condition}_{fold_1}_{idx}_rad_dist.png", bbox_inches='tight', pad_inches = 0)
                plt.close()
            
            data_simple[6] = fitting_parameters[1]
            data_simple[7] = fitting_parameters[0]
            data_simple[8] = (fitting_parameters[1] - fitting_parameters[0])/(2*fitting_parameters[2]) + (fitting_parameters[0] - background_0)/(H/scale)
            data_simple[2] = fitting_parameters[2]
            
            mC_image = imread(f"{dir_crop_images}mC_crop_{fold_1}_object_{idx}.tif")
            
            img_size = np.size(mC_image)
            img_shape = mC_image.shape
            
            fraction_min = 100/2
            fraction_max = 100/10
            
            portion_min_size = int(img_size/fraction_min)
            
            background_pixels = np.partition(mC_image.flatten(), portion_min_size)
            background = np.mean(background_pixels[:portion_min_size])
            
            intensity = np.mean(mC_image[int(np.floor(img_shape[0]/2 - r_/fraction_max)):int(np.ceil(img_shape[0]/2 + r_/fraction_max)), int(np.floor(img_shape[1]/2 - r_/fraction_max)):int(np.ceil(img_shape[1]/2 + r_/fraction_max))])
            
            data_simple[9] = intensity
            data_simple[10] = background
            data_simple[11] = intensity - background
            data_simple[12] = (intensity - background)/r_
            
            mC_conc_list += [data_simple[12]]
            radius_list += [r_*scale]
                        
            if len(data) == 0:
                data = copy.deepcopy(data_simple)
            else:
                if np.mean(data_simple) > 0:
                    data = np.vstack((data, data_simple))
                
            if len(img_labels) == 0:
                img_labels = np.array([fold_1])
            else:
                if np.mean(data_simple) > 0:
                    img_labels = np.vstack((img_labels, np.array([fold_1])))
    
    
    data_formated['Radius [um]'] += radius_list
    data_formated['Replicate'] += replicate_list
    data_formated['Object'] += object_list
    data_formated[f"mC concentration  - MOI {condition} [a.u.]"] += mC_conc_list
    for u in list_conditions:
        if u != condition:
            data_formated[f"mC concentration  - MOI {u} [a.u.]"]+= [np.nan]*len(radius_list)
    
    plt.figure(100)
    plt.scatter(data[:,2]*scale, data[:,12], label = f"MOI: {condition}")

df = pd.DataFrame(data_formated)
df.to_excel(writer)
writer.close() 

plt.figure(100)
plt.ylabel('mCherry concentration [a.u.]')
plt.xlabel('Radius [$\mu$m]')
plt.title('Effect of liposome size on infectivity')
plt.legend(prop={'size': 10})
plt.savefig(f"{directory_data}recap_lin-lin.svg", bbox_inches='tight', pad_inches = 0)
plt.close()