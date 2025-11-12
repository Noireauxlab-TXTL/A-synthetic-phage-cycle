# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 22:33:02 2024

@author: Ziane
"""

import os
import pandas as pd
import numpy as np
np.seterr(divide='ignore', invalid='ignore')
import matplotlib
import matplotlib.pyplot as plt
import math
import copy

plt.close('all')
plt.ioff()

dir_root = (((os.getcwd()).replace("\\", "/")).split('Script'))[0]
dir_source = f"{dir_root}Source/"
dir_data = f"{dir_root}Output/"

scale = 0.29 # pixel to um
background_0 = 500
H = 59/scale

matplotlib.rcParams.update({'font.size': 15})

process_flag = False # display plots radial profiles of fluorescence intensity with spherical fit, and radial profiles of concentration
adjust_flag = False # display plots of the details of spherical fits for each object

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
        weight_of_pos = 5
        weight_of_neg = 1
        sign = math.copysign(1, v)
        return sign*(weight_of_pos - weight_of_neg)/2 + (weight_of_pos + weight_of_neg)/2
    return sum(abs(weight_function(i)*i) for i in u)

illustrative_examples_list = [['TF-IN_only', 'with_PEG-PE', ['9_1', '9_2', '9_3', '9_4']], ['TF-IN_and_OUT', 'with_PEG-PE', ['1_1', '1_2', '1_3', '1_4']]]
def condition_is_in_sup(_exp_name): # find the objects used in figures S10-11
    res = False
    for i in range(len(illustrative_examples_list)):
        if illustrative_examples_list[i][0] in _exp_name and illustrative_examples_list[i][1] in _exp_name and any(u in _exp_name for u in illustrative_examples_list[i][2]):
            res = True
    return res

r_int_diff_screening_range = 6
excluded_center_range_pixels = 10
excluded_outer_range_fraction = 1/7
int_diff_peak_cutoff = 0.8
int_diff_peak_edge_cutoff = 0.15
background_radius_cutoff = 2
localization_edge_cutoff  = 0.8
center_cutoff = 0.2
# localization_edge_cutoff  = 0.6

writer = pd.ExcelWriter(f"{dir_data}Consolidated_data_localization.xlsx", engine = 'xlsxwriter')
data_consolidated = []
data_consolidated_detail = []
column_names = ['Radius [um]', 'Localization']

# start_time = time.time()
for condition_TF in [_ for _ in os.listdir(dir_source) if os.path.isdir(os.path.join(dir_source, _))]: # Processing of each condition
    dir_condition = f"{dir_source}{condition_TF}/" 
    for condition_PEG in [_ for _ in os.listdir(dir_condition) if os.path.isdir(os.path.join(dir_condition, _))]: # Processing of each sub condition
        dir_rad_dist = f"{dir_condition}{condition_PEG}/Radial intensity distributions/"
        
        data_concatenated = []
        data_replicated = []
        
        list_rad_dist = [_ for _ in os.listdir(dir_rad_dist) if _.endswith('.csv')]
            
        for object_rad_dist in list_rad_dist: # Processing of each object
            name_rad_dist = object_rad_dist[0:-4]
            object_name = f"{condition_TF} - {condition_PEG} - {name_rad_dist}"
            
            print('\n' +f"***** {object_name} *****")
            
            df_rad_dist = pd.read_csv(f"{dir_rad_dist}{object_rad_dist}")
            
            intensity = df_rad_dist['Normalized_Integrated_Intensity']
            radius = df_rad_dist['Radius_[pixels]']
            
            intensity_diff = np.diff(intensity)/np.diff(radius)/(np.nanmax(intensity) - np.nanmin(intensity))
            
            ## Find the main peak, where the fluorescence intensity shows the sharpest decrease
            radius_screening_range = [u for u in range(excluded_center_range_pixels,int(len(intensity_diff)*(1 - excluded_outer_range_fraction)))]
            min_intensity_diff = np.nanmin(intensity_diff[radius_screening_range])
            min_intensity_diff_idx = [u for u in range(len(intensity_diff)) if intensity_diff[u] == min_intensity_diff][0]
            
            ## Find all the possible peaks, where the fluorescence intensity's decrease is with (1-int_diff_peak_cutoff)% away from the sharpest decrease 
            rad_from_diff_idx = [u for u in range(len(intensity_diff)) if intensity_diff[u] <= int_diff_peak_cutoff*min_intensity_diff if u in radius_screening_range if np.diff(intensity_diff)[u] > 0 if np.diff(intensity_diff)[u-1] < 0 if u > min_intensity_diff_idx*0.5]
            
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
                aa = ([u for u in range(peak_edge_high_onset_idx[uu], min(peak_edge_high_onset_idx[uu] + r_int_diff_screening_range, len(intensity_diff)-1)) if np.diff(intensity_diff)[u]>0 if np.diff(intensity_diff)[u+1]<0])
                aa.append(min(peak_edge_high_onset_idx[uu] + r_int_diff_screening_range, len(intensity_diff)-1))
                rad_cutoff_high_idx[uu] = min(aa) + 1 #min(([u for u in range(peak_edge_high_onset_idx[uu], min(peak_edge_high_onset_idx[uu] + r_int_diff_screening_range, len(intensity_diff)-1)) if np.diff(intensity_diff)[u]>0 if np.diff(intensity_diff)[u+1]<0]).append(min(peak_edge_high_onset_idx[uu] + r_int_diff_screening_range, len(intensity_diff)-1)))+1
            
            if adjust_flag:                
                ## Plot and save the normalized derivative of the fluorescence intensity, with horizontal lines showing various levels, and vertical lines showing peaks and search boundaries
                plt.figure()
                plt.plot(radius[0:-1], intensity_diff/abs(min_intensity_diff), color = 'g')
                for uu in range(len(rad_from_diff_idx)):
                    plt.axvline(x = radius[rad_from_diff_idx[uu]], color = 'r', linestyle = ':') # peaks
                    plt.axvline(x = radius[rad_cutoff_low_idx[uu]], color = 'k', linestyle = ':') # lower search boundaries
                    plt.axvline(x = radius[rad_cutoff_high_idx[uu]], color = 'k', linestyle = ':') # higher search boundaries
                plt.axhline(y = -int_diff_peak_edge_cutoff, color = 'b', linestyle = ':')
                plt.axhline(y = -int_diff_peak_cutoff, color = 'r', linestyle = ':')
                plt.title(object_name)
                plt.xlabel('Radius [pixel]')
                plt.ylabel('Normalized fluo. int. variation [a.u.]')
                plt.ylim((-1,0))
                plt.savefig(f"{dir_data}{condition_TF}/{condition_PEG}/{object_name}_rad_search_diff_norm.svg", bbox_inches='tight', pad_inches = 0)
                plt.close()
                
                ## Plot and save the fluorescence intensity, with vertical lines showing inflection points and search boundaries    
                plt.figure()
                plt.plot(radius, intensity, color = 'g')
                for uu in range(len(rad_from_diff_idx)):
                    plt.axvline(x = radius[rad_from_diff_idx[uu] + 1], color = 'r', linestyle = ':')
                    plt.axvline(x = radius[rad_cutoff_low_idx[uu] + 1], color = 'k', linestyle = ':')
                    plt.axvline(x = radius[rad_cutoff_high_idx[uu] + 1], color = 'k', linestyle = ':')
                plt.title(object_name)
                plt.xlabel('Radius [pixel]')
                plt.ylabel('Fluo. int. [a.u.]')
                plt.savefig(f"{dir_data}{condition_TF}/{condition_PEG}/{object_name}_rad_search_intensity.svg", bbox_inches='tight', pad_inches = 0)
                plt.close()

            
            ## Screening range for the radius as the collection of the ranges around all the inflecition points found
            r_range = []
            for uu in range(len(rad_from_diff_idx)):
                for vv in range(rad_from_diff_idx[uu] + 1, rad_cutoff_high_idx[uu] + 1):
                    if vv not in r_range:
                        r_range.append(radius[vv])
           
            ## Screening range for the background as the range of intensity with background_radius_cutoff% radius away from the higher boundary of the peaks found
            I_o_min = min(intensity[max(rad_cutoff_high_idx) : min(int(round(max(rad_cutoff_high_idx)*background_radius_cutoff)), len(df_rad_dist) - 1)])
            I_o_max = intensity[max(rad_cutoff_high_idx)]
            I_o_range = increments(I_o_min, I_o_max)
            
            ## Screening range for the liposome intensity as the range of intensity between the minimum intensity over half the range of the highest boundary and the maximum intensity over half the range of the lowest boundary
            I_i_min = np.min(intensity[0:int(np.ceil(max(rad_cutoff_high_idx)/2))])
            I_i_max = np.max(intensity[0:int(np.ceil(min(rad_cutoff_low_idx)/2))])
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
            
            fluo_ref = [fluo_sphere(u, fitting_parameters) for u in radius]
            data = [(intensity[u] - fluo_ref[u] + abs(fluo_ref[u] - fitting_parameters[0]))/abs(fluo_ref[u] - fitting_parameters[0]) for u in range(max([v for v in range(len(df_rad_dist)) if radius[v] < fitting_parameters[2]]))] # exclude the last point within the liposome to avoid artifacts
            center_cutoff_idx = int(round(len(data)*center_cutoff))
            
            if process_flag or condition_is_in_sup(object_name): # Display the results of the object analysis
                # Plot fluorescence intensity vs radial position for each object, with the fitted uniform distribution and the fitting ranges used
                plt.figure()
                plt.ylabel('Fluorescence intensity [a.u.]')
                plt.xlabel('Radius [$\mu$m]')
                plt.title(object_name)
                plt.plot(radius*scale, intensity, color='g', label = 'Experiment')
                plt.plot(radius*scale, [fluo_sphere(u, fitting_parameters) for u in radius], color='r', label = 'Spherical fit')
                plt.xlim(left = 0, right=round(fitting_parameters[2]*scale,1)*1.4)
                plt.legend(prop={'size': 15}, loc='lower left')
                plt.savefig(f"{dir_data}{condition_TF}/{condition_PEG}/{object_name}_rad_dist.svg", bbox_inches='tight', pad_inches = 0)
                plt.close()
        
            concentration_alpha = [(intensity[u] - fitting_parameters[0])/(2*height(radius[u], fitting_parameters[2])) + (fitting_parameters[0] - background_0)/H for u in range(len(data))]
            
            interfacial_localization = max([concentration_alpha[u] for u in range(int(round(len(data)*0.5)), len(data))])/(sum([concentration_alpha[u]*radius[u] for u in range(int(round(len(data)*center_cutoff)))])/sum([radius[u] for u in range(int(round(len(data)*center_cutoff)))]))
            
            if process_flag or condition_is_in_sup(object_name): # Display the results of the object analysis
                # Plot concentration*alpha vs radial position for each object
                plt.figure()
                plt.ylabel('Concentration [a.u.]')
                plt.xlabel('Radius [$\mu$m]')
                plt.title(object_name)
                plt.plot(radius[0:len(data)]*scale, concentration_alpha, label = 'Inside liposome - loc. coeff.: ' + f'{interfacial_localization:.2f}')
                plt.ylim(bottom = 0, top = np.nanmax(concentration_alpha)*1.2)
                plt.axhline(y = (fitting_parameters[0] - background_0)/H, color = 'k', linestyle = ':', label = 'Outside liposome')
                plt.legend(prop={'size': 15}, loc='lower left')
                plt.savefig(f"{dir_data}{condition_TF}/{condition_PEG}/{object_name}_Concentration.svg", bbox_inches='tight', pad_inches = 0)
                plt.close()
            
            if condition_is_in_sup(object_name):
                writer_fig_SI = pd.ExcelWriter(f"{dir_data}{condition_TF}/{condition_PEG}/{object_name}_details.xlsx", engine = 'xlsxwriter')
                data_object = {'Radius [um]': (radius*scale).to_list(), 'Fluo. int. [a.u.]': intensity.to_list(), 'Spherical fit [a.u.]': [fluo_sphere(u, fitting_parameters) for u in radius], 'Concentration [a.u.]' : concentration_alpha  + [np.nan]*(len(radius) - len(concentration_alpha))}
                df_data_object = pd.DataFrame(data_object)
                df_data_object.to_excel(writer_fig_SI)
                writer_fig_SI.close()
            
            data_replicated.append([fitting_parameters[2]*scale, interfacial_localization])
                
        uuu = copy.deepcopy(np.array(data_replicated))
        if np.shape(data_concatenated)[0] == 0:
            data_concatenated = copy.deepcopy(uuu)
        else:
            radius_to_concat = copy.deepcopy(uuu[:,0][...,None])
            array_to_concat_vert = np.hstack((radius_to_concat,np.array([[float('nan')]*(np.shape(data_concatenated)[1] - 1)]*len(radius_to_concat))))
            data_to_add = copy.deepcopy(uuu[:,1:3])
            array_to_concat_hor = np.vstack((np.array([[float('nan')]*2]*np.shape(data_concatenated)[0]), data_to_add))
            data_concatenated = np.vstack((data_concatenated, array_to_concat_vert))
            data_concatenated = np.hstack((data_concatenated, array_to_concat_hor))
        
        data_summary = [f"{condition_TF} - {condition_PEG}", np.mean(np.array(data_replicated)[:,1]), np.std(np.array(data_replicated)[:,1]), len(np.array(data_replicated))]
        data_consolidated = data_consolidated + [data_summary]
        data_summary_detail = [f"{condition_TF} - {condition_PEG}", np.array(data_replicated)[:,1]]
        data_consolidated_detail = data_consolidated_detail + [data_summary_detail]
        df_concatenated = pd.DataFrame(data_concatenated, index=[u + 1 for u in range(len(data_concatenated))], columns=column_names)
        df_concatenated.to_excel(writer, sheet_name = f"{condition_TF} - {condition_PEG}")
   
samples_sizes = [len(u) for u in np.array(data_consolidated_detail, dtype=object)[:, 1]]
       
df_consolidated = pd.DataFrame(data_consolidated, index=[u for u in range(len(data_consolidated))], columns=['Condition', 'Localization - Average', 'Localization - StD', 'Sample size'])
df_consolidated.to_excel(writer, sheet_name = 'Summary')

df_consolidated_detail = pd.DataFrame(data_consolidated_detail, index=[u for u in range(len(data_consolidated_detail))], columns=['Condition', 'Localization - detail'])
df_consolidated_detail.to_excel(writer, sheet_name = 'Summary - detail')

writer.close()