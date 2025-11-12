# -*- coding: utf-8 -*-
"""
Created on Mon Jul 28 14:26:46 2025

@author: Ziane
"""

import os
from os import listdir

import skimage.io
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd


plt.close('all')
plt.ioff()


dir_root = (((os.getcwd()).replace("\\", "/")).split('Script'))[0]
dir_source = f"{dir_root}Source/"
dir_data = f"{dir_root}Output/"


def sigmoid(x, L, x0, k, b):
    return L / (1 + np.exp(-k * (x - x0))) + b

exp_param_percentile_dict = {}
time_collection_dict = {key: [] for key in listdir(dir_source)}
fraction_min = 100/20
fraction_max = 100/0.1

column_names_full = ['Time [min]', 'mC Replicate 1 [a.u.]', 'mC Replicate 2 [a.u.]', 'mC Average [a.u.]', 'mC Average fit [a.u.]']

# Process each stack, generate a plot for each stack, the average plot for each condition, store data of all the averages in a dictionary
for condition in listdir(dir_source):
    target_file = f"{dir_data}{condition}/data_{condition}.xlsx"
    writer = pd.ExcelWriter(target_file, engine = 'xlsxwriter')
    data_formated = {key: [] for key in column_names_full}
    image_path = f"{dir_source}{condition}/"
    
    if condition == 'T7GFP DNA':
        t_max = 450
        t_step = 10
    elif condition == 'T7mC DNA':
        t_max = 450
        t_step = 10
    elif condition == 'T7GFP phage':
        t_max = 400
        t_step = 15
    elif condition == 'T7mC phage':
        t_max = 800
        t_step = 15
                
    stacks = [f for f in listdir(image_path) if ('red' in f or 'green' in f or 'gfp' in f)]
    
    percentile_max_avg = [0]*int(t_max/t_step + 1)
    
    plt.figure(200)
    rep_idx = 0
                
    for k in stacks:
        replicate_tag = f"{condition} - {k.split('.')[0]}"
        print(replicate_tag)
        img_stack = skimage.io.imread(image_path + k) # load stack of images
        
        stack_percentile_kinetics = [0]*int(t_max/t_step + 1)
        img_size = np.size(img_stack[0])
        portion_min_size = int(img_size/fraction_min)
        portion_max_size = int(img_size/fraction_max)
        noise = np.min(img_stack)
        background_max = 0
        
        for l in range(len(stack_percentile_kinetics)):
            background_pixels = np.partition(img_stack[l].flatten(), portion_min_size)
            background_max = max(np.mean(background_pixels[:portion_min_size]), background_max)
        
        for l in range(len(stack_percentile_kinetics)):
            background_pixels = np.partition(img_stack[l].flatten(), portion_min_size)
            background = np.mean(background_pixels[:portion_min_size])
            
            bright_pixels = -np.partition(-img_stack[l].flatten(), portion_max_size)
            intensity = np.mean(bright_pixels[:portion_max_size])
            
            intensity_adjusted = (intensity - noise)*(background_max - noise)/(background - noise) - (background_max - noise)
            
            stack_percentile_kinetics[l] = intensity_adjusted
            percentile_max_avg[l] += intensity_adjusted
        
        # relevant channel
        if 'red' in replicate_tag:
            color = 'r'
            reporter = 'mC'
        elif ('green' in k or 'gfp' in replicate_tag):
            color = 'g'
            reporter = 'GFP'
        
        if 'DNA' in condition:
            mode = 'DNA'
        elif 'phage' in condition:
            mode = 'phage'
        
        t_ = [_*t_step for _ in range(len(stack_percentile_kinetics))]
        
        # Sigmoidal best fit
        p0 = [max(stack_percentile_kinetics) - min(stack_percentile_kinetics), np.median(t_), max(np.diff(stack_percentile_kinetics))/t_step/(max(stack_percentile_kinetics) - min(stack_percentile_kinetics))*4, min(stack_percentile_kinetics)]
        try:
            popt, pcov = curve_fit(sigmoid, t_, stack_percentile_kinetics, p0 = p0)
            t_half_list = [t_[u] for u in range(len(t_)) if (stack_percentile_kinetics[u] - popt[3])/popt[0] > 0.5]
            if len(t_half_list) > 0:
                t_half = min(t_half_list)
            else:
                t_half = []
            time_collection_dict[condition] += [t_half]
        except:
            print('oups')
        
        # Time evolution for each replicate
        plt.figure()
        plt.plot(t_, stack_percentile_kinetics, color = color)
        plt.plot(t_, [sigmoid(_, popt[0], popt[1], popt[2], popt[3]) for _ in t_], color = color, linestyle = '--')
        plt.xlabel('Time [min]')
        plt.ylabel('Gene expression [a.u.]')
        plt.xlim((0, t_max))
        plt.title(replicate_tag)
        plt.tight_layout(pad=1)
        plt.savefig(f"{dir_data}{condition}/{replicate_tag}.svg")
        if not plt.isinteractive():
            plt.close()
        
        rep_idx += 1

        data_formated[f"mC Replicate {rep_idx} [a.u.]"] = stack_percentile_kinetics
        
        plt.figure(200)
        label_submission = f"Replicate {rep_idx}"
        plt.scatter(t_, stack_percentile_kinetics, s = 3, label = label_submission)
        
    for l in range(len(stack_percentile_kinetics)):
        percentile_max_avg[l] = percentile_max_avg[l]/rep_idx # generate the average over all the stacks of a given condition
    
    t_ = [_*t_step for _ in range(len(percentile_max_avg))]
    
    # Sigmoidal best fit
    p0 = [max(percentile_max_avg) - min(percentile_max_avg), np.median(t_), max(np.diff(percentile_max_avg))/t_step/(max(percentile_max_avg) - min(percentile_max_avg))*4, min(percentile_max_avg)]
    popt, pcov = curve_fit(sigmoid, t_, percentile_max_avg, p0 = p0)
        
    plt.figure(200)
    plt.scatter([_*t_step for _ in range(len(percentile_max_avg))], percentile_max_avg, s = 50, c = color, label = 'Average', marker = '+')
    plt.plot([_*t_step for _ in range(len(percentile_max_avg))], [sigmoid(_*t_step, popt[0], popt[1], popt[2], popt[3]) for _ in range(len(percentile_max_avg))], color = color, label = 'Fit')
    plt.xlabel('Time [min]')
    plt.ylabel('Gene expression [a.u.]')
    plt.title(condition)
    plt.xlim((0, t_max))
    plt.legend()
    plt.tight_layout(pad=1)
    plt.savefig(f"{dir_data}{condition}/{condition}.svg", bbox_inches='tight', pad_inches = 0)
    if not plt.isinteractive():
        plt.close()
    
    data_formated['Time [min]'] = t_
    data_formated['mC Average [a.u.]'] = percentile_max_avg
    data_formated['mC Average fit [a.u.]'] = [sigmoid(_*t_step, popt[0], popt[1], popt[2], popt[3]) for _ in range(len(percentile_max_avg))]
    
    df = pd.DataFrame(data_formated)
    df.to_excel(writer)
    writer.close() 
                    
print("\n\nt_half replicates:\n")
for key, value in time_collection_dict.items():
    print(f"{key}: {value}")