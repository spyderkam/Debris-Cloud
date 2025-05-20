#!usr/bin/env python3

__author__ = "Kamyar Modjtahedzadeh"
__date__ = "May 14, 2025"

# This script calculates the portion of fragments inside cloud radii and prints the average.

nRuns = 50
nFrags = 1000
s_Lc = 0.001     # [m]
#m_Lc = 0.11     # [m]
#l_Lc = 0.115     # [m]

from gdmpidc import *
from numpy import linalg, mean
from time import time

## Inside and outside probabilities of s/m/l fragments.
in_ratios_s = []
#in_ratios_m = []
#in_ratios_l = []

start_time = time()

for _ in range(nRuns):
    ## Generate s/m/l clouds and fragments.
    smallCloud = Cloud(s_Lc, nFrags);  vec_r_s = smallCloud.sample_positions()
    #medCloud =  Cloud(m_Lc, nFrags);   vec_r_m = medCloud.sample_positions()
    #largeCloud = Cloud(l_Lc, nFrags);  vec_r_l = largeCloud.sample_positions()
        
    ## Inside and outside numbers of s/m/l fragments.
    inl_s = []
    #inl_m = []
    #inl_l = []

    fragment_data = {
        'small': [smallCloud, inl_s, in_ratios_s],
        #'medium': [medCloud, inl_m, in_ratios_m],
        #'large': [largeCloud, inl_l, in_ratios_l],
    }

    for i in range(nFrags):
        r_s = vec_r_s[i]
        #r_m = vec_r_m[i]
        #r_l = vec_r_l[i]
        
        for size_category in fragment_data.keys():
            fragment_data[size_category].insert(0, locals()[f"r_{size_category[0]}"])
            r, cloud, inside, _, = fragment_data[size_category]
        
            if linalg.norm(r) <= cloud.updated_radius(t=0):
                inside.append(r)
            fragment_data[size_category].pop(0)

    for size_category in fragment_data.keys():
        _, inl, in_ratios = fragment_data[size_category]
        in_ratios.append(len(inl)/cloud.nFrag)

for size_category in fragment_data.keys():
    _, _, in_ratios = fragment_data[size_category]
    print(f"{size_category.capitalize()} fragments (Lc = {locals()[size_category[0]+'_Lc']*100} cm): {mean(in_ratios)*100}% inside")
    
print(f"Process time: {round(time() - start_time, 2)} seconds")

