# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 21:16:18 2019

@author: Xu Jin
"""

from __future__ import absolute_import, division, print_function

import os, sys, urllib2, zipfile, arcpy
import numpy as np
import arcpy.sa

# workspace path, output file path
arcpy.env.workspace = sys.argv[1]
# output data projection
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference("WGS 1984 Web Mercator (Auxiliary Sphere)")
 


# DEFINING FUNCTIONS
# fetch and unzip
def fetch_unzip(url, out_dir):
    '''
    x = fetch_unzip(download_url, unzip_path)
    
    Function:
        fetch zipped file from url 
        save to out_dir
        unzip to out_dir
    
    Return:
        a list of .shp filename in the .zip
    '''
    # fatch data from url
    print("Fetching data from {0}".format(url))
    
    response = urllib2.urlopen(url)
    bin_content = response.read()
    response.close()
    
    out_file = out_dir + os.path.basename(url)
    with open(out_file, "wb") as outf:
        outf.write(bin_content)
    print("Fetch complete")
    
    # unzip data
    with zipfile.ZipFile(out_file, "r") as zip_obj:
        zip_obj.extractall(out_dir)
        namelist = zip_obj.namelist()
        shp_name = []
        
        for filename in namelist:
            print("Unzipping: {0}".format(filename)) 
            # recording shpfile name
            if filename[-3:] == "shp":
                shp_name.append(filename) 
        
        # return shp filename
        return shp_name[0]


# for litter index only    
#def unzip(in_file, out_dir):
#    '''
#    x = unzip(in_file, out_dir)
    
#    Function:
#        unzip in_file to out_dir
#        
#    Return:
#        a list of .shp filename in the .zip
#    '''    
#    # unzip data
#    with zipfile.ZipFile(in_file, "r") as zip_obj:
#       zip_obj.extractall(out_dir)
#        namelist = zip_obj.namelist()
#        shp_name = []
        
#        for filename in namelist:
#            print("Unzipping: {0}".format(filename)) 
#            # recording shpfile name
#            if filename[-3:] == "shp":
#                shp_name.append(filename) 
        
        # return shp filename
#        return shp_name[0]


# data visualization
def visual(in_lyr, in_mxd, field, til_obj, out_path):
    '''
    visual(input_lyr_path, in_mxd, value_field, title_obj, out_png_path)
    
    Function:
        modify the first layer in thefirst data frame
        output png
    
    Note:
        in_mxd = arcpy.mapping.MapDocument(mxd_path)
        mxd should be opened
    '''
    # first dataframe
    data_frames = arcpy.mapping.ListDataFrames(in_mxd)
    df = data_frames[0]
    # get first layer to modify
    lyr_list = arcpy.mapping.ListLayers(in_mxd)
    layer_modify = lyr_list[0]

    src_obj = arcpy.mapping.Layer(in_lyr)
    arcpy.mapping.UpdateLayer(df, layer_modify, src_obj)
    
    # value field
    layer_modify.symbology.valueField = field
    arcpy.RefreshActiveView()
    
    # change title and legend
    elems = arcpy.mapping.ListLayoutElements(in_mxd)  
    # finding title, and legend, only 1 text element/legend in mxd
    for e in elems:
        if "TEXT_ELEMENT" in e.type:
            title = e  
        elif "LEGEND_ELEMENT" in e.type:
            legend = e
    
    # changing title name
    title.text = til_obj    
    # change legend height
    legend.elementHeight = 7
    
    # refresh view
    arcpy.RefreshActiveView()
    
    # export
    out_file = out_path + til_obj + ".pdf"
    arcpy.mapping.ExportToPDF(in_mxd, out_file)
    
    print("Mapped")



    
# VARIABLES
# data created mannually in, the .zip file        
#indir = sys.argv[2] + "/"
# storing everthing generated from thhis script     
outdir = arcpy.env.workspace + "/"

# output projection
coordinate = arcpy.SpatialReference('WGS 1984 Web Mercator (Auxiliary Sphere)')

# maping files url on github.  add raw. before github, delete blob/
original_url = "https://raw.github.com/XuJinGabie/FA_GUS_8066_public/master/FA_original.zip"

tract_url = "http://data.phl.opendata.arcgis.com/datasets/8bc0786524a4486bb3cf0f9862ad0fbf_0.zip"
building_url = "http://data.phl.opendata.arcgis.com/datasets/ab9e89e1273f445bb265846c90b38a96_0.zip"
tree_url = "http://data.phl.opendata.arcgis.com/datasets/957f032f9c874327a1ad800abd887d17_0.zip"
light_url = "http://data.phl.opendata.arcgis.com/datasets/9059a7546b6a4d658bef9ce9c84e4b03_0.zip"

#litter_loc = indir + "litter_index_survey.zip"

# namelist for object name
objs = ["tract", "building", "tree", "light", "litter"]

# namelist for analyst object
analyst_objs = ["building", "tree", "light", "litter"]

# object for correlated with litter
corr_objs = ["building", "tree", "light"]

# grid filename
grid = "Grid.shp"

# report for correlation coefficient
corr_file = outdir + "correlation.txt"

# the mxd file is created mannually as blank file with map elements
map_name = "empty.mxd"

# dictionary for object-projected filename
proj_dic = {}

# dictionary for object-column name
field_dic = {}



# Data
tract_org = fetch_unzip(tract_url, outdir)
building_org = fetch_unzip(building_url, outdir)
tree_org = fetch_unzip(tree_url, outdir)
light_org = fetch_unzip(light_url, outdir)

#litter_org = unzip(litter_loc, outdir)
litter_org= fetch_unzip(original_url, outdir)


# turn building into point
building_point = "Building.shp"
arcpy.FeatureToPoint_management(building_org, building_point, "CENTROID")

# namelist for unprojected filenames
orgs = [tract_org, building_point, tree_org, light_org, litter_org]

del building_point




# PROCESS: projection
# projection
for i, org in enumerate(orgs):
    key = objs[i]
    value = "P_"+ org
    arcpy.Project_management(org, value, coordinate)
    proj_dic[key] = value

print(proj_dic)

# delete unprojected filename variables
del orgs

del tract_org
del building_org
del tree_org
del light_org
del litter_org





# PROCESS: create grid
# get extent
desc = arcpy.Describe(proj_dic["tract"])
# set xmin and ymin as the begining coordinate
begin_coord = str(desc.extent.lowerLeft)
# set 0 point and direction_coord is the angle for turning fishnet, no turn
direction_coord = str(desc.extent.XMin) + " " + str(desc.extent.YMax + 10)
# set xmax and ymax as opposit coordinate
end_coord = str(desc.extent.upperRight)

# create 1000X1000m grid
arcpy.CreateFishnet_management(grid, begin_coord, direction_coord,"1000","1000","0","0", end_coord, "NO_LABELS","#","POLYGON")

# delete temporaly varables
del desc
del begin_coord
del direction_coord
del end_coord



                               

# PROCESS: MAIN ANALYSIS
# count for how many time have calculated for point number 
count_freq = 0

# proj_dic[anl_obj] = filename of object for analyst
for anl_obj in analyst_objs:
    
    # Variables

    intersect_point = "Inter_"+ proj_dic[anl_obj]
    count_table = "Sum_"+ anl_obj + ".dbf"
    
    
    print("Analyzing: {0}".format(proj_dic[anl_obj]))
    
    # Intersect for grid id
    arcpy.Intersect_analysis([grid,proj_dic[anl_obj]], intersect_point, "ALL", "", "INPUT")
    
    if anl_obj == "litter":
        # average litter index
        arcpy.Statistics_analysis(intersect_point,count_table,[["rating","MEAN"]],"FID_Grid") 
        
        # average litter index column name
        #field_names.append("MEAN_FID")
        key = anl_obj
        value = "MEAN_ratin"
        field_dic[key] = value
        print("Data in MEAN_FID")
        
    else:
        # count number of points
        arcpy.Statistics_analysis(intersect_point,count_table,[["FID","COUNT"]],"FID_Grid") 
        
        # counted frequency for 1 time                
        count_freq += 1
        
        # first time counting frequency column name
        if count_freq == 1:
            #field_dic.append("FREQUENCY")
            key = anl_obj
            value = "FREQUENCY"
            field_dic[key] = value
            
            print("Data in FREQUENCY")
        
        # after first time counting frequency column name          
        else:
            #field_dic.append("FREQUENC_{0}".format(i))
            key = anl_obj
            value = "FREQUENC_{0}".format(count_freq-1)
            field_dic[key] = value

            print("Data in FREQUENC_{0}".format(count_freq-1))
            
    # Join to grid
    arcpy.JoinField_management(grid,"FID",count_table,"FID_Grid")
    
    print("Finish Analyzing {0}".format(proj_dic[anl_obj]))

print(field_dic)

del count_freq





# PROCESS: Correlation 
field_names = list(field_dic.values())
corr_table = arcpy.da.FeatureClassToNumPyArray(grid, field_names)
# litter field
litter_field = corr_table[field_dic["litter"]]

# result for stoting coefficient
corr_result = open(corr_file, "a")

# title for content
corr_result.write("Correlation Coefficient for Feature Number-Litter Index: \n")
corr_result.write("\n")

# content
for corr_obj in corr_objs:
    
    # object's field for corrlation
    obj_field = corr_table[field_dic[corr_obj]]
    # correlation: obj-litter
    obj_lit_corr = np.corrcoef(obj_field, litter_field)
    # correlation coefficient
    obj_lit_cef = obj_lit_corr[0][1]
    
    # write in txt
    corr_result.write("{0}-litter = {1} \n".format(corr_obj, obj_lit_cef))
    
corr_result.close()

del field_names
del corr_table

#building_field = corr_table[field_dic["building"]]
#tree_field = corr_table[field_dic["tree"]]
#light_field = corr_table[field_dic["light"]]

#building_lit_corr = np.corrcoef(building_field, litter_field)
#building_lit_cef = building_lit_corr[0][1]
#print("Correlation Coefficient for Building-Litter Index: {0}".format(building_lit_cef))

#tree_lit_corr = np.corrcoef(tree_field, litter_field)
#tree_lit_cef = tree_lit_corr[0][1]
#print("Correlation Coefficient for Tree-Litter Index: {0}".format(tree_lit_cef))

#light_lit_corr = np.corrcoef(light_field, litter_field)
#light_lit_cef = light_lit_corr[0][1]
#print("Correlation Coefficient for Light Pole-Litter Index: {0}".format(light_lit_cef))






# PROCESS: Map
# open mxd and data frames
mxd = arcpy.mapping.MapDocument(outdir + map_name)
# 1st dataframe
data_frames = arcpy.mapping.ListDataFrames(mxd)
df = data_frames[0]
# add grid to 1st data frame 1st layer
lyr_file = arcpy.mapping.Layer(grid)
arcpy.mapping.AddLayer(df, lyr_file)

# visualization

# input source .lyr, export filename, .lyr mannually create, natural break, excluded 0
for k in field_dic.keys():
    
    src_lyr = outdir + "/" + field_dic[k] + ".lyr"
    field = field_dic[k]
    export = arcpy.env.workspace + "/"
    
    visual(src_lyr, mxd, field, k, export)

    print("Finished maping {0}".format(k))

# save new file
mxd.saveACopy(arcpy.env.workspace + "/" + "new" + map_name)

del mxd
del data_frames
del lyr_file
del df



