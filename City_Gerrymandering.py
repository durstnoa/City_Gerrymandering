# -*- coding: utf-8 -*-f
"""
Created on Thu Feb 11 11:58:24 2021

@author: noahd
"""

# In[1]:

## loading packages
from shapely.ops import cascaded_union
import geopandas as gpd
#get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt

import pickle
import numpy as np
import math
import csv

from shapely.ops import polygonize
from shapely.geometry import Polygon
from shapely.geometry import LinearRing
from shapely.geometry import Point
import pandas as pd
from random import random
import os
import time
from time import process_time 


# In[3]:

# =============================================================================
# Loading data
# =============================================================================

# Setting working directory
os.chdir('D:\OneDrive - Michigan State University\Spatial Annexation')
os.getcwd()


# Getting city shapes for 2010
US_place_2010=gpd.read_file('D:/OneDrive - Michigan State University/Spatial Annexation/Automated Identification/shapefiles/US_place_2010.shp')
# Identifying cities, dropping census designated places
US_place_2010=US_place_2010.loc[(US_place_2010['CLASSFP10'] =='C1') | (US_place_2010['CLASSFP10'] == 'C2') | (US_place_2010['CLASSFP10'] == 'C3')| (US_place_2010['CLASSFP10'] == 'C4')| (US_place_2010['CLASSFP10'] == 'C5')| (US_place_2010['CLASSFP10'] == 'C6')| (US_place_2010['CLASSFP10'] == 'C7')| (US_place_2010['CLASSFP10'] == 'C8')| (US_place_2010['CLASSFP10'] == 'C9')]
US_place_2010_name= US_place_2010[['geometry', 'GISJOIN', 'NAME10', 'STATEFP10']].copy()
# saving only geometry and GISJOIN codes
US_place_2010 = US_place_2010[['geometry', 'GISJOIN']].copy()


# Getting city shapes for 2000
US_place_2000=gpd.read_file('D:/OneDrive - Michigan State University/Spatial Annexation/Automated Identification/shapefiles/US_place_2000.shp')
US_place_2000=US_place_2000.loc[(US_place_2000['FIPSCC'] =='C1') | (US_place_2000['FIPSCC'] == 'C2') | (US_place_2000['FIPSCC'] == 'C3')| (US_place_2000['FIPSCC'] == 'C4')| (US_place_2000['FIPSCC'] == 'C5')| (US_place_2000['FIPSCC'] == 'C6')| (US_place_2000['FIPSCC'] == 'C7')| (US_place_2000['FIPSCC'] == 'C8')| (US_place_2000['FIPSCC'] == 'C9')]
US_place_2000=US_place_2000[['GISJOIN', 'geometry']].copy()


Post_DF=US_place_2010.copy()
Pre_DF=US_place_2000.copy()

# Identifying cities that reported annexations to the BAS between 2000 and 2010
BAS_data=pd.read_csv('D:/OneDrive - Michigan State University/Spatial Annexation/Static analysis/BAS_2000_to_2010.csv')
Annexed_list=list(BAS_data['GISJOIN'].unique())

Post_DF_annexed=Post_DF[Post_DF['GISJOIN'].isin(Annexed_list)]


#%%

# =============================================================================
# # =============================================================================
# # # =============================================================================
# # # Defining primary functions
# # # =============================================================================
# # =============================================================================
# =============================================================================


## Identifying cities' longitude
def long2UTM(long):
    """Identifying cities' longitude for iterative reprojection of cities to appropriate UTM"""
    UTM=(((long + 180)/6)%60)+ 1
    return UTM

## Identifying cities' EPSG and UTM codes
def find_UTM(DF_preprocess):
    """Identifying cities' UTM projection"""
    DF_WGS84=DF_preprocess.copy()
    DF_WGS84=DF_WGS84.to_crs("EPSG:4326")
    exploded=DF_WGS84.reset_index(drop=True).explode()
    exploded = exploded.reset_index() 
    exploded=exploded[exploded['level_1']==0]
    x_list=[]
    for index, row in exploded.iterrows():
        #x=row[5].exterior.coords.xy[0][0]
        x=row['geometry'].exterior.coords.xy[0][0]
        x_list.append(x)
    utm_list=[]
    for i in x_list:
        utm=math.floor(long2UTM(i))
        if utm<10:
            utm = '0'+str(utm)
        epsg_code = 'epsg:326' + str(utm)
        utm_list.append([epsg_code, utm])
    utm_df=pd.DataFrame(utm_list)
    DF_WGS84=DF_WGS84.reset_index(drop=True)
    DF_WGS84.reset_index(drop=True)
    DF_WGS84_merged=DF_WGS84.merge(utm_df, left_index=True, right_index=True, how='outer')
    DF_WGS84_merged=DF_WGS84_merged.rename(columns={0: 'EPSG'})
    DF_WGS84_merged=DF_WGS84_merged.rename(columns={1: 'UTM'})
    return DF_WGS84_merged




state_codes = {
    '01':'AL',
    '02':'AK',
    '04':'AZ',
    '05':'AR',
    '06':'CA',
    '08':'CO',
    '09':'CT',
    '10':'DE',
    '11':'DC',
    '12':'FL',
    '13':'GA',
    '15':'HI',
    '16':'ID',
    '17':'IL',
    '18':'IN',
    '19':'IA',
    '20':'KS',
    '21':'KY',
    '22':'LA',
    '23':'ME',
    '24':'MD',
    '25':'MA',
    '26':'MI',
    '27':'MN',
    '28':'MS',
    '29':'MO',
    '30':'MT',
    '31':'NE',
    '32':'NV',
    '33':'NH',
    '34':'NJ',
    '35':'NM',
    '36':'NY',
    '37':'NC',
    '38':'ND',
    '39':'OH',
    '40':'OK',
    '41':'OR',
    '42':'PA',
    '44':'RI',
    '45':'SC',
    '46':'SD',
    '47':'TN',
    '48':'TX',
    '49':'UT',
    '50':'VT',
    '51':'VA',
    '53':'WA',
    '54':'WV',
    '55':'WI',
    '56':'WY'
}


def spatial_index(cities):
    """Creates spatial index to identify nearby cities"""
    cities_newIDX=cities.reset_index(drop=True)
    cities_sindex=cities_newIDX.sindex
    return cities_sindex, cities_newIDX


def find_shapes(city_test, boundary):
    """ Identifies most of the primary shapes of interest, returned as shapely objects"""
    GISJOIN=city_test['GISJOIN'].iloc[0]
    city=city_test['geometry'].iloc[0]
    city=city.buffer(0)
    box_min = city.minimum_rotated_rectangle
    x, y = box_min.exterior.coords.xy
    edge_length = (Point(x[0], y[0]).distance(Point(x[1], y[1])), Point(x[1], y[1]).distance(Point(x[2], y[2])))
    width=min(edge_length)/2
    city_buffer_reduce=city.buffer(width).buffer(-1.01*width)
    city_buffer_reduce=city_buffer_reduce.difference(boundary)
    city_buffer_temp=city.buffer(width/4)
    city_buffer_temp_2=city_buffer_temp.difference(boundary)
    city_buffer_reduce_2=city_buffer_reduce.buffer(1)
    city_buffer=city_buffer_temp_2.difference(city_buffer_reduce_2)
    
    city_buffer_temp_3=city_buffer_temp.difference(city)
    other_cities=city_buffer_temp_3.buffer(-1).intersection(boundary)
    return city, city_buffer, city_buffer_reduce, other_cities, GISJOIN

def nearby_cities(DF, cities_newIDX, cities_sindex, epsg):
    """Identifies nearby cities"""     
    
    city=DF['geometry'].iloc[0]
    city=city.buffer(0)
    box_min = city.minimum_rotated_rectangle
    x, y = box_min.exterior.coords.xy
    edge_length = (Point(x[0], y[0]).distance(Point(x[1], y[1])), Point(x[1], y[1]).distance(Point(x[2], y[2])))
    width=min(edge_length)/2
    DF['geometry']=DF.buffer(width)
    
    boxes=DF.geometry.bounds
    bounds=[]
    bounds.append(boxes['minx'].min())
    bounds.append(boxes['miny'].min())
    bounds.append(boxes['maxx'].max())
    bounds.append(boxes['maxy'].max())            
    place_candidate_idx = list(cities_sindex.intersection(bounds))
    place_candidates=cities_newIDX.loc[place_candidate_idx]
    place_candidates['Other_Place']=place_candidates[['GISJOIN']]
    place_candidates=place_candidates[['Other_Place', 'geometry']].copy()
    place_candidates=place_candidates.set_geometry('geometry')
    city_shapes=place_candidates['geometry']
    cities_buffer=city_shapes.buffer(.001)
    cities_buffer=cities_buffer.to_crs(epsg)
    boundary = gpd.GeoSeries(cascaded_union(cities_buffer))
    boundary=boundary.loc[0]
    return boundary

def underbounded_simple(city, city_buffer_reduce, boundary, GISJOIN, epsg):
    """" Identifies underbounded and unaffected areas (i.e., all unincorporated areas) 
         -- used only for cross-sectional analysis"""
    if boundary.type=='Polygon':
        cities_DF=gpd.GeoDataFrame([boundary])
    else:
        cities_DF=gpd.GeoDataFrame(boundary)
             
    if len(cities_DF)==0:
        cities_DF['geometry']=""
        cities_DF=cities_DF.set_geometry('geometry')
        cities_DF=cities_DF.set_crs(epsg)
    else:
        cities_DF=cities_DF.rename(columns={0:'geometry'})
        cities_DF=cities_DF.set_geometry('geometry')
        cities_DF=cities_DF.set_crs(epsg)
        cities_DF['GISJOIN']=GISJOIN
        cities_DF['geometry']=cities_DF.buffer(1)
    if city_buffer_reduce.type=='Polygon':
        city_buffer_reduce_DF=gpd.GeoDataFrame([city_buffer_reduce])
    else:
        city_buffer_reduce_DF=gpd.GeoDataFrame(city_buffer_reduce)
    
    city_buffer_reduce_DF=city_buffer_reduce_DF.rename(columns={0:'geometry'})
    city_buffer_reduce_DF=city_buffer_reduce_DF.set_geometry('geometry')
    city_buffer_reduce_DF=city_buffer_reduce_DF.set_crs(epsg)
    city_buffer_reduce_DF=city_buffer_reduce_DF.explode()
    city_buffer_reduce_DF=city_buffer_reduce_DF.reset_index()
    city_buffer_reduce_DF['index'] = city_buffer_reduce_DF.index
    city_buffer_reduce_DF['AreaAcres']=city_buffer_reduce_DF.area*0.000247105
    too_small=city_buffer_reduce_DF[city_buffer_reduce_DF['AreaAcres']<1].copy()
    too_small['Type']='Unaffected'
    
    city_buffer_reduce_DF=city_buffer_reduce_DF[city_buffer_reduce_DF['AreaAcres']>=1].copy()
    city_buffer_reduce_DF['Type']="Underbounded"


    results=pd.concat([city_buffer_reduce_DF, too_small])
    results['GISJOIN']=GISJOIN
    results=results[['GISJOIN', 'geometry', 'AreaAcres', 'Type']].copy()
    return  results

def cities(city, other_cities,city_buffer, Pre_DF, GISJOIN, epsg):
    """ Identifies cities and unincorporated areas blocked by other cities
        -- used for analysis of boundary changes over time"""

    city_DF=Pre_DF[Pre_DF['GISJOIN']==GISJOIN].copy()

    city_DF=city_DF.to_crs(epsg)
    city_DF['GISJOIN']=GISJOIN
    city_DF['Type']='City' 
    if other_cities.type=='Polygon':
        other_cities_DF=gpd.GeoDataFrame([other_cities])
    else:
        other_cities_DF=gpd.GeoDataFrame(other_cities)
    
    if len(other_cities_DF)==0:
        other_cities_DF['geometry']=""
    else:
        other_cities_DF=other_cities_DF.rename(columns={0:'geometry'})
        other_cities_DF=other_cities_DF.set_geometry('geometry')
        other_cities_DF=other_cities_DF.set_crs(epsg)
        other_cities_DF['GISJOIN']=GISJOIN
        other_cities_DF['Type']='Other_cities'        
        
    if city_buffer.type=='Polygon':
        city_buffer_DF=gpd.GeoDataFrame([city_buffer])
    else:
        city_buffer_DF=gpd.GeoDataFrame(city_buffer)
    
    if len(city_buffer_DF)==0:
        city_buffer_DF['geometry']=""
    else:
        city_buffer_DF=city_buffer_DF.rename(columns={0:'geometry'})
        city_buffer_DF=city_buffer_DF.set_geometry('geometry')
        city_buffer_DF=city_buffer_DF.set_crs(epsg)
        city_buffer_DF['GISJOIN']=GISJOIN
        city_buffer_DF['Type']='Unaffected' 
        
        # Identifying areas that don't touch the city
        
        city_buffer_DF_temp=city_buffer_DF.copy()
        city_buffer_DF_temp['join_index']=city_buffer_DF_temp.index
        
        if city.type=='Polygon':
            city_DF_temp=gpd.GeoDataFrame([city])
        else:
            city_DF_temp=gpd.GeoDataFrame(city)
        
        if len(city_DF)==0:
            city_DF_temp['geometry']=""
        else:
            city_DF_temp=city_DF_temp.rename(columns={0:'geometry'})
            city_DF_temp=city_DF_temp.set_geometry('geometry')
            city_DF_temp=city_DF_temp.set_crs(epsg)

        city_DF_temp['geometry']=city_DF_temp.buffer(5)
        city_buffer_DF_temp= gpd.overlay(city_buffer_DF_temp, city_DF_temp,  how='intersection', make_valid=True)
        city_buffer_DF_temp['Not_Blocked']='Not_Blocked'
        
        city_buffer_DF_temp=city_buffer_DF_temp[['Not_Blocked','join_index']].copy()
        
        city_buffer_DF=pd.merge(city_buffer_DF, city_buffer_DF_temp, left_index=True, right_on='join_index', how='outer')
        city_buffer_DF['Type']= np.where(city_buffer_DF['Not_Blocked']!='Not_Blocked', 'Blocked', city_buffer_DF['Type'])
        
    Cities=pd.concat([city_DF, other_cities_DF, city_buffer_DF])
    Cities=Cities.reset_index()
    Cities=Cities[['GISJOIN', 'geometry', 'Type']].copy()
    return Cities               

def cities_simple(city, other_cities,city_buffer, Pre_DF, GISJOIN, epsg):
    """ Identifies cities and unincorporated areas blocked by other cities
        -- used only for cross-sectional analysis"""

   
    if city.type=='Polygon':
        city_DF=gpd.GeoDataFrame([city])
    else:
        city_DF=gpd.GeoDataFrame(city)
    
    if len(city_DF)==0:
        city_DF['geometry']=""
    else:
        city_DF=city_DF.rename(columns={0:'geometry'})
        city_DF=city_DF.set_geometry('geometry')
        city_DF=city_DF.set_crs(epsg)
        city_DF['GISJOIN']=GISJOIN
        city_DF['Type']='City' 

    
    if other_cities.type=='Polygon':
        other_cities_DF=gpd.GeoDataFrame([other_cities])
    else:
        other_cities_DF=gpd.GeoDataFrame(other_cities)
    
    if len(other_cities_DF)==0:
        other_cities_DF['geometry']=""
    else:
        other_cities_DF=other_cities_DF.rename(columns={0:'geometry'})
        other_cities_DF=other_cities_DF.set_geometry('geometry')
        other_cities_DF=other_cities_DF.set_crs(epsg)
        other_cities_DF['GISJOIN']=GISJOIN
        other_cities_DF['Type']='Other_cities' 
        

    if city_buffer.type=='Polygon':
        city_buffer_DF=gpd.GeoDataFrame([city_buffer])
    else:
        city_buffer_DF=gpd.GeoDataFrame(city_buffer)
    
    if len(city_buffer_DF)==0:
        city_buffer_DF['geometry']=""
    else:
        city_buffer_DF=city_buffer_DF.rename(columns={0:'geometry'})
        city_buffer_DF=city_buffer_DF.set_geometry('geometry')
        city_buffer_DF=city_buffer_DF.set_crs(epsg)
        city_buffer_DF['GISJOIN']=GISJOIN
        city_buffer_DF['Type']='Unaffected' 
        
        # Identifying areas that don't touch the city
        
        city_buffer_DF_temp=city_buffer_DF.copy()
        city_buffer_DF_temp['join_index']=city_buffer_DF_temp.index
        
        city_DF_temp=city_DF.copy()
        city_DF_temp['geometry']=city_DF_temp.buffer(5)
        city_buffer_DF_temp= gpd.overlay(city_buffer_DF_temp, city_DF_temp,  how='intersection', make_valid=True)
        city_buffer_DF_temp['Not_Blocked']='Not_Blocked'
        
        city_buffer_DF_temp=city_buffer_DF_temp[['Not_Blocked','join_index']].copy()
        
        city_buffer_DF=pd.merge(city_buffer_DF, city_buffer_DF_temp, left_index=True, right_on='join_index', how='outer')
        city_buffer_DF['Type']= np.where(city_buffer_DF['Not_Blocked']!='Not_Blocked', 'Blocked', city_buffer_DF['Type'])
        
    Cities=pd.concat([city_DF, other_cities_DF, city_buffer_DF])
    Cities=Cities[['GISJOIN', 'geometry', 'Type']].copy()
    return Cities               

def annexed(DF, Pre_DF, epsg):
    """Takes Pre and Post-annexation city shapefile and classifies areas that were annexed"""
    Post=DF.iloc[[i]]
    
    GISJOIN=Post['GISJOIN'].iloc[0]
    Pre=Pre_DF[Pre_DF['GISJOIN']==GISJOIN].copy()
    Pre=Pre.to_crs(epsg)
    if len(Pre)>0:
        sym_diff = gpd.overlay(Pre, Post, how='symmetric_difference')
        
        sym_diff['geometry']=sym_diff.buffer(0)
        sym_diff['AreaAcres'] = sym_diff.area*0.000247105
        sym_diff=sym_diff[sym_diff['AreaAcres']>1].copy()
        Pre['geometry']=Pre.buffer(1)
        intersect = gpd.overlay(sym_diff, Pre, how='intersection')
        
        if len(intersect) == 0:
            diff = sym_diff.explode()
        else:
            intersect['geometry']=intersect.buffer(2)
            diff = gpd.overlay(sym_diff, intersect, how='difference').explode()
    
        
        diff['AreaAcres'] = diff.area*0.000247105
        diff = diff.reset_index(level=[0])
        diff = diff.rename(columns={'GISJOIN_1': 'GISJOIN'})
        diff = diff[['GISJOIN', 'geometry', 'AreaAcres']].copy()

#        diff=diff[diff['AreaAcres']>1].copy()
        diff=diff.reset_index(drop=True)
        diff['join_index'] = diff.index
        bounds = diff.copy()
        bounds['BoundaryGeometry'] = bounds.boundary
        bounds = bounds[['GISJOIN', 'AreaAcres', 'join_index', 'BoundaryGeometry']].copy()
        bounds = bounds.rename(columns={'BoundaryGeometry': 'geometry'})
        bounds = bounds.set_geometry('geometry')
         
    
        bounds['A']=bounds.length
        Pre_buffer=Pre.copy()
        Pre_buffer['geometry']=Pre_buffer.buffer(5)
    
        interior_line = gpd.overlay(bounds, Pre_buffer, how='intersection', make_valid=True, keep_geom_type=False)
    
        interior_line=interior_line.dissolve(by=['join_index'])
        interior_line=interior_line.reset_index()
        interior_line['A']= np.where(interior_line['A'].isnull(), 0, interior_line['A'])
        
        interior_line['C_pre'] = interior_line.length
        interior_line['C_pre']= np.where(interior_line['C_pre'].isnull(), 0, interior_line['C_pre'])
    
    
        interior_line['Shape_index'] = (interior_line['A'] - interior_line['C_pre']) /interior_line['A'] 
    
        interior_line = interior_line[['Shape_index', 'C_pre', 'geometry', 'join_index', 'A']].copy()
        interior_line['Geom_Type'] = interior_line['geometry'].geom_type
        interior_line = interior_line[interior_line.Geom_Type != 'Point']
        interior_line = interior_line[interior_line.Geom_Type != 'MultiPoint']
        
        poly_line = diff.merge(interior_line, on='join_index', how='left')
        poly_line['Shape_index']= np.where(poly_line['Shape_index'].isnull(), 1, poly_line['Shape_index'])
    
        poly_line.rename(columns={'geometry_x': 'geometry'}, inplace=True)
        poly_line = poly_line.set_geometry('geometry')
        poly_line.drop(['geometry_y'], axis=1, inplace=True)
        poly_line=poly_line.reset_index()
        
        poly_line['Type']="Annexed"
        poly_line['GISJOIN']=GISJOIN
        poly_line=poly_line[['GISJOIN', 'geometry', 'AreaAcres', 'Shape_index', 'Type', 'A', 'C_pre']].copy()
       
    else:
        poly_line=pd.DataFrame([])
    return  poly_line



def underbounded(city, city_buffer_reduce, boundary, GISJOIN, epsg):
    if boundary.type=='Polygon':
        cities_DF=gpd.GeoDataFrame([boundary])
    else:
        cities_DF=gpd.GeoDataFrame(boundary)
             
    if len(cities_DF)==0:
        cities_DF['geometry']=""
        cities_DF=cities_DF.set_geometry('geometry')
        cities_DF=cities_DF.set_crs(epsg)
    else:
        cities_DF=cities_DF.rename(columns={0:'geometry'})
        cities_DF=cities_DF.set_geometry('geometry')
        cities_DF=cities_DF.set_crs(epsg)
        cities_DF['GISJOIN']=GISJOIN
        cities_DF['geometry']=cities_DF.buffer(1)
    if city_buffer_reduce.type=='Polygon':
        city_buffer_reduce_DF=gpd.GeoDataFrame([city_buffer_reduce])
    else:
        city_buffer_reduce_DF=gpd.GeoDataFrame(city_buffer_reduce)
    
    city_buffer_reduce_DF=city_buffer_reduce_DF.rename(columns={0:'geometry'})
    city_buffer_reduce_DF=city_buffer_reduce_DF.set_geometry('geometry')
    city_buffer_reduce_DF=city_buffer_reduce_DF.set_crs(epsg)
    city_buffer_reduce_DF=city_buffer_reduce_DF.explode()
    city_buffer_reduce_DF=city_buffer_reduce_DF.reset_index()
    city_buffer_reduce_DF['join_index'] = city_buffer_reduce_DF.index
    city_buffer_reduce_DF['AreaAcres']=city_buffer_reduce_DF.area*0.000247105
    too_small=city_buffer_reduce_DF[city_buffer_reduce_DF['AreaAcres']<1].copy()
    too_small['Type']='Unaffected'
    
    city_buffer_reduce_DF=city_buffer_reduce_DF[city_buffer_reduce_DF['AreaAcres']>=1].copy()
    
    bounds=city_buffer_reduce_DF.copy()
    bounds['BoundaryGeometry'] = city_buffer_reduce_DF.boundary   
    bounds = bounds[['AreaAcres', 'BoundaryGeometry', 'join_index']].reset_index(drop=True)
    bounds = bounds.rename(columns={'BoundaryGeometry': 'geometry'})
    bounds = bounds.set_geometry('geometry')
    bounds['U']=bounds.length
    
    interior_line = gpd.overlay(bounds, cities_DF, how='intersection', make_valid=True, keep_geom_type=False)
    interior_line=interior_line.dissolve(by=['join_index'])
    interior_line=interior_line.reset_index()
    interior_line['C_post'] = interior_line.length
    interior_line['C_post']= np.where(interior_line['C_post'].isnull(), 0, interior_line['C_post'])

    interior_line.loc[interior_line['C_post'].isnull(), 'C_post'] = 0
    interior_line['Shape_index'] = (interior_line['C_post'] /(interior_line['U']*.5))-1 
    interior_line['Shape_index']= np.where(interior_line['Shape_index']<0, 0, interior_line['Shape_index'])
    interior_line = interior_line[['Shape_index', 'U', 'C_post', 'geometry', 'join_index']].copy()
    interior_line['Geom_Type'] = interior_line['geometry'].geom_type
    interior_line = interior_line[interior_line.Geom_Type != 'Point']
    interior_line = interior_line[interior_line.Geom_Type != 'MultiPoint']
    
    poly_line = city_buffer_reduce_DF.merge(interior_line, on='join_index', how='left')
    
    poly_line.rename(columns={'geometry_x': 'geometry'}, inplace=True)
    poly_line = poly_line.set_geometry('geometry')
    poly_line.drop(['geometry_y'], axis=1, inplace=True)
    
    poly_line['Type']="Underbounded"
    
    if city.type=='Polygon':
        city_DF_bounds=gpd.GeoDataFrame([city])
    else:
        city_DF_bounds=gpd.GeoDataFrame(city)
    
    city_DF_bounds=city_DF_bounds.rename(columns={0:'geometry'})
    city_DF_bounds=city_DF_bounds.set_geometry('geometry')
    city_DF_bounds=city_DF_bounds.set_crs(epsg)
    city_DF_bounds['geometry']=city_DF_bounds.buffer(5)
    city_DF_bounds['geometry']=city_DF_bounds.boundary
    
    if len(poly_line)==0:
        poly_line=poly_line.copy()    
    elif len(poly_line)==1 & poly_line['geometry'].iloc[0].is_empty: 
        poly_line=poly_line.copy()    
    else:
        test_intersect= gpd.overlay(city_DF_bounds, poly_line,  how='intersection', make_valid=True, keep_geom_type=False)
        test_intersect=test_intersect.dissolve(by='join_index')
        test_intersect=test_intersect.reset_index()
        test_intersect['Length_Touch']=test_intersect.length
        test_intersect=test_intersect[['join_index', 'Length_Touch']].copy()
        poly_line=pd.merge(poly_line, test_intersect, left_on='join_index', right_on='join_index', how='outer')
        poly_line['Length_Touch'] = poly_line['Length_Touch'].fillna(0)
        poly_line['Type']= np.where(poly_line['Length_Touch']==0, 'Blocked', poly_line['Type'])
        poly_line['Type']= np.where(poly_line['Length_Touch']==np.nan, 'Blocked', poly_line['Type'])
    
    poly_line=pd.concat([poly_line, too_small])
    poly_line=poly_line.reset_index()
    poly_line['GISJOIN']=GISJOIN
    poly_line=poly_line[['GISJOIN', 'geometry', 'AreaAcres', 'Shape_index', 'Type','U', 'C_post']].copy()
    return poly_line
#%%



#%%
# =============================================================================
# # =============================================================================
# # # =============================================================================
# # # Processing cities for 2000 to 2010 analysis
# # # =============================================================================
# # =============================================================================
# =============================================================================

# =============================================================================
# Reprojecting city boundaries to appropriate UTMs; storing as listed dataframes
# =============================================================================

#np.random.seed(123)
#Post_DF_annexed_selected=Post_DF_annexed.sample(1000)


Post_UTM=find_UTM(Post_DF_annexed)
#Post_UTM=find_UTM(Post_DF_annexed_selected)

print(Post_UTM.UTM.value_counts())

DF_list_2=[]
UTM_list=list(Post_UTM.UTM.unique())
UTM_list = [int(i) for i in UTM_list] 
UTM_list=[10, 11, 12, 13, 14, 15, 16, 17, 18]

for UTM in UTM_list:
    exec("Post_UTM_{utm}=Post_UTM[Post_UTM['UTM']=={utm}]".format(utm=UTM))
    exec("DF_list_2.append(Post_UTM_{utm})".format(utm=UTM))


Error_list=[]
error=0
k=0
City_Stats=[]

cities_sindex, cities_newIDX=spatial_index(Post_DF) 
sindex_crs=cities_newIDX.crs
for DF in DF_list_2:
    t1_start = process_time()  
    
    print("DF number", k)
    
    epsg=DF['EPSG'].iloc[0]
    epsg_2=DF.crs
    DF=DF.to_crs(epsg)

    cities_newIDX_reproject=cities_newIDX.to_crs(epsg)
    print("finished data preparation")
    GDF_list=[]
# =============================================================================
#     Looping over cities, one at a time to identify shapes
# =============================================================================
    if len(DF)>0:
        for i in range(len(DF)):
            #print(DF['GISJOIN'].iloc[[i]])

            GISJOIN_test=DF['GISJOIN'].iloc[[i]].iloc[0]
            Pre_test=Pre_DF[Pre_DF['GISJOIN']==GISJOIN_test].copy()
            
            city_test=DF.iloc[[i]].copy()
            city_test_project=city_test.to_crs(sindex_crs)
            boundary=nearby_cities(city_test_project, cities_newIDX_reproject, cities_sindex, epsg)
            
            if (i+1)%50==0:
                print('Processing metrics for %sth city'%(i+1))
            if len(Pre_test)>0:
                try:
                    city, city_buffer, city_buffer_reduce, other_cities, GISJOIN= find_shapes(city_test, boundary)
                    #print(GISJOIN)
                    Cities=cities(city, other_cities, city_buffer, Pre_DF, GISJOIN, epsg )
                    Cities=Cities.to_crs(epsg_2)
                    Unincorporated = underbounded(city, city_buffer_reduce, boundary, GISJOIN, epsg)
                    Unincorporated=Unincorporated.to_crs(epsg_2)
                    Annexed=annexed(DF, Pre_DF, epsg)
                    Annexed=Annexed.to_crs(epsg_2)
                    GDF_list.append(Unincorporated)
                    GDF_list.append(Cities)
                    GDF_list.append(Annexed)

                except Exception as e:
                      Error_list.append([GISJOIN, e])
                      error+=1
                      print("Error number", error)
                      pass
    k+=1
    t1_stop = process_time()  
    print("processed ", (i+1), " cities in ", (t1_stop-t1_start)/60, "minutes")
    Cities_per_minute=(i+1)/((t1_stop-t1_start)/60)
    print("processed ", Cities_per_minute, " cities per minute ")
# =============================================================================
#     Combining and saving results
# =============================================================================

    GDF=gpd.GeoDataFrame(pd.concat(GDF_list, ignore_index=True), crs=GDF_list[0].crs)
    # Identifying cases where the geography is empty
    
    GDF['Empty_Geo']=GDF.is_empty
    
    GDF=GDF[GDF['Empty_Geo']!=True].copy()
    UTM=epsg[-2:]
    exec('GDF[GDF["Empty_Geo"]!=True].to_pickle("Type_Results_2000_to_2010_UTM_{utm}.pkl")'.format(utm=UTM))
    


GDF_list_4=[]
for UTM in UTM_list:    
    exec('GDF=pd.read_pickle("D:\OneDrive - Michigan State University\Spatial Annexation\Type_Results_2000_to_2010_UTM_{utm}.pkl")'.format(utm=UTM))
    GDF_list_4.append(GDF)

GDF=gpd.GeoDataFrame(pd.concat(GDF_list_4, ignore_index=True), crs=GDF_list[0].crs)

GDF=GDF.drop_duplicates(keep='last')

GDF=GDF[(GDF['AreaAcres']>1) |(GDF['Type']!="Annexed")].copy()

GDF[GDF["Empty_Geo"]!=True].to_pickle("Type_Results_2000_to_2010_Final.pkl")

GDF_export=GDF[['GISJOIN', 'AreaAcres', 'Shape_index', 'Type']].copy()
GDF_export.to_csv("Type_Results_2000_to_2010_Final.csv")

#%%

GDF = pd.read_pickle("D:\OneDrive - Michigan State University\Spatial Annexation\Type_Results_2000_to_2010_Final.pkl")

#%%


# =============================================================================
# # =============================================================================
# # Dissolving for city-level statistics for 2000 to 2010 temporal analysis
# # =============================================================================
# =============================================================================



national_crs=US_place_2010_name.crs
US_place_2010_name=US_place_2010_name[['NAME10', 'GISJOIN']].copy()

GDF_dissolve=GDF.copy()
GDF_dissolve=GDF_dissolve.to_crs(national_crs)
GDF_dissolve=GDF_dissolve.dissolve(by=['Type', 'GISJOIN'])
GDF_dissolve=GDF_dissolve.reset_index()


GDF_dissolve['AreaAcres']=GDF_dissolve.area*0.000247105

GDF_dissolve_1=GDF_dissolve.pivot_table(index="GISJOIN", columns='Type', values='AreaAcres')
GDF_dissolve_1=GDF_dissolve_1.fillna(0)

column_list=(GDF_dissolve_1.columns)
GDF_dissolve_1['Denominator']= GDF_dissolve_1['City']
for column in column_list:
    GDF_dissolve_1[column]=GDF_dissolve_1[column]/GDF_dissolve_1['Denominator']

GDF_dissolve_1=pd.merge(GDF_dissolve_1, US_place_2010_name, left_index=True, right_on='GISJOIN', how='inner')
GDF_dissolve_2=GDF[['GISJOIN', 'Shape_index', 'Type']].copy()
GDF_dissolve_2=GDF_dissolve_2[(GDF_dissolve_2['Type']=='Underbounded') | (GDF_dissolve_2['Type']=='Annexed')].copy()
GDF_dissolve_2['Shape_index'] = GDF_dissolve_2['Shape_index'].astype(float)
GDF_dissolve_2=GDF_dissolve_2.pivot_table(index="GISJOIN", columns='Type', values='Shape_index', aggfunc='mean')
GDF_dissolve_2=GDF_dissolve_2.rename(columns={'Annexed':'Shape_Annex', 'Underbounded':'Shape_Under'})
GDF_dissolve_3=pd.merge(GDF_dissolve_1, GDF_dissolve_2, left_on='GISJOIN', right_on='GISJOIN', how='inner')

# Shape_Annex and Shape_Under are the average gerrymandering index for each city.

# Other variables (Annexed, etc.) are the ratio of the type's area to the area of the city itself

GDF_dissolve_3.to_pickle("Type_Results_2000_to_2010_Dissolve.pkl")
GDF_dissolve_3.to_csv("Type_Results_2000_to_2010_Dissolve.csv")
GDF_dissolve=GDF_dissolve_3.copy()

#%%
GDF_dissolve = pd.read_pickle("D:\OneDrive - Michigan State University\Spatial Annexation\Type_Results_2000_to_2010_Dissolve.pkl")
#%%
# Descriptive statistics
print("Dissolved statistics - city level")
print(GDF_dissolve.mean())
print("Number of areas by type")
print(GDF['Type'].value_counts())
print("Number of unincorporated islands and satellite annexations")
print(GDF[GDF['Shape_index']==1]['Type'].value_counts())
GDF['Shape_index'] = GDF['Shape_index'].astype(float)
print("Average gerrymandering index by type")
print(GDF.groupby('Type')['Shape_index'].mean().round(3))
GDF=GDF.to_crs(national_crs)
GDF['Area_Sq_Mi']=GDF.area*0.000000386102159
print("Size of area in square miles by type")
print(GDF.groupby('Type')['Area_Sq_Mi'].sum().round(0))

#%%
# =============================================================================
# # =============================================================================
# # Processing Houston for map
# # =============================================================================
# =============================================================================


#Post_DF_Houston=Post_DF[Post_DF['GISJOIN']=='G35002000'].copy()
Post_DF_Houston=Post_DF[Post_DF['GISJOIN']=='G48035000'].copy()
#Post_DF_Houston=Post_DF[Post_DF['GISJOIN']=='G20071125'].copy()
Post_DF_Houston_UTM=find_UTM(Post_DF_Houston)


DF_list_3=[]
UTM_list=list(Post_DF_Houston_UTM.UTM.unique())
UTM_list = [int(i) for i in UTM_list] 
UTM_list=[15]
#UTM_list=[14]
for UTM in UTM_list:
    exec("Post_DF_Houston_UTM_{utm}=Post_DF_Houston_UTM[Post_DF_Houston_UTM['UTM']=={utm}]".format(utm=UTM))
    exec("DF_list_3.append(Post_DF_Houston_UTM_{utm})".format(utm=UTM))
    
   
GDF_list=[]

Error_list=[]
error=0
k=0
City_Stats=[]

cities_sindex, cities_newIDX=spatial_index(Post_DF) 
sindex_crs=cities_newIDX.crs
for DF in DF_list_3:
    t1_start = process_time()  

    print("DF number", k)
    
    epsg=DF['EPSG'].iloc[0]
    epsg_2=DF.crs
    DF=DF.to_crs(epsg)

    cities_newIDX_reproject=cities_newIDX.to_crs(epsg)
    print("finished data preparation")


    if len(DF)>0:
        for i in range(len(DF)):
            
            city_test=DF.iloc[[i]].copy()
            city_test_project=city_test.to_crs(sindex_crs)
            boundary=nearby_cities(city_test_project, cities_newIDX_reproject, cities_sindex, epsg)
            if (i+1)%10==0:
                print('Processing metrics for %sth city'%(i+1))
            try:
                city, city_buffer, city_buffer_reduce, other_cities, GISJOIN= find_shapes(city_test, boundary)
                #print(GISJOIN)
                Cities=cities(city, other_cities, city_buffer, Pre_DF, GISJOIN, epsg )
                Cities=Cities.to_crs(epsg_2)
                Unincorporated = underbounded(city, city_buffer_reduce, boundary, GISJOIN, epsg)
                Unincorporated=Unincorporated.to_crs(epsg_2)

                Annexed=annexed(DF, Pre_DF, epsg)
                Annexed=Annexed.to_crs(epsg_2)
                GDF_list.append(Unincorporated)
                GDF_list.append(Cities)
                GDF_list.append(Annexed)

            except Exception as e:
                  Error_list.append([GISJOIN, e])
                  error+=1
                  print("Error number", error)
                  pass
    k+=1
    t1_stop = process_time()  
    print("processed ", (i+1), " cities in ", (t1_stop-t1_start)/60, "minutes")
    Cities_per_minute=(i+1)/((t1_stop-t1_start)/60)
    print("processed ", Cities_per_minute, " cities per minute ")


GDF_Houston=gpd.GeoDataFrame(pd.concat(GDF_list, ignore_index=True), crs=GDF_list[0].crs)

# Identifying cases where the geography is empty

GDF_Houston['Empty_Geo']=GDF_Houston.is_empty

GDF_Houston=GDF_Houston[GDF_Houston['Empty_Geo']!=True].copy()


GDF_Houston[GDF_Houston['Empty_Geo']!=True].to_pickle("Type_Results_2000_to_2010_Houston.pkl")


#%%

# =============================================================================
# Drawing map of Houston area with various shapes
# =============================================================================

GDF_Houston = pd.read_pickle("D:\OneDrive - Michigan State University\Spatial Annexation\Type_Results_2000_to_2010_Houston.pkl")

shape_list=[]
#['Underbounded', 'Blocked', 'Unaffected', 'City', 'Other_cities', 'Annexed']
color_list=['lightblue', 'indigo','royalblue' ,'white', 'cadetblue',  'darkviolet']

Houston_Type_list=list(GDF_Houston['Type'].unique())
for count, types in enumerate(Houston_Type_list):
    exec("Houston_{types}=GDF_Houston[GDF_Houston['Type']=='{types}'].copy()".format(types=types))
    exec("Houston_{types}['Color']=color_list[count]".format(types=types))
    exec("shape_list.append(Houston_{types})".format(types=types))

Houston_shapes=pd.concat(shape_list)
Houston_shapes=Houston_shapes.dissolve(by='Type')
Houston_shapes=Houston_shapes.reset_index()


first_plot=Houston_shapes.plot(color=Houston_shapes['Color'], legend=True, figsize=(12,12))
Houston_shapes[Houston_shapes['Type']=='City'].plot(color='white',ec='black', linewidth=1.5,
                    legend=True, figsize=(12,12), ax=first_plot)
import matplotlib.patches as mpatches

Underbounded_patch = mpatches.Patch(color='lightblue', label='Underbounded Areas')
Blocked_patch = mpatches.Patch(color='indigo', label='Blocked Areas')
Unaffected_patch = mpatches.Patch(color='royalblue', label='Unaffected Areas')
City_patch = mpatches.Patch(color='white',ec='black',  label='City of Houston, TX')
Other_cities_patch = mpatches.Patch(color='cadetblue', label='Other Cities')
Annexed_patch = mpatches.Patch(color='darkviolet', label='Annexed Areas')
plt.axis('off')
plt.legend(handles=[City_patch,Underbounded_patch,  Other_cities_patch,  Unaffected_patch, Blocked_patch, Annexed_patch], 
               bbox_to_anchor=(.85, -.05), ncol=2, prop={'size': 15}, 
               fancybox=True,  borderpad=1)
plt.savefig('Figure_1_B.png', dpi=600)
plt.show()
#%%


# =============================================================================
# Illustrating the creation of the concave hull for Houston
# =============================================================================

    
city_DF=gpd.GeoDataFrame(city)
city_DF=city_DF.rename(columns={0: 'geometry'})

box_min = city.minimum_rotated_rectangle
x, y = box_min.exterior.coords.xy
edge_length = (Point(x[0], y[0]).distance(Point(x[1], y[1])), Point(x[1], y[1]).distance(Point(x[2], y[2])))
width=min(edge_length)/2

city_buffer=city.buffer(width)
city_buffer_reduce=city.buffer(width).buffer(-1.01*width)


city_buffer_DF=gpd.GeoDataFrame([city_buffer])
city_buffer_DF=city_buffer_DF.rename(columns={0: 'geometry'})

city_buffer_reduce_DF=gpd.GeoDataFrame([city_buffer_reduce])
city_buffer_reduce_DF=city_buffer_reduce_DF.rename(columns={0: 'geometry'})

part_list=[city_DF, city_buffer_DF, city_buffer_reduce_DF ]

def draw_maps(part_list, n_rows, n_cols):
    fig=plt.figure(figsize=(7,7))
    titles=['City of Houston, TX', 'City Buffer (.5W)', "Concave Hull"]
    axis_range=part_list[1]
    axis_range=axis_range.set_crs(epsg)
    axis_range=axis_range['geometry']
    box=axis_range.bounds

    for i, part in enumerate(part_list):
        ax=fig.add_subplot(n_rows,n_cols,i+1)
        part=part.set_crs(epsg)
        part.plot(ax=ax, color='grey')
        ax.set_title(titles[i])
        ax.axis('off')
        plt.xlim(box['minx'].min(), box['maxx'].max())
        plt.ylim(box['miny'].min(), box['maxy'].max())
    
    fig.tight_layout()  # Improves appearance a bit.
    plt.savefig('Figure_1_A.png', dpi=600)
    plt.show()
    

draw_maps(part_list, 1, 3)

#%%


# =============================================================================
# # =============================================================================
# # # =============================================================================
# # # Joining city shapes from 2000 to 2010 temporal analysis to census blocks 
# # # =============================================================================
# # =============================================================================
# =============================================================================

Post_UTM_2=Post_UTM[['GISJOIN', 'UTM']].copy()

Selected_Areas=GDF[(GDF['Type']!='Blocked') & (GDF['Type']!='Other_cities') & (GDF['Type']!='City')].copy()
Selected_Areas['StateFIPS'] = Selected_Areas['GISJOIN'].str[1:3]

State_list=[]
state_list=list(Selected_Areas['StateFIPS'].unique())
for state in state_list:
    State=state_codes[state]
    State_list.append(State)

d_inverted = {v: k for k, v in state_codes.items()}
#state_list=list(state_codes.values())
#state_list_selected=[state_list[0]]+state_list[2:11]+state_list[12:]

start = time.time()
DF_list=[]
error_list=[]

# =============================================================================
# Looping over states
# =============================================================================

for State in State_list:
    print("Processing data for state: ", State)
    start_city_time = time.time()
    
    state_cities=Selected_Areas[Selected_Areas['StateFIPS']==d_inverted[State]].copy()
    
# =============================================================================
#     Getting census blocks
# =============================================================================
    
    exec('Blocks_Shapes_Data= gpd.read_file(r"D:/OneDrive - Michigan State University/Spatial Annexation/Static analysis/State blocks clean/{}_Blocks_Clean.shp")'.format(State))
    temp_crs=Blocks_Shapes_Data.crs
    state_cities=state_cities.to_crs(temp_crs)
    Blocks_Shapes_Data['Area_Pre']=Blocks_Shapes_Data.area
    Blocks_Shapes_Data['GISJOIN_Block']=Blocks_Shapes_Data['GISJOIN']    
    Blocks_Shapes_Data['geometry']=Blocks_Shapes_Data.buffer(0)
    state_cities['geometry']=state_cities.buffer(-.01) ## reducing by buffer of 1 centimeter to prevent non-noded intersections
    state_cities=state_cities.rename(columns={'GISJOIN': 'GISJOIN_City'})
    state_cities['Area_Part']=state_cities.area
# =============================================================================
#     Finding census blocks that intersect cities
# =============================================================================
    blocks_intersect=gpd.overlay(state_cities, Blocks_Shapes_Data, how='intersection')
    column_list=list(Blocks_Shapes_Data.columns)[2:-3]
    blocks_intersect['Area_Post']=blocks_intersect.area
    blocks_intersect['Area_Pct']=blocks_intersect['Area_Post']/blocks_intersect['Area_Pre']
    blocks_intersect=blocks_intersect[['GISJOIN_City','GISJOIN_Block', 'Type','Shape_index', 'Area_Part','Area_Pre', 'Area_Post', 'Area_Pct', 'geometry']+column_list].copy()
    blocks_intersect_temp=blocks_intersect[blocks_intersect['TotPop']>0].copy()
    blocks_intersect_temp=blocks_intersect_temp[blocks_intersect_temp['Area_Pct']>.99].copy()
    
    
# =============================================================================
#     Weighting population estimates by % of block that intersects with the city
# =============================================================================
    temp_array=blocks_intersect_temp['Area_Pct']*blocks_intersect_temp['TotPop']
    blocks_intersect_temp['TotPop']=temp_array
    temp_array=blocks_intersect_temp['Area_Pct']*blocks_intersect_temp['Urb_Area']
    blocks_intersect_temp['Urb_Area']=temp_array
    temp_array=blocks_intersect_temp['Area_Pct']*blocks_intersect_temp['Urb_Clust']
    blocks_intersect_temp['Urb_Clust']=temp_array
    temp_array=blocks_intersect_temp['Area_Pct']*blocks_intersect_temp['White']
    blocks_intersect_temp['White']=temp_array
    temp_array=blocks_intersect_temp['Area_Pct']*blocks_intersect_temp['Black']
    blocks_intersect_temp['Black']=temp_array
    temp_array=blocks_intersect_temp['Area_Pct']*blocks_intersect_temp['Hispanic']
    blocks_intersect_temp['Hispanic']=temp_array
    temp_array=blocks_intersect_temp['Area_Pct']*blocks_intersect_temp['Tot_Units']
    blocks_intersect_temp['Tot_Units']=temp_array
    temp_array=blocks_intersect_temp['Area_Pct']*blocks_intersect_temp['Owner']
    blocks_intersect_temp['Owner']=temp_array
    temp_array=blocks_intersect_temp['Area_Pct']*blocks_intersect_temp['TotalJobsR']
    blocks_intersect_temp['TotalJobsR']=temp_array
    temp_array=blocks_intersect_temp['Area_Pct']*blocks_intersect_temp['TotalUnder']
    blocks_intersect_temp['TotalUnder']=temp_array
    temp_array=blocks_intersect_temp['Area_Pct']*blocks_intersect_temp['TotalOver3']
    blocks_intersect_temp['TotalOver3']=temp_array
    temp_array=blocks_intersect_temp['Area_Pct']*blocks_intersect_temp['TotalJobsW']
    blocks_intersect_temp['TotalJobsW']=temp_array

    Final_estimates=blocks_intersect_temp.copy()
    Final_estimates['Pct_Urban']=(Final_estimates['Urb_Area']+Final_estimates['Urb_Clust'])/Final_estimates['TotPop']
    Final_estimates['Pct_White']=Final_estimates['White']/Final_estimates['TotPop']
    Final_estimates['Pct_Hisp']=Final_estimates['Hispanic']/Final_estimates['TotPop']
    Final_estimates['Pct_Black']=Final_estimates['Black']/Final_estimates['TotPop']
    Final_estimates['Pct_Owner']=Final_estimates['Owner']/Final_estimates['Tot_Units']
    Final_estimates['Pct_LowInc']=Final_estimates['TotalUnder']/Final_estimates['TotalJobsR']
    Final_estimates['Pct_HighInc']=Final_estimates['TotalOver3']/Final_estimates['TotalJobsR']
    column_list_final=['Pct_Urban', 'Pct_White', 'Pct_Hisp', 'Pct_Black', 'Pct_Owner', 'Pct_LowInc', 'Pct_HighInc', 'TotPop', 'Tot_Units', 'TotalJobsR', 'Shape_index','Area_Pct' ,'Area_Part', 'GISJOIN_City', 'GISJOIN_Block', 'Type', 'geometry']
    Final_estimates=Final_estimates[column_list_final].copy()
    stop_city_time = time.time()
    print(f'\nTime to complete processing: { (stop_city_time - start_city_time)/60:.2f} minutes\n')

    DF_list.append(Final_estimates)
    
end = time.time()

print(f'\nTime to complete processing: { (end - start)/60:.2f} minutes  \n')

# =============================================================================
# Saving results
# =============================================================================
DF_list_combined=pd.concat(DF_list)

State_ID_list=DF_list_combined['GISJOIN_City'].copy()

State_IDs = [z[1:3] for z in State_ID_list]
DF_list_combined['STATE_ID']=State_IDs

State_FIPS=[]
for state in State_IDs:
    State_FIPS.append(state_codes[state])
    
DF_list_combined['STATE']=State_FIPS

City_Parts_Data=DF_list_combined.copy()

City_Parts_Data=City_Parts_Data.drop_duplicates(keep='last')


City_Parts_Data['Sq_Mi']=City_Parts_Data.area*0.000000386102159
City_Parts_Data['Pop_Den']=City_Parts_Data['TotPop']/City_Parts_Data['Sq_Mi']

City_Parts_Data.to_pickle("City_Parts_00_to_10_Data.pkl")



City_Parts_Data.to_csv("City_Parts_00_to_10_Data.csv")
#%%

City_Parts_Data = pd.read_pickle("D:\OneDrive - Michigan State University\Spatial Annexation\City_Parts_00_to_10_Data.pkl")


#%%

# =============================================================================
# Statistics for census blocks
# =============================================================================

print(City_Parts_Data.groupby('Type').mean().T.round(2)) 



#%%

# =============================================================================
# # =============================================================================
# # # =============================================================================
# # # Processing cities for 2010 analysis
# # # =============================================================================
# # =============================================================================
# =============================================================================

Post_UTM_2=find_UTM(Post_DF)


print(Post_UTM_2.UTM.value_counts())

DF_list_5=[]

UTM_list=[10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

for UTM in UTM_list:
    exec("Post_UTM_2_{utm}=Post_UTM_2[Post_UTM_2['UTM']=={utm}]".format(utm=UTM))
    exec("DF_list_5.append(Post_UTM_2_{utm})".format(utm=UTM))

Error_list=[]
error=0
k=0
City_Stats=[]

cities_sindex, cities_newIDX=spatial_index(Post_DF) 
sindex_crs=cities_newIDX.crs
for DF in DF_list_5:
    t1_start = process_time()  
    
    print("DF number", k)
    
    epsg=DF['EPSG'].iloc[0]
    epsg_2=DF.crs
    DF=DF.to_crs(epsg)

    cities_newIDX_reproject=cities_newIDX.to_crs(epsg)
    print("finished data preparation")
    GDF_list=[]
# =============================================================================
#     Looping over cities, one at a time to identify shapes
# =============================================================================
    if len(DF)>0:
        for i in range(len(DF)):
            #print(DF['GISJOIN'].iloc[[i]])
            city_test=DF.iloc[[i]].copy()
            city_test_project=city_test.to_crs(sindex_crs)
            boundary=nearby_cities(city_test_project, cities_newIDX_reproject, cities_sindex, epsg)
            if (i+1)%50==0:
                print('Processing metrics for %sth city'%(i+1))
            try:
                city, city_buffer, city_buffer_reduce, other_cities, GISJOIN= find_shapes(city_test, boundary)
                #print(GISJOIN)
                Cities=cities_simple(city, other_cities, city_buffer, Pre_DF, GISJOIN, epsg )
                Cities=Cities.to_crs(epsg_2)
                Unincorporated = underbounded_simple(city, city_buffer_reduce, boundary, GISJOIN, epsg)
                Unincorporated=Unincorporated.to_crs(epsg_2)
                GDF_list.append(Unincorporated)
                GDF_list.append(Cities)

            except Exception as e:
                  Error_list.append([GISJOIN, e])
                  error+=1
                  print("Error number", error)
                  pass
    k+=1
    t1_stop = process_time()  
    print("processed ", (i+1), " cities in ", (t1_stop-t1_start)/60, "minutes")
    Cities_per_minute=(i+1)/((t1_stop-t1_start)/60)
    print("processed ", Cities_per_minute, " cities per minute ")
# =============================================================================
#     Combining and saving results
# =============================================================================

    GDF=gpd.GeoDataFrame(pd.concat(GDF_list, ignore_index=True), crs=GDF_list[0].crs)
    # Identifying cases where the geography is empty
    
    GDF['Empty_Geo']=GDF.is_empty
    
    GDF=GDF[GDF['Empty_Geo']!=True].copy()
    UTM=epsg[-2:]
    exec('GDF[GDF["Empty_Geo"]!=True].to_pickle("Type_Results_2010_UTM_{utm}.pkl")'.format(utm=UTM))

GDF_list_4=[]
for UTM in UTM_list:    
    exec('GDF=pd.read_pickle("D:\OneDrive - Michigan State University\Spatial Annexation\Type_Results_2010_UTM_{utm}.pkl")'.format(utm=UTM))
    GDF_list_4.append(GDF)

GDF=gpd.GeoDataFrame(pd.concat(GDF_list_4, ignore_index=True), crs=GDF_list_4[0].crs)

GDF=GDF.drop_duplicates(keep='last')

GDF[GDF["Empty_Geo"]!=True].to_pickle("Type_Results_2010_Final.pkl")


#%%

GDF=pd.read_pickle("Type_Results_2010_Final.pkl")

#%%

# =============================================================================
# # =============================================================================
# # Dissolving for city-level statistics for 2010 cross-sectional analysis
# # =============================================================================
# =============================================================================



national_crs=US_place_2010_name.crs
US_place_2010_name=US_place_2010_name[['NAME10', 'GISJOIN']].copy()

GDF_dissolve=GDF.copy()
GDF_dissolve=GDF_dissolve.to_crs(national_crs)
GDF_dissolve=GDF_dissolve.dissolve(by=['Type', 'GISJOIN'])
GDF_dissolve=GDF_dissolve.reset_index()


GDF_dissolve['AreaAcres']=GDF_dissolve.area*0.000247105

GDF_dissolve_1=GDF_dissolve.pivot_table(index="GISJOIN", columns='Type', values='AreaAcres')
GDF_dissolve_1=GDF_dissolve_1.fillna(0)
column_list=(GDF_dissolve_1.columns)
GDF_dissolve_1['Denominator']= GDF_dissolve_1['City']
for column in column_list:
    GDF_dissolve_1[column]=GDF_dissolve_1[column]/GDF_dissolve_1['Denominator']

GDF_dissolve_1=pd.merge(GDF_dissolve_1, US_place_2010_name, left_index=True, right_on='GISJOIN', how='inner')

GDF_dissolve=GDF_dissolve_1.copy()

GDF_dissolve.to_pickle("Type_Results_2010_Dissolve.pkl")
GDF_dissolve.to_csv("Type_Results_2010_Dissolve.csv")

#%%


GDF_dissolve = pd.read_pickle("D:\OneDrive - Michigan State University\Spatial Annexation\Type_Results_2010_Dissolve.pkl")
#%%

# =============================================================================
# Statistics for 2010 analysis
# =============================================================================

print("Dissolved statistics - city level")
print(GDF_dissolve.mean())

print("Number of areas by type")
print(GDF['Type'].value_counts())

GDF=GDF.to_crs(national_crs)
GDF['Area_Sq_Mi']=GDF.area*0.000000386102159

print("Size of area in square miles by type")
print(GDF.groupby('Type')['Area_Sq_Mi'].sum().round(0))


#%%
# =============================================================================
# # =============================================================================
# # # =============================================================================
# # # Creating map of underbounding by state
# # # =============================================================================
# # =============================================================================
# =============================================================================
epsg=US_place_2010.crs
GDF=GDF.to_crs(epsg)
GDF['Area']=GDF.area
GDF_wide=GDF.pivot_table(index="GISJOIN", columns='Type', values='Area')

GDF_wide['Underbounded']=GDF_wide['Underbounded'].fillna(0)
GDF_wide['City']=GDF_wide['City'].fillna(0)
GDF_wide['Other_cities']=GDF_wide['Other_cities'].fillna(0)
GDF_wide['Unaffected']=GDF_wide['Other_cities'].fillna(0)
GDF_wide['Pct_Under']=GDF_wide['Underbounded']/GDF_wide['City']
GDF_wide=GDF_wide.reset_index()
GDF_wide=GDF_wide[['GISJOIN', 'Pct_Under']].copy()

US_county_2010=gpd.read_file("US_county_2010.shp")
US_county_2010=US_county_2010[['GISJOIN', 'geometry']].copy()
US_county_2010=US_county_2010.rename(columns={'GISJOIN': 'GISJOIN_County', 'geometry': 'geometry_County'})
US_county_2010=US_county_2010.set_geometry('geometry_County')

US_state_2010=US_county_2010.copy()
US_state_2010['state']=US_state_2010['GISJOIN_County'].str[:3]
US_state_2010 = US_state_2010.dissolve(by='state')
US_state_2010=US_state_2010.reset_index()
US_state_2010=US_state_2010[US_state_2010['state']!="G02"].copy()
US_state_2010=US_state_2010[US_state_2010['state']!="G15"].copy()
US_state_2010=US_state_2010[US_state_2010['state']!="G72"].copy()

US_place_2010_GISJOIN=US_place_2010[['GISJOIN', 'geometry']].copy()

GDF_wide_centroid=pd.merge(US_place_2010_GISJOIN, GDF_wide, left_on='GISJOIN', right_on='GISJOIN', how='inner')

GDF_wide_centroid=GDF_wide_centroid.set_geometry('geometry')
GDF_wide_centroid=GDF_wide_centroid.copy()
GDF_wide_centroid['geometry']=GDF_wide_centroid['geometry'].centroid

State_intersection=gpd.sjoin(GDF_wide_centroid, US_state_2010, how='right', op='intersects')

State_Dissolve_1=State_intersection.groupby('state').Pct_Under.mean()
State_Dissolve_2=State_intersection.groupby('state').Pct_Under.count()

State_Dissolve_2=pd.DataFrame(State_Dissolve_2)

State_Dissolve_2=State_Dissolve_2.rename(columns={'Pct_Under': 'City_Count'})

US_state_2010_Pct_Under=pd.merge(US_state_2010, State_Dissolve_1, left_on="state", right_index=True, how='inner')
US_state_2010_Pct_Under=pd.merge(US_state_2010_Pct_Under, State_Dissolve_2, left_on="state", right_index=True, how='inner')
US_state_2010_Pct_Under=US_state_2010_Pct_Under[US_state_2010_Pct_Under['state']!='G11'].copy()

US_state_2010_Pct_Under.to_pickle("US_state_2010_Pct_Under.pkl")

#%%

US_state_2010_Pct_Under=pd.read_pickle("US_state_2010_Pct_Under.pkl")
#%%

fig = plt.figure()
ax=US_state_2010_Pct_Under.plot(column='Pct_Under', figsize=(12,12), legend=True, cmap="gist_heat_r", legend_kwds={'shrink': 0.3})
ax.axis('off')

plt.rcParams.update({'font.size': 13})
plt.show()


#%%

# =============================================================================
# # =============================================================================
# # # =============================================================================
# # # Joining city shapes from 2010 cross-sectional analysis to census blocks 
# # # =============================================================================
# # =============================================================================
# =============================================================================

Post_UTM_2=Post_UTM[['GISJOIN', 'UTM']].copy()


Selected_Areas=GDF[(GDF['Type']!='Blocked') & (GDF['Type']!='Other_cities') & (GDF['Type']!='City')].copy()
Selected_Areas['StateFIPS'] = Selected_Areas['GISJOIN'].str[1:3]


State_list=[]
state_list=list(Selected_Areas['StateFIPS'].unique())
for state in state_list:
    State=state_codes[state]
    State_list.append(State)

d_inverted = {v: k for k, v in state_codes.items()}
#state_list=list(state_codes.values())
#state_list_selected=[state_list[0]]+state_list[2:11]+state_list[12:]

start = time.time()
DF_list=[]
error_list=[]

# =============================================================================
# Looping over states
# =============================================================================

for State in State_list:
    print("Processing data for state: ", State)
    start_city_time = time.time()
    
    state_cities=Selected_Areas[Selected_Areas['StateFIPS']==d_inverted[State]].copy()
    
# =============================================================================
#     Getting census blocks
# =============================================================================
    
    exec('Blocks_Shapes_Data= gpd.read_file(r"D:/OneDrive - Michigan State University/Spatial Annexation/Static analysis/State blocks clean/{}_Blocks_Clean.shp")'.format(State))
    temp_crs=Blocks_Shapes_Data.crs
    state_cities=state_cities.to_crs(temp_crs)
    Blocks_Shapes_Data['Area_Pre']=Blocks_Shapes_Data.area
    Blocks_Shapes_Data['GISJOIN_Block']=Blocks_Shapes_Data['GISJOIN']    
    Blocks_Shapes_Data['geometry']=Blocks_Shapes_Data.buffer(0)
    state_cities['geometry']=state_cities.buffer(-.01) ## reducing by buffer of 1 centimeter to prevent non-noded intersections
    state_cities=state_cities.rename(columns={'GISJOIN': 'GISJOIN_City'})
    
# =============================================================================
#     Finding census blocks that intersect cities
# =============================================================================
    blocks_intersect=gpd.overlay(state_cities, Blocks_Shapes_Data, how='intersection')
    column_list=list(Blocks_Shapes_Data.columns)[2:-3]
    blocks_intersect['Area_Post']=blocks_intersect.area
    blocks_intersect['Area_Pct']=blocks_intersect['Area_Post']/blocks_intersect['Area_Pre']
    blocks_intersect=blocks_intersect[['GISJOIN_City','GISJOIN_Block', 'Type', 'Area_Pre', 'Area_Post', 'Area_Pct', 'geometry']+column_list].copy()
    blocks_intersect_temp=blocks_intersect[blocks_intersect['TotPop']>0].copy()
    blocks_intersect_temp=blocks_intersect_temp[blocks_intersect_temp['Area_Pct']>.99].copy()
    
    
# =============================================================================
#     Weighting population estimates by % of block that intersects with the city
# =============================================================================
    temp_array=blocks_intersect_temp['Area_Pct']*blocks_intersect_temp['TotPop']
    blocks_intersect_temp['TotPop']=temp_array
    temp_array=blocks_intersect_temp['Area_Pct']*blocks_intersect_temp['Urb_Area']
    blocks_intersect_temp['Urb_Area']=temp_array
    temp_array=blocks_intersect_temp['Area_Pct']*blocks_intersect_temp['Urb_Clust']
    blocks_intersect_temp['Urb_Clust']=temp_array
    temp_array=blocks_intersect_temp['Area_Pct']*blocks_intersect_temp['White']
    blocks_intersect_temp['White']=temp_array
    temp_array=blocks_intersect_temp['Area_Pct']*blocks_intersect_temp['Black']
    blocks_intersect_temp['Black']=temp_array
    temp_array=blocks_intersect_temp['Area_Pct']*blocks_intersect_temp['Hispanic']
    blocks_intersect_temp['Hispanic']=temp_array
    temp_array=blocks_intersect_temp['Area_Pct']*blocks_intersect_temp['Tot_Units']
    blocks_intersect_temp['Tot_Units']=temp_array
    temp_array=blocks_intersect_temp['Area_Pct']*blocks_intersect_temp['Owner']
    blocks_intersect_temp['Owner']=temp_array
    temp_array=blocks_intersect_temp['Area_Pct']*blocks_intersect_temp['TotalJobsR']
    blocks_intersect_temp['TotalJobsR']=temp_array
    temp_array=blocks_intersect_temp['Area_Pct']*blocks_intersect_temp['TotalUnder']
    blocks_intersect_temp['TotalUnder']=temp_array
    temp_array=blocks_intersect_temp['Area_Pct']*blocks_intersect_temp['TotalOver3']
    blocks_intersect_temp['TotalOver3']=temp_array
    temp_array=blocks_intersect_temp['Area_Pct']*blocks_intersect_temp['TotalJobsW']
    blocks_intersect_temp['TotalJobsW']=temp_array

    Final_estimates=blocks_intersect_temp.copy()
    Final_estimates['Pct_Urban']=(Final_estimates['Urb_Area']+Final_estimates['Urb_Clust'])/Final_estimates['TotPop']
    Final_estimates['Pct_White']=Final_estimates['White']/Final_estimates['TotPop']
    Final_estimates['Pct_Hisp']=Final_estimates['Hispanic']/Final_estimates['TotPop']
    Final_estimates['Pct_Black']=Final_estimates['Black']/Final_estimates['TotPop']
    Final_estimates['Pct_Owner']=Final_estimates['Owner']/Final_estimates['Tot_Units']
    Final_estimates['Pct_LowInc']=Final_estimates['TotalUnder']/Final_estimates['TotalJobsR']
    Final_estimates['Pct_HighInc']=Final_estimates['TotalOver3']/Final_estimates['TotalJobsR']
    column_list_final=['Pct_Urban', 'Pct_White', 'Pct_Hisp', 'Pct_Black', 'Pct_Owner', 'Pct_LowInc', 'Pct_HighInc', 'TotPop', 'Tot_Units', 'TotalJobsR', 'Area_Pct' ,'GISJOIN_City', 'GISJOIN_Block', 'Type', 'geometry']
    Final_estimates=Final_estimates[column_list_final].copy()
    stop_city_time = time.time()
    print(f'\nTime to complete processing: { (stop_city_time - start_city_time)/60:.2f} minutes\n')

    DF_list.append(Final_estimates)
    
end = time.time()

print(f'\nTime to complete processing: { (end - start)/60:.2f} minutes  \n')

# =============================================================================
# Saving results
# =============================================================================
DF_list_combined=pd.concat(DF_list)

State_ID_list=DF_list_combined['GISJOIN_City'].copy()

State_IDs = [z[1:3] for z in State_ID_list]
DF_list_combined['STATE_ID']=State_IDs

State_FIPS=[]
for state in State_IDs:
    State_FIPS.append(state_codes[state])
    
DF_list_combined['STATE']=State_FIPS

City_Parts_Data=DF_list_combined.copy()

City_Parts_Data=City_Parts_Data.drop_duplicates(keep='last')

City_Parts_Data['Sq_Mi']=City_Parts_Data.area*0.000000386102159
City_Parts_Data['Pop_Den']=City_Parts_Data['TotPop']/City_Parts_Data['Sq_Mi']


City_Parts_Data.to_pickle("City_Parts_10_Data.pkl")

City_Parts_Data_Export=City_Parts_Data[['Pct_Urban', 'Pct_White', 'Pct_Hisp', 'Pct_Black', 'Pct_Owner', 'Pct_LowInc', 'Pct_HighInc', 'TotPop', 'Tot_Units', 'TotalJobsR', 'Area_Pct' ,'Pop_Den', 'Sq_Mi' , 'GISJOIN_City', 'GISJOIN_Block', 'Type']].copy()

City_Parts_Data_Export.to_csv("City_Parts_10_Data.csv")

#%%
City_Parts_Data = pd.read_pickle("City_Parts_10_Data.pkl")



#%%

# =============================================================================
# Calculating nationwide statistics for 2010 cross-sectional analysis
# =============================================================================

Total_Fringe_Duplicated_Dropped=City_Parts_Data[['TotPop', 'Pct_White', 'Pct_Hisp', 'Pct_Black', 'Pct_Owner', 'Pct_LowInc',  'Pct_HighInc', 'Tot_Units', 'TotalJobsR' ,'GISJOIN_City', 'GISJOIN_Block', 'Type']].copy()
Total_Fringe_Duplicated_Dropped['Tot_White']=Total_Fringe_Duplicated_Dropped['TotPop']*Total_Fringe_Duplicated_Dropped['Pct_White']
Total_Fringe_Duplicated_Dropped['Tot_Hisp']=Total_Fringe_Duplicated_Dropped['TotPop']*Total_Fringe_Duplicated_Dropped['Pct_Hisp']
Total_Fringe_Duplicated_Dropped['Tot_Black']=Total_Fringe_Duplicated_Dropped['TotPop']*Total_Fringe_Duplicated_Dropped['Pct_Black']
Total_Fringe_Duplicated_Dropped['Tot_Owner']=Total_Fringe_Duplicated_Dropped['Tot_Units']*Total_Fringe_Duplicated_Dropped['Pct_Owner']
Total_Fringe_Duplicated_Dropped['Tot_HighInc']=Total_Fringe_Duplicated_Dropped['TotalJobsR']*Total_Fringe_Duplicated_Dropped['Pct_HighInc']

Total_Fringe_Duplicated_Dropped=Total_Fringe_Duplicated_Dropped.drop_duplicates(subset=['GISJOIN_Block'])

Total_Fringe=Total_Fringe_Duplicated_Dropped.groupby('Type').sum().round(0)

Total_Fringe=Total_Fringe[[list(Total_Fringe.columns)[0]]+list(Total_Fringe.columns)[-7:]].copy()

Total_Fringe['Pct_White']=Total_Fringe['Tot_White']/Total_Fringe['TotPop']
Total_Fringe['Pct_White']=Total_Fringe['Tot_White']/Total_Fringe['TotPop']
Total_Fringe['PctBlack']=Total_Fringe['Tot_Black']/Total_Fringe['TotPop']
Total_Fringe['Pct_Hisp']=Total_Fringe['Tot_Hisp']/Total_Fringe['TotPop']
Total_Fringe['Pct_Owner']=Total_Fringe['Tot_Owner']/Total_Fringe['Tot_Units']
Total_Fringe['Pct_HighInc']=Total_Fringe['Tot_HighInc']/Total_Fringe['TotalJobsR']
#%%
print("Nationwide statistics for 2010 analysis")

print(Total_Fringe.T.round(2))

#%%
Total_Fringe_Duplicated_Dropped['dissolve']=Total_Fringe_Duplicated_Dropped['GISJOIN_City']+Total_Fringe_Duplicated_Dropped['Type']
Total_Fringe_Join=Total_Fringe_Duplicated_Dropped[['GISJOIN_City', 'Type', 'dissolve']].copy()
Total_Fringe_Join= Total_Fringe_Join.drop_duplicates(keep='last')
Total_Fringe_Duplicated_Dropped['dissolve']=Total_Fringe_Duplicated_Dropped['GISJOIN_City']+Total_Fringe_Duplicated_Dropped['Type']
Total_Fringe_Cities=Total_Fringe_Duplicated_Dropped.groupby('dissolve').sum().round(0)
Total_Fringe_Cities=Total_Fringe_Cities[[list(Total_Fringe_Cities.columns)[0]]+list(Total_Fringe_Cities.columns)[-7:]].copy()

Total_Fringe_Cities=pd.merge(Total_Fringe_Cities, Total_Fringe_Join, left_on='dissolve', right_on='dissolve', how='left')

Total_Fringe_Nation=Total_Fringe_Cities.groupby('Type')['TotPop', 'Tot_White'].sum()

Total_Fringe_Cities['Pct_White']=Total_Fringe_Cities['Tot_White']/Total_Fringe_Cities['TotPop']
Total_Fringe_Cities['Pct_Black']=Total_Fringe_Cities['Tot_Black']/Total_Fringe_Cities['TotPop']
Total_Fringe_Cities['Pct_Hisp']=Total_Fringe_Cities['Tot_Hisp']/Total_Fringe_Cities['TotPop']
Total_Fringe_Cities['Pct_Owner']=Total_Fringe_Cities['Tot_Owner']/Total_Fringe_Cities['Tot_Units']
Total_Fringe_Cities['Pct_HighInc']=Total_Fringe_Cities['Tot_HighInc']/Total_Fringe_Cities['TotalJobsR']

Total_Fringe_Cities=Total_Fringe_Cities.pivot_table(index="GISJOIN_City", columns='Type', values=list(Total_Fringe_Cities.columns))



Total_Fringe_Cities=Total_Fringe_Cities.fillna(0)
Total_Fringe_Cities.columns = [f'{i}{j}' for i, j in Total_Fringe_Cities.columns]

#Importing place shapefile
US_place_2010=gpd.read_file('D:/OneDrive - Michigan State University/Spatial Annexation/Automated Identification/shapefiles/US_place_2010.shp')

# Identifying cities, dropping census designated places
US_place_2010=US_place_2010.loc[(US_place_2010['CLASSFP10'] =='C1') | (US_place_2010['CLASSFP10'] == 'C2') | (US_place_2010['CLASSFP10'] == 'C3')| (US_place_2010['CLASSFP10'] == 'C4')| (US_place_2010['CLASSFP10'] == 'C5')| (US_place_2010['CLASSFP10'] == 'C6')| (US_place_2010['CLASSFP10'] == 'C7')| (US_place_2010['CLASSFP10'] == 'C8')| (US_place_2010['CLASSFP10'] == 'C9')]

US_place_2010_name= US_place_2010[['geometry', 'GISJOIN', 'NAME10', 'STATEFP10']].copy()

US_place_2010_name=US_place_2010_name.replace({"STATEFP10": state_codes})


Total_Fringe_Cities=pd.merge(Total_Fringe_Cities,US_place_2010_name, left_on='GISJOIN_City', right_on='GISJOIN', how='inner')

Total_Fringe_Cities=Total_Fringe_Cities.rename(columns={'NAME10': 'Name', 'STATEFP10': 'State'})

Total_Fringe_Cities['Diff_Pct_White']=Total_Fringe_Cities['Pct_WhiteUnderbounded']-Total_Fringe_Cities['Pct_WhiteUnaffected']
Total_Fringe_Cities['Diff_Pct_Black']=Total_Fringe_Cities['Pct_BlackUnderbounded']-Total_Fringe_Cities['Pct_BlackUnaffected']
Total_Fringe_Cities['Diff_Pct_Hisp']=Total_Fringe_Cities['Pct_HispUnderbounded']-Total_Fringe_Cities['Pct_HispUnaffected']

Total_Fringe_Cities=Total_Fringe_Cities[[ 'Name', 'State', 'TotPopUnderbounded' ,'Tot_WhiteUnderbounded', 'TotPopUnaffected' ,'Tot_WhiteUnaffected', 'Diff_Pct_White', 'Diff_Pct_Black', 'Diff_Pct_Hisp']].copy()
#Total_Fringe_Cities=Total_Fringe_Cities[[ 'Name', 'TotPopUnderbounded' ,'TotPopUnaffected', 'Diff_Pct_White', 'Diff_Pct_Black', 'Diff_Pct_Hisp']].copy()
#%%

print("Statistics for 50 largest metros - 2010 analysis")

print(Total_Fringe_Cities.nlargest(50, ['TotPopUnderbounded']).round(2))
#%%
print("Nationwide statistics for 2010 analysis")

print(Total_Fringe_Cities[['TotPopUnderbounded' ,'Tot_WhiteUnderbounded', 'TotPopUnaffected' ,'Tot_WhiteUnaffected']].nlargest(50, ['TotPopUnderbounded']).round(2))


#%%

# =============================================================================
# Temporal analysis - 1980 to 2019
# =============================================================================

# Identifying crosswalk ID from 1980 to 2019

Place_Pop_80_to_10=pd.read_csv('D:/OneDrive - Michigan State University/Spatial Annexation/Automated Identification/shapefiles/Place_Pop_80_to_10.csv')
Place_Pop_80_to_10['CDP_80']=Place_Pop_80_to_10['NAME1980'].str.contains('CDP')
Place_Pop_80_to_10=Place_Pop_80_to_10.rename(columns={'GJOIN2012':'GISJOIN_10'})
Place_Pop_80_to_10=Place_Pop_80_to_10.rename(columns={'GJOIN1980':'GISJOIN_80'})
Place_Pop_80_to_10=Place_Pop_80_to_10[['NHGISCODE', 'GISJOIN_10' ,'GISJOIN_80', 'CDP_80']].copy()


#Importing 1980 and 2019 place shapefiles
US_place_1980=gpd.read_file('D:/OneDrive - Michigan State University/Spatial Annexation/Automated Identification/shapefiles/US_place_1980.shp')
US_place_1980['Area_80']=(US_place_1980.area*.0000003861021585).round(5)
working_crs=US_place_1980.crs
US_place_1980=US_place_1980.rename(columns={'geometry':'geometry_80'})
US_place_1980=US_place_1980.rename(columns={'GISJOIN':'GISJOIN_80'})
US_place_1980=US_place_1980[['GISJOIN_80', 'geometry_80', 'Area_80']].copy()


US_place_2010=gpd.read_file('D:/OneDrive - Michigan State University/Spatial Annexation/Automated Identification/shapefiles/US_place_2010.shp')
US_place_2010['Area_10']=(US_place_2010.area*.0000003861021585).round(5)
US_place_2010=US_place_2010.rename(columns={'geometry':'geometry_10'})
US_place_2010=US_place_2010.rename(columns={'GISJOIN':'GISJOIN_10'})
US_place_2010['CDP_10']=US_place_2010['NAME10'].str.contains('CDP')
US_place_2010=US_place_2010[['GISJOIN_10', 'CDP_10', 'geometry_10', 'Area_10']].copy()

# Dropping CDPs or any place not found in both years
Place_80_to_10=pd.merge(Place_Pop_80_to_10, US_place_1980, left_on='GISJOIN_80', right_on='GISJOIN_80', how='outer')
Place_80_to_10=pd.merge(Place_80_to_10, US_place_2010, left_on='GISJOIN_10', right_on='GISJOIN_10', how='outer')
Place_80_to_10=Place_80_to_10[(Place_80_to_10['CDP_10']!=True)].copy()
Place_80_to_10=Place_80_to_10[(Place_80_to_10['CDP_80']!=True)].copy()
Place_80_to_10=Place_80_to_10.dropna()


Place_80_to_10=Place_80_to_10.set_geometry('geometry_80')
Place_80_to_10=Place_80_to_10.set_crs(working_crs)


Place_Pop_80_to_10=pd.read_csv('D:/OneDrive - Michigan State University/Spatial Annexation/Automated Identification/shapefiles/Place_Pop_80_to_10.csv')

US_city_1980=gpd.read_file('D:/OneDrive - Michigan State University/Spatial Annexation/Automated Identification/shapefiles/US_place_1980.shp')
US_city_1980=pd.merge(Place_Pop_80_to_10, US_city_1980, left_on='GJOIN1980', right_on='GISJOIN', how='inner')

US_city_1980['CDP_80']=US_city_1980['NAME1980'].str.contains('CDP')
US_city_1980=US_city_1980[US_city_1980['CDP_80']!='CDP'].copy()
US_city_1980=US_city_1980[['GISJOIN', 'geometry']].copy()
US_city_1980=US_city_1980.set_geometry('geometry')


# Measuring change in size of city
Place_80_to_10['Area_Pct_Chg']=100*((Place_80_to_10['Area_10']/Place_80_to_10['Area_80'])-1)


Place_80_to_10_Growing=Place_80_to_10[Place_80_to_10['Area_Pct_Chg']>10].copy()
Place_80_to_10_Growing['state']=Place_80_to_10_Growing['NHGISCODE'].str[:3]

Place_80_to_10_Growing=Place_80_to_10_Growing[Place_80_to_10_Growing['state']!="G02"].copy()
Place_80_to_10_Growing=Place_80_to_10_Growing[Place_80_to_10_Growing['state']!="G15"].copy()
Place_80_to_10_Growing=Place_80_to_10_Growing[Place_80_to_10_Growing['state']!="G72"].copy()



Place_80_to_10_Growing=Place_80_to_10_Growing[['GISJOIN_80', 'geometry_80']].copy()

Place_80_to_10_Growing=Place_80_to_10_Growing.rename(columns={'GISJOIN_80': 'GISJOIN', 'geometry_80': 'geometry'})


Place_80_to_10_Growing=Place_80_to_10_Growing.set_geometry('geometry')
#%%

# =============================================================================
# # =============================================================================
# # # =============================================================================
# # # Processing cities for 1980 to 2010 analysis
# # # =============================================================================
# # =============================================================================
# =============================================================================

Post_UTM_2=find_UTM(Place_80_to_19_Growing)


print(Post_UTM_2.UTM.value_counts())

DF_list_5=[]

UTM_list=list(Post_UTM_2['UTM'].unique())

for UTM in UTM_list:
    exec("Post_UTM_2_{utm}=Post_UTM_2[Post_UTM_2['UTM']=={utm}]".format(utm=UTM))
    exec("DF_list_5.append(Post_UTM_2_{utm})".format(utm=UTM))

Error_list=[]
error=0
k=0
City_Stats=[]

cities_sindex, cities_newIDX=spatial_index(US_city_1980) 
sindex_crs=cities_newIDX.crs
for DF in DF_list_5:
    t1_start = process_time()  
    
    print("DF number", k)
    
    epsg=DF['EPSG'].iloc[0]
    epsg_2=DF.crs
    DF=DF.to_crs(epsg)

    cities_newIDX_reproject=cities_newIDX.to_crs(epsg)
    print("finished data preparation")
    GDF_list=[]
# =============================================================================
#     Looping over cities, one at a time to identify shapes
# =============================================================================
    if len(DF)>0:
        for i in range(len(DF)):
            #print(DF['GISJOIN'].iloc[[i]])
            city_test=DF.iloc[[i]].copy()
            city_test_project=city_test.to_crs(sindex_crs)
            boundary=nearby_cities(city_test_project, cities_newIDX_reproject, cities_sindex, epsg)
            if (i+1)%50==0:
                print('Processing metrics for %sth city'%(i+1))
            try:
                city, city_buffer, city_buffer_reduce, other_cities, GISJOIN= find_shapes(city_test, boundary)
                #print(GISJOIN)
                Cities=cities_simple(city, other_cities, city_buffer, Pre_DF, GISJOIN, epsg )
                Cities=Cities.to_crs(epsg_2)
                Unincorporated = underbounded_simple(city, city_buffer_reduce, boundary, GISJOIN, epsg)
                Unincorporated=Unincorporated.to_crs(epsg_2)
                GDF_list.append(Unincorporated)
                GDF_list.append(Cities)

            except Exception as e:
                  Error_list.append([GISJOIN, e])
                  error+=1
                  print("Error number", error)
                  pass
    k+=1
    t1_stop = process_time()  
    print("processed ", (i+1), " cities in ", (t1_stop-t1_start)/60, "minutes")
    Cities_per_minute=(i+1)/((t1_stop-t1_start)/60)
    print("processed ", Cities_per_minute, " cities per minute ")
# =============================================================================
#     Combining and saving results
# =============================================================================

    GDF=gpd.GeoDataFrame(pd.concat(GDF_list, ignore_index=True), crs=GDF_list[0].crs)
    # Identifying cases where the geography is empty
    
    GDF['Empty_Geo']=GDF.is_empty
    
    GDF=GDF[GDF['Empty_Geo']!=True].copy()
    UTM=epsg[-2:]
    exec('GDF[GDF["Empty_Geo"]!=True].to_pickle("Type_Results_1980_UTM_{utm}.pkl")'.format(utm=UTM))

GDF_list_4=[]
for UTM in UTM_list:    
    exec('GDF=pd.read_pickle("D:\OneDrive - Michigan State University\Spatial Annexation\Type_Results_1980_UTM_{utm}.pkl")'.format(utm=UTM))
    GDF_list_4.append(GDF)

GDF=gpd.GeoDataFrame(pd.concat(GDF_list_4, ignore_index=True), crs=GDF_list_4[0].crs)

GDF=GDF.drop_duplicates(keep='last')

GDF[GDF["Empty_Geo"]!=True].to_pickle("Type_Results_1980_Final.pkl")


#%%


Type_Results_80_Final = pd.read_pickle("Type_Results_1980_Final.pkl")

Type_Results_10_Final = pd.read_pickle("Type_Results_2010_Final.pkl")
#%%
US_place_2010=gpd.read_file('D:/OneDrive - Michigan State University/Spatial Annexation/Automated Identification/shapefiles/US_place_2010.shp')

# Identifying cities, dropping census designated places
US_place_2010=US_place_2010.loc[(US_place_2010['CLASSFP10'] =='C1') | (US_place_2010['CLASSFP10'] == 'C2') | (US_place_2010['CLASSFP10'] == 'C3')| (US_place_2010['CLASSFP10'] == 'C4')| (US_place_2010['CLASSFP10'] == 'C5')| (US_place_2010['CLASSFP10'] == 'C6')| (US_place_2010['CLASSFP10'] == 'C7')| (US_place_2010['CLASSFP10'] == 'C8')| (US_place_2010['CLASSFP10'] == 'C9')]

US_place_2010_name= US_place_2010[['geometry', 'GISJOIN', 'NAME10', 'STATEFP10']].copy()

# saving only geometry and GISJOIN codes
US_place_2010 = US_place_2010[['geometry', 'GISJOIN']].copy()


working_crs=US_place_2010.crs

#%%
def Clean_results(Type_Results_80_Final, Type_Results_10_Final):
    Type_Results_10_Final_Dissovle=Type_Results_10_Final.copy()
    Type_Results_10_Final_Dissovle=Type_Results_10_Final_Dissovle.to_crs(working_crs)
    Type_Results_10_Final_Dissovle['Area_SqMi']=Type_Results_10_Final_Dissovle.area*0.000000386102159
    Type_Results_10_Final_Dissovle=Type_Results_10_Final_Dissovle.groupby(['GISJOIN', 'Type'])['Area_SqMi'].sum()
    Type_Results_10_Final_Dissovle=Type_Results_10_Final_Dissovle.reset_index()
    Type_Results_10_Final_Dissovle=Type_Results_10_Final_Dissovle.pivot_table(index="GISJOIN", columns='Type', values='Area_SqMi')
    Type_Results_10_Final_Dissovle=pd.merge(Type_Results_10_Final_Dissovle, Place_Pop_80_to_10, left_on='GISJOIN', right_on='GJOIN2010', how='inner')
    Type_Results_10_Final_Dissovle=Type_Results_10_Final_Dissovle[list(Type_Results_10_Final_Dissovle.columns)[:5]+['NHGISCODE','GJOIN1980', 'GJOIN2010']].copy()
    Type_Results_10_Final_Dissovle=pd.merge(Type_Results_10_Final_Dissovle, Place_80_to_10_Growing, left_on='GJOIN1980', right_on='GISJOIN', how='inner')
    Type_Results_10_Final_Dissovle=Type_Results_10_Final_Dissovle.fillna(0)
    column_list=list(Type_Results_10_Final_Dissovle.columns)[:5]
    Type_Results_10_Final_Dissovle['Total_Area']= Type_Results_10_Final_Dissovle.sum(axis=1)
    
    
    
    Type_Results_80_Final_Dissovle=Type_Results_80_Final.copy()
    Type_Results_80_Final_Dissovle=Type_Results_80_Final_Dissovle.to_crs(working_crs)
    Type_Results_80_Final_Dissovle['Area_SqMi']=Type_Results_80_Final_Dissovle.area*0.000000386102159
    Type_Results_80_Final_Dissovle=Type_Results_80_Final_Dissovle.groupby(['GISJOIN', 'Type'])['Area_SqMi'].sum()
    Type_Results_80_Final_Dissovle=Type_Results_80_Final_Dissovle.reset_index()
    Type_Results_80_Final_Dissovle=Type_Results_80_Final_Dissovle.pivot_table(index="GISJOIN", columns='Type', values='Area_SqMi')
    
    Type_Results_80_Final_Dissovle=pd.merge(Type_Results_80_Final_Dissovle, Place_Pop_80_to_10, left_on='GISJOIN', right_on='GJOIN1980', how='inner')
    Type_Results_80_Final_Dissovle=Type_Results_80_Final_Dissovle[list(Type_Results_80_Final_Dissovle.columns)[:5]+['NHGISCODE', 'GJOIN1980', 'GJOIN2010']].copy()
    Type_Results_80_Final_Dissovle=pd.merge(Type_Results_80_Final_Dissovle, Place_80_to_10_Growing, left_on='GJOIN1980', right_on='GISJOIN', how='inner')
    Type_Results_80_Final_Dissovle=Type_Results_80_Final_Dissovle.fillna(0)
    column_list=list(Type_Results_80_Final_Dissovle.columns)[:5]
    Type_Results_80_Final_Dissovle['Total_Area']= Type_Results_80_Final_Dissovle.sum(axis=1)
    return Type_Results_80_Final_Dissovle, Type_Results_10_Final_Dissovle


Type_Results_80_Final_Dissovle, Type_Results_10_Final_Dissovle=Clean_results(Type_Results_80_Final, Type_Results_10_Final)
#%%
print("1980 Size of Underbounded Ares in SqMi",(Type_Results_80_Final_Dissovle['Underbounded']).mean())
print("2010 Size of Underbounded Ares in SqMi",(Type_Results_10_Final_Dissovle['Underbounded']).mean())
print("1980 Ratio of Underbounded to City Area",(Type_Results_80_Final_Dissovle['Underbounded']/Type_Results_80_Final_Dissovle['City']).mean())
print("2010 Ratio of Underbounded to City Area",(Type_Results_10_Final_Dissovle['Underbounded']/Type_Results_10_Final_Dissovle['City']).mean())

#%%

# =============================================================================
# # =============================================================================
# # # =============================================================================
# # # Historical analysis of five cities
# # # =============================================================================
# # =============================================================================
# =============================================================================

# Setting working directory
os.chdir('D:\OneDrive - Michigan State University\Spatial Annexation')
os.getcwd()

new_crs='ESRI:102003'
#%%


#Importing Columbus data
Columbus=gpd.read_file('D:\OneDrive - Michigan State University\Spatial Annexation\City Historical Annexations\Columbus, OH\Corporate_Boundary_Annexations.shp')
Columbus['ACCEPTANCE'] = pd.to_datetime(Columbus['ACCEPTANCE'])
Columbus['year']= Columbus['ACCEPTANCE'].dt.year
Columbus['City']='Columbus, OH'
Columbus=Columbus[['year', 'City', 'geometry']].copy()
Columbus=Columbus.to_crs(new_crs)

#Importing Albuqurque data
Albuquerque=gpd.read_file('D:\OneDrive - Michigan State University\Spatial Annexation\City Historical Annexations\Albequerque, NM\\annexations.shp')
Albuquerque['ORDINANCED'] = pd.to_datetime(Albuquerque['ORDINANCED'])
Albuquerque['year']= Albuquerque['ORDINANCED'].dt.year
Albuquerque['City']='Albuquerque, NM'
Albuquerque=Albuquerque[['year', 'City', 'geometry']].copy()
Albuquerque=Albuquerque.to_crs(new_crs)

#Importing Huntsville data
Huntsville=gpd.read_file('D:\OneDrive - Michigan State University\Spatial Annexation\City Historical Annexations\Huntsville, AL\\Boundary_City_Annexation_poly.shp')
Huntsville['Eff_Date'] = pd.to_datetime(Huntsville['Eff_Date'])
Huntsville['year']= Huntsville['Eff_Date'].dt.year
Huntsville['City']='Huntsville, AL'
Huntsville=Huntsville[['year', 'City', 'geometry']].copy()
Huntsville=Huntsville.to_crs(new_crs)

#Importing San Jose data
San_Jose=gpd.read_file('D:\OneDrive - Michigan State University\Spatial Annexation\City Historical Annexations\San Jose, CA\\city-of-san-jose-annexations.shp')
San_Jose=San_Jose[San_Jose['CURRENT_ST']=='Annexed'].copy()
San_Jose['ANNEX_DATE'] = pd.to_datetime(San_Jose['ANNEX_DATE'])
San_Jose['year']= San_Jose['ANNEX_DATE'].dt.year
San_Jose['City']='San Jose, CA'
San_Jose=San_Jose[['year', 'City', 'geometry']].copy()
San_Jose=San_Jose.to_crs(new_crs)

#Importing Aurora data
Aurora=gpd.read_file('D:\OneDrive - Michigan State University\Spatial Annexation\City Historical Annexations\Aurora, CO\\Annexations.shp')
Aurora['ANNEXATI_2'] = pd.to_datetime(Aurora['ANNEXATI_2'])
Aurora['year']= Aurora['ANNEXATI_2'].dt.year
Aurora['City']='Aurora, CO'
Aurora=Aurora[['year', 'City', 'geometry']].copy()
Aurora=Aurora.to_crs(new_crs)

#%%

# =============================================================================
# Processing data for five cities
# =============================================================================

City_DF_List=[Albuquerque, Aurora,  Columbus, Huntsville, San_Jose]
Five_City_List=['G35002000', 'G08004000', 'G39018000', 'G01037000', 'G06068000']

Counter_List=[]
Year_List=[]
Avoided_List=[]
Avoided_Ratio_List=[]
City_Area_List=[]
Avoided_Area_List=[]
City_Name_List=[]
Geometry_List=[]
GISJOIN_List=[]
new_part=None
Error_list=[]
for city_count, City_Selected in enumerate(City_DF_List):
    try:
        city_selected=City_Selected.iloc[0]['City']
        print(city_selected)
        City_Selected['decade']=City_Selected['year']*.1
        City_Selected['decade']=City_Selected['decade'].apply(np.ceil)
        City_Selected['decade']=City_Selected['decade']*10
        City_Selected_Years = City_Selected.dissolve(by='decade')
        City_Selected_Years=City_Selected_Years.sort_values(by='year')
        City_Selected_Years=City_Selected_Years.reset_index()
        City_Selected_Years['year']=City_Selected_Years['decade']
        counter=0
    
        for index, row in City_Selected_Years.iterrows():  
            
            year=row['year']
            print(year)
            if counter==0:
                new_part=row['geometry']
                old_city=new_part
            else:
                new_part=row['geometry']
                old_city = old_city.buffer(1).union(new_part)
            box_min = old_city.minimum_rotated_rectangle
            x, y = box_min.exterior.coords.xy
            edge_length = (Point(x[0], y[0]).distance(Point(x[1], y[1])), Point(x[1], y[1]).distance(Point(x[2], y[2])))
            width=min(edge_length)/2
            city_buffer=old_city.buffer(width)
            city_buffer_reduce=city_buffer.buffer(-1*width)
            avoided=city_buffer_reduce.buffer(1).difference(old_city)
            Year_List.append(year)
            Avoided_Ratio=avoided.area/(old_city.area)
            Avoided_Ratio_List.append(Avoided_Ratio)
            City_Area_List.append(old_city.area)
            Avoided_Area_List.append(avoided.area)
            City_Name_List.append(city_selected)
            Counter_List.append(counter)
            Geometry_List.append(new_part)
            GISJOIN=Five_City_List[city_count]
            GISJOIN_List.append(GISJOIN)
            counter+=1
            print(counter)
    except:
        Error_list.append(city_selected)
        pass


    Historical_results = pd.DataFrame(
        {'year': Year_List,
         'City_Area': City_Area_List,
         'Under_Area': Avoided_Area_List,
         'Under_Ratio': Avoided_Ratio_List,
         'City': City_Name_List,
         'Counter': Counter_List,
         'geometry': Geometry_List,
         'GISJOIN': GISJOIN_List
           })

Historical_results=Historical_results[Historical_results['year']<2021].copy()
AvoideHistorical_resultsd_Areas = Historical_results.sort_values(['City', 'year'])
Historical_results['City_Area']=Historical_results['City_Area']*.000000386
Historical_results['City_Area']=Historical_results['City_Area'].round(3)
Historical_results['Under_Area']=Historical_results['Under_Area']*.000000386
Historical_results['Under_Area']=Historical_results['Under_Area'].round(3)


Historical_results=Historical_results.set_geometry('geometry')
#%%
Historical_results['Decade_recent']=Historical_results['year']
condition=Historical_results.year>1980
Historical_results.loc[condition, 'Decade_recent']=2020
condition=Historical_results.year<1990
Historical_results.loc[condition, 'Decade_recent']=1980
condition=Historical_results.year<1950
Historical_results.loc[condition, 'Decade_recent']=1940

new_crs='ESRI:102003'
Historical_results=Historical_results.set_crs(new_crs)
#%%
Historical_results.to_pickle("Historical_results_spatial.pkl")

#%%
Historical_results_spatial=pd.read_pickle("Historical_results_spatial.pkl")
#%%

# =============================================================================
# # =============================================================================
# # # =============================================================================
# # # Creating maps of five cities
# # # =============================================================================
# # =============================================================================
# =============================================================================
new_crs='ESRI:102003'
Five_City_List=['G35002000', 'G08004000', 'G39018000', 'G01037000', 'G06068000']

def draw_maps(Historical_results_spatial, Five_City_List, GDF, n_rows, n_cols):
    fig=plt.figure(figsize=(12,12))

    for i, GISJOIN in enumerate(Five_City_List):
        city=GDF[(GDF['GISJOIN']==GISJOIN) & (GDF['Type']=='City')].copy()
        city=city.dissolve(by='GISJOIN')
        city=city.to_crs(new_crs)
        city=city.iloc[0].geometry.buffer(0)
    
        box_min = city.minimum_rotated_rectangle
        x, y = box_min.exterior.coords.xy
        edge_length = (Point(x[0], y[0]).distance(Point(x[1], y[1])), Point(x[1], y[1]).distance(Point(x[2], y[2])))
        width=min(edge_length)/2
        city_buffer_reduce=city.buffer(width).buffer(-.98*width)
        city_buffer_reduce=pd.DataFrame([city_buffer_reduce])
        city_buffer_reduce=city_buffer_reduce.rename(columns={0:'geomtery'})
        city_buffer_reduce=city_buffer_reduce.set_geometry('geomtery')
        city_buffer_reduce=city_buffer_reduce.set_crs(new_crs)
        ax=fig.add_subplot(n_rows,n_cols,i+1)
        
        concave_plot=city_buffer_reduce.plot(ax=ax,color='none',ec='darkgrey', linewidth=2, figsize=(10,10))
        
        city_plot=Historical_results_spatial[(Historical_results_spatial['GISJOIN']==GISJOIN)].plot(column='Decade_recent', cmap='autumn_r', aspect=1, ax=concave_plot)
        
        other_cities=GDF[(GDF['GISJOIN']==GISJOIN) & (GDF['Type']=='Other_cities') |
                         (GDF['GISJOIN']==GISJOIN) & (GDF['Type']=='Blocked')].copy()
        other_cities=other_cities.to_crs(new_crs)
        other_cities=gpd.overlay(other_cities, city_buffer_reduce, how='intersection')
        other_cities.plot(color='black', ax=city_plot, aspect=1, figsize=(10,10))
        ax.axis('off')
        City_name=Historical_results_spatial[(Historical_results_spatial['GISJOIN']==GISJOIN)]['City'].iloc[0]
        ax.set_title(City_name)
        
    #fig.add_subplot(n_rows,n_cols,i+1)
    
    import matplotlib.patches as mpatches
    early = mpatches.Patch(color='yellow',edgecolor='yellow', label='Before 1940')
    mid = mpatches.Patch(color='orange',edgecolor='orange', label='1940 to 1980')
    late = mpatches.Patch(color='red', edgecolor='red',label='1980 to 2020')
    other_cities = mpatches.Patch(color='black', label='Other Cities (2010)')
    concave_hull = mpatches.Patch(color='none', ec='darkgrey',linewidth=2, label='Underbounded Areas')
    ax.axis('off')
    fig.add_subplot(n_rows,n_cols,6)
    plt.axis('off')
    plt.legend(handles=[early, mid, late, other_cities, concave_hull], 
               loc='center',ncol=1, prop={'size': 15}, 
               fancybox=True,  borderpad=1)
    fig.tight_layout()  
    plt.savefig('Five_city_maps.png', dpi=600)
    plt.show()
    

draw_maps(Historical_results_spatial, Five_City_List,GDF,2, 3)

##%
#%%


Five_cities=City_Parts_Data[(City_Parts_Data['GISJOIN_City']=='G35002000')|(City_Parts_Data['GISJOIN_City']=='G08004000')|
                         (City_Parts_Data['GISJOIN_City']=='G39018000')|(City_Parts_Data['GISJOIN_City']=='G01037000')|
                         (City_Parts_Data['GISJOIN_City']=='G06068000')].copy()

print(Five_cities.groupby('STATE').mean().T.round(2))
#%%
City_Parts_Data[(City_Parts_Data['GISJOIN_City']=='G35002000')].plot(column='Pct_Hisp', legend=True, figsize=(10,10))

#%%
# =============================================================================
# Adjusting estimates to correct for other neighboring cities
# =============================================================================

GDF=pd.read_pickle("Type_Results_2010_Final.pkl")

national_crs=US_place_2010_name.crs
US_place_2010_name=US_place_2010_name[['NAME10', 'GISJOIN']].copy()

GDF_dissolve=GDF.copy()
GDF_dissolve=GDF_dissolve.to_crs(national_crs)
GDF_dissolve=GDF_dissolve.dissolve(by=['Type', 'GISJOIN'])
GDF_dissolve=GDF_dissolve.reset_index()
GDF_dissolve['Sq_Mi']=GDF_dissolve.area*0.000000386102159
GDF_dissolve_1=GDF_dissolve.pivot_table(index="GISJOIN", columns='Type', values='Sq_Mi')
GDF_dissolve_1=GDF_dissolve_1.fillna(0)
column_list=(GDF_dissolve_1.columns)
GDF_dissolve_1=pd.merge(GDF_dissolve_1, US_place_2010_name, left_index=True, right_on='GISJOIN', how='inner')
GDF_dissolve=GDF_dissolve_1.copy()


Under_SqMi_List=[]
Five_City_List=['G35002000', 'G08004000', 'G39018000', 'G01037000', 'G06068000']
for city in Five_City_List:
    under=GDF_dissolve[GDF_dissolve['GISJOIN']==city]['Underbounded'].iloc[0]
    Under_SqMi_List.append(under)

Under_Area_Correct=pd.DataFrame(Under_SqMi_List)
Under_Area_Correct=Under_Area_Correct.rename(columns={0:'Under_Correct'})
#
City_Name=pd.DataFrame(list(Historical_results['City'].unique()))
City_Name=City_Name.rename(columns={0:'City'})
Under_Wrong=Historical_results[Historical_results['year']==2010]['Under_Area'].copy()
Under_Wrong=Under_Wrong.reset_index(drop=True)
Under_Wrong=pd.DataFrame(Under_Wrong)
Under_Wrong=Under_Wrong.rename(columns={'Under_Area':'Under_Wrong'})
Merge_1=pd.merge(Under_Area_Correct,City_Name, left_index=True, right_index=True, how='inner')
Merge_2=pd.merge(Merge_1,Under_Wrong, left_index=True, right_index=True, how='inner')
Merge_2['Adjustment']=Merge_2['Under_Correct']/Merge_2['Under_Wrong']
Historical_results=pd.merge(Historical_results, Merge_2, left_on='City', right_on='City', how='inner')
Historical_results['Under_Area']=Historical_results['Under_Area']*Historical_results['Adjustment']
Historical_results['Under_Ratio']=Historical_results['Under_Ratio']*Historical_results['Adjustment']

Historical_results.to_pickle("Historical_results.pkl")
#%%

Historical_results=pd.read_pickle("Historical_results.pkl")
#%%

# =============================================================================
# Drawing plots
# =============================================================================

plt.rcParams.update({'font.size': 15})
plt.rcParams.update({'axes.titlesize': 17})
from itertools import cycle
colors=[ 'red','darkorange', 'green', 'blue', 'indigo']
city_list=Historical_results['City'].unique()
#samples=['City_Area', 'Under_Area', 'Under_Ratio']
samples=['City_Area', 'Under_Ratio']
#labels=['City Area in Square Miles', 'Underbounded Area in Square Miles', 'Ratio of Underbounded Area to City Area']
labels=['City Area in Square Miles', 'Ratio of Underbounded Area to City Area']
year_label=['', '', 'Year']
def draw_maps(df, samples, n_rows, n_cols, labels):
    fig=plt.figure(figsize=(14,6))

    for j, sample in enumerate(samples):
        ax=fig.add_subplot(n_rows,n_cols,j+1)
        lines = ["-","--","-.",":", "-.."]
        linecycler = cycle(lines)
        for i, city_selected in enumerate(city_list):
            years=Historical_results[Historical_results['City']==city_selected]['year'].copy()
            variable=Historical_results[Historical_results['City']==city_selected][sample].copy()
        
            plt.plot(years, variable,  next(linecycler), label = city_selected,color=colors[i], linewidth=2.5)
            plt.legend(loc="upper left")
            plt.xlabel(year_label[j])
            plt.ylabel(labels[j])
            plt.xticks(np.arange(1820, 2021, 40))
        
    fig.tight_layout()  
    plt.savefig('Figure_4.png', dpi=600)
    plt.show()
    
draw_maps(Historical_results, samples, 1, 2, labels)


#%%

# =============================================================================
# # =============================================================================
# # Creating plots for methodological appendix
# # =============================================================================
# =============================================================================
Post_UTM=find_UTM(Post_DF_annexed)
#Post_UTM=find_UTM(Post_DF_annexed_selected)

print(Post_UTM.UTM.value_counts())

DF_list_2=[]
UTM_list=list(Post_UTM.UTM.unique())
UTM_list = [int(i) for i in UTM_list] 
UTM_list=[10, 11, 12, 13, 14, 15, 16, 17, 18]

for UTM in UTM_list:
    exec("Post_UTM_{utm}=Post_UTM[Post_UTM['UTM']=={utm}]".format(utm=UTM))
    exec("DF_list_2.append(Post_UTM_{utm})".format(utm=UTM))

#%%

DF_list_1=[Post_UTM_13[(Post_UTM_13['GISJOIN']=="G35002000")].copy()]

#%%

# =============================================================================
# Creating maps of Houston to illustrate method
# =============================================================================

for DF in DF_list_1:
    epsg=DF['EPSG'].iloc[0]
    DF=DF.to_crs(epsg)
    US_place_2010_projected=US_place_2010.to_crs(epsg)
    US_place_2010_newIDX=US_place_2010_projected.reset_index(drop=True)
    US_place_2010_sindex=US_place_2010_newIDX.sindex
    boxes=DF.geometry.bounds
    bounds=[]
    bounds.append(boxes['minx'].min())
    bounds.append(boxes['miny'].min())
    bounds.append(boxes['maxx'].max())
    bounds.append(boxes['maxy'].max())            
    place_candidate_idx = list(US_place_2010_sindex.intersection(bounds))
    place_candidates=US_place_2010_newIDX.loc[place_candidate_idx]
    place_candidates['Other_Place']=place_candidates[['GISJOIN']]
    place_candidates=place_candidates[['Other_Place', 'geometry']].copy()
    place_candidates=place_candidates.set_geometry('geometry')
    city_shapes=place_candidates['geometry']
    
    cities_buffer=city_shapes.buffer(.001)
    boundary = gpd.GeoSeries(cascaded_union(cities_buffer))
    boundary=boundary.loc[0]

    print("finished data preparation")
    if len(DF)>0:
        for i in range(len(DF)):
            
            city_test=DF.iloc[[i]]
            city, city_buffer, city_buffer_reduce, other_cities, GISJOIN= find_shapes(city_test, boundary)
            #GISJOIN=city_test['GISJOIN'].iloc[0]
            #city=city_test['geometry'].iloc[0]
            #city=city.buffer(0)

            box_min = city.minimum_rotated_rectangle
            x, y = box_min.exterior.coords.xy

            edge_length = (Point(x[0], y[0]).distance(Point(x[1], y[1])), Point(x[1], y[1]).distance(Point(x[2], y[2])))
            width=min(edge_length)/2
            city_buffer=city.buffer(width)
            
            city_buffer_reduce=city.buffer(width).buffer(-1*width)
            city_fringe=city.buffer(width/4)
            city_fringe_2=city_fringe.difference(boundary)
            
            unaffected=city_fringe.difference(city)
            unaffected=unaffected.difference(other_cities)

            fringe_other=city_buffer_reduce.difference(city_fringe_2)
            other_city=fringe_other.difference(city.buffer(.01))
 
            avoided=city_buffer_reduce.difference(city)
            avoided=avoided.difference(boundary)
            unaffected=unaffected.difference(avoided)


#%%
city_DF=gpd.GeoDataFrame(city)
city_DF=city_DF.rename(columns={0: 'geometry'})

x, y = box_min.exterior.coords.xy
box_min_clean = Polygon(zip(list(x), list(y)))
box_min_DF = gpd.GeoDataFrame(index=[0], crs=epsg, geometry=[box_min_clean])       
box_min_DF=box_min_DF.rename(columns={0: 'geometry'})

city_buffer_DF=gpd.GeoDataFrame([city_buffer])
city_buffer_DF=city_buffer_DF.rename(columns={0: 'geometry'})

city_buffer_reduce_DF=gpd.GeoDataFrame([city_buffer_reduce])
city_buffer_reduce_DF=city_buffer_reduce_DF.rename(columns={0: 'geometry'})

city_fringe_DF=gpd.GeoDataFrame([city_fringe])
city_fringe_DF=city_fringe_DF.rename(columns={0: 'geometry'})


unaffected_DF=gpd.GeoDataFrame([unaffected]).T
unaffected_DF=unaffected_DF.rename(columns={0: 'geometry'})

other_city_DF=gpd.GeoDataFrame(other_cities)
other_city_DF=other_city_DF.rename(columns={0: 'geometry'})

avoided_DF=gpd.GeoDataFrame(avoided)
avoided_DF=avoided_DF.rename(columns={0: 'geometry'})

part_list=[city_DF, box_min_DF, city_buffer_DF,city_buffer_reduce_DF, 
           city_fringe_DF,  other_city_DF,  unaffected_DF,  avoided_DF ]

def draw_maps(part_list, n_rows, n_cols):
    fig=plt.figure(figsize=(12,9))
    titles=['(1) City of Albuquerque', '(2) Minimum Rotated Rectangle', 
            '(3) City Buffer = .5W','(4) Concave Hull', '(5) City Buffer = .125W', 
            '(6) Other Cities', '(7) Unaffected Areas',  "(8) Underbounded Areas"]
    axis_range=part_list[2]
    axis_range=axis_range.set_crs(epsg)
    axis_range=axis_range['geometry']
    box=axis_range.bounds

    for i, part in enumerate(part_list):
        ax=fig.add_subplot(n_rows,n_cols,i+1)
        part=part.set_crs(epsg)
        part.plot(ax=ax, color='black')
        ax.set_title(titles[i])
        ax.axis('off')
        plt.xlim(box['minx'].min(), box['maxx'].max())
        plt.ylim(box['miny'].min(), box['maxy'].max())
    
    fig.tight_layout()  # Improves appearance a bit.
    plt.savefig('Figure_Method.png', dpi=600)
    plt.show()
    

draw_maps(part_list, 3, 3)
