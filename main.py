#!/usr/bin/env python
# coding: utf-8

# # Data Visualization Exam

# ## Initialization

# ### Imports

# In[ ]:

import geopandas as gpd
import json
import numpy as np
import os
import pandas as pd
import plotly.graph_objs as go
import plotly.express as px
import re
from shapely.wkt import loads
import shutil as sh
import time


# ### Parameters

# In[ ]:

dataFolderName = 'data'
geoJsonFolder = dataFolderName+'/geoJson/'
figureOutputFolder = 'exported_figures'
dataFileName = dataFolderName + '/DCSC_RACLI_01092021113430630.csv' # for data loading (salaries)

default_font_family = "Bahnschrift Light"
colors_palette = ['#003a2b','#249e89','#f5f5f5','#d86e58','#6a0000']

# for animations
animation_duration_frame = 1000         # millisecs
animation_duration_transition = 100     # millisecs

exportFigure = False    # set to true if you want to export the figure



########################### START OF PROGRAM ###########################
startTime = time.time_ns()



# ## Data acquisition
# CSV file creation from raw data and loading the CSV file into the program.

# In[ ]:

print("\n\n## DATA ACQUISITION ##")

dataFolderName = "data"
fileName = dataFolderName + '/DCSC_RACLI_01092021113430630.csv'
df = pd.read_csv(fileName)      # load data from CSV to program
df.head() # data loaded


# ## Data parsing
# In[ ]:

print("\n\n## DATA PARSING ##")
print("DATASET")
print('\tNumber of rows: ' + str(len(df)) + ' rows')
df.drop_duplicates()
print('\tNumber of rows without duplicates: ' + str(len(df)) + ' rows')

# are values reasonable?
print('\tMin value is ' + str(df['Value'].min()) + ' €/h')
print('\tMax value is ' + str(df['Value'].max()) + ' €/h')
print('\tMean value is ' + str(df['Value'].mean()) + ' €/h')


# ### Trasform the data
# In[ ]:

print("\n\n\t### DATA TRANSFORMING ##")
df['Territorio'] = df['Territorio'].str.replace(' / ','/')


# #### Translation
# In[ ]:

# rename sectors in english
it_sec_names = tuple(
    df.query('`Ateco 2007`!="TOTALE" & `ATECO_2007`>="A" & `ATECO_2007`<="Z"')['Ateco 2007'] # Column 'ATECO_2007' has a letter for macro sectors (it is different from 'Ateco 2007'!!
    .drop_duplicates().sort_values().reset_index(drop=True)
)

en_sec_names = (
    "Other service activities",
    "Arts, sports, entertainment and recreation", 
    "Accommodation and food service activities",
    "Financial and insurance activities",
    "Real estate activities",
    "Manufacturing activities",
    "Professional, scientific and technical activities",
    "Wholesale and retail trade, repair of motor vehicles and motorcycles",
    "Constructions",
    "Extraction of minerals from quarries and mines",
    "Water supply sewerage, waste management and sanitation activities",
    "Supply of electricity, gas, steam and air conditioning",
    "Education",
    "Rental, travel agencies, business support services", 
    "Health and social work",
    "Information and communication services", 
    "Transport and storage"
)

df["Ateco 2007 IT"] = df["Ateco 2007"].tolist()  # save the column in Italian

for i in range(0, len(it_sec_names)):
    df.loc[df['Ateco 2007']==it_sec_names[i],"Ateco 2007"] = en_sec_names[i]


# Modify long names to fit two lines instead of one (when used in charts)
df["Ateco 2007 BR"] = df["Ateco 2007"].tolist()  # save the column with long names
long_names = (
    "Financial and insurance activities",
    "Supply of electricity, gas, steam and air conditioning",
    "Extraction of minerals from quarries and mines",
    "Information and communication services",
    "Professional, scientific and technical activities",
    "Arts, sports, entertainment and recreation",
)
br_names = (
    "Financial and insurance  "       +"<br>"+  "activities",
    "Supply of electricity, gas,  "   +"<br>"+  "steam and air conditioning",
    "Extraction of minerals  "        +"<br>"+  "from quarries and mines",
    "Information and  "               +"<br>"+  "      communication services",
    "Professional, scientific  "      +"<br>"+  "and technical activities",
    "Arts, sports, entertainment  "   +"<br>"+  "and recreation",
)
for i in range(0, len(long_names)):
    df.loc[df['Ateco 2007']==long_names[i],"Ateco 2007 BR"] = br_names[i]


print(df.head())


# ## Data filtering
# In[ ]:

print("\n\n## DATA FILTERING ##")

df2 = df.copy()

# unique data
del df2['TIPO_DATO7'] # always the same (HOUWAG_ENTEMP_AV_MI)
del df2['Tipo dato']  # always the same (Retribuzione lorda oraria per ora retribuita delle posizioni lavorative dipendenti in euro (media).)

# ridondance of information
df2 = df2.drop(['SEXISTAT1', 'ETA1_A','PROFILO_PROF','CLLVT','Seleziona periodo'], axis=1)
# del df2['ATECO_2007']

df2 = df2[df2['Flag Codes'] != 'c'] # delete incomplete data

del df2['Flags']
del df2['Flag Codes']

df2.head()


# ## Data mining

# In[ ]:

print("\n\n## DATA MINING ##")


# granularity of sectors exists only for entire Italy (no territorial granularity)
df_sectors = df2.query('`Ateco 2007`!="TOTALE"')

# choose granularity of sectors
df_sectors = df_sectors.query('`ATECO_2007`>="A" & `ATECO_2007`<="Z"').drop(['Territorio', 'ATECO_2007'], axis=1)

# choose granularity of territory
df_territory = df2.query('`Ateco 2007`=="TOTALE"').drop(['Ateco 2007'], axis=1)


# In[ ]:

print("\n\n## Question 1: Where do people earn more? ##")

# ## Question 1

# ### Utility functions

# In[ ]:


# Utility functions to read data from csv and shape files, remove useless columns frome the dataframe and transfrom the data for the program

def loadDataFromCSV(forceUpdate=False):
    '''
    Load data about salaries into the program.
    Returns a Pandas Dataframe.
    If the parameter forceUpdate is set to True, this function will
    reload the dataframe from the file even if it was already loaded
    (to be used when suspecting the data are chenged on the file).
    '''

    df_ref = {}

    def closureFun(forceUpdate=False):

        if (forceUpdate==True):
            df_ref.clear() # clear the df
        
        if (len(df_ref)==0):

            df = pd.read_csv(dataFileName).drop_duplicates()

            # Transform the data and remove useless columns
            df['Territorio'] = df['Territorio'].str.replace(' / ','/')
            df = df.drop('TIPO_DATO7', axis=1) # always the same (HOUWAG_ENTEMP_AV_MI)
            df = df.drop('Tipo dato', axis=1)  # always the same (Retribuzione lorda oraria per ora retribuita delle posizioni lavorative dipendenti in euro (media).)
            df = df.drop(['SEXISTAT1', 'ETA1_A','PROFILO_PROF','CLLVT','Seleziona periodo'], axis=1)  # ridondance of information
            df = df[df['Flag Codes'] != 'c'].drop(['Flags','Flag Codes'], axis=1) # delete incomplete data and drop columns with corresponding flag ('c' is the flag for hidden data)

            # Transform data for consistency with datasets of geocoords
            df.loc[df['Territorio']=="Forlì-Cesena", "Territorio"] = "Forli'-Cesena"

            # Save the dataframe
            df_ref[0] = df

        return df_ref[0]
    
    return closureFun

loadDataFromCSV = loadDataFromCSV() # use the closure
    


def getDataAboutTerritory():
    '''
    Returns data about salaries in territories (data about sectors are excluded).
    '''
    return loadDataFromCSV().query('`Ateco 2007`=="TOTALE"').drop(['Ateco 2007', 'ATECO_2007'], axis=1)


def getDataAboutProvinces():
    '''
    Returns data about salaries in provinces (data about sectors, regions, entire Italy are excluded).
    '''
    df_territory = getDataAboutTerritory()
    years = loadDataFromCSV()['TIME'].drop_duplicates()

    # Note: this column is present also in geo-data and can be used to join the datasets
    df_territory["TerritorioAnno"] = df_territory["Territorio"] + df_territory['TIME'].astype(str)
    return df_territory[df_territory['ITTER107'].str.contains('.{5}')].drop('ITTER107', axis=1)   # for provinces, 'ITTER107' code is 5 chars long


def getDataAboutProvincesInDictHavingYearsAsKey(years=-1):
    '''
    Returns data about salaries in provinces (data about sectors, regions, entire Italy are excluded),
    organized in a dictionary having years (the parameters) as keys.
    Params: years, e.g.: years=range(2014,2018).
    If the parameter years is not specified, all the years are considered.
    '''
    dataProvinces = getDataAboutProvinces()
    if(years==-1):
        years = dataProvinces['TIME'].drop_duplicates()

    return {year: dataProvinces.query(f'TIME=={year}').drop_duplicates() for year in years}


def getProvinceSalaryvalue(year=-1):        # TODO: take a list as input parameter
    '''
    Returns a Pandas Dataframe with three columns: one for Province names ("Territorio"), the second for the
    year ("TIME") and the third for the corresponding salary value ("Value"); column names are the ones inside
    the brackets ("Territorio", "TIME", "Value").
    Returned data refer to the year which is given as parameter.
    If the year parameter is not specified, also the column 'TIME' is returned, with the corresponding year
    '''
    df_years = getDataAboutProvincesInDictHavingYearsAsKey([year]) if year!=-1                                                              else getDataAboutProvincesInDictHavingYearsAsKey()
    
    years = sorted(df_years.keys())

    df_years = {year: df_years[year].query("Sesso=='totale' & `Classe di età`=='totale' & `Qualifica contrattuale`=='totale' & `Classe di dipendenti`=='totale'")
                                    .drop(['Sesso', 'Classe di età', 'Qualifica contrattuale', 'Classe di dipendenti'], axis=1)
                    for year in years}

    # Categorization of Salary values (grouping in categories)
    valueCountedData = {year: np.floor(df_years[year]["Value"]).astype(int).value_counts() for year in years}

    # NOTE: This part should be part of data transforming? But ranges should adapt to the context?

    salaryCategoryBorders = range(11,20,2)   # same category subdivion for all years
    for year in years:
        oldCategory=0
        df = df_years[year]
        for category in salaryCategoryBorders:
            numberProvinceInThisCategory = sum([valueCountedData[year][key] for key in np.intersect1d(valueCountedData[year].keys().tolist(), range(oldCategory,category))])
            df.loc[(oldCategory<=df['Value']) & ( (df['Value']<category) | (df['Value']>=salaryCategoryBorders[-1]) ), "SalaryCategory"] =                                   (f"{oldCategory} ≤ " if oldCategory >= salaryCategoryBorders[0] else "       ")                                                                              + ".."                                                                                                                                                       + (f" < {category}"  if category < salaryCategoryBorders[-1] else "       ")                                                                                 + f"\t  €/h\t({numberProvinceInThisCategory} province{'s' if numberProvinceInThisCategory>2 else ''})"
            oldCategory = category
        
        # sort (needed to respect the range-scale in plots if categorization is used)
        df.sort_values(by=['Value'], ascending=True, inplace=True)
        
        df_years[year] = df

    df = pd.concat(tuple(df_years[year] for year in years))

    # sort (needed to respect the range-scale in plots if categorization is used)
    #   Sort (first) ascending wrt 'TIME' (oldest first) then descending wrt 'Value'
    df['Value'] = -df['Value']  # invert sign, so 'Value' can be sorted descending
    df.sort_values(by=['TIME', 'Value'], ascending=True, inplace=True)
    df['Value'] = -df['Value']  # restore the correct sign
    
    return  df


def categorization(salaryCategoryBorders = range(11,20,2)): # TODO: code duplication with getProvinceSalaryvalue()
    '''
    Returns a Pandas Dataframe with three columns: one for the year, the second for the
    salary category and the third for the corresponding number of provinces where people
    earn as much as declared in the category.
    You can specifiy the range for the categories as parameter.
    '''

    # TODO : REFACTORING (code duplication with the previous function)

    df_years = getDataAboutProvincesInDictHavingYearsAsKey()
    years = sorted(df_years.keys())
    df_years = {year: df_years[year].query("Sesso=='totale' & `Classe di età`=='totale' & `Qualifica contrattuale`=='totale' & `Classe di dipendenti`=='totale'")
                                    .drop(['Sesso', 'Classe di età', 'Qualifica contrattuale', 'Classe di dipendenti'], axis=1)
                    for year in years}

    # Categorization of Salary values (grouping in categories)
    valueCountedData = {year: np.floor(df_years[year]["Value"]).astype(int).value_counts() for year in years}

    df_toReturn = pd.DataFrame(columns=['Year', 'Gross salary  [€/h]', '#Provinces'])
    columnNames = tuple(df_toReturn.columns)

    for year in years:
        oldCategory=0
        for category in salaryCategoryBorders:
            numberProvinceInThisCategory = sum([valueCountedData[year][key] for key in np.intersect1d(valueCountedData[year].keys().tolist(), range(oldCategory,category))])
            categoryStr = (f"{oldCategory} ≤ " if oldCategory >= salaryCategoryBorders[0] else "        ") + ".."                                + (f" < {category}"  if category < salaryCategoryBorders[-1] else "        ")
            oldCategory = category
            df_toReturn = df_toReturn.append({columnNames[0]: year, columnNames[1]: categoryStr, columnNames[2]: numberProvinceInThisCategory}, ignore_index=True)


        
        # sort (needed to respect the range-scale in plots if categorization is used)
        df_toReturn.sort_values(by=[columnNames[0]], ascending=True, inplace=True)
    
    return  df_toReturn


def avgSalary(territory='Italia', year=-1):
    '''
    Returns the average salary value in a given territory for a given year (parameters).
    The default value for the territory is entire Italy.
    If the year is not specified, the average value is computed over all the years which
    are available.
    '''
    query = f"Territorio=='Italia' & Sesso=='totale' & `Classe di età`=='totale' & `Qualifica contrattuale`=='totale' & `Classe di dipendenti`=='totale'"               + (f" & `TIME=={year}" if year!=-1 else "")
    return round(100*getDataAboutTerritory().query(query)['Value'].mean())/100  # round(100*..)/100 is used to have two decimal digits


# Utility functions for geo-data
def readGeoDataToDictHavingYearAsKey():
    '''
    Import data Geo-data (coordinates) and returns the dictionary having as key
    the year and as values the dataframe with geodata loaded from shape files.
    '''
    map_df = {} # dictionary, year as key
    map_df[2014] = gpd.read_file(f'{dataFolderName}/province_shapes/Prov01012014_g/Prov01012014_g_WGS84.shp')
    map_df[2014]['DEN_PCM'] = map_df[2014]['DEN_PROV']  # duplicate this column to make the dataframe compliant with those of subsequent years 
    map_df[2014].loc[ map_df[2014].DEN_PCM=="Forlì-Cesena","DEN_PCM" ] = "Forli'-Cesena"

    for year in range(2015,2018):
        fp = f'{dataFolderName}/province_shapes/ProvCM01012017_g/ProvCM01012017_g_WGS84.shp' # data updated to 1st Jan 2017 work for our purposes
        map_df[year] = gpd.read_file(fp) #reading the file stored in variable fp
        map_df[year].loc[ map_df[year].DEN_PCM=="Aosta","DEN_PCM" ] = "Valle d'Aosta/Vallée d'Aoste"
        map_df[year].loc[ map_df[year].DEN_PCM=="Massa Carrara","DEN_PCM" ] = "Massa-Carrara"
        map_df[year].loc[ map_df[year].DEN_PCM=="Bolzano","DEN_PCM" ] = "Bolzano/Bozen"

    # Note: territories coords change over the year, hence we save the year near the territory names
    for year in map_df.keys():
        map_df[year]["TerritorioAnno"] = map_df[year]["DEN_PCM"] + str(year)
        map_df[year] = map_df[year][['DEN_PCM','TerritorioAnno','geometry']]
    
    return map_df


# Function to convert (project) coordinates to latitude/longitude
def convertCrsToLatLong(inputGeopandasDf, inplace=False):
    '''
    Convert the geo-coordinates of the iunput GeoPandas Dataframe to EPSG:4326 (latitude and longitude)
    and returns a new GeoPandas dataframe having the data in the new coordinates system.
    You can specify the parameter inplace=True if you want to change the coordinate system "inplace",
    i.e., directly in the input GeoPandas Dataframe.
    '''
    outputGeopandasDf = inputGeopandasDf.set_geometry("geometry") # The original geometry column is replaced with "geometry" (if it was different).
    outputGeopandasDf = outputGeopandasDf.to_crs("EPSG:4326", inplace=inplace)
    return outputGeopandasDf
    

def createGeoJsonFromFile(geoJsonFolder, shapeDataDictYears, convertCrsToLatLongFlag=True):
    '''
    Creates GeoJson files in the folder whose path is specified as parameter as string,
    from the given dictionary having years as keys and the corresponding shape file data
    (GeoPandas dataframe) as values.
    The parameter shapeDataDictYears can also be the shape file data directly, i.e. the
    value of onlyh one record of a dictionary.
    Specify the parameter convertCrsToLatLongFlag=False if you do NOT want to convert the
    geo-coordinate system to EPSG:4326; default is True.
    Returns a dictionary having as key the years (the same as the input dictionary) and
    the corresponding GeoJson data as values.
    '''
    geoJsonData = {}
    if not os.path.exists(geoJsonFolder):
        os.makedirs(geoJsonFolder)              # TODO : check for issues (everything correct? Warning: '"writeGeoJson" is not accessed', as if os.makedirs was never used)

    isInputShapeDataAsDict = type(shapeDataDictYears) is dict # true id a dictionary is given as input parameter
    if(not isInputShapeDataAsDict):
        shapeDataDictYears = {'': shapeDataDictYears}    # converted to dict to use the same code

    for year in shapeDataDictYears.keys():
        if(convertCrsToLatLongFlag):
            shapeDataDictYears[year] = convertCrsToLatLong(shapeDataDictYears[year])
        geoJsonPathThisYear = geoJsonFolder+str(year)+'.json'
        shapeDataDictYears[year].to_file(geoJsonPathThisYear, driver="GeoJSON")
        with open(geoJsonPathThisYear, encoding="utf-8") as geofile:
            geoJsonData[year] = json.load(geofile)    
    
    return geoJsonData if(isInputShapeDataAsDict)                        else geoJsonData[[v for v in shapeDataDictYears.keys()][0]]


def loadDataMultipleYears(provinceNames=[], years=[], compress=-1, simplify=-1):
    '''
    Returns the GeoJson data and the dataframe of provinces (only with territories, economic sectors
    excluded) for all the years. Additionally, the dataframe of provinces is improved with percentages
    of hourly average gross salary of the province w.r.t. the national hourly average gross salary
    (named "Salary wrt. national average [%]") and with the column named "Value2" with salary values
    follwed by the measure unit (€/h), for printing purposes.
    The two dataframes (geoJsonData, df_province) have to be unpacked.
    This function can be used to rapidly load both geo-data and data about salaries in provinces, over
    all the years (province granularity only).
    If the parameter provinceNames is specified, only data about the desired provinces will be loaded
    (a list is expected).
    If the parameter years (a list is expected) is specified, only data about selected years will be
    returned.
    If the parameter 'compress' is specified and set to a positive value, then a compressed version of the GeoJson
    data will be provided. The compression is given by rounding the precision of the geo coordinates to
    the specified number of decimal digits.
    Similarly, you can specify a tolerance value for the parameter 'simplify'.
    See: https://geopandas.org/docs/user_guide/geometric_manipulations.html#GeoSeries.simplify
    '''

    # Read geo-data
    map_df = readGeoDataToDictHavingYearAsKey() # dictionary, year as key

    # Load data about salaries for each province
    df_province = getProvinceSalaryvalue()

    if(len(years)>0):   # filter according to years
        df_province = df_province.query(' | '.join({f"(TIME=={year})" for year in years}))
        map_df = {year: map_df[year] for year in years}
    else:
        years = map_df.keys()

    if(len(provinceNames)>0):
        df_province = df_province.query(' | '.join({f'(Territorio=="{provinceName}")' for provinceName in provinceNames}))
        map_df = {year: map_df[year].query(' | '.join({f'(DEN_PCM=="{provinceName}")' for provinceName in provinceNames})) for year in years}
    
    # Compute percentage of salaries in each province wrt. the national average value and add to the dataframe
    for year in years:
        # Percentage increment (I) for a province wrt. national average value (A):      V = A + I/100*A ,   V=value in the province, ==> I = 100(V/A-1)  [%]
        nationalAvgSalary = df_province.query(f'TIME=={year}')['Value'].mean()
        df_province.loc[df_province.TIME==year, "Salary wrt. national average [%]"] = round(100*(df_province.loc[df_province.TIME==year,'Value']/nationalAvgSalary-1), 2)

    # Format the percentage salary
    df_province["Salary wrt. national average [%]"] = df_province["Salary wrt. national average [%]"].map('{:+.2f} %'.format)   # Same as: # df_province["Salary wrt. national average [%]"] = df_province["Salary wrt. national average [%]"].map(lambda val: ('+' if val>0 else '') + str(val))

    # Format the value
    df_province["Value2"] = df_province["Value"].map('{:.2f} €/h'.format)

    # Union over years of geodata and conversion of coordinates
    geoData = pd.concat(tuple(convertCrsToLatLong(map_df[year]) for year in years))

    # Compression of geo data (from: https://gis.stackexchange.com/a/321531)
    if compress>=0:
        # Round coordinates to the specified number of decimal digits. Topology may not be preserved
        simpledec = re.compile(r"\d*\.\d+")
        geoData['geometry'] = geoData['geometry'].apply(lambda x: loads(re.sub(simpledec, lambda match: f"{float(match.group()):.{compress}f}", x.wkt)))                                                     .simplify(0) # 0 means no tolerance
    if simplify>0:
        geoData['geometry'] = geoData['geometry'].simplify(simplify)
    
    # Create GeoJson from SHP dataframe (union over years of shp files)
    geoJsonData = createGeoJsonFromFile(geoJsonFolder, geoData)
    return geoJsonData, df_province


def createFigure(dataframe, geoJsonData):
    '''
    Create the figure for answering to this question.
    Parameters:
     -  dataframe   :  the dataframe with values to use
     -  geoJsonData :  the geoJson dataframe associated with the given dataframe
     If data refer to more than one year, an animation over the years is shown.
     Returns: the created figure.
    '''

    # consider only a subset of columns (less size)
    dataframe = dataframe[['TIME', 'TerritorioAnno', 'Territorio', 'Value2', 'SalaryCategory', 'Salary wrt. national average [%]']]

    showAnimationFlag = len(dataframe['TIME'].drop_duplicates()) > 1        # in order to show the animation over the years, data about more than one year must be available

    if showAnimationFlag:
        # Keep only the salary category (drop out the number of provinces belonging to it)
        dataframe['SalaryCategory2'] = tuple( aMatch[0] for aMatch in re.findall(r"((\s*[0-9]*\s*([≤][ ])?[.]{2}([ ][<])?\s*[0-9]*)(\s*€\/h))", ''.join(dataframe['SalaryCategory'].tolist()) ) )
        salaryCategories = dataframe['SalaryCategory2'].drop_duplicates().sort_values().tolist()
        salaryCategories = tuple([salaryCategories[-1]] + salaryCategories[:-1])
        if len(salaryCategories)!=len(colors_palette):
            raise Exception('Number of colors is different than the number of categories')

    colorLabel = 'SalaryCategory2' if showAnimationFlag else 'SalaryCategory'
    fig = px.choropleth(
        # title=None, title='Salaries in private companies',
        data_frame=dataframe, 
        geojson=geoJsonData, 
        locations='TerritorioAnno',               # name of dataframe column
        hover_name='Territorio',
        hover_data={'TIME':False, 'Value2':True, colorLabel:False, 'Territorio': False, 'TerritorioAnno': False, 'Salary wrt. national average [%]': True},          # TODO: improve this (see "hovertemplate")
        featureidkey='properties.TerritorioAnno', # path to field in GeoJSON feature object with which to match the values passed in to locations
        color=colorLabel,
        color_discrete_sequence=colors_palette,   # for discrete scale of colors
        center={"lat": 42, "lon": 13},
        projection='mercator',
        labels={
            colorLabel: '<br><br>Average hourly gross salary                              ',    # <br> and spaces used for layout
            'Value2': 'Avg salary', "Salary wrt. national average [%]": "Percentage"},
        animation_frame="TIME" if showAnimationFlag else None,
    )
    fig.update_traces(marker=dict(opacity=1, line=dict(color='black', width=0.1)))      # TODO: look for "hovertemplate, https://plotly.com/python/reference/choropleth/#choropleth-hovertemplate"
    fig.update_layout(        
        hoverlabel=dict(font_family=default_font_family),
        plot_bgcolor='white',
        font=dict(color='dimgray', family=default_font_family),
        # margin={"r":0,"t":0,"l":0,"b":0},
        title_font_family=default_font_family,
        legend_itemsizing='trace'               # Determines if the legend items symbols scale with their corresponding "trace" attributes or remain "constant" independent of the symbol size on the graph. # TODO: NOT working
    )

    # Add text annotation outside the map
    fig.add_annotation(
        dict(font=dict(color="dimgray",size=10),
            x=1.225 if showAnimationFlag else 1.67,
            y=0.45,
            showarrow=False,
            text='Percentages on hover labels refer to the <br>national average gross salary for the year.',
            textangle=0,
            align='left',
            xref="paper",
            yref="paper"
        )
    )

    fig.update_geos(showcountries=False, showcoastlines=False, showland=False, fitbounds="locations")

    if showAnimationFlag:
        fig.layout.updatemenus[0].buttons[0].args[1]['frame']['duration'] = animation_duration_frame
        fig.layout.updatemenus[0].buttons[0].args[1]['transition']['duration'] = animation_duration_transition

    return fig


# ### Plot maps

# In[ ]:



# Load both geodata and data about salaries (only for the desired province), for each year

salaryCategories = range(11,20,2)
df_year_salaryCategory_nProvs = categorization(salaryCategories)
df_year_salaryCategory_nProvs.iloc[:,1] = "    " + df_year_salaryCategory_nProvs.iloc[:,1] + "   " # update the text
df_year_salaryCategory_nProvs = {year: df_year_salaryCategory_nProvs.query(f"{df_year_salaryCategory_nProvs.columns[0]}=={year}").iloc[:,1:3].sort_values(by=[df_year_salaryCategory_nProvs.columns[1]], ascending=True) for year in df_year_salaryCategory_nProvs.iloc[:,0].drop_duplicates()}  # create a dictionary {year: [salaryCategory, numberOfProvinces]}
print(df_year_salaryCategory_nProvs[2014])

years = years = df_year_salaryCategory_nProvs.keys()
max_x_val = max(df_year_salaryCategory_nProvs[year].iloc[:,1].max() for year in years)


if exportFigure:
    figureOutputFolder_this = figureOutputFolder + '/question1'
    if os.path.exists(figureOutputFolder_this): # remove old data
        sh.rmtree(figureOutputFolder_this)
    os.makedirs(figureOutputFolder_this)
    

# Load data and compressed geodata about salaries (for animation over years)
simplifyTolerance=0.01       # TODO : should be a parameter?
geoJsonData, df_province = loadDataMultipleYears(simplify=simplifyTolerance) # simplify geodata (to save memory space)
# decimalDigitsCompression = 2 # TODO : should be a parameter?
# geoJsonData, df_province = loadDataMultipleYears(compress=decimalDigitsCompression) # simplify geodata (to save memory space)
fig = createFigure(df_province, geoJsonData)

# Export the figure
if exportFigure:
    fig.write_html(f"{figureOutputFolder_this}/geoMapSlider.html")


# Re-load data and load uncompressed geodata
geoJsonData, df_province = loadDataMultipleYears()
    
for year in years:
    geoJsonData, df_province = loadDataMultipleYears(years=[year])

    print(f"\n\nYear: {year}")
    maxSalary = max(df_province['Value'])
    minSalary = min(df_province['Value'])
    best_province  = df_province.query(f"Value=={maxSalary}")
    worst_province = df_province.query(f"Value=={minSalary}")
    if(len(best_province)>1 or len(worst_province)>1):
        print("WARNING: query returned more than one result, only the first result is showed")
    
    print(f"\tBest province{'s' if len(best_province)>1 else ''}:\t{str(best_province.to_dict('records'))}")
    print(f"\tWorst province{'s' if len(worst_province)>1 else ''}:\t{str(worst_province.to_dict('records'))}")
    print(f"\tItalian average gross salary: %.2f €/h" % df_province['Value'].mean())

    # Choropleth by categories
    fig = createFigure(df_province, geoJsonData)

    # Export the figure
    if exportFigure:
        fig.write_image(f"{figureOutputFolder_this}/geoMap{year}.svg")

    

    # Barplot

    fig = px.bar(
        data_frame = df_year_salaryCategory_nProvs[year],
        x = df_year_salaryCategory_nProvs[year].columns[1],
        y = df_year_salaryCategory_nProvs[year].columns[0],
        orientation = 'h',  # horizontal bar chart
        text=df_year_salaryCategory_nProvs[year].columns[1],
        height=300,
        width=450,
        # log_x=True  # logarithmic scale
    )

    fig.update_traces(
        marker_color=[color for color in reversed(colors_palette)],
        marker_line_color='black', marker_line_width=1, opacity=1,
        # texttemplate='%{text:d} ', textposition='inside'
    )

    fig.update_layout(
        hoverlabel=dict(font_family=default_font_family),
        title_text=f'{year}',
        yaxis_title=df_year_salaryCategory_nProvs[year].columns[0],
        xaxis_title="Number of provinces",
        xaxis=dict(showline=True, showticklabels=True, ticks='outside',
            linecolor='rgb(204, 204, 204)', linewidth=2, dtick = 10,
            range = [max_x_val, 0]),  # reversed xaxis
            # range = [2, 0]),  # reversed xaxis if log xaxis
        yaxis=dict( showgrid=False, showline=False, side='right'),              # yaxis on the right side
        paper_bgcolor='white',
        plot_bgcolor='white',
        title_font_family=default_font_family,
        font=dict(family=default_font_family),
        showlegend=False,
        hovermode=False
    )
    fig.update_xaxes(title_font_family=default_font_family)
    fig.update_yaxes(title_font_family=default_font_family)


    # Export the figure
    if exportFigure:
        fig.write_image(f"{figureOutputFolder_this}/legend_barChartSectors{year}.svg")


if exportFigure:
    del figureOutputFolder_this


# ## Question 2

# ### Filter dataset for q2

# In[ ]:

print("\n\n## Question 2: Has the gender gap decreased over time? ##")

df_sex = df_territory.query('Sesso!="totale" & Territorio=="Italia"')[['Sesso','TIME','Value']]


# ### Plot line chart

# In[ ]:


labels = ['Male','Female','Gap']
colors = ['#5b8592','#fbb4b9','#171717']
first_year = df_sex.TIME.min()
last_year = df_sex.TIME.max()

x_year = np.arange(first_year,last_year+1)
x_data = np.vstack((x_year,)*3)

df_sex.sort_values(by='TIME')
df_mal = df_sex.query('Sesso=="maschi"')['Value'].to_list()
df_fem = df_sex.query('Sesso=="femmine"')['Value'].to_list()
gap = tuple(round(df_mal[i] - df_fem[i],2) for i in range(0,last_year-first_year+1))
gapPercentageWrtFirtYear = [str(round((gapValThisYear-gap[0])/gap[0]*100, 2)) + ' %' for gapValThisYear in gap]
gapPercentageWrtFirtYear[0] = 'ref.'
y_data = np.array([df_mal,df_fem,gap])

fig = go.Figure()

annotations = []

for i in range(0, len(labels)):
    fig.add_trace(go.Scatter(x=x_data[i], y=y_data[i], mode='lines',
        name=labels[i], line=dict(color=colors[i]), connectgaps=True))
    # endpoints
    if i==0: # male
        fig.add_trace(
            go.Scatter(x=x_data[i], y=y_data[i],
            mode='markers+text', marker=dict(color=colors[i]),
            text=y_data[i] , textposition="top center")
        )
    elif i==1: # female
        fig.add_trace(
            go.Scatter(x=x_data[i], y=y_data[i],
            mode='markers+text', marker=dict(color=colors[i]),
            text=y_data[i] , textposition="bottom center")
        )
    else:
        fig.add_trace(
            go.Scatter(x=x_data[i], y=y_data[i],
            mode='markers+text', marker=dict(color=colors[i]),
            text=y_data[i], textposition="bottom center")
        )
        
    
    # Name of lines
    annotations.append(dict(text=labels[i],showarrow=False,
        xref='x', x=x_data[i,3]+0.05, y=y_data[i,3], xanchor='left', yanchor='middle', 
        font=dict(family=default_font_family,size=16,color=colors[i])))

fig.update_layout(annotations=annotations)

fig.update_layout(
    xaxis_title="Year",
    yaxis_title="Gross salary [€/h]",
    xaxis=dict( showgrid=False,showline=True, showticklabels=True, ticks='outside',
        linecolor='rgb(204, 204, 204)', linewidth=2, dtick = 1),
    yaxis=dict( showgrid=False,showline=True, showticklabels=True, ticks='outside', 
        linecolor='rgb(204, 204, 204)', linewidth=2, dtick = 2,
        # gridcolor = "lightgrey",
        range = [0, round(max(y_data[0])+1)]),
    showlegend=False,
    plot_bgcolor='white',
    title_font_family=default_font_family,
    font=dict(family=default_font_family,size=10,color="grey"),
    width=800, height=500
)

fig.update_xaxes(title_font_family=default_font_family)
fig.update_yaxes(title_font_family=default_font_family)

# Add shapes
color_shp = colors_palette[2]

fig.update_layout(        
    hoverlabel=dict(font_family=default_font_family),
    shapes=[
        # male-female
        dict(
            type="path",
            path=" M "  +str(x_data[0][0])+","+str(y_data[0][0])+
                    " L"+str(x_data[0][1])+","+str(y_data[0][1])+
                    " L"+str(x_data[0][2])+","+str(y_data[0][2])+
                    " L"+str(x_data[0][3])+","+str(y_data[0][3])+
                    " L"+str(x_data[1][3])+","+str(y_data[1][3])+
                    " L"+str(x_data[1][2])+","+str(y_data[1][2])+
                    " L"+str(x_data[1][1])+","+str(y_data[1][1])+
                    " L"+str(x_data[1][0])+","+str(y_data[1][0])+" Z",
            fillcolor=color_shp,
            line_width=0, 
            layer="below"
        ),

        # gap
        dict(
            type="path",
            path=" M "  +str(x_data[2][0])+","+str(y_data[2][0])+
                    " L"+str(x_data[2][1])+","+str(y_data[2][1])+
                    " L"+str(x_data[2][2])+","+str(y_data[2][2])+
                    " L"+str(x_data[2][3])+","+str(y_data[2][3])+
                    " L"+str(x_data[2][3])+",0"
                    " L"+str(x_data[2][2])+",0"
                    " L"+str(x_data[2][1])+",0"
                    " L"+str(x_data[2][0])+",0 Z",
            fillcolor=color_shp,
            line_width=0, 
            layer="below"
        ),

        # horizontal line to emphasize the gap decreased
        dict(
            type="line",
            x0=first_year, y0=gap[0], x1=last_year, y1=gap[0],
            line=dict(color="gray", width=2, dash="dash"),
            opacity=1, layer="above"
        ),
    ]
)


# Export images
if exportFigure:
    figureOutputFolder_this = figureOutputFolder + '/question2'
    if os.path.exists(figureOutputFolder_this): # remove old data
        sh.rmtree(figureOutputFolder_this)
    os.makedirs(figureOutputFolder_this)
    fig.write_image(f"{figureOutputFolder_this}/genderGapLine.svg")
    del figureOutputFolder_this

print("Percentage increment of gender gap wrt. 2014: " + str(gapPercentageWrtFirtYear))


# ## Question 3

# ### Filter dataset for q3

# In[ ]:

print("\n\n## Question 3: What are the most profitable sectors? ##")

df_sectors_tot = df_sectors.query('Sesso=="totale" & `Classe di età`=="totale" &  \
                                  `Classe di dipendenti`=="totale" & `Qualifica contrattuale`=="totale"'
                            )[['Ateco 2007','Ateco 2007 BR','TIME','Value']]


# ### Plot horizontal bar chart for sectors

# In[ ]:

howMany=5
val_x_axis = max(df_sectors_tot['Value'])

if exportFigure:
    figureOutputFolder_this = figureOutputFolder + '/question3'
    if os.path.exists(figureOutputFolder_this): # remove old data
        sh.rmtree(figureOutputFolder_this)
    os.makedirs(figureOutputFolder_this)

for year in range(2014,2018,1):
    tmp = df_sectors_tot.query(f'TIME=={year}').sort_values(by='Value')
    print(f"\n\tYear {year}: \n\t\t" + "\n\t\t".join(reversed((tmp['Ateco 2007'] + "\t( " + tmp['Value'].astype(str) + " €/h )").tolist()) ) )

    fig = px.bar(tmp.tail(howMany), x="Value", y="Ateco 2007 BR", text="Value")

    fig.update_traces(
        texttemplate='%{text:.2f} ', textposition='inside',
        marker_color= colors_palette[1], opacity=0.8
    )
    fig.update_traces(texttemplate='%{text:.2f} ', textposition='inside')
    fig.update_traces(marker_color= colors_palette[1], opacity=0.8)
    fig.update_layout(
        hoverlabel=dict(font_family=default_font_family),
        #title_text=f'{year}',
        yaxis_title=None,
        xaxis_title="Gross salary [€/h]",
        xaxis=dict(showline=True, showticklabels=True, ticks='outside',
            linecolor='rgb(204, 204, 204)', linewidth=2, dtick = 5,
            range = [0, val_x_axis]),
        yaxis=dict( showgrid=False, showline=False, ticksuffix='  '),
        paper_bgcolor='white',
        plot_bgcolor='white',
        font=dict(family=default_font_family,size=12,color="grey"),
        title_font_family=default_font_family,
        showlegend=False,
        width=800, height=350
    )
    fig.update_xaxes(title_font_family=default_font_family)
    fig.update_yaxes(title_font_family=default_font_family)

    avg = round(np.average(tmp["Value"]),2)
    fig.add_shape(
        type="line",
        x0=avg, y0=-0.5, x1=avg, y1=4.5,
        line=dict(color="grey",width=2),
        opacity=0.5, layer="below"
    )
    fig.add_annotation(
        x=avg, y=4.9,
        text="Average<br>{avg}",
        font=dict(family=default_font_family,size=12,color="grey"),
        showarrow=False
    )

    # Export images
    if exportFigure:
        fig.write_image(f"{figureOutputFolder_this}/barChartSectors{year}.svg")


# ### Plot with slider

# In[ ]:

df_new = pd.DataFrame(columns=['Ateco 2007','TIME','Value'])

for year in years:
    tmp = df_sectors_tot.query(f'TIME=={year}').sort_values(by='Value',ascending=False)
    df_new = df_new.append(tmp.head(howMany))

    if year==2014:
        avg_2014 = round(np.average(tmp["Value"]),2)

df_new = df_new.sort_values(by=['TIME','Value']).reset_index(drop=True)

fig = px.bar(df_new, x="Value", y="Ateco 2007", text="Value",
    animation_frame="TIME", range_x=[0,df_new['Value'].max()*1.1],
    color_discrete_sequence=[colors_palette[1]]*howMany
)

fig.update_traces(texttemplate='%{text:.2f} ', textposition='inside')

fig.update_layout(        
    hoverlabel=dict(font_family=default_font_family),
    #title_text=f'{year}',
    yaxis_title=None,
    xaxis_title="Gross salary [€/h]",
    xaxis=dict(showline=True, showticklabels=True, ticks='outside',
                linecolor='rgb(204, 204, 204)', linewidth=2, dtick = 5,
                range = [0, val_x_axis]),
    yaxis=dict( showgrid=False, showline=False, ticksuffix='  '),
    paper_bgcolor='white',
    plot_bgcolor='white',
    title_font_family=default_font_family,
    font=dict(family=default_font_family,size=12,color="grey"),
    showlegend=False,
    annotations = [
        dict(
            x=avg_2014, y=4.9,
            text=f"Average<br>{avg_2014}",
            font=dict(family=default_font_family,size=12,color="grey"),
            showarrow=False
        )
    ],
    shapes = [
        dict(
            type="line",
            x0=avg_2014, y0=-0.5, x1=avg_2014, y1=4.5,
            line=dict(color="grey",width=2),
            opacity=0.5, layer="below"
        )
    ]
)
fig.update_xaxes(title_font_family=default_font_family)
fig.update_yaxes(title_font_family=default_font_family)

for step in fig.layout.sliders[0].steps:
    step["args"][1]["frame"]["redraw"] = True

for k in range(len(fig.frames)):
    tmp = df_sectors_tot.query(f'TIME=={2014+k}')
    avg = round(np.average(tmp["Value"]),2)

    annotation = [
        dict(
            x=avg, y=4.9,
            text=f"Average<br>{avg}",
            font=dict(family=default_font_family,size=12,color="grey"),
            showarrow=False
        )
    ]
    shape = [
        dict(
            type="line",
            x0=avg, y0=-0.5, x1=avg, y1=4.5,
            line=dict(color="grey",width=2),
            opacity=0.5, layer="below"
        )
    ]

    fig.frames[k]['layout'].update(annotations=annotation,shapes=shape)

fig.layout.updatemenus[0].buttons[0].args[1]['frame']['duration'] = animation_duration_frame

if exportFigure: 
    fig.write_html(f"{figureOutputFolder_this}/barChartSectors.html")
    del figureOutputFolder_this



endTime = time.time_ns()
print("\n\nTime elapsed: " + str((endTime-startTime)/1000000) + " ms")
print("\nEND OF PROGRAM")