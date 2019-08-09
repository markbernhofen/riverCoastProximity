########################################################################
# analysis.py
#
#   This code takes all the files produced in the riverBuffer.py and
#   burnCoastline.py codes and other input files and produces a CSV
#   file that reports for each country: The count of urban population
#   within 1km of a river or coast, the count of rural population within
#   1km of a river or coast, the count of urban population not within
#   1km of a river or coast, the count of rural population not within
#   1 km of a river or coast, and the total population for the country.
#
#   NOTE: can also be a distance other than 1km but that is dependent on
#   the river buffer that was chosen in riverBuffer.py
#
#   NOTE: only argument required when running code is to specify which
#   'mega' continent the analysis should be carried out on.
#
#   INPUTS: river buffer: 'buffer.tif'; population: 'pop.tif'; urban density:
#   'density.tif'; national identifier: 'nationID.tif'; accumulation file:
#   'acc.tif'
#
#   OUTPUTS: CSV containing results: 'results.csv'
#
# Author: Mark Bernhofen (2019) cn13mvb@leeds.ac.uk
########################################################################


import argparse
import os
import sys
import numpy as np
import pandas as pd
from osgeo import gdal
import csv

class analysis:

    def __init__(self, inputDir, overwrite, continent):
        self.inputDir = inputDir
        if not os.path.exists(self.inputDir):
            print("Input directory doesn't exist")
            sys.exit(-1)
        self.overwrite = overwrite
        self.continent = str(continent)

        # Define input and output filenames
        self.bufferFile = os.path.join(self.inputDir, "buffer.tif")
        self.popFile = os.path.join(self.inputDir, "pop.tif")
        self.densityfile = os.path.join(self.inputDir, "density.tif")
        self.nationID = os.path.join(self.inputDir, "nationID.tif")
        self.acc = os.path.join(self.inputDir, "acc.tif")
        self.resultsFile = os.path.join(self.inputDir, "results.csv")

        # List of NationIDs for each continent

        self.Africa = [108, 174, 262, 231, 232, 404, 450, 454, 480, 174, 508, 638, 646, 690,
                  706, 728, 800, 834, 894, 716, 24, 120, 140, 148, 178, 180, 226, 266, 678,
                  12, 818, 434, 504, 729, 788, 732, 72, 426, 516, 710, 748, 204, 854, 132, 384,
                  270, 288, 324, 624, 430, 466, 478, 562, 566, 654, 686, 694, 768, 724, 620, 887]
        self.Eurasia = [156, 344, 446, 158, 408, 392, 496, 410, 398, 417, 762, 795, 860, 4, 50, 64,
                   356, 364, 462, 524, 586, 144, 96, 116, 360, 418, 458, 104, 608, 702, 764,
                   626, 704, 51, 31, 48, 196, 268, 368, 376, 400, 414, 422, 512, 634, 682, 275,
                   760, 792, 784, 887, 112, 100, 203, 348, 616, 498, 642, 643, 703, 804, 830, 208,
                   233, 234, 246, 352, 372, 833, 428, 440, 578, 752, 826, 925, 8, 20, 70, 191, 292,
                   300, 336, 380, 470, 499, 620, 674, 688, 705, 724, 807, 40, 56, 250, 276, 438,
                   442, 492, 528, 756, 999, 832, 831, 818, 90, 36, 520, 598, 798, 296, 316, 580,
                   583, 584, 585]
        self.America = [660, 28, 533, 44, 52, 92, 535, 136, 192, 531, 212, 214, 308, 312, 332,
                   388, 474, 500, 630, 659, 662, 670, 534, 780, 796, 850, 84, 188, 222, 320, 340,
                   484, 558, 591, 32, 68, 76, 152, 170, 218, 238, 254, 328, 600, 604, 740, 858,
                   862, 60, 124, 304, 666, 840, 652, 663]
        self.Australia = [36, 554, 242, 540, 598, 90, 548, 316, 296, 584, 583, 520, 580, 585, 16, 184,
                     258, 570, 882, 772, 776, 798, 876, 574]

    def mainProcess(self):

        # Load Rasters

        gdalBuffer = gdal.Open(self.bufferFile)
        gdalPop = gdal.Open(self.popFile)
        gdalDensity = gdal.Open(self.densityfile)
        gdalNationID = gdal.Open(self.nationID)
        gdalAcc = gdal.Open(self.acc)

        # Whack them into an array
        buffer = np.array(gdalBuffer.GetRasterBand(1).ReadAsArray())
        pop = np.array(gdalPop.GetRasterBand(1).ReadAsArray())
        density = np.array(gdalDensity.GetRasterBand(1).ReadAsArray())
        nationID = np.array(gdalNationID.GetRasterBand(1).ReadAsArray())
        acc = np.array(gdalAcc.GetRasterBand(1).ReadAsArray())

        # Flatten original arrays
        buffer = buffer.flatten()
        pop = pop.flatten()
        density = density.flatten()
        nationID = nationID.flatten()
        acc = acc.flatten()

        # Stack the array
        combArray = np.column_stack((buffer, pop, density, nationID, acc))

        # Load the stacked array into a data frame
        dataframe = pd.DataFrame(combArray)
        dataframe.columns = ['buffer', 'pop', 'density', 'nationID', 'acc']

        # Interested only in region we clipped for analysis
        dataframe = dataframe[dataframe.acc > 0]

        #### get data
        # Loop through countries
        # what country?
        if self.continent == "Africa":
            countries = self.Africa
        if self.continent == "Eurasia":
            countries = self.Eurasia
        if self.continent == "America":
            countries = self.America
        if self.continent == "Australia":
            countries = self.Australia

        with open(self.resultsFile, 'w+') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=',')

            csvwriter.writerow(['Country ID', 'Urban Pop < 1km', 'Rural Pop < 1km', 'Urban Pop > 1km',
                                'Rural Pop > 1km', 'Total Population'])

            for i in countries:
                country = i
                # Urban / Rural cut-off according to European Commission. Urban clusters are those
                # with at least 300 inhabitants per kilometre squared.
                rural = 300
                dataframe1 = dataframe[dataframe.nationID==country]
                # debug to get rid of no-data cells
                dataframe1 = dataframe1[dataframe1['pop'] > 0]
                # what is the total population?
                totalPop = dataframe1['pop'].sum()
                # Urban population within 1km of river or coast
                dataframe2 = dataframe1[dataframe1.density>rural]
                uPopWithin1km = dataframe2['buffer'].sum()
                # Rural population within 1km of river or coast
                dataframe3 = dataframe1[dataframe1.density<=rural]
                rPopWithin1km = dataframe3['buffer'].sum()
                # Urban population not within 1km of river or coast
                uPop = dataframe2['pop'].sum()
                uPopOutside1km = uPop - uPopWithin1km
                # Rural population not within 1km of river or coast
                rPop = dataframe3['pop'].sum()
                rPopOutside1km = rPop - rPopWithin1km

                col1= country
                col2 = uPopWithin1km
                col3 = rPopWithin1km
                col4 = uPopOutside1km
                col5 = rPopOutside1km
                col6 = totalPop

                csvwriter.writerow([col1, col2, col3, col4, col5, col6])

##################################################################################################
# Main
#
# analysis.py -i 'inputDirectory' -n 'continent

if __name__ == "__main__":
    version_num = int(gdal.VersionInfo('VERSION_NUM'))
    if version_num < 1800:  # because of GetGeoTransform(can_return_null)
        print('ERROR: Python bindings of GDAL 1.8.0 or later required')
        sys.exit(-1)

    parser = argparse.ArgumentParser(description='Run analysis')
    apg_input = parser.add_argument_group('Input')
    apg_input.add_argument("-i", "--inputDir", nargs='?',
                           help="Path to directory containing the 5 required input files. "
                                "Directory should contain the following files: 'acc.tif', 'buffer.tif',"
                                "'density.tif', 'pop.tif', 'nationID.tif'")
    apg_input.add_argument("-o", "--overwrite", action='store_true',
                           help="Will overwrite any existing files found in output directory")
    apg_input.add_argument("-n", "--continent", nargs='?',
                           help="Name of the 'mega' continent on which analysis is taking place.")
    options = parser.parse_args()

    A = analysis(inputDir=options.inputDir, overwrite=options.overwrite, continent=options.continent)

    if not os.path.isfile(A.resultsFile) or options.overwrite:
        A.mainProcess()
    else:
        print('Results file exists', A.resultsFile)

##################################################################################################
