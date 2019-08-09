######################################################################
# burnCoastline.py
#
#   This code takes an upstream accumulation file (name='acc.tif') that
#   has already been assigned a no-data value (in QGIS) and burns a user
#   specified 'coastline' value onto the the accumulation file.
#
#   NOTE: make sure the coastline burn value you choose is greater than
#   the minimum river accumulating area threshold that you plan on using
#
#   INPUTS: accumulation file 'acc.tif'
#   OUTPUTS: accumulation file (with coastline burn) 'acc.tif'
#
# Author: Mark Bernhofen (2019) cn13mvb@leeds.ac.uk
######################################################################


import argparse
import os
import sys
from datetime import datetime

from osgeo import gdal

class burnCoastline:

    def __init__(self, inputDir, outputDir, overwrite, noData, coastlineVal):
        self.inputDir = inputDir
        if not os.path.exists(self.inputDir):
            print("Input directory doesn't exist")
            sys.exit(-1)
        self.outputDir = outputDir
        if not os.path.exists(self.outputDir):
            print("Output directory does not exist. Creating directory...")
            os.makedirs(self.outputDir)
        self.overwrite = overwrite
        self.noData = noData
        self.coastlineVal = coastlineVal

        # Define input filenames and then read the data
        self.accFile = os.path.join(self.inputDir, 'acc.tif')
        self.acc = self.readData(self.accFile, saveInfo=True)
        self.newAccFile = os.path.join(self.outputDir, 'acc.tif')

    def readData(self, fileName, saveInfo=False):
        fileHandle = gdal.Open(fileName)
        if fileHandle is None:
            print('Error: data file no data: ', fileName)
            sys.exit(-1)

        self.rasterXSize = fileHandle.RasterXSize
        self.rasterYSize = fileHandle.RasterYSize

        if saveInfo:
            self.dataInfo = fileHandle

        band = fileHandle.GetRasterBand(1)
        data = band.ReadAsArray(0, 0, fileHandle.RasterXSize, fileHandle.RasterYSize)

        fileHandle = None
        return data

    def mainProcess(self):

        # Create new accumulation array
        self.accNew = self.acc

        # Loop through the accumulation file
        for r in range(self.rasterYSize):
            for c in range(self.rasterXSize):
                if self.acc[r,c] > 0:
                    neighbors = self.neighborhood(r, c)
                    for idx, n in enumerate(neighbors):
                        if n is not None:
                            neighbor = (n[0], n[1], idx)
                            self.neighborProcess(neighbor, r, c)

    def neighborhood(self, row, col):
        '''
        This function returns a vector containing the coordinates of the current pixel's 8
        neighbors. The neighboring elements are labelled as follows:
        5  6  7
        4  x  0
        3  2  1
        '''
        neighbors = [[row, col + 1], [row + 1, col + 1], [row + 1, col], [row + 1, col - 1], [row, col - 1],
                     [row - 1, col - 1], [row - 1, col], [row - 1, col + 1]]
        for i in range(8):
            neighborR, neighborC = neighbors[i]
            if (neighborR < 0) or (neighborR > self.rasterYSize - 1):
                neighbors[i] = None
            elif (neighborC < 0) or (neighborC > self.rasterXSize - 1):
                neighbors[i] = None

        return neighbors

    def neighborProcess(self, neighbor, r, c):
        '''
        This function checks if the neighboring cell is a no-data cell (ocean). If it is,
        then the current cell is taken as a coastal cell.
        '''
        nRow = neighbor[0]
        nCol = neighbor[1]
        cRow = r
        cCol = c

        if self.acc[nRow, nCol] == self.noData:
            self.accNew[cRow, cCol] = self.coastlineVal

        return

    def writeTifOutput(self):
        print(str(datetime.now()), 'Saving GeoTiffs...')

        # Write text file with layer information
        fLoc = os.path.join(self.outputDir, "read_me.txt")
        f = open(fLoc, "w+")
        f.write("River coastline burn value " + str(self.coastlineVal))
        f.close()

        driver = gdal.GetDriverByName("GTiff")
        dst = driver.Create(self.newAccFile, self.rasterXSize, self.rasterYSize, 1, gdal.GDT_Int32)
        band = dst.GetRasterBand(1)
        band.WriteArray(self.accNew, 0, 0)

        projection = self.dataInfo.GetProjection()
        geotransform = self.dataInfo.GetGeoTransform()

        dst.SetGeoTransform(geotransform)
        dst.SetProjection(projection)
        dst.FlushCache()

        dst = None
        self.dataInfo = None

        return

##################################################################################################
# Main
#
# burnCoastline.py -i 'inputDirectory' -d 'outputDirectory' -n 'noData' -c 'coastlineBurn'

if __name__ == "__main__":
    version_num = int(gdal.VersionInfo('VERSION_NUM'))
    if version_num < 1800:  # because of GetGeoTransform(can_return_null)
        print('ERROR: Python bindings of GDAL 1.8.0 or later required')
        sys.exit(-1)

    parser = argparse.ArgumentParser(description='Burn coastline into accumulation file')
    apg_input = parser.add_argument_group('Input')
    apg_input.add_argument("-i", "--inputDir", nargs='?',
                           help="Path to directory containing the 2 required input files. "
                                "Directory should contain an accumulation file with the name 'acc'")
    apg_input.add_argument("-d", "--outputDir", nargs='?',
                           help="Path to desired output directory. If no directory exists, one will be created.")
    apg_input.add_argument("-o", "--overwrite", action='store_true',
                           help="Will overwrite any existing files found in output directory")
    apg_input.add_argument("-n", "--noData", default=-9999, nargs='?', type=int,
                           help="What is the noData value of the input accumulation dataset")
    apg_input.add_argument("-c", "--coastlineBurn", default=9999, nargs='?', type=float,
                           help="What coastline value do you want to burn into the accumulation file ")
    options = parser.parse_args()

    B = burnCoastline(inputDir=options.inputDir, outputDir=options.outputDir, overwrite=options.overwrite,
             noData=options.noData, coastlineVal=options.coastlineBurn)

    if not os.path.isfile(B.newAccFile) or options.overwrite:
        B.mainProcess()
        B.writeTifOutput()
    else:
        print('Coastline burned accumulation file exists', B.newAccFile)

##################################################################################################