########################################################################
# riverBuffer.py
#
#   This code takes a river accumulation file, 'acc.tif' (from burnCoastline.py)
#   and a population file 'pop.tif' and carries out a neighborhood analysis
#   to determine which population cells are within a specified distance from
#   the river.
#
#   INPUTS:  accumulation file: 'acc.tif'; population file: 'pop.tif'
#   OUTPUTS: population within river buffer file: 'buffer.tif'
#
#   NOTE: user needs to specify the river accumulation threshold (in km^2) and
#   the chosen buffer distance from the river (in km)
#
# Author: Mark Bernhofen (2019) cn13mvb@leeds.ac.uk
########################################################################


import argparse
import os
import sys
from datetime import datetime

import numpy as np
from osgeo import gdal
import math

class riverBuffer:

    def __init__(self, inputDir, outputDir, overwrite, accThresh, bufferDistance):
        self.inputDir = inputDir
        if not os.path.exists(self.inputDir):
            print("Input directory doesn't exist")
            sys.exit(-1)
        self.outputDir = outputDir
        if not os.path.exists(self.outputDir):
            print("Output directory does not exist. Creating directory...")
            os.makedirs(self.outputDir)
        self.overwrite = overwrite
        self.accThresh = accThresh
        self.bufferDistance = bufferDistance

        # Define input filenames and then read in the data
        self.accFile = os.path.join(self.inputDir, 'acc.tif')
        self.acc = self.readData(self.accFile, saveInfo=True)
        self.popFile = os.path.join(self.inputDir, 'pop.tif')
        self.pop = self.readData(self.popFile).astype(np.float32)
        self.bufferFile = os.path.join(self.outputDir, "buffer.tif")

        # Get raster info
        self.ulx, self.xres, self.xskew, self.uly, self.yskew, self.yres = self.dataInfo.GetGeoTransform()

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

        # Extract river from the accumulation file
        self.river = np.where(self.acc >= int(self.accThresh), 1, 0)
        self.pop = np.where(self.pop>0, self.pop, 0)

        # Create buffered population array
        self.buffer = np.full(self.acc.shape, 0, dtype=np.float32)

        # Loop through all the river cells
        for r in range(self.rasterYSize):
            for c in range(self.rasterXSize):
                if self.river[r, c] == 1:
                    # If river cell overlaps with population cell, copy population onto buffer
                    if self.pop[r, c] > 0:
                        pop = self.pop[r, c]
                        self.buffer[r, c] = pop
                    for m in range(int(self.bufferDistance)+1):
                        for n in range(int(self.bufferDistance)+1):
                            neighbors = self.neighborhood(r, c, m, n)
                            for idx, n in enumerate(neighbors):
                                pixel = (n[0], n[1], idx)
                                if pixel[0]<0:
                                    continue
                                if pixel[1]<0:
                                    continue
                                if pixel[0] > self.rasterYSize - 1:
                                    continue
                                if pixel[1] > self.rasterXSize - 1:
                                    continue
                                self.neighborProcess(pixel, r, c)

    def neighborhood(self, row, col, m, n):
        '''
        This function returns a vector containing the coordinates of the current pixel's 8
        neighbors. The neighboring elements are labelled as follows:
        5  6  7
        4  x  0
        3  2  1
        '''
        neighbors = [[row + m, col + n], [row - m, col - n], [row + m, col - n], [row - m, col + n]]

        return neighbors

    def neighborProcess(self, neighbor, riverRow, riverCol):
        row = neighbor[0]
        col = neighbor[1]
        if row < 0:
            return
        if col < 0:
            return

        cRow = riverRow
        cCol = riverCol

        distance = self.haversineDistance(cRow, cCol, row, col)
        if distance < int(self.bufferDistance) and distance > 0:
            if self.buffer[row, col] > 0:
                return
            if self.pop[row, col] > 0:
                pop = self.pop[row, col]
                self.buffer[row, col] = pop
                return
        else:
            return

    def haversineDistance(self, cRow, cCol, row, col):
        '''
        This function uses the haversine formula to calculate the 2D distance between the current cell
        and the river cell.
        '''

        coastCoord = self.coordinateConverter(cRow, cCol)
        neighborCoord = self.coordinateConverter(row, col)

        R = 6372.8 # Earth's radius in kilometers
        lat1, lon1 = coastCoord
        lat2, lon2 = neighborCoord

        phi1, phi2  = math.radians(lat1), math.radians(lat2)
        dphi        = math.radians(lat2 - lat1)
        dlamda      = math.radians(lon2 - lon1)

        a = math.sin(dphi/2)**2 + \
            math.cos(phi1)*math.cos(phi2)*math.sin(dlamda/2)**2

        distance = 2*R*math.atan2(math.sqrt(a), math.sqrt(1 - a))

        return distance

    def coordinateConverter(self, r, c):
        '''
        This function converts cell array coordinates into geographic coordinates
        '''

        row = r
        col = c

        x = self.ulx + (col*self.xres)
        y = self.uly - (row*self.yres)
        coordinate = (x, y)

        return coordinate

    def writeTifOutput(self):
        print(str(datetime.now()), 'Saving GeoTiffs...')

        # Write text file with layer information
        fLoc = os.path.join(self.outputDir, "read_me.txt")
        f = open(fLoc, "w+")
        f.write("River accumulation threshold " + str(self.accThresh))
        f.write("Buffer distance (kms)  " + str(self.bufferDistance))
        f.close()

        driver = gdal.GetDriverByName("GTiff")
        dst = driver.Create(self.bufferFile, self.rasterXSize, self.rasterYSize, 1, gdal.GDT_Float32)
        band = dst.GetRasterBand(1)
        band.WriteArray(self.buffer, 0, 0)

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
# riverBuffer.py -i 'inputDirectory' -d 'outputDirectory' -a 'accumulationThreshold' -b 'bufferDistance'

if __name__ == "__main__":
    version_num = int(gdal.VersionInfo('VERSION_NUM'))
    if version_num < 1800:  # because of GetGeoTransform(can_return_null)
        print('ERROR: Python bindings of GDAL 1.8.0 or later required')
        sys.exit(-1)

    parser = argparse.ArgumentParser(description='Generate river/coast buffer')
    apg_input = parser.add_argument_group('Input')
    apg_input.add_argument("-i", "--inputDir", nargs='?',
                           help="Path to directory containing the 2 required input files. "
                                "Directory should contain accumulation and population files with names"
                                "of 'acc' and 'pop' respectively.")
    apg_input.add_argument("-d", "--outputDir", nargs='?',
                           help="Path to desired output directory. If no directory exists, one will be created.")
    apg_input.add_argument("-o", "--overwrite", action='store_true',
                           help="Will overwrite any existing files found in output directory")
    apg_input.add_argument("-a", "--threshold", default=25, nargs='?',
                           help="River upstream drainage area threshold to consider [kilometres squared]")
    apg_input.add_argument("-b", "--bufferDistance", default=1, nargs='?',
                           help="Buffer distance to calculate [kilometres]")
    options = parser.parse_args()

    B = riverBuffer(inputDir=options.inputDir, outputDir=options.outputDir, overwrite=options.overwrite,
             accThresh=options.threshold, bufferDistance=options.bufferDistance)

    if not os.path.isfile(B.bufferFile) or options.overwrite:
        B.mainProcess()
        B.writeTifOutput()
    else:
        print('Buffer file exists', B.bufferFile)

##################################################################################################






