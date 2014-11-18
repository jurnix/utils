#!/bin/sh
"""
Thai met

toNetcdf_1.py converts Thai met data from Excel files to netCDF format.
It tries to follow the convention specified under Climate Forecast. 
Check the link below to get more information:

http://cf-pcmdi.llnl.gov/

Met stations are sparsed in different locations in Thailand. So they do not 
form a grid, the most common shape for NetCDF files. For this reason,
code stations have become an index of each station. In consequence,
'time' and 'code stations' are a dimension. Otherwise, using coordinates 
would suppose an inefficent way of handle the data (empty spaces, ...).

Note: this script doesn't work with scipy 0.9 version. 
There is a bug when writing data.

https://github.com/scipy/scipy/blob/master/scipy/io/netcdf.py
Types of data on netCDF
TYPEMAP = { NC_BYTE:    ('b', 1),
            NC_CHAR:    ('c', 1),
            NC_SHORT:   ('h', 2),
            NC_INT:     ('i', 4),
            NC_FLOAT:   ('f', 4),
            NC_DOUBLE:  ('d', 8)}
            
Version 0.2 - Moved netCDF to version 4

@author: Albert Jornet Puig
"""

import xlrd
import logging
import unittest
import numpy as np
from datetime import date
#from scipy.io import netcdf
import netCDF4 as netcdf

logging.basicConfig()
log = logging.getLogger('toNetcdf')
log.setLevel(logging.INFO)

RELATED_DIMS = ('time','stncode', )

EMPTY = -999

STARTDATE = [1951, 1, 1]
ENDATE = [2012,12,31 ]
TIMEUNITS = 'Days since 1951-01-01 00:00:00'

NVARS = ["maxtmp", "mintmp", "rain", "avgrh"]
NDIMS = ["stn_name", "stncode", "latitude", "longitude", "time"]

XLS_CELL_TYPE = dict({0: 'Empty', 
                   1: 'Text', 
                   2:'Number', 
                   3: 'Date', 
                   4:'Boolean', 
                   5:'Error', 
                   6:'Blank'})

# save dates into int 
reverse_date = {}    
    
class TestNetCDFValues(unittest.TestCase):
    """ Test data has been well written in NetCDF data file. """
    filename = 'data.nc'
    cases = [
    [455201, 1985, 5, 29, 32.2, 25.9, 41.7, 84],
    [455601, 1992, 12, 25, 21.4, 16.9, 4.9, 88]
    ]
    def setUp(self):
        self.nc = netcdf.Dataset(self.filename, 'r')

    def test_first(self):
        ind = 0
        
        values = self.cases[ind]
        codes = self.nc.variables['stncode']
        rain = self.nc.variables['rain']
        avg = self.nc.variables['avghr']
        mintmp = self.nc.variables['mintmp']
        maxtmp = self.nc.variables['maxtmp']
        
        nday = (date(values[1], values[2], values[3]) - \
                date(STARTDATE[0], STARTDATE[1], STARTDATE[2])).days
        idx_code = list(codes[:]).index(values[0])
        log.debug("day= %s" % nday)
        log.debug(rain[nday, :])

        # test data        
        self.assertEqual(values[4], maxtmp[nday, idx_code])
        self.assertEqual(values[5], mintmp[nday, idx_code])
        self.assertEqual(values[6], rain[nday, idx_code])
        self.assertEqual(values[7], avg[nday, idx_code])
    
    def test_second(self):
        ind = 1
        
        values = self.cases[ind]
        codes = self.nc.variables['stncode']
        avghr = self.nc.variables['avghr']
        rain = self.nc.variables['rain']
        mintmp = self.nc.variables['mintmp']
        maxtmp = self.nc.variables['maxtmp']
        
        nday = (date(values[1], values[2], values[3]) - \
                date(STARTDATE[0], STARTDATE[1], STARTDATE[2])).days
        idx_code = list(codes[:]).index(values[0])
        
        # test data        
        self.assertEqual(values[4], maxtmp[nday, idx_code])
        self.assertEqual(values[5], mintmp[nday, idx_code])
        self.assertEqual(values[6], rain[nday, idx_code])
        self.assertEqual(values[7], avghr[nday, idx_code])


def clean_string(s):
    ''' Remove non desired chars '''
    tmp = s.replace('\t','')
    tmp = tmp.replace('*','')
    tmp = tmp.replace(' ','_')
    tmp = tmp.replace('\n','')
    return tmp
    

class MetData(object):
    ''' Load Thai met data into memory '''
    metstation = []
    
    def __init__(self, filename):
        print "MetData::init: Load file %s" % filename
        self.filename = filename
        # load data into memory
        f = open(filename,'r')
        for line in f:            
            d = {}
            code, name, lat, lon =  line.split('\t')[:4]
            d['stncode'] = code
            d['stnname'] = name
            d['latitude'] = lat
            d['longitude'] = lon
            self.metstation.append(d)
            
        f.close()
        print "MetData::init: OK"

    def searchByCode(self, code):
        #print "code to search: ", code, type(code)
        return self.__get_index('stncode', code)
    
    def __get_index(self, attr, value):        
        ''' Search in list of dictionary '''
        try:
            return next(d for (index, d) in enumerate(self.metstation) \
                    if d[attr].endswith(value))
        except StopIteration:
            return None
            
    def getSubset(self, codes):
        ''' It returns a new list with all data related with each code given '''
        output = []
        for code in codes:
            entry = self.searchByCode(code)
            output.append(entry.copy())
            
        return output

class xlsSheetParser:
    """ Iterate each row found on the specified spreadSheet """
    def __init__(self, filename, nsheet):
        """
        Keyword arguments:
            filename -- xls file to read
            nsheet -- number of worksheet to read from file
        """
        self.filename = filename        
        self.nsheet = nsheet
        # open
        book = xlrd.open_workbook(self.filename)
        # point to xls worksheet
        self.curr_sheet = book.sheet_by_index(self.nsheet)
        #indices
        self.total_rows = self.curr_sheet.nrows - 1
        self.total_cols = self.curr_sheet.ncols - 1
        self.current_row = 0
        
    def __parse_row(self):
        current_col = -1
        output = []
        
        while current_col < self.total_cols:
            current_col += 1
            # Cell Types: 0=Empty, 1=Text, 2=Number, 
            #3=Date, 4=Boolean, 5=Error, 6=Blank
            #cell_type = self.curr_sheet.cell_type(self.current_row, current_col)
            cell_value = self.curr_sheet.cell_value(self.current_row, current_col)
            #print '	', XLS_CELL_TYPE[cell_type], ':', cell_value
            if isinstance(cell_value, basestring):
                cell_value = clean_string(cell_value)
                
            if not cell_value:
                cell_value = EMPTY
                
            output.append(cell_value)
            
        return output

    def __iter__(self):
        return self
    
    def next(self):
        if self.current_row >= self.total_rows:
            raise StopIteration
        else:
            self.current_row += 1            
            return self.__parse_row()
            
    def getColumns(self, keys = [], init = -1):
        ''' 
        Returns data from specified cols
        returns: dict
        '''
        output = []
        
        header_pos = []
        self.current_row = init
        
        # save index position of list
        for header in self:
            for i, key in enumerate(keys):
                if not key in header:
                    raise Exception("key %s not found on %s" % (key, header))
                header_pos.append(i)
            break
        
        for row in self:
            d = []
            for pos in header_pos:
                #d[keys[pos]] = row[pos]
                #print row[pos]
                d.append(row[pos])

            output.append(d)
                
        return output
        
    def getSetColumns(self, keys = [], init = -1):
        ''' 
        Returns data from specified cols with non repeated elements
        returns: dict
        '''
        return frozenset( self.getColumns(keys, init) )


def discover_files():
    """ Returns a list with all xls with its relative path """
    import fnmatch
    import os

    matches = []
    for root, dirnames, filenames in os.walk('.'):
        for filename in fnmatch.filter(filenames, '*.xls'):
            matches.append(os.path.join(root, filename))
    return matches

def get_coordinates_from_id(database, code):
    ''' from MetData clase search data from the given code '''          
    str_c = str(code).split(".")[0]
    latAndLong = database.searchByCode(str_c)
    if latAndLong:
        return float(latAndLong['latitude']), float(latAndLong['longitude'])
    else:
        #if not code == '360201' and not code == '429301':
        #    print latAndLong
        #    raise Exception("Station's coordinates with id %s not found" % code)
        #else:
        print "ignored code", code, "NOT FOUND"
        
    return EMPTY, EMPTY
            
def get_zone(string):
    """ Extracts from input the zone. E.g. ./XX/filename """
    return string.split('/')[1]

def createNcStrlen(nc, length):
    ''' '''
    nc.createDimension('name_strlen', length) 

def createNcZoneStrlen(nc):
    ''' '''
    nc.createDimension('zone_strlen', 2)

def createNcDate(nc, ndays):
    ''' Set date in netcdf as dimension '''   
    nc.createDimension('time', ndays) 
    dates = nc.createVariable('time', 'i', ('time', ))
    dates.units = TIMEUNITS
    dates.long_name = 'time'
    dates[:] = np.arange(ndays)
    
def createNcCodeStations(nc, nstations, l_code_stations):
    ''' Set code stations as dimension '''
    nc.createDimension('stncode', nstations) # station identifier
    codes = nc.createVariable('stncode', 'i', ('stncode', ))
    codes.units = 'str'
    codes.standard_name = 'station_codes'
    codes.long_name = 'station codes'
    codes[:] = l_code_stations

def createNcStations(nc, l_names):
    ''' Stations Name '''
    names = nc.createVariable('stnname', 'c', ('stncode', 'name_strlen'))
    names.units = 'str'
    names.standard_name = 'stations_names'
    names.long_name = 'met stations name'

    for i, name in enumerate(l_names):
        for j, letter in enumerate(name):
            # copy letter by letter
            names[i,j] = letter     
        # fill remaining chars
        #print len(names[i]), len(name)
        #print "total", (len(names[i])-len(name))
        #names[i, len(name):] = '' * (len(names[i])-len(name))
    
def createNcZones(nc, l_zones):
    ''' set zones '''
    zones = nc.createVariable('zone', 'c', ('stncode', 'zone_strlen' ))
    zones.units = 'str'
    zones.standard_name = 'region'
    zones.long_name = 'regions areas of Thailand'
    for i, name in enumerate(l_zones):
        for j, letter in enumerate(name):
            # copy letter by letter
            zones[i,j] = letter
    #zones[:] = l_zones
    
def createNcLat(nc, l_lats):
    ''' set latitude '''
    lats = nc.createVariable('latitude', 'd', ('stncode', ))
    lats.units = 'degrees_north'
    lats.long_name = 'latitude'
    lats.standard_name = 'latitude'
    lats[:] = l_lats
    
def createNcLong(nc, l_longs):
    ''' set longitude '''
    longs = nc.createVariable('longitude', 'd', ('stncode', ))
    longs.units = 'degrees_east'
    longs.long_name = 'longitude'
    longs.standard_name = 'longitude'
    longs[:] = l_longs
    
def createNcMaxTemp(nc):
    ''' set for netcdf the max temperature '''
    maxtmps = nc.createVariable('maxtmp', 'd', RELATED_DIMS)
    maxtmps.units = 'degrees'
    maxtmps.standard_name = "air_temperature"
    maxtmps.long_name = "maximum temperature in degress celsius"
    maxtmps.missval = '-999'
    return maxtmps
    
def createNcMinTemp(nc):
    ''' set for netcdf the min temp '''
    mintmps = nc.createVariable('mintmp', 'd', RELATED_DIMS)
    mintmps.units = 'degrees'    
    mintmps.standard_name = "air_temperature"
    mintmps.long_name = "minimum temperature in degress celsius"
    mintmps.missval = '-999'
    return mintmps

def createNcRain(nc):
    ''' set for netcdf rain measurements '''
    rains = nc.createVariable('rain', 'd', RELATED_DIMS)    
    rains.units = 'mm/day' # not as CF convention
    rains.standard_name  = 'precipitation_amount'    
    rains.long_name  = 'relative humidity'
    rains.missval = '-999'
    return rains

def createNcAvgHum(nc):
    ''' set for netcdf average relative humidity '''
    avghrs = nc.createVariable('avghr', 'i', RELATED_DIMS)
    avghrs.units = '%' # not as CF convention
    avghrs.standard_name  = 'relative_humidity'
    avghrs.long_name  = 'average relative humidity'
    avghrs.missval = '-999'
    return avghrs
    
def main():
    """ """
    global reverse_date
    xlsSheet = 0
    output = "data.nc"
    metdata = MetData("met_stations_coords.txt")
    #NZONES = 6# cc, ee, ne, nn, se, sw
    
    files = sorted(discover_files())
    log.info("files=%s" % files)
    nstations = 0    
    found_zones = set([])
    all_data = set([])

    ## -------------------------------------------------------
    ## COLLECT ALL DATA FROM XLS, COORDINATES AND ZONES
    ## -------------------------------------------------------

    # collect all name and code stations and region codes
    for filename in files:
        tmp_zone = get_zone(filename)
        
        found_zones.add(tmp_zone)
        xls = xlsSheetParser(filename, xlsSheet)
        cols = xls.getColumns(['stn_name','stncode'])
        
        # add zone
        [x.append(tmp_zone) for x in cols]
        # from list to set -identify unique names-
        tmp_data = [tuple(x) for x in cols]
        tmp_data = frozenset(tmp_data)        
        # join results
        all_data = all_data | tmp_data

    # convert to list    
    all_data = list(all_data)
    all_data = [list(x) for x in all_data]
    nstations = len(all_data)    
    
    log.debug("nstations= %d" % nstations)
    log.debug("all_data=%s" % all_data)
    # collect all latitudes and longitudes
    for station in all_data:
        str_c = str(int(station[1]))
        latt, longtt = get_coordinates_from_id(metdata, str_c)
        station.append(latt)
        station.append(longtt)
    
    # number of days
    d0 = date(ENDATE[0], ENDATE[1], ENDATE[2])
    d1 = date(STARTDATE[0], STARTDATE[1], STARTDATE[2])
    delta = d0 - d1
    totalDays = delta.days + 1
    
    ## ---------------------------------------
    ## GENERATE NETCDF 
    ##
    ## CREATE DIMENSIONS AND VARIABLES
    ## ---------------------------------------
    nc = netcdf.Dataset(output, 'w')
    nc.history = 'Release 0.2'
    nc.conventions = 'CF-1.0'
    nc.title = '''
    Thailand data of temp, rain and average relative humidity from 1951 to 2012
    obtained from meteorological stations
    '''    
    
    l_names  = np.array([x[0] for x in all_data])
    l_codes  = np.array([int(x[1]) for x in all_data], dtype=np.int)
    l_zones  = np.array([x[2] for x in all_data])
    l_lats   = np.array([x[3] for x in all_data])
    l_longs  = np.array([x[4] for x in all_data])
    
    # set dimensions and variables
    createNcDate(nc, totalDays)
    createNcCodeStations(nc, nstations, l_codes)
    # metadata
    createNcStrlen(nc, 40)
    createNcZoneStrlen(nc)
    createNcStations(nc, l_names)
    createNcZones(nc, l_zones)
    createNcLat(nc, l_lats)
    createNcLong(nc, l_longs)
    # values to save
    maxtmps = createNcMaxTemp(nc)
    mintmps = createNcMinTemp(nc)
    rains = createNcRain(nc)
    avghrs = createNcAvgHum(nc)
    
    # DEBUG    
    log.debug( "Names = %s" % l_names)
    log.debug("Codes = %s" % l_codes)
    log.debug("Zones = %s" % l_zones)
    log.debug("Lats = %s" % l_lats)
    log.debug("Longs = %s" % l_longs)
    
    #nc.flush()
    
    ## ---------------------------------------
    ## FILL NETCDF WITH DATA
    ## ---------------------------------------
    rain_data = np.zeros(shape=(totalDays, len(l_codes)), dtype = np.float)
    avg_data = np.zeros(shape=(totalDays, len(l_codes)), dtype = np.int)
    
    # save data to netCDF        
    for i, filename in enumerate(files):  
        log.info("%d/%d filename %s"% (i, len(files), filename))
        itexls = xlsSheetParser(filename, 0)        
        # ignore header
        for row in itexls:
            break
        
        #matr = np.array(ndmin=2)
        #j = -1
        # get rows data
        for row in itexls:
            #print row
            # splitted into real meaning
            rname, rcode = row[0], row[1]
            r_max, rmin, r_rain, r_avg = row[5], row[6], row[7], row[8]
            try:
                r_date = date(int(row[2]), int(row[3]), int(row[4]))
                
            except ValueError:
                log.error( "Not a valid date %d/%d/%d with data %s" % \
                        (int(row[2]), int(row[3]), int(row[4]), str(row[5:])))
                continue
            
            # format
            str_code = str(int(rcode))
            latt, longtt = get_coordinates_from_id(metdata, str_code)
            nday = (r_date - date(STARTDATE[0], STARTDATE[1], STARTDATE[2])).days
            reverse_date[nday] = r_date

            code_ind = list(l_codes).index(rcode)
            log.debug("%s to %d position=%d (%d) value=%f" % 
                        (r_date, nday, code_ind, rcode, r_rain))
            
            # write data           
            mintmps[nday, code_ind] = rmin
            maxtmps[nday, code_ind] = r_max
            rain_data[nday, code_ind] = r_rain
            avg_data[nday, code_ind] = r_avg
            
        log.info("-----------------------------------------------------------------")
       
    for i, line in enumerate(rain_data):
        if not i in reverse_date:
            log.debug("value not found ", i, line)
            continue
        log.debug( "%s(%d) %s" % (reverse_date[i], i, line))
            
    rains[:,:] = rain_data
    avghrs[:,:] = avg_data
    #nc.flush()
    
    nc.close()
    log.info("Done")
    
if __name__ == "__main__":
    main()
    unittest.main()
