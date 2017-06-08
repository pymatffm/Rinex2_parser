#!/usr/bin/env python
#*-* coding: utf-8 *-*
from __builtin__ import file
import os.path
import sys
import datetime
import pandas as pd
import numpy as np
import scipy as sp
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from math import isnan 
import platform
import pylab
import rinex3_parser
import constants as cs
import vtec_converter
import rinex3_p_parser
import rinex2_nav
import dcb_calculation
import gps_sow_time
import tec_calculations
from matplotlib.ticker import MaxNLocator
from collections import Counter
from mpl_toolkits.basemap import Basemap
from obspy.imaging.mopad_wrapper import beach

#Global functions
def floatornan(x):
    if x == '' or x[-1] == ' ':
        return np.nan #np.NaN
    else:
        return float(x)

def digitorzero(x):
    if x == ' ' or x == '':
        return 0
    else:
        return int(x)

def padline(l, n=16):
    x = len(l)
    x_ = n * ((x + n - 1) / n)
    test = l + ' ' * (x_ - x)
    return l + ' ' * (x_ - x)


class RINEX_2_Formatting:
        def __init__(self, filename, elevation_file, igs_station, azimuth_file, selected_nav_rinex, selected_dcb_p1c1, selected_dcb_p1p2, selected_dcb_p2c2):  

            workspacePath = os.getcwd()
            self.pathList = workspacePath.split(os.sep)
            
            # DCB calculations: Preliminary results of DCB C1-C2 in a GPS system, KRASUSKI
            # Slide 3 of Biases in GNSS Analysis, Schaer
            # C1-C2: (P1-P2) + (P1-C1) - (P2-C2)
            d = dcb_calculation.DCB_generation(selected_dcb_p1c1)
            dcb_p1c1 = d.dcb_dictionary
            
            e = dcb_calculation.DCB_generation(selected_dcb_p1p2)
            dcb_p1p2 = e.dcb_dictionary
            
            f = dcb_calculation.DCB_generation(selected_dcb_p2c2)
            dcb_p2c2 = f.dcb_dictionary
            
            self.final_dict = self.merge_dcb_values(dcb_p1c1, dcb_p1p2, dcb_p2c2)
                        
            with open(filename, 'r+') as RINEXfile:
                version_line = padline(RINEXfile.readline(), 80)[0:9]
                version_line = float(version_line)
                path, file = os.path.split(selected_nav_rinex)
                nav_type = file[-1:]
                if nav_type == 'g':
                    gnss_sys = 'R'
                elif nav_type == 'n':
                    gnss_sys = 'G'
                
                if version_line > 2.12:
                    rinex3_p_parser.RINEX3_p_Files(selected_nav_rinex, filename)

                else:
                    #gnss_sys = 'E'
                    a = vtec_converter.VTEC_Calculating(elevation_file, igs_station, azimuth_file, gnss_sys)
                    elevation_array = a.elevation_array
                    azimuth_array = a.revised_azimuth_array
                    self.mf_dataframe = a.inv_mf_dataframe
                    self.read_header(RINEXfile)
                    self.read_data(RINEXfile, selected_nav_rinex)
                    tec_calculations.Calculate_VTEC(RINEXfile, elevation_array)

					
        def merge_dcb_values(self, dcb_p1c1, dcb_p1p2, dcb_p2c2):
            # -B_(P1-P2) + B_(P1-C1) - B_(P2-C2)
            # Stefan Schaer 24.10.2016         
            intermediate_dict = {}
            intermediate_dict.update(dcb_p1p2)
            intermediate_dict.update(dcb_p1c1)
            
            for i in dcb_p1c1.keys():
                try:
                    addition = -dcb_p1p2[i] + dcb_p1c1[i]
                    intermediate_dict[i] = addition
                except KeyError:
                    continue
                
                final_dict = {}
                final_dict.update(dcb_p2c2)
                final_dict.update(intermediate_dict)
                
                for i in dcb_p2c2.keys():
                    try:
                        subtraction = intermediate_dict[i] - dcb_p2c2[i]
                        final_dict[i] = subtraction
                    except KeyError:
                        continue

            return final_dict
					
					
        def read_header(self, RINEXfile):
            version_line = padline(RINEXfile.readline(), 80)

            self.comment = ""
            self.observ = ""
            count = 0
            while True:
                line = padline(RINEXfile.readline(), 80)
                label = line[60:80].rstrip()
                if label == cs.END_OF_HEADER_2:
                    break
                if label == cs.APPROX_POSITION_XYZ or label == cs.APPROX_POSITION_XYZ_2 or label == cs.APPROX_POSITION_XYZ_3:
                    self.X_rx = float(line[:14])
                    self.Y_rx = float(line[15:28])  
                    self.Z_rx = float(line[29:42]) 
                if label == cs.COMMENT:
                    self.comment += line[:60] + '\n'
                if label == cs.MARKER_NAME_2:
                    self.marker_name = line[:60].rstrip()
                    temp = self.marker_name[:4]
                    if temp == '1427':
                        self.marker_name = 'UNSW'
                    if self.marker_name == '':
                        self.marker_name = cs.UNKNOWN
                if label == cs.TYPES_OF_OBSERV:
                    if self.observ == "":
                        self.observ = line[:60]
                        n_obs = int(self.observ[0:6])
                        self.obs_types = []
                        if n_obs < 10:
                            for i in range(0, n_obs):
                                self.obs_types.append(self.observ[10+6*i:12+6*i])
                        else:
                            for i in range(0, 9):
                                self.obs_types.append(self.observ[10+6*i:12+6*i])
                            count += 1
                    elif self.observ != "" and count < 2:
                        self.observ = line #[6:60]
                        if n_obs < 19:
                            for i in range(0, n_obs - 9):
                                self.obs_types.append(self.observ[10+6*i:12+6*i])
                        else:
                            for i in range(0, 9):
                                self.obs_types.append(self.observ[10+6*i:12+6*i])
                            count += 1
                    else:
                        self.observ = line
                        for i in range(0, n_obs - 18):
                            self.obs_types.append(self.observ[10+6*i:12+6*i])


        def read_epoch_header(self, RINEXfile):
            epoch_hdr = RINEXfile.readline()
            if epoch_hdr == '':
                return None

            year = int(epoch_hdr[1:3])
            if year >= 80:
                year += 1900
            else:
                year += 2000

            month = int(epoch_hdr[4:6])
            day = int(epoch_hdr[7:9])
            hour = int(epoch_hdr[10:12])
            minute = int(epoch_hdr[13:15])
            second = int(epoch_hdr[15:18])
            microsecond = int(epoch_hdr[19:25]) # Discard the least significant digits (use microseconds only). 
            date_USA = datetime.date(year, month, day)
            temp_date = str(date_USA)
            list_date = temp_date.split('-')
            seq = (list_date[2], list_date[1], list_date[0])
            self.date = '-'.join(seq)
            epoch_time = int((hour * 3600) + (minute * 60) + second)
            print "{0}:{1}:{2}".format(hour, minute, second)

            n_sats = int(epoch_hdr[29:32])
            sats = []
            for i in range(0, n_sats):
                if ((i % 12) == 0) and (i > 0):
                    epoch_hdr = RINEXfile.readline()
                sats.append(epoch_hdr[(32+(i%12)*3):(35+(i%12)*3)])
                
            hour_fraction = float(hour+(minute/60.0))
            julian_day = self.julian_day_conversion(year, month, day, hour_fraction)

            d = gps_sow_time.GPS_Time(year, month, day, hour_fraction)
            sow = d.sow

            sow = sow + second

            return epoch_time, sats, sow


        def read_obs(self, RINEXfile, n_sat, sat_map):
            obs = np.empty((cs.TOTAL_SATS, len(self.obs_types)), dtype=np.float64) * np.NaN
            lli = np.zeros((cs.TOTAL_SATS, len(self.obs_types)), dtype=np.uint8)
            signal_strength = np.zeros((cs.TOTAL_SATS, len(self.obs_types)), dtype=np.uint8)

            for i in range(n_sat):      
                obs_line = ''.join(RINEXfile.readline().ljust(80) for i in range((len(self.obs_types) + 4) / 5))
                for j in range(len(self.obs_types)):
                    obs_block = obs_line[(16*j):(16*(j+1))]
                    obs[sat_map[i], j] = floatornan(obs_block[1:14])
        
            return obs, lli, signal_strength


        def read_data_chunk(self, RINEXfile, selected_nav_rinex):
            obss = np.empty((cs.CHUNK_SIZE, cs.TOTAL_SATS, len(self.obs_types)), dtype=np.float64) * np.NaN
            llis = np.zeros((cs.CHUNK_SIZE, cs.TOTAL_SATS, len(self.obs_types)), dtype=np.uint8)
            signal_strengths = np.zeros((cs.CHUNK_SIZE, cs.TOTAL_SATS, len(self.obs_types)), dtype=np.uint8)
            epochs = np.zeros(cs.CHUNK_SIZE, dtype='datetime64[us]')
            epochs_times = []
            sow_list = np.zeros(cs.CHUNK_SIZE, dtype=np.float64) * np.NaN

            i = 0
            while True:
                hdr = self.read_epoch_header(RINEXfile)
                if hdr is None:
                    break
                epoch_time, sats, sow = hdr           
                sow_list[i] = sow
                epochs_times.append(epoch_time)
                epochs[i] = epoch_time
                sat_map = np.ones(len(sats)) * -1
                for n, self.sat in enumerate(sats):
                    if self.sat[0] == 'R':
                        sat_map[n] = int(self.sat[1:]) - 1
                    elif self.sat[0] == 'G':
                        sat_map[n] = int(self.sat[1:]) - 1
                    elif self.sat[0] == 'E':
                        sat_map[n] = int(self.sat[1:]) - 1
                obss[i], llis[i], signal_strengths[i] = self.read_obs(RINEXfile, len(sats), sat_map)
                i += 1
                if i >= cs.CHUNK_SIZE:
                    break
                

            #The below code should be commented out, if not interested in the RINEX XYZ coords for producing a map
            ### START
            rx_X = self.X_rx
            rx_Y = self.Y_rx
            rx_Z = self.Z_rx
            
            a = rinex2_nav.RINEX_2_Nav(selected_nav_rinex, sow_list, epochs_times, rx_X, rx_Y, rx_Z)   
            #### END

            return obss[:i], llis[:i], signal_strengths[:i], epochs[:i]


        def read_data(self, RINEXfile, selected_nav_rinex):
            obs_data_chunks = []

            obss, _, _, epochs= self.read_data_chunk(RINEXfile, selected_nav_rinex)
            epochs = epochs.astype(np.int64)
            epochs = np.divide(epochs, float(3600.000))

            if obss.shape[0] == 0:
                print "obss.shape[0] appears false."

            if self.sat[0] == 'R':
                obs_data_chunks.append(pd.Panel(
                    np.rollaxis(obss, 1, 0),
                    items=['R%02d' % d for d in range(1, 33)],
                    major_axis=epochs,
                    minor_axis=self.obs_types
                ).dropna(axis=0, how='all'))#.dropna(axis=2, how='all'))
            elif self.sat[0] == 'G':
                obs_data_chunks.append(pd.Panel(
                    np.rollaxis(obss, 1, 0),
                    items=['G%02d' % d for d in range(1, 33)],
                    major_axis=epochs,
                    minor_axis=self.obs_types
                ).dropna(axis=0, how='all'))#.dropna(axis=2, how='all'))
            elif self.sat[0] == 'E':
                obs_data_chunks.append(pd.Panel(
                    np.rollaxis(obss, 1, 0),
                    items=['E%02d' % d for d in range(1, 33)],
                    major_axis=epochs,
                    minor_axis=self.obs_types
                ).dropna(axis=0, how='all').dropna(axis=2, how='all'))

            self.obs_data_chunks_dataframe = obs_data_chunks[0]
            odc_dimensions = self.obs_data_chunks_dataframe.shape
            self.itemSize = odc_dimensions[0]

            self.windows_or_mac()


        def windows_or_mac(self):
            if sys.platform.startswith('win'):
                self.winpath = r'<your_path>'
                if not os.path.exists(self.winpath):
                    os.makedirs(self.winpath)
            elif sys.platform.startswith('darwin'):
                workspacePath = os.getcwd()
                pathList = workspacePath.split(os.sep)
                if 'iMac' in pathList:
                    self.STECpath = r'<your_path>'
                    self.VTECpath = r'<your_path>'
                    self.vZOOMpath = r'<your_path>'
                    if not os.path.exists(self.STECpath):
                        os.makedirs(self.STECpath)
                    if not os.path.exists(self.VTECpath):
                        os.makedirs(self.VTECpath)
                    if not os.path.exists(self.vZOOMpath):
                        os.makedirs(self.vZOOMpath)
                else:
                    self.STECpath = r'<your_path>'
                    self.VTECpath = r'<your_path>'
                    self.vZOOMpath = r'<your_path>'
                    if not os.path.exists(self.STECpath):
                        os.makedirs(self.STECpath)
                    if not os.path.exists(self.VTECpath):
                        os.makedirs(self.VTECpath)
                    if not os.path.exists(self.vZOOMpath):
                        os.makedirs(self.vZOOMpath)



if __name__ == '__main__':
    RINEX_filename = 'kaik3180_GAL.16o'
    RINEX_navigation = 'brdm3180.16p' 
    elevation_file = 'kaik3180_GPS.ele'
    azimuth_file = 'kaik3180_GPS.azi'
    dcb_file_p1c1 = 'P1C11611.DCB'
    dcb_file_p1p2 = 'P1P21611_ALL.DCB'
    dcb_file_p2c2 = 'P2C21611_RINEX.DCB'
    
    igs_station = elevation_file[:4].upper()
    workspacePath = os.getcwd()
    pathList = workspacePath.split(os.sep)
    
    if 'iMac' in pathList:
        path_to_be_used = '<your_path>' 

    else:
        path_to_be_used = '<your_path>' 
        
    selected_rinex = path_to_be_used + RINEX_filename
    selected_nav_rinex = path_to_be_used + RINEX_navigation
    selected_ele = path_to_be_used + elevation_file
    selected_azi = path_to_be_used + azimuth_file
    selected_dcb_p1c1 = path_to_be_used + dcb_file_p1c1
    selected_dcb_p1p2 = path_to_be_used + dcb_file_p1p2
    selected_dcb_p2c2 = path_to_be_used + dcb_file_p2c2

    r = RINEX_2_Formatting(selected_rinex, selected_ele, igs_station, selected_azi, selected_nav_rinex, selected_dcb_p1c1, selected_dcb_p1p2, selected_dcb_p2c2)


