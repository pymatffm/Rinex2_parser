#!/usr/bin/env python
#*-* coding: utf-8 *-*

# TODO: split header constants in "if" statement so that it is automated chosen
#Header constants
END_OF_HEADER = "END OF HEADER\n"
END_OF_HEADER_2 = "END OF HEADER"
END_OF_HEADER_3 = "END OF HEADER       "#\n"
END_OF_HEADER_4 = "END OF HEADER\r\n"
COMMENT = "COMMENT"
INTERVAL = "INTERVAL\n"
TIME_OF_FIRST_OBS = "TIME OF FIRST OBS\n"
MARKER_NAME = "MARKER NAME\n"
MARKER_NAME_2 = "MARKER NAME"
MARKER_NAME_3 = "MARKER NAME\r\n"
UNKNOWN = "UNKNOWN"
APPROX_POSITION_XYZ = 'APPROX POSITION XYZ\n'
APPROX_POSITION_XYZ_2 = 'APPROX POSITION XYZ '
APPROX_POSITION_XYZ_3 = 'APPROX POSITION XYZ'
TYPES_OF_OBSERV = "# / TYPES OF OBSERV"
TOTAL_SATS = 32
OBSERVE_TYPES_3 = "SYS / # / OBS TYPES\n"
OBSERVE_TYPES_4 = "SYS / # / OBS TYPES\r"
DATE = "TIME OF FIRST OBS"
PHASE_SHIFT = "SYS / PHASE SHIFT\n"

#Mapping function constants
SHELL_HEIGHT = 300  #UNITS: km
RADIUS_EARTH = 6371.0 #UNITS: km
RADIUS_RX = 6371    #UNITS: km
MF_COEFF = RADIUS_RX / (SHELL_HEIGHT + RADIUS_EARTH)

#Other constants
CHUNK_SIZE = 15000
LIGHT_V = 299792458                                         # Vacuum speed of light, m/s
RINEX_NAVIGATION_PARAMETERS = 29 
RINEX_GLONASS_NAVIGATION_PARAMETERS = 15
NAVIGATION_PARAM_LINES = 8
GLONASS_NAVIGATION_PARAM_LINES = 4                                 # Gravitational constant times mass of Earth, m**3 / s**2 
OMEGAE_DOT = 7.2921151467 * (10**-5)                        # Mean angular velocity of Earth, rad / s
INVERSE_OF_FLATTENING = 298.257223563
SEMI_MAJOR_AXIS = 6378137
GM = 3.986004418 * (10**14)
C_20 = -1.08263 * (10**-3)

#Array constants
SEQUENCE_SIZE = 143
SEQUENCE_SIZE_GPS = 20
SEQUENCE_SIZE_GALILEO = 143 #max values
SEQUENCE_SIZE_GLONASS = 60 #max values
SEQUENCE_SIZE_BEIDOU = 26

NAVIGATION_PARAMS_GPS = 29 #? 
NAVIGATION_PARAMS_GALILEO = 28
NAVIGATION_PARAMS_GLONASS = 15
NAVIGATION_PARAMS_BEIDOU = 29

NAVI_PARAMS_LINES_GPS = 8 
NAVI_PARAMS_LINES_GALILEO = 8 
NAVI_PARAMS_LINES_GLONASS = 4 
NAVI_PARAMS_LINES_BEIDOU = 8 

#GPS related constants
L1 = 154*(10.23*10**6)                                      # L1 frequency Hz
L2 = 120*(10.23*10**6)                                      # L2 frequency Hz
GPS_LAMBDA_1 = LIGHT_V / L1                                 # Wavelength on L1:  0.190293672798 m
GPS_LAMBDA_2 = LIGHT_V / L2                                 # Wavelength on L2:  0.244210213425 m
COEFF_GPS = (1 - (L1/L2)**2)                                 # Ionosphere QC coefficient as from Kai Borre

TECU_L1 = (40.3 * (10**16)) / ((L1)**2)                     
TECU_L2 = (40.3 * (10**16)) / ((L2)**2)
TECU_L1L2 = TECU_L2 - TECU_L1

#GALILEO related constants (to be extended depending on Signals used)
E1 = 1575.42*(10**6)  
E5a = 1176.45*(10**6)
E5b = 1207.14*(10**6)
E5 = 1191.795*(10**6)
E6 = 1278.75*(10**6)

TECU_E1 = (40.3 * (10**16)) / ((E1)**2)
TECU_E5a = (40.3 * (10**16)) / ((E5a)**2)
TECU_E5b = (40.3 * (10**16)) / ((E5b)**2)
TECU_E5 = (40.3 * (10**16)) / ((E5)**2)

TECU_E1E5   = TECU_E5 - TECU_E1
TECU_E1E5a = TECU_E5a - TECU_E1
TECU_E1E5b = TECU_E5b - TECU_E1
TECU_E5bE5 = TECU_E5 - TECU_E5b
TECU_E5bE5a = TECU_E5a - TECU_E5b

GALILEO_LAMBDA_E1  = LIGHT_V / E1                           # L1
GALILEO_LAMBDA_E5a = LIGHT_V / E5a                          # L5
GALILEO_LAMBDA_E5b = LIGHT_V / E5b                          # L7
GALILEO_LAMBDA_E5  = LIGHT_V / E5                           # L8

COEFF_GALILEO_GAMMA = (1 - (E1/E5)**2)                        # Ionosphere QC coefficient
COEFF_GALILEO_DELTA = (1 - (E1/E5a)**2)
COEFF_GALILEO_EPSILON = (1 - (E1/E5b)**2)
COEFF_GALILEO_IOTA = (1 - (E5b/E5)**2)
COEFF_GALILEO_KAPPA = (1 - (E5b/E5a)**2)

#BeiDou related constants 
B1 = 1561.098*(10**6)       #E2
B2 = 1207.14*(10**6)        #E5b
B3 = 1268.52*(10**6)        #E6

TECU_B1 = (40.3 * (10**16)) / ((B1)**2)
TECU_B2 = (40.3 * (10**16)) / ((B2)**2)
TECU_B3 = (40.3 * (10**16)) / ((B3)**2)

BEIDOU_LAMBDA_B1 = LIGHT_V / B1
BEIDOU_LAMBDA_B2 = LIGHT_V / B2
BEIDOU_LAMBDA_B3 = LIGHT_V / B3

TECU_B1B2 = TECU_B2 - TECU_B1
TECU_B1B3 = TECU_B3 - TECU_B1

COEFF_BEIDOU_THETA = (1 - (B1/B2)**2)   
COEFF_BEIDOU_OMEGA = (1 - (B1/B3)**2) 


#GLONASS related constants
SEMI_MAJOR_AXIS_GLONASS = 6378.137   # / km
GM_GLONASS = 398600.4418             # / km**3 / s**2
TECU_L1_glo = (40.3 * (10**16)) / ((L1)**2)
TECU_L2_glo = (40.3 * (10**16)) / ((L2)**2)
TECU_L1L2_glo = TECU_L2_glo - TECU_L1_glo

def glonass_frequencies(glonassSV): 
    channelDict = {'R01': 1, 'R02': -4, 'R03': 5, 'R04': 6, 'R05': 1, 'R06': -4, 'R07': 5, 'R08': 6,
                   'R09': -2, 'R10': -7, 'R11': 0, 'R12': -1, 'R13': -2, 'R14': -7, 'R15': 0, 'R16': -1,
                   'R17': -6, 'R18': -3, 'R19': 3, 'R20': 2, 'R21': 4, 'R22': -3, 'R23': 3, 'R24': 2,
                   }

    channelNumber = channelDict[glonassSV]
    centre_L1 = 1602*10**6
    centre_L2 = 1246*10**6
    deltaL1 = 562.5*10**3
    deltaL2 = 437.5*10**3
    L1 = float(centre_L1 + channelNumber*deltaL1)
    L2 = float(centre_L2 + channelNumber*deltaL2)
    glonassL1L2 = [L1, L2]

    TECU_L1_glo = (40.3 * (10**16)) / ((L1)**2)
    TECU_L2_glo = (40.3 * (10**16)) / ((L2)**2)
    TECU_L1L2_glo = TECU_L2_glo - TECU_L1_glo

    BETA = (L1/L2)**2
    COEFF_GLONASS = (1 - BETA)

    return glonassL1L2, COEFF_GLONASS, TECU_L1L2_glo
