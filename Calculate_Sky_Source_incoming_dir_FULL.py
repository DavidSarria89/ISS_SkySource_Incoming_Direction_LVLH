from skyfield.api import load
import numpy as np
import datetime
import pyproj
from math import acos, sin, asin
from skyfield.positionlib import position_from_radec
from satellite_coordinates import satellite_coordinates
import sys

###### check python version, required >= 3.6
if sys.version_info[0] != 3 or sys.version_info[1] < 6:
    print("This script requires Python version >= 3.6")
    sys.exit(1)
######

## INPUT
# datetime to calculate ISS's position (latitude longitude altitude)
input_datetime = datetime.datetime(year=2018, month=4, day=26,
                                   hour=2, minute=3, second=15,
                                   microsecond=575003)

# Specification of the source position in the sky (equatorial coordinate system, J2000)
ra_hours = 5.378/ 15.0  # or ra_hours = 47.0 / 60.0 + 30.0 / 3600.0
# One hour of right ascension (1h) is 15°. Since 24x15°=360°, there are 24h of right ascension around the celestial equator.
dec_degrees = -20.392

##

rad_to_deg = 180.0 / np.pi
deg_to_rad = np.pi / 180.0
indices = {"f107s": 72.6, "f107": 71.5, "Ap": 9.5}

def gps_to_ecef_pyproj(lat, lon, alt):
    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    x, y, z = pyproj.transform(lla, ecef, lon, lat, alt, radians=False)
    return x, y, z

def ecef_to_gps_pyproj(X, Y, Z):
    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    lon, lat, alt = pyproj.transform(ecef,lla, X, Y, Z, radians=False)
    return lat, lon, alt


ts = load.timescale(builtin=True)

to = ts.utc(input_datetime.year, input_datetime.month, input_datetime.day, input_datetime.hour,
            input_datetime.minute, input_datetime.second + input_datetime.microsecond / 1.0e6)


## Getting Source's direction vector from Equatorial (J2000) to ITRF/ECEF
earth = 399  # NAIF code for the Earth center of mass

# distance is not important here
Source = position_from_radec(ra_hours, dec_degrees, t=to, center=earth, distance=0.000045, epoch=ts.J2000)
# input ra is in hours, dec in degrees

# ECEF and ITRF frames are the same thing
dir_source_itrf = np.array(Source.itrf_xyz().km)  # this is pointing from Earth center to space

dir_source_itrf = - dir_source_itrf / np.linalg.norm(dir_source_itrf)  # normalized, pointing from space to Earth
print('\n','Source direction unit vector from space to Earth in ECEF/ITRF frame: ')
print(dir_source_itrf,'\n') # ITRF = ECEF

# alternative method (give similar results)
# subpoint = Source.subpoint()
# print('Latitude:', subpoint.latitude.degrees)
# print('Longitude:', subpoint.longitude.degrees)
# print('elevation:', subpoint.elevation.km)  # elevation is not important here
# x, y, z = gps_to_ecef_pyproj(subpoint.latitude.degrees, subpoint.longitude.degrees, subpoint.elevation.km)
# dir_vect = np.array([x, y, z])
# dir_vect = dir_vect / np.linalg.norm(dir_vect)
# print(dir_vect)

iss_coordinates = satellite_coordinates(name='ISS')
lon, lat, alt, v_vec = iss_coordinates.get_satellite_coordinates(input_datetime)
X, Y, Z = gps_to_ecef_pyproj(lat, lon, alt*1000.0)

print('ISS Position in ECEF/ITRF frame (meter): ')
print(X,Y,Z,'\n')
X_ISS=X
Y_ISS=Y
Z_ISS=Z
print('ISS Position in Geographic(geodetic) frame (lat, lon, alt ;; degrees, degrees, km): ')
lat2, lon2, alt2 = ecef_to_gps_pyproj(X, Y, Z)
print(lat, lon, alt,'\n')

## finding minimum altitude along the source's light of sight
alt_list=[]
min_alt=1e12;
time_vect=np.linspace(0.0,20000.0,20000)

idx=0
for ii in range(len(time_vect)):

    dr=500 # m

    dt = time_vect[ii]
    X_p = X_ISS - dt*dir_source_itrf[0]*dr
    Y_p = Y_ISS - dt*dir_source_itrf[1]*dr
    Z_p = Z_ISS - dt*dir_source_itrf[2]*dr

    lat, lon, alt = ecef_to_gps_pyproj(X_p, Y_p, Z_p)
    
    if alt<min_alt:
        min_alt=alt
    
print('Minimum altitude crossed along source\'s light of sight (km): ')
print(min_alt/1000.0,'\n')

## Building rotation matrix from ECEF(ITRF) to ISS's LVLH

LVLH_baseZ = np.array([-1.0 * X, -1.0 * Y, -1.0 * Z])
LVLH_baseZ = LVLH_baseZ / np.linalg.norm(LVLH_baseZ)

v_vec = np.array(v_vec)
v_vec = v_vec / np.linalg.norm(v_vec)
# Y: Z cross velocity
LVLH_baseY = np.cross(LVLH_baseZ, v_vec)
LVLH_baseY = LVLH_baseY / np.linalg.norm(LVLH_baseY)
# X: Y cross Z
LVLH_baseX = np.cross(LVLH_baseY, LVLH_baseZ)
LVLH_baseX = LVLH_baseX / np.linalg.norm(LVLH_baseX)

rotation_mat = np.array([LVLH_baseX, LVLH_baseY, LVLH_baseZ])

print('rotation matrix from ECEF to LVLH: ')
print(rotation_mat,'\n')

##
dir_Source_ecef = rotation_mat.dot(dir_source_itrf)
print('Source direction unit vector in LVLH frame: ')
print(dir_Source_ecef,'\n')
