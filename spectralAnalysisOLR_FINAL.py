import cftime as nc4
import netCDF4 as nc
import matplotlib
from matplotlib import pyplot as plt
import random as random
import netCDF4
from numpy.ma import MaskError, anomalies
from pylab import *
import math
from math import *
import warnings
import plotly.graph_objects as go
import matplotlib.dates as dates
#from matplotlib.dates import num2date, date2num 
from scipy import signal
warnings.simplefilter('ignore')
import numpy as np
from numpy.fft import fftfreq, rfft, ifft, fft
import scipy as sp
from scipy.fftpack import *
from scipy.signal import butter, filtfilt, freqs, signaltools, windows 
import pandas as pd
from netCDF4 import Dataset, num2date
import datetime
from datetime import timedelta
#from cftime import num2date, date2num #this isn't working for some reason
import sympy 
from sympy import fourier_series, pi
from sympy.abc import Y, x 
from scipy.stats import chisquare, spearmanr, pearsonr
import xarray as xr

#extract data from netCDF file
fn = 'E:/Omen Laptop Backup (Lucas Dimoveo)/reccs/olr.day.mean.nc'
ds = nc.Dataset(fn)

time = ds['time']
time_reduced = time[1800:] #loop find gap in array

years_preFix = 17153/365
print(years_preFix)

day_difference = 17153-15352
year_difference = day_difference/365
print("Difference between original timeseries and time reduced is", day_difference, "days, ", year_difference, "years") #1801 days, 4.934246575342466 years
#.934246575342466  years is 11.210946618140695819 months, or 341.00000000000011369 days into the year. 
julian_to_gregorian_preFix = datetime.date(1974, 6, 1) + datetime.timedelta(days=17153)
print(julian_to_gregorian_preFix)
julian_to_gregorian_Fix = julian_to_gregorian_preFix - datetime.timedelta(days=15352) 
print(julian_to_gregorian_Fix)
startDate = julian_to_gregorian_Fix

#create a loop that creates an array of datetime for range in 0,15352
time_reduced_dates = []
for t in range(0,15352):
  arr = np.array([startDate + datetime.timedelta(days=t)])
  time_reduced_dates.append(arr)
  # time_reduced_dates += arr

olr = np.array(ds['olr'][:])
olr_reduced = olr[1800:,:,:]


olr_average = []
for t in range(0,15352): 
  average = np.average(olr_reduced[t, 35:37, 75:79])
  olr_average.append(average)

olr_reduced_averageOLR = np.array(olr_average)

Fs = 1
Ts = 1.0/Fs
t = arange(0,15352,Ts) 
n = len(olr_reduced_averageOLR)
#print(n)
window = signal.tukey(n, alpha=.10)
k = arange(n)
frq = k/Ts # two sides frequency range
frq = frq[range(math.floor(n/2))] 
olr_reduced_averageOLR = olr_reduced_averageOLR - np.average(olr_reduced_averageOLR)
#print("Average value for olr_reduced_averageOLR is", np.average(olr_reduced_averageOLR)) 
olr_reduced_averageOLR_tapered = window * olr_reduced_averageOLR
tapered_fft = np.fft.fft(olr_reduced_averageOLR_tapered)
Y = tapered_fft/n
Y= Y[range(math.floor(n/2))]

years = n / 365 #n = number of inputs, or 15152 days


olr_mean = np.mean(olr_reduced_averageOLR_tapered)
olr_variance = np.var(olr_reduced_averageOLR_tapered)
olr_std = np.std(olr_reduced_averageOLR_tapered)
#OLR time series has mean:  -0.32251191996630024 a variance:  692.6303354766661 and a standard deviation:  26.317871028574217


#isolate DJF
year_days = 365
years_Array = []
years_Array = [i for i in range((42))]


winter_dates_array = []
winter_OLR_array = np.array([])
start_Date = time_reduced_dates[208] #find type, convert to int
end_Date = time_reduced_dates[299] #find type, convert to int

#label your variables

for i in range(0,42):
  year_day_iterator = year_days * i
  arrayTime = time_reduced_dates[(208 + year_day_iterator) :(299 + year_day_iterator)]
  #winter_dates_array.append(arrayTime)
  winter_dates_array = np.append(winter_dates_array, arrayTime)


for i in range(0,42):
  year_day_iterator = year_days * i
  arrayOLR = olr_reduced_averageOLR_tapered[(208 + year_day_iterator) :(299 + year_day_iterator)]
  winter_OLR_array = np.append(winter_OLR_array, arrayOLR)

print("shape of winter OLR array is: ",shape(winter_OLR_array))
##Calculations for DJF##

winter_mean = np.mean(winter_OLR_array)
winter_var = np.var(winter_OLR_array)
winter_std = np.std(winter_OLR_array)
print("winter mean is: ", winter_mean, "winter standard deviation is: ", winter_std, "winter variance is: ", winter_var)

winter_mean_array_line = np.full(42, winter_mean)
winter_std_array_line = np.full(42, winter_std)

#find djf mean for each year
djf_mean = []
for i in range(0,42):
  year_day_iterator = year_days * i
  meanOLR = np.mean(olr_reduced_averageOLR_tapered[(208 + year_day_iterator) :(299 + year_day_iterator)])
  djf_mean.append(meanOLR)

#find djf std for all years
djf_std = []
for i in range(0,42):
  year_day_iterator = year_days * i
  arrayOLR = np.std(olr_reduced_averageOLR_tapered[(208 + year_day_iterator) :(299 + year_day_iterator)])
  stdOLR = np.std(arrayOLR)
  djf_std.append(stdOLR)

winter_negative_std_array = [(-1) * i for i in winter_std_array_line]

djf_fft = np.fft.fft(winter_OLR_array)
djf_fft_mean = np.mean(djf_fft) 
djf_fft_variance = np.var(djf_fft)
djf_fft_std = np.std(djf_fft)

#find djf mean of each season. You can identify el nino years for years when OLR values are below one standard deviation

djf_fft_season_chunks = []#[djf_fft[i:i + daysInSeason] for i in range(0, len(djf_fft) ,daysInSeason)]
djf_fft_mean_array = []


print("shape of djf_fft is: ",shape(djf_fft))
djf_fft_season_chunks = np.reshape(djf_fft, (42,91))
winter_season_chunks = np.reshape(winter_OLR_array, (42,91))


winter_mean_array = []
for i in winter_season_chunks:
  ave = np.mean(i)
  winter_mean_array = np.append(winter_mean_array, ave)
djf_mean_powerSpectrum = abs(winter_mean_array**2) #fix this #square magnitude of complex values in array, then test array to see if it is real or complex valued

for i in djf_fft_season_chunks:
  ave = np.mean(i)
  djf_fft_mean_array = np.append(djf_fft_mean_array, ave)
djf_powerSpectrum = abs(djf_fft_mean_array**2) #square magnitude of complex values in array

m = len(djf_fft_season_chunks)
k = arange(m)
frq = k/Ts # two sides frequency range
frq = frq[range(math.floor(m/2))] 
L = djf_fft/m
L= L[range(math.floor(m/2))]

plt.title('DJF OLR Power Spectrum from 1979 to 2020')
plt.stem(frq,abs(L),'r')
plt.xlabel('Freq (Cycles per Day)')
plt.ylabel('|fft^2|')
plt.tight_layout()
plt.show()
plt.clf()
plt.cla()
plt.close('all')

winter_negative_std = winter_std * (-1)

winter_std_array_POSITIVE_line = [i for i in winter_std_array_line]

seventyNine_to_twentyTwenty = np.array(range(1979,2021))
#plt.plot(djf_fft_std)
#plt.plot(djf_fft_mean_array)
plt.plot(seventyNine_to_twentyTwenty,winter_mean_array, 'o')
plt.plot(seventyNine_to_twentyTwenty, winter_mean_array_line)
plt.plot(seventyNine_to_twentyTwenty, winter_negative_std_array)
plt.plot(seventyNine_to_twentyTwenty,winter_std_array_POSITIVE_line)
plt.title("Winter OLR Values")
plt.xlabel("DJF Time from 1979 to 2020")
plt.ylabel("OLR Values")
plt.legend("OLR values", "Standard Deviation and Mean" )
plt.show()


print("djf_fft_season_chunks shape is: ", shape(djf_fft_season_chunks))


print("djf OLR has mean: ",winter_mean, "a variance: " , winter_var, "and a standard deviation: ", winter_std )
print("djf fft has mean: ", djf_fft_mean , "a variance: ", djf_fft_variance , "and a standard deviation: ", djf_fft_std)
#djf OLR has mean:  -5.819333865211153 a variance:  987.902855423144 and a standard deviation:  31.43092196266511
#djf fft has mean:  (-2.855400776991091-2.498611818057495e-16j) a variance:  92972.68937367505 and a standard deviation:  304.91423281584457

#plt.title("DJF OLR Histogram")
#plt.hist(winter_OLR_array, bins=100, density=True, edgecolor='black')
#plt.xlabel("DJF OLR Values W/m^2")
#plt.ylabel("DJF OLR Amounts")
#plt.show()
#plt.clf()
#plt.cla()
#plt.close('all')

djf_fft_ELNINO = []
djf_fft_LANINA = []

djf_fft_ELNINO = np.where(djf_fft_mean_array < djf_fft_std, djf_fft_mean_array, np.append(djf_fft_ELNINO, djf_fft_mean_array))
djf_fft_LANINA = np.where(djf_fft_mean_array > djf_fft_std, djf_fft_mean_array, np.append(djf_fft_LANINA, djf_fft_mean_array))

elNino_powerSpectrum = abs(djf_fft_ELNINO**2)
laNina_powerSPectrum = abs(djf_fft_LANINA**2)

winter_std_Positive = winter_std 

djf_anomalous_Range = [winter_mean - winter_std, winter_mean + winter_std_Positive]
djf_fft_anomalous_Range = [djf_fft_mean - djf_fft_std, djf_fft_mean + djf_fft_std]

elNinoLine = djf_fft_mean - djf_fft_std

djf_warm = winter_mean_array > djf_anomalous_Range[1] 
djf_cold = winter_mean_array < djf_anomalous_Range[0]
djf_middle = np.logical_and(winter_mean_array < winter_OLR_array[1], winter_mean_array > winter_OLR_array[0])


elNinoOLR_djf = winter_mean_array[djf_cold]
laNinaOLR_djf = winter_mean_array[djf_warm]
print("elNinoOLR_djf", elNinoOLR_djf) #elNinoOLR_djf [-39.18686349 -59.05102422 -54.3518338  -44.67970544 -38.89612018 -39.26798391 -49.01100544 -70.32523673]
print("laNinaOLR_djf", laNinaOLR_djf) #laNinaOLR_djf [18.49448235 15.33016901 17.52971138 12.88266343 18.76898143 12.6443124313.27458627 19.07623274 14.88367806 16.05585991 13.99969332]

#el nino = index 3,7,12, 15, 18, 23, 30, 36 ... 1982, 1986, 1991, 1994, 1997, 2001, 2004, 2015
#la nina = index is more difficult to figure out ... 4 (1983), 20 (1999)
F = djf_fft_season_chunks[36]/m
F=F[range(math.floor(m/2))]

plt.title('DJF Mean OLR Power Spectrum 2015 (El Nino)') #'DJF Mean OLR Power Spectrum 2015 (El Nino)'
plt.stem(frq,abs(F),'b')
plt.xlabel('Freq (Cycles per Day)')
plt.ylabel('|fft^2|')
plt.tight_layout()
plt.show()

#write as comprehensions instead of for loops ^^^^  
#learn nested comprehensions 









