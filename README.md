# spectralanalysisRECCS by Lucas Dimoveo

#This was a project done under the mentorship of Stefan Tulich, PhD, for the RECCS program June-August 2021.

#Apologies for the messiness of the code, as this was my first project. Updates and comments will be improved in the coming weeks

#the netCDF4 file used can be found here: https://psl.noaa.gov/data/gridded/data.interp_OLR.html

#the netCDF4 file has the variable for olr (Outgoing Longwave Radiation) from the Nino 3.4 region of the Pacific Ocean from 1970 to 2021.
#the code takes the data from the olr variable and does the following:
      -- gets statistical data such as the mean, standard deviation, and variance
      -- identifies El Nino, La Nina, and Neutral years 
