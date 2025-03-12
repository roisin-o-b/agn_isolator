import numpy as np
from astropy.table import Table
from astropy.io import ascii

def AverageExposures():

    counts_file = input("Enter the counts file name, including full path: ")

    # Assigns data from the counts file to an array
    counts_table = ascii.read(counts_file)

    # New observation table
    # Each observation refers to the average of all exposures in one night
    observations = Table(names = ('Decimal Year', \
        'Star 1 average', 'Star 1 error', \
        'Star 2 average', 'Star 2 error', \
        'Star 3 average', 'Star 3 error', \
        'Target average', 'Target error'))

    # List of observation dates
    obs_dates = np.unique(counts_table['Decimal Year'])

    for date in obs_dates:
        # Finds all exposures on each observation night
        night_exposures = counts_table['Decimal Year'] == date

        # Calculates the error-weighted averages of all exposures in each night
        # i.e. observations
        star1_avg = np.average(counts_table[night_exposures]['Star 1 counts'], \
            weights= 1/counts_table[night_exposures]['Star 1 error']**2)
        star2_avg = np.average(counts_table[night_exposures]['Star 2 counts'], \
            weights= 1/counts_table[night_exposures]['Star 2 error']**2)
        star3_avg = np.average(counts_table[night_exposures]['Star 3 counts'], \
            weights= 1/counts_table[night_exposures]['Star 3 error']**2)
        target_avg = np.average(counts_table[night_exposures]['Target counts'], \
            weights= 1/counts_table[night_exposures]['Target error']**2)

        # Calculates the errors for each observation
        star1_error = np.sqrt(
            1/np.sum(1/counts_table[night_exposures]['Star 1 error']**2)
            )
        star2_error = np.sqrt(
            1/np.sum(1/counts_table[night_exposures]['Star 2 error']**2)
            )
        star3_error = np.sqrt(
            1/np.sum(1/counts_table[night_exposures]['Star 3 error']**2)
            )
        target_error = np.sqrt(
            1/np.sum(1/counts_table[night_exposures]['Target error']**2)
            )

        # Adds the values for each date to the observations table
        observations.add_row([date, star1_avg, star1_error, star2_avg, star2_error,
                          star3_avg, star3_error, target_avg, target_error])

    ascii.write(observations, 'observations.csv', format = 'csv', overwrite = True)

# TODO change code to update rather than overwrite file
