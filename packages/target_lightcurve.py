import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import ascii

def target_lightcurve(observations):

    # Requests user input
    # STELLA zero point is 15.87
    zeropoint = input("Enter the telescope zero point: ")

    # Converts counts to logarithmic magnitudes scale
    log_counts = Table(data = (-2.5*np.log10(observations['Star 1 average']),
                               -2.5*np.log10(observations['Star 2 average']),
                               -2.5*np.log10(observations['Star 3 average']),
                               -2.5*np.log10(observations['Target average'])),
                       names = ('Star 1 log', 'Star 2 log',
                                'Star 3 log', 'Target log'))

    # Calculates the magnitude errors
    log_counts_errors = Table(data = (1.09*(observations['Star 1 error']/
                                      observations['Star 1 average']),
                                1.09*(observations['Star 2 error']/
                                      observations['Star 2 average']),
                                1.09*(observations['Star 3 error']/
                                      observations['Star 3 average']),
                                1.09*(observations['Target error']/
                                      observations['Target average']),),
                     names = ('Star 1 log error', 'Star 2 log error',
                                'Star 3 log error', 'Target log error'))

    # Finds diffential magnitude with respect to star 3
    # Corrects for zero point
    magnitudes = Table(data = (observations['Decimal_year'],
                              log_counts['Star 1 mag'] -
                              log_counts['Star 3 mag'] + zeropoint,
                              log_counts['Star 2 mag'] -
                              log_counts['Star 3 mag'] + zeropoint,
                              log_counts['Star 3 mag'] -
                              log_counts['Star 3 mag'] + zeropoint,
                              log_counts['Target mag'] -
                              log_counts['Star 3 mag'] + zeropoint),
                      names = ('Decimal Year',
                               'Star 1 mag', 'Star 2 mag', 'Star 3 mag',
                               'Target mag'))

    mags_error = Table(data=(np.sqrt(log_counts_errors['Star 1 log error']**2 +
                                     log_counts_errors['Star 3 log error']**2),
                              np.sqrt(log_counts_errors['Star 2 log error']**2 +
                                      log_counts_errors['Star 3 log error']**2),
                              np.sqrt(log_counts_errors['Star 3 log error']**2 +
                                      log_counts_errors['Star 3 log error']**2),
                              np.sqrt(log_counts_errors['Target log error']**2 +
                                      log_counts_errors['Star 3 log error']**2)),
                        names = ('Star 1 mag error', 'Star 2 mag error',
                                 'Star 3 mag error', 'Target mag error'))

    # Save magnitude results to file
    ascii.write(magnitudes, 'magnitudes.csv', format = 'csv', overwrite = True)
    ascii.write(mags_error, 'magnitude_errors.csv', format = 'csv', overwrite = True)

    # Plot the target lightcurve in magnitudes
    plt.gca().invert_yaxis()
    #plt.ylim(17.5, 16.3)
    # TODO add band name eg u'-band
    plt.ylabel('Magnitude')
    plt.xlabel('Decimal Year')
    plt.errorbar(observations['Decimal Year'], magnitudes['Target mag'],
                 yerr = mags_error['Target mag error'], linestyle='None',
                 marker='.', color='purple')
    plt.title('Target lightcurve')
    plt.ticklabel_format(useOffset=False)
    plt.show()
    plt.savefig('target_lightcurve')
    plt.close()
