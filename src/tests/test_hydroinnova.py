import numpy as np
import pandas as pd
from crnpy import crnpy

def hydroinnova_example_mean_value():
    col_names = 'RecordNum,Date Time(UTC),PTB110_mb,P4_mb,P1_mb,T1_C,RH1,T_CS215,RH_CS215,Vbat,N1Cts,N2Cts,N3Cts,N4Cts,N5Cts,N6Cts,N7Cts,N8Cts,N1ETsec,N3ETsec,N5ETsec,N7ETsec,N1T(C),N1RH,N5T(C),N5RH,GpsUTC,LatDec,LongDec,Alt,Qual,NumSats,HDOP,Speed_kmh,COG,SpeedQuality,strDate'.split(
        ',')
    df = pd.read_csv(
        'https://raw.githubusercontent.com/soilwater/crnpy/main/docs/examples/rover/gypsum_transect_01_may_2018.KSU',
        skiprows=20, names=col_names)
    df['LongDec'] = df['LongDec'] * -1  # Raw data is in absolute values

    # Parse timestamps and set as index
    df['timestamp'] = pd.to_datetime(df['Date Time(UTC)'])
    df.set_index(df['timestamp'], inplace=True)

    # Remove rows with missing coordinates
    df['LatDec'].replace(0.0, np.nan, inplace=True)
    df['LongDec'].replace(0.0, np.nan, inplace=True)
    df.dropna(axis=0, subset=['LatDec', 'LongDec'], inplace=True)
    df.reset_index(drop=True, inplace=True)
    df['x'], df['y'] = crnpy.latlon_to_utm(df['LatDec'], df['LongDec'], 14, missing_values=0.0)

    # Define columns names
    counts_colums = ['N1Cts', 'N2Cts', 'N3Cts','N4Cts', 'N5Cts', 'N6Cts', 'N7Cts', 'N8Cts']
    cont_times_col = ['N1ETsec', 'N1ETsec', 'N3ETsec','N3ETsec', 'N5ETsec', 'N5ETsec', 'N7ETsec', 'N7ETsec']

    # Normalize counts to counts/min
    df[counts_colums] = \
        crnpy.normalize_counts(df[counts_colums], \
                               count_time=60, count_times=df[cont_times_col])

    # Compute total neutron counts
    df['total_counts'] = crnpy.compute_total_raw_counts(df[counts_colums])

    # Find stations with cutoff rigidity similar to estimated by lat,lon
    crnpy.find_neutron_detectors(crnpy.cutoff_rigidity(df['LatDec'][0], df['LongDec'][0]),
                                 start_date=df.iloc[0]['timestamp'], end_date=df.iloc[-1]['timestamp'])

    # Download data for one of the similar stations and add to df
    incoming_neutrons = crnpy.get_incoming_neutron_flux(start_date=df.iloc[0]['timestamp'],
                                                        end_date=df.iloc[-1]['timestamp'], station="NEWK",
                                                        utc_offset=-5, expand_window=2)
    # Interpolate incomming flux hourly to measured timestamps (~ every minute)
    df['incoming_flux'] = crnpy.interpolate_incoming_flux(incoming_neutrons, timestamps=df['timestamp']).values

    # Fill NaN values in atmospheric data
    df[['PTB110_mb', 'RH_CS215', 'T_CS215']] = crnpy.fill_missing_atm(df[['PTB110_mb', 'RH_CS215', 'T_CS215']])

    # Correct count by atmospheric variables and incoming flux
    df['corrected_counts'] = crnpy.atm_correction(df['total_counts'], pressure=df['PTB110_mb'], humidity=df['RH_CS215'],
                                                  temp=df['T_CS215'],
                                                  Pref=df['PTB110_mb'].mean(), Aref=0, L=130,
                                                  incoming_neutrons=df['incoming_flux'])

    # Smooth variable
    df['corrected_smoothed'] = crnpy.smooth_2d(df['x'],
                                               df['y'],
                                               df['corrected_counts'],
                                               buffer=800, method='median', rnd=True)

    # Estimate Soil Columetric Water Content
    df['VWC'] = crnpy.counts_to_vwc(df['corrected_smoothed'], N0=550, bulk_density=1.3, Wlat=0.03, Wsoc=0.01)

    # Drop VWC NaN values before interpolating values
    df = df.dropna(subset=['VWC'])

    # Interpolate variable
    X_pred, Y_pred, Z_pred = crnpy.interpolate_2d(df['x'],
                                                  df['y'],
                                                  df['VWC'],
                                                  dx=250, dy=250, method='cubic')

    return X_pred, Y_pred, Z_pred


def test_rover():
    X_pred, Y_pred, Z_pred = hydroinnova_example_mean_value()
    assert np.nanmean(Z_pred) > 0.155 and np.nanmean(Z_pred) < 0.16, f"Rover survey test failed. Mean value is {np.nanmean(Z_pred)}, expected value is between 0.155 and 0.16."