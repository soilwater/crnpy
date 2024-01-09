import numpy as np
import pandas as pd
import crnpy


def hydroinnova_example_mean_value():
    # Load sample dataset
    filepath = "https://raw.githubusercontent.com/soilwater/crnpy/main/docs/examples/rover/gypsum_transect_01_may_2018.KSU"
    col_names = 'RecordNum,Date Time(UTC),barometric_pressure_Avg,P4_mb,P1_mb,T1_C,RH1,air_temperature_Avg,relative_humidity_Avg,Vbat,N1Cts,N2Cts,N3Cts,N4Cts,N5Cts,N6Cts,N7Cts,N8Cts,N1ETsec,N3ETsec,N5ETsec,N7ETsec,N1T(C),N1RH,N5T(C),N5RH,GpsUTC,LatDec,LongDec,Alt,Qual,NumSats,HDOP,Speed_kmh,COG,SpeedQuality,strDate'.split(
        ',')

    df = pd.read_csv(filepath, skiprows=20, names=col_names)
    df['LongDec'] = df['LongDec'] * -1  # Raw data is in absolute values

    # Parse timestamps and set as index
    df['timestamp'] = pd.to_datetime(df['Date Time(UTC)'])

    # Remove rows with missing coordinates
    df['LatDec'].replace(0.0, np.nan, inplace=True)
    df['LongDec'].replace(0.0, np.nan, inplace=True)
    df.dropna(axis=0, subset=['LatDec', 'LongDec'], inplace=True)
    df.reset_index(drop=True, inplace=True)

    # Convert Lat and Lon to X and Y
    df['x'], df['y'], zone_letter, zone_number = crnpy.latlon_to_utm(df['LatDec'], df['LongDec'])
    df['x'], df['y'] = crnpy.rover_centered_coordinates(df['x'], df['y'])  # Estimate the center of the observation

    # Define columns names
    counts_colums = ['N1Cts', 'N2Cts', 'N3Cts', 'N4Cts', 'N5Cts', 'N6Cts', 'N7Cts', 'N8Cts']
    cont_times_col = ['N1ETsec', 'N1ETsec', 'N3ETsec', 'N3ETsec', 'N5ETsec', 'N5ETsec', 'N7ETsec', 'N7ETsec']

    # Compute total neutron counts
    df['total_raw_counts'] = crnpy.total_raw_counts(df[counts_colums])

    # Fill NaN values in atmospheric data
    df[['barometric_pressure_Avg', 'relative_humidity_Avg', 'air_temperature_Avg']] = df[
        ['barometric_pressure_Avg', 'relative_humidity_Avg', 'air_temperature_Avg']].interpolate(
        method='pchip', limit=24, limit_direction='both')

    # Compute the pressure correction factor
    df['fp'] = crnpy.correction_pressure(pressure=df['barometric_pressure_Avg'],
                                         Pref=df['barometric_pressure_Avg'].mean(), L=130)

    # Estimate the absolute humidity and compute the vapor pressure correction factor
    df['abs_humidity'] = crnpy.abs_humidity(df['relative_humidity_Avg'], df['air_temperature_Avg'])
    df['fw'] = crnpy.correction_humidity(df['abs_humidity'], Aref=0)

    # Apply correction factors
    df['total_corrected_neutrons'] = df['total_raw_counts'] * df['fw'] / (df['fp'])
    # Smooth variable
    df['corrected_neutrons_smoothed'] = crnpy.spatial_average(df['x'],
                                                              df['y'],
                                                              df['total_corrected_neutrons'],
                                                              buffer=800, method='median', rnd=True)

    # Estimate lattice water (%) based on texture, and convert it to decimal.
    lattice_water = crnpy.lattice_water(clay_content=0.35)

    # Estimate Soil Columetric Water Content
    df['VWC'] = crnpy.counts_to_vwc(df['corrected_neutrons_smoothed'], N0=550, bulk_density=1.3, Wlat=lattice_water, Wsoc=0.01)

    # Drop VWC NaN values before interpolating values
    df = df.dropna(subset=['VWC'])

    # Interpolate variable using IDW (https://en.wikipedia.org/wiki/Inverse_distance_weighting)
    X_pred, Y_pred, Z_pred = crnpy.interpolate_2d(df['x'],
                                                  df['y'],
                                                  df['VWC'],
                                                  dx=250, dy=250, method='idw')
    return X_pred, Y_pred, Z_pred


def test_rover():
    print("\nRunning Rover survey test...")
    X_pred, Y_pred, Z_pred = hydroinnova_example_mean_value()
    expected_min = 0.16
    expected_max = 0.17
    print(f"Rover survey test passed. Mean value is {np.nanmean(Z_pred)}, expected value is between {expected_min} and {expected_max}")
    assert np.nanmean(Z_pred) > expected_min and np.nanmean(Z_pred) < expected_max, f"Rover survey test failed. Mean value is {np.nanmean(Z_pred)}, expected value is between {expected_min} and {expected_max}"