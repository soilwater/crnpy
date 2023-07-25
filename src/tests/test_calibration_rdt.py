import pandas as pd
import numpy as np
import crnpy
from scipy.optimize import root


# content of test_calibration_rdt.py
def calibration_example():
    # Load the soil samples data and the CRNP dataset using pandas
    url_soil_samples = "https://raw.githubusercontent.com/soilwater/crnpy/main/docs/examples/calibration/soil_data.csv"
    df_soil = pd.read_csv(url_soil_samples)

    # Define start and end of field survey calibration
    calibration_start = pd.to_datetime("2021-10-22 08:00")
    calibration_end = pd.to_datetime("2021-10-22 16:00")

    # Load the station data
    url_station_data = "https://raw.githubusercontent.com/soilwater/crnpy/main/docs/examples/calibration/station_data.csv"
    df_station = pd.read_csv(url_station_data, skiprows=[0, 2, 3])

    #  Parse dates (you can also use the option `parse_dates=['TIMESTAMP]` in pd.read_csv()
    df_station['TIMESTAMP'] = pd.to_datetime(df_station['TIMESTAMP'], format='%Y-%m-%d %H:%M:%S')

    # Select station data only until the calibration date
    # Define date in which the probe was deployed in the field (i.e., first record)
    deployment_date = df_station['TIMESTAMP'].iloc[0]

    # Filter station data from the first record to the end of the field survey calibration
    # This is important since we are considering the incoming flux on the first day as the reference value
    idx_period = (df_station['TIMESTAMP'] >= deployment_date) & (df_station['TIMESTAMP'] <= calibration_end)
    df_station = df_station[idx_period]

    # Compute total neutron counts by adding the counts from both probe detectors
    df_station['total_raw_counts'] = crnpy.total_raw_counts(df_station[['counts_1_Tot', 'counts_2_Tot']],
                                                                    nan_strategy='average')

    # Atmospheric corrections

    # Fill NaN values in atmospheric data
    df_station[['barometric_pressure_Avg', 'relative_humidity_Avg', 'air_temperature_Avg']] = df_station[
        ['barometric_pressure_Avg', 'relative_humidity_Avg', 'air_temperature_Avg']].interpolate(method='pchip',
                                                                                                 limit=24,
                                                                                                 limit_direction='both')

    # Calculate absolute humidity
    df_station['abs_humidity'] = crnpy.abs_humidity(df_station['relative_humidity_Avg'],
                                                             df_station['air_temperature_Avg'])

    # Compute correction factor for atmospheric pressure
    # Reference atmospheric pressure for the location is 976 Pa
    # Using an atmospheric attentuation coefficient of 130 g/cm2
    df_station['fp'] = crnpy.correction_pressure(pressure=df_station['barometric_pressure_Avg'],
                                                 Pref=976, L=130)

    # Compute correction factor for air humidity
    df_station['fw'] = crnpy.correction_humidity(abs_humidity=df_station['abs_humidity'],
                                                 Aref=0)

    # Incoming neutron flux correction


    df_station['total_corrected_neutrons'] = df_station['total_raw_counts'] * df_station['fw'] / (
                df_station['fp'])

    # Compute the weights of each sample for the field average
    nrad_weights = crnpy.nrad_weight(df_station['abs_humidity'].mean(), df_soil['theta_v'],
                                     df_soil['distance_from_station'],
                                     (df_soil['bottom_depth'] + df_soil['top_depth']) / 2,
                                     rhob=df_soil['bulk_density'].mean())

    # Apply distance weights to volumetric water content and bulk density
    field_theta_v = np.sum(df_soil['theta_v'] * nrad_weights)
    field_bulk_density = np.sum(df_soil['bulk_density'] * nrad_weights)

    # Determine the mean corrected counts during the calibration survey
    idx_cal_period = (df_station['TIMESTAMP'] >= calibration_start) & (df_station['TIMESTAMP'] <= calibration_end)
    mean_cal_counts = df_station.loc[idx_cal_period, 'total_corrected_neutrons'].mean()

    print(f"Mean volumetric Water content during calibration survey: {round(field_theta_v, 3)}")
    print(f"Mean corrected counts during calibration: {round(mean_cal_counts)} counts")

    # Define the function for which we want to find the roots
    VWC_func = lambda N0: crnpy.counts_to_vwc(mean_cal_counts, N0, bulk_density=field_bulk_density, Wlat=0.03,
                                              Wsoc=0.01) - field_theta_v

    # Make an initial guess for N0
    N0_initial_guess = 1000

    # Find the root
    sol = int(root(VWC_func, N0_initial_guess).x)

    # Print the solution
    print(f"The solved value for N0 is: {sol}")

    return sol

def test_calibration():
    solved_N0 = calibration_example()
    assert solved_N0 > 2500 and solved_N0 < 2800, f"Calibration test failed. Solved N0: {solved_N0}, expected value between 2500 and 2800."