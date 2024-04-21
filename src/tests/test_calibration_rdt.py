import pandas as pd
import numpy as np
import crnpy
from scipy.optimize import root


# content of test_calibration_rdt.py
def calibration_example():
    # Load the soil samples data and the CRNP dataset using pandas
    url_soil_samples = "https://raw.githubusercontent.com/soilwater/crnpy/main/docs/examples/calibration/soil_data.csv"
    df_soil = pd.read_csv(url_soil_samples)

    df_soil['ID'] = df_soil['latitude'].astype(str) +'_'+ df_soil['longitude'].astype(str)

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
    df_station['total_raw_counts'] = crnpy.total_raw_counts(df_station[['counts_1_Tot', 'counts_2_Tot']])

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
    field_theta_v_1, w = crnpy.nrad_weight(df_station['abs_humidity'].mean(), df_soil['theta_v'],
                                         df_soil['distance_from_station'],
                                         (df_soil['bottom_depth'] + df_soil['top_depth']) / 2,
                                         rhob=df_soil['bulk_density'].mean(), method="Kohli_2015")

    field_theta_v_2, w = crnpy.nrad_weight(df_station['abs_humidity'].mean(), df_soil['theta_v'],
                                         df_soil['distance_from_station'],
                                         (df_soil['bottom_depth'] + df_soil['top_depth']) / 2,
                                         rhob=df_soil['bulk_density'].mean(),
                                         p=df_station['barometric_pressure_Avg'].mean(),
                                         profiles=df_soil['ID'],
                                         method="Schron_2017")

    # Apply distance weights to volumetric water content and bulk density
    field_bulk_density = np.sum(df_soil['bulk_density'].mean())

    # Determine the mean corrected counts during the calibration survey
    idx_cal_period = (df_station['TIMESTAMP'] >= calibration_start) & (df_station['TIMESTAMP'] <= calibration_end)
    mean_cal_counts = df_station.loc[idx_cal_period, 'total_corrected_neutrons'].mean()

    print(f"Mean volumetric Water content during calibration survey: {round(field_theta_v_2, 3)}")
    print(f"Mean corrected counts during calibration: {round(mean_cal_counts)} counts")

    # Define the function for which we want to find the roots
    VWC_func = lambda N0: crnpy.counts_to_vwc(mean_cal_counts, N0, bulk_density=field_bulk_density, Wlat=0.03,
                                              Wsoc=0.01) - field_theta_v_1

    # Make an initial guess for N0
    N0_initial_guess = 1000

    # Find the root
    sol = int(root(VWC_func, N0_initial_guess).x)

    # Print the solution
    print(f"The solved value for N0 is: {sol}")

    # Compute close and far soil moisture to compare both weight methods
    close_sm = df_soil[(df_soil['distance_from_station'] < 10) & (df_soil['bottom_depth'] < 10)]['theta_v'].mean()
    far_sm = df_soil[(df_soil['distance_from_station'] > 10) & (df_soil['bottom_depth'] < 10)]['theta_v'].mean()
    return sol, field_theta_v_1, field_theta_v_2, close_sm, far_sm

def test_calibration():
    solved_N0, field_sm_Kohli, field_sm_Schron, close_sm, far_sm  = calibration_example()
    print("The expected value for N0 is between 2500 and 2800")

    assert solved_N0 > 2500 and solved_N0 < 2800, f"Calibration test failed. Solved N0: {solved_N0}, expected value between 2500 and 2800."
    print("The expected value for field soil moisture using Kohli et al. (2015) or Schron et al. (2017) is between 0.25 and 0.30")
    assert field_sm_Kohli > 0.25 and field_sm_Kohli < 0.30, f"Calibration test failed. Field soil moisture using Kohli et al. (2015): {field_sm_Kohli}, expected value between 0.25 and 0.30."
    assert field_sm_Schron > 0.2 and field_sm_Schron < 0.30, f"Calibration test failed. Field soil moisture using Schron et al. (2017): {field_sm_Schron}, expected value between 0.25 and 0.30."

    print(f"Close soil moisture: {close_sm}")
    print(f"Far soil moisture: {far_sm}")
    print(f"Field soil moisture using Kohli et al. (2015): {field_sm_Kohli}")
    print(f"Field soil moisture using Schron et al. (2017): {field_sm_Schron}")

    # Schron et al. (2017) method should more influenced by close soil moisture samples according to CRNP sensitivity assumptions

    dif_Kohli = abs(field_sm_Kohli - close_sm)
    dif_Schron = abs(field_sm_Schron - close_sm)

    assert dif_Schron < dif_Kohli, f"Calibration test failed. Schron et al. (2017) method should be more influenced by close soil moisture samples according to CRNP sensitivity assumptions."