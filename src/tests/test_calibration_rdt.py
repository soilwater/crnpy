import pandas as pd
import numpy as np
from crnpy import crnpy
from scipy.optimize import root



# content of test_calibration_rdt.py
def calibration_example():
    df_soil = pd.read_csv(
        "https://raw.githubusercontent.com/soilwater/crnpy/main/docs/examples/calibration/soil_data.csv")
    df_station = pd.read_csv(
        "https://raw.githubusercontent.com/soilwater/crnpy/main/docs/examples/calibration/station_data.csv",
        skiprows=[0, 2, 3])

    #  Parse dates
    df_station['TIMESTAMP'] = pd.to_datetime(df_station['TIMESTAMP'])

    # Filter data matching the sampling date
    df_station_calib = df_station[(df_station['TIMESTAMP'] > pd.to_datetime("2021-10-22 08:00")) & (
                df_station['TIMESTAMP'] < pd.to_datetime("2021-10-22 16:00"))].copy()
    df_station_calib['abs_h'] = crnpy.estimate_abs_humidity(df_station_calib['relative_humidity_Avg'],
                                                            df_station_calib['air_temperature_Avg'])
    nrad_weights = crnpy.nrad_weight(df_station_calib['abs_h'].mean(), df_soil['theta_v'],
                                     df_soil['distance_from_station'],
                                     (df_soil['bottom_depth'] + df_soil['top_depth']) / 2, rhob=df_soil['bulk_density'])

    field_theta_v = np.sum(df_soil['theta_v'] * nrad_weights)
    field_bulk_density = np.sum(df_soil['bulk_density'] * nrad_weights)
    print(f"Field Volumetric Water content: {round(field_theta_v, 3)}")
    df = df_station_calib.copy()
    # Set timestamp as index
    df.set_index(df['TIMESTAMP'], inplace=True)

    df['total_counts'] = crnpy.compute_total_raw_counts(df[['counts_1_Tot', 'counts_2_Tot']], nan_strategy='average')
    # Find stations with cutoff rigidity similar to estimated by lat,lon,
    # filtering the time window from experiment setup to the end of the calibration
    crnpy.find_neutron_detectors(crnpy.cutoff_rigidity(39.1, -96.6), start_date = df_station['TIMESTAMP'].iloc[0], end_date = df['TIMESTAMP'].iloc[-1])

    #Download data for one of the similar stations and add to df
    incoming_neutrons = crnpy.get_incoming_neutron_flux(df_station['TIMESTAMP'].iloc[0], df['TIMESTAMP'].iloc[-1], station="DRBS", utc_offset=-5)
    df['incoming_flux']=crnpy.interpolate_incoming_flux(incoming_neutrons, timestamps=df['TIMESTAMP'])
    ref_incoming_flux = incoming_neutrons.iloc[0]
    # Fill NaN values in atmospheric data
    df[['pressure', 'RH', 'T']] = crnpy.fill_missing_atm(
        df[['barometric_pressure_Avg', 'relative_humidity_Avg', 'air_temperature_Avg']])
    # Correct count by atmospheric variables and incoming flux
    df['total_counts'] = crnpy.fill_counts(df['total_counts'])
    df['corrected'] = crnpy.atm_correction(df['total_counts'], pressure=df['pressure'], humidity=df['RH'], temp=df['T'],
                                           Pref=976, Aref=0, L=130, incoming_neutrons=df['incoming_flux'],
                                           incoming_Ref=ref_incoming_flux).dropna()
    print(f"Mean corrected neutron count during sampling: {df['corrected'].mean().round()}")

    # Define the function for which we want to find the roots
    def func(N0):
        return crnpy.counts_to_vwc(df['corrected'].mean(), N0, bulk_density=field_bulk_density, Wlat=0.03,
                                   Wsoc=0.01) - field_theta_v

    # Make an initial guess for N0
    N0_initial_guess = 1000

    # Find the root
    sol = root(func, N0_initial_guess)

    # Print the solution
    print(f"The solved value for N0 is: {sol.x.round()}")
    return sol.x.round().item()


def test_calibration():
    solved_N0 = calibration_example()
    assert solved_N0 > 2500 and solved_N0 < 2800, f"Calibration test failed. Solved N0: {solved_N0}, expected value between 2500 and 2800."