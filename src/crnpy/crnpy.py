# crnpy/crnpy.py
"""
`crnpy` is a Python package for processing cosmic ray neutron data.

 Created by Joaquin Peraza and Andres Patrignani.
"""

import crnpy.data as data
import io
import numbers
import numpy as np
import pandas as pd
import requests
import sys
import time
import utm
import warnings

from scipy.interpolate import griddata
from scipy.signal import savgol_filter
from scipy.special import erfcinv

# Define python version
python_version = (3, 8)  # tuple of (major, minor) version requirement
python_version_str = str(python_version[0]) + "." + str(python_version[1])

# produce an error message if the python version is less than required
if sys.version_info < python_version:
    msg = "Module only runs on python version >= %s" % python_version_str
    raise Exception(msg)


def remove_incomplete_intervals(df, timestamp_col, integration_time, remove_first=False):
    """Function that removes rows with incomplete integration intervals.
    
    Args:
        df (pandas.DataFrame): Pandas Dataframe with data from stationary or roving CRNP devices.
        timestamp_col (str): Name of the column with timestamps in datetime format.
        integration_time (int): Duration of the neutron counting interval in seconds. Typical values are 60 seconds and 3600 seconds.
        remove_first (bool, optional): Remove first row. Default is False.
        
    Returns:
        (pandas.DataFrame): 
    """

    # Check format of timestamp column
    if not pd.api.types.is_datetime64_any_dtype(df[timestamp_col]):
        raise TypeError('timestamp_col must be datetime64. Use `pd.to_datetime()` to fix this issue.')

    # Check if differences in timestamps are below or above the provided integration time
    idx_delta = df[timestamp_col].diff().dt.total_seconds() != integration_time

    # The first row has no preceding interval (NaT difference), so keep it unless remove_first is requested
    idx_delta.iloc[0] = bool(remove_first)

    # Select rows that meet the specified integration time (copy so later column assignments never act on a view)
    df = df[~idx_delta].copy()
    df.reset_index(drop=True, inplace=True)

    # Notify user about the number of rows that have been removed
    print(f"Removed a total of {sum(idx_delta)} rows.")

    return df


def fill_missing_timestamps(df, timestamp_col='timestamp', freq='h', round_timestamp=True, verbose=False):
    """Helper function to fill rows with missing timestamps in datetime record. Rows are filled with NaN values.

     Args:
         df (pandas.DataFrame): Pandas DataFrame.
         timestamp_col (str, optional): Column with the timestamp. Must be in datetime format. Default column name is 'timestamp'.
         freq (str, optional): Timestamp frequency using pandas offset aliases: 'h' for hourly, 'min' for minute,
             or e.g. '3h' for a 3 hour frequency. Default is 'h'.
         round_timestamp (bool, optional): Whether to round timestamps to the nearest frequency. Default is True.
         verbose (bool, optional): Prints the missing timestamps added to the DataFrame.

     Returns:
         (pandas.DataFrame): DataFrame with filled missing timestamps.

     """

    # Check format of timestamp column
    if not pd.api.types.is_datetime64_any_dtype(df[timestamp_col]):
        raise TypeError('timestamp_col must be datetime64. Use `pd.to_datetime()` to fix this issue.')

    # Round timestamps to nearest frequency. This step must precede the filling of rows.
    if round_timestamp:
        df[timestamp_col] = df[timestamp_col].dt.round(freq)

    # Fill in rows with missing timestamps
    start_date = df[timestamp_col].iloc[0]
    end_date = df[timestamp_col].iloc[-1]
    date_range = pd.date_range(start_date, end_date, freq=freq)
    counter = 0
    for date in date_range:
        if date not in df[timestamp_col].values:
            if verbose:
                print('Adding missing date:', date)
            new_line = pd.DataFrame({timestamp_col: date}, index=[-1])  # By default fills columns with np.nan
            df = pd.concat([df, new_line])
            counter += 1

    df.sort_values(by=timestamp_col, inplace=True)
    df.reset_index(drop=True, inplace=True)

    # Notify user about the number of rows that have been removed
    print(f"Added a total of {counter} missing timestamps.")

    return df


def total_raw_counts(counts):
    """Compute the sum of uncorrected neutron counts for all detectors.

    Args:
        counts (pandas.DataFrame): Dataframe containing only the columns with neutron counts.

    Returns:
        (pandas.DataFrame): Dataframe with the sum of uncorrected neutron counts for all detectors.
    """

    # With more than one detector, fill missing values of a detector with the mean of the other detectors in that row
    if counts.shape[1] > 1:
        counts = counts.apply(lambda x: x.fillna(counts.mean(axis=1)), axis=0)

    # Compute sum of counts
    total_raw_counts = counts.sum(axis=1)

    # Replace zeros with NaN
    total_raw_counts = total_raw_counts.replace(0, np.nan)

    return total_raw_counts


def is_outlier(x, method, window=11, min_val=None, max_val=None):
    """Function that tests whether values are outliers using a range check and/or a dispersion-based method.

    Args:
        x (pd.Series): Series with the variable to test, typically neutron counts.
        method (str): Outlier detection method. One of: range, iqr, moviqr, zscore, movzscore, modified_zscore, and scaled_mad.
            The range check defined by `min_val` and `max_val` is always applied in addition to the selected method.
        window (int, optional): Window size for the moving central tendency. Default is 11.
        min_val (int or float): Minimum value for a reading to be considered valid. Default is None.
        max_val(int or float): Maximum value for a reading to be considered valid. Default is None.

    Returns:
        (pandas.Series): Boolean indicating outliers.

    References:
        Iglewicz, B. and Hoaglin, D.C., 1993. How to detect and handle outliers (Vol. 16). Asq Press.
    """

    if not isinstance(x, pd.Series):
        raise TypeError('x must be of type pandas.Series')

    # Separate this method to allow usage together with other methods below
    if isinstance(min_val, numbers.Number) and isinstance(max_val, numbers.Number):
        idx_range_outliers = (x < min_val) | (x > max_val)
    else:
        idx_range_outliers = pd.Series(False, index=x.index)

    # Apply other methods in addition to a range check
    if method == 'range':
        if not (isinstance(min_val, numbers.Number) and isinstance(max_val, numbers.Number)):
            raise ValueError("The 'range' method requires min_val and max_val.")
        idx_outliers = idx_range_outliers

    elif method == 'iqr':
        q1 = x.quantile(0.25)
        q3 = x.quantile(0.75)
        iqr = q3 - q1
        high_fence = q3 + (1.5 * iqr)
        low_fence = q1 - (1.5 * iqr)
        idx_outliers = (x < low_fence) | (x > high_fence)

    elif method == 'moviqr':
        q1 = x.rolling(window, center=True).quantile(0.25)
        q3 = x.rolling(window, center=True).quantile(0.75)
        iqr = q3 - q1
        ub = q3 + (1.5 * iqr)  # Upper boundary
        lb = q1 - (1.5 * iqr)  # Lower boundary
        idx_outliers = (x < lb) | (x > ub)

    elif method == 'zscore':
        zscore = (x - x.mean()) / x.std()
        idx_outliers = (zscore < -3) | (zscore > 3)

    elif method == 'movzscore':
        movmean = x.rolling(window=window, center=True).mean()
        movstd = x.rolling(window=window, center=True).std()
        movzscore = (x - movmean) / movstd
        idx_outliers = (movzscore < -3) | (movzscore > 3)

    elif method == 'modified_zscore':
        # Compute median absolute difference
        movmedian = x.rolling(window, center=True).median()
        abs_diff = np.abs(x - movmedian)
        mad = abs_diff.rolling(window, center=True).median()

        # Compute modified z-score
        modified_z_score = 0.6745 * abs_diff / mad
        idx_outliers = (modified_z_score < -3.5) | (modified_z_score > 3.5)

    elif method == 'scaled_mad':
        # Returns true for elements more than three scaled MAD from the median (as in MATLAB isoutlier).
        c = -1 / (np.sqrt(2) * erfcinv(3 / 2))
        median = np.nanmedian(x)
        mad = c * np.nanmedian(np.abs(x - median))
        idx_outliers = np.abs(x - median) > 3 * mad

    else:
        raise ValueError('Outlier detection method not found.')

    return idx_outliers | idx_range_outliers


def correction_pressure(pressure, Pref, L):
    r"""Correction factor for atmospheric pressure.

    This function corrects neutron counts for atmospheric pressure using the method described in Andreasen et al. (2017).
    The correction is performed using the following equation:

    $$
    C_{corrected} = \frac{C_{raw}}{fp}
    $$

    where:

    - Ccorrected: corrected neutron counts
    - Craw: raw neutron counts
    - fp: pressure correction factor

    $$
    fp = e^{\frac{P_{ref} - P}{L}}
    $$

    where:

    - P: atmospheric pressure
    - Pref: reference atmospheric pressure
    - L: mass attenuation length for high-energy neutrons.


    Args:
        pressure (list or array): Atmospheric pressure readings, in the same units as Pref and L (e.g. hPa).
        Pref (float): Reference atmospheric pressure. The long-term average pressure at the site is recommended (Zreda et al., 2012).
        L (float): Mass attenuation length for high-energy neutrons, in the same units as the pressure (hPa) or in g cm-2.
            It varies from about 128 g cm-2 at high latitudes to 142 g cm-2 at the equator (Zreda et al., 2012).

    Returns:
        (list): fp pressure correction factor.

    References:
        Zreda, M., Shuttleworth, W. J., Zeng, X., Zweck, C., Desilets, D., Franz, T., and Rosolem, R.: COSMOS: the COsmic-ray Soil Moisture Observing System, Hydrol. Earth Syst. Sci., 16, 4079–4099, https://doi.org/10.5194/hess-16-4079-2012, 2012.

        M. Andreasen, K.H. Jensen, D. Desilets, T.E. Franz, M. Zreda, H.R. Bogena, and M.C. Looms. 2017. Status and perspectives on the cosmic-ray neutron method for soil moisture estimation and other environmental science applications. Vadose Zone J. 16(8). doi:10.2136/vzj2017.04.0086
    """

    # Compute pressure correction factor
    fp = np.exp((Pref - pressure) / L)  # Zreda et al. (2012) Eq. 5

    return fp


def correction_humidity(abs_humidity, Aref):
    r"""Correction factor for absolute humidity.

    This function corrects neutron counts for absolute humidity using the method described in Rosolem et al. (2013) and Andreasen et al. (2017). The correction is performed using the following equation:

    $$
    C_{corrected} = C_{raw} \cdot f_w
    $$

    where:

    - Ccorrected: corrected neutron counts
    - Craw: raw neutron counts
    - fw: absolute humidity correction factor

    $$
    f_w = 1 + 0.0054(A - A_{ref})
    $$

    where:

    - A: absolute humidity
    - Aref: reference absolute humidity

    Args:
        abs_humidity (list or array): Absolute humidity readings (g/m^3). See `abs_humidity()`.
        Aref (float): Reference absolute humidity (g/m^3). The value on the day of the instrument calibration is recommended.

    Returns:
        (list): fw correction factor.

    References:
        Rosolem, R., W. J. Shuttleworth, M. Zreda, T. E. Franz, X. Zeng, and S. A. Kurc, 2013: The Effect of Atmospheric Water Vapor on Neutron Count in the Cosmic-Ray Soil Moisture Observing System. J. Hydrometeor., 14, 1659–1671, https://doi.org/10.1175/JHM-D-12-0120.1.

        M. Andreasen, K.H. Jensen, D. Desilets, T.E. Franz, M. Zreda, H.R. Bogena, and M.C. Looms. 2017. Status and perspectives on the cosmic-ray neutron method for soil moisture estimation and other environmental science applications. Vadose Zone J. 16(8). doi:10.2136/vzj2017.04.0086
    """
    A = abs_humidity
    fw = 1 + 0.0054 * (A - Aref)  # Rosolem et al. (2013); Andreasen et al. (2017) Eq. 3
    return fw


def correction_incoming_flux(incoming_neutrons, incoming_Ref=None, fill_na=None, Rc_method=None, Rc_site=None,
                             site_atmdepth=None, Rc_ref=None, ref_atmdepth=None):

    r"""Correction factor for incoming neutron flux.

    This function corrects neutron counts for incoming neutron flux using the method described in Andreasen et al. (2017). The correction is performed using the following equation:

    $$
    C_{corrected} = \frac{C_{raw}}{f_i}
    $$

    where:

    - Ccorrected: corrected neutron counts
    - Craw: raw neutron counts
    - fi: incoming neutron flux correction factor

    $$
    f_i = \frac{I}{I_{ref}}
    $$

    where:

    - I: incoming neutron flux
    - Iref: reference incoming neutron flux

    Args:
        incoming_neutrons (list or array): Incoming neutron flux readings.
        incoming_Ref (float): Reference incoming neutron flux. Baseline incoming neutron flux.
        fill_na (float): Value to fill missing data. If None, missing data remains as NaN.
        Rc_method (str): Optional to correct for differences in cutoff rigidity between the site and the reference station. Possible values are 'McJannetandDesilets2023' or 'Hawdonetal2014'. If None, no correction is performed.
        Rc_site (float): Cutoff rigidity at the monitoring site.
        site_atmdepth (float): Atmospheric depth at the monitoring site.
        Rc_ref (float): Cutoff rigidity at the reference station.
        ref_atmdepth (float): Atmospheric depth at the reference station.

    Returns:
        (list): fi correction factor.

    References:
        Hawdon, A., D. McJannet, and J. Wallace (2014), Calibration and correction procedures for cosmic-ray neutron soil moisture probes located across Australia, Water Resour. Res., 50, 5029–5043, doi:10.1002/2013WR015138.

        M. Andreasen, K.H. Jensen, D. Desilets, T.E. Franz, M. Zreda, H.R. Bogena, and M.C. Looms. 2017. Status and perspectives on the cosmic-ray neutron method for soil moisture estimation and other environmental science applications. Vadose Zone J. 16(8). doi:10.2136/vzj2017.04.0086

        McJannet, D. L., & Desilets, D. (2023). Incoming neutron flux corrections for cosmic-ray soil and snow sensors using the global neutron monitor network. Water Resources Research, 59, e2022WR033889. https://doi.org/10.1029/2022WR033889
    """
    if incoming_Ref is None and not isinstance(incoming_neutrons, type(None)):
        incoming_Ref = incoming_neutrons[0]
        warnings.warn('Reference incoming neutron flux not provided. Using first value of incoming neutron flux.')
    fi = incoming_neutrons / incoming_Ref

    if Rc_method is not None:
        if Rc_ref is None:
            raise ValueError('Reference cutoff rigidity not provided.')
        if Rc_site is None:
            raise ValueError('Site cutoff rigidity not provided.')

        if Rc_method == 'McJannetandDesilets2023':
            tau = location_factor(site_atmdepth, Rc_site, ref_atmdepth, Rc_ref)
            # McJannet and Desilets (2023) Eq. 10 give the multiplicative factor [tau*I/Iref + 1 - tau]^-1.
            # This library divides counts by fi (C_corrected = C_raw / fi), so fi is the bracket itself.
            fi = tau * fi + 1 - tau

        elif Rc_method == 'Hawdonetal2014':
            Rc_corr = -0.075 * (Rc_site - Rc_ref) + 1.0
            fi = (fi - 1.0) * Rc_corr + 1.0

        else:
            raise ValueError(
                'Cutoff rigidity method not found. Valid options are: McJannetandDesilets2023, Hawdonetal2014.')

    if fill_na is not None:
        fi.fillna(fill_na, inplace=True)  # Use a value of 1 for days without data

    return fi


def get_incoming_neutron_flux(start_date, end_date, station, utc_offset=0, expand_window=0, verbose=False):
    """Function to retrieve neutron flux from the Neutron Monitor Database.

    Args:
        start_date (datetime): Start date of the time series.
        end_date (datetime): End date of the time series.
        station (str): Neutron Monitor station to retrieve data from.
        utc_offset (int): UTC offset in hours. Default is 0.
        expand_window (int): Number of hours to expand the time window to retrieve extra data. Default is 0.
        verbose (bool): Print information about the request. Default is False.

    Returns:
        (pandas.DataFrame): Neutron flux in counts per hour and timestamps.

    Note:
        Data retrieved via NMDB are the property of the individual data providers and are free for non-commercial use
        within the restrictions imposed by the providers. Please acknowledge the NMDB database (www.nmdb.eu), founded
        under the European Union's FP7 programme (contract no. 213007), and the PIs of the individual neutron monitors.
        The acknowledgement text is printed when `verbose=True`.

    References:
        Documentation available:https://www.nmdb.eu/nest/help.php#howto
    """

    # Example: get_incoming_neutron_flux(station='IRKT',start_date='2020-04-10 11:00:00',end_date='2020-06-18 17:00:00')
    # Template url = 'http://nest.nmdb.eu/draw_graph.php?formchk=1&stations[]=KERG&output=ascii&tabchoice=revori&dtype=corr_for_efficiency&date_choice=bydate&start_year=2009&start_month=09&start_day=01&start_hour=00&start_min=00&end_year=2009&end_month=09&end_day=05&end_hour=23&end_min=59&yunits=0'

    # Expand the time window by `expand_window` hours to ensure extra observations are included in the request.
    start_date -= pd.Timedelta(hours=expand_window)
    end_date += pd.Timedelta(hours=expand_window)

    # Convert local time to UTC
    start_date = start_date - pd.Timedelta(hours=utc_offset)
    end_date = end_date - pd.Timedelta(hours=utc_offset)
    root = 'http://www.nmdb.eu/nest/draw_graph.php?'
    url_par = ['formchk=1',
               'stations[]=' + station,
               'output=ascii',
               'tabchoice=revori',
               'dtype=corr_for_efficiency',
               'tresolution=' + str(60),
               'date_choice=bydate',
               'start_year=' + str(start_date.year),
               'start_month=' + str(start_date.month),
               'start_day=' + str(start_date.day),
               'start_hour=' + str(start_date.hour),
               'start_min=' + str(start_date.minute),
               'end_year=' + str(end_date.year),
               'end_month=' + str(end_date.month),
               'end_day=' + str(end_date.day),
               'end_hour=' + str(end_date.hour),
               'end_min=' + str(end_date.minute),
               'yunits=0']

    url = root + '&'.join(url_par)

    if verbose:
        print(f"Retrieving data from {url}")

    # The NMDB draw_graph.php endpoint occasionally returns the HTML page without the ASCII data block,
    # even when data are available, and appears to rate-limit rapid requests. Try at most 3 times, pausing
    # between attempts, when the 'RCORR_E' marker is missing.
    n_attempts = 3
    wait_seconds = 5
    for attempt in range(n_attempts):
        r = requests.get(url).content.decode('utf-8')
        if "RCORR_E" in r:
            break
        if attempt < n_attempts - 1:
            print(f"Request to NMDB unsuccessful, retrying in {wait_seconds} seconds "
                  f"(attempt {attempt + 1} of {n_attempts})...")
            time.sleep(wait_seconds)

    # If the data block is still missing after the retries, no data were returned for the request. The NMDB
    # server appears to limit the number of requests per minute, so waiting before trying again usually helps.
    if "RCORR_E" not in r:
        print("Could not retrieve data from NMDB after "
              f"{n_attempts} attempts. The server may be limiting requests; please wait a minute and try again.")
        if verbose:
            print(f"Error retrieving data from {url}")
        return None

    start = r.find("RCORR_E\n") + 8
    end = r.find('\n</code></pre><br>Total') - 1
    s = r[start:end]
    s2 = ''.join([row.replace(';', ',') for row in s])
    try:
        df_flux = pd.read_csv(io.StringIO(s2), names=['timestamp', 'counts'])
    except:
        if verbose:
            print(f"Error retrieving data from {url}")
        return None

    # Check if all values from selected detector are NaN. If yes, warn the user
    if df_flux['counts'].isna().all():
        warnings.warn('Data for selected neutron detectors appears to be unavailable for the selected period')

    # Convert timestamp to datetime and apply UTC offset
    df_flux['timestamp'] = pd.to_datetime(df_flux['timestamp'])
    df_flux['timestamp'] = df_flux['timestamp'] + pd.Timedelta(hours=utc_offset)

    # Print acknowledgement to inform users about restrictions and to acknowledge the NMDB database
    acknowledgement = """Data retrieved via NMDB are the property of the individual data providers. These data are free for non commercial
use to within the restriction imposed by the providers. If you use such data for your research or applications, please acknowledge
the origin by a sentence like 'We acknowledge the NMDB database (www.nmdb.eu) founded under the European Union's FP7 programme 
(contract no. 213007), and the PIs of individual neutron monitors at: IGY Jungfraujoch 
(Physikalisches Institut, University of Bern, Switzerland)"""
    if verbose:
        print(acknowledgement)

    return df_flux


def get_reference_neutron_flux(station, date=pd.to_datetime("2011-05-01")):
    """Function to retrieve reference neutron flux from the Neutron Monitor Database. Default date is 2011-05-01, following previous studies (Zreda et a., 2012, Hawdon et al., 2014, Bogena et al., 2022).

    Args:
        station (str): Neutron Monitor station to retrieve data from.
        date (datetime): Date of the reference neutron flux. Default is 2011-05-01.

    Returns:
        (float): Reference neutron flux in counts per hour.

    References:
        Zreda, M., Shuttleworth, W. J., Zeng, X., Zweck, C., Desilets, D., Franz, T., & Rosolem, R. (2012). COSMOS: The cosmic-ray soil moisture observing system. Hydrology and Earth System Sciences, 16(11), 4079-4099.

        Hawdon, A., McJannet, D., & Wallace, J. (2014). Calibration and correction procedures for cosmic‐ray neutron soil moisture probes located across Australia. Water Resources Research, 50(6), 5029-5043.

        Bogena, H. R., Schrön, M., Jakobi, J., Ney, P., Zacharias, S., Andreasen, M., ... & Vereecken, H. (2022). COSMOS-Europe: a European network of cosmic-ray neutron soil moisture sensors. Earth System Science Data, 14(3), 1125-1151.

"""

    # Get flux for 2011-05-01
    df_flux = get_incoming_neutron_flux(date, date + pd.Timedelta(hours=24), station=station)
    if df_flux is None:
        warnings.warn(f"Reference neutron flux for {station} not available. Returning NaN.")
    else:
        return df_flux['counts'].median()


def smooth_1d(values, window=5, order=3, method='moving_median'):
    """Use a Savitzky-Golay filter to smooth the signal of corrected neutron counts or another one-dimensional array (e.g. computed volumetric water content).

    Args:
        values (pd.DataFrame or pd.Series): Dataframe containing the values to smooth.
        window (int): Window size for the Savitzky-Golay filter. Default is 5.
        method (str): Method to use for smoothing the data. Default is 'moving_median'.
            Options are 'moving_average', 'moving_median' and 'savitzky_golay'.
        order (int): Order of the Savitzky-Golay filter. Default is 3.

    Returns:
        (pd.DataFrame): DataFrame with smoothed values.

    References:
        Franz, T.E., Wahbi, A., Zhang, J., Vreugdenhil, M., Heng, L., Dercon, G., Strauss, P., Brocca, L. and Wagner, W., 2020.
        Practical data products from cosmic-ray neutron sensing for hydrological applications. Frontiers in Water, 2, p.9.
        doi.org/10.3389/frwa.2020.00009

        Savitzky, A., & Golay, M. J. (1964). Smoothing and differentiation of data by simplified least squares procedures.
        Analytical chemistry, 36(8), 1627-1639.
    """

    if not isinstance(values, pd.Series) and not isinstance(values, pd.DataFrame):
        raise ValueError('Input must be a pandas Series or DataFrame')

    if method == 'moving_average':
        corrected_counts = values.rolling(window=window, center=True, min_periods=1).mean()
    elif method == 'moving_median':
        corrected_counts = values.rolling(window=window, center=True, min_periods=1).median()

    elif method == 'savitzky_golay':
        if np.any(values.isna().values):
            print('Dataframe contains NaN values. Please remove NaN values before smoothing the data.')

        if type(values) == pd.core.series.Series:
            filtered = savgol_filter(values, window, order)
            corrected_counts = pd.DataFrame(filtered, columns=['smoothed'], index=values.index)
        elif type(values) == pd.core.frame.DataFrame:
            corrected_counts = values.copy()
            for col in values.columns:
                corrected_counts[col] = savgol_filter(values[col], window, order)
    else:
        raise ValueError(
            'Invalid method. Please select a valid filtering method., options are: moving_average, moving_median, savitzky_golay')
    corrected_counts = corrected_counts.ffill(limit=window).bfill(limit=window).copy()
    return corrected_counts


def correction_bwe(counts, bwe, r2_N0=0.0053):
    """Function to correct for biomass effects in neutron counts.
    following the approach described in Baatz et al., 2015.

    Args:
        counts (array or pd.Series or pd.DataFrame): Array of epithermal neutron counts.
        bwe (float): Biomass water equivalent kg m-2.
        r2_N0 (float): Ratio r2/N0 of Baatz et al. (2015), i.e. the fractional reduction in neutron counts per kg m-2 of
            biomass water equivalent. Default is 0.0053 (r2 = 6.4 cph per kg m-2 BWE and N0 = 1210 cph, about 0.5% per kg m-2).

    Returns:
        (array or pd.Series or pd.DataFrame): Array of corrected neutron counts for biomass effects.

    References:
        Baatz, R., H. R. Bogena, H.-J. Hendricks Franssen, J. A. Huisman, C. Montzka, and H. Vereecken (2015),
        An empirical vegetation correction for soil water content quantification using cosmic ray probes,
        Water Resour. Res., 51, 2030–2046, doi:10.1002/2014WR016443.
    """

    return counts / (1 - bwe * r2_N0)


def biomass_to_bwe(biomass_dry, biomass_fresh, fWE=0.494):
    """Function to convert biomass to biomass water equivalent.

    Args:
        biomass_dry (array or pd.Series or pd.DataFrame): Above ground dry biomass in kg m-2.
        biomass_fresh (array or pd.Series or pd.DataFrame): Above ground fresh biomass in kg m-2.
        fWE (float): Stoichiometric ratio of H2O to organic carbon molecules in the plant (assuming this is mostly cellulose)
            Default is 0.494 (Wahbi & Avery, 2018).

    Returns:
        (array or pd.Series or pd.DataFrame): Biomass water equivalent in kg m-2.

    References:
        Wahbi, A., Avery, W. (2018). In Situ Destructive Sampling. In:
        Cosmic Ray Neutron Sensing: Estimation of Agricultural Crop Biomass Water Equivalent.
        Springer, Cham. https://doi.org/10.1007/978-3-319-69539-6_2
    """
    return (biomass_fresh - biomass_dry) + fWE * biomass_dry


def correction_road(counts, theta_N, road_width, road_distance=0.0, theta_road=0.12, p0=0.42, p1=0.5, p2=1.06, p3=4,
                    p4=0.16, p5=0.39, p6=0.94, p7=1.10, p8=2.70, p9=0.01):
    """Function to correct for road effects in neutron counts.
    following the approach described in Schrön et al., 2018. The parameters p0 to p9 of the correction function
    default to the values of Table 1 in Schrön et al. (2018): p0 and p1 for the geometry term, p2 to p5 for the
    moisture term, and p6 to p9 for the distance term.

    Args:
        counts (array or pd.Series or pd.DataFrame): Array of epithermal neutron counts.
        theta_N (float): Volumetric water content of the soil estimated from the uncorrected neutron counts.
        road_width (float): Width of the road in m.
        road_distance (float): Distance of the road from the sensor in m. Default is 0.0.
        theta_road (float): Volumetric water content of the road. Default is 0.12.

    Returns:
        (array or pd.Series or pd.DataFrame): Array of corrected neutron counts for road effects.

    References:
        Schrön, M., Rosolem, R., Köhli, M., Piussi, L., Schröter, I., Iwema, J., et al. (2018). Cosmic-ray neutron rover surveys
        of field soil moisture and the influence of roads. Water Resources Research, 54, 6441–6459.
        https://doi.org/10.1029/2017WR021719
    """
    road_width = np.asarray(road_width, dtype=float)
    F1 = p0 * (1 - np.exp(-p1 * road_width))
    # Moisture term F2' of Schrön et al. (2018), Eq. 6 and Table 1, which uses the soil moisture
    # estimated from the uncorrected counts instead of prior knowledge of the field soil moisture.
    F2 = p2 - p3 * theta_road - (p4 + theta_road) / (p5 + theta_N)
    # For a road width of zero F1 = 0 and no correction is applied (Schrön et al., 2018, requirement 4);
    # the width is replaced by 1 m inside F3 only to avoid raising zero to a negative power.
    width_safe = np.where(road_width > 0, road_width, 1.0)
    F3 = p6 * np.exp(-p7 * (width_safe ** -p8) * road_distance ** 4) + (1 - p6) * np.exp(-p9 * road_distance)

    C_roads = 1 + F1 * F2 * F3

    corrected_counts = counts / C_roads

    return corrected_counts


def counts_to_vwc(counts, N0, Wlat, Wsoc, bulk_density, a0=0.0808, a1=0.372, a2=0.115):
    r"""Function to convert corrected and filtered neutron counts into volumetric water content.

    This method implements soil moisture estimation using the non-linear relationship between neutron count and soil water content
    of Desilets et al. (2010), extended with the lattice water and soil organic matter terms and the conversion to volumetric
    units as in Eq. 7 of Hawdon et al. (2014):

    $\theta_v = \left( \frac{a_0}{N/N_0 - a_1} - a_2 - W_{lat} - W_{soc} \right) \rho_b$

    Args:
        counts (array or pd.Series or pd.DataFrame): Array of corrected and filtered neutron counts.
        N0 (float): Device-specific neutron calibration constant (neutron intensity above dry soil).
        Wlat (float): Gravimetric lattice water content in g of water per g of dry soil.
        Wsoc (float): Water equivalent of soil organic matter in g of water per g of dry soil.
        bulk_density (float): Soil dry bulk density in g cm-3.
        a0 (float): Parameter given in Desilets et al., 2010. Default is 0.0808.
        a1 (float): Parameter given in Desilets et al., 2010. Default is 0.372.
        a2 (float): Parameter given in Desilets et al., 2010. Default is 0.115.

    Returns:
        (array or pd.Series or pd.DataFrame): Volumetric water content in m3 m-3.

    References:
        Desilets, D., M. Zreda, and T.P.A. Ferré. 2010. Nature’s neutron probe:
        Land surface hydrology at an elusive scale with cosmic rays. Water Resour. Res. 46:W11505.
        doi.org/10.1029/2009WR008726

        Hawdon, A., D. McJannet, and J. Wallace (2014), Calibration and correction procedures for cosmic-ray neutron
        soil moisture probes located across Australia, Water Resour. Res., 50, 5029–5043, doi:10.1002/2013WR015138.
    """

    # Convert neutron counts into vwc
    vwc = (a0 / (counts / N0 - a1) - a2 - Wlat - Wsoc) * bulk_density
    return vwc


def sensing_depth(vwc, pressure, p_ref, bulk_density, Wlat, dist=None, method='Schron_2017'):
    """Function that computes the estimated sensing depth of the cosmic-ray neutron probe.
    The function offers several methods to compute the depth at which 86 % of the neutrons
    probe the soil profile.

    Args:
        vwc (array or pd.Series or pd.DataFrame): Estimated volumetric water content for each timestamp.
        pressure (array or pd.Series or pd.DataFrame): Atmospheric pressure in hPa for each timestamp.
        p_ref (float): Reference pressure in hPa.
        bulk_density (float): Soil bulk density.
        Wlat (float): Lattice water content.
        method (str): Method to compute the sensing depth. Options are 'Schron_2017' or 'Franz_2012'.
        dist (list or array): List of radial distances at which to estimate the sensing depth. Only used for the 'Schron_2017' method.

    Returns:
        (array or pd.Series or pd.DataFrame): Estimated sensing depth in cm.

    References:
        Franz, T.E., Zreda, M., Ferre, T.P.A., Rosolem, R., Zweck, C., Stillman, S., Zeng, X. and Shuttleworth, W.J., 2012.
        Measurement depth of the cosmic ray soil moisture probe affected by hydrogen from various sources.
        Water Resources Research, 48(8). doi.org/10.1029/2012WR011871

        Schrön, M., Köhli, M., Scheiffele, L., Iwema, J., Bogena, H. R., Lv, L., et al. (2017).
        Improving calibration and validation of cosmic-ray neutron sensors in the light of spatial sensitivity.
        Hydrol. Earth Syst. Sci. 21, 5009–5030. doi.org/10.5194/hess-21-5009-2017
    """

    # Determine sensing depth (D86)
    if method == 'Schron_2017':
        # See Appendix A of Schrön et al. (2017)
        Fp = 0.4922 / (0.86 - np.exp(-1 * pressure / p_ref))
        results = []
        for d in dist:
            # Compute r_star
            r_start = d / Fp

            # Compute soil depth that accounts for 86% of the neutron flux
            # theta is the total volumetric water equivalent, theta_sm + theta_lw (Schrön et al., 2017, Eq. 2),
            # so gravimetric lattice water is converted with the bulk density as in Franz et al. (2012, Eq. 3)
            D86 = 1 / bulk_density * (8.321 + 0.14249 * (0.96655 + np.exp(-0.01 * r_start)) * (20 + (bulk_density * Wlat + vwc)) / (
                    0.0429 + (bulk_density * Wlat + vwc)))

            results.append(D86)

    elif method == 'Franz_2012':
        results = 5.8 / (bulk_density * Wlat + vwc + 0.0829)
    else:
        raise ValueError('Method not recognized. Please select either "Schron_2017" or "Franz_2012".')
    return results


def abs_humidity(relative_humidity, temp):
    """
    Compute the absolute humidity (mass of water vapor per volume of air) in g m^-3 from relative humidity (%) and
    air temperature (Celsius). The saturation vapor pressure follows Eq. 3.8 of Campbell and Norman (1998).

    Args:
        relative_humidity (float): relative humidity (%)
        temp (float): temperature (Celsius)

    Returns:
        (float): Absolute humidity (g m^-3)

    References:
        Campbell, G. S., & Norman, J. M. (1998). An Introduction to Environmental Biophysics, 2nd ed. Springer, New York.
    """

    ### Atmospheric water vapor factor
    # Saturation vapor pressure
    e_sat = 0.611 * np.exp(17.502 * temp / (
            temp + 240.97)) * 1000  # in Pascals Eq. 3.8 p.41 Environmental Biophysics (Campbell and Norman)

    # Vapor pressure Pascals
    Pw = e_sat * relative_humidity / 100

    # Absolute humidity (g/m^3)
    C = 2.16679  # g K/J;
    abs_h = C * Pw / (temp + 273.15)
    return abs_h


def nrad_weight(h, theta, distances, depth, profiles=None, rhob=1.4, p=None, Hveg=0, tol=0.01):
    """Function to compute distance weights corresponding to each soil sample following the revised
    footprint weighting functions of Schrön et al. (2017).

    Args:
        h (float): Air humidity from 0.1 to 50 g/m^3. A single value for the calibration period is expected (e.g. mean humidity over the calibration hours). If an array is provided its mean is used.
        theta (np.array or pd.Series): Soil Moisture for each sample (0.02 - 0.50 m^3/m^3)
        distances (np.array or pd.Series): Distances from the location of each sample to the origin (0.5 - 600 m)
        depth (np.array or pd.Series): Depths for each sample (cm)
        profiles (np.array or pd.Series): Soil profiles ID for each sample. Required.
        rhob (np.array or pd.Series): Bulk density in g/cm^3
        p (float): Atmospheric pressure in hPa. A single value for the calibration period is expected. Required.
        Hveg (np.array or pd.Series): Vegetation height in m.
        tol (float): Tolerance for the iterative solution. Default is 0.01.

    Returns:
        theta_new (float): Weighted soil moisture values.
        weights (list): [theta_P, r_stars, Wrs] with the vertically averaged soil moisture, the scaled distance and the horizontal weight of each profile.

    References:
        Köhli, M., Schrön, M., Zreda, M., Schmidt, U., Dietrich, P., and Zacharias, S. (2015).
        Footprint characteristics revised for field-scale soil moisture monitoring with cosmic-ray
        neutrons. Water Resour. Res. 51, 5772–5790. doi:10.1002/2015WR017169

        Schrön, M., Köhli, M., Scheiffele, L., Iwema, J., Bogena, H. R., Lv, L.,
        Martini, E., Baroni, G., Rosolem, R., Weimar, J., Mai, J., Cuntz, M., Rebmann, C.,
        Oswald, S. E., Dietrich, P., Schmidt, U., and Zacharias, S.: Improving calibration and
        validation of cosmic-ray neutron sensors in the light of spatial sensitivity,
        Hydrol. Earth Syst. Sci., 21, 5009–5030, https://doi.org/10.5194/hess-21-5009-2017, 2017.
    """

    # Horizontal distance weights According to Eq. 6 and Table A1 in Schrön et al. (2017)
    # Method for calculating the horizontal distance weights from 0 to 1m
    def WrX(r, x, y):
        x00 = 3.7
        a00 = 8735;
        a01 = 22.689;
        a02 = 11720;
        a03 = 0.00978;
        a04 = 9306;
        a05 = 0.003632
        a10 = 2.7925e-2;
        a11 = 6.6577;
        a12 = 0.028544;
        a13 = 0.002455;
        a14 = 6.851e-5;
        a15 = 12.2755
        a20 = 247970;
        a21 = 23.289;
        a22 = 374655;
        a23 = 0.00191;
        a24 = 258552
        a30 = 5.4818e-2;
        a31 = 21.032;
        a32 = 0.6373;
        a33 = 0.0791;
        a34 = 5.425e-4

        x0 = x00
        A0 = (a00 * (1 + a03 * x) * np.exp(-a01 * y) + a02 * (1 + a05 * x) - a04 * y)
        A1 = ((-a10 + a14 * x) * np.exp(-a11 * y / (1 + a15 * y)) + a12) * (1 + x * a13)
        A2 = (a20 * (1 + a23 * x) * np.exp(-a21 * y) + a22 - a24 * y)
        A3 = a30 * np.exp(-a31 * y) + a32 - a33 * y + a34 * x

        return ((A0 * (np.exp(-A1 * r)) + A2 * np.exp(-A3 * r)) * (1 - np.exp(-x0 * r)))

    # Method for calculating the horizontal distance weights from 1 to 50m
    def WrA(r, x, y):
        a00 = 8735;
        a01 = 22.689;
        a02 = 11720;
        a03 = 0.00978;
        a04 = 9306;
        a05 = 0.003632
        a10 = 2.7925e-2;
        a11 = 6.6577;
        a12 = 0.028544;
        a13 = 0.002455;
        a14 = 6.851e-5;
        a15 = 12.2755
        a20 = 247970;
        a21 = 23.289;
        a22 = 374655;
        a23 = 0.00191;
        a24 = 258552
        a30 = 5.4818e-2;
        a31 = 21.032;
        a32 = 0.6373;
        a33 = 0.0791;
        a34 = 5.425e-4

        A0 = (a00 * (1 + a03 * x) * np.exp(-a01 * y) + a02 * (1 + a05 * x) - a04 * y)
        A1 = ((-a10 + a14 * x) * np.exp(-a11 * y / (1 + a15 * y)) + a12) * (1 + x * a13)
        A2 = (a20 * (1 + a23 * x) * np.exp(-a21 * y) + a22 - a24 * y)
        A3 = a30 * np.exp(-a31 * y) + a32 - a33 * y + a34 * x

        return A0 * np.exp(-A1 * r) + A2 * np.exp(-A3 * r)

    # Method for calculating the horizontal distance weights from 50 to 600m
    def WrB(r, x, y):
        b00 = 39006;
        b01 = 15002337;
        b02 = 2009.24;
        b03 = 0.01181;
        b04 = 3.146;
        b05 = 16.7417;
        b06 = 3727
        b10 = 6.031e-5;
        b11 = 98.5;
        b12 = 0.0013826
        b20 = 11747;
        b21 = 55.033;
        b22 = 4521;
        b23 = 0.01998;
        b24 = 0.00604;
        b25 = 3347.4;
        b26 = 0.00475
        b30 = 1.543e-2;
        b31 = 13.29;
        b32 = 1.807e-2;
        b33 = 0.0011;
        b34 = 8.81e-5;
        b35 = 0.0405;
        b36 = 26.74

        B0 = (b00 - b01 / (b02 * y + x - 0.13)) * (b03 - y) * np.exp(-b04 * y) - b05 * x * y + b06
        B1 = b10 * (x + b11) + b12 * y
        B2 = (b20 * (1 - b26 * x) * np.exp(-b21 * y * (1 - x * b24)) + b22 - b25 * y) * (2 + x * b23)
        B3 = ((-b30 + b34 * x) * np.exp(-b31 * y / (1 + b35 * x + b36 * y)) + b32) * (2 + x * b33)

        return B0 * np.exp(-B1 * r) + B2 * np.exp(-B3 * r)

    # Wrapper method for calculating the horizontal distance weights
    def Wr(r, x, y):
        if r <= 1:
            return WrX(r, x, y)
        elif r <= 50:
            return WrA(r, x, y)
        elif r <= 600:
            return WrB(r, x, y)
        else:
            raise ValueError("r must be between 0 and 600m")

    def rscaled(r, p, y, Hveg = 0):
        Fp = 0.4922 / (0.86 - np.exp(-p / 1013.25))
        Fveg = 1 - 0.17 * (1 - np.exp(-0.41 * Hveg)) * (1 + np.exp(-9.25 * y))
        return r / Fp / Fveg

    if profiles is None:
        raise ValueError("Profile ID's must be provided")

    if p is None:
        raise ValueError("Atmospheric pressure must be provided")

    # Rename variables to be consistent with the revised paper
    r = distances
    theta_ = np.mean(theta) # Start with the mean value of theta as initial guess
    bd = np.mean(rhob) # Neutrons are impacted by the bulk density across the whole area and not just the sample area. https://github.com/soilwater/crnpy/issues/9#issuecomment-2003813777

    # The Schrön et al. (2017) procedure (Sect. 3, step 4) evaluates the horizontal weighting
    # function with a single site air humidity for the calibration period, W_r(h, <theta>, p, Hveg).
    # A scalar (e.g. mean humidity over the calibration hours) is expected. If an array is
    # provided, its mean is used.
    h_ = float(np.mean(h))
    if np.ndim(h) > 0 and (np.nanmax(h) - np.nanmin(h)) > 1e-6:
        warnings.warn("A single air humidity value is expected for the calibration period. "
                      "The mean of the provided values was used.")

    # Vertical distance weights functions
    def D86(r, bd, y):
        return 1 / bd * (8.321 + 0.14249 * (0.96655 + np.exp(-0.01 * r)) * (20 + y) / (0.0429 + y))

    def Wd(d, r, bd, y):
        return np.exp(-2 * d / D86(r, bd, y))

    step = 0
    diff = 1
    while diff > tol:
        step += 1
        print(f"Step {step}, diff = {diff}", end="\r")
        # Calculate the scaled distance and D86
        r = rscaled(distances, p, theta_, Hveg)

        # Calculate the vertical average for each profile
        P = np.unique(profiles)
        theta_P = []
        r_stars = []
        for i in range(len(P)):
            profile = P[i]
            idx = profiles == profile
            depths_P = depth[idx]
            r_P = r[idx]
            theta_Pi = theta[idx]
            # Calculate the vertical distance weights. Penetration depth uses the field-average
            # soil moisture of the current iteration, D_P = D86(<theta>, r*_P) (Schrön et al., 2017, step 2)
            Wd_P = Wd(depths_P, r_P, bd, theta_)
            # Calculate the vertical average of theta
            theta_P_i = np.sum(Wd_P * theta_Pi) / np.sum(Wd_P)
            theta_P.append(theta_P_i)
            r_stars.append(np.mean(r_P))

        # Calculate the horizontal distance weights, W_r(r*, x=air humidity, y=field-average soil moisture)
        # following Schrön et al. (2017, step 4). Hveg enters only through the scaled distance r*.
        Wrs = np.array([Wr(r_star, h_, theta_) for r_star in r_stars])
        theta_new = np.sum(Wrs * theta_P) / np.sum(Wrs)
        diff = np.abs(theta_new - theta_)
        theta_ = theta_new

    print(f"Solution converged after {step} steps, the average soil moisture is {theta_new}")

    return theta_new, [theta_P, r_stars, Wrs]


def exp_filter(sm, T=1):
    """Exponential filter to estimate soil moisture in the rootzone from surface observations.

    Args:
        sm (list or array): Soil moisture in mm of water for the top layer of the soil profile.
        T (float): Characteristic time length in the same units as the measurement interval.

    Returns:
        sm_subsurface (list or array): Subsurface soil moisture in the same units as the input.

    References:
        Albergel, C., Rüdiger, C., Pellarin, T., Calvet, J.C., Fritz, N., Froissard, F., Suquia, D., Petitpa, A., Piguet, B. and Martin, E., 2008.
        From near-surface to root-zone soil moisture using an exponential filter: an assessment of the method based on in-situ observations and model
        simulations. Hydrology and Earth System Sciences, 12(6), pp.1323-1337.

        Franz, T.E., Wahbi, A., Zhang, J., Vreugdenhil, M., Heng, L., Dercon, G., Strauss, P., Brocca, L. and Wagner, W., 2020.
        Practical data products from cosmic-ray neutron sensing for hydrological applications. Frontiers in Water, 2, p.9.

        Rossini, P. and Patrignani, A., 2021. Predicting rootzone soil moisture from surface observations in cropland using an exponential filter.
        Soil Science Society of America Journal.
    """

    # Parameters
    sm = np.asarray(sm, dtype=float)
    sm_min = np.nanmin(sm)
    sm_max = np.nanmax(sm)
    ms = (sm - sm_min) / (sm_max - sm_min)

    # Pre-allocate soil water index array and recursive constant K
    SWI = np.ones_like(ms) * np.nan
    K = np.ones_like(ms) * np.nan

    # Initial conditions (Albergel et al., 2008): K = 1 and SWI equal to the first available observation
    valid = np.where(~np.isnan(ms))[0]
    if valid.size == 0:
        return SWI
    last = valid[0]
    SWI[last] = ms[last]
    K[last] = 1

    # Recursive formulation, Eqs. 4 and 6 in Albergel et al. (2008). The gain uses the time elapsed since the
    # last available observation, so missing values (NaN) do not interrupt the recursion.
    for n in range(last + 1, len(SWI)):
        if np.isnan(ms[n]):
            continue
        t_delta = n - last
        K[n] = K[last] / (K[last] + np.exp(-t_delta / T))
        SWI[n] = SWI[last] + K[n] * (ms[n] - SWI[last])
        last = n

    # Rootzone storage
    sm_subsurface = SWI * (sm_max - sm_min) + sm_min

    return sm_subsurface


def cutoff_rigidity(lat, lon):
    """Function to estimate the approximate cutoff rigidity for any point on Earth by interpolating the
    world grid of calculated vertical cutoff rigidities for epoch 1995.0 of Smart and Shea (2008), tabulated
    every 5 degrees in latitude and 15 degrees in longitude. The returned value can be used to select the appropriate
    neutron monitor station to estimate the cosmic-ray neutron intensity at the location of interest.

    Args:
        lat (float): Geographic latitude in decimal degrees. Value in range -90 to 90
        lon (float): Geographic longitude in decimal degrees. Values in range from 0 to 360.
            Typical negative longitudes in the west hemisphere will fall in the range 180 to 360.

    Returns:
        (float): Cutoff rigidity in GV. Compared with 102 neutron monitors the mean absolute error is about 0.2 GV,
            with the largest deviations (about 1.5 GV) in South America.

    Examples:
        Estimate the cutoff rigidity for Newark, NJ, US

        >>> zq = cutoff_rigidity(39.68, -75.75)
        >>> print(zq)
        2.25 GV (Value from NMDB is 2.40 GV)

    References:
        Smart, D. F., & Shea, M. A. (2008). World grid of calculated cosmic ray vertical cutoff rigidities
        for epoch 1995.0. Proceedings of the 30th International Cosmic Ray Conference (Mérida), 1, 733-736.
    """
    xq = lon
    yq = lat

    if xq < 0:
        xq = xq + 360  # west longitudes expressed on the 0-360 east longitude grid
    Z = np.array(data.cutoff_rigidity)
    x = np.linspace(0, 360, Z.shape[1])
    y = np.linspace(90, -90, Z.shape[0])
    X, Y = np.meshgrid(x, y)
    points = np.array((X.flatten(), Y.flatten())).T
    values = Z.flatten()
    zq = griddata(points, values, (xq, yq))

    return np.round(zq, 2)

def atmospheric_depth(elevation, latitude):
    """Function to estimate the atmospheric depth for any point on Earth according to McJannet and Desilets, 2023

    This function is required in the calculation of the location-dependent reference correction proposed by McJannet and Desilets, 2023.

    Args:
        elevation (float): Elevation in meters above sea level.
        latitude (float): Geographic latitude in decimal degrees. Value in range -90 to 90

    Returns:
        (float): Atmospheric depth in g/cm2

    References:
        Atmosphere, U. S. (1976). US standard atmosphere. National Oceanic and Atmospheric Administration.

        McJannet, D. L., & Desilets, D. (2023). Incoming Neutron Flux Corrections for Cosmic‐Ray Soil and Snow Sensors Using the Global Neutron Monitor Network. Water Resources Research, 59(4), e2022WR033889.
    """

    density_of_rock = 2670  # Density of rock in kg/m3
    air_pressure_sea_level = 1013.25  # Air pressure at sea level in hPa
    air_molar_mass = 0.0289644  # Air molar mass in kg/mol
    universal_gas_constant = 8.31432  # Universal gas constant in J/(mol*K)
    reference_temperature = 288.15  # Reference temperature Kelvin
    temperature_lapse_rate = -0.0065  # Temperature lapse rate in K/m

    # Gravity at sea-level calculation
    gravity_sea_level = 9.780327 * (
            1 + 0.0053024 * np.sin(np.radians(latitude)) ** 2 - 0.0000058 * np.sin(2 * np.radians(latitude)) ** 2)
    # Free air correction
    free_air = -3.086 * 10 ** -6 * elevation
    # Bouguer correction
    bouguer_corr = 4.193 * 10 ** -10 * density_of_rock * elevation
    # Total gravity
    gravity = gravity_sea_level + free_air + bouguer_corr

    # Air pressure calculation
    reference_air_pressure = air_pressure_sea_level * (
            1 + temperature_lapse_rate / reference_temperature * elevation) ** ((-gravity * air_molar_mass) / (
            universal_gas_constant * temperature_lapse_rate))

    # Atmospheric depth calculation
    atmospheric_depth = (10 * reference_air_pressure) / gravity
    return atmospheric_depth


def location_factor(site_atmospheric_depth, site_Rc, reference_atmospheric_depth, reference_Rc):
    """
    Function to estimate the location factor between two sites according to McJannet and Desilets, 2023.


    Args:
        site_atmospheric_depth (float): Atmospheric depth at the site in g/cm2. Can be estimated using the function `atmospheric_depth()`
        site_Rc (float): Cutoff rigidity at the site in GV. Can be estimated using the function `cutoff_rigidity()`
        reference_atmospheric_depth (float): Atmospheric depth at the reference location in g/cm2.
        reference_Rc (float): Cutoff rigidity at the reference location in GV.

    Returns:
        (float): Location-dependent correction factor.

    References:
        McJannet, D. L., & Desilets, D. (2023). Incoming Neutron Flux Corrections for Cosmic‐Ray Soil and Snow Sensors Using the Global Neutron Monitor Network. Water Resources Research, 59(4), e2022WR033889.

    """

    # Renamed variables based on the provided table
    c0 = -0.0009  # from C39
    c1 = 1.7699  # from C40
    c2 = 0.0064  # from C41
    c3 = 1.8855  # from C42
    c4 = 0.000013  # from C43
    c5 = -1.2237  # from C44
    epsilon = 1  # from C45

    # Translated formula with renamed variables from McJannet and Desilets, 2023
    tau_new = epsilon * (c0 * reference_atmospheric_depth + c1) * (
            1 - np.exp(
        -(c2 * reference_atmospheric_depth + c3) * reference_Rc ** (c4 * reference_atmospheric_depth + c5)))

    norm_factor = 1 / tau_new

    # Calculate the result using the provided parameters
    tau = epsilon * norm_factor * (c0 * site_atmospheric_depth + c1) * (
            1 - np.exp(-(c2 * site_atmospheric_depth + c3) * site_Rc ** (c4 * site_atmospheric_depth + c5)))
    return tau



def find_neutron_monitor(Rc, start_date=None, end_date=None, verbose=False):
    """Search for potential reference neutron monitoring stations based on cutoff rigidity.

    Args:
        Rc (float): Cutoff rigidity in GV. Values in range 1.0 to 3.0 GV.
        start_date (datetime): Start date for the period of interest.
        end_date (datetime): End date for the period of interest.
        verbose (bool): If True, print a expanded output of the incoming neutron flux data.

    Returns:
        (list): List of top five stations with closes cutoff rigidity.
            User needs to select station according to site altitude.

    Examples:
        >>> from crnpy import crnpy
        >>> Rc = 2.40 # 2.40 Newark, NJ, US
        >>> crnpy.find_neutron_monitor(Rc)
        Select a station with an altitude similar to that of your location. For more information go to: 'https://www.nmdb.eu/nest/help.php#helpstations

        Your cutoff rigidity is 2.4 GV
                STID                          NAME     R  Altitude_m
        40   NEWK                        Newark  2.40          50
        33   MOSC                        Moscow  2.43         200
        27   KIEL                          Kiel  2.36          54
        28  KIEL2                        KielRT  2.36          54
        31   MCRL  Mobile Cosmic Ray Laboratory  2.46         200
        32   MGDN                       Magadan  2.10         220
        42   NVBK                   Novosibirsk  2.91         163
        26   KGSN                      Kingston  1.88          65
        9    CLMX                        Climax  3.00        3400
        57   YKTK                       Yakutsk  1.65         105

    References:
        https://www.nmdb.eu/nest/help.php#helpstations
    """

    # Load file with list of neutron monitoring stations
    stations = pd.DataFrame(data.neutron_detectors, columns=["STID", "NAME", "R", "Altitude_m"])

    # Sort stations by closest cutoff rigidity
    idx_R = (stations['R'] - Rc).abs().argsort()

    if start_date is not None and end_date is not None:
        stations["Period available"] = False
        for i in range(10):
            station = stations.iloc[idx_R[i]]["STID"]
            try:
                if get_incoming_neutron_flux(start_date, end_date, station, verbose=verbose) is not None:
                    stations.iloc[idx_R[i], -1] = True
            except:
                pass

        if sum(stations["Period available"] == True) == 0:
            print("No stations available for the selected period!")
            result = stations.reindex(idx_R).head(10).rename_axis(None)
        else:
            stations = stations[stations["Period available"] == True]
            idx_R = (stations['R'] - Rc).abs().argsort()
            result = stations.iloc[idx_R.iloc[:10]]
    else:
        result = stations.reindex(idx_R).head(10).rename_axis(None)

    # Print results
    print('')
    print(
        """Select a station with an altitude similar to that of your location. For more information go to: 'https://www.nmdb.eu/nest/help.php#helpstations""")
    print('')
    print(f"Your cutoff rigidity is {Rc} GV")
    print(result)
    return result


def interpolate_incoming_flux(nmdb_timestamps, nmdb_counts, crnp_timestamps):
    """Function to interpolate incoming neutron flux to match the timestamps of the observations.

    Args:
        nmdb_timestamps (pd.Series or np.array): Series or array of timestamps in datetime format from the NMDB
        nmdb_counts (pd.Series or np.array): Series or array of incoming neutron flux counts from the NMDB
        crnp_timestamps (pd.Series or np.array): Series or array of timestamps in datetime format from the CRNP device

    Returns:
        (np.array): Incoming neutron flux matched to each CRNP timestamp. Same length as crnp_timestamps. Periods without
            NMDB data remain NaN; see the `fill_na` option of `correction_incoming_flux()`.
    """
    # Create a DataFrame from nmdb timestamps and counts
    df_nmdb = pd.DataFrame({'timestamp': nmdb_timestamps, 'counts': nmdb_counts})


    # Set the Timestamp column as the index
    df_nmdb.set_index('timestamp', inplace=True)

    # Reindex the DataFrame to the timestamps from the CRNP device using the nearest method
    # This will match each CRNP timestamp with the nearest NMDB timestamp
    interpolated_flux = df_nmdb.reindex(crnp_timestamps, method='nearest')['counts'].values

    return interpolated_flux


def lattice_water(clay_content, total_carbon=None):
    r"""Estimate the gravimetric lattice water content of the soil, i.e. the mass of water bound in the
    lattice of clay minerals per mass of dry soil, from pedotransfer functions.

    ![img1](img/lattice_water_simple.png) | ![img2](img/lattice_water_multiple.png)
    :-------------------------:|:-------------------------:
    $\omega_{lat} = 0.097 * clay(\%)$ | $\omega_{lat} = -0.028 + 0.077 * clay(\%) + 0.459 * carbon(\%)$
    Linear regression [gravimetric lattice water (% by mass) as a function of clay (% by mass)] done with data from Kansas State University - Soil Water Processes Lab. |  Multiple linear regression [gravimetric lattice water (% by mass) as a function of clay (% by mass) and soil carbon (% by mass)] done with data from Soil Water Processes Lab.

    All quantities are on a mass basis (gravimetric), not a volume basis.

    Args:
        clay_content (float): Clay content of the soil in percent by mass (g of clay per 100 g of dry soil).
        total_carbon (float, optional): Total carbon content of the soil in percent by mass.
            If None, the lattice water is estimated from the clay content only.

    Returns:
        (float): Gravimetric lattice water content in percent by mass (g of water per 100 g of dry soil).
            Divide by 100 to obtain the g/g fraction expected as `Wlat` by `counts_to_vwc()` and `sensing_depth()`.
            The clay-only regression has no intercept, so passing the clay content as a fraction returns the
            lattice water directly as a g/g fraction; the regression with carbon must be used with inputs in percent.
    """
    if total_carbon is None:
        lattice_water = 0.097 * clay_content
    else:
        lattice_water = -0.028 + 0.077 * clay_content + 0.459 * total_carbon
    return lattice_water


def latlon_to_utm(lat, lon, utm_zone_number=None, utm_zone_letter=None):
    """Convert geographic coordinates (lat, lon) to projected coordinates in the Universal Transverse Mercator (UTM) system.

    Function only applies to non-polar coordinates.
    If further functionality is required, consider using the utm module. See references for more information.

    ![UTM zones](https://upload.wikimedia.org/wikipedia/commons/thumb/b/b7/Universal_Transverse_Mercator_zones.svg/1920px-Universal_Transverse_Mercator_zones.svg.png)
    UTM zones on an equirectangular world map with irregular zones in red and New York City's zone highlighted. See [UTM zones](https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system#UTM_zones) for a full description.


    Args:
        lat (float, array): Latitude in decimal degrees.
        lon (float, array): Longitude in decimal degrees.
        utm_zone_number (int): UTM zone number. If None, the zone number is automatically calculated.
        utm_zone_letter (str): UTM zone letter. If None, the zone letter is automatically calculated.

    Returns:
        (float, float, int, str): Tuple of easting, northing, zone number and zone letter. First element is easting, second is northing, third is zone number and fourth is zone letter.

    References:
         Code adapted from utm module created by Tobias Bieniek (Github username: Turbo87)
         [https://github.com/Turbo87/utm](https://github.com/Turbo87/utm)

         [https://www.maptools.com/tutorials/grid_zone_details#](https://www.maptools.com/tutorials/grid_zone_details#)
    """
    # utm module requires numpy arrays
    if not isinstance(lat, np.ndarray):
        lat = np.array(lat)
    if not isinstance(lon, np.ndarray):
        lon = np.array(lon)

    if utm_zone_number is None or utm_zone_letter is None:
        easting, northing, zone_number, zone_letter = utm.from_latlon(lat, lon)
    else:
        easting, northing, zone_number, zone_letter = utm.from_latlon(lat, lon, utm_zone_number, utm_zone_letter)

    return easting, northing, zone_number, zone_letter


def euclidean_distance(px, py, x, y):
    """Function that computes the Euclidean distance between one point
    in space and one or more points.

    Args:
        px (float): x projected coordinate of the point.
        py (float): y projected coordinate of the point.
        x (list, ndarray, pandas.series): vector of x projected coordinates.
        y (list, ndarray, pandas.series): vector of y projected coordinates.

    Returns:
        (ndarray): Numpy array of distances from the point (px,py) to all the points in x and y vectors.
    """
    d = np.sqrt((px - x) ** 2 + (py - y) ** 2)
    return d


def spatial_average(x, y, z, buffer=100, min_neighbours=3, method='mean', rnd=False):
    """Moving buffer filter to smooth georeferenced two-dimensional data.

    Args:
        x (list or array): UTM x coordinates in meters.
        y (list or array): UTM y coordinates in meters.
        z (list or array): Values to be smoothed.
        buffer (float): Radial buffer distance in meters.
        min_neighbours (int): Minimum number of neighbours to consider for the smoothing.
        method (str): One of 'mean' or 'median'.
        rnd (bool): Boolean to round the final result. Useful in case of z representing neutron counts.

    Returns:
        (array): Smoothed version of z with the same dimension as z.
    """

    # Convert input data to Numpy arrays
    if (type(x) is not np.ndarray) or (type(y) is not np.ndarray):
        try:
            x = np.array(x)
            y = np.array(y)
        except Exception:
            raise ValueError("Input values cannot be converted to Numpy arrays.")

    if len(x) != len(y):
        raise ValueError(f"The number of x and y must be equal. Input x has {len(x)} values and y has {len(y)} values.")

    # Compute distances
    N = len(x)
    z_smooth = np.array([])
    for k in range(N):
        px = x[k]
        py = y[k]

        distances = euclidean_distance(px, py, x, y)
        idx_within_buffer = distances <= buffer

        if np.isnan(z[k]):
            z_new_val = np.nan
        elif len(distances[idx_within_buffer]) > min_neighbours:
            if method == 'mean':
                z_new_val = np.nanmean(z[idx_within_buffer])
            elif method == 'median':
                z_new_val = np.nanmedian(z[idx_within_buffer])
            else:
                raise ValueError(f"Method {method} does not exist. Provide either 'mean' or 'median'.")
        else:
            z_new_val = z[k]  # If there are not enough neighbours, keep the original value

        # Append smoothed value to array
        z_smooth = np.append(z_smooth, z_new_val)

    if rnd:
        z_smooth = np.round(z_smooth, 0)

    return z_smooth


def idw(x, y, z, X_pred, Y_pred, neighborhood=1000, p=1):
    """Function to interpolate data using inverse distance weight.

    Args:
        x (list or array): UTM x coordinates in meters.
        y (list or array): UTM y coordinates in meters.
        z (list or array): Values to be interpolated.
        X_pred (list or array): UTM x coordinates where z values need to be predicted.
        Y_pred (list or array): UTM y coordinates where z values need to be predicted.
        neighborhood (float): Only points within this radius in meters are considered for the interpolation.
        p (int): Exponent of the inverse distance weight formula. Typically, p=1 or p=2.

    Returns:
        (array): Interpolated values.

    References:
        [https://en.wikipedia.org/wiki/Inverse_distance_weighting](https://en.wikipedia.org/wiki/Inverse_distance_weighting)


    """

    # Flatten arrays to handle 1D and 2D arrays with the same code
    s = X_pred.shape  # Save shape
    X_pred = X_pred.flatten()
    Y_pred = Y_pred.flatten()

    # Pre-allocate output array
    Z_pred = np.full_like(X_pred, np.nan)

    for n in range(X_pred.size):
        # Distance between current and observed points
        d = euclidean_distance(X_pred[n], Y_pred[n], x, y)

        # Select points within neighborhood only for interpolation
        idx_neighbors = d < neighborhood

        # A prediction point that coincides with an observation takes the observed value (exact interpolation)
        if np.any(d == 0):
            Z_pred[n] = np.mean(z[d == 0])
        else:
            # Compute interpolated value at point of interest
            Z_pred[n] = np.sum(z[idx_neighbors] / d[idx_neighbors] ** p) / np.sum(1 / d[idx_neighbors] ** p)

    return np.reshape(Z_pred, s)


def interpolate_2d(x, y, z, dx=100, dy=100, method='cubic', neighborhood=1000):
    """Function for interpolating irregular spatial data into a regular grid.

    Args:
        x (list or array): UTM x coordinates in meters.
        y (list or array): UTM y coordinates in meters.
        z (list or array): Values to be interpolated.
        dx (float): Pixel width in meters.
        dy (float): Pixel height in meters.
        method (str): Interpolation method. One of 'cubic', 'linear', 'nearest', or 'idw'.
        neighborhood (float): Only points within this radius in meters are considered for the interpolation.

    Returns:
        x_pred (array): 2D array with x coordinates.
        y_pred (array): 2D array with y coordinates.
        z_pred (array): 2D array with interpolated values.

    References:
        [https://soilwater.github.io/pynotes-agriscience/notebooks/interpolation.html](https://soilwater.github.io/pynotes-agriscience/notebooks/interpolation.html)
    """

    # Drop NaN values in x y and z
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    idx_nan = np.isnan(x) | np.isnan(y) | np.isnan(z)

    if idx_nan.any():
        print(
            f"WARNING: {np.isnan(x).sum()}, {np.isnan(y).sum()}, and {np.isnan(z).sum()} NaN values were dropped from x, y, and z.")

    x = x[~idx_nan]
    y = y[~idx_nan]
    z = z[~idx_nan]

    # Create 2D grid for interpolation
    Nx = round((np.max(x) - np.min(x)) / dx) + 1
    Ny = round((np.max(y) - np.min(y)) / dy) + 1
    X_vec = np.linspace(np.min(x), np.max(x), Nx)
    Y_vec = np.linspace(np.min(y), np.max(y), Ny)
    X_pred, Y_pred = np.meshgrid(X_vec, Y_vec)

    if method in ['linear', 'nearest', 'cubic']:
        points = list(zip(x, y))
        Z_pred = griddata(points, z, (X_pred, Y_pred), method=method)

    elif method == 'idw':
        Z_pred = idw(x, y, z, X_pred, Y_pred, neighborhood)

    else:
        raise ValueError(f"Method {method} does not exist. Provide either 'cubic', 'linear', 'nearest', or 'idw'.")

    return X_pred, Y_pred, Z_pred


def rover_centered_coordinates(x, y):
    """Function to estimate the intermediate locations between two points, assuming the measurements were taken at a constant speed.

    Args:
        x (array): x coordinates.
        y (array): y coordinates.

    Returns:
        x_est (array): Estimated x coordinates.
        y_est (array): Estimated y coordinates.
    """

    # Make it datatype agnostic
    if (isinstance(x, pd.Series)):
        x = x.values
    if (isinstance(y, pd.Series)):
        y = y.values

    # Do the average of the two points
    x_est = (x[1:] + x[:-1]) / 2
    y_est = (y[1:] + y[:-1]) / 2

    # Add the first point to match the length of the original array
    x_est = np.insert(x_est, 0, x[0])
    y_est = np.insert(y_est, 0, y[0])

    return x_est, y_est


def uncertainty_counts(raw_counts, metric="std", fp=1, fw=1, fi=1):
    r"""Function to estimate the uncertainty of raw counts.

    Measurements of proportional neutron detector systems are governed by counting statistics that follow a Poissonian probability distribution (Zreda et al., 2012).
    The expected uncertainty in the neutron count rate $N$ is defined by the standard deviation $ \sqrt{N} $ (Jakobi et al., 2020).
    The CV% can be expressed as $ N^{-1/2} $

    Args:
        raw_counts (array): Raw neutron counts.
        metric (str): Either 'std' or 'cv' for standard deviation or coefficient of variation.
        fp (float): Pressure correction factor.
        fw (float): Humidity correction factor.
        fi (float): Incoming neutron flux correction factor.

    Returns:
        uncertainty (float): Uncertainty of raw counts.

    References:
        Jakobi J, Huisman JA, Schrön M, Fiedler J, Brogi C, Vereecken H and Bogena HR (2020) Error Estimation for Soil Moisture Measurements With
        Cosmic Ray Neutron Sensing and Implications for Rover Surveys. Front. Water 2:10. doi: 10.3389/frwa.2020.00010

        Zreda, M., Shuttleworth, W. J., Zeng, X., Zweck, C., Desilets, D., Franz, T., and Rosolem, R.: COSMOS: the COsmic-ray Soil Moisture Observing System,
        Hydrol. Earth Syst. Sci., 16, 4079–4099, https://doi.org/10.5194/hess-16-4079-2012, 2012.

    """

    s = fw / (fp * fi)
    if metric == "std":
        uncertainty = np.sqrt(raw_counts) * s
    elif metric == "cv":
        # sigma_N = s*sqrt(N) and N_cor = s*N (Jakobi et al., 2020, Eq. 7), so the coefficient of
        # variation sigma_N / N_cor = 1/sqrt(N) does not depend on the correction factors.
        uncertainty = 1 / np.sqrt(raw_counts)
    else:
        raise ValueError(f"Metric {metric} does not exist. Provide either 'std' or 'cv' for standard deviation or coefficient of variation.")
    return uncertainty


def uncertainty_vwc(raw_counts, N0, bulk_density, fp=1, fw=1, fi=1, a0=0.0808, a1=0.372, a2=0.115):
    r"""Function to estimate the uncertainty propagated to volumetric water content.

    The uncertainty of the volumetric water content is estimated by propagating the uncertainty of the raw counts.
    Following Eq. 10 in Jakobi et al. (2020), the uncertainty of the volumetric water content can be expressed as:
    $$
    \sigma_{\theta_g}(N) = \sigma_N \frac{a_0 N_0}{(N_{cor} - a_1 N_0)^4} \sqrt{(N_{cor} - a_1 N_0)^4 + 8 \sigma_N^2 (N_{cor} - a_1 N_0)^2 + 15 \sigma_N^4}
    $$

    Args:
        raw_counts (array): Raw neutron counts.
        N0 (float): Calibration parameter N0.
        bulk_density (float): Bulk density in g cm-3.
        fp (float): Pressure correction factor.
        fw (float): Humidity correction factor.
        fi (float): Incoming neutron flux correction factor.

    Returns:
        sigma_VWC (float): Uncertainty in terms of volumetric water content.

    References:
        Jakobi J, Huisman JA, Schrön M, Fiedler J, Brogi C, Vereecken H and Bogena HR (2020) Error Estimation for Soil Moisture Measurements With
        Cosmic Ray Neutron Sensing and Implications for Rover Surveys. Front. Water 2:10. doi: 10.3389/frwa.2020.00010
    """

    Ncorr = raw_counts * fw / (fp * fi)
    sigma_N = uncertainty_counts(raw_counts, metric="std", fp=fp, fw=fw, fi=fi)
    sigma_GWC = sigma_N * ((a0 * N0) / ((Ncorr - a1 * N0) ** 4)) * np.sqrt(
        (Ncorr - a1 * N0) ** 4 + 8 * sigma_N ** 2 * (Ncorr - a1 * N0) ** 2 + 15 * sigma_N ** 4)
    sigma_VWC = sigma_GWC * bulk_density

    return sigma_VWC
