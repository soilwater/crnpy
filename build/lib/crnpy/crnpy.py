# crnpy/crnpy.py
"""
`crnpy` is a Python package for processing cosmic ray neutron data.

 Created by Joaquin Peraza and Andres Patrignani.
"""

# Import modules
import sys
import warnings
import numpy as np
import pandas as pd
import requests
import io, datetime, os
from scipy.signal import savgol_filter
from scipy.interpolate import griddata
import crnpy.data as data

# Define python version
python_version = (3, 7)  # tuple of (major, minor) version requirement
python_version_str = str(python_version[0]) + "." + str(python_version[1])

# produce an error message if the python version is less than required
if sys.version_info < python_version:
    msg = "Module only runs on python version >= %s" % python_version_str
    raise Exception(msg)


def format_dates_df(df, col='timestamp', format='%Y-%m-%d %H:%M:%S', freq='H', round_time=True):
    """Helper function to change the format and round timestamps.

     Args:
         df (pandas.DataFrame): DataFrame with timestamp in the index.
         col (str, optional): Column with the timestamp. Default is 'timestamp'.
         format (str, optional): Format of the timestamp. Default is '%Y-%m-%d %H:%M:%S'.
         freq (str, optional): Rounding interval. 'H' for hourly, 'M' for minute, or None. Default is 'H'.
         round_time (bool, optional): Whether to round timestamps to the nearest frequency. Default is True.

     Returns:
         (pandas.DataFrame): DataFrame with formatted timestamps and rounded time.

     Examples:
            >>> from crnpy import crnpy
            >>> import pandas as pd
            >>> df = pd.DataFrame({'timestamp':['2020-01-01 00:00:00','2020-01-01 00:30:00','2020-01-01 01:00:00']})
            >>> df = crnpy.format_dates_df(df, col='timestamp', format='%Y-%m-%d %H:%M:%S', freq='H', round_time=True)
            >>> df
                            timestamp
            0   2020-01-01 00:00:00
            1   2020-01-01 00:00:00
            2   2020-01-01 01:00:00
     """

    # Change format of timestamp if needed
    if df[col].dtype != 'datetime64[ns]':
        df[col] = pd.to_datetime(df[col], format=format)

    # Round timestamps to nearest frequency
    if round_time:
        df[col] = df[col].dt.round(freq)

    # Fill in rows with missing timestamps
    start_date = df[col].iloc[0]
    end_date = df[col].iloc[-1]
    date_range = pd.date_range(start_date, end_date, freq=freq)
    for date in date_range:
        if date not in df[col].values:
            print('Adding missing date:',date)
            new_line = pd.DataFrame({col:date}, index=[-1]) # By default fills columns with np.nan
            source = pd.concat([df,new_line])

    df.sort_values(by=col, inplace=True)
    df.reset_index(drop=True, inplace=True)
    df.set_index(col, inplace=True)
    return df

def count_time(counts=None, timestamp_col=None):
    """Approximate counting time.

    Args:
        counts (pandas.DataFrame): DataFrame with neutron counts, might have DateTimeIndex.
        timestamp_col (pandas.Series): Series with timestamps. If counts has a DateTimeIndex, timestamp_col is not needed.

    Returns:
        (pandas.Series): Series with the approximate counting time for each observation.

    Examples:
        Using `count_time` in a console environment:

        >>> df = pd.DataFrame(...)
        >>> count_time(timestamp_col=df['timestamp'])
        0   3600.0
        1   3600.0
        2   3600.0
    """
    if timestamp_col is not None:
        if not isinstance(timestamp_col, pd.Series):
            raise TypeError('timestamp_col must be a pandas Series.')
        if timestamp_col.dtype != 'datetime64[ns]':
            raise TypeError('timestamp_col must be a pandas Series with datetime64[ns] dtype.')

        count_time = timestamp_col.diff().dt.total_seconds()
        return count_time

    if counts is not None:
        if not isinstance(counts, pd.DataFrame):
            raise TypeError('counts must be a pandas DataFrame.')

        if not counts.index.dtype == 'datetime64[ns]':
            raise TypeError('counts must have a DateTimeIndex.')

        count_time = counts.index.to_series().diff().dt.total_seconds()
        return count_time

    raise TypeError('Either counts or timestamp_col must be provided.')

def fill_counts(counts, count_times=None, timestamp_col=None, expected_time=False, threshold=0.25, limit=3):
    """Fill missing neutron counts. Observation periods shorter than threshold are discarded (replaced with NaN).
    
    Args:
        counts (pandas.DataFrame): DataFrame with neutron counts, might have DateTimeIndex.
        count_time (pandas.Series or pandas.DataFrame): Counting time in seconds. If a DataFrame is provided, it must have the same number of columns as `counts`.
        timestamp_col (pandas.Series): Series with timestamps. If counts has a DateTimeIndex, timestamp_col is not needed.
        expected_time (int): Expected counting time in seconds. If not provided, it is calculated as the median of the counting times.
        threshold (float): Minimum fraction of the neutron integration time. Default is 0.25.

    Returns:
        (pandas.DataFrame): DataFrame with linearly interpolated neutron counts.

    Examples:
        Using `fill_counts` in a console environment:

        >>> counts = pd.DataFrame({'counts':[100,105,98,102], count_time:[3600,200,3600,3600]})
        >>> fill_counts(counts, count_time=count_time, expected_time=3600, threshold=0.25)
        0   100.0
        1   NaN
        2   98.0
        3   102.0
    """

    counts=counts.copy()

    if type(counts.index) == pd.core.indexes.datetimes.DatetimeIndex and isinstance(count_times, type(None)):
        print("No count time columns provided. Using timestamp index to compute count time.")
        count_times = count_time(timestamp_col=counts.index.to_series())

    elif not isinstance(timestamp_col, type(None)) and isinstance(count_times, type(None)):
        if timestamp_col.dtype != 'datetime64[ns]':
            if len(timestamp_col) != len(counts):
                raise ValueError('Timestamp column length does not match number of readings.')
            print("No count time columns provided. Using timestamp column to compute count time.")
            count_times = count_time(timestamp_col=timestamp_col)
        else:
            raise TypeError('Timestamp column must be a pandas Series with datetime64[ns] dtype.')

    if type(counts.index) != pd.core.indexes.datetimes.DatetimeIndex and isinstance(count_times, type(None)) and isinstance(timestamp_col, type(None)):
        raise ValueError('Count_times must be provided, or timestamp_col must be provided, or counts must have a DatetimeIndex.')

    if len(counts) != len(count_times):
        raise ValueError('Count times length does not match number of readings.')

    if expected_time is False:
        expected_time = count_times.median()
        print('Using median count time as expected count time:', expected_time)

    # Replace values below threshold with NaN
    time_threshold = round(expected_time * threshold)

    if type(count_times) == pd.core.frame.DataFrame:
        if len(count_times.columns) == 1:
            idx_nan = count_times[count_times < time_threshold].index
            counts.loc[idx_nan] = np.nan
        else:
            for i in range(len(count_times.columns)):
                idx_nan = count_times[count_times.iloc[:,i] < time_threshold].index
                counts.iloc[:,i].loc[idx_nan] = np.nan
    elif type(count_times) == pd.core.series.Series:
        idx_nan = count_times[count_times < time_threshold].index
        counts.loc[idx_nan] = np.nan
    elif type(count_times) == np.ndarray:
        idx_nan = np.where(count_times < time_threshold)
        counts.loc[idx_nan] = np.nan

    # Fill missing values with linear interpolation and round to nearest integer
    counts = counts.interpolate(method='linear', limit=limit, limit_direction='both').round()
    return counts

def normalize_counts(counts, count_time=3600, count_times=None, timestamp_col=None):
    """Normalize neutron counts to the desired counting time.
    
    Args:
        counts (pandas.DataFrame): Dataframe containing only the columns with neutron counts.
        count_time (int): Count time in seconds for normalization. Default is 3600 seconds.
        count_times (pandas.Series or pandas.DataFrame): Counting time in seconds. If a DataFrame is provided, it must have the same number of columns as counts.
        timestamp_col (pandas.Series): Timestamp column, used to calculate count time if count_times is not provided, it must have the same number of rows as counts.
        
    Returns:
        (pandas.DataFrame): Normalized neutron counts.

    """

    if count_times is None and type(counts.index) == pd.core.indexes.datetimes.DatetimeIndex:
        print("No count_times columns provided. Using timestamp index to compute count time.")
        count_times = counts.index.to_series().diff().dt.total_seconds()

    elif count_times is None and not isinstance(timestamp_col, type(None)):
        if len(timestamp_col) != len(counts):
            raise ValueError('Timestamp column length does not match number of readings.')
        print("No count_times columns provided. Using timestamp column to compute count time.")
        if timestamp_col.dtype != 'datetime64[ns]':
            raise TypeError('Timestamp column must be a pandas Series with datetime64[ns] dtype.')
        count_times = count_time(timestamp_col=timestamp_col)


    if isinstance(count_times, type(None)):
        raise ValueError('Count time must be provided, or `timestamp_col` must be provided, or counts must have a DatetimeIndex.')

    if len(counts) != len(count_times):
        raise ValueError('Count times length does not match number of readings.')

    #Normalize counts rounded to integer
    if type(count_times) == pd.core.series.Series or len(count_times.columns) == 1:
        normalized_counts = counts.div(count_times, axis=0).mul(count_time).round()
        return normalized_counts
    else:
        normalized_counts = counts.copy()
        count_times = count_times.copy()
        for i in range(len(count_times.columns)):
            normalized_counts[normalized_counts.columns[i]] = normalized_counts.iloc[:,i].div(count_times.iloc[:,i], axis=0).mul(count_time).round()
        return normalized_counts


def compute_total_raw_counts(counts, nan_strategy=None, timestamp_col=None):
    """Compute the sum of uncorrected neutron counts for all detectors.

    Args:
        counts (pandas.DataFrame): Dataframe containing only the columns with neutron counts.
        nan_strategy (str): Strategy to use for NaN values. Options are 'interpolate', 'average', or None. Default is None.

    Returns:
        (pandas.DataFrame): Dataframe with the sum of uncorrected neutron counts for all detectors.
    """

    counts=counts.copy()

    if counts.isnull().values.any():
        if nan_strategy is None:
            raise ValueError('NaN values found. Please fill missing values or provide a strategy. See documentation for more information.')
        elif nan_strategy == 'interpolate':
            print('NaN values found. Interpolating missing values using fill_counts().')
            if not isinstance(timestamp_col, type(None)):
                if type(timestamp_col) != pd.core.series.Series:
                    if timestamp_col.dtype != 'datetime64[ns]':
                        raise TypeError('Timestamp column must be a pandas Series with datetime64[ns] dtype.')
                counts = fill_counts(counts, timestamp_col)
            elif type(counts.index) != pd.core.indexes.datetimes.DatetimeIndex:
                raise ValueError('`timestamp_col` must be provided if `counts` does not have a DatetimeIndex.')
            counts = fill_counts(counts)
        elif nan_strategy == 'average':
            if len(counts.columns) == 1:
                raise ValueError('Only one detector found. Cannot use average strategy.')
            print('NaN values found. Replacing missing values with average of other detectors before summing.')
            counts = counts.apply(lambda x: x.fillna(counts.mean(axis=1)),axis=0)
        else:
            raise ValueError('Invalid strategy.')

    #Compute sum of counts
    total_raw_counts = counts.sum(axis=1)
    # Replace zeros with NaN
    total_raw_counts = total_raw_counts.replace(0, np.nan)
    return total_raw_counts


def drop_outliers(raw_counts, window=5, store_outliers=False, min_counts=None, max_counts=None):
    """Computation of a moving modified Z-score based on the median absolute difference.
    
    Args:
        raw_counts (pandas.DataFrame): Dataframe containing only the columns with neutron counts.
        window (int): Window size for the moving median. Default is 11.
        store_outliers (bool): If True, store the outliers in a new column. Default is False.
        min_counts (int): Minimum number of counts for a reading to be considered valid. Default is None.
        max_counts (int): Maximum number of counts for a reading to be considered valid. Default is None.

    Returns:
        (pandas.DataFrame): Dataframe without outliers.
        or
        (pandas.DataFrame, pandas.DataFrame): Dataframe without outliers and dataframe with outliers.

    References:
        Iglewicz, B. and Hoaglin, D.C., 1993. How to detect and handle outliers (Vol. 16). Asq Press.
    """

    if min_counts is not None:
        lower_count = np.sum(raw_counts < min_counts)
        if lower_count > len(raw_counts) * 0.25:
            print(f"WARNING: Discarded {lower_count} counts below {min_counts}. This is more than 25% of the total number of readings. Consider increasing the minimum counts threshold.")
        else:
            print(f"Discarded counts below {min_counts}: {lower_count}")
        raw_counts = raw_counts[raw_counts >= min_counts]
    if max_counts is not None:
        upper_count = np.sum(raw_counts > max_counts)
        print(f"Discarded counts above {max_counts}: {upper_count}")
        raw_counts = raw_counts[raw_counts <= max_counts]

    # Compute median absolute difference
    median = raw_counts.rolling(window, center=True).median()
    diff = np.abs(raw_counts - median)
    mad = diff.rolling(window, center=True).median()

    # Compute modified Z-score
    modified_z_score = 0.6745 * diff / mad
    outliers = raw_counts[modified_z_score > 3.5]
    # Drop outliers
    raw_counts = raw_counts[modified_z_score < 3.5]

    if store_outliers:
        return raw_counts, outliers
    print(f"Discarded {len(outliers)} outliers using modified Z-score.")
    return raw_counts


def fill_missing_atm(cols_atm, limit=24):
    """Fill missing values in atmospheric variables. Gap filling is performed using a
    piecewise cubic Hermite interpolating polynomial (pchip method) that is restricted to intervals
    of missing data with a limited number of values and surrounded by valid observations.
    There is no interpolation at the end of the time series.

    Args:
        col_atm (pandas.Series or pandas.DataFrame): Atmospheric variables to fill.
        limit (int): Maximum number of consecutive missing values to interpolate. Default is 24.

    Returns:
        (pandas.DataFrame): Atmospheric variables with filled missing values using a piecewise cubic Hermite polynomial.

    References:
        https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.interpolate.html
    """

    # Fill missing values in atmospheric variables
    return cols_atm.interpolate(method='pchip', limit=limit, limit_direction='both')

def atm_correction(raw_counts, pressure, humidity, temp, Pref, Aref, L, incoming_neutrons=None, incoming_Ref=None):
    r"""Correct neutron counts for atmospheric factors and incoming neutron flux.

    This function corrects neutron counts for atmospheric pressure, and absolute humidity using the method described in Zreda et al. (2012) and Anderson et al. (2017). The correction is performed using the following equation:

    $$
    C_{corrected} = \frac{C_{raw} \cdot f_w}{f_p \cdot f_i}
    $$

    - fp: pressure correction factor
    - fw: abosolute humidity correction factor
    - fi: incoming neutron flux correction factor
    
    Args:
        raw_counts (list or array): Neutron counts to correct.
        pressure (list or array): Atmospheric pressure readings.
        humidity (list or array): Atmospheric humidity readings in %.
        temp (list or array): Atmospheric temperature readings in Celsius.
        Pref (float): Reference atmospheric pressure in millibars (mbar).
        Aref (float): Reference absolute humidity.
        L (float): Atmospheric attenuation coefficient.
        incoming_neutrons (list or array): Incoming neutron flux. Default is None.
        incoming_Ref (float): Reference incoming neutron flux. Default is None.
        
    Returns:
        (numpy.array): Total neutron counts corrected by atmospheric conditions.
        
    References:
        Zreda, M., Shuttleworth, W. J., Zeng, X., Zweck, C., Desilets, D., Franz, T., et al. (2012).
        COSMOS: the cosmic-ray soil moisture observing system. Hydrol. Earth Syst. Sci. 16, 4079–4099.
        doi: 10.5194/hess-16-4079-2012

        Andreasen, M., Jensen, K.H., Desilets, D., Franz, T.E., Zreda, M., Bogena, H.R. and Looms, M.C., 2017.
        Status and perspectives on the cosmic‐ray neutron method for soil moisture estimation and other
        environmental science applications. Vadose zone journal, 16(8), pp.1-11. doi.org/10.2136/vzj2017.04.0086
    """

    ### Barometric pressure factor
    fp = np.exp((Pref - pressure) / L) # Zreda et al. 2017 Eq 5.

    ### Atmospheric water vapor factor
    # Saturation vapor pressure
    e_sat = 0.611 * np.exp(17.502 * temp / (temp + 240.97)) * 1000 # in Pascals Eq. 3.8 p.41 Environmental Biophysics (Campbell and Norman)

    # Vapor pressure Pascals
    Pw = e_sat * humidity/100

    # Absolute humidity (g/m^3)
    C = 2.16679 # g K/J;
    A = C * Pw / (temp + 273.15)
    fw = 1 + 0.0054*(A - Aref) # Zreda et al. 2017 Eq 6.

    ### Incoming neutron flux factor
    if incoming_neutrons is None:
        fi = 1
        warnings.warn("Ignoring incoming neutron flux correction factor (using value fi=1)")
    else:
        if incoming_Ref is None and not isinstance(incoming_neutrons, type(None)):
            incoming_Ref = incoming_neutrons[0]
            warnings.warn('Reference incoming neutron flux not provided. Using first value of incoming neutron flux.')

        fi = incoming_neutrons/incoming_Ref
        fi.fillna(1.0, inplace=True) # Use a value of 1 for days without data

    # Apply correction factors
    return np.round((raw_counts*fw)/(fp*fi))



def get_incoming_neutron_flux(start_date, end_date, station, utc_offset=0, expand_window = 0,  verbose=False):
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

    References:
        Documentation available:https://www.nmdb.eu/nest/help.php#howto
    """

    # Example: get_incoming_flux(station='IRKT',start_date='2020-04-10 11:00:00',end_date='2020-06-18 17:00:00')
    # Template url = 'http://nest.nmdb.eu/draw_graph.php?formchk=1&stations[]=KERG&output=ascii&tabchoice=revori&dtype=corr_for_efficiency&date_choice=bydate&start_year=2009&start_month=09&start_day=01&start_hour=00&start_min=00&end_year=2009&end_month=09&end_day=05&end_hour=23&end_min=59&yunits=0'


    # Expand the time window by 1 hour to ensure an extra observation is included in the request.
    start_date -= pd.Timedelta(hours=expand_window)
    end_date += pd.Timedelta(hours=expand_window)


    # Convert local time to UTC
    start_date = start_date - datetime.timedelta(hours=utc_offset)
    end_date = end_date - datetime.timedelta(hours=utc_offset)
    date_format = '%Y-%m-%d %H:%M:%S'
    root = 'http://www.nmdb.eu/nest/draw_graph.php?'
    url_par = [ 'formchk=1',
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

    if verbose > 0:
        print(f"Retrieving data from {url}")

    r = requests.get(url).content.decode('utf-8')

    # Subtract 1 hour to restore the last date included in the request.
    end_date -= pd.Timedelta('1H')
    start = r.find("RCORR_E\n") + 8
    end = r.find('\n</code></pre><br>Total') - 1
    s = r[start:end]
    s2 = ''.join([row.replace(';',',') for row in s])
    try:
        df_flux = pd.read_csv(io.StringIO(s2), names=['timestamp','counts'])
    except:
        if verbose > -1:
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
    #print(acknowledgement)

    return df_flux.set_index('timestamp')

def interpolate_incoming_flux(df_flux, timestamps):
    """Function to interpolate incoming neutron flux to match the timestamps of the observations.

    Args:
        df_flux (pd.DataFrame): Dataframe returned by get_incoming_flux method.
        timestamps (pd.series or pd.DataFrame or pd.DatetimeIndex): Timestamps to interpolate the incoming neutron flux.

    Returns:
        (pd.DataFrame): Dataframe containing interpolated incoming neutron flux.
    """

    # Check that index is datetime
    if not isinstance(df_flux.index, pd.DatetimeIndex):
        raise ValueError('Index of df_flux must be datetime')

    for timestamp in timestamps:
        if timestamp not in df_flux.index:
            df_flux.loc[timestamp] = np.nan
    df_flux = df_flux.sort_index()

    # Interpolate nan values
    df_flux = df_flux['counts'].interpolate(method='nearest', limit_direction='both')

    # Retur only the values for the selected timestamps
    return df_flux.loc[timestamps]


def smooth_1d(corrected_counts, window=5, order=3, method='moving_median'):
    """Use a Savitzky-Golay filter to smooth the signal of corrected neutron counts or another one-dimensional array (e.g. computed volumetric water content).

    Args:
        corrected_counts (pd.DataFrame): Dataframe containing the corrected neutron counts.
        window (int): Window size for the Savitzky-Golay filter. Default is 5.
        method (str): Method to use for smoothing the data. Default is 'moving_median'.
            Options are 'moving_average', 'moving_median' and 'savitzky_golay'.
        order (int): Order of the Savitzky-Golay filter. Default is 3.

    Returns:
        (pd.DataFrame): DataFrame with smoothed neutron counts.

    References:
        Franz, T.E., Wahbi, A., Zhang, J., Vreugdenhil, M., Heng, L., Dercon, G., Strauss, P., Brocca, L. and Wagner, W., 2020.
        Practical data products from cosmic-ray neutron sensing for hydrological applications. Frontiers in Water, 2, p.9.
        doi.org/10.3389/frwa.2020.00009
    """

    if method == 'moving_average':
        corrected_counts = corrected_counts.rolling(window=window, center=True, min_periods=1).mean()
    elif method == 'moving_median':
        corrected_counts = corrected_counts.rolling(window=window, center=True, min_periods=1).median()

    elif method == 'savitzky_golay':
        if corrected_counts.isna().any():
            print('Dataframe contains NaN values. Please remove NaN values before smoothing the data.')

        if type(corrected_counts) == pd.core.series.Series:
            filtered = np.round(savgol_filter(corrected_counts,window,order))
            corrected_counts = pd.DataFrame(filtered,columns=['counts'], index=corrected_counts.index)
        elif type(corrected_counts) == pd.core.frame.DataFrame:
            for col in corrected_counts.columns:
                corrected_counts[col] = np.round(savgol_filter(corrected_counts[col],window,order))
    else:
        raise ValueError('Invalid method. Please select a valid filtering method., options are: moving_average, moving_median, savitzky_golay')
    corrected_counts = corrected_counts.ffill(limit=window).bfill(limit=window).copy()
    return corrected_counts


def bwe_correction(counts, bwe, r2_N0=0.05):
    """Function to correct for biomass effects in neutron counts.
    following the approach described in Baatz et al., 2015.

    Args:
        counts (array or pd.Series or pd.DataFrame): Array of ephithermal neutron counts.
        bwe (float): Biomass water equivalent kg m-2.
        r2_N0 (float): Ratio of the neutron counts reduction (counts kg-1) to the neutron calibration constant (N0). Default is 0.05 (Baatz et al., 2015).

    Returns:
        (array or pd.Series or pd.DataFrame): Array of corrected neutron counts for biomass effects.

    References:
        Baatz, R., H. R. Bogena, H.-J. Hendricks Franssen, J. A. Huisman, C. Montzka, and H. Vereecken (2015),
        An empiricalvegetation correction for soil water content quantification using cosmic ray probes,
        Water Resour. Res., 51, 2030–2046, doi:10.1002/ 2014WR016443.
    """

    return counts/(1 - bwe*r2_N0)

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


def road_correction(counts, theta_N, road_width, road_distance=0.0, theta_road=0.12, p0=0.42, p1=0.5, p2=1.06, p3=4, p4=0.16, p6=0.94, p7=1.10, p8=2.70, p9=0.01):
    """Function to correct for road effects in neutron counts.
    following the approach described in Schrön et al., 2018.

    Args:
        counts (array or pd.Series or pd.DataFrame): Array of ephithermal neutron counts.
        theta_N (float): Volumetric water content of the soil estimated from the uncorrected neutron counts.
        road_width (float): Width of the road in m.
        road_distance (float): Distance of the road from the sensor in m. Default is 0.0.
        theta_road (float): Volumetric water content of the road. Default is 0.12.
        p0-p9 (float): Parameters of the correction function. Default values are from Schrön et al., 2018.

    Returns:
        (array or pd.Series or pd.DataFrame): Array of corrected neutron counts for road effects.

    References:
        Schrön,M.,Rosolem,R.,Köhli,M., Piussi,L.,Schröter,I.,Iwema,J.,etal. (2018).Cosmic-ray neutron rover surveys
        of field soil moisture and the influence of roads.WaterResources Research,54,6441–6459.
        https://doi. org/10.1029/2017WR021719
    """
    F1 = p0 * (1-np.exp(-p1*road_width))
    F2 = -p2 - p3 * theta_road - ((p4 + theta_road) / (theta_N))
    F3 = p6 * np.exp(-p7 * (road_width ** -p8) * road_distance ** 4) + (1 - p6) * np.exp(-p9 * road_distance)

    C_roads = 1 + F1 * F2 * F3

    corrected_counts = counts / C_roads

    return corrected_counts

def counts_to_vwc(counts, N0, Wlat, Wsoc ,bulk_density, a0=0.0808,a1=0.372,a2=0.115):
    r"""Function to convert corrected and filtered neutron counts into volumetric water content.

    This method implements soil moisture estimation using the non-linear relationship between neutron count and soil volumetric water content following the approach described in Desilets et al., 2010.

    $\theta(N) =\frac{a_0}{(\frac{N}{N_0}) - a_1} - a_2 $

    Args:
        counts (array or pd.Series or pd.DataFrame): Array of corrected and filtered neutron counts.
        N0 (float): Device-specific neutron calibration constant.
        Wlat (float): Lattice water content.
        Wsoc (float): Soil organic carbon content.
        bulk_density (float): Soil bulk density.
        a0 (float): Parameter given in Zreda et al., 2012. Default is 0.0808.
        a1 (float): Parameter given in Zreda et al., 2012. Default is 0.372.
        a2 (float): Parameter given in Zreda et al., 2012. Default is 0.115.
        
    Returns:
        (array or pd.Series or pd.DataFrame): Volumetric water content in m3 m-3.
        
    References:
        Desilets, D., M. Zreda, and T.P.A. Ferré. 2010. Nature’s neutron probe:
        Land surface hydrology at an elusive scale with cosmic rays. Water Resour. Res. 46:W11505.
        doi.org/10.1029/2009WR008726
    """

    # Convert neutron counts into vwc
    vwc = (a0 / (counts/N0-a1) - a2 - Wlat - Wsoc) * bulk_density
    return vwc



def sensing_depth(vwc, pressure, p_ref, bulk_density, Wlat, dist, method='Schron_2017'):
    # Convert docstring to google format
    """Function that computes the estimated sensing depth of the cosmic-ray neutron probe.
    The function offers several methods to compute the depth at which 86 % of the neutrons
    probes the soil profile.

    Args:
        vwc (array or pd.Series or pd.DataFrame): Estimated volumetric water content for each timestamp.
        pressure (array or pd.Series or pd.DataFrame): Atmospheric pressure in hPa for each timestamp.
        p_ref (float): Reference pressure in hPa.
        bulk_density (float): Soil bulk density.
        Wlat (float): Lattice water content.
        method (str): Method to compute the sensing depth. Options are 'Schron_2017' or 'Franz_2012'.
        dist (list or array): List of radial distances at which to estimate the sensing depth. Only used for the 'Schron_2017' method.

    Returns:
        (array or pd.Series or pd.DataFrame): Estimated sensing depth in m.

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
        Fp = 0.4922 / (0.86 - np.exp(-1 * pressure / p_ref));
        Fveg = 0
        results = []
        for d in dist:
            # Compute r_star
            r_start = d/Fp

            # Compute soil depth that accounts for 86% of the neutron flux
            D86 = 1/ bulk_density * (8.321+0.14249*(0.96655 + np.exp(-0.01*r_start))*(20+(Wlat+vwc)) / (0.0429+(Wlat+vwc)))
            results.append(D86)

    elif method == 'Franz_2012':
        results = [5.8/(bulk_density*Wlat+vwc+0.0829)]

    return results

def estimate_abs_humidity(RH, temp):
    """
    Compute the actual vapor pressure (e) in g m^-3 using RH (%) and current temperature (c) observations.

    Args:
        RH (float): relative humidity (%)
        temp (float): temperature (Celsius)

    Returns:
        float: actual vapor pressure (g m^-3)
    """

    ### Atmospheric water vapor factor
    # Saturation vapor pressure
    e_sat = 0.611 * np.exp(17.502 * temp / (
                temp + 240.97)) * 1000  # in Pascals Eq. 3.8 p.41 Environmental Biophysics (Campbell and Norman)

    # Vapor pressure Pascals
    Pw = e_sat * RH / 100

    # Absolute humidity (g/m^3)
    C = 2.16679  # g K/J;
    abs_h = C * Pw / (temp + 273.15)
    return abs_h


def nrad_weight(h,theta,distances,depth,rhob=1.4):
    """Function to compute distance weights corresponding to each soil sample.

    Args:
        h (float): Air Humidity  from 0.1  to 50 in g/m^3. When h=0, the function will skip the distance weighting.
        theta (array or pd.Series or pd.DataFrame): Soil Moisture for each sample (0.02 - 0.50 m^3/m^3)
        distances (array or pd.Series or pd.DataFrame): Distances from the location of each sample to the origin (0.5 - 600 m)
        depth (array or pd.Series or pd.DataFrame): Depths for each sample (m)
        rhob (float): Bulk density in g/cm^3

    Returns:
        (array or pd.Series or pd.DataFrame): Distance weights for each sample.

    References:
        Köhli, M., Schrön, M., Zreda, M., Schmidt, U., Dietrich, P., and Zacharias, S. (2015).
        Footprint characteristics revised for field-scale soil moisture monitoring with cosmic-ray
        neutrons. Water Resour. Res. 51, 5772–5790. doi:10.1002/2015WR017169
    """

    # Table A1. Parameters for Fi and D86
    p10 = 8735;       p11 = 17.1758; p12 = 11720;      p13 = 0.00978;   p14 = 7045;      p15 = 0.003632;
    p20 = 2.7925e-2;  p21 = 5.0399;  p22 = 2.8544e-2;  p23 = 0.002455;  p24 = 6.851e-5;  p25 = 9.2926;
    p30 = 247970;     p31 = 17.63;   p32 = 374655;     p33 = 0.00191;   p34 = 195725;
    p40 = 5.4818e-2;  p41 = 15.921;  p42 = 0.6373;     p43 = 5.99e-2;   p44 = 5.425e-4;
    p50 = 1383702;    p51 = 4.156;   p52 = 5325;       p53 = 0.00238;   p54 = 0.0156;    p55 = 0.130;     p56 = 1521;
    p60 = 6.031e-5;   p61 = 98.5;    p62 = 1.0466e-3;
    p70 = 11747;      p71 = 41.66;   p72 = 4521;       p73 = 0.01998;   p74 = 0.00604;   p75 = 2534;      p76 = 0.00475;
    p80 = 1.543e-2;   p81 = 10.06;   p82 = 1.807e-2;   p83 = 0.0011;    p84 = 8.81e-5;   p85 = 0.0405;    p86 = 20.24;
    p90 = 8.321;      p91 = 0.14249; p92 = 0.96655;    p93 = 26.42;     p94 = 0.0567;


    # Numerical determination of the penetration depth (86%) (Eq. 8)
    D86 = 1/rhob*(p90+p91*(p92+np.exp(-1*distances/100))*(p93+theta)/(p94+theta))

    # Depth weights (Eq. 7)
    Wd = np.exp(-2*depth/D86)

    if h == 0:
        W = 1 # skip distance weighting

    elif (h >= 0.1) and (h<= 50):
        # Functions for Fi (Appendix A in Köhli et al., 2015)
        F1 = p10*(1+p13*h)*np.exp(-p11*theta)+p12*(1+p15*h)-p14*theta
        F2 = ((-p20+p24*h)*np.exp(-p21*theta/(1+p25*theta))+p22)*(1+h*p23)
        F3 = (p30*(1+p33*h)*np.exp(-p31*theta)+p32-p34*theta)
        F4 = p40*np.exp(-p41*theta)+p42-p43*theta+p44*h
        F5 = p50*(0.02-1/p55/(h-p55+p56*theta))*(p54-theta)*np.exp(-p51*(theta-p54))+p52*(0.7-h*theta*p53)
        F6 = p60*(h+p61)+p62*theta
        F7 = (p70*(1-p76*h)*np.exp(-p71*theta*(1-h*p74))+p72-p75*theta)*(2+h*p73)
        F8 = ((-p80+p84*h)*np.exp(-p81*theta/(1+p85*h+p86*theta))+p82)*(2+h*p83)

        # Distance weights (Eq. 3)
        W = np.ones_like(distances)*np.nan
        for i in range(len(distances)):
            if (distances[i]<=50) and (distances[i]>0.5):
                W[i]=F1[i]*(np.exp(-F2[i]*distances[i]))+F3[i]*np.exp(-F4[i]*distances[i])

            elif (distances[i]>50) and (distances[i]<600):
                W[i]=F5[i]*(np.exp(-F6[i]*distances[i]))+F7[i]*np.exp(-F8[i]*distances[i])

            else:
                raise ValueError('Input distances are not valid.')

    else:
        raise ValueError('Air humidity values are out of range.')


    # Combined and normalized weights
    weights = Wd*W/np.nansum(Wd*W)

    return weights



def storage(sm,T=1,Z_surface=150,Z_subsurface=1000):
    """Exponential filter to estimate soil moisture in the rootzone from surface observtions.

    Args:
        sm (list or array): Soil moisture in mm of water.
        T (float): Characteristic time length in the same units as the measurement interval.
        Z_surface (float): Depth of surface layer in mm. This should be an intermediate value according to the
            sensing depth computed using the D86 method.
        Z_subsurface (float): Depth of subsurface layer in mm.

    Returns:
        (tuple): tuple containing:
            - **Surface soil water storage** (*array*): Surface soil water storage in mm of water.
            - **Subsurface soil water storage** (*array*): Subsurface soil water storage in mm of water.

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
    t_delta = 1
    sm_min = np.min(sm)
    sm_max = np.max(sm)
    ms = (sm-sm_min)/(sm_max-sm_min)

    # Pre-allocate soil water index array
    SWI = np.ones_like(ms)*np.nan
    K = np.ones_like(ms)*np.nan

    # Initial conditions
    SWI[0] = ms[0]
    K[0] = 1

    # Values from 2 to N
    for n in range(1,len(SWI)):
        if ~np.isnan(ms[n]) & ~np.isnan(ms[n-1]):
            K[n] = K[n-1] / (K[n-1] + np.exp(-t_delta/T))
            SWI[n] = SWI[n-1] + K[n]*(ms[n] - SWI[n-1])
        else:
            continue

    # Surface storage
    storage_surface = sm*Z_surface

    # Rootzone storage
    storage_subsurface = (SWI*(sm_max-sm_min) + sm_min)*Z_subsurface

    return storage_surface, storage_subsurface


def cutoff_rigidity(lat,lon):
    """Function to estimate the approximate cutoff rigidity for any point on Earth according to the
    tabulated data of Smart and Shea, 2019. Values are approximations so that users have an idea of
    what neutron detectors from the Neutron Monitor Database (NMD).

    Args:
        lat (float): Geographic latitude in decimal degrees. Value in range -90 to 90
        lon (float): Geographic longitude in decimal degrees. Values in range from 0 to 360.
            Typical negative longitudes in the west hemisphere will fall in the range 180 to 360.

    Returns:
        (float): Cutoff rigidity in GV. Error is about +/- 0.3 GV

    Examples:
        Estimate the cutoff rigidity for Newark, NJ, US

        >>> zq = cutoff_rigidity(39.68, -75.75)
        >>> print(zq)
        2.52 GV (Value from NMD is 2.40 GV)

    References:
        Smart, D. & Shea, Matthew. (2001). Geomagnetic Cutoff Rigidity Computer Program:
        Theory, Software Description and Example. NASA STI/Recon Technical Report N.

        Shea, M. A., & Smart, D. F. (2019, July). Re-examination of the First Five Ground-Level Events.
        In International Cosmic Ray Conference (ICRC2019) (Vol. 36, p. 1149).
    """
    xq = lon
    yq = lat

    if xq < 0:
        xq = xq*-1 + 180
    Z = np.array(data.cutoff_rigidity)
    x = np.linspace(0, 360, Z.shape[1])
    y = np.linspace(90, -90, Z.shape[0])
    X, Y = np.meshgrid(x, y)
    points = np.array( (X.flatten(), Y.flatten()) ).T
    values = Z.flatten()
    zq = griddata(points, values, (xq,yq))

    return np.round(zq,2)


def find_neutron_detectors(Rc, start_date=None, end_date=None):
    """Search for potential reference neutron monitoring stations based on cutoff rigidity.
    
    Args:
        Rc (float): Cutoff rigidity in GV. Values in range 1.0 to 3.0 GV.
        start_date (datetime): Start date for the period of interest.   
        end_date (datetime): End date for the period of interest.
        
    Returns:
        (list): List of top five stations with closes cutoff rigidity.
            User needs to select station according to site altitude.
            
    Examples:
        >>> from crnpy import crnpy
        >>> Rc = 2.40 # 2.40 Newark, NJ, US
        >>> crnpy.find_neutron_detectors(Rc)
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
    stations = pd.DataFrame(data.neutron_detectors, columns=["STID","NAME","R","Altitude_m"])

    # Sort stations by closest cutoff rigidity
    idx_R = (stations['R'] - Rc).abs().argsort()

    if start_date is not None and end_date is not None:
        stations["Period available"] = False
        for i in range(10):
            station = stations.iloc[idx_R[i]]["STID"]
            try:
                if get_incoming_neutron_flux(start_date, end_date, station, verbose=-1) is not None:
                    stations.iloc[idx_R[i],-1] = True
            except:
                pass
        if sum(stations["Period available"] == True) == 0:
            print("No stations available for the selected period!")
        else:
            stations = stations[stations["Period available"] == True]
            idx_R = (stations['R'] - Rc).abs().argsort()
            result = stations.iloc[idx_R.iloc[:10]]
    else:
        result = stations.reindex(idx_R).head(10).rename_axis(None)

    # Print results
    print('')
    print("""Select a station with an altitude similar to that of your location. For more information go to: 'https://www.nmdb.eu/nest/help.php#helpstations""")
    print('')
    print(f"Your cutoff rigidity is {Rc} GV")
    print(result)
    return result


def estimate_lattice_water(clay_content, total_carbon=None):
    r"""Estimate the amount of water in the lattice of clay minerals.

    ![img1](img/lattice_water_simple.png) | ![img2](img/lattice_water_multiple.png)
    :-------------------------:|:-------------------------:
    $\omega_{lat} = 0.097 * clay(\%)$ | $\omega_{lat} = -0.028 + 0.077 * clay(\%) + 0.459 * carbon(\%)$
    Linear regression [lattice water (%) as a function of clay (%)] done with data from Soil Water Processes Lab and Dong and Ochsner (2018) |  Multiple linear regression [lattice water (%) as a function of clay (%) and soil carbon (%)] done with data from Soil Water Processes Lab.

    Args:
        clay_content (float): Clay content in the soil in percent.
        total_carbon (float, optional): Total carbon content in the soil in percent.
            If None, the amount of water is estimated based on clay content only.

    Returns:
        (float): Amount of water in the lattice of clay minerals in percent

    References:
        Dong, J., & Ochsner, T. E. (2018). Soil texture often exerts a stronger influence than precipitation
         on mesoscale soil moisture patterns. Water Resources Research, 54, 2199– 2211.
         https://doi.org/10.1002/2017WR021692

    """
    if total_carbon is None:
        lattice_water = 0.097 * clay_content
    else:
        lattice_water = -0.028 + 0.077 * clay_content + 0.459 * total_carbon
    return lattice_water


def latlon_to_utm(lat, lon, utm_zone_number, missing_values=None):
    """Convert geographic coordinates (lat, lon) to projected coordinates (utm) using the Military Grid Reference System.

    Function only applies to non-polar coordinates.
    If further functionality is required, consider using the utm module. See references for more information.

    ![UTM zones](https://upload.wikimedia.org/wikipedia/commons/thumb/b/b7/Universal_Transverse_Mercator_zones.svg/1920px-Universal_Transverse_Mercator_zones.svg.png)
    UTM zones on an equirectangular world map with irregular zones in red and New York City's zone highlighted. See [UTM zones](https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system#UTM_zones) for a full description.


    Args:
        lat (float, array): Latitude in decimal degrees.
        lon (float, array): Longitude in decimal degrees.
        utm_zone_number (int): Universal Transverse Mercator (UTM) zone.

    Returns:
        (float, float): Tuple of easting and northing coordinates in meters. First element is easting, second is northing.

    References:
         Code adapted from utm module created by Tobias Bieniek (Github username: Turbo87)
         [https://github.com/Turbo87/utm](https://github.com/Turbo87/utm)

         [https://www.maptools.com/tutorials/grid_zone_details#](https://www.maptools.com/tutorials/grid_zone_details#)
    """


    # Define constants
    R = 6_378_137  # Earth's radius at the Equator in meters

    # Convert input data to Numpy arrays
    if (type(lat) is not np.ndarray) or (type(lon) is not np.ndarray):
        try:
            lat = np.array(lat)
            lon = np.array(lon)
        except:
            raise "Input values cannot be converted to Numpy arrays."

    # Check latitude range
    if np.any(lat < -80) | np.any(lat > 84):
        raise "One or more latitude values exceed the range -80 to 84"

    # Check longitude range
    if np.any(lon < -180) | np.any(lon > 180):
        raise "One or more longitude values exceed the range -180 to 180"

    # Constants
    K0 = 0.9996
    E = 0.00669438
    E_P2 = E / (1 - E)

    M1 = (1 - E / 4 - 3 * E ** 2 / 64 - 5 * E ** 3 / 256)
    M2 = (3 * E / 8 + 3 * E ** 2 / 32 + 45 * E ** 3 / 1024)
    M3 = (15 * E ** 2 / 256 + 45 * E ** 3 / 1024)
    M4 = (35 * E ** 3 / 3072)

    # Trigonometric operations
    lat_rad = np.radians(lat)
    lon_rad = np.radians(lon)

    lat_sin = np.sin(lat_rad)
    lat_cos = np.cos(lat_rad)
    lat_tan = lat_sin / lat_cos
    lat_tan2 = lat_tan * lat_tan
    lat_tan4 = lat_tan2 * lat_tan2

    # Find central meridian.
    central_lon = (utm_zone_number * 6 - 180) - 3  # Zones are every 6 degrees.
    central_lon_rad = np.radians(central_lon)

    n = R / np.sqrt(1 - E * lat_sin ** 2)
    c = E_P2 * lat_cos ** 2

    with np.errstate(divide='ignore', invalid='ignore'):
        a = lat_cos * (np.remainder(((lon_rad - central_lon_rad) + np.pi), (2 * np.pi)) - np.pi)
    m = R * (M1 * lat_rad - M2 * np.sin(2 * lat_rad) + M3 * np.sin(4 * lat_rad) - M4 * np.sin(6 * lat_rad))

    easting = K0 * n * (a + a ** 3 / 6 * (1 - lat_tan2 + c) + a ** 5 / 120 * (
                5 - 18 * lat_tan2 + lat_tan4 + 72 * c - 58 * E_P2)) + 500_000
    northing = K0 * (m + n * lat_tan * (
                a ** 2 / 2 + a ** 4 / 24 * (5 - lat_tan2 + 9 * c + 4 * c ** 2) + a ** 6 / 720 * (
                    61 - 58 * lat_tan2 + lat_tan4 + 600 * c - 330 * E_P2)))

    if np.any(lat < 0):
        northing += 10_000_000

    return easting, northing

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


def smooth_2d(x, y, z, buffer=100, min_neighbours=3, method='mean', rnd=False):
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
        except:
            raise "Input values cannot be converted to Numpy arrays."

    if len(x) != len(y):
        raise f"The number of x and y must be equal. Input x has {len(x)} values and y has {len(y)} values."

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
                raise f"Method {method} does not exist. Provide either 'mean' or 'median'."
        else:
            z_new_val = z[k] # If there are not enough neighbours, keep the original value

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
    idx_nan = np.isnan(x) | np.isnan(y) | np.isnan(z)
    x = x[~idx_nan]
    y = y[~idx_nan]
    z = z[~idx_nan]

    if idx_nan.any():
        print(f"WARNING: {np.isnan(x).sum()}, {np.isnan(y).sum()}, and {np.isnan(z).sum()} NaN values were dropped from x, y, and z.")

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
        raise f"Method {method} does not exist. Provide either 'cubic', 'linear', 'nearest', or 'idw'."

    return X_pred, Y_pred, Z_pred



