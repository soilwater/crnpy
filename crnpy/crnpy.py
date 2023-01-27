#!/usr/bin/env python3
# Created by Andres Patrignani Nov-2021

# Import modules
import sys
import warnings
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from scipy.interpolate import griddata, pchip_interpolate
import requests
import io, datetime, os



# Define python version
python_version = (3, 7)  # tuple of (major, minor) version requirement
python_version_str = str(python_version[0]) + "." + str(python_version[1])

# produce an error message if the python version is less than required
if sys.version_info < python_version:
    msg = "Module only runs on python version >= %s" % python_version_str
    raise Exception(msg)


def format_dates_df(df, col='timestamp', format='%Y-%m-%d %H:%M:%S', freq='H', round_time=True):
    """Helper function to change the format and round timestamps
    
    Parameters
    ----------
    round_time : str or None
        String denoting the rounding interval. 'H'=hourly, 'M'=minute, or None

    Returns
    -------
    DataFrame with formatted timestamps and rounded time.
    """
    
    # Change format of timestamp
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
      
def count_time(df):
    """Approximate counting time

    Parameters
    ----------
    count_cols : DataFrame
        Dataframe containing only the columns with neutron counts and timestamp in the index.

    Returns
    -------
    DataFrame
        Dataframe with the approximate counting time for each detector.
    """

    # Check that index is a timestamp
    if type(df.index) != pd.core.indexes.datetimes.DatetimeIndex:
        raise ValueError('Index must be a timestamp.')

    # Calculate time difference between rows
    count_time = df.index.to_series().diff().dt.total_seconds()

    return count_time

def fill_counts(df, count_times=None, expected_time=False, threshold=0.25, limit=3):
    """Fill missing neutron counts. Periods shorter than threshold are replaced with NaN.

    Parameters
    ----------
    df : DataFrame
        DataFrame with neutron counts, might have count_time column(s).

    count_time : pd.Series or pd.DataFrame
        Counting time in seconds. If a DataFrame is provided, it must have the same number of columns as df.
    count_times_normal : int
        Expected counting time in seconds. Default is 3600 seconds.
    threshold : float
        Minimum fraction of the neutron integration time. Default is 0.25.

    Returns
    -------
    DataFrame
        with linearly interpolated neutron counts.
    """
    df=df.copy()

    if type(df.index) == pd.core.indexes.datetimes.DatetimeIndex and isinstance(count_times, type(None)):
        warnings.warn("No count time columns provided. Using timestamp index to compute count time.")
        count_times = df.index.to_series().diff().dt.total_seconds()

    if type(df.index) != pd.core.indexes.datetimes.DatetimeIndex and isinstance(count_times, type(None)):
        raise ValueError('Index must be a timestamp or count times must be provided.')

    if len(df) != len(count_times):
        raise ValueError('Count times length does not match number of readings.')

    if expected_time is False:
        expected_time = count_times.median()
        print('Using median count time as expected count time:', expected_time)

    # Replace values below threshold with NaN
    time_threshold = round(expected_time * threshold)

    if type(count_times) == pd.core.frame.DataFrame:
        if len(count_times.columns) == 1:
            idx_nan = count_times[count_times < time_threshold].index
            df.loc[idx_nan] = np.nan
        else:
            for i in range(len(count_times.columns)):
                idx_nan = count_times[count_times.iloc[:,i] < time_threshold].index
                df.iloc[:,i].loc[idx_nan] = np.nan
    elif type(count_times) == pd.core.series.Series:
        idx_nan = count_times[count_times < time_threshold].index
        df.loc[idx_nan] = np.nan
    elif type(count_times) == np.ndarray:
        idx_nan = np.where(count_times < time_threshold)
        df.loc[idx_nan] = np.nan

    # Fill missing values with linear interpolation and round to nearest integer
    df = df.interpolate(method='linear', limit=limit, limit_direction='both').round()
    return df

def normalize_counts(df, count_time=3600, count_times=None):
    """Normalize neutron counts to the desired counting time.

    Parameters
    ----------
    df : DataFrame
        Dataframe containing only the columns with neutron counts.
    count_times : pd.Series or pd.DataFrame
        Counting time in seconds. If a DataFrame is provided, it must have the same number of columns as df.
    count_time : int
        Count time in seconds for normalization. Default is 3600 seconds.

    Returns
    -------
    DataFrame
        Normalized neutron counts.
    """

    if count_times is None and type(df.index) == pd.core.indexes.datetimes.DatetimeIndex:
        print("No count_times columns provided. Using timestamp index to compute count time.")
        count_times = df.index.to_series().diff().dt.total_seconds()

    if isinstance(count_times, type(None)):
        raise ValueError('Count time must be provided or index must be a timestamp.')

    if len(df) != len(count_times):
        raise ValueError('Count times length does not match number of readings.')


    #Normalize counts rounded to integer
    if type(count_times) == pd.core.series.Series or len(count_times.columns) == 1:
        df_r = df.div(count_times, axis=0).mul(count_time).round()
        return df_r
    else:
        df = df.copy()
        count_times = count_times.copy()
        for i in range(len(count_times.columns)):
            df.iloc[:,i] = df.iloc[:,i].div(count_times.iloc[:,i], axis=0).mul(count_time).round()
        return df



def compute_total_raw_counts(df, nan_strategy=None):
    """Compute the sum of uncorrected neutron counts for all detectors.
    
    Parameters
    ----------
    df : DataFrame
        Dataframe containing only the columns with neutron counts.
    
    Returns
    -------
    DataFrame
        Dataframe with the sum of uncorrected neutron counts for all detectors.
    """
    df=df.copy()

    if df.isnull().values.any():
        if nan_strategy is None:
            raise ValueError('NaN values found. Please fill missing values or provide a strategy. See documentation for more information.')
        elif nan_strategy == 'interpolate':
            print('NaN values found. Interpolating missing values using fill_counts().')
            if type(df.index) != pd.core.indexes.datetimes.DatetimeIndex:
                raise ValueError('Index must be a timestamp to use interpolation strategy.')
            df = fill_counts(df)
        elif nan_strategy == 'average':
            if len(df.columns) == 1:
                raise ValueError('Only one detector found. Cannot use average strategy.')
            print('NaN values found. Replacing missing values with average of other detectors before summing.')
            df = df.apply(lambda x: x.fillna(df.mean(axis=1)),axis=0)
        else:
            raise ValueError('Invalid strategy.')

    #Compute sum of counts
    df['total_raw_counts'] = df.sum(axis=1)
    # Replace zeros with NaN
    df['total_raw_counts'] = df['total_raw_counts'].replace(0, np.nan)
    return df['total_raw_counts']


def drop_outlier(df, window=11, store_outliers=False, min_counts=None, max_counts=None):
    """Computation of a moving modified Z-score based on the median absolute difference.
    
    References
    ----------
    Iglewicz, B. and Hoaglin, D.C., 1993. How to detect and handle outliers (Vol. 16). Asq Press.

    Parameters
    ----------
    cols_counts : DataFrame
        Dataframe containing only the columns with neutron counts.
    window : int
        Window size for the moving median. Default is 11.
    store_outliers : bool
        If True, store the outliers in a new column. Default is False.
    min_counts : int
        Minimum number of counts for a reading to be considered valid. Default is None.
    max_counts : int
        Maximum number of counts for a reading to be considered valid. Default is None.

    Returns
    -------
    DataFrame
        Dataframe without outliers.
    """
    if min_counts is not None:
        lower_count = np.sum(df < min_counts)
        print(f"Discarded counts below {min_counts}: {lower_count}")
        df = df[df >= min_counts]
    if max_counts is not None:
        upper_count = np.sum(df > max_counts)
        print(f"Discarded counts above {max_counts}: {upper_count}")
        df = df[df <= max_counts]

    # Compute median absolute difference
    median = df.rolling(window, center=True).median()
    diff = np.abs(df - median)
    mad = diff.rolling(window, center=True).median()

    # Compute modified Z-score
    modified_z_score = 0.6745 * diff / mad
    outliers = df[modified_z_score > 3.5]
    # Drop outliers
    df = df[modified_z_score < 3.5]

    if store_outliers:
        return df, outliers
    print(f"Discarded {len(outliers)} outliers using modified Z-score.")
    return df


def fill_missing_atm(cols_atm, limit=24):
    """Fill missing values in atmospheric variables. Gap filling is performed using a 
    piecewise cubic Hermite interpolating polynomial (pchip method) that is restricted to intervals
    of missing data with a limited number of values and surrounded by valid observations.
    There is no interpolation at the end of the time series.
    
    Parameters
    ----------
    col_atm : pd.Series or pd.DataFrame
        Atmospheric variables to fill.
    limit : int
        Maximum number of consecutive missing values to interpolate. Default is 24.

    Returns
    -------
    DataFrame
        Atmospheric variables with filled missing values using a piecewise cubic Hermite polynomial.
    """

    # Fill missing values in atmospheric variables
    return cols_atm.interpolate(method='pchip', limit=limit, limit_direction='both')

def atm_correction(counts, pressure, humidity, temp, Pref, Aref, L, incoming_neutrons=None, incoming_Ref=None):
    """Correct neutron counts for atmospheric factors and incoming neutron flux.
    
    Parameters
    ----------
    counts : list or array
        Neutron counts to correct.
    pressure : list or array
        Atmospheric pressure readings.
    humidity : list or array
        Atmospheric humidity readings in %.
    temp : list or array
        Atmospheric temperature readings in Celsius.
    Pref : float
        Reference atmospheric pressure in millibars (mbar).
    Aref : float
        Reference absolute humidity.
    L : float
        Atmospheric attenuation coefficient.
    incoming_neutrons : list or array
        Incoming neutron flux. Default is None.
    incoming_Ref : float
        Reference incoming neutron flux. Default is None.
    
    Returns
    -------
    Total neutron counts corrected by atmospheric conditions.
    
    References:
    Zreda, M., Shuttleworth, W. J., Zeng, X., Zweck, C., Desilets, D., Franz, T., et al. (2012). 
    COSMOS: the cosmic-ray soil moisture observing system. Hydrol. Earth Syst. Sci. 16, 4079–4099.
    doi: 10.5194/hess-16-4079-2012
    
    Hawdon, A., McJannet, D. and Wallace, J., 2014. Calibration and correction procedures for cosmic‐ray 
    neutron soil moisture probes located across Australia. Water Resources Research, 50(6), pp.5029-5043.
    
    Andreasen, M., Jensen, K.H., Desilets, D., Franz, T.E., Zreda, M., Bogena, H.R. and Looms, M.C., 2017. 
    Status and perspectives on the cosmic‐ray neutron method for soil moisture estimation and other 
    environmental science applications. Vadose zone journal, 16(8), pp.1-11. doi.org/10.2136/vzj2017.04.0086
    """
       
    ### Barometric pressure factor
    fp = np.exp((pressure - Pref) / L)

    ### Atmospheric water vapor factor
    # Saturation vapor pressure
    e_sat = 0.611 * np.exp(17.502 * temp / (temp + 240.97)) * 1000 # in Pascals Eq. 3.8 p.41 Environmental Biophysics (Campbell and Norman)
    
    # Vapor pressure Pascals
    Pw = e_sat * humidity/100;
    
    # Absolute humidity (g/m^3)
    C = 2.16679 # g K/J;
    A = C * Pw / (temp + 273.15)
    fw = 1 + 0.0054*(A - Aref);

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
    return np.round(counts*(fp*fw/fi))



def get_incoming_neutron_flux(timestamps, station, utc_offset=0, verbose=False):
    """Function to retrieve neutron flux from the Neutron Monitor Database.

    Parameters
    ----------
    timestamps : pd.series or pd.DataFrame or pd.DatetimeIndex
        Timestamps to retrieve neutron flux.
    station : str
        Neutron Monitor station to retrieve data from.
    utc_offset : int
        UTC offset in hours. Default is 0.

    Keyword arguments:
    utc_offset -- UTC offset in hours. Default is 0.
    
    Returns:
    Neutron flux in counts per hour and timestamps.
    
    References:
    Documentation available:https://www.nmdb.eu/nest/help.php#howto
    
    """

    # Example: get_incoming_flux(station='IRKT',start_date='2020-04-10 11:00:00',end_date='2020-06-18 17:00:00')
    # Template url = 'http://nest.nmdb.eu/draw_graph.php?formchk=1&stations[]=KERG&output=ascii&tabchoice=revori&dtype=corr_for_efficiency&date_choice=bydate&start_year=2009&start_month=09&start_day=01&start_hour=00&start_min=00&end_year=2009&end_month=09&end_day=05&end_hour=23&end_min=59&yunits=0'


    # Convert timestamps to datetime objects
    start_date = timestamps.min()
    end_date = timestamps.max()
    
    # Add 1 hour to ensure the last date is included in the request.
    end_date += pd.Timedelta(hours=1)

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
        return None, None

    if(len(timestamps) != len(df_flux) and verbose > -1):
        print('Warning: The number of timestamps does not match the number of neutron flux values.')
        print('Check time resolution, or try the interpolate_incoming_neutron_flux() method.')

    # Check if all values from selected detector are NaN. If yes, warn the user
    if df_flux['counts'].isna().all():
        warnings.warn('Data for selected neutron detectors appears to be unavailable for the selected period')
    

    # Print acknowledgement to inform users about restrictions and to acknowledge the NMDB database
    acknowledgement = """Data retrieved via NMDB are the property of the individual data providers. These data are free for non commercial
use to within the restriction imposed by the providers. If you use such data for your research or applications, please acknowledge
the origin by a sentence like 'We acknowledge the NMDB database (www.nmdb.eu) founded under the European Union's FP7 programme 
(contract no. 213007), and the PIs of individual neutron monitors at: IGY Jungfraujoch 
(Physikalisches Institut, University of Bern, Switzerland)"""
    #print(acknowledgement)

    return df_flux

def interpolate_incoming_flux(df_flux, timestamps):
    """Function to interpolate incoming neutron flux.

    Parameters
    ----------
    df_flux : pd.DataFrame
        Dataframe returned by get_incoming_flux method.
    timestamps : pd.series or pd.DataFrame or pd.DatetimeIndex
        Timestamps to interpolate the incoming neutron flux.

    Returns
    -------
    pd.DataFrame
        Dataframe containing interpolated incoming neutron flux.
    """

    # Add timestamps with nan values to the dataframe
    df_flux['timestamp'] = pd.to_datetime(df_flux['timestamp'])
    df_flux = df_flux.set_index('timestamp')
    for timestamp in timestamps:
        if timestamp not in df_flux.index:
            df_flux.loc[timestamp] = np.nan
    df_flux = df_flux.sort_index()

    # Interpolate nan values
    df_flux = df_flux['counts'].interpolate(method='pchip')

    # Retur only the values for the selected timestamps
    return df_flux.loc[timestamps]


def smooth_counts(df,window=11,order=3):
    """Use a Savitzky-Golay filter to smooth the signal of corrected neutron counts.
    
    Parameters
    ----------
    df : DataFrame
        Dataframe containing the corrected neutron counts.
    window : int
        Window size for the Savitzky-Golay filter. Default is 11.
    order : int
        Order of the Savitzky-Golay filter. Default is 3.

    Returns
    -------
    DataFrame with smoothed neutron counts.


    References:
    Franz, T.E., Wahbi, A., Zhang, J., Vreugdenhil, M., Heng, L., Dercon, G., Strauss, P., Brocca, L. and Wagner, W., 2020. Practical data products from cosmic-ray neutron sensing for hydrological applications. Frontiers in Water, 2, p.9. doi.org/10.3389/frwa.2020.00009
    """
    if df.isna().all().any():
        raise ValueError('Dataframe contains NaN values. Please remove NaN values before smoothing the data.')

    if type(df) == pd.core.series.Series:
        filtered = np.round(savgol_filter(df,window,order))
        df = pd.DataFrame(filtered,columns=['counts'], index=df.index)
    elif type(df) == pd.core.frame.DataFrame:
        for col in df.columns:
            df[col] = np.round(savgol_filter(df[col],window,order))

    return df

def counts_to_vwc(counts, N0, Wlat, Wsoc ,bulk_density, a0=0.0808,a1=0.372,a2=0.115):
    """Function to convert corrected and filtered neutron counts into volumetric water content
    following the approach described in Desilets et al., 2010.

    Parameters
    ----------
    counts : array or pd.Series or pd.DataFrame
        Array of corrected and filtered neutron counts.
    N0 : float
        Device-specific neutron calibration constant.
    Wlat : float
        Lattice water content.
    Wsoc : float
        Soil organic carbon content.
    bulk_density : float
        Soil bulk density.

    Keyword arguments:

    a0, a1, a2 -- Parameters given in Zreda et al., 2012.
    
    References: 
    Desilets, D., M. Zreda, and T.P.A. Ferré. 2010. Nature’s neutron probe: 
    Land surface hydrology at an elusive scale with cosmic rays. Water Resour. Res. 46:W11505.
    doi.org/10.1029/2009WR008726
    """
    
    # Convert neutron counts into vwc
    vwc = (a0 / (counts/N0-a1) - a2 - Wlat - Wsoc) * bulk_density
    return vwc


def fill_missing_vwc(df,limit=24):
    """Fill missing values in volumetric water content using a piecewise cubic Hermite 
    interpolating polynomial (pchip method)

    Keyword arguments:
    
    df -- This is the main DataFrame with tabular data to correct
    limit -- Maximum number of consecutive missing values to fill
    """
    
    # Interpolate rows with missing volumetric water content.
    df.interpolate(method='pchip', limit_area='inside', limit=limit, inplace=True)
    return df


def sensing_depth(df,par,method='Schron_2017',dist=[1]):
    """Function that computes the estimated sensing depth of the cosmic-ray neutron probe.
    The function offers several methods to compute the depth at which 86 % of the neutrons
    probes the soil profile.
    
    Keyword arguments:
    
    df -- This is the main DataFrame with tabular data to correct.
    par -- User-defined arguments for the particular experiment or field
    method -- 'Schron_2017' or 'Franz_2012'
    dist -- List of radial distances at which to estimate the sensing depth. Only used for the 'Schron_2017' method.
    
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
        Fp = 0.4922 / (0.86 - np.exp(-1*df['pressure']/par['Pref']));
        Fveg = 0
        
        for d in dist:
            # Compute r_star
            r_start = d/Fp 
            
            # Compute soil depth that accounts for 86% of the neutron flux
            col_name = 'sensing_depth_D86_' + str(d) + 'm'
            df[col_name] = 1/par['bulk_density'] * (8.321+0.14249*(0.96655 + np.exp(-0.01*r_start))*(20+(par['Wlat']+df['vwc'])) / (0.0429+(par['Wlat']+df['vwc'])))

    elif method == 'Franz_2012':
        df['sensing_depth_D86'] = 5.8/(par['bulk_density']*par['Wlat']+ df['vwc']+0.0829)
        
    return df
    
    
def nrad_weight(h,theta,distances,depth,rhob=1.4):
    """Function to compute distance weights corresponding to each soil sample.
    
    Keyword arguments:
    
    h -- Air Humidity  from 0.1  to 50    in g/m^3. When h=0, the function will skip the distance weighting.
    theta -- Soil Moisture for each sample (0.02 - 0.50 m^3/m^3)
    distances -- Distances from the location of each sample to the origin (0.5 - 600 m)
    depth -- Depths for each sample (m)
    rhob -- Bulk density in g/cm^3

    Example: 
    W = nrad_weight(5,np.array([0.25,0.25,0.25]),np.array([5,10,150]),np.array([0.05,0.05,0.05]),rhob=1.4)
    
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
    
    Parameters
    ----------
    df : This is the main DataFrame with tabular data to correct.
    Zsurface : Depth of surface layer in mm. This should be an intermediate value according to the 
                sensing depth computed using the D86 method.
    T : Characteristic time length in the same units as the measurement interval    
    
    Returns
    -------
    Surface and subsurface soil water storage in mm of water.
    
    References
    ----------
- Albergel, C., Rüdiger, C., Pellarin, T., Calvet, J.C., Fritz, N., Froissard, F., Suquia, D., Petitpa, A., Piguet, B. and Martin, E., 2008. From near-surface to root-zone soil moisture using an exponential filter: an assessment of the method based on in-situ observations and model simulations. Hydrology and Earth System Sciences, 12(6), pp.1323-1337.
    
    - Franz, T.E., Wahbi, A., Zhang, J., Vreugdenhil, M., Heng, L., Dercon, G., Strauss, P., Brocca, L. and Wagner, W., 2020. Practical data products from cosmic-ray neutron sensing for hydrological applications. Frontiers in Water, 2, p.9.
    
    - Rossini, P. and Patrignani, A., 2021. Predicting rootzone soil moisture from surface observations in cropland using an exponential filter. Soil Science Society of America Journal.
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
    
    Inputs:
        - Geographic latitude in decimal degrees. Value in range -90 to 90
        - Geographic longitude in decimal degrees. Values in range from 0 to 360. 
          Typical negative longitudes in the west hemisphere will fall in the range 180 to 360.
          
    Outputs:
        - Cutoff rigidity in GV. Error is about +/- 0.3 GV
        
    Example:
        Estimate the cutoff rigidity for Newark, NJ, US
        zq = cutoff_rigidity(39.68, -75.75)
        
        2.52 GV (Value from NMD is 2.40 GV)
    """
    xq = lon
    yq = lat

    if xq < 0:
        xq = xq*-1 + 180
    this_dir, this_filename = os.path.split(__file__)
    DATA_PATH = os.path.join(this_dir, "global_cutoff_rigidity_2015.csv")
    Z = np.loadtxt(open(DATA_PATH, "rb"), delimiter=",", skiprows=3)
    x = np.linspace(0, 360, Z.shape[1])
    y = np.linspace(90, -90, Z.shape[0])
    X, Y = np.meshgrid(x, y)

    points = np.array( (X.flatten(), Y.flatten()) ).T
    values = Z.flatten()
    zq = griddata(points, values, (xq,yq))
   
    return np.round(zq,2)


def find_neutron_detectors(Rc, timestamps=None):
    """
    Inputs:
        - Rc = Cutoff rigidity in GV
        
    Outputs:
        -List of top five stations with closes cutoff rigidity. 
         User needs to select station according to site altitude.
         
    Example:
        Rc = 2.40 # 2.40 Newark, NJ, US
        suggested_neutron_detectors(Rc)
        
         
    Source: https://www.nmdb.eu/nest/help.php#helpstations. Last accessed 20-Dec-2021
    """
    # Load file with list of neutron monitoring stations
    this_dir, this_filename = os.path.split(__file__)
    DATA_PATH = os.path.join(this_dir, "global_neutron_detectors.csv")
    stations = pd.read_csv(DATA_PATH, skiprows=1)

    # Sort stations by closest cutoff rigidity
    idx_R = (stations['R'] - Rc).abs().argsort()

    if timestamps is not None:
        stations["Period available"] = False
        for i in range(10):
            station = stations.iloc[idx_R[i]]["STID"]
            try:
                get_incoming_neutron_flux(timestamps, station, verbose=-1)
                stations.iloc[idx_R[i],-1] = True
            except:
                pass
        if sum(stations["Period available"] == True) == 0:
            print("No stations available for the selected period!")
        else:
            stations = stations[stations["Period available"] == True]
            idx_R = (stations['R'] - Rc).abs().argsort()
            result = stations.iloc[idx_R[:10]]
    else:
        result = stations.reindex(idx_R).head(10).rename_axis(None)

    # Print results
    print('')
    print("""Select a station with an altitude similar to that of your location. For more information go to: 'https://www.nmdb.eu/nest/help.php#helpstations""")
    print('')
    print(f"Your cutoff rigidity is {Rc} GV")
    print(result)
    