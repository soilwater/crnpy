#!/usr/bin/env python3
# Created by Andres Patrignani Nov-2021

# Import modules
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from scipy.interpolate import griddata
import warnings
import requests
import io
import sys

# Define python version
python_version = (3, 7)  # tuple of (major, minor) version requirement
python_version_str = str(python_version[0]) + "." + str(python_version[1])

# produce an error message if the python version is less than required
if sys.version_info < python_version:
    msg = "Module only runs on python version >= %s" % python_version_str
    raise Exception(msg)
    

def read_data(filename,variables,skip_rows=0,round_time=None):
    """Function that reads a .csv file with raw neutron counts and associated atmospheric variables.

    Parameters
    ----------
    filename: str
        Name of comma-separated value file with raw sensor data (.csv, .txt, .dat)
    skip_rows: int
        Positive integer representing the number of lines to skip (e.g. comments)
    variables:
        Dictionary of lists with the column position for each pre-defined variable (1-index)
    par: dict
        User-defined parameters for specific field and device
    round_time: str or None
        String denoting the rounding of the timestamp. 'H'=hourly, 'M'=minute, or None

    Returns
    -------
    DataFrame
        Pandas DataFrame with the selected variables.
    """

    # Extract column names and rename if necessary (eg. more than one detector)
    col_names = []
    col_positions = []
    for var in variables.keys():
        for k,col in enumerate(variables[var]):
            col_positions.append(col-1) # Variable positions are 1-index
            if len(variables[var])>1:
                col_names.append(var + '_' + str(k+1))
            else:
                col_names.append(var)

    # Read data from .CSV file
    df = pd.read_csv(filename, skiprows=skip_rows, usecols=col_positions)

    # Re-order column names and assign name to column
    df.columns = [x for _,x in sorted(zip(col_positions, col_names))]
    df['timestamp'] = pd.to_datetime(df['timestamp'], format='%Y/%m/%d %H:%M:%S')

    # Round timestamps to nearest hour
    if round_time == 'H':
        df['timestamp'] = df['timestamp'].dt.round('H')
    elif round_time == 'M':
        df['timestamp'] = df['timestamp'].dt.round('M')

    return df


def fill_counts_v2(source, cols_count, cols_time, count_time=3600, threshold=0.25):
    """Fill missing neutron counts. Periods shorter than threshold are replaced with NaN.

    Parameters
    ----------
    threshold : float
        Minimum fraction of the neutron integration time. Default is 0.25.

    Returns
    -------
    DataFrame
        with linearly interpolated neutron counts.
    """

    # Find columns with count values and count time
    time_threshold = round(count_time*threshold)

    for k in range(len(cols_count)):
        idx_nan = source[cols_time[k]] < time_threshold
        source.loc[idx_nan,cols_count[k]] = np.nan
        source.loc[idx_nan,cols_time[k]] = np.nan

        idx_fill = source[cols_time[k]] > time_threshold
        source.loc[idx_fill,cols_count[k]] = (source.loc[idx_fill,cols_count[k]]*count_time/source.loc[idx_fill,cols_time[k]]).round()
        source.loc[idx_fill,cols_time[k]] = count_time

    return source

       
        
    
def fill_counts(threshold=0.25):
    """Fill missing neutron counts. Periods shorter than threshold are replaced with NaN.

    Parameters
    ----------
    threshold : float
        Minimum fraction of the neutron integration time. Default is 0.25.

    Returns
    -------
    DataFrame
        with column for neutron counts for individual detectors with filled counts using linear interpolation. 
    """

    # Find columns with count values and count time
    cols_count = [col for col in df.columns if "counts_raw" in col]
    cols_time = [col for col in df.columns if "counts_time" in col]
    time_threshold = par['count_time']*threshold

    if len(cols_count) == len(cols_time):
        try:
            for k in range(len(cols_count)):
                idx_nan = df[cols_time[k]] < time_threshold
                df.loc[idx_nan,cols_count[k]] = np.nan
                df.loc[idx_nan,cols_time[k]] = np.nan

                idx_fill = (df[cols_time[k]] >= time_threshold) & (df[cols_time[k]] < self.par['count_time'])
                df.loc[idx_fill,cols_count[k]] = (df.loc[idx_fill,cols_count[k]] * self.par['count_time'] / df.loc[idx_fill,cols_time[k]]).round()
                df.loc[idx_fill,cols_time[k]] = self.par['count_time']

            return df

        except ValueError as error:
            logger.error(error)
            raise
    else:
        raise ValueError("Number of column counts and number of column times don't match")


def compute_total_counts(source, cols_counts):
    """Compute the sum of uncorrected neutron counts for all detectors.
    
    Parameters
    ----------
    df : This is the main DataFrame with tabular data to correct.
    
    Returns:
    Sum of neutron counts for all detectors.
    """
    
    df['counts_total_raw'] = df[cols_count].sum(skipna=True, axis=1, min_count=1)
    return df


def fill_missing_atm(df,limit=24):
    """Fill missing values in atmospheric variables. Gap filling is performed using a 
    piecewise cubic Hermite interpolating polynomial (pchip method) that is restricted to intervals
    of missing data with a limited number of values and surrounded by valid observations.
    There is no interpolation at the end of the time series.
    
    Keyword arguments:
    df -- This is the main DataFrame with tabular data to correct.
    limit -- Positive integer denoting the largest gap in missing data that will be interpolated.
    
    Returns:
    Atmospheric varisbales with filled missing values using a piecewise cubic Hermite polynomial.
    """
    
    df['temperature'].interpolate(method='pchip', limit_area='inside', limit=limit, inplace=True)
    df['humidity'].interpolate(method='pchip', limit_area='inside', limit=limit, inplace=True)
    df['pressure'].interpolate(method='pchip', limit_area='inside', limit=limit, inplace=True)
    return df


def atm_correction(df,par):
    """Correct neutron counts for atmospheric factors and incoming neutron flux.
    
    Keyword arguments:
    df -- This is the main DataFrame with tabular data to correct.
    par -- User-defined parameters for the particular experiment or field
    
    Returns:
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
    fp = np.exp((df['pressure'] - par['Pref'])/par['L']); # Pressure factor (fp) [Eq. 1] Hawdon et al. 2014
    
    ### Atmospheric water vapor factor
    # Saturation vapor pressure 
    e_sat = 0.611*np.exp(17.502*df['temperature'] / (df['temperature']+240.97))*1000 # in Pascals Eq. 3.8 p.41 Environmental Biophysics (Campbell and Norman)
    
    # Vapor pressure Pascals
    Pw = e_sat * df['humidity']/100; 
    
    # Absolute humidity (g/m^3)
    C = 2.16679 # g K/J;
    A = C * Pw / (df['temperature']+273.15) 
    fw = 1 + 0.0054*(A - par['Aref']);

    ### Incoming neutron flux factor
    if 'incoming_counts' not in df.columns:
        fi = 1
        warnings.warn("Ignoring incoming neutron flux correction factor (using value fi=1)")
    else:
        fi = df['incoming_counts']/df['incoming_counts'].iloc[0] # Use first value as reference
        fi.fillna(1.0, inplace=True) # Use a value of 1 for days without data 

    # Apply correction factors
    df['counts_total_corrected'] = (df['counts_total_raw']*(fp*fw/fi)).round(0)
    
    return df


def get_incoming_neutron_flux(df,par):
    """Function to retrieve neutron flux from the Neutron Monitor Database.
    
    Keyword arguments:
    par -- User-defined parameters for the particular experiment or field
    
    Returns:
    Neutron flux in counts per hour
    
    References:
    Documentation available:https://www.nmdb.eu/nest/help.php#howto
    
    """

    # Example: get_incoming_flux(station='IRKT',start_date='2020-04-10 11:00:00',end_date='2020-06-18 17:00:00')
    # Template url = 'http://nest.nmdb.eu/draw_graph.php?formchk=1&stations[]=KERG&output=ascii&tabchoice=revori&dtype=corr_for_efficiency&date_choice=bydate&start_year=2009&start_month=09&start_day=01&start_hour=00&start_min=00&end_year=2009&end_month=09&end_day=05&end_hour=23&end_min=59&yunits=0'

    start_date = df['timestamp'].iloc[0]
    end_date = df['timestamp'].iloc[-1]
    
    # Add 1 hour to ensure the last date is included in the request.
    end_date += pd.Timedelta('1H') 

    # Convert local time to UTC
    if par['utc_offset']>0:
        start_date -= pd.Timedelta( str(par['utc_offset']) + 'H')
        end_date -= pd.Timedelta( str(par['utc_offset']) + 'H')
    else:
        start_date += pd.Timedelta( str(par['utc_offset']) + 'H')
        end_date += pd.Timedelta( str(par['utc_offset']) + 'H')

    date_format = '%Y-%m-%d %H:%M:%S'
    root = 'http://www.nmdb.eu/nest/draw_graph.php?'
    url_par = [ 'formchk=1',
                'stations[]=' + par['nmdb_station'],
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
                'yunits=0' ]

    url = root + '&'.join(url_par)
    r = requests.get(url).content.decode('utf-8')

    # Subtract 1 hour to restore the last date included in the request.
    end_date -= pd.Timedelta('1H') 

    start = r.find('\n' + start_date.strftime(format=date_format))
    end = r.find('\n' + end_date.strftime(format=date_format)) - 1
    s = r[start:end]
    s = r[start:end]
    s2 = ''.join([row.replace(';',',') for row in s])
    df_flux = pd.read_csv(io.StringIO(s2), names=['timestamp','counts'])
    df['incoming_counts'] = df_flux['counts']
    
    # Check if all values from selected detector are NaN. If yes, warn the user
    if df['incoming_counts'].isna().all():
        warnings.warn('Data for selected neutron detectors appears to be unavailable for the selected period')
    

    # Print acknowledgement to inform users about restrictions and to acknowledge the NMDB database
    acknowledgement = """Data retrieved via NMDB are the property of the individual data providers. These data are free for non commercial
use to within the restriction imposed by the providers. If you use such data for your research or applications, please acknowledge
the origin by a sentence like 'We acknowledge the NMDB database (www.nmdb.eu) founded under the European Union's FP7 programme 
(contract no. 213007), and the PIs of individual neutron monitors at: IGY Jungfraujoch 
(Physikalisches Institut, University of Bern, Switzerland)"""
    #print(acknowledgement)

    return df


def smooth_counts(df,window=11,order=3):
    """Use a Savitzky-Golay filter to smooth the signal of corrected neutron counts.
    
    Keyword arguments:
    
    df -- This is the main DataFrame with tabular data to correct
    window --  Window of the filter, typically odd integer [Default is 11]
    order -- Order of the polynomial. [Default is 3]
    
    References:
    Franz, T.E., Wahbi, A., Zhang, J., Vreugdenhil, M., Heng, L., Dercon, G., Strauss, P., Brocca, L. and Wagner, W., 2020. Practical data products from cosmic-ray neutron sensing for hydrological applications. Frontiers in Water, 2, p.9. doi.org/10.3389/frwa.2020.00009
    """
    
    # Smooth data
    mode = 'interp'
    df['counts_total_filtered'] = np.round(savgol_filter(df['counts_total_corrected'], window, order, mode=mode))

    return df


def counts_to_vwc(a0=0.0808,a1=0.372,a2=0.115):
    """Function to convert corrected and filtered neutron counts into volumetric water content
    following the approach described in Desilets et al., 2010.
    
    Keyword arguments:
    
    df -- This is the main DataFrame with tabular data to correct
    a0, a1, a2 -- Parameters given in Zreda et al., 2012.
    
    References: 
    Desilets, D., M. Zreda, and T.P.A. Ferré. 2010. Nature’s neutron probe: 
    Land surface hydrology at an elusive scale with cosmic rays. Water Resour. Res. 46:W11505.
    doi.org/10.1029/2009WR008726
    """
    
    # Convert neutron counts into vwc
    vwc = (a0 / (df['counts_total_filtered']/par['N0']-a1) - a2 - par['Wlat'] - par['Wsoc'])*par['bulk_density']
    return df


def fill_missing_vwc(df,limit=24):
    """Fill missing values in volumetric water content using a piecewise cubic Hermite 
    interpolating polynomial (pchip method)

    Keyword arguments:
    
    df -- This is the main DataFrame with tabular data to correct
    limit -- Maximum number of consecutive missing values to fill
    """
    
    # Interpolate rows with missing volumetric water content.
    df['vwc'].interpolate(method='pchip', limit_area='inside', limit=limit, inplace=True)
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

    
def storage(df,T=1,Z_surface=150,Z_subsurface=1000):
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
    sm_min = df['vwc'].min()
    sm_max = df['vwc'].max()
    ms = (df['vwc']-sm_min)/(sm_max-sm_min)
    
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
    df['storage_surface'] = df['vwc']*Z_surface
    
    # Rootzone storage
    df['storage_subsurface'] = (SWI*(sm_max-sm_min) + sm_min)*Z_subsurface
    
    return df
    
    
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

    Z = np.loadtxt(open("global_cutoff_rigidity_2015.csv", "rb"), delimiter=",", skiprows=3)
    x = np.linspace(0, 360, Z.shape[1])
    y = np.linspace(90, -90, Z.shape[0])
    X, Y = np.meshgrid(x, y)

    points = np.array( (X.flatten(), Y.flatten()) ).T
    values = Z.flatten()
    zq = griddata(points, values, (xq,yq))
   
    return np.round(zq,2)


def suggested_neutron_detectors(Rc):
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
    stations = pd.read_csv("global_neutron_detectors.csv", skiprows=1)

    # Sort stations by closest cutoff rigidity
    idx_R = (stations['R'] - Rc).abs().argsort()
    
    # Print results
    print(stations.reindex(idx_R).head().rename_axis(None))
    print("""Select a station with an altitude similar to that of your location. 
    For more information go to: 'https://www.nmdb.eu/nest/help.php#helpstations""")