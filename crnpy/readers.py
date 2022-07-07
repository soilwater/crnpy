
def read_data(self,filename,skip_rows,variables,round_time=None):
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
    df: Object
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