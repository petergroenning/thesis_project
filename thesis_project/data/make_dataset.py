import pandas as pd
from thesis_project.data import DEPTHS_A_STRING, DEPTHS_B_STRING

DATA_PATH_RAW = 'data/raw/'
DATA_PATH_PROCESSED = 'data/processed/'
DATA_NAME = 'Damvarmelager temperatur målinger 2023-09.csv'

DEPTHS = {**DEPTHS_A_STRING, **DEPTHS_B_STRING}

def parse_data(data_path: str, data_name: str, columns: list = ['EVENT_TIME', 'POINT_NAME', 'VALUE_CUR']) -> pd.DataFrame:
    '''
    Parse raw data from csv file to pandas DataFrame

    Output format (example):
    | EVENT_TIME | POINT_NAME_X | POINT_NAME_Y | ... | POINT_NAME_Z |
    |------------|--------------|--------------|-----|--------------|
    | 2023-09-01 | 20.0         | 21.0         | ... | 22.0         |
    '''
    data = pd.read_csv(data_path + data_name, sep=';', header=0, parse_dates=True, usecols=columns)
    # Rename POINT_NAME colum (remove path)
    data.POINT_NAME = data.POINT_NAME.apply(lambda x: '_'.join(x.split('/')[5:]) if 'GRADER' not in x else x.split('/')[5])
    # Expand based on POINT_NAME
    data = data.pivot(index='EVENT_TIME', columns='POINT_NAME', values='VALUE_CUR')
    # Drop that one weird row (2023-09-04 14:00:00.000)
    data = data[~data.index.isin(['2023-09-04 14:00:00.000'])]
    # Drop zero column (not used)
    data = data[~(data == 0).all(axis = 1)]
    
  
    # Get Columns with temps (sorted):
    depths = sorted(DEPTHS.items(), key=lambda x:x[1])
    columns = [x[0] for x in depths]
    columns += ['FTDVTIL_M3H', 'FTDVFRA_M3H']
    data = data[columns]

    # Get t columns
    data.index = pd.to_datetime(data.index)
    data['t'] = (data.index - data.index[0]).total_seconds() / 3600


    depths = pd.DataFrame(depths, columns = ['POINT_NAME', 'DEPTH'])

    # rename columns
    data.rename(columns={'FTDVTIL_M3H': 'U_IN', 'FTDVFRA_M3H':'U_OUT'}, inplace=True)
    return data, depths

def add_depths(data_path:str, data_name:str, columns: list = ['EVENT_TIME', 'POINT_NAME', 'VALUE_CUR']) -> pd.DataFrame:
    '''
    Add depths to data DataFrame

    Output format (example):
    | EVENT_TIME | POINT_NAME_X | POINT_NAME_Y | ... | POINT_NAME_Z | DEPTH |
    |------------|--------------|--------------|-----|--------------|-------|
    | 2023-09-01 | 20.0         | 21.0         | ... | 22.0         | 0.25  |
    '''
    data = pd.read_csv(data_path + data_name, sep=';', header=0, parse_dates=True, usecols=columns)
    data['DEPTHS'] = data.POINT_NAME.apply(lambda x: DEPTHS.get(x.split('/')[5], None))

    return data

if __name__ == '__main__':
    data, depths = parse_data(DATA_PATH_RAW, DATA_NAME)
    data.to_csv(DATA_PATH_PROCESSED + 'data.csv')
    depths.to_csv(DATA_PATH_PROCESSED + 'depths.csv')
    # data = add_depths(DATA_PATH_RAW, DATA_NAME)
    # data.to_csv(DATA_PATH_PROCESSED + 'data_depths.csv')
    
