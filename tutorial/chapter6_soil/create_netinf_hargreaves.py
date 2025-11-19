"""
Creates net infiltration using climate data, 
i.e. precipiation + evapotranspiration using Hargreaves based on temperature 

download data from: TERENO (Eifel-Rur), Climate/Runoff/Water Quality station Rollesbroich, Germany
(RO_AKRW_003.2007-01-01T00_00_00.2015-01-01T00_00_00.nc)
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr


def calculate_et0_hargreaves(dates, tmin, tmax, latitude_deg):
    """
    Calculates daily reference evapotranspiration (ET0) using the Hargreaves equation.
    
    Parameters:
        dates (array-like): Array of datetime objects (used to extract day of year).
        tmin (array-like): Minimum daily temperature (°C).
        tmax (array-like): Maximum daily temperature (°C).
        latitude_deg (float): Latitude in degrees (positive for North, negative for South).
        
    Returns:
        et0 (array-like): Estimated daily evapotranspiration (mm/day).
    """
    latitude_rad = np.radians(latitude_deg)
    doy = pd.to_datetime(dates).dayofyear  # Get day of year (1–365 or 366) from dates
    G_sc = 0.0820  # Solar constant [MJ m^-2 min^-1]
    dr = 1 + 0.033 * np.cos(2 * np.pi * doy / 365)  # Inverse relative distance Earth–Sun (dimensionless)
    delta = 0.409 * np.sin(2 * np.pi * doy / 365 - 1.39)  # Solar declination [radians]
    omega_s = np.arccos(-np.tan(latitude_rad) * np.tan(delta))  # Sunset hour angle [radians]
    Ra = (24 * 60 / np.pi) * G_sc * dr * (
        omega_s * np.sin(latitude_rad) * np.sin(delta) +
        np.cos(latitude_rad) * np.cos(delta) * np.sin(omega_s)
    )  # Extraterrestrial radiation (Ra) in MJ m^-2 day^-1
    tmean = (tmax + tmin) / 2  # Mean daily temperature
    et0 = 0.0023 * (Ra / 2.45) * (tmean + 17.8) * np.sqrt(tmax - tmin)  # Hargreaves equation to compute ET0 [mm/day]; 1 mm/day of evapotranspiration≈2.45 MJ/m2/day
    return et0


def estimate_et0(ds, start_date = 0, end_date = -1):
    """
    Estimate daily reference evapotranspiration (ET₀) using the Hargreaves method for each station in a dataset.

    Parameters:
        ds (xarray.Dataset): Dataset containing meteorological variables (NetCDF format, e.g. Tereno database).
        start_date (int or str): Optional index or date to start the analysis (default is 0).
        end_date (int or str): Optional index or date to end the analysis (default is -1, i.e. till the end).

    Returns:
        results (list of pandas.DataFrame): A list of DataFrames, one per station, each containing:
            - 'date': date of observation
            - 'ET0_mm_day': estimated evapotranspiration (mm/day)
            - 'station': station identifier
    """

    station_ids = ds['station'].values
    results = []

    for i, station in enumerate(station_ids):  # Loop over stations

        # Extract air temperature for current station
        temp = ds['AirTemperature_2m'].isel(station = i).to_series()
        temp.name = 'temp'
        temp = temp.dropna()

        # Extract time and make DataFrame
        df = temp.to_frame()[start_date:end_date]
        df['date'] = df.index.date
        daily_tmin = df.groupby('date')['temp'].min()
        daily_tmax = df.groupby('date')['temp'].max()

        # Get latitude for the station
        lat = ds['lat'].isel(station = i).values.item()

        # Align dates
        common_dates = daily_tmin.index.intersection(daily_tmax.index)
        tmin = daily_tmin.loc[common_dates].values
        tmax = daily_tmax.loc[common_dates].values
        dates = pd.to_datetime(common_dates)

        # Calculate ET₀
        et0 = calculate_et0_hargreaves(dates, tmin, tmax, lat)
        df_et0 = pd.DataFrame({
            'date': common_dates,
            'ET0_mm_day': et0,
            'station': station.item() if hasattr(station, "item") else station
        })

        results.append(df_et0)

    return results


def get_daily_precip(ds, start_date = 0, end_date = -1):
    """
    Calculate daily precipitation from cumulative 10-min values.

    Parameters:
    - ds: xarray.Dataset
    - start_date: optional start index for slicing
    - end_date: optional end index for slicing

    Returns:
    - List of DataFrames per station with daily precipitation (mm/day)
    """
    station_ids = ds['station'].values
    results = []

    for i, station in enumerate(station_ids):

        # Extract the cumulative precipitation series for station i
        precip = ds['PrecipitationAmount_Cum10min'].isel(station = i).to_series()
        precip.name = 'precip'
        precip = precip.dropna()

        # Check units
        ds['PrecipitationAmount_Cum10min'].attrs.get('units', 'No units attribute found')
        # print(f"Units for station {station}: {units}")

        # Calculate incremental 10-min precipitation by differencing
        precip_diff = precip.diff()
        precip_diff = precip_diff.clip(lower = 0)  # Remove negative resets (e.g., sensor reset)

        # Build DataFrame
        df = precip_diff.to_frame()[start_date:end_date]
        df['date'] = df.index.date

        # Sum increments per day
        daily_precip = df.groupby('date')['precip'].sum()

        # Store in output DataFrame
        df_out = pd.DataFrame({
            'date': daily_precip.index,
            'Precip_mm_day': daily_precip,
            'station': station.item() if hasattr(station, "item") else station
        })

        results.append(df_out)

    return results


def get_netinf(ds, start_date = 0, end_date = -1):
    """
    Calculate net infiltration for the first station in the dataset.

    Parameters:
    - ds: xarray.Dataset
        Dataset containing climate data (precipitation, temperature, etc.).
    - start_date: int, optional
        Index for start of time range (default: 0 = start of dataset).
    - end_date: int, optional
        Index for end of time range (default: -1 = end of dataset).

    Returns:
    - df: pandas.DataFrame
        DataFrame with daily date, station ID, precipitation, ET₀,
        and computed net infiltration (precip - ET₀), for one station.
    """
    et_df = estimate_et0(ds, start_date, end_date)[0]
    pr_df = get_daily_precip(ds, start_date, end_date)[0]
    et_df = et_df.reset_index(drop = True)
    pr_df = pr_df.reset_index(drop = True)
    df = pd.merge(pr_df, et_df, on = ['date', 'station'], how = 'inner')
    df['NetInfiltration_mm_day'] = df['Precip_mm_day'] - df['ET0_mm_day']
    return df


if __name__ == '__main__':

    # export net infiltration from
    # TERENO (Eifel-Rur), Climate/Runoff/Water Quality station Rollesbroich, Germany
    ds = xr.open_dataset("RO_AKRW_003.2007-01-01T00_00_00.2015-01-01T00_00_00.nc")

    # Plot
    start_date_str = '2013-05-01'
    end_date_str = '2013-08-01'

    res_et0 = estimate_et0(ds, start_date = start_date_str, end_date = end_date_str)[0]  # only one station
    print("Cumulative evapotranspiration", np.sum(res_et0['ET0_mm_day']), "mm = l/m2", len(res_et0['ET0_mm_day']))
    res_prec = get_daily_precip(ds, start_date = start_date_str, end_date = end_date_str)[0]  # only one station (start_date = '2012-01-01', end_date = '2012-03-01')
    print("Cumulative precipitation", np.sum(res_prec['Precip_mm_day']), "mm = l/m2", len(res_prec['Precip_mm_day']))
    # res_netinf = get_netinf(ds, start_date = start_date_str, end_date = end_date_str)  # only one station
    res_netinf = get_netinf(ds)  # only one station
    print("Net Infiltration", np.sum(res_netinf['NetInfiltration_mm_day']), "mm = l/m2")

    plt.figure(figsize = (12, 6))
    plt.bar(res_et0['date'], res_et0['ET0_mm_day'], color = 'skyblue', edgecolor = 'k', width = 0.8)
    plt.ylabel('ET0 (mm/day)')
    # plt.bar(res_prec['date'], res_prec['Precip_mm_day'], color = 'skyblue', edgecolor = 'k', width = 0.8)
    # plt.ylabel('Preciptiation (mm/day)')
    # plt.bar(res_netinf['date'], res_netinf['NetInfiltration_mm_day'], color = 'skyblue', edgecolor = 'k', width = 0.8)
    # plt.ylabel('Net Infiltration (mm/day)')
    plt.xlabel('Date')
    plt.show()

    print("write csv")
    res_netinf[['date', 'station', 'NetInfiltration_mm_day']].to_csv('RO_AKRW_003.2007-01-01T00_00_00.2015-01-01T00_00_00_net_infiltration.csv', index = False)
