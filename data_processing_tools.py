import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
from os.path import exists as file_exists
from scipy.signal import sosfiltfiltv


def load_raw_data(data_file_main):
    """
    data_file_main should be everthing up to _daq_BlockX.txt. All blocks will be 
    loaded, as well as the relevant pitch file.
    """
    print('Beginning to load raw data...')

    df_raw_daq = pd.read_csv(data_file_main+'_daq_Block1.txt', delimiter='\t', 
        header=None, names=raw_data_columns())
    i = 2
    while file_exists(data_file_main+'_daq_Block'+str(i)+'.txt'):
        df_raw_daq = pd.concat((df_raw_daq, pd.read_csv(
            data_file_main+'_daq_Block'+str(i)+'.txt', delimiter='\t', 
            header=None, names=raw_data_columns())),
            axis=0, ignore_index=True
        )
        i = i + 1

    df_raw_pitch = pd.read_csv(data_file_main+'_pitch_Block1.txt', 
        delimiter='\t', header=None, names=pitch_data_columns())
        
    print('Raw data loaded.')
    return df_raw_daq, df_raw_pitch

def process_raw_data(df_raw_daq, df_raw_pitch, offset_env_file, hw_calib_file):
    """
    Based on einlesen_MoWiTO_azimuthHW.m.
    """
    # Starting clock time for the data set
    t0 = df_raw_daq.epoch_time[0]
    Fs1 = 5000 # Raw data sampling frequency [Hz]

    # Deal with blade pitch data
    df_pitch = pd.DataFrame()
    for b in range(1,4):
        # What exactly is this conversion factor (30/171)?
        df_pitch['pitch_'+str(b)] = df_raw_pitch['pitch_'+str(b)]/171*30 
        df_pitch['time_'+str(b)] = df_raw_pitch['time_'+str(b)]-t0
    
    # Load offsets data
    with open(offset_env_file+'.csv', newline='') as f:
        reader = csv.reader(f)
        offsets_env = [float(row[0]) for row in reader]
        
    off_Mby1 = offsets_env[0]
    off_Mby2 = offsets_env[1]
    off_Mby3 = offsets_env[2]
    off_TwrFA = offsets_env[3]
    off_TwrSS = offsets_env[4]
    off_Prandtl = offsets_env[5]
    
    temperature = offsets_env[6] # deg C
    ambientpressure = offsets_env[7] # hPa
    humidity = offsets_env[8] # as percentage
    density = air_density(temperature, humidity, ambientpressure*10)
    
    # Build up processed dataframe
    df_daq = pd.DataFrame()
    df_daq['Torque'] = 2*df_raw_daq.U_Torque
    df_daq['Mby_1'] = (df_raw_daq.U_Mby_1 - off_Mby1) / -0.2448
    df_daq['Mby_2'] = (df_raw_daq.U_Mby_2 - off_Mby2) / -0.2722
    df_daq['Mby_3'] = (df_raw_daq.U_Mby_3 - off_Mby3) / -0.2586
    df_daq['Welle'] = df_raw_daq.U_Welle
    df_daq['Twr_FA'] = (df_raw_daq.U_Twr_FA - off_TwrFA) / -0.0337
    df_daq['Twr_SS'] = (df_raw_daq.U_Twr_SS - off_TwrSS) / -0.0263

    # Apply calibration curve to raw hotwire data
    # (assuming taking the mean over the coefficients is appropriate)
    with open(hw_calib_file+'.csv', newline='') as f: 
        reader = csv.reader(f) 
        hw_calib_data = [row for row in reader] 
    p_HW1 = np.flip(np.array(hw_calib_data, float).mean(axis=0))
    df_daq['hw_v1'] = p_HW1[0] + p_HW1[1]*df_raw_daq.HW1 + \
        p_HW1[2]*df_raw_daq.HW1**2 + p_HW1[3]*df_raw_daq.HW1**3 + \
        p_HW1[4]*df_raw_daq.HW1**4
    df_daq['RPM'], df_daq['azimuth'] = construct_RPM_and_azimuth(
        df_raw_daq.count_.to_numpy(copy=True), 
        df_raw_daq.countzero.to_numpy(copy=True), 
        df_raw_daq.saved_data.to_numpy(copy=True), Fs1)

    # Power, thrust parts
    diameter = 1.8 # m
    PowerWind = 0.5 * density * df_daq.hw_v1**3 * np.pi * (diameter/2)**2
    ThrustRef = 0.5 * density * df_daq.hw_v1**2 * np.pi * (diameter/2)**2
    df_daq['omega'] = df_daq.RPM * 2*np.pi/60
    df_daq['Power'] = df_daq.omega * df_daq.Torque
    df_daq['Cp'] = df_daq.Power/PowerWind
    df_daq['Ct'] = df_daq.Twr_FA/ThrustRef

    # Time
    df_daq['time'] = np.linspace(0, len(df_raw_daq)/Fs1, len(df_raw_daq), 
        endpoint=False)

    # High-level channels
    df_daq['saved_data'] = df_raw_daq.saved_data
    
    print('Raw data processed.')
    return df_daq, df_pitch

def downsample_data(df_daq, df_pitch, apply_filtering=True):
    """
    Downsample from high-frequency data (5000 Hz) to low frequency (100 Hz). 
    In doing so, extract the saved_data, and filter the rotor speed and 
    torque if apply_filtering is true.

    Note that resampler.mean() ignores NaNs in mean, as desired.
    """
    Fs1 = 5000 # Raw data sampling frequency [Hz]
    Fs2 = 100 # Higher-level data sampling freqency (to downsample to) [Hz]
    
    df_ds = pd.DataFrame()
    df_ds['time'] = np.linspace(0, len(df_daq)/Fs1, 
        round(len(df_daq)/(Fs1/Fs2)), endpoint=False)

    # Set up appropriate time index for resampling
    df_daq.time = pd.to_timedelta(df_daq.time, unit='S')
    df_daq.set_index('time', inplace=True)
    
    # Blade pitch angles
    for b in range(1,4):
        pitch = df_pitch['pitch_'+str(b)].values
        if len(df_ds) >= len(pitch): # pad pitch with NaNs
            df_ds['pitch_'+str(b)] = np.concatenate((pitch, 
                np.nan*np.ones(len(df_ds) - len(pitch))))
        else: # truncate pitch
            df_ds['pitch_'+str(b)] = pitch[:len(df_ds)]
            
    
    # Filtered signals
    if apply_filtering:
        use_saved_gains = False # Was creating numerical issues
        with open('filtpar_bw_lp_Fp30_1_Fs36_80_o55_SOS.csv', newline='') as f: 
            reader = csv.reader(f) 
            SOS = [row for row in reader]
        SOS = np.array(SOS, float)
        with open('filtpar_bw_lp_Fp12c5_1_Fs14c5_80_o67_SOS.csv', newline='') \
            as f: 
            reader = csv.reader(f) 
            SOS2 = [row for row in reader]
        SOS2 = np.array(SOS2, float)
        if use_saved_gains:
            with open('filtpar_bw_lp_Fp30_1_Fs36_80_o55_G.csv', newline='') \
                as f: 
                reader = csv.reader(f) 
                G = [row for row in reader]
            G = np.array(G, float)
            with open('filtpar_bw_lp_Fp12c5_1_Fs14c5_80_o67_G.csv', newline='')\
                as f: 
                reader = csv.reader(f) 
                G2 = [row for row in reader]
            G2 = np.array(G2, float)
        else: # Compute gains directly
            G = 1/(np.sum(SOS[:,0:3], axis=1)/np.sum(SOS[:,3:6], axis=1))
            G2 = 1/(np.sum(SOS2[:,0:3], axis=1)/np.sum(SOS2[:,3:6], axis=1))

        RPM_uf = df_daq.RPM.interpolate(method='linear', limit_direction='both')
        RPM_f = sosfiltfilt(SOS2, RPM_uf.values) * np.prod(G2)**2
        
        Torque_uf = df_daq.Torque.interpolate(method='linear', 
            limit_direction='both')
        Torque_f = sosfiltfilt(SOS, Torque_uf.values) * np.prod(G)**2
        Power_f = RPM_f * (2*np.pi/60) * Torque_f

        df_ds['RPM'] = pd.Series(RPM_f, index=df_daq.index).\
            resample(str(1/Fs2)+'S').mean().tolist()
        df_ds['Torque'] = pd.Series(Torque_f, index=df_daq.index).\
            resample(str(1/Fs2)+'S').mean().tolist()
        df_ds['Power'] = pd.Series(Power_f, index=df_daq.index).\
            resample(str(1/Fs2)+'S').mean().tolist()
    else: # Don't filter
        df_ds['RPM'] = df_daq.RPM.resample(str(1/Fs2)+'S').mean().tolist() 
        df_ds['Torque'] = df_daq.Torque.resample(str(1/Fs2)+'S').mean().tolist()
        df_ds['Power'] = df_daq.Power.resample(str(1/Fs2)+'S').mean().tolist()

    # Take mean over most of the other columns
    df_ds[resampling_columns()] = df_daq[resampling_columns()].\
        resample(str(1/Fs2)+'S').mean().to_numpy()

    # Extract special channel data. Just name s0,...,sX for now.
    S = df_daq.saved_data.to_numpy(copy=True).reshape(len(df_ds),round(Fs1/Fs2))
    n_keep = 20
    S = S[:,:n_keep]
    labels = ['s'+str(s) for s in range(n_keep)]
    df_ds[labels] = S

    print('Data downsampled.')

    return df_ds

def save_test_data_to_pickle(df_ds, save_path=None, time_range=None, 
    rename_channels=None, save_channels=None):
    """
    Keep only some channels of interest and the relevant test times.
    time_range should be a two-tuple of (start_time, stop_time) (inclusive)
    rename_channels should be a dictionary mapping original channel names to 
        new ones.
    save_channels should be a list of the channels which are to be kept (in the
        new names specified by rename channels).

    If save_path is None, data will not be saved.
    """

    if time_range is not None:
        df_ds = df_ds[(df_ds.time >= time_range[0]) & 
                      (df_ds.time <= time_range[1])]
    
    if rename_channels is not None:
        df_ds.rename(columns=rename_channels, inplace=True)

    if save_channels is not None:
        df_ds = df_ds[save_channels]

    if save_path is not None: 
        df_ds.to_pickle(save_path)
        print('Data written to pickle.')
    else:
        print('Data not saved.')

    return df_ds

def run_all_steps_mpc(data_file, offset_file, hw_calib_file, save_file=None, 
    time_range=None, plot_time_range=True):
    """
    Run through each step.
    """
    df_raw_daq, df_raw_pitch = load_raw_data(data_file)
    df_daq, df_pitch = process_raw_data(df_raw_daq, df_raw_pitch, offset_file,
        hw_calib_file)
    df_ds = downsample_data(df_daq, df_pitch, apply_filtering=True)
    if plot_time_range:
        plot_time_window(df_ds, time_range=time_range,  title=data_file)
    df_mpc = save_test_data_to_pickle(df_ds, save_path=save_file, 
        rename_channels=mpc_columns())

    print('All steps run.')

    return df_mpc

def plot_time_window(df, time_range=None, title=None):
    """
    Plot wind speed, blade pitch angle, and rotor speed across all times in
    df, but indicate window selected in time_range.
    """

    fig, ax = plt.subplots(3,1,sharex=True)
    fig.set_size_inches(10, 10)

    # y axis limits for plots
    WSlims = [0, 15]
    BPlims = [-5, 20]
    RSlims = [400, 600]

    ax[0].plot(df.time.values, df.hw_v1.values, color='C0', zorder=1)
    ax[1].plot(df.time.values, df.pitch_1.values, color='red', zorder=1)
    ax[2].plot(df.time.values, df.RPM.values, color='black', zorder=1)

    if time_range is not None:
        fill_color = 'gray'
        fill_alpha = 0.5
        for a, l in zip(ax, [WSlims, BPlims, RSlims]):
            a.fill_between([df.time.iloc[0], time_range[0]], 
                           [l[1], l[1]], [l[0], l[0]],
                           facecolor=fill_color, edgecolor=None, 
                           alpha=fill_alpha, zorder=2)
            a.fill_between([time_range[1], df.time.iloc[-1]], 
                           [l[1], l[1]], [l[0], l[0]],
                           facecolor=fill_color, edgecolor=None, 
                           alpha=fill_alpha, zorder=2)

    # Plot aesthetics
    ax[0].set_title(title)
    ax[0].set_ylim(WSlims)
    ax[0].set_ylabel('HW wind speed [m/s]')
    ax[1].set_ylim(BPlims)
    ax[1].set_ylabel('Blade 1 pitch [deg]')
    ax[2].set_ylim(RSlims)
    ax[2].set_ylabel('Rotor speed [rpm]')
    ax[2].set_xlim([df.time.iloc[0], df.time.iloc[-1]])
    ax[2].set_xlabel('Time [s]')

    return ax

def construct_RPM_and_azimuth(count, countzero, time, Fs1):
    """
    Compute rotor/generator speed in rpm and azimuth position from raw data.
    Inputs should be 1D numpy arrays.
    """
    count_mask = (count == 9999)
    count[count_mask] = np.nan
    time[count_mask] = np.nan

    RPM = (count/time) * (1000/1440 * 60) # Not sure about unit conversion here.

    omega2 = pd.Series(RPM*(2*np.pi/60)).interpolate(
        method='linear', limit_direction='both')
    omega2 = omega2.values*(360/(Fs1*2*np.pi))
    
    azimuth = np.zeros(len(countzero))

    for i in range(1,len(countzero)):
        if countzero[i] == 1 and countzero[i-1] == 0:
            pass # azimuth already zero here. 
        else:
            azimuth[i] = azimuth[i-1]+omega2[i]

    azimuth = azimuth - 86
    azimuth[azimuth < 0] = azimuth[azimuth < 0] + 360
    
    return RPM, azimuth

def air_density(t, hr, p):
    """
    AIR_DENSITY calculates density of air
    %  Usage :[ro] = air_density(t,hr,p)
    %  Inputs:   t = ambient temperature (ºC)
    %           hr = relative humidity [%]
    %            p = ambient pressure [Pa]  (1000 mb = 1e5 Pa)
    %  Output:  ro = air density [kg/m3]

    %
    %  Refs:
    % 1)'Equation for the Determination of the Density of Moist Air' P. Giacomo 
    %    Metrologia 18, 33-40 (1982)
    % 2)'Equation for the Determination of the Density of Moist Air' R. S. Davis
    %    Metrologia 29, 67-70 (1992)
    %
    % ver 1.0   06/10/2006    Jose Luis Prego Borges (Sensor & System Group, 
    %    Universitat Politecnica de Catalunya)
    % ver 1.1   05-Feb-2007   Richard Signell (rsignell@usgs.gov)  Vectorized 
    """
    # -------------------------------------------------------------------------
    T0 = 273.16        # Triple point of water (aprox. 0C)
    T = T0 + t         # Ambient temperature in Kelvin

    #-------------------------------------------------------------------------
    #-------------------------------------------------------------------------
    # 1) Coefficients values

    R =  8.314510            # Molar ideal gas constant   [J/(mol.ºK)]
    Mv = 18.015e-3        # Molar mass of water vapour [kg/mol]
    Ma = 28.9635e-3       # Molar mass of dry air      [kg/mol]

    A =  1.2378847e-5    # [K^-2]
    B = -1.9121316e-2    # [K^-1]
    C = 33.93711047         #
    D = -6.3431645e3     # [K]
    
    a0 =  1.58123e-6;      # [K/Pa]
    a1 = -2.9331e-8;       # [1/Pa]
    a2 =  1.1043e-10;      # [1/(K.Pa)]
    b0 =  5.707e-6;        # [K/Pa]
    b1 = -2.051e-8;        # [1/Pa]
    c0 =  1.9898e-4;       # [K/Pa]
    c1 = -2.376e-6;        # [1/Pa]
    d =  1.83e-11;         # [K^2/Pa^2]
    e = -0.765e-8;         # [K^2/Pa^2]

    #-------------------------------------------------------------------------
    # 2) Calculation of the saturation vapour pressure at ambient temperature, 
    # in [Pa]
    psv = np.exp(A*(T**2) + B*T + C + D/T)  # [Pa]


    #-------------------------------------------------------------------------
    # 3) Calculation of the enhancement factor at ambient temperature and 
    # pressure
    fpt = 1.00062 + (3.14e-8)*p + (5.6e-7)*(t**2)


    #-------------------------------------------------------------------------
    # 4) Calculation of the mole fraction of water vapour
    xv = hr*fpt*psv*(1/p)*(1e-2)


    #-------------------------------------------------------------------------
    # 5) Calculation of the compressibility factor of air
    Z = 1 - ((p/T)*(a0 + a1*t + a2*(t**2) + (b0+b1*t)*xv + (c0+c1*t)*(xv**2)))+\
        ((p**2/T**2)*(d + e*(xv**2)))


    #-------------------------------------------------------------------------
    # 6) Final calculation of the air density in [kg/m^3]
    ro = (p*Ma/(Z*R*T))*(1 - xv*(1-Mv/Ma))

    return ro

def raw_data_columns():
    """
    Return list of column names for the raw daq data.
    """

    columns = [
        'U_Mby_1',
        'U_Mby_2',
        'U_Mby_3',
        'U_Welle',
        'unknown_4',
        'U_Twr_FA',
        'U_Twr_SS',
        'U_Torque',
        'unknown_8',
        'unknown_9',
        'unknown_10',
        'HW1',
        'unknown_12',
        'HW2',
        'HW3',
        'unknown_15',
        'unknown_16',
        'unknown_17',
        'unknown_18',
        'unknown_19',
        'count_',
        'countzero',
        'saved_data',
        'epoch_time'
    ]

    return columns

def pitch_data_columns():
    """
    Return list of column names for the raw pitch data.
    """

    columns = [
        'time_1',
        'pitch_1',
        'time_2',
        'pitch_2',
        'time_3',
        'pitch_3',
    ]
    
    return columns

def resampling_columns():
    """
    Return list of columns for direct mean resampling in downsample_data.
    """

    columns = [
        'Mby_1',
        'Mby_2',
        'Mby_3',
        'Welle',
        'Twr_FA',
        'Twr_SS',
        'hw_v1',
    ]

    return columns

def mpc_columns():
    """
    Return dict to assign special data columns.
    """

    col_names = {
        's1':'control_code',
        's10':'pitch_command',
        's11':'HW_wind_speed',
        's12':'feedback_only',
        's13':'solve_status',
        's14':'solve_time_ms',
        's15':'iterations',
    }

    return col_names

if __name__ == "__main__":
    # Perform example data processing.
    example_data = 'MPC_Juan_N20_FB_240'
    main_data_dir = '/Volumes/MS_HD2_WD/MoWiTO_testdata/2021_December/Dec_9'
    data_file = main_data_dir+'/'+example_data
    offset_file = main_data_dir+'/offsets_temperature_09122021_2030'
    hw_calib_file = main_data_dir+'/hw_calib_parameters_09122021_1600_2'

    df_raw_daq, df_raw_pitch = load_raw_data(data_file)
    df_daq, df_pitch = process_raw_data(df_raw_daq, df_raw_pitch, offset_file,
        hw_calib_file)
    df_ds = downsample_data(df_daq, df_pitch, apply_filtering=True)

    rough_time_range = [645., 745.]
    plot_time_window(df_ds, time_range=rough_time_range,  title=example_data)
    
    df_mpc = save_test_data_to_pickle(df_ds, save_path=None, 
        rename_channels=mpc_columns())

    print(df_mpc)
    
    plt.show()
    # Extract the relevant data
    
