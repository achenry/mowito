% clearvars -except switch_number 
tic

%% Eingaben

% MoWiTO files
%%%%%%%%%%%%%%%
main_data_dir = '../2021_December/Dec_9';

% load particular case results
switch switch_number % needs to be defined before running script.
    case 1
        file = fullfile(main_data_dir, 'MPC_Juan_N10_FB_240_daq');
        pitchfile = fullfile(main_data_dir, 'MPC_Juan_N10_FB_240_pitch_Block1.txt');
        blocks = 1;
        hwcalibfile = fullfile(main_data_dir, 'hw_calib_parameters_09122021_1600_2.mat'); 
        offset_env_file = fullfile(main_data_dir, 'offsets_temperature_09122021_2030.mat');
        numberhw = 1; % (1-4)
        %save_file = '..\data_store\ExperimentData.mat';
    case 2
        file = fullfile(main_data_dir, 'MPC_Juan_N15_FB_240_daq');
        pitchfile = fullfile(main_data_dir, 'MPC_Juan_N15_FB_240_pitch_Block1.txt');
        blocks = 4;
        hwcalibfile = fullfile(main_data_dir, 'hw_calib_parameters_09122021_1600_2.mat'); 
        offset_env_file = fullfile(main_data_dir, 'offsets_temperature_09122021_2030.mat');
        numberhw = 1; % (1-4)
        %save_file = '..\data_store\ExperimentData.mat';
    case 3
        file = fullfile(main_data_dir, 'MPC_Juan_N20_FB_240_daq');
        pitchfile = fullfile(main_data_dir, 'MPC_Juan_N20_FB_240_pitch_Block1.txt');
        blocks = 4;
        hwcalibfile = fullfile(main_data_dir, 'hw_calib_parameters_09122021_1600_2.mat'); 
        offset_env_file = fullfile(main_data_dir, 'offsets_temperature_09122021_2030.mat');
        numberhw = 1; % (1-4)
        %save_file = '..\data_store\ExperimentData.mat';
    case 4
        file = fullfile(main_data_dir, 'MPC_Juan_N10_FF_240_daq');
        pitchfile = fullfile(main_data_dir, 'MPC_Juan_N10_FF_240_pitch_Block1.txt');
        blocks = 2;
        hwcalibfile = fullfile(main_data_dir, 'hw_calib_parameters_09122021_1600_2.mat'); 
        offset_env_file = fullfile(main_data_dir, 'offsets_temperature_09122021_2030.mat');
        numberhw = 1; % (1-4)
        %save_file = '..\data_store\ExperimentData.mat';
    case 5
        file = fullfile(main_data_dir, 'MPC_Juan_N15_FF_240_daq');
        pitchfile = fullfile(main_data_dir, 'MPC_Juan_N15_FF_240_pitch_Block1.txt');
        blocks = 2;
        hwcalibfile = fullfile(main_data_dir, 'hw_calib_parameters_09122021_1600_2.mat'); 
        offset_env_file = fullfile(main_data_dir, 'offsets_temperature_09122021_2030.mat');
        numberhw = 1; % (1-4)
        %save_file = '..\data_store\ExperimentData.mat';
    case 6
        file = fullfile(main_data_dir, 'MPC_Juan_N20_FF_240_daq');
        pitchfile = fullfile(main_data_dir, 'MPC_Juan_N20_FF_240_pitch_Block1.txt');
        blocks = 2;
        hwcalibfile = fullfile(main_data_dir, 'hw_calib_parameters_09122021_1600_2.mat'); 
        offset_env_file = fullfile(main_data_dir, 'offsets_temperature_09122021_2030.mat');
        numberhw = 1; % (1-4)
        %save_file = '..\data_store\ExperimentData.mat';
end


%% Ende Eingaben

load(offset_env_file);
temperature = offsets_env(7,1); % in �C
ambientpressure = offsets_env(8,1); % in hPa
humidity = offsets_env(9,1); % in %
density = f_air_density2(temperature,humidity,ambientpressure*10); % in kg/m^3

% in Hz
Fs1 = 5000;
Fs2 = 100;

data = ReadBlocks(file,blocks);
pitch = load(pitchfile);

einlesen_MoWiTO_azimuthHW(data,Fs1,Fs2,1,hwcalibfile,offset_env_file,density,temperature,numberhw);

% Umrechnen Pitchfile; Spalte 2, 4, 6 sind Winkel der Blatter 1,2,3
% und Spalten 1, 3, 5 die dazugeh�rigen Zeiten; Zeitschritt ca. 10 ms aber
% nicht En�htzeit, daher separate Vektoren

pitch(:,1) = pitch(:,1) - data(1,24);
pitch(:,3) = pitch(:,3) - data(1,24);
pitch(:,5) = pitch(:,5) - data(1,24);
pitch(:,2) = pitch(:,2)./171.*30;
pitch(:,4) = pitch(:,4)./171.*30;
pitch(:,6) = pitch(:,6)./171.*30;

% filtering torque and rpm data and calculating filtered power
load(fullfile(main_data_dir,'filtpar_bw_lp_Fp30_1_Fs36_80_o55'))
load(fullfile(main_data_dir,'filtpar_bw_lp_Fp12c5_1_Fs14c5_80_o67'))
RPM_res = resample(nonan(RPM),5000,100);
RPM_res_f = filtfilt(SOS2,G2,RPM_res(:,1));
Torque_f = filtfilt(SOS,G,Torque(:,1));
Power_f = (RPM_res_f.*2*pi/60) .* Torque_f;

% Save subset of processed data so that it's easy to access later
if 0 %exist('save_file', 'var')
    trigger = data(:,15);
    save(save_file,...
         'RPM_res_f', 'hw_v1', 'pitch', 'Torque_f', 't1', 'trigger');
end
%plot_measurements

figure(1)

subplot(4,1,1)
plot(t1,Mby_1)
hold on
plot(t1,Mby_2)
hold on
plot(t1,Mby_3)
xlabel('t [s]')
ylabel('Mby [Nm]')
grid on

subplot(4,1,2)
plot(t1,smooth(Twr_FA,10))
xlabel('t [s]')
ylabel('thrust [N]')
grid on

subplot(4,1,3)
plot(t1,Torque_f)
xlabel('t [s]')
ylabel('toruqe [Nm]')
ylim([-1, 5]);
grid on

subplot(4,1,4)
plot(t1,RPM_res_f)
xlabel('t [s]')
ylabel('rpm [1/min]')
ylim([0 600]);
grid on

%subplot(3,1,3)
%plot(t2,smooth(nonan(RPM),1))
%xlabel('t [s]')
%ylabel('rpm [1/min]')
%grid on

vars = {'ambientpressure','blocks','dir','file','humidity',...
    'hwcalibfile','numberhw','offset_env_file','offsets_env'...
    'pitchfile','temperature','Welle'};
clear(vars{:},'vars')

% Aufrufen Skript LDA Einlesen und Synchronisieren
% LDA_universal_read_and_reference

toc