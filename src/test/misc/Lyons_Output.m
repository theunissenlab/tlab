function out = Lyons_Output(data,agcf)
% Returns the output pof Lyon's cochlear model on waveform data, with adaptive gain control controlled by agc (1 = yes, 0 = no)

addpath('/auto/fhome/pgill/From_Natalia/AuditoryToolbox');
addpath('/auto/fhome/pgill/From_Natalia/matlab');
addpath('/home/pgill/MATLAB/plasticity/Auditory_Toolbox');
start_dir = pwd;
%cd('/auto/fhome/pgill/From_Natalia/AuditoryToolbox');

%  Parameters
sr = 32000; % Sampling rate in Hz
df = 32; % Decimation factor is Ampsamprate / wavesamprate
low_freq = 250; % Frequency of the lowest filter
high_freq = 8000; % Frequency of the highest filter
earQ = 8; % Quality factor of each filter
stepfactor = 0.25; % 1/stepfactor is approximately the number of filters per bandwidth of one filter
differ = 1; % Differential gain control on
% agcf is input to this function
taufactor = 3; % time constant of gain control

out = LyonPassiveEar_new_mod(data,sr,df,low_freq,high_freq,earQ,stepfactor,differ,agcf,taufactor);

cd(start_dir);