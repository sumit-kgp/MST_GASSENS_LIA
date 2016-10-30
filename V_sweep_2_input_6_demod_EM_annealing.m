function data = amplitude_sweeper(port, api_level)
% AMPLITUDE_SWEEPER Perform a voltage amplitude sweep using ziDAQ's sweep module
% Perform a voltage amplitude sweep and gather demodulator data from two inputs.

clear all;
clear ziDAQ;

% Check ziDAQ's ziAutoConnect (in the Utils/ subfolder) is in the path
if exist('ziAutoConnect','file') ~= 2
    fprintf('Please configure your path using the ziDAQ function ziAddPath().\n')
    fprintf('This can be found in the API subfolder of your LabOne installation.\n');
    fprintf('On Windows this is typically:\n');
    fprintf('C:\\Program Files\\Zurich Instruments\\LabOne\\API\\MATLAB2012\\\n');
    return
end

% open a connection to a Zurich Instruments server
if exist('port','var') && exist('api_level','var')
    ziAutoConnect(port, api_level);
elseif exist('port','var')
    ziAutoConnect(port);
else
    ziAutoConnect();
end

% get device name (e.g. 'dev234')
device = ziAutoDetect();

% get the device type and its options (in order to set correct device-specific configuration)
devtype = ziDAQ('getByte',[ '/' device '/features/devtype' ] );
options = ziDAQ('getByte',[ '/' device '/features/options' ] );

fprintf('Device type ''%s'' with options ''%s''.\n',devtype,regexprep(options,'\n','|'));

format shortg;
time = clock;
time(6) = [];
chip_name = 'Chip15apr1';
device_name = 'EM42R';
R_s = 50; %[Ohm]
measurement_name = sprintf('%d-%02d-%02d_%02d-%02d_Psweep_%s_%s',time,chip_name,device_name);
measurement_name_plot = sprintf('%d-%02d-%02d_%02d-%02d_%s_%s',time,chip_name,device_name);
mkdir(chip_name,measurement_name);
folder_path = sprintf('%s/%s',chip_name,measurement_name);

sigout_frequency = 1333; %[Hz]
out_c1_range = 10; % signal output range [V]
start_amplitude = 0.1; %Factor of the range!
stop_amplitude = 0.45; %Factor of the range! %burning around 5.5
sweep_samplecount = 40;
bandwidth = 3; %[Hz]
demod_idx = [1,2,3,4,5,6]; % 1-based indexing
data = run_amplitude_sweep(device,devtype,options,sigout_frequency,start_amplitude,stop_amplitude,sweep_samplecount,bandwidth,demod_idx,out_c1_range);
%test_example(data,sigout_amplitude,sweep_samplecount);

% Process and plot the acquired data
table = process_data(demod_idx,data,out_c1_range,folder_path,measurement_name_plot);

% Wrtie data to file
write_data(folder_path,measurement_name,R_s,sigout_frequency,start_amplitude,stop_amplitude,sweep_samplecount,bandwidth,table,out_c1_range);

% Plot
plot_data_V_sweep(table,R_s,folder_path,measurement_name_plot);

end

function data = run_amplitude_sweep(device,devtype,options,sigout_frequency,start_amplitude,stop_amplitude,sweep_samplecount,bandwidth,demod_idx,out_c1_range)

out_c1 = '0'; % signal output channel
in_c1 = '0'; % signal input channel
in_c2 = '1'; % signal input channel
osc_c1 = '0'; % oscillator
for d = 1:length(demod_idx)
  demod_c(d) = num2str(demod_idx(d)-1,'%d');
end

in_c1_range = 2.0; %Vpk [V]
in_c2_range = 1.2; %Vpk [V]
tc = 0.07; % [s]
demod_rate = 13e3;

% define the output mixer channel based on the device type and its options
if strfind(devtype,'UHF') & isempty(strfind(options,'MF'))
    out_mixer_c = '3';
elseif strfind(devtype,'HF2') & isempty(strfind(options,'MF'))
    out_mixer_c = '6';
else
% instrument with multi-frequency option
    out_mixer_c = '0';
end

% create a base configuration: disable all outputs, demods and scopes
ziDAQ('setDouble',['/' device '/demods/*/rate'], 0.0);
ziDAQ('setInt',['/' device '/demods/*/trigger'], 0);
ziDAQ('setInt',['/' device '/sigouts/*/enables/*'], 0);
if strfind(devtype,'UHF')
% if the device is a UHF additionally disable all demodulators
    ziDAQ('setInt',['/' device '/demods/*/enable'], 0);
    ziDAQ('setInt',['/' device '/scopes/*/enable'], 0);
elseif strfind(devtype,'HF2')
    ziDAQ('setInt',['/' device '/scopes/*/trigchannel'],-1)
end

%% configure the device ready for this experiment
% Input channel 1
% Select input impedance
ziDAQ('setInt',['/' device '/sigins/' in_c1 '/imp50'], 1);
% Select if AC coupling is OFF(0) or ON(1)
ziDAQ('setInt',['/' device '/sigins/' in_c1 '/ac'], 0);
% Set input range [Vpk]
ziDAQ('setDouble',['/' device '/sigins/' in_c1 '/range'], in_c1_range);

% Input channel 2
% Select input impedance
ziDAQ('setInt',['/' device '/sigins/' in_c2 '/imp50'], 0);
% Select if AC coupling is OFF(0) or ON(1)
ziDAQ('setInt',['/' device '/sigins/' in_c2 '/ac'], 0);
% Set input range [Vpk]
ziDAQ('setDouble',['/' device '/sigins/' in_c2 '/range'], in_c2_range);

% Output channel 1
% Select the range of the Signal Output
ziDAQ('setDouble',['/' device '/sigouts/' out_c1 '/range'], out_c1_range);
% Set the Signal Output Add, Out, Sync amplitudes to zero
ziDAQ('setDouble',['/' device '/sigouts/' out_c1 '/amplitudes/*'], 0);
% Set the Signal Output amplitude to sigout_amplitude?
ziDAQ('setDouble',['/' device '/sigouts/' out_c1 '/amplitudes/' out_mixer_c], 0);
% Enable the Signal Output?
ziDAQ('setDouble',['/' device '/sigouts/' out_c1 '/enables/' out_mixer_c], 1);
% Set the Signal Output OFF(0) or ON(1)
ziDAQ('setInt',['/' device '/sigouts/' out_c1 '/on'], 1);

if strfind(devtype,'HF2')
    % Select if the Signal Input of channel 1 is differential (1) or not (0)
    ziDAQ('setInt',['/' device '/sigins/' in_c1 '/diff'], 0);
    % Select if the Signal Input of channel 2 is differential (1) or not (0)
    ziDAQ('setInt',['/' device '/sigins/' in_c2 '/diff'], 1);
    % Select if the Add port of Signal Output of channel 1 is enabled (1) or not (0)
    ziDAQ('setInt',['/' device '/sigouts/' out_c1 '/add'], 0);
end

ziDAQ('setInt',['/' device '/demods/*/oscselect'], 0);
ziDAQ('setDouble',['/' device '/demods/*/phaseshift'], 0);
ziDAQ('setInt',['/' device '/demods/*/order'], 4);
for d = demod_c
  ziDAQ('setDouble',['/' device '/demods/' d '/rate'], demod_rate);
  if strfind(devtype,'UHF')
    ziDAQ('setInt',['/' device '/demods/' d '/enable'], 1);
  end
end
ziDAQ('setInt',['/' device '/demods/' demod_c(1) '/harmonic'], 1);
ziDAQ('setInt',['/' device '/demods/' demod_c(2) '/harmonic'], 2);
ziDAQ('setInt',['/' device '/demods/' demod_c(3) '/harmonic'], 3);
ziDAQ('setInt',['/' device '/demods/' demod_c(4) '/harmonic'], 1);
ziDAQ('setInt',['/' device '/demods/' demod_c(5) '/harmonic'], 2);
ziDAQ('setInt',['/' device '/demods/' demod_c(6) '/harmonic'], 3);

if strfind(options,'MF')
% HF2IS and HF2LI multi-frequency option do not support the node oscselect.
    ziDAQ('setInt',['/' device '/demods/*/oscselect'], str2double(osc_c1));
    ziDAQ('setInt',['/' device '/demods/*/adcselect'], str2double(in_c1));
    ziDAQ('setInt',['/' device '/demods/*/adcselect'], str2double(in_c2));
end

ziDAQ('setDouble',['/' device '/demods/*/timeconstant'], tc);
ziDAQ('setDouble',['/' device '/oscs/' osc_c1 '/freq'], sigout_frequency); % [Hz]

%% Sweeper settings
% Create a thread for the sweeper 
timeout = 500 % milliseconds
h = ziDAQ('sweep', timeout);
% Device on which sweeping will be performed
ziDAQ('set', h, 'sweep/device', device);
% Sweeping setting is the amplitude of the output signal
ziDAQ('set', h, 'sweep/gridnode', ['sigouts/' out_c1 '/amplitudes/' out_mixer_c]);
%ziDAQ('set', h, 'sweep/gridnode', ['sigouts/' out_c '/amplitutes/']);
% Set start amplitude
ziDAQ('set', h, 'sweep/start', start_amplitude);
% Set stop amplitude
ziDAQ('set', h, 'sweep/stop', stop_amplitude);
% sweep_samplecount measurement points (for sweep_samplecount different
% amplitudes input signal Parameters will be recorded)
ziDAQ('set', h, 'sweep/samplecount', sweep_samplecount);
% Single sweep 
ziDAQ('set', h, 'sweep/loopcount', 1);
% Logarithmic sweep mode
ziDAQ('set', h, 'sweep/xmapping', 0);
% Binary scan type -1, Bidirectional -2, Unidirectional - 0
ziDAQ('set', h, 'sweep/scan', 2);
% The is no user system involved in the setup, so it's needed to wait only
% filter response to settle 
ziDAQ('set', h, 'sweep/settling/time', 0);
% Wait 15 time constants for filter response to settle
ziDAQ('set', h, 'sweep/settling/tc', 10);
% Minimum time to record data is 50 time constants  
ziDAQ('set', h, 'sweep/averaging/tc', 10);
% Minimal number of samples that we want to record is 100
ziDAQ('set', h, 'sweep/averaging/sample', 5);
% Bandwidth control for each measurement
% Automatic:
% ziDAQ('set', h, 'sweep/bandwidthcontrol', 2);
% Fixed: set sweep/bandwidthcontrol to 1 and specify a bandwidth via 'sweep/bandwidth'
ziDAQ('set', h, 'sweep/bandwidthcontrol', 1);
ziDAQ('set', h, 'sweep/bandwidth', bandwidth);
% Manual: set  sweep/bandwidthcontrol to 2. sweep/bandwidth must also be set to a value > 0 although it is ignored. 
% Otherwise Auto control is automatically chosen (for backwards compatibility reasons).
% ziDAQ('set', h, 'sweep/bandwidthcontrol', 2);
% ziDAQ('set', h, 'sweep/bandwidth', bandwidth);

% Subscribe to the node from which data will be recorded
for d = demod_c
  ziDAQ('subscribe', h, ['/' device '/demods/' d '/sample']);
end

% Start sweeping
ziDAQ('execute', h);

data = [];

% Waiting untill sweeping is not finished
while ~ziDAQ('finished', h)
    pause(1);
    fprintf('Sweep progress %0.0f%%\n', ziDAQ('progress', h) * 100);
end

% Set the Signal Output OFF(0)
ziDAQ('setInt',['/' device '/sigouts/' out_c1 '/on'], 0);

% now read the data. This command can also be executed during the waiting.
tmp = ziDAQ('read', h);

% unsubscribe from the node; stop filling the data from that node to the internal buffer in the server
ziDAQ('unsubscribe', h, ['/' device '/demods/*/sample']);

% aqcuire the read data
if isfield(tmp,device)
    tmp = tmp.(device);
    if isfield(tmp,'demods')
        if ~isempty(tmp.demods(demod_idx(1)).sample)
            if ~isempty(tmp.demods(demod_idx(2)).sample)
                if ~isempty(tmp.demods(demod_idx(3)).sample)
                    if ~isempty(tmp.demods(demod_idx(4)).sample)
                        if ~isempty(tmp.demods(demod_idx(5)).sample)
                            if ~isempty(tmp.demods(demod_idx(6)).sample)
                                % As several sweeps may be returned, a cell array is used.          
                                data = tmp;
                            end
                        end
                    end
                end
            end
        end
    end
end

end

function table = process_data(demod_idx,data,out_c1_range,folder_path,measurement_name_plot)
% Process the acquired data
% As several sweeps may be returned, a cell array is used.
% In this case we pick the first sweep result by {1}.    
for d = demod_idx
    mydata(d) = data.demods(d).sample{1};
end
for d = demod_idx
    r(:,d) = mydata(d).r;
end
for d = demod_idx
    theta(:,d) = mydata(d).phase;
end
% Amplitude values at which measurement points were taken
for d = demod_idx
    amplitudes(:,d) = (mydata(d).grid)*out_c1_range;
end
% if isequal(amplitudes(:,demod_idx(1)),amplitudes(:,demod_idx(2))) 
    amplitudes = amplitudes(:,demod_idx(1));
% else
%     error('Error. The amplitude sweep values do not match.')
% end
% Plot the final result
close all;
for d = demod_idx
    plot_data(d, amplitudes, r(:,d), theta(:,d), '-',folder_path,measurement_name_plot)
end

% Put the data in a table
table = [amplitudes r theta];
table( :, ~any(table,1) ) = [];

end

function plot_data(demodulator, amplitudes, r, theta, style,folder_path,measurement_name_plot)
% Plot data
paper_size = [16 12];
figure(); clf;
title_str = sprintf('Demodulator %d',demodulator);
subplot(2,1,1)
s = plot(amplitudes, r, style);
set(s, 'LineWidth', 1.5)
set(s, 'Color', 'black');
grid on
title(title_str)
xlabel('Vs [Vpk]')
ylabel('r [Vrms]')
subplot(2,1,2)
s = plot(amplitudes, theta*180/pi, style);
set(s, 'LineWidth', 1.5)
set(s, 'Color', 'black');
grid on
xlabel('V_s [V]')
ylabel('Phase [deg]')
figure_name = sprintf('%s/%s_demod%d',folder_path,measurement_name_plot,demodulator);
set(gcf,'paperunits','centimeters','papersize',paper_size,'paperposition',[0,0,paper_size(1),paper_size(2)])
print('-dpdf',figure_name)
print('-dpng',figure_name)
end

function plot_data_V_sweep(data,R_s,folder_path,measurement_name_plot)
% Plot data
paper_size = [8 6];

deltaT_func = inline('(R_em-R_0)./(2e-3*R_0)','R_em','R_0');
I = data(:,2)./R_s;          %[mA]
R_em = (data(:,3))./I;       %[Ohm]
P_em = (data(:,3)).*I;       %[mW]
% % De-amplifier at voltage across EM:
% R_em = ((data(:,3)).*10)./I;       %[Ohm]
% P_em = ((data(:,3)).*10).*I;       %[mW]

figure()
plot(I,R_em,'.','MarkerSize',13)
xlabel('Current [A]','FontSize',12)
ylabel('Resistance [\Omega]','FontSize',12)
figure_name = sprintf('%s/%s_R_vs_I',folder_path,measurement_name_plot);
set(gcf,'paperunits','centimeters','papersize',paper_size,'paperposition',[0,0,paper_size(1),paper_size(2)])
print('-dpdf',figure_name)
print('-dpng',figure_name)

figure()
plot(P_em,R_em,'.','MarkerSize',13)
xlabel('Power [W]','FontSize',12)
ylabel('Resistance [\Omega]','FontSize',12)
figure_name = sprintf('%s/%s_R_vs_P',folder_path,measurement_name_plot);
set(gcf,'paperunits','centimeters','papersize',paper_size,'paperposition',[0,0,paper_size(1),paper_size(2)])
print('-dpdf',figure_name)
print('-dpng',figure_name)
end

function write_data(folder_path,measurement_name,R_s,sigout_frequency,start_amplitude,stop_amplitude,sweep_samplecount,bandwidth,table,out_c1_range)
% Write data to file
file_path_name = sprintf('%s/%s.txt',folder_path,measurement_name);
fileID = fopen(file_path_name,'w');
fprintf(fileID,'%%%s\n',measurement_name);
fprintf(fileID,'%%No thermal release tape between chip and chuck. Kapton tape fix.\n');
% fprintf(fileID,'%%De-amplifier at voltage across EM.\n');
fprintf(fileID,'%%R_s = %d Ohm (Input impedance)\n',R_s);
fprintf(fileID,'%%Frequency of driving signal = %.2f Hz\n',sigout_frequency);
fprintf(fileID,'%%Start amplitude = %.3f V \n%%Stop amplitude = %.3f V \n%%Measurement points = %d\n%%Bandwidth = %.3f Hz\n',...
    start_amplitude*out_c1_range,stop_amplitude*out_c1_range,sweep_samplecount,bandwidth);
fprintf(fileID,'%%\n%%%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
    'Vs [V]','r_Rs1 [V_rms]','r_Rs2 [V_rms]','r_Rs3 [V_rms]','r_h1 [V_rms]','r_h2 [V_rms]','r_h3 [V_rms]',...
    'theta_Rs1 [°]', 'theta_Rs2 [°]', 'theta_Rs3 [°]','theta_h1 [°]','theta_h2 [°]','theta_h3 [°]');
dlmwrite(file_path_name,table,'-append','delimiter','\t');
fclose(fileID);
end

function test_data(data,amplitude,sweep_samplecount)
% TEST_EXAMPLE run some simple tests on the data returned by the example
assert(iscell(data.demods.sample),'Assertion failed: iscell(data.demods.sample).');
assert(length(data.demods.sample) == 1,'Assertion failed: length(data.demods.sample) == 1.');
assert(length(data.demods.sample{1}.r) == sweep_samplecount,'Assertion failed: length(data.demods.sample) == sweep_samplecount.');
end

% Local variables:
% matlab-indent-level: 4
% matlab-indent-function-body: nil
% End:
