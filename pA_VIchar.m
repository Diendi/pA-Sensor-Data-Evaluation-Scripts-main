%pA UI characteristic plotting
%reads a test parameter file (applied U-t), and reads a given pA sensor
%data file in chunks according to the timings given in the former. Raw data
%chunks are then plotted individually, filtered and plotted again, and 
%a U/I pair at the end of each voltage step is taken to plot a UI
%characteristic over the whole test duration.
%
%Matthias Musil
%FH Wels RING
%21.02.2024

clear;
close all;

warning('off','all')

%% setup
disp('Setup...')
global figcnt;

figcnt = 0;

% test number
test_nr = 125;

fft = 'fft'; % 'fft' to plot fft, '' to not plot fft
fft_lim = [0 500]; %Hz, FFT frequency limits, use '' for no limits

step_duration_known = 1;        %if step duration is given by test parameter file and known, set this to 1
                                %otherwise, the full sensor log file is
                                %loaded and iterated through first to find
                                %timings of individual steps

fn_test = 'vitestlongv2.csv';       %test parameter file
%folder_sensor = 'sensordaten/';
folder_sensor = '';
fn_sensor = strcat('test_',num2str(test_nr),'.csv');   %sensor log file
test_name = strtok(fn_sensor,'.'); %test_name to store figures, removes .csv from filename
test_name = strcat(test_name,''); %optional renaming for testing of different filters
test_name_tex = strcat('test\_',num2str(test_nr));
fn_mpd = strcat(num2str(test_nr),'/dcFile1.csv'); %MPD pd measurement file

fs = 12.5e3;        %sensor sampling rate

%depecrated - not used anymore
%n_filt = 1;       %number of averages for movmean filter
%f_lp = 10;          %lowpass fc
%filt_type = 'auto';
%filt_steep = 0.95; 

% filter characteristics given by firFiltObj which is loaded here:
%load filter object from lowpass_characterization.m
disp('Loading filter...')
load fir_filt_10_100.mat %load FIR filter with fpass = 10 and fstop = 100 Hz
%get pass and stopband frequencies
filter_info = measure(firFiltObj)


%characteristic point search method
char_method = 'fft'; %fft for FFT, mean for mean value
%characteristic point search range
t_UIpair = 70;      %s before end of step where U/I pair should be taken
dt_UIpair = 60;      %s of data which should be averaged to take the U/I pair


%plotting y axis limits for current in voltage in nA / kV
ylims = [-15 15; 0 80]; %nA; kV
%limits for step plots
ylims_step = [-2 10; 25 65]; %nA; kV


%% read test params
disp('Reading test parameters...')
test_params = readmatrix(fn_test);

test_voltage = test_params(:,1);    %kV
test_duration = test_params(:,2);   %s
test_ramp = test_params(:,3);       %kV/s

%% read files

if ~exist(test_name, 'dir') %create test_name folder if it doesnt exist
   mkdir(test_name)
end

%pA sensor data
disp('Reading sensor data...')
sensordata=readmatrix(strcat(folder_sensor,fn_sensor));
v_sensor = sensordata(:,3); 
i_sensor = sensordata(:,2); 
t_sensor = sensordata(:,1);
%determine offset
start_time = 2  ;  %duration of start in seconds, when input is 0
n_samples = round(start_time*fs); %number of samples in given duration
disp('Sensor offsets:')
I_os = mean(i_sensor(1:n_samples))
V_os = mean(v_sensor(1:n_samples))
disp('Removing offsets...')
i_sensor = i_sensor-I_os;
v_sensor = v_sensor-V_os;
%MPD q data
disp('Reading MPD data...')
data_mpd = readmatrix(fn_mpd);
t_mpd = data_mpd(:,1);
q_mpd = data_mpd(:,2);
%remove any q smaller 1 pC to reduce data in RAM and plot
disp('Removing any q smaller 1 pC...')
k_rem = find(abs(q_mpd)<1e-12);
t_mpd(k_rem)=[];
q_mpd(k_rem)=[];
data_mpd = [t_mpd, q_mpd]; %update array 

%% plot full log data, raw and filtered

%plot unfiltered data
disp('Plotting raw data...')
figcnt = figcnt+1;
figure(figcnt)
plotPASenseData(sensordata,ylims,data_mpd);
titlestr = test_name_tex + " overview: raw data";
title(titlestr);
%save figure
saveas(figure(figcnt),strcat(test_name,'/',test_name,'-full'),'png');

%filter and plot same overview
disp('Filtering sensor data...')

%filter signals
v_sensor_filt = filtSignal(v_sensor,firFiltObj);
i_sensor_filt = filtSignal(i_sensor,firFiltObj);


%not using this one. as 'lowpass_characterization.m' shows, this filter is
%quite ill-designed for this purpose. will use specifically designed filter
%instead
% v_sensor_filt = lowpass(movmean(v_sensor,n_filt),f_lp,fs,ImpulseResponse=filt_type,Steepness=filt_steep);
% [i_sensor_filt, filt_obj]= lowpass(movmean(i_sensor,n_filt),f_lp,fs,ImpulseResponse=filt_type,Steepness=filt_steep);
% 
% save filter filt_obj
% save('filter.mat','filt_obj');

disp('Plotting filtered sensor data...')
figcnt=figcnt+1;
figure(figcnt)
plotPASenseData([t_sensor, i_sensor_filt, v_sensor_filt],ylims,data_mpd);
titlestr = test_name_tex + " overview: FIR filtered, F_{pass}=" + filter_info.Fpass + " Hz, F_{stop}="+filter_info.Fstop;
title(titlestr)
saveas(figure(figcnt),strcat(test_name,'/',test_name,'-full-filtered'),'png');


%% find step timings if unknown

if(step_duration_known==0)

    
    t_sensor = sensordata(:,1);
    v_sensor = lowpass(sensordata(:,3),0.1,fs,ImpulseResponse="iir",Steepness=0.8); %this filter specifically like this is needed for finding timings
    
    ptr = 1; %log file data pointer
    n_data = length(v_sensor); %number of datapoints
    for k = 1 : length(test_voltage)-1 %iterate through voltage steps
    
        %start time of this voltage step
        t1_step = t_sensor(ptr)
        %find step durations by time when logged voltage exceeds 3% of the
        %step delta to the next step

        test_delta_V(k+1) = (test_voltage(k+1)-test_voltage(k)); 
        if k == length(test_voltage)-1 %last step before ramping down to 0, part where it ramps down is not interesting
            while(v_sensor(ptr) > (test_voltage(k)+test_delta_V(k+1)*0.005)*1000 && ptr < n_data) %scaled from kV to V
                ptr = ptr+1;
            end
            test_duration(k) = t_sensor(ptr)-t1_step;
        else
            if k == 1 %first step, could get issues with filter edge effect, need to have a larger tolerance here
                tol = 0.15;
            else
                tol = 0.04;
            end
            while(v_sensor(ptr) < (test_voltage(k)+test_delta_V(k+1)*tol)*1000) %scaled from kV to V
                ptr = ptr+1;
            end
            t2_step = t_sensor(ptr)
            test_duration(k) = t2_step-t1_step;
            while(v_sensor(ptr) < (test_voltage(k)+test_delta_V(k+1)*0.97)*1000)
                ptr = ptr+1;
            end
            t3_step = t_sensor(ptr)
            test_ramp(k+1)=test_delta_V(k+1)/(t3_step-t2_step);
            %cont here
        end
    end

end

%% find index range of individual voltage steps based on test params

disp('Finding index range of voltage steps...')
ptr = 2; %log file position pointer - start from 2 because first value might have an error

for k = 1 : length(test_duration) %add -1 here if last step in test file is 0V

%calculate voltage difference from previous to this step
if k == 1 %first step is always 0 V
    delta_V = 0;
else 
    delta_V = test_voltage(k)-test_voltage(k-1);
end
dur = test_duration(k); %set duration of this step at voltage
ramp = test_ramp(k); %ramp in kV/s from previous to this step voltage
step_time = dur + abs(delta_V/ramp);  %full time of this step

n_samples = round(step_time*fs); %number of samples in given duration

%set range for this step
test_range(k,:) = [ptr, ptr+n_samples];

ptr = ptr + n_samples; %increment pointer

end

%% find characteristic values and plot individual voltage ranges
disp('Finding characteristic I/V and plotting steps...')
%doesnt plot when using parfor?
for k=1:length(test_duration)
disp('Step: ')
disp(k)

%check if index exceeds max data length, e.g. happens when test is
%prematurely ended
    if test_range(k,2) > length(t_sensor)
            test_range(k,2) = length(t_sensor);
    end

    t_range = t_sensor(test_range(k,1):test_range(k,2));
    i_range = i_sensor_filt(test_range(k,1):test_range(k,2));
    v_range = v_sensor_filt(test_range(k,1):test_range(k,2));

    %find I/V pair for VJ characteristic based on samples
    k_start = length(t_range)-t_UIpair*fs;
    if(k_start < 0)
        k_start = 1;
    end
    k_stop = k_start + dt_UIpair*fs;
    if(k_stop > length(t_range))
        k_stop = length(t_range)-2*fs;
    end

    if strcmp(char_method,'mean')
        I_char(k) = mean(i_range(k_start:k_stop));
        V_char(k) = mean(v_range(k_start:k_stop));
    elseif strcmp(char_method,'fft')
         %calculate FFT in given range
        [f_fft_i i_fft] = getFFT(fs,i_range(k_start:k_stop));
        [f_fft_v v_fft] = getFFT(fs,v_range(k_start:k_stop));
        %take DC component as characteristic value
        I_char(k) = i_fft(1);
        V_char(k) = v_fft(1);
    end

    %first step V = 0, overwrite
    if k==1
        V_char(k)=0;
    end


    %plot this range
    titlestr = test_name_tex + ": step "+test_voltage(k)+" kV";
    fig_k = figcnt+k;
    figure(fig_k)    
    if strcmp(fft, 'fft')
        fn = strcat(test_name,'/',test_name,'-',num2str(test_voltage(k)),'kV+fft.png');
        plotPASenseData([t_range, i_range, v_range],ylims_step,data_mpd,'fft',fft_lim);
    else
        fn = strcat(test_name,'/',test_name,'-',num2str(test_voltage(k)),'kV.png');
        plotPASenseData([t_range, i_range, v_range],ylims_step,data_mpd);
    end
    title(titlestr)
    saveas(figure(fig_k),fn,'png');
end

figcnt = figcnt+length(test_duration);

%% plot V/I char
disp('Plotting V/I characteristic...')
figcnt=figcnt+1;
figure(figcnt)
clf
hold on
grid on
grid minor
plot(V_char/1e3, I_char/1e-9,'*-')
titlestr = test_name_tex + ": V-I characteristic";
title(titlestr)
xlabel('V / kV')
ylabel('I / nA')

hold off

%save figure
fn = strcat(test_name,'/',test_name,'-VI-char-',char_method,'.png');
saveas(figure(figcnt),fn,'png');

%save characteristic in .mat file for use in comparative plots
disp('Saving VI characteristic data...')

char_data_folder = 'VI-characteristic data';
if ~exist(char_data_folder, 'dir') %create folder if doesn't exist
   mkdir(char_data_folder)
end
data_name = char_data_folder + "\" + test_name + "-VI-"+char_method;
data_char = [V_char', I_char'];
save(data_name, "data_char")

%% done / cleaning up
disp('Done!')

%% functions

%filter signal with supplied filter object, compensates for filter delays
%and edge effects
%acc. to matlab docu: https://de.mathworks.com/help/signal/ref/lowpass.html
%"fir" — The function designs a minimum-order, linear-phase, finite impulse response (FIR) filter. 
% To compensate for the delay, the function appends to the input signal N/2 zeros, where N is the 
% filter order. The function then filters the signal and removes the first N/2 samples of the output.

function [filt_data] = filtSignal(data,filt_obj)

    N = length(filt_obj.Numerator); %find filter order

    data(end+1:end+round(N/2)) = 0; %append zeros
    filt_data = filt_obj(data); %filter
    filt_data = filt_data(round(N/2)+1:end); %remove first N/2 outputs    

end

%not used in this version of the script
%plots filtered and unfiltered data next to each other
function [] = plotData(t,i,v,i_filt, v_filt, t_mpd, pd_mpd, I_char, V_char, titlestr,fft)
    global figcnt;
 

    %scale to nA and kV, and nC
    i = i/1e-9;
    v = v/1e3;
    i_filt = i_filt/1e-9;
    v_filt = v_filt/1e3;
    pd_mpd = pd_mpd/1e-9;

    %plotting y limits
    imin_lim = (min(i)-abs(min(i)*0.05));
    imax_lim = max(i)*1.05;
    vmin_lim = (min(v)-abs(min(v)*0.05));
    vmax_lim = max(v)*1.05;
    if(I_char == 0 || V_char == 0)
        ifiltmin_lim = (min(i_filt)-abs(min(i_filt)*0.01));
        ifiltmax_lim = max(i_filt)*1.01;
        vfiltmin_lim = (min(v_filt)-abs(min(v_filt)*0.01));
        vfiltmax_lim = max(v_filt)*1.01;
    else
        ifiltmin_lim = I_char * 0.7 /1e-9;
        ifiltmax_lim = I_char * 1.5 /1e-9;
        vfiltmin_lim = V_char * 0.99 /1e3;
        vfiltmax_lim = V_char * 1.01 /1e3;
    end

    figcnt = figcnt +1;
    figure(figcnt)
    clf
    hold on   
    ax1 = subplot(3,2,1);   
    plot(t,i)
    title(titlestr,'Interpreter', 'none')
    ax1.TitleFontSizeMultiplier = 0.8;
    grid on
    grid minor
    ylim([imin_lim imax_lim])
    ylabel('i / nA')
    xlabel('t / s')
    ax2 = subplot(3,2,3);
    plot(t,v)
    grid on
    grid minor
    ylim([vmin_lim vmax_lim])
    ylabel('v / kV')
    xlabel('t / s') 

    ax3 = subplot(3,2,2);
    plot(t,i_filt)
    title(strcat(titlestr,' filtered'),'Interpreter', 'none')
    ax3.TitleFontSizeMultiplier = 0.8;
    grid on
    grid minor
    ylim([ifiltmin_lim ifiltmax_lim])
    ylabel('i / nA')
    xlabel('t / s')
    ax4 = subplot(3,2,4);
    plot(t,v_filt)
    grid on
    grid minor
    ylim([vfiltmin_lim vfiltmax_lim])
    ylabel('v / kV')
    xlabel('t / s') 
    ax5 = subplot(3,2,6);
    plot(t_mpd,pd_mpd,'.')
    grid on
    grid minor
    %ylim([(min(v)-abs(min(v)*0.05)) max(v)*1.05])
    ylabel('Q / nC')
    xlabel('t / s')    


    linkaxes([ax1, ax2, ax3, ax4, ax5],'x');  
    hold off

    if strcmp(fft, 'fft') %plot fft
        figcnt = figcnt +1;
        figure(figcnt)
        fs = 1/ (t(2)-t(1));
        plotfft(fs,i,'lin');
        grid on
        grid minor
        ylabel('i / dBnA')

    end

end