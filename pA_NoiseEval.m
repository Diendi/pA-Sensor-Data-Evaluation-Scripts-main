%plots and compares pA-sensor noise measurement data
%
%Matthias Musil
%FH Wels Ring
%03.02.2024

%in this folder compares previous noise measurements from 1.2.2024 to newer
%ones, including battery, from 4.3.2024

clear
close all

global figcnt;
figcnt = 0;

colors = distinguishable_colors(5);
%% setup

folder = 'noise'; %folder that contains noise meausrements
fn_base = 'noise'; %base name of measurement files
tests_to_plot = [5, 2, 4, 6, 7]; %number of measurements to plot
tests_descriptions =  {'0 V battery 1'; '0 V PS 1';'0 V PS 2';'0 V PS 3';'0 V battery 2'}; %description of tests
test_length = 14; %seconds
test_fs = 12500; %sps
N = test_length*test_fs;
%% load and plot time data

figcnt = figcnt + 1;
figure(figcnt);
hold on
grid on
grid minor

for k = 1:length(tests_to_plot)
    fn = folder + "/" + fn_base + "_" + tests_to_plot(k)+ ".csv";
    sensordata = readmatrix(fn);
    if length(sensordata(:,1))>N
        sensordata = sensordata(1:N,:);
    end
    i(:,k) = sensordata(:,2)-mean(sensordata(:,2)); %removing offset also
    t(:,k) = sensordata(:,1);
    fs(k) = 1/(t(2,k)-t(1,k));

    %find and remove offset
    I_os = mean(i(:,k));
    i(:,k) = i(:,k)-I_os;

    %calculate noise RMS
    noise_rms(k)=sqrt(mean(abs(i(:,k)/1e-9).^2));
    disp("Noise rms of "+tests_to_plot(k)+": "+noise_rms(k)+" nA")

    plot(t(:,k),i(:,k)/1e-9,'Color',colors(k,:),'DisplayName',""+tests_descriptions{k}+": "+noise_rms(k)+" nA RMS")
    
end
xlabel('t / s')
ylabel('i / nA')
title('Comparison of noise')
legend show

%save plot
plot_folder = folder + "/plots2";
if ~exist(plot_folder, 'dir') %create test_name folder if it doesnt exist
   mkdir(plot_folder)
end
saveas(figure(figcnt),strcat(plot_folder,'/noise-comparison'),'png');

%% plot FFTs

figcnt = figcnt + 1;
figure(figcnt);
hold on

for k = 1:length(tests_to_plot)    

    plotfft(fs(k),i(:,k)/1e-9,'lin','Color',colors(k,:),'DisplayName',tests_descriptions{k});
    
end
grid on
grid minor
xlabel('f / Hz')
ylabel('i / dBnA')
title('Comparison of noise - FFT')
legend show
set(gcf,'Position',[100 100 1120 630])
%save full view
saveas(figure(figcnt),strcat(plot_folder,'/noise-comparison-fft'),'png');
%zoom to different view and save also
xlim([0 1100])
saveas(figure(figcnt),strcat(plot_folder,'/noise-comparison-fft-zoom'),'png');