%plots together VI characteristics of differet pA sensor tests
%
%Matthias Musil
%FH Wels RING
%03.02.2024

%constants
global figcnt;
char_data_folder = 'VI-characteristic data'; %folder containing data
char_method = 'fft'; %how VI char was found, fft for FFT, mean for mean value

%enter test numbers to plot in array below
%script will fetch according data and plot them together
tests_to_plot = [59, 60, 61];
test_description = {'clean','lightly polluted','heavily polluted'};
figcnt = figcnt+1;
figure(figcnt)
clf
hold on
grid on
grid minor
for k = 1:length(tests_to_plot)
    
    data_name = char_data_folder + "\test_" + tests_to_plot(k) + "-VI-"+char_method+".mat";
    load(data_name);
    V_char = data_char(:,1);
    I_char = data_char(:,2);
    plot(V_char/1e3, I_char/1e-9,'*-','DisplayName',strcat('test\_',num2str(tests_to_plot(k)),': ',test_description{k}))    
end
titlestr = 'V-I characteristics of multiple tests';
title(titlestr)
xlabel('V / kV')
ylabel('I / nA')
legend('Location','northwest')
legend show
hold off    

char_plots_folder = 'VI-characteristics plots';
if ~exist(char_plots_folder, 'dir') %create folder if doesn't exist
   mkdir(char_plots_folder)
end
fn = char_plots_folder + "\" + "VI-char-tests";
saveas(figure(figcnt),fn,'png');
