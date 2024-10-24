%plots specified pA Sensor data
%plots current and voltage in separate graphs
%if MPD data is passed, it will also plot q in the voltage graph
%if desired, performs and plots the fft besides the other graphs
%enter limits in format [imin imax; vmin vmax]

%TODO: have this function use plotPASenseData instead of direct code
function [sensordata, data_mpd, fftdata] = plotPASenseDataFromFile(fn,ylims,fn_mpd, fft)


sensordata=readmatrix(fn);
vsense = sensordata(:,3)/1e3; 
isense = sensordata(:,2)/1e-9;
tsense = sensordata(:,1);

clf
hold on



switch nargin
    case 4 %include fft 

        if strcmp(fft,'fft')
            
            if ~strcmp(fn_mpd,'')
                data_mpd = readmatrix(fn_mpd);
                t_mpd = data_mpd(:,1);
                q_mpd = data_mpd(:,2)/1e-12;
                %remove any q smaller 1 pC to reduce data in plot
                k_rem = find(abs(q_mpd)<1);
                t_mpd(k_rem)=[];
                q_mpd(k_rem)=[];
                data_mpd = [t_mpd, q_mpd*1e-12]; %updata array for return
                %converted back to actual value of C to give out raw data
            end
    
            fs = 1/(tsense(2)-tsense(1));
            %limit data for FFT(opt)
            %n = round(20*fs);
            n = length(isense);

            ax1 = subplot(2,2,1);
            plot(tsense,isense)
            grid on
            grid minor
            ylabel('i / nA')
            xlabel('t / s')
            if ~strcmp(ylims,'') 
                ylim(ylims(1,:)) 
            end
            ax2 = subplot(2,2,3);
            yyaxis left
            plot(tsense,vsense)
            grid on
            grid minor
            ylabel('v / kV')
            if ~strcmp(ylims,'') 
                ylim(ylims(2,:))
            end
            if ~strcmp(fn_mpd,'')
                yyaxis right
                plot(t_mpd,q_mpd,'.')
                ylabel('q / pC')
                xlabel('t / s')
                align_yyaxis_zero(ax2)
            end
            ax3 = subplot(2,2,[2 4]);
            fftdata = plotfft(fs,isense(1:n),'lin');
            grid on
            grid minor
            ylabel('i / dBnA')
            
            linkaxes([ax1, ax2],'x');
            xlim([0 tsense(end)]) %limit time axis                        
        end



    case 3 %include mpd data

        if ~strcmp(fn_mpd,'')
         
            data_mpd = readmatrix(fn_mpd);
            t_mpd = data_mpd(:,1);
            q_mpd = data_mpd(:,2)/1e-12; %convert to pC
    
            %remove any q smaller 1 pC to reduce data in plot
            k_rem = find(abs(q_mpd)<1);
            t_mpd(k_rem)=[];
            q_mpd(k_rem)=[];
            data_mpd = [t_mpd, q_mpd]; %updata array for return

            ax1 = subplot(2,1,1);
            plot(tsense,isense)
            grid on
            grid minor
            ylabel('i / nA')
            xlabel('t / s')
            if ~strcmp(ylims,'') 
                ylim(ylims(1,:)) 
            end
            ax2 = subplot(2,1,2);
            yyaxis left
            plot(tsense,vsense)
            grid on
            grid minor
            ylabel('v / kV')
            if ~strcmp(ylims,'') 
                ylim(ylims(2,:))
            end
            yyaxis right
            plot(t_mpd,q_mpd,'.')
            ylabel('q / pC')
            xlabel('t / s')            
            align_yyaxis_zero(ax2)

            linkaxes([ax1, ax2],'x');

            xlim([0 tsense(end)]) %limit time axis

            fftdata=0;
        end
    case 2 %just sensor data with limits
        ax1 = subplot(2,1,1);
        plot(tsense,isense)
        grid on
        grid minor
        ylabel('i / nA')
        xlabel('t / s')
        ylim(ylims(1,:))
        ax2 = subplot(2,1,2);
        plot(tsense,vsense)
        grid on
        grid minor
        ylabel('v / kV')
        xlabel('t / s') 
        ylim(ylims(2,:))
        linkaxes([ax1, ax2],'x');
        xlim([0 tsense(end)]) %limit time axis
        fftdata=0;
    case 1 %just sensor data no limits
        ax1 = subplot(2,1,1);
        plot(tsense,isense)
        grid on
        grid minor
        ylabel('i / nA')
        xlabel('t / s')
        ax2 = subplot(2,1,2);
        plot(tsense,vsense)
        grid on
        grid minor
        ylabel('v / kV')
        xlabel('t / s')        
        linkaxes([ax1, ax2],'x');
        xlim([0 tsense(end)]) %limit time axis
        fftdata=0;
end





end

