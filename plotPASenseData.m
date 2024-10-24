%plots specified pA Sensor data
%plots current and voltage in separate graphs
%if MPD data is passed, it will also plot q in the voltage graph
%if desired, performs and plots the fft besides the other graphs
%enter limits in format [imin imax; vmin vmax]
function [fftdata] = plotPASenseData(data_sensor,ylims,data_mpd, fft,fftlim)

vsense = data_sensor(:,3)/1e3; 
isense = data_sensor(:,2)/1e-9;
tsense = data_sensor(:,1);

clf
hold on



switch nargin
    case 5 %include fft with frequency limits

        if strcmp(fft,'fft')
            
            if ~strcmp(data_mpd,'')
                t_mpd = data_mpd(:,1);
                q_mpd = data_mpd(:,2)/1e-12;
                %remove any q smaller 1 pC to reduce data in plot
                k_rem = find(abs(q_mpd)<1);
                t_mpd(k_rem)=[];
                q_mpd(k_rem)=[];
            end
    
            fs = 1/(tsense(2)-tsense(1));
            %limit data for FFT(opt)
            %n = round(20*fs);
            n = length(isense);

            ax1 = subplot(3,1,1);
            plot(tsense,isense)
            grid on
            grid minor
            ylabel('i / nA')
            xlabel('t / s')
            if ~strcmp(ylims,'') 
                ylim(ylims(1,:)) 
            end
            xlim([tsense(1) tsense(end)]) %limit time axis 
            ax2 = subplot(3,1,2);
            if ~strcmp(data_mpd,'')
                yyaxis left
            end
            plot(tsense,vsense)
            grid on
            grid minor
            ylabel('v / kV')
            if ~strcmp(ylims,'') 
                ylim(ylims(2,:))
            end
            if ~strcmp(data_mpd,'')
                yyaxis right
                plot(t_mpd,q_mpd,'.')
                ylabel('q / pC')
                xlabel('t / s')
                align_yyaxis_zero(ax2)
            end
            xlim([tsense(1) tsense(end)]) %limit time axis 
            ax3 = subplot(3,1,3);
            fftdata = plotfft(fs,isense(1:n),'log');
            grid on
            grid minor
            ylabel('i / dBnA')
            xlabel(' f / Hz')
            if ~strcmp(fftlim,'')
                xlim(fftlim)
            end
            linkaxes([ax1, ax2],'x');
            set(gcf,'Position',[100 100 560 630])
            subplot(3,1,1)                 
        end

    case 4 %include fft 

        if strcmp(fft,'fft')
            
            if ~strcmp(data_mpd,'')
                t_mpd = data_mpd(:,1);
                q_mpd = data_mpd(:,2)/1e-12;
                %remove any q smaller 1 pC to reduce data in plot
                k_rem = find(abs(q_mpd)<1);
                t_mpd(k_rem)=[];
                q_mpd(k_rem)=[];
            end
    
            fs = 1/(tsense(2)-tsense(1));
            %limit data for FFT(opt)
            %n = round(20*fs);
            n = length(isense);

            ax1 = subplot(3,1,1);
            plot(tsense,isense)
            grid on
            grid minor
            ylabel('i / nA')
            xlabel('t / s')
            if ~strcmp(ylims,'') 
                ylim(ylims(1,:)) 
            end
            xlim([tsense(1) tsense(end)]) %limit time axis 
            ax2 = subplot(3,1,2);
            if ~strcmp(data_mpd,'')
                yyaxis left
            end
            plot(tsense,vsense)
            grid on
            grid minor
            ylabel('v / kV')
            if ~strcmp(ylims,'') 
                ylim(ylims(2,:))
            end
            if ~strcmp(data_mpd,'')
                yyaxis right
                plot(t_mpd,q_mpd,'.')
                ylabel('q / pC')
                xlabel('t / s')
                align_yyaxis_zero(ax2)
            end
            xlim([tsense(1) tsense(end)]) %limit time axis 
            ax3 = subplot(3,1,3);
            fftdata = plotfft(fs,isense(1:n),'log');
            grid on
            grid minor
            ylabel('i / dBnA')
            xlabel(' f / Hz')
            linkaxes([ax1, ax2],'x');
            set(gcf,'Position',[100 100 560 630])
            subplot(3,1,1)                 
        end



    case 3 %include mpd data

        if ~strcmp(data_mpd,'')
         
            t_mpd = data_mpd(:,1);
            q_mpd = data_mpd(:,2)/1e-12; %convert to pC
    
            %remove any q smaller 1 pC to reduce data in plot
            k_rem = find(abs(q_mpd)<1);
            t_mpd(k_rem)=[];
            q_mpd(k_rem)=[];

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

            xlim([tsense(1) tsense(end)]) %limit time axis 

            fftdata=0;
            set(gcf,'Position',[100 100 560 630])
            subplot(2,1,1)
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
        xlim([tsense(1) tsense(end)]) %limit time axis
        fftdata=0;
        set(gcf,'Position',[100 100 560 630])
        subplot(2,1,1)
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
        xlim([tsense(1) tsense(end)]) %limit time axis 
        fftdata=0;
        set(gcf,'Position',[100 100 560 630])
        subplot(2,1,1)
end





end

