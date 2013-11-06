function [Comodulogram] = make_comodulogram_plv(data,Fs,PhaseLow,PhaseHigh,AmpLow,AmpHigh,bad_times,Ampwidth)
%%used PLV method to calculate PAC. No statistics (but you can modify to do, use Tort script as a reference),
%will exclude noisy segments of data if you specify "badtimes"
%inputs:
%data= data (from 1 channel)you want to analyze as a row vector
%Fs = sampling rate
%PhaseLow = Low frequency value for phase in Hz (ex. 4). 
%PhaseHigh = high value of phase frequency vector in Hz (ex. 45).
%bad_times = vector of data to exclude (b/c of artifacts etc), in units of sample points, leave empty if none. Note removing the bad 
%segments of data is done after filtering to avoid edge artifacts.
%Ampwidth = the width of each band in the amplitude component of the PAC (the phase component is 1 by default)
%5 for amp and 1 for phase is what Canolty 2006 used (ex. 5)

%outputs
%Comodulogram = Comodulogram of PLV values for each phase and
%ampltide frequency examines (matrix phase x amp)


%initialize variables
PhaseFreqVector=[PhaseLow:PhaseHigh];
AmpFreqVector=[AmpLow:Ampwidth:AmpHigh];
Comodulogram = NaN(length(PhaseFreqVector),length(AmpFreqVector));
    %loop through phase frequencies
    phaseCount=1;
    for p=PhaseLow:PhaseHigh
        ampCount=1;
        %filter data at phase frequencies with a bandwidth of 1. 
        filtered1=eegfilt(data,Fs,p,[]);
        filtered1=eegfilt(filtered1,Fs,[],p+1);
        %take angle of filtered data
        slow =  angle(hilbert(filtered1));
        %loop through amplitude frequencies
        for a=AmpLow:Ampwidth:AmpHigh
            %take amplitude of filtered data for amplitude frequencies
            fast =   abs(hilbert(eegfilt(data,Fs,a,a+Ampwidth-1)));
            %filter amplitude envelope at low frequency with bandwidth of 1
            fast = eegfilt(fast,Fs,p,[]);
            fast = eegfilt(fast,Fs,[],p+1);
            %take angle of amplitude envelope filtered signal.
            fast = angle(hilbert(fast));
           
         
            x = slow;
            y = fast;
            
            %remove time points with noise
            if (~isempty(bad_times))
                x(bad_times)=[];
                y(bad_times)=[];
            end
            
            %calculate the plv value as the difference in phase angle between low freq filtered signal and amplitude envelope
            %filtered signal. Note that the difference is multipled by i and taken to the exponent b/c values are circular. 
            %then sum plv values over time points 
            plvs=exp(i*(y(1)-x(1)));
            for ee=2:length(x)
                plvs=plvs+exp(i*(y(ee)-x(ee)));
            end
            %take magnitide of plv (vector length), and divide by number of
            %time points (i.e. average plv)
            PLV=abs(plvs)/length(x);
            %put in correct place in matrix
            Comodulogram(phaseCount,ampCount)=PLV;
            clear PLV
            
            ampCount=ampCount+1;
        end
%        fprintf('\n*****\n\n\n\nfinished about %2.0f percent of electrode %d\n\n\n\n*****\n', phaseCount/gencnt*100, e)
        phaseCount=phaseCount+1;
    end
   
end