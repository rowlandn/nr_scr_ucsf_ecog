
function [Comodulogram,pComodulogram,zComodulogram] = pac_art_reject_surr(data,srate,PhaseFreqVector,AmpFreqVector,bad_times,skip,ch_plot) 
%%uses Tort method to calculate PAC. Will do statistics if you specify a
%%variable "skip", will exclude noisy segments of data if you specify "band
%%times"
%inputs:
%data= data (from 1 channel)you want to analyze as a row vector
%srate = sampling rate
%PhaseFreqVector = vector of frequences to analyze for phase in Hz ex. PhaseFreqVector=[4:2:45];
%AmpFreqVector = vector of frequences to analyze for amplitude in Hz ex. AmpFreqVector=[50:4:150];
%bad_times = vector of data to exclude (b/c of artifacts etc), in units of sample points, leave empty if none
%Note that these are excluded after filtering to avoid edge artifacts.
%skip = vector of random time offsets to run for surrogated. Length of skip
%is number of surrogates. Recommend getting skip using "make skip list". leave skip empty if you don't want to run statistics

%outputs
%Comodulogram = Comodulogram of modulation indexes for each phase and
%ampltide frequency examines (matrix phase x amp)

%pComodulogram = matrix same dimensions of Comodulogram which has the
%p-value for each data point. Note that p value is calculated as percent
%surrogates with larger ampltidue than the real data. It is possible to get
%p values of 0. Still a work in progress.
%
%zComodulogram = matrix same dimensions of Comodulogram which has the
%zscore for each data point. Note that z value is calculated as a real
%value - mean of surrogates/ sd of surrogates. This may not be proper if
%the surrogate distrubition is not normal. Still a work in progress. 
%


%initialize variables. 
if ~isempty(skip)
numsurrogate = length(skip);
end

%data_length = length(signal);

data_length = length(data);

Comodulogram=single(zeros(length(PhaseFreqVector),length(AmpFreqVector)));
pComodulogram=single(zeros(length(PhaseFreqVector),length(AmpFreqVector)));
zComodulogram=single(zeros(length(PhaseFreqVector),length(AmpFreqVector)));

AmpFreqTransformed = zeros(length(AmpFreqVector), data_length-length(bad_times));
PhaseFreqTransformed = zeros(length(PhaseFreqVector), data_length-length(bad_times));
    
    
    %% For comodulation calculation (only has to be calculated once)
    nbin = 18;
    position=zeros(1,nbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
    winsize = 2*pi/nbin;
    for j=1:nbin
        position(j) = -pi+(j-1)*winsize;
    end
    
    'CPU filtering'
    
    
    
    PhaseFreq_BandWidth=2;
    AmpFreq_BandWidth=4;

    %loop through all frequencies for amplitude, filter at
    %desired frequencies, and store the amplitude of the signal
    for ii=1:length(AmpFreqVector)
        Af1 = AmpFreqVector(ii);
        Af2=Af1+AmpFreq_BandWidth;
        AmpFreq=eegfilt_nicki(data,srate,Af1,Af2); % just filtering
        temp = hilbert(AmpFreq);
        temp(bad_times) = [];
        AmpFreqTransformed(ii, :) = abs(temp); % getting the amplitude envelope
        clear temp
    end
    
    %loop through all frequencies for phase, filter at
    %desired frequencies, and store the phase of the signal
    for jj=1:length(PhaseFreqVector)
        Pf1 = PhaseFreqVector(jj);
        Pf2 = Pf1 + PhaseFreq_BandWidth;
        %PhaseFreq=eegfilt_nicki(data,srate,Pf1,Pf2); % this is just filtering
        PhaseFreq=eegfilt_nicki(data,srate,Pf1,Pf2); % this is just filtering
        temp = hilbert(PhaseFreq);
        temp(bad_times) = [];
        PhaseFreqTransformed(jj, :) = angle(temp); % this is getting the phase time series
    end
    
  
    %% Do comodulation calculation
    'Comodulation loop'
    
    counter1=0;
    for ii=1:length(PhaseFreqVector)
        ii
        counter1=counter1+1;
        ii
        counter2=0;
        for jj=1:length(AmpFreqVector)
            counter2=counter2+1;
            jj
            %main modulation index calculation and save modulation index
            %value which reflects the degree of PAC for each phase and
            %ampltude frequency
            [MI,MeanAmp]=ModIndex_v2(PhaseFreqTransformed(ii, :), AmpFreqTransformed(jj, :), position);
            Comodulogram(counter1,counter2)=MI;
            
            %run surrogates (statistics) if desired
            if ~isempty(skip)
                surrogate_m=zeros(numsurrogate,1);
                
                
                for s=1:numsurrogate
                    Amp_surr =[AmpFreqTransformed(jj,skip(s):end) AmpFreqTransformed(jj,1:skip(s)-1)];
                    [MI_S,MeanAmp_S]=ModIndex_v2(PhaseFreqTransformed(ii,:), Amp_surr, position);
                    MI_surr(s) = MI_S;
                    
                end
                high_surr = length(find(abs(MI_surr)>=(abs(Comodulogram(counter1,1)))));
                %calculate statistics
                pComodulogram(counter1,1) = high_surr./numsurrogate;
                zComodulogram(counter1,1) = (Comodulogram(counter1,1) - mean(abs(MI_surr)))./std(abs(MI_surr));
            end
        end
    end
    
    
%     %%plot figure if desired
%     figure
%     h = subplot(1,1,1);
%             cmax=max(abs(Comodulogram(:)))*.75;
%             cmin=0;
%             tempmat=double(squeeze(Comodulogram(:,:,ch_plot)));
%             pcolor(PhaseFreqVector',AmpFreqVector,tempmat');
%             shading interp;
%             caxis([cmin cmax]);
%             colorbar;
%             hold on;
%              title(['ch',int2str(ch_plot)]);
%             hold on;
            
            


    
   