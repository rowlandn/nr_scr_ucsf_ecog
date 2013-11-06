
function Comodulogram_rest_all_chan(ecog,name,M1_ch,Fs, surrogate)

load(name)
name = strrep(name,'_Com_chan.mat','');
name = strrep(name,'_ecog.mat','');

%% Define the Amplitude- and Phase- Frequencies
% PhaseFreqVector=[4:2:40];
% AmpFreqVector=[50:4:200];
% Comodulogram = nan*zeros(24,98,6);
PhaseFreqVector=[4:2:50];
AmpFreqVector=[10:4:400];

PhaseFreq_BandWidth=2;
AmpFreq_BandWidth=4;

srate=Fs;
% if  length(ecog.contact_pair)>6
%     M1_ch1=M1_ch1;
%     M1_ch2=M1_ch2;
% end

%% For comodulation calculation (only has to be calculated once)
nbin = 18;
position=zeros(1,nbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
winsize = 2*pi/nbin;
for j=1:nbin
    position(j) = -pi+(j-1)*winsize;
end

%% Do the remontage
for i = 1: length(ecog.contact_pair)
    if i<5
        ecog.contact_pair(i).remontaged_ecog_signal = ecog.contact_pair(i).raw_ecog_signal - ecog.contact_pair(i+1).raw_ecog_signal;
    else
        ecog.contact_pair(i).remontaged_ecog_signal = ecog.contact_pair(i).raw_ecog_signal;
    end
end
% ecog.contact_pair(i).remontaged_ecog_signal = ecog.contact_pair(i).remontaged_ecog_signal
%%Compute COM for each ecog pair or contact
for chan =1:length(ecog.contact_pair)
    %% Remove zeros added at the end of the ecog signal
    chan
         lfp = ecog.contact_pair(chan).remontaged_ecog_signal;
        
         WINDOW = 1024;
         NOVERLAP = 512;
         NFFT = 1024;
         % remove artifacts
         [n1_b, n1_a]=butter(3,2*[118 122]/Fs,'stop'); %120hz
         lfp=filtfilt(n1_b, n1_a, lfp);
         [n1_b, n1_a]=butter(3,2*[178 182]/Fs,'stop'); %180hz
         lfp=filtfilt(n1_b, n1_a, lfp);
         if ~isempty(strfind(name,'DBSON'))
             %remove stim artifact
             [psd,f] = pwelch(lfp,WINDOW ,NOVERLAP,NFFT,Fs);
             freq = find(f>= 125 & f<= 175);
             [v,p] = max(psd(freq));
             f_max = p + freq(1)-1;
             f_stim = f(f_max);
             [n1_b, n1_a]=butter(3,2*[f_stim-3 f_stim+3]/Fs,'stop'); %180hz
             lfp=filtfilt(n1_b, n1_a, lfp);
         end

%     lfp = ecog.contact_pair(chan).raw_ecog_signal;
%     lfp = lfp(find(lfp~=0));
    lfp = lfp-mean(lfp);
            
%     if length(lfp)>=60000 % file must have the same duration 1 min
%         lfp=lfp(:,1:60000);
%     elseif length(lfp)>=30000
if length(lfp)>=30000
        lfp=lfp(:,1:30000);
    end
%     if length(lfp)>=30000
%         lfp=lfp(:,1:30000);
%     end
    data_length = length(lfp);
    
    if ~isempty(lfp)
        
        %% Do filtering and Hilbert transform on CPU
        
        'CPU filtering'
        tic
        AmpFreqTransformed = zeros(length(AmpFreqVector), data_length);
        PhaseFreqTransformed = zeros(length(PhaseFreqVector), data_length);
        
        for ii=1:length(AmpFreqVector)
            Af1 = AmpFreqVector(ii)-AmpFreq_BandWidth/2;
            Af2=AmpFreqVector(ii)+AmpFreq_BandWidth/2;
            AmpFreq=eegfilt(lfp,srate,Af1,Af2); % just filtering
            AmpFreqTransformed(ii, :) = abs(hilbert(AmpFreq)); % getting the amplitude envelope
        end
        
        for jj=1:length(PhaseFreqVector)
            Pf1 = PhaseFreqVector(jj) - PhaseFreq_BandWidth/2;
            Pf2 = PhaseFreqVector(jj) + PhaseFreq_BandWidth/2;
            PhaseFreq=eegfilt(lfp,srate,Pf1,Pf2); % this is just filtering
            PhaseFreqTransformed(jj, :) = angle(hilbert(PhaseFreq)); % this is getting the phase time series
        end
        % PhaseFreqTransformed = PhaseFreqTransformed(:,500: end-500);
        % AmpFreqTransformed = AmpFreqTransformed(:,500: end-500);
        toc
        
        % clear 'PhaseFreqTransformed' 'AmpFreqTransformed'
        %% Do comodulation calculation
        'Comodulation loop'
        
        counter1=0;
        for ii=1:length(PhaseFreqVector)
            counter1=counter1+1;
            ii
            counter2=0;
            for jj=1:length(AmpFreqVector)
                counter2=counter2+1;
                
                [MI,MeanAmp]=ModIndex_v2(PhaseFreqTransformed(ii, :), AmpFreqTransformed(jj, :), position);
                Comodulogram(counter1,counter2,chan)=MI;
                x = 10:20:720;
                [val,pos]=max(MeanAmp);
                Comodulogram_phase(counter1,counter2,chan) = x(pos);
                if surrogate ==1
                    if chan == M1_ch-1 || chan == M1_ch
                        numpoints=size(AmpFreqTransformed,2);
                        numsurrogate=200; %% number of surrogate values to compare to actual value
                        minskip=Fs; %% time lag must be at least this big
                        maxskip=numpoints-Fs; %% time lag must be smaller than this
                        skip=ceil(numpoints.*rand(numsurrogate*2,1));
                        skip((skip>maxskip))=[];
                        skip(skip<minskip)=[];
                        skip=skip(1:numsurrogate,1); % creates vector with of time lags "tau" (the skip values) used for surrogate MIs
                        surrogate_m=zeros(numsurrogate,1);
                        
                        for s=1:numsurrogate
                            Amp_surr =[AmpFreqTransformed(jj,skip(s):end) AmpFreqTransformed(jj,1:skip(s)-1)];
                            [MI_S,MeanAmp_S]=ModIndex_v2(PhaseFreqTransformed(ii,:), Amp_surr, position);
                            MI_surr(s) = MI_S;
                        end
                        
                        % fit gaussian to surrogate data, uses normfit.m from MATLAB Statistics toolbox
                        [surrogate_mean,surrogate_std]=normfit(MI_surr);
                        Mean_surr(counter1,counter2,chan)=surrogate_mean;
                        Std_surr(counter1,counter2,chan)=surrogate_std;
                        Comodulogram_surr(counter1,counter2,chan)=(abs(MI)-surrogate_mean)/surrogate_std;
                        p_surr(counter1,counter2,chan)=  prctile(MI_surr,99);
                        Comodulogram_surr(counter1,counter2,chan)=MI;
                        if MI<prctile(MI_surr,99)
                            Comodulogram_surr(counter1,counter2,chan)=0;
                        end
                        
                    else
                        Comodulogram_surr(counter1,counter2,chan)=MI;
                        Mean_surr(counter1,counter2,chan)=nan;
                        Std_surr(counter1,counter2,chan)=nan;
                        p_surr(counter1,counter2,chan)=nan;
                        %                     Comodulogram_surr(counter1,counter2,chan)=(abs(Comodulogram(counter1,counter2,chan))-Mean_surr(counter1,counter2,chan))./Std_surr(counter1,counter2,chan);
                    end
                else
                    Comodulogram_surr(counter1,counter2,chan)=MI;
                    Mean_surr(counter1,counter2,chan)=nan;
                    Std_surr(counter1,counter2,chan)=nan;
                    p_surr(counter1,counter2,chan)=nan;
                    %                     Comodulogram_surr(counter1,counter2,chan)=(abs(Comodulogram(counter1,counter2,chan))-Mean_surr(counter1,counter2,chan))./Std_surr(counter1,counter2,chan);
                end
            end
        end
        toc
    end
end
% if length(ecog.contact_pair)<=6
save([name '_Com_chan'],'ecog', 'Comodulogram_surr','Comodulogram','Comodulogram_phase','p_surr','Mean_surr','Std_surr','M1_ch','PhaseFreqVector','PhaseFreq_BandWidth','AmpFreqVector','AmpFreq_BandWidth');
% else
%     save([name '_Com_chan'],'ecog', 'Comodulogram_surr','Comodulogram','Comodulogram_phase','p_surr','Mean_surr','Std_surr','M1_ch1','M1_ch2','PhaseFreqVector','PhaseFreq_BandWidth','AmpFreqVector','AmpFreq_BandWidth');
% end

%% Graph comodulogram
Comodulogram_surr(1,1,:)=0.00001;
Clim2 = max(max(max(Comodulogram(:,:,M1_ch))));
Clim1 = min(min(min(Comodulogram(:,:,M1_ch))));

necog=length(ecog.contact_pair);
if necog <=6
    length_row = necog;
elseif necog <=32
    length_row = 14;
elseif necog >= 64;
    length_row = 16;
end

row=1;
for ii = 1:length_row:necog
    hf1 = figure;
    hold on
    for i = 1:length_row
        chan=ii+i-1;
        
        if length_row <=6
            subplot(3,3,i)
        else
            subplot(4,4,i)
        end
        
        C=squeeze(Comodulogram(:,:,chan));
        contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,C',30,'lines','none')
        set(gca,'fontsize',14)
        ylabel('Amplitude Frequency (Hz)')
        xlabel('Phase Frequency (Hz)')
        colorbar
        caxis([Clim1 Clim2])
        if chan==1
            title([name num2str(chan) ]);
        elseif chan==M1_ch
            title(['M1=' num2str(chan) ],'FontWeight','b');
        else
            title( num2str(chan) );
        end
    end
    row = row;
    saveas(gcf,[name '_Com_raw' num2str(ii)],'fig');
    row = row+1;
end
close all

%
% for chan = 1:length(ecog.contact_pair)
%     if chan==1
%         figure;
%         hold on
%     end
%     if length(ecog.contact_pair)<=6
%         subplot(3,3,chan)
%         C=squeeze(Comodulogram(:,:,chan));
%         contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,C',30,'lines','none')
%         set(gca,'fontsize',14)
%         ylabel('Amplitude Frequency (Hz)')
%         xlabel('Phase Frequency (Hz)')
%         colorbar
%         caxis([Clim1 Clim2])
%         if chan==1
%             title([name num2str(chan) ]);
%         elseif chan==M1_ch
%             title(['M1=' num2str(chan) ],'FontWeight','b');
%         else
%             title( num2str(chan) );
%         end
%         if chan>=5
%             saveas(gcf,[name '_Com_chan'],'fig');
%         end
%     else
%         if chan<15
%             subplot(3,5,chan)
%             C=squeeze(Comodulogram_surr(:,:,chan));
%             contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,C',30,'lines','none')
%             set(gca,'fontsize',14)
%             ylabel('Amplitude Frequency (Hz)')
%             xlabel('Phase Frequency (Hz)')
%             colorbar
%             caxis([Clim1 Clim2])
%             if chan==1
%                 title([name num2str(chan) ]);
%             elseif chan==M1_ch1
%                 title(['M1=' num2str(chan) ],'FontWeight','b');
%             else
%                 title( num2str(chan) );
%             end
%             if chan==14
%                 saveas(gcf,[name '_Com_raw1'],'fig');
%             end
%         else
%             if chan==15
%                 figure;
%                 hold on
%             end
%             subplot(3,5,chan-14)
%             C=squeeze(Comodulogram_surr(:,:,chan));
%             contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,C',30,'lines','none')
%             set(gca,'fontsize',14)
%             ylabel('Amplitude Frequency (Hz)')
%             xlabel('Phase Frequency (Hz)')
%             colorbar
%             caxis([Clim1 Clim2])
%             if chan==15
%                 title([name num2str(chan) ]);
%             elseif chan==M1_ch2
%                 title(['M1=' num2str(chan) ],'FontWeight','b');
%             else
%                 title( num2str(chan) );
%             end
%
%             if chan==length(ecog.contact_pair)
%                 saveas(gcf,[name '_Com_raw2'],'fig');
%             end
%         end
%     end
% end
%
%
% close all
