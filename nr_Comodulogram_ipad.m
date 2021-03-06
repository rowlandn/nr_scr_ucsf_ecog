
function Comodulogram_ipad(ecog,aux,emg,name,M1,Fs,surr,trials_ok) %(name,surr,M1) %

% load(name)
% 
% name = strrep(name,'_ecog.mat','');
name = strrep(name,'_ecog.mat','');

srate=Fs;
epoch = 2*Fs; % 2 sec duration
WINDOW = 1024;
NOVERLAP = 512;
NFFT = 1024;
lfp = ecog.contact_pair(M1).remontaged_ecog_signal;
lfp = lfp-mean(lfp);

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
data_length = length(lfp);


ecog.rest_time = int32(ecog.rest_time);
ecog.prep_time = int32(ecog.prep_time);
ecog.go_time = int32(ecog.go_time);
ecog.touch_time = int32(ecog.touch_time);
ecog.move_time = int32(ecog.active_time);
ecog.end_time = int32(ecog.trial_end_time);

ecog.rest_time=[ecog.trial_beg_time(1) ecog.rest_time];
%% Define the Amplitude- and Phase- Frequencies

PhaseFreqVector=[4:2:50];
AmpFreqVector=[50:4:200];

PhaseFreq_BandWidth=2;
AmpFreq_BandWidth=4;

%% For comodulation calculation (only has to be calculated once)
nbin = 18;
position=zeros(1,nbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
winsize = 2*pi/nbin;
for j=1:nbin
    position(j) = -pi+(j-1)*winsize;
end
%% Do filtering and Hilbert transform on CPU

'CPU filtering'
tic
Comodulogram=single(zeros(length(PhaseFreqVector),length(AmpFreqVector)));
AmpFreqTransformed = zeros(length(AmpFreqVector), data_length);
PhaseFreqTransformed = zeros(length(PhaseFreqVector), data_length);

for ii=1:length(AmpFreqVector)
    Af1 = AmpFreqVector(ii)-AmpFreq_BandWidth/2;
    Af2=AmpFreqVector(ii)+AmpFreq_BandWidth/2;
    AmpFreq=eegfilt(lfp,srate,Af1,Af2); % just filtering
    AmpFreqTransformed(ii, :) = abs(hilbert(AmpFreq)); % getting the amplitude envelope
end

for jj=1:length(PhaseFreqVector)
    Pf1 = PhaseFreqVector(jj)- PhaseFreq_BandWidth/2;
    Pf2 = PhaseFreqVector(jj) + PhaseFreq_BandWidth/2;
    PhaseFreq=eegfilt(lfp,srate,Pf1,Pf2); % this is just filtering
    PhaseFreqTransformed(jj, :) = angle(hilbert(PhaseFreq)); % this is getting the phase time series
end
% PhaseFreqTransformed = PhaseFreqTransformed(:,500: end-500);
% AmpFreqTransformed = AmpFreqTransformed(:,500: end-500);
toc

%% Do comodulation calculation
'Comodulation loop'
psd_all=[]
for tt = 1:4 % 1 == all trials, 2 == rest, 3== preparation, 4== mvt
    period_t =[];

    if tt ~= 1 % analyse COM for all trials
        if tt == 2
            for nn =trials_ok
                period = 'rest';
%                 period_t = [period_t ecog.prep_time(nn)- epoch:ecog.prep_time(nn)];
                period_t = [period_t ecog.rest_time(nn)+500:ecog.prep_time(nn)-500];

            end
        elseif tt== 3 %prep period
            period_t =[];
            for nn =trials_ok
                period = 'prep';
                t = (ecog.go_time(nn)-ecog.prep_time(nn))/2;
%                 period_t = [period_t ecog.prep_time(nn)+ (t - epoch/2): ecog.prep_time(nn)+ (t + epoch/2)];
                period_t = [period_t ecog.prep_time(nn)+ 500: ecog.go_time(nn)-500];

            end
        elseif tt== 4 %mvt period
            period_t =[];
            for nn =trials_ok
%                 period_t = [period_t ecog.move_time(nn): ecog.move_time(nn)+  epoch];
                period_t = [period_t ecog.move_time(nn): ecog.rest_time(nn+1)-500];

                period = 'mvt';
            end
        end
        Phase_period = PhaseFreqTransformed(:, period_t);
        Amp_period = AmpFreqTransformed(:, period_t);
    else
        period_t=[1:length(ecog.contact_pair(M1).remontaged_ecog_signal)];
        Phase_period = PhaseFreqTransformed(:,ecog.rest_time(1):ecog.end_time(end));
        Amp_period = AmpFreqTransformed(:,ecog.rest_time(1):ecog.end_time(end));
        period = 'all';
    end
    
    %Compute PSD
    [psd,f] = pwelch(lfp(period_t), WINDOW, NOVERLAP, NFFT, Fs);
    psd_all = [psd_all; psd'];
    %Compute mean accel
    accel_all(tt)=mean(abs(aux.chan(2).raw(period_t(1:end-1))));
    emg_all(tt)=mean(abs(emg.chan(1).raw(period_t(1:end-1))));
    
    counter1=0;
    for ii=1:length(PhaseFreqVector)
        counter1=counter1+1;
        ii
        counter2=0;
        for jj=1:length(AmpFreqVector)
            counter2=counter2+1;
            
            [MI,MeanAmp]=ModIndex_v2(Phase_period(ii, :), Amp_period(jj, :), position);
            Comodulogram(counter1,counter2,tt)=MI;
            x = 10:20:720;
            [val,pos]=max(MeanAmp);
            Comodulogram_phase(counter1,counter2,tt) = x(pos);
           
%             pac_value = abs(sum(exp(1i * (Phase(ii, :) - PhaseAmp(jj, :))), 'double')) / length(Phase); % pac calculation
%             pac(counter1,counter2,tt) = pac_value;
            
            if surr==0
                Comodulogram_surr(counter1,counter2,tt)=MI;
                Mean_surr(counter1,counter2,tt)=nan;
                Std_surr(counter1,counter2,tt)=nan;
                p_surr(counter1,counter2,tt)=  nan;
            else
                numpoints=size(Amp_period,2);
                numsurrogate=200; %% number of surrogate values to compare to actual value
                minskip=2*Fs; %% time lag must be at least this big
                maxskip=numpoints-Fs; %% time lag must be smaller than this
                skip=ceil(numpoints.*rand(numsurrogate*2,1));
                skip((skip>maxskip))=[];
                skip(skip<minskip)=[];
                skip=skip(1:numsurrogate,1); % creates vector with of time lags "tau" (the skip values) used for surrogate MIs
                surrogate_m=zeros(numsurrogate,1);
                
                for s=1:numsurrogate
                    Amp_surr =[Amp_period(jj,skip(s):end) Amp_period(jj,1:skip(s)-1)];
                    [MI_S,MeanAmp_S]=ModIndex_v2(Phase_period(ii,:), Amp_surr, position);
                    MI_surr(s) = MI_S;
                end
                
                % fit gaussian to surrogate data, uses normfit.m from MATLAB Statistics toolbox
                [surrogate_mean,surrogate_std]=normfit(MI_surr);
                Mean_surr(counter1,counter2,tt)=surrogate_mean;
                Std_surr(counter1,counter2,tt)=surrogate_std;
                
                p_surr(counter1,counter2,tt)=  prctile(MI_surr,99);
                Comodulogram_surr(counter1,counter2,tt)=MI;
                if MI<prctile(MI_surr,99)
                    Comodulogram_surr(counter1,counter2,tt)=0;
                end
            end
        end
    end
end

%% Graph comodulogram
Comodulogram_surr(1,1,:)=0.00001;
Clim2 = max(max(max(Comodulogram_surr)));
Clim1 = min(min(min(Comodulogram_surr)));

figure
hold on
for tt=1:4
    subplot(3,2,tt)
    C=squeeze(Comodulogram_surr(:,:,tt));
    contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,C',30,'lines','none')
    set(gca,'fontsize',9)
    ylabel('Amplitude Frequency (Hz)')
    xlabel('Phase Frequency (Hz)')
    colorbar
    caxis([Clim1 Clim2])
    if tt==1
        title(['all  ch=' num2str(M1)])
    elseif tt==1
        title(['rest  ch=' num2str(M1)])
    elseif tt==1
        title(['prep  ch=' num2str(M1)])
    elseif tt==1
        title(['mvt  ch=' num2str(M1)])
    end
    mean_COM_B(tt)=nanmean(nanmean(Comodulogram_surr(:,:,tt)));
end
subplot(3,2,5)
plot(f,log10(psd_all(2:4,:)))
xlim([1 200])
ylim([-1 3])
title('psd')
subplot(3,2,6)
title('accel/emg')
plot(accel_all(2:4))
hold on
plot(emg_all(2:4)/1000,'r')

%% save
save([name '_Com'], 'p_surr','Comodulogram_surr','Mean_surr','Std_surr','Comodulogram','Comodulogram_phase','lfp',...
    'MeanAmp','M1','PhaseFreqVector','PhaseFreq_BandWidth','AmpFreqVector','AmpFreq_BandWidth','psd_all','f','accel_all','emg_all');
saveas(gcf,[name '_Com'],'fig');

close all


