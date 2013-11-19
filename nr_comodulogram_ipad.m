function nr_comodulogram_ipad(sum_filedirnm,condition,chand,surr) 

load(sum_filedirnm)

if ~isempty(strfind(sum_filedirnm,'ps_'))
    sum_filedirnm_ps = strfind(sum_filedirnm,'ps_');
    sbj = sum_filedirnm(sum_filedirnm_ps(1):sum_filedirnm_ps(1)+9);
elseif ~isempty(strfind(sum_filedirnm,'ec_'))
    sum_filedirnm_ec = strfind(sum_filedirnm,'ec_');
    sbj = sum_filedirnm(sum_filedirnm_ec(1):sum_filedirnm_ec(1)+9);
end

if chand == 'M1_ch1'
    chan = M1_ch1;
elseif chand == 'M1_ch2'
    if exist('M1_ch2')
        chan = M1_ch2;
    else
        chan = M1_ch1;
    end
elseif chand == 'S1_ch1'
    if exist('M1_ch2')
        chan = M1_ch1-2;
    else
        chan = M1_ch1-1;
    end
elseif chand == 'S1_ch2'
    if exist('M1_ch2')
        chan = M1_ch2-2;
    else
        chan = M1_ch1-1;
    end
elseif chand == 'P1_ch1'
    if M1_ch1 < 5
        if exist('M1_ch2')
            chan = M1_ch1+2;
        else
            chan = M1_ch1+1;
        end
    else
        if exist('M1_ch2')
            chan = M1_ch1+2;
        else
            chan = M1_ch1;   % remember, here you are setting P1 to M1
        end
    end
elseif chand == 'P1_ch2'
    if M1_ch1 < 5
        if exist('M1_ch2')
            chan = M1_ch2+2;
        else
            chan = M1_ch1+1;
        end
    else
        if exist('M1_ch2')
            chan = M1_ch2+2;
        else
            chan = M1_ch1;   % remember, here you are setting P1 to M1
        end
    end
end




srate=Fs;
epoch = 2*Fs; % 2 sec duration
WINDOW = 1024;
NOVERLAP = 512;
NFFT = 1024;
lfp = ecog.contact_pair(chan).remontaged_ecog_signal;
lfp = lfp-mean(lfp);

% remove artifacts
[n1_b, n1_a]=butter(3,2*[118 122]/Fs,'stop'); %120hz
lfp=filtfilt(n1_b, n1_a, lfp);
[n1_b, n1_a]=butter(3,2*[178 182]/Fs,'stop'); %180hz
lfp=filtfilt(n1_b, n1_a, lfp);
if ~isempty(strfind(sbj,'DBSON'))
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
    AmpFreq=nr_eegfilt(lfp,srate,Af1,Af2); % just filtering
    AmpFreqTransformed(ii, :) = abs(hilbert(AmpFreq)); % getting the amplitude envelope
end

for jj=1:length(PhaseFreqVector)
    Pf1 = PhaseFreqVector(jj)- PhaseFreq_BandWidth/2;
    Pf2 = PhaseFreqVector(jj) + PhaseFreq_BandWidth/2;
    PhaseFreq=nr_eegfilt(lfp,srate,Pf1,Pf2); % this is just filtering
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
        period_t=[1:length(ecog.contact_pair(chan).remontaged_ecog_signal)];
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
            
            [MI,MeanAmp]=nr_ModIndex_v2(Phase_period(ii, :), Amp_period(jj, :), position);
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

% % %% Graph comodulogram
% % Comodulogram_surr(1,1,:)=0.00001;
% % Clim2 = max(max(max(Comodulogram_surr)));
% % Clim1 = min(min(min(Comodulogram_surr)));
% % 
% % figure
% % hold on
% % for tt=1:4
% %     subplot(3,2,tt)
% %     C=squeeze(Comodulogram_surr(:,:,tt));
% %     contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,C',30,'lines','none')
% %     set(gca,'fontsize',9)
% %     ylabel('Amplitude Frequency (Hz)')
% %     xlabel('Phase Frequency (Hz)')
% %     colorbar
% %     caxis([Clim1 Clim2])
% %     if tt==1
% %         title(['all  ch=' num2str(chan)])
% %     elseif tt==1
% %         title(['rest  ch=' num2str(chan)])
% %     elseif tt==1
% %         title(['prep  ch=' num2str(chan)])
% %     elseif tt==1
% %         title(['mvt  ch=' num2str(chan)])
% %     end
% %     mean_COM_B(tt)=nanmean(nanmean(Comodulogram_surr(:,:,tt)));
% % end
% % subplot(3,2,5)
% % plot(f,log10(psd_all(2:4,:)))
% % xlim([1 200])
% % ylim([-1 3])
% % title('psd')
% % subplot(3,2,6)
% % title('accel/emg')
% % plot(accel_all(2:4))
% % hold on
% % plot(emg_all(2:4)/1000,'r')

%% save


str_eval1 = [chand,'.p_surr = p_surr;'];
eval(str_eval1)
str_eval2 = [chand,'.Comodulogram_surr = Comodulogram_surr;'];
eval(str_eval2)
str_eval3 = [chand,'.Mean_surr = Mean_surr;'];
eval(str_eval3)
str_eval4 = [chand,'.Std_surr = Std_surr;'];
eval(str_eval4)
str_eval5 = [chand,'.Comodulogram = Comodulogram;'];
eval(str_eval5)
str_eval6 = [chand,'.Comodulogram_phase = Comodulogram_phase;'];
eval(str_eval6)
str_eval7 = [chand,'.lfp = lfp;'];
eval(str_eval7)
str_eval8 = [chand,'.MeanAmp = MeanAmp;'];
eval(str_eval8)
str_eval9 = [chand,'.chan = chan;'];
eval(str_eval9)
str_eval0 = [chand,'.PhaseFreqVector = PhaseFreqVector;'];
eval(str_eval0)
str_eval11 = [chand,'.PhaseFreq_BandWidth = PhaseFreq_BandWidth;'];
eval(str_eval11)
str_eval12 = [chand,'.AmpFreqVector = AmpFreqVector;'];
eval(str_eval12)
str_eval13 = [chand,'.AmpFreq_BandWidth = AmpFreq_BandWidth;'];
eval(str_eval13)
str_eval14 = [chand,'.psd_all = psd_all;'];
eval(str_eval14)
str_eval15 = [chand,'.f = f;'];
eval(str_eval15)
str_eval16 = [chand,'.accel_all = accel_all;'];
eval(str_eval16)
str_eval17 = [chand,'.emg_all = emg_all;'];
eval(str_eval17)

mfnm = mfilename;

if ~exist(['com_',sbj,'.mat'],'file')
    allt = struct('prelead',[],'postlead',[]);
    tres = struct('prelead',[],'postlead',[]);
    if strcmp(mfnm,'nr_comodulogram_ipad') == 1
        allt_eval = ['allt.',condition,'.',chand,' = ',chand,';'];
        eval(allt_eval)
        save(['com_',sbj],'allt','tres')
    elseif strcmp(mfnm,'nr_comodulogram_ipad_time') == 1
        tres_eval = ['tres.',condition,'.',chand,' = ',chand,';'];
        eval(tres_eval)
        save(['com_',sbj],'allt','tres')
    end
else
    load(['com_',sbj])
    if strcmp(mfnm,'nr_comodulogram_ipad') == 1
        allt_eval = ['allt.',condition,'.',chand,' = ',chand,';'];
        eval(allt_eval)
        save(['com_',sbj],'allt','tres')
    elseif strcmp(mfnm,'nr_comodulogram_ipad_time') == 1
        tres_eval = ['tres.',condition,'.',chand,' = ',chand,';'];
        eval(tres_eval)
        save(['com_',sbj],'allt','tres')
    end
end
    
    






