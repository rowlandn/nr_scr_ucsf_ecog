function nr_comodulogram_ipad_time(sum_filedirnm,condition,chand)


load(sum_filedirnm)

if ~isempty(strfind(sum_filedirnm,'ps_'))
    sum_filedirnm_ps = strfind(sum_filedirnm,'ps_');
    sbj = sum_filedirnm(sum_filedirnm_ps(1):sum_filedirnm_ps(1)+9);
elseif ~isempty(strfind(sum_filedirnm,'ec_'))
    sum_filedirnm_ec = strfind(sum_filedirnm,'ec_');
    sbj = sum_filedirnm(sum_filedirnm_ec(1):sum_filedirnm_ec(1)+9);
end

%name = strrep(name,'_ecog.mat','');
%% General
srate=Fs;

epoch = 1*Fs; % duration of each period
ecog.move_time=ecog.active_time;
ecog.move_off_time=ecog.rest_time;


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


%% data
if chan<5
    lfp = ecog.contact_pair(chan).raw_ecog_signal - ecog.contact_pair(chan+1).raw_ecog_signal;
else
    lfp = ecog.contact_pair(chan).raw_ecog_signal;
end

st_r = ecog.trial_beg_time;
st_p = ecog.prep_time;
st_a = ecog.active_time;
st_e = ecog.rest_time+500; % add 500 to be sure that patient is resting.
st_off = ecog.trial_end_time;

%% subdivide file into rest-prep-mvt or xs periods

st_a = fix(st_a);
st_r = fix(st_r);
if st_a(1)==0
    st_a(1)=1;
end
if st_r(1)==0
    st_r(1)=1;
end

stim=zeros(length(lfp),1);

for k=1:length(st_r)
    stim(st_r(k):st_p(k))=1;
end
for k=1:length(st_r)   
    stim(st_p(k)+1:st_a(k))=2;
    stim(st_a(k)+1:st_e(k))=3;
end
% stim(1:st_r(1))=0;
% stim(st_off(end):end)=0;
% delete bad trials
bad_trials = setdiff(1:size(Time_events,1),trials_ok);
for i= 1:length(bad_trials)
    x=bad_trials(i);
    if x< size(Time_events,1)
    stim(st_r(x):st_r(x+1))=0;
    else
        stim(st_r(x):length(stim))=0;
    end
end
ns=find((stim-[2; stim(1:(end-1))])~=0);
pts=[];
for k = 1:length(ns)
    k
    if k<length(ns)
        n = floor((ns(k+1)-ns(k))/epoch);
        t = [1:n]*epoch;
        xx = [ns(k) ns(k)+t];
        pts = [pts xx];
    else
        n = floor((length(lfp)-ns(k))/epoch);
        t = [1:n]*epoch;
        xx = [ns(k) ns(k)+t];
        pts = [pts xx];        
    end
end
pts =pts';
pts=[pts pts stim(pts(:,1))];
pts = pts(find(pts(:,3)~=0),:);
good=find(pts(2:end,1)-pts(1:end-1,1)>=epoch);
pts=pts(good,:);


%% Define the Amplitude- and Phase- Frequencies

PhaseFreqVector=[4:2:50];%[12:2:30];
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
AmpFreqTransformed = zeros(length(AmpFreqVector), length(lfp));
PhaseFreqTransformed = zeros(length(PhaseFreqVector), length(lfp));

for ii=1:length(AmpFreqVector)
    Af1 = AmpFreqVector(ii)-AmpFreq_BandWidth/2;
    Af2=AmpFreqVector(ii)+AmpFreq_BandWidth/2;
    AmpFreq=nr_eegfilt(lfp,srate,Af1,Af2); % just filtering
    AmpFreqTransformed(ii, :) = abs(hilbert(AmpFreq)); % getting the amplitude envelope
end

for jj=1:length(PhaseFreqVector)
    Pf1 = PhaseFreqVector(jj) - PhaseFreq_BandWidth/2;
    Pf2 = PhaseFreqVector(jj) + PhaseFreq_BandWidth/2;
    PhaseFreq=nr_eegfilt(lfp,srate,Pf1,Pf2); % this is just filtering
    PhaseFreqTransformed(jj, :) = angle(hilbert(PhaseFreq)); % this is getting the phase time series
end
toc

%% Do comodulation calculation
'Comodulation loop'

MI_means = [];
Phase_means = [];
counter0=0;
for kk = 1:length(pts(:,1))
    aa = pts(kk,1);
    counter0=counter0+1;
    
    counter1=0;
    for ii=1:length(PhaseFreqVector)
        counter1=counter1+1;
        ii
        counter2=0;
        for jj=1:length(AmpFreqVector)
            counter2=counter2+1;
            
            [MI,MeanAmp]=nr_ModIndex_v2(PhaseFreqTransformed(ii, aa:aa+epoch), AmpFreqTransformed(jj, aa:aa+epoch), position);
            Comodulogram(counter0,counter1,counter2)=MI;
            x = 10:20:720;
            [val,pos]=max(MeanAmp);
            Comodulogram_phase(counter0,counter1,counter2) = x(pos);
            Comodulogram_surr(counter0,counter1,counter2)=MI;
            Mean_surr(counter0,counter1,counter2)=nan;
            Std_surr(counter0,counter1,counter2)=nan;
%                             numpoints=size(AmpFreqTransformed_down,2);
%                             numsurrogate=50; %% number of surrogate values to compare to actual value
%                             minskip=Fs; %% time lag must be at least this big
%                             maxskip=numpoints-Fs; %% time lag must be smaller than this
%                             skip=ceil(numpoints.*rand(numsurrogate*2,1));
%                             skip((skip>maxskip))=[];
%                             skip(skip<minskip)=[];
%                             skip=skip(1:numsurrogate,1); % creates vector with of time lags "tau" (the skip values) used for surrogate MIs
%                             surrogate_m=zeros(numsurrogate,1);
%             
%                             for s=1:numsurrogate
%                                 Amp_surr =[AmpFreqTransformed_down(jj,skip(s):end) AmpFreqTransformed_down(jj,1:skip(s)-1)];
%                                 [MI_S,MeanAmp_S]=ModIndex_v2(PhaseFreqTransformed_down(ii,:), Amp_surr, position);
%                                 MI_surr(s) = MI_S;
%                             end
%             
%                             % fit gaussian to surrogate data, uses normfit.m from MATLAB Statistics toolbox
%                             [surrogate_mean,surrogate_std]=normfit(MI_surr);
%                             Mean_surr(counter1,counter2)=surrogate_mean;
%                             Std_surr(counter1,counter2)=surrogate_std;
%                             Comodulogram_surr(counter1,counter2)=(abs(Comodulogram(counter1,counter2))-surrogate_mean)/surrogate_std;
%             
        end
    end
    toc
    MI_means = [MI_means; mean(mean(Comodulogram_surr(kk,4:end,:)))];
    Phase_means = [Phase_means; mean(mean(Comodulogram_phase(kk,4:end,:)))];
end
    %% Graph comodulogram
figure
plot([1:length(MI_means)],MI_means,'k*')
x=find(pts(:,3)==2);
hold on
plot(x,MI_means(x),'r*')
x=find(pts(:,3)==3);
plot(x,MI_means(x),'g*')
%saveas(gcf,[sum_filedirnm '_Com_time'],'fig');
%% test
[P,ANOVATAB,STATS] = kruskalwallis(MI_means,pts(:,3));
title(P)
%saveas(gcf,[sum_filedirnm '_KW_time'],'fig');
%save(['com_time_' sbj], 'Comodulogram_surr','Mean_surr','Std_surr','Comodulogram','Comodulogram_phase','lfp','MeanAmp','chan','PhaseFreqVector','PhaseFreq_BandWidth','AmpFreqVector','AmpFreq_BandWidth', 'MI_means','Phase_means','pts','P');

str_eval1 = [chand,'.Comodulogram_surr = Comodulogram_surr;'];
eval(str_eval1)
str_eval2 = [chand,'.Mean_surr = Mean_surr;'];
eval(str_eval2)
str_eval3 = [chand,'.Std_surr = Std_surr;'];
eval(str_eval3)
str_eval4 = [chand,'.Comodulogram = Comodulogram;'];
eval(str_eval4)
str_eval5 = [chand,'.Comodulogram_phase = Comodulogram_phase;'];
eval(str_eval5)
str_eval6 = [chand,'.lfp = lfp;'];
eval(str_eval6)
str_eval7 = [chand,'.MeanAmp = MeanAmp;'];
eval(str_eval7)
str_eval8 = [chand,'.chan = chan;'];
eval(str_eval8)
str_eval9 = [chand,'.PhaseFreqVector = PhaseFreqVector;'];
eval(str_eval9)
str_eval10 = [chand,'.PhaseFreq_BandWidth = PhaseFreq_BandWidth;'];
eval(str_eval10)
str_eval11 = [chand,'.AmpFreqVector = AmpFreqVector;'];
eval(str_eval11)
str_eval12 = [chand,'.AmpFreq_BandWidth = AmpFreq_BandWidth;'];
eval(str_eval12)
str_eval13 = [chand,'.MI_means = MI_means;'];
eval(str_eval13)
str_eval14 = [chand,'.Phase_means = Phase_means;'];
eval(str_eval14)
str_eval15 = [chand,'.pts = pts;'];
eval(str_eval15)
str_eval16 = [chand,'.P = P;'];
eval(str_eval16)

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