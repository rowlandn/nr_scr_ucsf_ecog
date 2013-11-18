
function Comodulogram_ipad_time(name,M1)






load(name)

%name = strrep(name,'_ecog.mat','');
%% General
srate=Fs;

epoch = 1*Fs; % duration of ecah period
ecog.move_time=ecog.active_time;
ecog.move_off_time=ecog.rest_time;

%% data
if M1<5
    lfp = ecog.contact_pair(M1).raw_ecog_signal - ecog.contact_pair(M1+1).raw_ecog_signal;
else
    lfp = ecog.contact_pair(M1).raw_ecog_signal;
end

st_r = ecog.trial_beg_time;
st_p = ecog.prep_time;
st_a = ecog.active_time;
st_e = ecog.rest_time+500; % add 500 to be sure that patient is resting.
st_off = ecog.trial_end_time;

%% subdivise file into rest-prep-mvt or xs periods

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
    AmpFreq=eegfilt(lfp,srate,Af1,Af2); % just filtering
    AmpFreqTransformed(ii, :) = abs(hilbert(AmpFreq)); % getting the amplitude envelope
end

for jj=1:length(PhaseFreqVector)
    Pf1 = PhaseFreqVector(jj) - PhaseFreq_BandWidth/2;
    Pf2 = PhaseFreqVector(jj) + PhaseFreq_BandWidth/2;
    PhaseFreq=eegfilt(lfp,srate,Pf1,Pf2); % this is just filtering
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
            
            [MI,MeanAmp]=ModIndex_v2(PhaseFreqTransformed(ii, aa:aa+epoch), AmpFreqTransformed(jj, aa:aa+epoch), position);
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
saveas(gcf,[name '_Com_time'],'fig');
%% test
[P,ANOVATAB,STATS] = kruskalwallis(MI_means,pts(:,3));
title(P)
saveas(gcf,[name '_KW_time'],'fig');
save([name '_Com_time'], 'Comodulogram_surr','Mean_surr','Std_surr','Comodulogram','Comodulogram_phase','lfp','MeanAmp','M1','PhaseFreqVector','PhaseFreq_BandWidth','AmpFreqVector','AmpFreq_BandWidth', 'MI_means','Phase_means','pts');


