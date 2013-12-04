function nr_ipad_htkmatref(htkmatconv_dir,cond,chan1,chan2,chan3,chan4,chan5,chan6,chan7,chan8,...
                           M1_ch1,M1_ch2,n_trials,json_filename,htkmatref_filedir)
                           
                       
                       
%% Parse subject ID and append 'sum' to filename
if ~isempty(strfind(htkmatconv_dir,'ps_'))
    htkmatconv_dir_ps = strfind(htkmatconv_dir,'ps_');
    sbj = htkmatconv_dir(htkmatconv_dir_ps(1):htkmatconv_dir_ps(1)+9);
elseif ~isempty(strfind(htkmatconv_dir,'ec_'))
    htkmatconv_dir_ec = strfind(htkmatconv_dir,'ec_');
    sbj = htkmatconv_dir(htkmatconv_dir_ec(1):htkmatconv_dir_ec(1)+9);
end

htkmatconv_dir_slash = strfind(htkmatconv_dir,'/');
find_dir = htkmatconv_dir_slash(end);
htkmatconv_sum_filename = ['dat_',sbj,'_ipad_',cond];
                                           
%% load the data and store processed Aux chan data into structure array
cd(htkmatconv_dir)

ecog = struct('chan',{});
aux = struct('chan',{});
aux(1).chan(1).raw=[];
aux(1).chan(2).raw=[];
aux(1).chan(3).raw=[];
emg = struct('chan',{});
X = struct('chan',{});
signal = struct('chan',{});

gain=1e6;
Fs=1000; % sampling freq after downsampling

files = dir('*.htk');

for abc =1 : length(files)
    if ~isempty(strfind(files(abc).name,'htk'))
        % load the data and convert in mat.file
        file_name=files(abc).name;
        [d,fs,dt,tc,t]=nr_readhtk(file_name);
       
        %store the downsampled data in structure arrays
        if ~isempty(strfind(file_name,'ipad.htk'))
            aux(1).chan(1).raw=resample(d,2^10,(5^5)*8);
        elseif ~isempty(strfind(file_name,'accel.htk'))
            aux(1).chan(2).raw=resample(d,2^10,(5^5)*8);
        elseif ~isempty(strfind(file_name,'trig.htk'))
            aux(1).chan(3).raw=resample(d,2^10,(5^5)*8);
        elseif ~isempty(strfind(file_name,'emg.htk'))
            emg(1).chan(1).raw=resample(d,2^10,(5^5)*8);
        elseif ~isempty(strfind(file_name,'emg2.htk'))
            emg(1).chan(2).raw=resample(d,2^10,(5^5)*8);
        elseif ~isempty(strfind(file_name,'signal.htk'))
            signal(1).chan(1).raw=resample(d,2^10,(5^5)*8);
        else            
            ecog(1).contact_pair(tc).raw_ecog_signal=resample(d,2^10,5^5)*gain;
        end
        abc
        length(ecog(1).contact_pair)
    end
end

% Include this for epilepsy cases where there is no emg
if isempty(emg)
    emg(1).chan = [];
end


C2 = 1e6;   % constant to convert ecog/lfp channel voltage level from V->microV to match GL4k system
Glfp = 25000;   % LFP channel gain
if ~isempty(signal)
    ecog(1).contact_pair(29).raw_ecog_signal=signal(1).chan(1).raw*C2/Glfp;
end

%% notch filter around 60Hz, 120Hz and 180Hz
% butterworth notch filter - model order, [low/(Fs/2) high/(Fs/2)]
[n1_b, n1_a]=butter(3,2*[57 63]/Fs,'stop'); %60hz
[n2_b, n2_a]=butter(3,2*[117 123]/Fs,'stop'); %120hz
[n3_b, n3_a]=butter(3,2*[177 183]/Fs,'stop'); %180hz
[n4_b, n4_a]=butter(3,2*[237 243]/Fs,'stop'); %240hz
[n5_b, n5_a]=butter(3,2*[297 303]/Fs,'stop'); %300hz
[n6_b, n6_a]=butter(3,2*[357 363]/Fs,'stop'); %360hz
%     [n7_b, n7_a]=butter(3,2*[88 94]/Fs,'stop'); %180hz
for k=1:length(ecog(1).contact_pair)
    ecog(1).contact_pair(k).raw_ecog_signal=filtfilt(n1_b, n1_a, ecog(1).contact_pair(k).raw_ecog_signal); %notch out at 60
    ecog(1).contact_pair(k).raw_ecog_signal=filtfilt(n2_b, n2_a, ecog(1).contact_pair(k).raw_ecog_signal); %notch out at 120
    ecog(1).contact_pair(k).raw_ecog_signal=filtfilt(n3_b, n3_a, ecog(1).contact_pair(k).raw_ecog_signal); %notch out at 180
    ecog(1).contact_pair(k).raw_ecog_signal=filtfilt(n4_b, n4_a, ecog(1).contact_pair(k).raw_ecog_signal); %notch out at 60
    ecog(1).contact_pair(k).raw_ecog_signal=filtfilt(n5_b, n5_a, ecog(1).contact_pair(k).raw_ecog_signal); %notch out at 120
    ecog(1).contact_pair(k).raw_ecog_signal=filtfilt(n6_b, n6_a, ecog(1).contact_pair(k).raw_ecog_signal); %notch out at 180
    %          ecog(1).contact_pair(k).raw_ecog_signal=filtfilt(n7_b, n7_a, ecog(1).contact_pair(k).raw_ecog_signal); %notch out at 180
end
%% remove DC offset
for k=1:length(ecog(1).contact_pair)
    ecog(1).contact_pair(k).raw_ecog_signal=ecog(1).contact_pair(k).raw_ecog_signal-mean(ecog(1).contact_pair(k).raw_ecog_signal); %notch out at 60
end
for k=1:length(aux.chan)
    aux.chan(k).raw=aux.chan(k).raw-mean(aux.chan(k).raw);
end
% %% contact to remove
% 
% unused=input('unused contacts');
% used = setdiff(1:length(ecog.contact_pair),unused);
% for i = 1:length(used)
%     j=used(i);
%     X(1).contact_pair(i).raw_ecog_signal=ecog.contact_pair(j).raw_ecog_signal;
% end
% ecog=X;

%% detection of bad electrodes
% figure;hold on
% for i = 1: 14
%     subplot(3,5,i)
%     plot(ecog.contact_pair(i).raw_ecog_signal)
%     title(num2str(i))
% %     ylim([-0.005 0.005])
% end
% figure;hold on
% for i = 15: 28
%     subplot(3,5,i-14)
%     plot(ecog.contact_pair(i).raw_ecog_signal)
%     title(num2str(i))
% %     ylim([-25 25])
% end
% bad=input('bad contacts');
% 
% close all
%% common reference
if length(ecog.contact_pair) <=30
     length_CAR = 28;
elseif length(ecog.contact_pair) <=64 
    length_CAR = 64;
else
    length_CAR = 256;
end
data = nan*ones(length_CAR,length(ecog.contact_pair(1).raw_ecog_signal));
for i = 1: length_CAR
    data(i,:) = ecog.contact_pair(i).raw_ecog_signal';
end
car=nanmean(data)/length_CAR;
for i = 1: length_CAR
    ecog.contact_pair(i).remontaged_ecog_signal = ecog.contact_pair(i).raw_ecog_signal-car;
%     ecog.contact_pair(i).remontaged_ecog_signal=ecog.contact_pair(i).raw_ecog_signal-ecog.contact_pair(i+1).raw_ecog_signal;
end
% ecog.contact_pair(28).remontaged_ecog_signal=ecog.contact_pair(28).raw_ecog_signal;
if size(ecog.contact_pair,2) ==29
    ecog.contact_pair(29).remontaged_ecog_signal=ecog.contact_pair(29).raw_ecog_signal;
elseif size(ecog.contact_pair,2) >28 && size(ecog.contact_pair,2) <=32
    for ii = 29:size(ecog.contact_pair,2)-1
        ecog.contact_pair(ii).remontaged_ecog_signal=ecog.contact_pair(ii).raw_ecog_signal-ecog.contact_pair(ii+1).raw_ecog_signal;
    end
end
%% Label channels

[auxchan_lbl] = nr_htk_auxchan_lbl(cond, chan1, chan2, chan3, chan4, chan5, chan6, chan7, chan8)

%% load ipad data

%if ~isempty(strfind(filename,'_ipad')) %evt == 1
    
%     [d,fs,dt,tc,t]=readhtk('trig.htk');
%     chan_trig=d;
%     % find active and rest periods using the raw signal
%     MARGIN = 0.005; % is a very wide range that handles large fluctuations in task voltage. ParseIntraOpEvents uses MARGIN = 0.05
%     
% %     thresh = max(chan_trig);
% %     ACTIVE = [thresh-MARGIN thresh+MARGIN]; %Determine the threshold
% %     in1 = inrange(chan_trig,ACTIVE);
%     thresh =  4;%mean(chan_trig(find(chan_trig~=0)))
%     in1=chan_trig>thresh;
%     inds = find(in1);
%     [pos,n] = evFindGroups(inds,1,1000); %find active period of 1s minimum
%     active_time = inds(pos(1,:))/(fs/1000);
%     rest_time = inds(pos(2,:))/(fs/1000);
%     
%     % save rest/active timestamps in ecog structure array
%     ecog(1).rest_time=int32(rest_time);
%     ecog(1).active_time=int32(active_time);
    
    % find the begining of each trial using the raw signal
    [d,fs,dt,tc,t]=nr_readhtk('ipad.htk');
    chan_ipad=d;
    plot(d)
    hold on
    thresh =  input('thresh');
    START =  15;% input('start');
    inds = find(chan_ipad(START:end-15)>=thresh);
    [pos,n] = nr_evFindGroups(inds,500,1);
    plot(inds(pos(1,1:end-1))+START,thresh, '*r')
    plot(inds(pos(1,2:end))+START,thresh-0.2*thresh, '*k')
    
    ipad_ON_time = (inds(pos(1,1:end-1))+START)/(fs/1000);
    ipad_OFF_time = (inds(pos(1,2:end))+START)/(fs/1000);
    DELETE =  input('delete');
    GOOD = setdiff([1:length(ipad_ON_time)],DELETE)
    ecog(1).ipad_ON_time=int32(ipad_ON_time(GOOD));
    ecog(1).ipad_OFF_time= int32(ipad_OFF_time(GOOD));
    
    
%% load ecog traces into matrix for eegplot

for i = 1:size(ecog.contact_pair,2)
    ecog_traces(i,:) = ecog.contact_pair(i).remontaged_ecog_signal;
end


%% Parse json file for timestamps

%[description, timestamp] = nr_import_from_ipad_02(json_filename);


n_steps = description{2};
n_steps = str2num(n_steps(length(n_steps)))-1;

% n_trials = description{end-1};
% n_trials = str2num(n_trials(length(n_trials)))+1;
% 
% n_trials_test = ipad.data.description{end-1}
% n_trials_test(length(n_trials_test))

%n_trials = 15;
% Sound
[sound]=strfind(description,'Sound');
sound = cellfun('isempty',sound);
sound_ON_idx = find(sound==0);
sound_ON = timestamp(sound_ON_idx)';
%assignin('base','sound_ON',sound_ON)

% trial start
[trial]=strfind(description,'Trial');
trial = cellfun('isempty',trial);
trial_ON_idx = find(trial==0);
trial_ON = timestamp(trial_ON_idx)';
%assigning('base','trial_ON',trial_ON)

% Rest epoch Beg Fixation point ON red dot
[rest]=strfind(description,'Scene: 0 - Rest');
rest = cellfun('isempty',rest);
rest_ON_idx = find(rest==0);
rest_ON = timestamp(rest_ON_idx)';
%assigning('base','rest_ON',rest_ON)

% Rest epoch End Fixation point OFF red dot
[rest]=strfind(description,'Timer: 0');
rest = cellfun('isempty',rest);
rest_OFF_idx = find(rest==0);
rest_OFF = timestamp(rest_OFF_idx)';

% Rest epoch Error Fixation error by mvt
[rest]=strfind(description,'Incorrect touch: 0');
rest = cellfun('isempty',rest);
rest_error_idx = find(rest==0);
rest_error = timestamp(rest_error_idx)';

% Preparation epoch Beg ON Cue ON blue dot
[prep]=strfind(description,'Scene: 1 - Preparation');
prep = cellfun('isempty',prep);
prep_ON_idx = find(prep==0);
prep_ON = timestamp(prep_ON_idx)';

% Preparation epoch End Cue OFF blue dot
[prep]=strfind(description,'Timer: 1');
prep = cellfun('isempty',prep);
prep_OFF_idx = find(prep==0);
prep_OFF = timestamp(prep_OFF_idx)';

% Preparation epoch error by mvt
[prep]=strfind(description,'Incorrect touch: 1');
prep = cellfun('isempty',prep);
prep_error_idx = find(prep==0);
prep_error = timestamp(prep_error_idx)';

% Movements
target_ON = nan*ones(n_steps,n_trials)';
target_idx = nan*ones(n_steps,n_trials)';

touch = nan*ones(n_steps,n_trials)';
notouch = nan*ones(n_steps,n_trials)';
incorrect_touch = nan*ones(n_steps,n_trials)';
touch = nan*ones(n_steps,n_trials)';

% Find 'Scene : movement'
[target]=strfind(description,'Movement');
target = cellfun('isempty',target);
target = find(target==0);
target = target(2:end);

for i = 1:n_steps
    %i
    % Target1 ON
    target_idx(:,i) = target(i:n_steps:end);
    target_ON(:,i) = timestamp(target(i:n_steps:end));
end

% Find the touch time
touch_idx = (target_idx) +1;
touch_ON = timestamp(touch_idx);

% Estimate ipad beeps
ecog.ipad_ON_time(end+1) = ecog.ipad_OFF_time(end);
ipad_beeps = ecog.ipad_ON_time';
ipad_sound_diff = ipad_beeps(1:length(sound_ON)) - int32(sound_ON);
mean_ipad_sound_diff = mean(ipad_sound_diff);
ipad_sound_shift = sound_ON + mean_ipad_sound_diff;
ipad_prep_ON_shift = prep_ON + mean_ipad_sound_diff;
ipad_prep_OFF_shift = prep_OFF + mean_ipad_sound_diff;
ipad_prep_error_shift = prep_error + mean_ipad_sound_diff;
ipad_rest_error_shift = rest_error + mean_ipad_sound_diff;
ipad_touch_ON_shift = touch_ON + mean_ipad_sound_diff;

ipad.data.description = description;
ipad.data.timestamps = timestamp;
ipad.events.prep_on = prep_ON;
ipad.events.prep_off = prep_OFF;
ipad.events.prep_error = prep_error;
ipad.events.rest_on = rest_ON;
ipad.events.rest_off = rest_OFF;
ipad.events.rest_error = rest_error;
ipad.events.sound_on = sound_ON;
ipad.events.touch_ON = touch_ON;
ipad.events_sc.factors.vec = ipad_sound_diff;
ipad.events_sc.factors.mean = mean_ipad_sound_diff;
ipad.events_sc.sound_ON_sc = ipad_sound_shift;
ipad.events_sc.prep_ON_sc = ipad_prep_ON_shift;
ipad.events_sc.prep_OFF_sc = ipad_prep_OFF_shift;
ipad.events_sc.prep_error_sc = ipad_prep_error_shift;
ipad.events_sc.rest_error_sc = ipad_rest_error_shift;
ipad.events_sc.touch_ON_sc = ipad_touch_ON_shift;

%% save data

cd(htkmatref_filedir)
save(htkmatconv_sum_filename,'auxchan_lbl','M1_ch1','M1_ch2','n_trials','aux','emg','ecog','Fs','thresh','ipad','GOOD','ecog_traces');


