function nr_ipad_aomatref(aomatconv_filename,cond,emg1,emg2,emg3,emg4,chan1,chan2,chan3,...
                          M1_ch1,n_trials,json_filename,json_mat_conv_dir,aomatref_filedir)

% converts mat file from AlphaOmega Mapfile converter program to a format
% that can be run on existing code for ecogPSD analysis and new ipad
% analyses

% Created by SAS (12/15/2008)pl
% Edited by SAS (5/29/09) to process EMG/accel data.
% 
% Input:
%       1) "fname.mat" = output file form Mapfileconverter
% Output:
%       1) "fname_ecog.mat" = ecog/LFP data


%% load data
load(aomatconv_filename);

%% Parse subject ID and append 'sum' to filename
if ~isempty(strfind(aomatconv_filename,'ps_'))
    aomatconv_filename_ps = strfind(aomatconv_filename,'ps_');
    sbj = aomatconv_filename(aomatconv_filename_ps(1):aomatconv_filename_ps(1)+9);
elseif ~isempty(strfind(aomatconv_filename,'ec_'))
    aomatconv_filename_ec = strfind(aomatconv_filename,'ec_');
    sbj = aomatconv_filename(aomatconv_filename_ec(1):aomatconv_filename_ec(1)+9);
end

aomatconv_filename_slash = strfind(aomatconv_filename,'/');
find_dir = aomatconv_filename_slash(end);
aomatconv_sum_filename = ['sum_',sbj,'_',aomatconv_filename(find_dir+12:end-4)]; 

% make sure sampling frequencies are equal for ecog/LFP/task recordings, then
% use CECOG_1_KHz varible as the sampling rate
% if isequal(CECOG_1_KHz,CECOG_2_KHz,CECOG_3_KHz,CECOG_4_KHz,CECOG_5_KHz, CLFP_KHz,CAn_In__1_KHz,CAn_In__2_KHz,CAn_In__3_KHz,CEMG_1_KHz,CEMG_2_KHz,CEMG_3_KHz)
if exist('CECOG1_KHz')
    Fs = CECOG1_KHz;
    Fs = Fs*1000;
else
    Fs = CECOG_1_KHz;
    Fs = Fs*1000; 
end
% multiply by 1000 for KHz->Hz conversion
% else
%     error('Sampling frequencies are different for ecog/LFP data');
% end

%% Label channels

[auxchan_lbl] = nr_aom_auxchan_lbl(cond,emg1,emg2,emg3,emg4,chan1,chan2,chan3);

%% force ecog/lfp data to be the same lengths for montage-ing

minlength=min([length(CECOG_1) length(CECOG_2) length(CECOG_3)...
    length(CECOG_4) length(CECOG_5) length(CLFP)]);
CECOG_1  = CECOG_1(1:minlength);
CECOG_2  = CECOG_2(1:minlength);
CECOG_3  = CECOG_3(1:minlength);
CECOG_4  = CECOG_4(1:minlength);
CECOG_5  = CECOG_5(1:minlength);
CLFP   = CLFP(1:minlength);

%% ecog/LFP voltage calibration
C1 = 5/32768; % constant for calibration
% note: value 32768 is equal to 5V in the Alpha Omega system
C2 = 1e6;   % constant to convert ecog/lfp channel voltage level from V->microV to match GL4k system
Gecog = 7000;   % ecog channel gain
Glfp = 25000;   % LFP channel gain

CAn_In__1 = CAn_In__1*C1;
CAn_In__2 = CAn_In__2*C1;
CAn_In__3 = CAn_In__3*C1;

% CEMG_1 = CEMG_1*C1;
% CEMG_2 = CEMG_2*C1;
% CEMG_3 = CEMG_3*C1;
% 
CECOG_1  = CECOG_1*C1*C2/Gecog;
CECOG_2  = CECOG_2*C1*C2/Gecog;
CECOG_3  = CECOG_3*C1*C2/Gecog;
CECOG_4  = CECOG_4*C1*C2/Gecog;
CECOG_5  = CECOG_5*C1*C2/Gecog;
CLFP   = CLFP*C1*C2/Glfp;

%% store processed ecog/LFP data into ecog structure array

% frq1 = menu('Select one','resample1K','resample2K','no');

d = round(Fs/1000);
if int32(Fs) == 1502
    frq1 = 2;
    frq2 = 3;
    d=Fs/1000;
else
    frq1=1;
    frq2=d;
    d=Fs/1000;
    
end

% ecog = struct('contact_pair',{},'rest_time',{},'active_time',{});
% evt = menu('Select one','LFP','No LFP');

ecog(1).contact_pair(1).raw_ecog_signal=resample(CECOG_1,frq1,frq2);
ecog(1).contact_pair(2).raw_ecog_signal=resample(CECOG_2,frq1,frq2);
ecog(1).contact_pair(3).raw_ecog_signal=resample(CECOG_3,frq1,frq2);
ecog(1).contact_pair(4).raw_ecog_signal=resample(CECOG_4,frq1,frq2);
ecog(1).contact_pair(5).raw_ecog_signal=resample(CECOG_5,frq1,frq2);

if ~isempty(strfind(aomatconv_filename,'lfp')) || ~isempty(strfind(aomatconv_filename,'LFP'))% evt == 1
    ecog(1).contact_pair(6).raw_ecog_signal=resample(CLFP,frq1,frq2);
end
Fs = 1000;

%% store processed Aux chan data into aux structure array

aux = struct('chan',{});
aux(1).chan(1).raw=resample(CAn_In__1,frq1,frq2);
aux(1).chan(2).raw=resample(CAn_In__2,frq1,frq2);
aux(1).chan(3).raw=resample(CAn_In__3,frq1,frq2);
aux.lbl = auxchan_lbl.prelead;

%% store processed EMG chan data into aux structure array

emg = struct('chan',{});
emg(1).chan(1).raw=resample(CEMG_1,frq1,frq2);
emg(1).chan(2).raw=resample(CEMG_2,frq1,frq2);
emg(1).chan(3).raw=resample(CEMG_3,frq1,frq2);
emg(1).chan(4).raw=resample(CEMG_4,frq1,frq2);

%% remove zeros add when map data are converted in mat file by 'mapconverter'


%% Do the remontage
for i = 1: length(ecog.contact_pair)
    if i<5
        ecog.contact_pair(i).remontaged_ecog_signal = ecog.contact_pair(i).raw_ecog_signal - ecog.contact_pair(i+1).raw_ecog_signal;
    else
        ecog.contact_pair(i).remontaged_ecog_signal = ecog.contact_pair(i).raw_ecog_signal;
    end
end

%% grab rest/active timestamps from task button

% ask user for aux channel with task button voltage
% evt = menu('Select one','rest/active data','rest only');

if ~isempty(strfind(aomatconv_filename,'ipad')) %evt == 1
    
    % find active and rest periods using the raw signal
    MARGIN = 0.5; % is a very wide range that handles large fluctuations in task voltage. ParseIntraOpEvents uses MARGIN = 0.05
    thresh = max(CAn_In__3); 
    ACTIVE = [thresh-MARGIN thresh+MARGIN]; %Determine the threshold
    in1 = nr_inrange(CAn_In__3,ACTIVE);
    inds = find(in1);
    [pos,n] = nr_evFindGroups(inds,1,1000*d); %find active period of 1s minimum
    active_time = inds(pos(1,:))/d;
    rest_time = inds(pos(2,:))/d;
    
    %assignin('base','thresh',thresh)
    %%assignin('base','active',active)
    %assignin('base','in1',in1)
    %%assignin('base','inds',inds)
    %assignin('base','active_time',active_time)
    %assignin('base','rest_time',rest_time)
        
    % save rest/active timestamps in ecog structure array
    ecog(1).rest_time=int32(rest_time);
    ecog(1).active_time=int32(active_time);
    
    % find the begining of each trial using the raw signal
    CAn_In__1=CAn_In__1-mean(CAn_In__1);
%     thresh = max(CAn_In__1);
%     MARGIN = thresh/5; % 0.3 is a very wide range that handles large fluctuations in task voltage. ParseIntraOpEvents uses MARGIN = 0.05
%     ACTIVE = [thresh-MARGIN thresh+MARGIN];
%     in1 = inrange(CAn_In__1,ACTIVE);
%     inds = find(in1);
%     thresh = max(CAn_In__1)-3*std(CAn_In__1);
    
    close ;plot(CAn_In__1)
    hold on
    thresh = input('thresh = ');
    
    %thresh =  0.02 %input('thresh');
    
    %assignin('base','CAn_In__1',CAn_In__1)
    %assignin('base','thresh',thresh)
    
     
    inds = find(abs(CAn_In__1)>thresh);
    
    %assignin('base','inds',inds)
    
    [pos,n] = nr_evFindGroups(inds,5*d,2);
    ipad_ON_time = inds(pos(1,1:end-1))/d;
    ipad_OFF_time = inds(pos(1,2:end))/d;
    % delete false detection
    diff_ipad_ON=  ipad_ON_time(2:end)-ipad_ON_time(1:end-1);
    ok_ON = find(diff_ipad_ON>5000);
    
    %%assignin('base','evFindGroups',evFindGroups)
    %assignin('base','d',d)
    %assignin('base','thresh',thresh)
    %%assignin('base','active',active)
    %assignin('base','in1',in1)
    %%assignin('base','inds',inds)
    %assignin('base','active_time',active_time)
    %assignin('base','rest_time',rest_time)
    
    
    
    
    
    %assignin('base','ipad_ON_time',ipad_ON_time)
    %assignin('base','diff_ipad_ON',diff_ipad_ON)
    %assignin('base','ok_ON',ok_ON)
    
    diff_ipad_OFF = ipad_OFF_time(2:end)-ipad_OFF_time(1:end-1);
    %assignin('base','diff_ipad_OFF',diff_ipad_OFF)
    %ok_OFF = find(diff_ipad_OFF>5000);
    % check the timing
    plot(ipad_ON_time*d, thresh,'*r')
    plot(ipad_OFF_time*d, thresh-0.2*thresh,'*k')
     
%     ecog(1).ipad_ON_time=int32(ipad_ON_time);
%     ecog(1).ipad_OFF_time= int32(ipad_OFF_time);
    
    DELETE =  input('delete');
    GOOD = setdiff([1:length(ipad_ON_time)],DELETE)
    ecog(1).ipad_ON_time=int32(ipad_ON_time(GOOD));
    ecog(1).ipad_OFF_time= int32(ipad_OFF_time(GOOD));
    
    %assignin('base','ok_ON',ok_ON)
    %assignin('base','ipad_ON_time',ipad_ON_time)
    

else
    ecog(1).rest_time=[];
    ecog(1).active_time=[];
    ecog(1).ipad_ON_time=[];
    ecog(1).ipad_OFF_time=[];
    
end

%% Parse json file for timestamps

[description, timestamp] = nr_import_from_ipad_02(json_filename);

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

cd(aomatref_filedir)
save(aomatconv_sum_filename,'auxchan_lbl','M1_ch1','n_trials','aux','emg','ecog','Fs','thresh','ipad');