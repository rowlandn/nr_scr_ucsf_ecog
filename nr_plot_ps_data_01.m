function nr_plot_ps_data_01(sum_file)

%% TO DO List
% be able to add 2 M1 channels if they exist

%% Load file
load(sum_file)

%% Scale traces 

if length(emg.chan) == 4

    % EMG Chan 1
    emgchan_1 = emg.chan(1).raw;
    std_emgchan_1 = std(emgchan_1)*5;
    emgchan_1_plus_std = emgchan_1 + std_emgchan_1;
    emgchan_1_plus_std_norm_1 = emgchan_1_plus_std/(std_emgchan_1*2);
    signal1_sc = emgchan_1_plus_std_norm_1;

    % EMG Chan 2
    emgchan_2 = emg.chan(2).raw;
    std_emgchan_2 = std(emgchan_2)*5;
    emgchan_2_plus_std = emgchan_2 + std_emgchan_2;
    emgchan_2_plus_std_norm_1 = emgchan_2_plus_std/(std_emgchan_2*2);
    signal2_sc = emgchan_2_plus_std_norm_1;
    
    % EMG Chan 3
    emgchan_3 = emg.chan(3).raw;
    std_emgchan_3 = std(emgchan_3)*5;
    emgchan_3_plus_std = emgchan_3 + std_emgchan_3;
    emgchan_3_plus_std_norm_1 = emgchan_3_plus_std/(std_emgchan_3*2);
    signal3_sc = emgchan_3_plus_std_norm_1;

    % EMG Chan 4
    emgchan_4 = emg.chan(4).raw;
    std_emgchan_4 = std(emgchan_4)*5;
    emgchan_4_plus_std = emgchan_4 + std_emgchan_4;
    emgchan_4_plus_std_norm_1 = emgchan_4_plus_std/(std_emgchan_4*2);
    signal4_sc = emgchan_4_plus_std_norm_1;
    
    % Aux Chan 1
    auxchan_1 = aux.chan(1).raw;
    auxchan_1_begin = auxchan_1(100);
    auxchan_1_end = auxchan_1(end-100);
    auxchan_1(1:100) = linspace(auxchan_1_begin,auxchan_1_begin,100);
    auxchan_1(end-100:end) = linspace(auxchan_1_end,auxchan_1_end,101);
    auxchan_1_min = min(auxchan_1);
    auxchan_1_zero = auxchan_1-auxchan_1_min;
    auxchan_1_max_zero = max(auxchan_1_zero);
    auxchan_1_one_zero = auxchan_1_zero/auxchan_1_max_zero;
    signal5_sc = auxchan_1_one_zero;

    % Aux Chan 2
    auxchan_2 = aux.chan(2).raw;
    auxchan_2_begin = auxchan_2(100);
    auxchan_2_end = auxchan_2(end-100);
    auxchan_2(1:100) = linspace(auxchan_2_begin,auxchan_2_begin,100);
    auxchan_2(end-100:end) = linspace(auxchan_2_end,auxchan_2_end,101);
    auxchan_2_min = min(auxchan_2);
    auxchan_2_zero = auxchan_2-auxchan_2_min;
    auxchan_2_max_zero = max(auxchan_2_zero);
    auxchan_2_one_zero = auxchan_2_zero/auxchan_2_max_zero;
    signal6_sc = auxchan_2_one_zero;

    % Aux Chan 3
    auxchan_3 = aux.chan(3).raw;
    auxchan_3_begin = auxchan_3(100);
    auxchan_3_end = auxchan_3(end-100);
    auxchan_3(1:100) = linspace(auxchan_3_begin,auxchan_3_begin,100);
    auxchan_3(end-100:end) = linspace(auxchan_3_end,auxchan_3_end,101);
    auxchan_3_min = min(auxchan_3);
    auxchan_3_zero = auxchan_3-auxchan_3_min;
    auxchan_3_max_zero = max(auxchan_3_zero);
    auxchan_3_one_zero = auxchan_3_zero/auxchan_3_max_zero;
    signal7_sc = auxchan_3_one_zero;
    
elseif length(emg.chan) == 2
    
    % EMG Chan 1
    emgchan_1 = emg.chan(1).raw;
    std_emgchan_1 = std(emgchan_1)*5;
    emgchan_1_plus_std = emgchan_1 + std_emgchan_1;
    emgchan_1_plus_std_norm_1 = emgchan_1_plus_std/(std_emgchan_1*2);
    signal1_sc = emgchan_1_plus_std_norm_1;

    % EMG Chan 2
    emgchan_2 = emg.chan(2).raw;
    std_emgchan_2 = std(emgchan_2)*5;
    emgchan_2_plus_std = emgchan_2 + std_emgchan_2;
    emgchan_2_plus_std_norm_1 = emgchan_2_plus_std/(std_emgchan_2*2);
    signal2_sc = emgchan_2_plus_std_norm_1;
    
    % Aux Chan 1
    auxchan_1 = aux.chan(1).raw;
    auxchan_1_begin = auxchan_1(100);
    auxchan_1_end = auxchan_1(end-100);
    auxchan_1(1:100) = linspace(auxchan_1_begin,auxchan_1_begin,100);
    auxchan_1(end-100:end) = linspace(auxchan_1_end,auxchan_1_end,101);
    auxchan_1_min = min(auxchan_1);
    auxchan_1_zero = auxchan_1-auxchan_1_min;
    auxchan_1_max_zero = max(auxchan_1_zero);
    auxchan_1_one_zero = auxchan_1_zero/auxchan_1_max_zero;
    signal3_sc = auxchan_1_one_zero;

    % Aux Chan 2
    auxchan_2 = aux.chan(2).raw;
    auxchan_2_begin = auxchan_2(100);
    auxchan_2_end = auxchan_2(end-100);
    auxchan_2(1:100) = linspace(auxchan_2_begin,auxchan_2_begin,100);
    auxchan_2(end-100:end) = linspace(auxchan_2_end,auxchan_2_end,101);
    auxchan_2_min = min(auxchan_2);
    auxchan_2_zero = auxchan_2-auxchan_2_min;
    auxchan_2_max_zero = max(auxchan_2_zero);
    auxchan_2_one_zero = auxchan_2_zero/auxchan_2_max_zero;
    signal4_sc = auxchan_2_one_zero;

    % Aux Chan 3
    auxchan_3 = aux.chan(3).raw;
    auxchan_3_begin = auxchan_3(100);
    auxchan_3_end = auxchan_3(end-100);
    auxchan_3(1:100) = linspace(auxchan_3_begin,auxchan_3_begin,100);
    auxchan_3(end-100:end) = linspace(auxchan_3_end,auxchan_3_end,101);
    auxchan_3_min = min(auxchan_3);
    auxchan_3_zero = auxchan_3-auxchan_3_min;
    auxchan_3_max_zero = max(auxchan_3_zero);
    auxchan_3_one_zero = auxchan_3_zero/auxchan_3_max_zero;
    signal5_sc = auxchan_3_one_zero;
    
elseif length(emg.chan) == 1
    
    % EMG Chan 1
    emgchan_1 = emg.chan(1).raw;
    std_emgchan_1 = std(emgchan_1)*5;
    emgchan_1_plus_std = emgchan_1 + std_emgchan_1;
    emgchan_1_plus_std_norm_1 = emgchan_1_plus_std/(std_emgchan_1*2);
    signal1_sc = emgchan_1_plus_std_norm_1;
    
    % Aux Chan 1
    auxchan_1 = aux.chan(1).raw;
    auxchan_1_begin = auxchan_1(100);
    auxchan_1_end = auxchan_1(end-100);
    auxchan_1(1:100) = linspace(auxchan_1_begin,auxchan_1_begin,100);
    auxchan_1(end-100:end) = linspace(auxchan_1_end,auxchan_1_end,101);
    auxchan_1_min = min(auxchan_1);
    auxchan_1_zero = auxchan_1-auxchan_1_min;
    auxchan_1_max_zero = max(auxchan_1_zero);
    auxchan_1_one_zero = auxchan_1_zero/auxchan_1_max_zero;
    signal2_sc = auxchan_1_one_zero;

    % Aux Chan 2
    auxchan_2 = aux.chan(2).raw;
    auxchan_2_begin = auxchan_2(100);
    auxchan_2_end = auxchan_2(end-100);
    auxchan_2(1:100) = linspace(auxchan_2_begin,auxchan_2_begin,100);
    auxchan_2(end-100:end) = linspace(auxchan_2_end,auxchan_2_end,101);
    auxchan_2_min = min(auxchan_2);
    auxchan_2_zero = auxchan_2-auxchan_2_min;
    auxchan_2_max_zero = max(auxchan_2_zero);
    auxchan_2_one_zero = auxchan_2_zero/auxchan_2_max_zero;
    signal3_sc = auxchan_2_one_zero;

    % Aux Chan 3
    auxchan_3 = aux.chan(3).raw;
    auxchan_3_begin = auxchan_3(100);
    auxchan_3_end = auxchan_3(end-100);
    auxchan_3(1:100) = linspace(auxchan_3_begin,auxchan_3_begin,100);
    auxchan_3(end-100:end) = linspace(auxchan_3_end,auxchan_3_end,101);
    auxchan_3_min = min(auxchan_3);
    auxchan_3_zero = auxchan_3-auxchan_3_min;
    auxchan_3_max_zero = max(auxchan_3_zero);
    auxchan_3_one_zero = auxchan_3_zero/auxchan_3_max_zero;
    signal4_sc = auxchan_3_one_zero;
    
elseif length(emg.chan) == 0
    
    % Aux Chan 1
    auxchan_1 = aux.chan(1).raw;
    auxchan_1_begin = auxchan_1(100);
    auxchan_1_end = auxchan_1(end-100);
    auxchan_1(1:100) = linspace(auxchan_1_begin,auxchan_1_begin,100);
    auxchan_1(end-100:end) = linspace(auxchan_1_end,auxchan_1_end,101);
    auxchan_1_min = min(auxchan_1);
    auxchan_1_zero = auxchan_1-auxchan_1_min;
    auxchan_1_max_zero = max(auxchan_1_zero);
    auxchan_1_one_zero = auxchan_1_zero/auxchan_1_max_zero;
    signal2_sc = auxchan_1_one_zero;

    if isempty(aux.chan(2).raw) % needed b/c ec_ep_0027 does not have aux chan 2
        signal3_sc = 0;
    else
        % Aux Chan 2
        auxchan_2 = aux.chan(2).raw;
        auxchan_2_begin = auxchan_2(100);
        auxchan_2_end = auxchan_2(end-100);
        auxchan_2(1:100) = linspace(auxchan_2_begin,auxchan_2_begin,100);
        auxchan_2(end-100:end) = linspace(auxchan_2_end,auxchan_2_end,101);
        auxchan_2_min = min(auxchan_2);
        auxchan_2_zero = auxchan_2-auxchan_2_min;
        auxchan_2_max_zero = max(auxchan_2_zero);
        auxchan_2_one_zero = auxchan_2_zero/auxchan_2_max_zero;
        signal3_sc = auxchan_2_one_zero;
    end

    % Aux Chan 3
    auxchan_3 = aux.chan(3).raw;
    auxchan_3_begin = auxchan_3(100);
    auxchan_3_end = auxchan_3(end-100);
    auxchan_3(1:100) = linspace(auxchan_3_begin,auxchan_3_begin,100);
    auxchan_3(end-100:end) = linspace(auxchan_3_end,auxchan_3_end,101);
    auxchan_3_min = min(auxchan_3);
    auxchan_3_zero = auxchan_3-auxchan_3_min;
    auxchan_3_max_zero = max(auxchan_3_zero);
    auxchan_3_one_zero = auxchan_3_zero/auxchan_3_max_zero;
    signal4_sc = auxchan_3_one_zero;
    
end

% ECOG M1_ch1
ecogchan_M1_ch1 = ecog.contact_pair(M1_ch1).remontaged_ecog_signal;
std_ecogchan_M1_ch1 = std(ecogchan_M1_ch1)*5;
ecogchan_M1_ch1_plus_std = ecogchan_M1_ch1 + std_ecogchan_M1_ch1;
ecogchan_M1_ch1_plus_std_norm_1 = ecogchan_M1_ch1_plus_std/(std_ecogchan_M1_ch1*2);

% ECOG M1_ch2
if exist('M1_ch2')
    ecogchan_M1_ch2 = ecog.contact_pair(M1_ch2).remontaged_ecog_signal;
    std_ecogchan_M1_ch2 = std(ecogchan_M1_ch2)*5;
    ecogchan_M1_ch2_plus_std = ecogchan_M1_ch2 + std_ecogchan_M1_ch2;
    ecogchan_M1_ch2_plus_std_norm_1 = ecogchan_M1_ch2_plus_std/(std_ecogchan_M1_ch2*2);
end
    

%% Plot channels

if length(emg.chan) == 4
    
    figure
    plot(auxchan_1_one_zero); hold on
    plot(emgchan_1_plus_std_norm_1,'y')
    plot(emgchan_2_plus_std_norm_1,'k')
    plot(emgchan_3_plus_std_norm_1,'c')
    plot(emgchan_4_plus_std_norm_1,'m')
    plot(auxchan_2_one_zero+1); hold on
    plot(auxchan_3_one_zero+2,'r')
    plot(ecogchan_M1_ch1_plus_std_norm_1+3,'k')  
    xlabel('ms')
    
    
elseif length(emg.chan) == 2
    
    figure
    %plot(auxchan_1_one_zero); hold on
    plot(emgchan_1_plus_std_norm_1,'k'); hold on
    plot(emgchan_2_plus_std_norm_1,'color',[0.5 0.5 0.5])
    plot(auxchan_2_one_zero+1,'color',[0.5 0.5 0.5])
    plot(auxchan_3_one_zero+2,'k','LineWidth',2)
    plot(ecogchan_M1_ch1_plus_std_norm_1+3,'k')
    if exist('M1_ch2')
        plot(ecogchan_M1_ch2_plus_std_norm_1+4,'color',[0.5 0.5 0.5])
    end
    xlabel('ms')
        
elseif length(emg.chan) == 1
    figure
    %plot(auxchan_1_one_zero); hold on
    plot(emgchan_1_plus_std_norm_1,'k'); hold on
    plot(auxchan_2_one_zero+1,'color',[0.5 0.5 0.5])
    plot(auxchan_3_one_zero+2,'k','LineWidth',2)
    plot(ecogchan_M1_ch1_plus_std_norm_1+3,'k')
    if exist('M1_ch2')
        plot(ecogchan_M1_ch2_plus_std_norm_1+4,'color',[0.5 0.5 0.5])
    end
    xlabel('ms')
    
elseif length(emg.chan) ==0
    figure
    plot(auxchan_1_one_zero); hold on
    if isempty(aux.chan(2).raw) % needed b/c ec_ep_0027 has no aux chan2
    else
        plot(auxchan_2_one_zero+1,'m'); hold on
    end
    plot(auxchan_3_one_zero+2,'r')
    plot(ecogchan_M1_ch1_plus_std_norm_1+3,'k')
    if exist('M1_ch2')
        plot(ecogchan_M1_ch2_plus_std_norm_1+4,'k')
    end
    xlabel('ms')
end



% plot TMPREJ epochs in red
for i = 1:size(TMPREJ,1)
    TMPREJ_vec = int32(linspace(TMPREJ(i,1),TMPREJ(i,2),TMPREJ(i,2)-TMPREJ(i,1)));
    plot(TMPREJ_vec,ecogchan_M1_ch1_plus_std_norm_1(int32(TMPREJ_vec))+3,'r')
    if exist('M1_ch2')
        plot(TMPREJ_vec,ecogchan_M1_ch2_plus_std_norm_1(int32(TMPREJ_vec))+4,'r')
    end
end

% Get y-limits
ylim = get(gca,'YLim');


for i = 1:length(ipad.events_sc.sound_ON_sc)
    plot([ipad.events_sc.sound_ON_sc(i) ipad.events_sc.sound_ON_sc(i)],ylim,'r','LineWidth',2) % here I am treating the beep as the start of rest
end

% Plot original ipad beep detection at resampled rate as an internal check
length_ecog_ipad_ON_time = length(ecog.ipad_ON_time);
mat_length_ecog_ipad_ON_time = linspace(0,length_ecog_ipad_ON_time-1,length_ecog_ipad_ON_time);
for i = 1:length_ecog_ipad_ON_time
    text(double(ecog.ipad_ON_time(i)),2,num2str(mat_length_ecog_ipad_ON_time(i)))
end


for i = 1:length(ipad.events_sc.prep_ON_sc)
    plot([ipad.events_sc.prep_ON_sc(i) ipad.events_sc.prep_ON_sc(i)],ylim,'b','LineWidth',2)
end

% prep_OFF is the same as movement 'on'
for i = 1:length(ipad.events_sc.prep_OFF_sc)
    plot([ipad.events_sc.prep_OFF_sc(i) ipad.events_sc.prep_OFF_sc(i)],ylim,'g','LineWidth',2)
end

set(gca,'YLim',[-2 5.5])

% plot results of Determine Events
for i = 1:length(ecog.active_time)
plot([ecog.active_time(i):50:ecog.rest_time(i)],2.5,'color',[0.5 0.5 0.5],'Marker','.','MarkerSize',8)
end

% plot results from 2-second epochs
for i = 1:length(rest.contact_pair(M1_ch1).epoch)
    plot(rest.contact_pair(M1_ch1).epoch(i).time(1):50:rest.contact_pair(M1_ch1).epoch(i).time(2),2.8,'r.')
end

for i = 1:length(prep.contact_pair(M1_ch1).epoch)
    plot(prep.contact_pair(M1_ch1).epoch(i).time(1):50:prep.contact_pair(M1_ch1).epoch(i).time(2),2.7,'b.')
end

for i = 1:length(active.contact_pair(M1_ch1).epoch)
    plot(active.contact_pair(M1_ch1).epoch(i).time(1):50:active.contact_pair(M1_ch1).epoch(i).time(2),2.6,'g.')
end



wichAxis = 'x';
minscroll = 1;
maxscroll = size(ecog_traces,2);
initialValue = minscroll;
width = size(ecog_traces,2)/2;
% x = min : max;
% 
% bar(x,x);

pos = get(gca,'position');
new_pos=[pos(1) pos(2)-0.1 pos(3) 0.05];

charLabel = @(axisHandle, farglist)( ...
    set(axisHandle, 'XTickLabel',( ...
    char(get(axisHandle,'XTick') + 96) ...
    )') ...
    );

panScrollZoom(gca, wichAxis, minscroll, maxscroll, initialValue, width, new_pos, charLabel, {})

%function scrollplotdemo
%%%%% Generate and plot data
% x=0:1e-2:2*pi;
% y=sin(x);
% dx=2;
% % if exist('M1_ch2')
% %         plot(ecogchan_M1_ch2_plus_std_norm_1+4,'color',[0.5 0.5 0.5])
% %     end
% 
% 
% % % dx is the width of the axis 'window'
% a=gca;
% %p=plot(x,y);
% %%%%% Set appropriate axis limits and settings
% set(gcf,'doublebuffer','on');
% % This avoids flickering when updating the axis
% set(a,'xlim',[0 dx]);
% %set(a,'ylim',[min(y) max(y)]);
% %%%%% Generate constants for use in uicontrol initialization
% pos=get(a,'position');
% Newpos=[pos(1) pos(2)-0.1 pos(3) 0.05];
% % This will create a slider which is just underneath the axis
% % but still leaves room for the axis labels above the slider
% xmax = length(ecogchan_M1_ch2_plus_std_norm_1);
% 
% %xmax=max();
% S=['set(gca,''xlim'',get(gcbo,''value'')+[0 ' num2str(dx) '])'];
% % Setting up callback string to modify XLim of axis (gca)
% % based on the position of the slider (gcbo)
% %%%%% Creating Uicontrol
% h=uicontrol('style','slider',...
%       'units','normalized','position',Newpos,...
%       'callback',S,'min',0,'max',xmax-dx);


