function nr_spect_plot(spec_data,chan,phase,fig_han,sp1,sp2,sp3)


%% parse subject and load all data

if ~isempty(strfind(spec_data,'ps_'))
    spec_data_ps = strfind(spec_data,'ps_');
    sbj = spec_data(spec_data_ps(1):spec_data_ps(1)+9);
elseif ~isempty(strfind(spec_data,'ec_'))
    spec_data_ec = strfind(spec_data,'ec_');
    sbj = spec_data(spec_data_ec(1):spec_data_ec(1)+9);
end

load(spec_data)
load_sum_eval = ['load sum_',sbj,'_ipad_postlead'];
eval(load_sum_eval)

%% define variables

if chan == 'M1_ch1'
    chan = M1_ch1;
elseif chan == 'M1_ch2'
    if exist('M1_ch2')
        chan = M1_ch2;
    else
        chan = M1_ch1;
    end
elseif chan == 'S1_ch1'
    if exist('M1_ch2')
        chan = M1_ch1-2;
    else
        chan = M1_ch1-1;
    end
elseif chan == 'S1_ch2'
    if exist('M1_ch2')
        chan = M1_ch2-2;
    else
        chan = M1_ch1-1;
    end
elseif chan == 'P1_ch1'
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
elseif chan == 'P1_ch2'
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


% if chan == 'M1_ch1'
%     chan = M1_ch1;
% elseif chan == 'M1_ch2'
%     if exist('M1_ch2')
%         chan = M1_ch2;
%     else
%         chan = M1_ch1;
%     end
% elseif chan == 'S1_ch1'
%     if exist('M1_ch2')
%         chan = M1_ch1-2;
%     else
%         chan = M1_ch1-1;
%     end
% elseif chan == 'S1_ch2'
%     if exist('M1_ch2')
%         chan = M1_ch2-2;
%     else
%         chan = M1_ch1-1;
%     end
% end

% if isempty(get(0,'Children'))
%     fig_han = figure;
% else
%     %fig_han = figure;
%     figure(fig_han)
% end

figure(fig_han)

%% create plot

if phase == 'rp'
    A2plot = A2plot_rp;
elseif phase == 'mn'
    A2plot = A2plot_on;
elseif phase == 'mf'
    A2plot = A2plot_off;
end




if n_ecog == 5 & phase == 'rp'
    val1 = min(min(min(A2plot(1:100,:,:))));
    val2 = max(max(max(A2plot(1:100,:,:))));
    clims1 = [val1 val2];
elseif n_ecog == 28 & phase == 'rp'
    ff = find(faxis==125);
    val1 = min(min(min(A2plot(1:ff,:,:))));
    val2 = max(max(max(A2plot(1:ff,:,:))));
    clims1 = [val1 val2];
elseif n_ecog == 5 & phase == 'mn' | n_ecog == 5 & phase == 'mf'
    val1 = min(min(min(A2plot(1:100,:,:))));
    val2 = max(max(max(A2plot(1:100,:,:))));
    clims1 = [val1 val2];
elseif n_ecog == 28 & phase == 'mn' | n_ecog == 28 & phase == 'mf'
    ff = find(faxis==125);
    val1 = min(min(min(A2plot(1:ff,:,:))));
    val2 = max(max(max(A2plot(1:ff,:,:))));
    clims1 = [0.8 1.3];
end
    

%assignin('base','A2plot2',A2plot)

%clims1
subplot(sp1,sp2,sp3);
hold(gca,'on');
tmp1 = A2plot(1:100,:,chan); 
faxis_new = faxis(1:100);
imagesc(taxis,faxis_new,tmp1,clims1);
plot([0 0],ylim,'k:');
hold(gca,'off');
set(gca,'YDir','normal');
set(gca,'Xlim',[0-PRE POST]);
set(gca,'Ylim',[0 200]);
% xlabel('time (sec)');
% ylabel('frequency (Hz)');
title([sbj(1:2),'\',sbj(3:5),'\',sbj(6:end)])
   
    
   


