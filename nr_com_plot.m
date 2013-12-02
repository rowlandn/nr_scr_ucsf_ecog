function nr_com_plot(com_file,cond,chans,lbl,sp1,sp2,sp3,sp4,sp5,sp6,sp7,sp8,sp9,sp10,sp11,sp12)

load(com_file)
% 
% if ~isempty(strfind(com_file,'ps_'))
%     com_file_ps = strfind(com_file,'ps_');
%     sbj = com_file(com_file_ps(1):com_file_ps(1)+9);
% elseif ~isempty(strfind(com_file,'ec_'))
%     com_file_ec = strfind(com_file,'ec_');
%     sbj = com_file(com_file_ec(1):com_file_ec(1)+9);
% end


%% Graph comodulogram

Comodulogram_surr_eval = ['Comodulogram_surr = allt.',cond,'.',chans,'.Comodulogram_surr;'];
eval(Comodulogram_surr_eval)
PhaseFreqVector_eval = ['PhaseFreqVector = allt.',cond,'.',chans,'.PhaseFreqVector;'];
eval(PhaseFreqVector_eval)
PhaseFreq_BandWidth_eval = ['PhaseFreq_BandWidth = allt.',cond,'.',chans,'.PhaseFreq_BandWidth;'];
eval(PhaseFreq_BandWidth_eval)
AmpFreqVector_eval = ['AmpFreqVector = allt.',cond,'.',chans,'.AmpFreqVector;'];
eval(AmpFreqVector_eval)
AmpFreq_BandWidth_eval = ['AmpFreq_BandWidth = allt.',cond,'.',chans,'.AmpFreq_BandWidth;'];
eval(AmpFreq_BandWidth_eval)
MI_means_eval = ['MI_means = tres.',cond,'.',chans,'.MI_means;'];
eval(MI_means_eval)
pts_eval = ['pts = tres.',cond,'.',chans,'.pts;'];
eval(pts_eval)
P_eval = ['P = tres.',cond,'.',chans,'.P;'];
eval(P_eval)

Comodulogram_surr(1,1,:)=0.00001;
Clim2 = max(max(max(Comodulogram_surr)));
Clim1 = min(min(min(Comodulogram_surr)));

subplot(sp1,sp2,sp3)
C=squeeze(Comodulogram_surr(:,:,2));
contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,C',30,'lines','none')
set(gca,'fontsize',9)
title([chans(1:2),' ',chans(4:6)])
% ylabel('Amplitude Frequency (Hz)')
% xlabel('Phase Frequency (Hz)')
% colorbar
caxis([Clim1 Clim2])
if strcmp(lbl,'off')
    set(gca,'XTickLabel','')
end

mean_COM_B(1)=nanmean(nanmean(Comodulogram_surr(:,:,2)));

subplot(sp4,sp5,sp6)
C=squeeze(Comodulogram_surr(:,:,3));
contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,C',30,'lines','none')
set(gca,'fontsize',9)
% ylabel('Amplitude Frequency (Hz)')
% xlabel('Phase Frequency (Hz)')
% colorbar
caxis([Clim1 Clim2])
if strcmp(lbl,'off')
    set(gca,'XTickLabel','')
    set(gca,'YTickLabel','')
elseif strcmp(lbl,'on')
    set(gca,'YTickLabel','')    
end
mean_COM_B(2)=nanmean(nanmean(Comodulogram_surr(:,:,3)));

subplot(sp7,sp8,sp9)
C=squeeze(Comodulogram_surr(:,:,4));
contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,C',30,'lines','none')
set(gca,'fontsize',9)
% ylabel('Amplitude Frequency (Hz)')
% xlabel('Phase Frequency (Hz)')
% colorbar
caxis([Clim1 Clim2])
if strcmp(lbl,'off')
    set(gca,'XTickLabel','')
    set(gca,'YTickLabel','')
elseif strcmp(lbl,'on')
    set(gca,'YTickLabel','')    
end
mean_COM_B(3)=nanmean(nanmean(Comodulogram_surr(:,:,4)));

subplot(sp10,sp11,sp12)
plot([1:length(MI_means)],MI_means,'r*','MarkerSize',4)
x=find(pts(:,3)==2);
hold on
plot(x,MI_means(x),'b*','MarkerSize',4)
x=find(pts(:,3)==3);
plot(x,MI_means(x),'g*','MarkerSize',4)
set(gca,'XLim',[0 length(MI_means)])
title(P,'FontSize',9)
set(gca,'fontsize',9)
if strcmp(lbl,'off')
    set(gca,'XTickLabel','')
    set(gca,'YTickLabel','')
end

assignin('base','mean_COM_B',mean_COM_B)



