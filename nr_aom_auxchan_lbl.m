function [auxchan_lbl] = nr_aom_auxchan_lbl(cond,emg1,emg2,emg3,emg4,chan1,chan2,chan3)

% This function organizes auxiliary channel data labels into a structure that can
% be saved along with analog data files



for i = 1:7
    if i < 5
        header_eval = ['auxchan_lbl.',cond,'.emg',num2str(i),' = emg',num2str(i),';'];
        eval(header_eval)
    else
        header_eval = ['auxchan_lbl.',cond,'.chan',num2str(i-4),' = chan',num2str(i-4),';'];
        eval(header_eval)
    end
end

% cd(dir)
% 
% save_eval = ['save(''',sbj,'_auxchan_lbl'',''auxchan_lbl'');'];
% eval(save_eval)
