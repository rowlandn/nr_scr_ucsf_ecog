function [auxchan_lbl] = nr_htk_auxchan_lbl(cond, chan1, chan2, chan3, chan4, chan5, chan6, chan7, chan8)

% This function organizes auxiliary channel data labels into a structure that can
% be saved along with analog data files

%% TO DO
% for ps_pd_0100 dlt hooked up before ecr and the script reversed the order
% so go and fix this (the sort needs to be perhaps a separate step

for i = 1:8
    header_eval = ['auxchan_con{',num2str(i),',1} = chan',num2str(i),';'];
    eval(header_eval)
end


ipd_idx = find(ismember(auxchan_con,'ipd'));
auxchan_lblchan1_str = ['auxchan_lbl.',cond,'.chan1 = auxchan_con{ipd_idx};'];
eval(auxchan_lblchan1_str)

acl_idx = find(ismember(auxchan_con,'acl'));
auxchan_lblchan2_str = ['auxchan_lbl.',cond,'.chan2 = auxchan_con{acl_idx};'];
eval(auxchan_lblchan2_str)

trg_idx = find(ismember(auxchan_con,'trg'));
auxchan_lblchan3_str = ['auxchan_lbl.',cond,'.chan3 = auxchan_con{trg_idx};'];
eval(auxchan_lblchan3_str)

if ~isempty(strfind(auxchan_con,'ecr'))
    ecr_idx = find(ismember(auxchan_con,'ecr'));
end

if ~isempty(strfind(auxchan_con,'fcr'))
    fcr_idx = find(ismember(auxchan_con,'fcr'));
end

if ~isempty(strfind(auxchan_con,'bcp'))
    bcp_idx = find(ismember(auxchan_con,'bcp'));
end

if ~isempty(strfind(auxchan_con,'dlt'))
    dlt_idx = find(ismember(auxchan_con,'dlt'));
end

emg_idx{1,1} = ecr_idx;
emg_idx{2,1} = fcr_idx;
emg_idx{3,1} = bcp_idx;
emg_idx{4,1} = dlt_idx;

find_empty = cellfun('isempty',emg_idx);
find_non_empty = find(not(find_empty));
length_non_empty = length(find_non_empty);

for i = 1:length_non_empty
    emg_idx_sort_eval = ['emg_idx_sort(',num2str(i),',1) = emg_idx{find_non_empty(',num2str(i),')};'];
    eval(emg_idx_sort_eval)
    emg_idx_find_non_empty_eval = ['emg_idx_sort(',num2str(i),',2) = find_non_empty(',num2str(i),');'];
    eval(emg_idx_find_non_empty_eval)
end
emg_idx_sort_sort = sortrows(emg_idx_sort,1);

for i = 1:length_non_empty
    emg_eval = ['auxchan_lbl.',cond,'.emg',num2str(i),' = auxchan_con{cell2mat(emg_idx(emg_idx_sort_sort(',num2str(i),',2)))};'];
    eval(emg_eval)
end
    
% cd(dir)
% 
% save_eval = ['save(''',sbj,'_auxchan'',''auxchan_con'',''auxchan_lbl'');'];
% eval(save_eval)
% 
