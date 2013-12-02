function nr_Determine_trials_ok(dat_file)

load(dat_file)

clear trials_ok*

%% Determine trials_ok_emg

for i = 1:length(active.contact_pair(M1_ch1).epoch)
    if active.contact_pair(M1_ch1).epoch(i).time(1) < ecog.go_time(i)
        trials_ok_emg(i) = 0;
    else
        trials_ok_emg(i) = i;
    end
end

%% Determine trials_ok_TMPREJ
    
for i = 1:size(TMPREJ,1)
    for j = 1:n_trials
        TMPREJ_array = int32(linspace(TMPREJ(i,1),TMPREJ(i,2),TMPREJ(i,2)-TMPREJ(i,1)));
        rest_array = int32(linspace(double(rest.contact_pair(M1_ch1).epoch(j).time(1)),double(rest.contact_pair(M1_ch1).epoch(j).time(2)),...
            double(rest.contact_pair(M1_ch1).epoch(j).time(2)-rest.contact_pair(M1_ch1).epoch(j).time(1))));
        prep_array = int32(linspace(double(prep.contact_pair(M1_ch1).epoch(j).time(1)),double(prep.contact_pair(M1_ch1).epoch(j).time(2)),...
            double(prep.contact_pair(M1_ch1).epoch(j).time(2)-prep.contact_pair(M1_ch1).epoch(j).time(1))));
        active_array = int32(linspace(double(active.contact_pair(M1_ch1).epoch(j).time(1)),double(active.contact_pair(M1_ch1).epoch(j).time(2)),...
            double(active.contact_pair(M1_ch1).epoch(j).time(2)-active.contact_pair(M1_ch1).epoch(j).time(1))));
        
        if ~isempty(find(ismember(rest_array,TMPREJ_array)))
            trials_ok_TMPREJ(1,j,i) = j;
        else
            trials_ok_TMPREJ(1,j,i) = 0;
        end
        
        if ~isempty(find(ismember(prep_array,TMPREJ_array)))
            trials_ok_TMPREJ(2,j,i) = j;
        else
            trials_ok_TMPREJ(2,j,i) = 0;
        end
        
        if ~isempty(find(ismember(active_array,TMPREJ_array)))
            trials_ok_TMPREJ(3,j,i) = j;
        else
            trials_ok_TMPREJ(3,j,i) = 0;
        end
    end
end
            
[trials_ok_TMPREJ_row,trials_ok_TMPREJ_col] = find(trials_ok_TMPREJ);
trials_ok_TMPREJ = unique(rem(trials_ok_TMPREJ_col,n_trials));

save(dat_file,'trials_ok_emg','trials_ok_TMPREJ')

