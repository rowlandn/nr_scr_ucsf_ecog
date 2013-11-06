function nr_auxchan_dat(sbj, cond, chan1, chan2, chan3, chan4, chan5, chan6, chan7, chan8, dir)

% This function organizes auxiliary channel data labels into a structure that can
% be saved along with analog data files

header = [sbj, '_auxchan_dat'];

header.sbj = sbj
header.cond.chan1 = chan1;
header.cond.chan2 = chan2;
header.cond.chan3 = chan3;
header.cond.chan4 = chan4;
header.cond.chan5 = chan5;
header.cond.chan6 = chan6;
header.cond.chan7 = chan7;
header.cond.chan8 = chan8;




% cd(dir)
% 
% save(header)
