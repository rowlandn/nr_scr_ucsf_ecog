function [json_mat_conv_file] = nr_import_from_ipad_01(json_filename,json_mat_conv_dir)
%% Rich Stoner, 2011
% import json file from responseApp ipad app into matlab
% clear all
% clc

% set this filename
% json_filename = 'Rich.stoner-exp-2011-09-13-14-52-32.json';
% [json_filename, json_dir] =uigetfile('*.json','Ipad files (*.json)'); % choose abf file from CD or proper directory

% Parse subject ID and append 'jsn' to filename
if ~isempty(strfind(json_filename,'ps'))
    json_filename_ps = strfind(json_filename,'ps');
    sbj = json_filename(json_filename_ps(1):json_filename_ps(1)+9);
elseif ~isempty(strfind(json_filename,'ec'))
    json_filename_ec = strfind(json_filename,'ec');
    sbj = json_filename(json_filename_ec(1):json_filename_ec(1)+9);
end

json_filename_slash = strfind(json_filename,'/');
find_dir = json_filename_slash(end);
json_mat_conv_file = [sbj,'_jsn_',json_filename(find_dir+12:end-5)]; 

% Read json file
json_str = fileread(json_filename);

% parsed using third party library
[dataset] = nr_parse_json(json_str);
array_list = dataset{1,1};

% distribute into individual arrays
description = {};   % cell array
timestamp = [];     % standard array of doubles

for i = 1:length(array_list)

%     uncomment to display 
%     disp(sprintf('%s\t%s ms', array_list{1,i}{1,1}, array_list{1,i}{1,2}));

    description{i} = array_list{1,i}{1,1};
    timestamp(i) = str2double(array_list{1,i}{1,2});
end

% cd to output directory and save
cd(json_mat_conv_dir)
save(json_mat_conv_file, 'description','timestamp')
