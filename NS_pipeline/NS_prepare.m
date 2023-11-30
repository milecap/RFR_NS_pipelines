%% NS prepare pipeline
%This main file will call the the neccessary routines to generate the NS
%results of the publication: "Stimulus-Induced Rotary Saturation imaging of
%visually evoked neuroelectric response: preliminary results and data
%analysis" (currently under review)

% Author: Milena Capiglioni, University of Bern
% Contact: milena.capiglioni@extern.insel.ch
% Last update: Nov.2023

% If you use/modify this code for your future publication, please cite the
% corresponding article



%%
clear global; clc; clear; close all;

% Add utils folder to the path
addpath(genpath([pwd filesep 'utils' filesep]));

%% Is it an only patient?
prompt = 'Do you want to analyse more than one subject? 0=NO 1=YES:  ';
mSub = input(prompt);

%------------ Select Patient / patients to analyze
if mSub == 0
    %select single patient
    work_dir_aux = strcat(uigetdir(path,'Select subject folder'),filesep);
    parts = strsplit(work_dir_aux,filesep);
    work_dir_aux = strcat(strjoin(parts(1:end-2),filesep),filesep);
    subjects.name = parts{end-1}; %Im sure there was a shorte way to do this
else
    %Select folder were all the patients are
    work_dir_aux = strcat(uigetdir(path,'Select folder with subjects'),filesep);
    cd(work_dir_aux);
    files = dir([work_dir_aux]);
    dirFlag = contains({files.name},{'Sub'})&[files.isdir];
    subjects = files(dirFlag);
end

cutoff = 0.1;

%% ------------ Run the processing for each subject------------------------------

for sub = 1:length(subjects)

    disp(['Starting analysis of subject ' , num2str(sub), ' of ', num2str(length(subjects)) ' : ',subjects(sub).name]);

    work_dir = strcat(work_dir_aux,subjects(sub).name,filesep);
    cd(work_dir);

    folders = dir([work_dir filesep 'data' filesep]);
    dirFlag = ismember({folders.name},{'noStim','VisStim'})&[folders.isdir];
    func_folders = folders(dirFlag);
    if isempty(func_folders)
        disp(['Subject', subjects(sub).name, ' does not have any sequence of interest']);
        continue;
    end

    %Define output directory, create if necessary
    out_dir_total = [work_dir 'NS_cof_' num2str(cutoff) filesep];
    if ~exist(out_dir_total,'dir')
        mkdir(out_dir_total);
    elseif ~isempty(dir([out_dir_total 'NS_output*.mat']))
        continue;
    end

    % Select data and run pipeline
    func_seq = {'noStim','VisStim'};

    for f = 1:2 %func_seq
        func_sloff = [work_dir 'data' filesep func_seq{f} filesep 'SLoff' filesep 'func' filesep 'rdata.nii'];
        func_slon = [work_dir 'data' filesep func_seq{f} filesep 'SLon' filesep 'func' filesep 'rdata.nii'];
        seg_data = [work_dir 'data' filesep func_seq{1} filesep 'SLoff' filesep 'anat' filesep 'rT1w_norm_seg.nii']; % segmentation only comes from noStim series

        NS(func_sloff,func_slon,seg_data,out_dir_total,cutoff,func_seq{f}) %Run pipeline
        clear func_slon func_sloff seg_data
    end
end

%% Perform statistics over the postprocessed data (Figures 5)

% This function requires user interaction!
NS_statistical(work_dir_aux,'NS_cof_0.1','V1');
