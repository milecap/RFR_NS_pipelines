% Author: Milena Capiglioni, University of Bern
% Contact: milena.capiglioni@extern.insel.ch
% Last update: Nov.2023

% If you use/modify this code for your future publication, please cite the
% corresponding article:  "Stimulus-Induced Rotary Saturation imaging of
%visually evoked neuroelectric response: preliminary results and data
%analysis" (currently under review)


function RFR_visualize(Subjects_folder,RFR_folder)
% RFR visualization

show_im = 1; % show images
save_im = 0; % save images

if nargin == 0
    % No arguments provided, interactively select folder with subjects and
    % RFR output folder
    Subjects_folder = uigetdir('Select folder with subjects');
    RFR_folder = input('Write name of RFR processing output on the form RFR_cof_n_dd_mm_yy: ');
    ROI = questdlg('Which ROI do you want to analyse', 'ROI','V1','G_subcallosal','S_circular_insula_ant','V1');
else
    % Check the number of arguments provided
    if nargin < 2
        error('Not enough input arguments. Provide both Subjects_folder and RFR_folder.');
    end
end

% Overlay over high or low res anatomical image
hr_anat = 0;

%% Define patient folders

cd(Subjects_folder);
files = dir(Subjects_folder);
dirFlag = contains({files.name},{'Sub_'})&[files.isdir];
subjects = files(dirFlag);

output_dir = [Subjects_folder RFR_folder filesep];

if ~isfolder(output_dir)
    mkdir(output_dir);
else
    if save_im == 1
        prompt = 'output folder already exists! if you continue, you might rewrite content, do you want to continue? 0=NO 1=YES:  ';
        answer = input(prompt);
        if answer == 0
            error('Output folder already exists');
        end
    end
end


%% ------------ Show the output of the RFR postprocessing
for sub = 1:length(subjects)

    work_dir = strcat(Subjects_folder,subjects(sub).name,filesep);
    current_subject = subjects(sub).name;

    folders = dir(work_dir);
    dirFlag = contains({folders.name},{RFR_folder})&[folders.isdir];
    output_folder = folders(dirFlag);
    if isempty(output_folder)
        disp(['Subject', current_subject, ' does not have any sequence of interest']);
        continue;
    end

    current_func = output_folder(1).name;
    data_dir = [work_dir 'data' filesep];
    results_dir = [work_dir current_func filesep];


    % Read .nii for each analized case
    seg_path_noStim = [data_dir 'noStim/SLoff/anat/rT1w_norm_seg.nii'];
    seg_image = load_nii(seg_path_noStim);
    seg_image = seg_image.img;
    if hr_anat
        anatoff_image = load_nii([data_dir 'noStim/SLon/anat/rT1w_norm.nii']);
        anatoff_image = double(anatoff_image.img);
        anaton_image = load_nii([data_dir 'noStim/SLon/anat/rT1w_norm.nii']);
        anaton_image = double(anaton_image.img);
    else
        anatoff_image = load_nii([data_dir 'noStim/SLoff/func/data.nii']);
        anatoff_image = squeeze(anatoff_image.img(:,:,:,1));
        anaton_image = load_nii([data_dir 'noStim/SLon/func/data.nii']);
        anaton_image = squeeze(anaton_image.img(:,:,:,1));
    end

    for slice = 1:size(seg_image,3)
        mask(:,:,slice) = logical(imrotate(seg_image(:,:,slice),90));
        anatoff(:,:,slice) = imrotate(anatoff_image(:,:,slice),90);
        anaton(:,:,slice) = imrotate(anaton_image(:,:,slice),90);
    end
    mask_mos = slices2mosaic(mask); 
    anatoff_mos = slices2mosaic(anatoff);
    anaton_mos = slices2mosaic(anaton);

    % Read result actuvation maps
    txt ={'_0.nii','_RFR.nii'};


    c = 2; % after RFR, c=1 before RFR

    act_maps = dir([results_dir]);
    dirFlag = logical(contains({act_maps.name},txt{c}) .* ~contains({act_maps.name},'seg'));
    act_maps = act_maps(dirFlag);

    %%%%%%%%%%%%% Figure 3 paper MRM %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% SL on and SL off maps for each volunteer %%%%%%%%%%%%%
    fig = figure(2);
    tlo = tiledlayout(1,length(act_maps)/2,'TileSpacing','compact','Padding','compact');

    %%%%%%%%%%% Create SL off maps %%%%%%%%%%%
    func_str = load_nii([act_maps(1).folder filesep act_maps(1).name]);
    map_1 = func_str.img; %VisStim

    func_str = load_nii([act_maps(1).folder filesep act_maps(2).name]);
    map_2 = func_str.img; %NoStim

    map_sloff = map_1 - map_2; %SL off -> VisStim - NoStim

    % rotate image and make mosaic
    for slice = 1:size(map_sloff,3)
        map_sloff(:,:,slice) = imrotate(map_sloff(:,:,slice),90);
    end
    act_sloff = slices2mosaic(map_sloff);

    mask_mos(mask_mos == 0) = NaN;

    %%%%%%%%%%% Create SL on maps %%%%%%%%%%%
    func_str = load_nii([act_maps(1).folder filesep act_maps(3).name]);
    map_3 = func_str.img; %VisStim

    func_str = load_nii([act_maps(1).folder filesep act_maps(4).name]);
    map_4 = func_str.img; %NoStim

    map_slon = map_3 - map_4; %SL on -> VisStim - NoStim

    % rotate image and make mosaic
    for slice = 1:size(map_slon,3)
        map_slon(:,:,slice) = imrotate(map_slon(:,:,slice),90);
    end
    act_slon = slices2mosaic(map_slon);


    %%%%%%%%%%% Plot maps and normalize scale %%%%%%%%%%%
    ax = nexttile(2); % SL on
    [min_s,max_s,norm,norm_min,norm_max] = im_overlay_RFR(anaton_mos,act_slon,mask_mos,-1,1,[],[],[],ax);
    title([act_maps(3).name '-' act_maps(4).name] ,'interpreter','none')

    ax = nexttile(1);% SL off
    [min_s,max_s,norm,norm_min,norm_max] = im_overlay_RFR(anatoff_mos,act_sloff,mask_mos,min_s,max_s,norm,norm_min,norm_max,ax);
    title([act_maps(1).name '-' act_maps(2).name] ,'interpreter','none')
    
    set(figure(2),'Position',[684 374 771 360]);
    sgtitle([current_subject ' ' current_func],'interpreter','none')

    if save_im == 1
        saveas(fig,[output_dir 'RFR_' current_subject '.svg'])
    end

    waitforbuttonpress;
end



