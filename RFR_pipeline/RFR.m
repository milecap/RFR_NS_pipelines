% This function runs the regression, filtering and rectification steps to
% the functional data series

% Author: Milena Capiglioni, University of Bern
% Contact: milena.capiglioni@extern.insel.ch
% Last update: Nov.2023

% Feel free to use the software and report any bugs!
% If you use/modify this code for your future publication, please cite the
% corresponding article:  "Stimulus-Induced Rotary Saturation imaging of
% visually evoked neuroelectric response: preliminary results and data
% analysis" (currently under review)

function RFR(func_path,seg_path,out_dir,cutoff,name)

func_data = load_nii(func_path);
func_data = double(func_data.img);
seg_data = load_nii(seg_path);
seg_data = seg_data.img;

% Check size match
if ~isequal(size(func_data,1:3), size(seg_data))
    error('Mask size and functional size differ, check data')
end

[xdim,ydim,zdim,tdim] = size(func_data);
TR = 2*zdim*139.55/1000;

% initial ft amplitude
ft_amp_0 = mean(func_data,4);

if strcmp(name,'')
    name = 'contrast';
end

niftiwrite(ft_amp_0,[out_dir name '_0.nii'])

%% Create mask of mean roi
seg_mean = zeros(size(func_data));
seg_aux = zeros(size(seg_data));
seg_list = unique(seg_data);

for s = 1:length(seg_list)
    for t = 1:tdim
        func_t = func_data(:,:,:,t);
        seg_mean(:,:,:,t) = seg_mean(:,:,:,t) + mean(func_t(seg_data == seg_list(s)))*(seg_data == seg_list(s));
    end
end

%% Apply regression
stats = struct([]);
resid = zeros(size(func_data));
regressors = zeros(size(seg_mean,4),3);

for x = 1:xdim
    for y = 1:ydim
        for z = 1:zdim
            % regression
            regressor(:,1) = squeeze(seg_mean(x,y,z,:));
            regressor(:,2) = linspace(0,1,size(seg_mean,4))';
            [~,~,stats] = glmfit(regressor,squeeze(func_data(x,y,z,:)));
            resid(x,y,z,:) = stats.resid;
        end
    end
end

% ft amplitude
ft_amp_R = mean(resid,4);
disp('Finished regression');

% save output nii
niftiwrite(ft_amp_R,[out_dir name '_R.nii'])

%% Apply high pass filter

m_hp = zeros(size(func_data));

parfor x = 1:xdim
    for y = 1:ydim
        for z = 1:zdim
            m_hp(x,y,z,:)= highpass(squeeze(resid(x,y,z,:)),cutoff,1/TR);
        end
    end
end

% ft amplitude
ft_amp_RF = mean(m_hp,4);
disp('Finished filtering');

% save output nii
niftiwrite(ft_amp_RF,[out_dir name '_RF.nii'])

%% Rectify signal

m_hp_rect= zeros(size(func_data));

parfor x = 1:xdim
    for y = 1:ydim
        for z = 1:zdim
            m_hp_rect(x,y,z,:) = m_hp(x,y,z,:) - mean(squeeze(m_hp(x,y,z,:)));
            m_hp_rect(x,y,z,:) = abs(m_hp_rect(x,y,z,:));
        end
    end
end

% ft amplitude
ft_amp_RFR = mean(m_hp_rect,4);
disp('Finished Rectifying');

% save output nii
niftiwrite(ft_amp_RFR,[out_dir name '_RFR.nii'])

%%
cd(out_dir)
save(['RFR_output' name])

%% Average in every roi

seg_power_std = zeros(length(seg_list),tdim);
seg_power_mean = zeros(length(seg_list),tdim);
seg_power_map = zeros(size(seg_data));

for s = 1:length(seg_list)
    power_f = squeeze(ft_amp_RFR(:,:,:));

    seg_power_mean(s) = mean(power_f(seg_data == seg_list(s)));
    seg_power_std(s) = std(power_f(seg_data == seg_list(s)));
    seg_power_map(:,:,:) = seg_power_map(:,:,:) + mean(power_f(seg_data == seg_list(s)))*(seg_data == seg_list(s));
end

% save output nii
niftiwrite(seg_power_map,[out_dir name '_seg_RFR.nii'])

end