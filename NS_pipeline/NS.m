% This function runs the Normalize substraction to
% the functional data series (+ high pass filtering)

% Author: Milena Capiglioni, University of Bern
% Contact: milena.capiglioni@extern.insel.ch
% Last update: Nov.2023

% Feel free to use the software and report any bugs!
% If you use/modify this code for your future publication, please cite the
% corresponding article:  "Stimulus-Induced Rotary Saturation imaging of
% visually evoked neuroelectric response: preliminary results and data
% analysis" (currently under review)


function NS(func_off,func_on,seg_off,out_dir,cutoff,name)

func_off = load_nii(func_off);
func_off = double(func_off.img);

func_on = load_nii(func_on);
func_on = double(func_on.img);

seg_off = load_nii(seg_off);
seg_off = double(seg_off.img);

% Check size match
if ~isequal(size(func_off), size(func_on))
    error('Functional sizes differ, check data')
end

[xdim,ydim,zdim,tdim] = size(func_off);
TR = 2*zdim*139.55/1000; %This is hardcoded to the TR of the sequence, could be read from the dicoms but no dicom data can be shared.

% initial ft amplitude
% ft_amp_0 = std(func_off,4);

if strcmp(name,'')
    name = 'contrast';
end

%% Apply deconvolution

r_dec = zeros(size(func_off));
resid = zeros(size(func_off));


for x = 1:xdim
    for y = 1:ydim
        for z = 1:zdim

            if seg_off(x,y,z)~=0 && squeeze(func_off(x,y,z,1))~=0

                % Deconvolution rest
                [~,r_dec(x,y,z,:)]  = deconv(squeeze(func_on(x,y,z,:)),squeeze(func_off(x,y,z,:)));

                % Claus post-processing map
                resid(x,y,z,:) = abs(squeeze(r_dec(x,y,z,:)) - mean(squeeze(r_dec(x,y,z,:))));

            end
        end
    end
end

% ft amplitude
dec_map = mean(resid,4);
disp('Finished deconvolution');

niftiwrite(dec_map,[out_dir name '_NS.nii'])

%% Apply high pass filter

m_hp = zeros(size(func_off));
m_hp_rec = zeros(size(func_off));

parfor x = 1:xdim
    for y = 1:ydim
        for z = 1:zdim

            if seg_off(x,y,z)~=0 && squeeze(func_off(x,y,z,1))~=0

                m_hp(x,y,z,:) = highpass(squeeze(r_dec(x,y,z,:)),cutoff,1/TR);
                m_hp_rec(x,y,z,:) = abs(squeeze(m_hp(x,y,z,:)) - mean(squeeze(m_hp(x,y,z,:))));
            end
        end
    end
end

% ft amplitude
dec_map_hp = mean(m_hp_rec,4);
disp('Finished filtering');

niftiwrite(dec_map_hp,[out_dir name '_hp.nii'])

%%
cd(out_dir)
save(['NS_output' name])

end