% Author: Milena Capiglioni, University of Bern
% Contact: milena.capiglioni@extern.insel.ch
% Last update: Nov.2023

% If you use/modify this code for your future publication, please cite the
% corresponding article:  "Stimulus-Induced Rotary Saturation imaging of
%visually evoked neuroelectric response: preliminary results and data
%analysis" (currently under review)


function NS_statistical(Subjects_folder,NS_folder,ROI)
% NS statistical analysis

show_im = 1; % show images
save_im = 1; % save images

if nargin == 0
    % No arguments provided, interactively select folder with subjects and
    % RFR output folder
    Subjects_folder = uigetdir('Select folder with subjects');
    NS_folder = input('Write name of NS processing output on the form NS_cof_n: ');
    ROI = questdlg('Which ROI do you want to analyse', 'ROI','V1','G_subcallosal','S_circular_insula_ant','V1');
else
    % Check the number of arguments provided
    if nargin < 2
        error('Not enough input arguments. Provide both Subjects_folder and NS_folder.');
    end
end

%% Define work directories

cd(Subjects_folder);
files = dir(Subjects_folder);
dirFlag = contains({files.name},{'Sub_'})&[files.isdir];
subjects = files(dirFlag);

output_dir = [Subjects_folder NS_folder filesep];

if ~isfolder(output_dir)
    mkdir(output_dir);
else
    if save_im
        prompt = 'output folder already exists! if you continue, you might rewrite content, do you want to continue? 0=NO 1=YES:  ';
        answer = input(prompt);
        if answer == 0
            error('Output folder already exists');
        end
    end
end


%% ------------ Run statistical analysis for each patient
for sub = 1:length(subjects)

    work_dir = strcat(Subjects_folder,subjects(sub).name,filesep);
    current_subject = subjects(sub).name;

    folders = dir;
    dirFlag = contains({folders.name},{NS_folder})&[folders.isdir];
    func_folders = folders(dirFlag);
    if isempty(func_folders)
        disp(['Subject', current_subject, ' does not have any sequence of interest']);
        continue;
    end

    current_dir = [work_dir func_folders.name filesep];

    % Read full data matrix
    noStim = load([current_dir 'NS_outputnoStim.mat']);
    VisStim = load([current_dir 'NS_outputVisStim.mat']);

    % Define mask in V1ROI
    mask_noStim = zeros(size(noStim.seg_off));
    mask_VisStim = zeros(size(VisStim.seg_off));
    switch ROI
        case 'V1'
            mask_noStim(noStim.seg_off == 11143) = 1;% V1 left side
            mask_noStim(noStim.seg_off == 12143) = 1;% V1 rigth side
            mask_VisStim(VisStim.seg_off == 11143) = 1;% V1 left side
            mask_VisStim(VisStim.seg_off == 12143) = 1;% V1 rigth side
        case 'G_subcallosal'
            mask_noStim(noStim.seg_off == 11132) = 1;% G_subcallosal left side
            mask_noStim(noStim.seg_off == 12132) = 1;% G_subcallosal rigth side
            mask_VisStim(VisStim.seg_off == 11132) = 1;% G_subcallosal left side
            mask_VisStim(VisStim.seg_off == 12132) = 1;% G_subcallosal rigth side
        case 'S_circular_insula_ant'
            mask_noStim(noStim.seg_off == 11148) = 1;% S_circular_insula_ant left side
            mask_noStim(noStim.seg_off == 12148) = 1;% S_circular_insula_ant rigth side
            mask_VisStim(VisStim.seg_off == 11148) = 1;% S_circular_insula_ant left side
            mask_VisStim(VisStim.seg_off == 12148) = 1;% S_circular_insula_ant rigth side
    end
    mask_noStim(mask_noStim == 0) = NaN;
    mask_noStim_mos = slices2mosaic(mask_noStim);
    mask_VisStim(mask_VisStim == 0) = NaN;
    mask_VisStim_mos = slices2mosaic(mask_VisStim);

    % (Optional)
    %         if show_im
    %             figure(10)
    %             subplot(1,2,1)
    %             imshow(slices2mosaic(squeeze(noStim.func_data(:,:,:,1))),[]); hold on;
    %             green = cat(3, zeros(size(mask_noStim_mos)), ones(size(mask_noStim_mos)), zeros(size(mask_noStim_mos)));
    %             auxim = imshow(green);
    %             hold off;
    %             set(auxim,'AlphaData', 0.75*mask_noStim_mos)
    %
    %             subplot(1,2,2)
    %             imshow(slices2mosaic(squeeze(VisStim.func_data(:,:,:,1))),[]); hold on;
    %             green = cat(3, zeros(size(mask_VisStim_mos)), ones(size(mask_VisStim_mos)), zeros(size(mask_VisStim_mos)));
    %             auxim = imshow(green);
    %             hold off;
    %             set(auxim,'AlphaData', 0.75*mask_VisStim_mos)
    %         end

    anat = noStim.func_off(:,:,:,1);

    std_array = mean(noStim.resid,4,'omitnan');
    mean_std_dec(sub,1) = mean(std_array(~isnan(mask_noStim)));
    std_std_dec(sub,1) = std(std_array(~isnan(mask_noStim)));
    std_array = mean(VisStim.resid,4,'omitnan');
    mean_std_dec(sub,2) = mean(std_array(~isnan(mask_VisStim)));
    std_std_dec(sub,2) = std(std_array(~isnan(mask_VisStim)));

    std_array = mean(noStim.m_hp_rec,4,'omitnan');
    mean_std_dec_hp(sub,1) = mean(std_array(~isnan(mask_noStim)));
    std_std_dec_hp(sub,1) = std(std_array(~isnan(mask_noStim)));
    std_array = mean(VisStim.m_hp_rec,4,'omitnan');
    mean_std_dec_hp(sub,2) = mean(std_array(~isnan(mask_VisStim)));
    std_std_dec_hp(sub,2) = std(std_array(~isnan(mask_VisStim)));



    for t = 1:size(noStim.func_off,4)
        data_t = noStim.func_off(:,:,:,t);
        V1_mean_noStim_0(sub,1,t) = mean(data_t(~isnan(mask_noStim)));
        V1_std_noStim_0(sub,1,t) = std(data_t(~isnan(mask_noStim)));
        data_t = noStim.func_on(:,:,:,t);
        V1_mean_noStim_0(sub,2,t) = mean(data_t(~isnan(mask_noStim)));
        V1_std_noStim_0(sub,2,t) = std(data_t(~isnan(mask_noStim)));

        data_t = VisStim.func_off(:,:,:,t);
        V1_mean_VisStim_0(sub,1,t) = mean(data_t(~isnan(mask_VisStim)));
        V1_std_VisStim_0(sub,1,t) = std(data_t(~isnan(mask_VisStim)));
        data_t = VisStim.func_on(:,:,:,t);
        V1_mean_VisStim_0(sub,2,t) = mean(data_t(~isnan(mask_VisStim)));
        V1_std_VisStim_0(sub,2,t) = std(data_t(~isnan(mask_VisStim)));

        data_t = noStim.resid(:,:,:,t);
        V1_mean_noStim_dec(sub,t) = mean(data_t(~isnan(mask_noStim)));
        V1_std_noStim_dec(sub,t) = std(data_t(~isnan(mask_noStim)));

        data_t = VisStim.resid(:,:,:,t);
        V1_mean_VisStim_dec(sub,t) = mean(data_t(~isnan(mask_VisStim)));
        V1_std_VisStim_dec(sub,t) = std(data_t(~isnan(mask_VisStim)));

        data_t = noStim.m_hp_rec(:,:,:,t);
        V1_mean_noStim_dec_hp(sub,t) = mean(data_t(~isnan(mask_noStim)));
        V1_std_noStim_dec_hp(sub,t) = std(data_t(~isnan(mask_noStim)));

        data_t = VisStim.m_hp_rec(:,:,:,t);
        V1_mean_VisStim_dec_hp(sub,t) = mean(data_t(~isnan(mask_VisStim)));
        V1_std_VisStim_dec_hp(sub,t) = std(data_t(~isnan(mask_VisStim)));
    end
    clear data_t

    for step = 1:2

        if step == 1 % dec
            noStim_map = noStim.dec_map;
            VisStim_map = VisStim.dec_map;
            txt_step{1} = 'dec';
        else % dec + hp
            noStim_map = noStim.dec_map_hp;
            VisStim_map = VisStim.dec_map_hp;
            txt_step{2} = 'hp';
        end

        % NoStim in V1ROI
        act_map_offres = squeeze(noStim_map(:,:,:));
        values_offres = act_map_offres(mask_noStim == 1);
        mean_offres(sub,step) = mean(values_offres);
        std_offres(sub,step) = std(values_offres);

        % VisStim in V1ROI
        act_map_onres = squeeze(VisStim_map(:,:,:));
        values_onres = act_map_onres(mask_VisStim == 1);
        mean_onres(sub,step) = mean(values_onres);
        std_onres(sub,step) = std(values_onres);

        % Perform t-test in ROI
        [h(sub,step),p(sub,step),~,~] = ttest(values_onres,values_offres,'Tail','right');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% tstatistics per slice
        %         for slice = 1:6
        %             % NoStim in V1ROI
        %             act_map_offres = squeeze(noStim_map(:,:,slice));
        %             values_offres = act_map_offres(mask_noStim(:,:,slice) == 1);
        %             m_values_offres(sub,step,slice) = mean(values_offres);
        %             v_offres{sub,step,slice} = act_map_offres(mask_noStim(:,:,slice) == 1);
        %             %                 quantity_offres(sub,step,slice) = length(act_map_offres(mask_noStim(:,:,slice) == 1));
        %
        %             % VisStim in V1ROI
        %             act_map_onres = squeeze(VisStim_map(:,:,slice));
        %             values_onres = act_map_onres(mask_VisStim(:,:,slice) == 1);
        %             m_values_onres(sub,step,slice) = mean(values_onres);
        %             v_onres{sub,step,slice} = act_map_onres(mask_VisStim(:,:,slice) == 1);
        %             %                 quantity_onres(sub,step,slice) = length(act_map_onres(mask_VisStim(:,:,slice) == 1));
        %
        %             %                 [hslice(sub,step,slice),pslice(sub,step,slice),~,~] = ttest(values_onres,values_offres,'Tail','right');
        %         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end


    %%%%%%%%%%%%%%%%%%%%%%% Figure paper
    if show_im
        col = lines(4);
        f=figure;

        ax(1) = nexttile;
        for slice = 1:size(anat,3)
            anat(:,:,slice) = imrotate(anat(:,:,slice),90);
            mask(:,:,slice) = imrotate(mask_noStim(:,:,slice),90);
        end
        mask_mos = slices2mosaic(mask);
        imshow(slices2mosaic(anat),[]); hold on;
        red = cat(3, ones(size(mask_mos)), zeros(size(mask_mos)), zeros(size(mask_mos)));
        auxim = imshow(red);
        hold off;
        set(auxim,'AlphaData', mask_mos)

        ax(2) = nexttile;
        samples = size(V1_mean_noStim_0,3);
        hold all;
        [x2,inBetween] = error_area(squeeze(V1_mean_noStim_0(sub,1,:))',squeeze(V1_std_noStim_0(sub,1,:))',1:samples);
        fill(x2, inBetween, col(1,:),'EdgeColor',col(1,:),'EdgeAlpha',0.3,'FaceAlpha',0.2,'HandleVisibility','off');
        plot(1:samples,squeeze(V1_mean_noStim_0(sub,1,:)),'color',col(1,:),'DisplayName','SL off');

        [x2,inBetween] = error_area(squeeze(V1_mean_VisStim_0(sub,1,:))',squeeze(V1_std_VisStim_0(sub,1,:))',samples + 1: 2*samples);
        fill(x2, inBetween, col(1,:),'EdgeColor',col(1,:),'EdgeAlpha',0.3,'FaceAlpha',0.2,'HandleVisibility','off');
        plot(samples + 1: 2*samples,squeeze(V1_mean_VisStim_0(sub,1,:)),'color',col(1,:),'HandleVisibility','off');

        [x2,inBetween] = error_area(squeeze(V1_mean_noStim_0(sub,2,:))',squeeze(V1_std_noStim_0(sub,2,:))',1:samples);
        fill(x2, inBetween, col(2,:),'EdgeColor',col(2,:),'EdgeAlpha',0.3,'FaceAlpha',0.2,'HandleVisibility','off');
        plot(1:samples,squeeze(V1_mean_noStim_0(sub,2,:)),'color',col(2,:),'DisplayName','SL on');

        [x2,inBetween] = error_area(squeeze(V1_mean_VisStim_0(sub,2,:))',squeeze(V1_std_VisStim_0(sub,2,:))',samples + 1: 2*samples);
        fill(x2, inBetween, col(2,:),'EdgeColor',col(2,:),'EdgeAlpha',0.3,'FaceAlpha',0.2,'HandleVisibility','off');
        plot(samples + 1: 2*samples,squeeze(V1_mean_VisStim_0(sub,2,:)),'color',col(2,:),'HandleVisibility','off');
        grid on; box on;
        xline(samples,'k--');
        xlim([0 2*samples])
        xlabel('Sample number');
        ylabel('Amplitude (a.u.)');
        title('Before RFR')

        ax(3) = nexttile;
        hold all;
        [x2,inBetween] = error_area(squeeze(V1_mean_noStim_dec(sub,:)),squeeze(V1_std_noStim_dec(sub,:)),1:samples);
        fill(x2, inBetween, col(3,:),'EdgeColor',col(3,:),'EdgeAlpha',0.3,'FaceAlpha',0.2,'HandleVisibility','off');
        plot(1:samples,squeeze(V1_mean_noStim_dec(sub,:))','color',col(3,:),'HandleVisibility','off');

        [x2,inBetween] = error_area(squeeze(V1_mean_VisStim_dec(sub,:)),squeeze(V1_std_VisStim_dec(sub,:)),samples + 1: 2*samples);
        fill(x2, inBetween, col(3,:),'EdgeColor',col(3,:),'EdgeAlpha',0.3,'FaceAlpha',0.2,'HandleVisibility','off');
        plot(samples + 1: 2*samples,squeeze(V1_mean_VisStim_dec(sub,:))','color',col(3,:),'HandleVisibility','off');
        grid on; box on;
        xlim([0 2*samples])
        xlabel('Sample number');
        ylabel('Amplitude (a.u.)');
        xline(samples,'k--');
        title('After dec')

        ax(4) = nexttile;
        hold all;
        [x2,inBetween] = error_area(squeeze(V1_mean_noStim_dec_hp(sub,:)),squeeze(V1_std_noStim_dec_hp(sub,:)),1:samples);
        fill(x2, inBetween, col(4,:),'EdgeColor',col(4,:),'EdgeAlpha',0.3,'FaceAlpha',0.2); plot(1:samples,squeeze(V1_mean_noStim_dec_hp(sub,:))','color',col(4,:));
        [x2,inBetween] = error_area(squeeze(V1_mean_VisStim_dec_hp(sub,:)),squeeze(V1_std_VisStim_dec_hp(sub,:)),samples + 1: 2*samples);
        fill(x2, inBetween, col(4,:),'EdgeColor',col(4,:),'EdgeAlpha',0.3,'FaceAlpha',0.2); plot(samples + 1: 2*samples,squeeze(V1_mean_VisStim_dec_hp(sub,:))','color',col(4,:));
        grid on; box on;
        xline(samples,'k--');
        xlim([0 2*samples])
        xlabel('Sample number');
        ylabel('Amplitude (a.u.)');
        title('After hp')

        sgtitle([current_subject ' - ' NS_folder],'interpreter','none');

        linkaxes([ax(3),ax(4)],'y')


        f.Position = [ 523         427        1268         284];

%         waitforbuttonpress;

        if save_im
            saveas(f,[output_dir 'timecourse_' current_subject '.svg'])
        end

        close(f);
    end
    %%%%%%%%%%%%%%%%%%end figure paper
end


cd(Subjects_folder)
for step = 1:2
    SequenceName  = [txt_step{step}];
    SubjectName = {subjects.name}';
    OnresMean = squeeze(mean_onres(:,step));
    OnresStd = squeeze(std_onres(:,step));
    OffresMean = squeeze(mean_offres(:,step));
    OffresStd = squeeze(std_offres(:,step));
    ttest_result = squeeze(h(:,step));
    pvalue = squeeze(p(:,step));

    T2 = table(SubjectName,OnresMean,OnresStd,OffresMean,OffresStd,ttest_result,pvalue);
    writetable(T2,[output_dir NS_folder ROI '_allslices_allsubjects_' SequenceName '.xls'],'WriteVariableNames',1);
end


% figure(100)
% subplot(2,1,1)
% x = 1:length(subjects);
% data = [mean_std_dec(:,1) mean_std_dec(:,2)];
% b = bar(data);
% xlim([0 13]);grid on;
% xlabel('Subject');ylabel('Amplitude (a.u.)');title('Before RFR')
% legend('noStim - SL off','VisStim - SL off','noStim - SL on','VisStim - SL on')
% 
% subplot(2,1,2)
% data = [mean_std_dec_hp(:,1) mean_std_dec_hp(:,2)];
% bar(data);
% xlim([0 13]); grid on;
% xlabel('Subject');ylabel('Amplitude (a.u.)'); title('After RFR')
% legend('noStim - SL off','VisStim - SL off','noStim - SL on','VisStim - SL on')
% sgtitle(NS_folder)


%%
function [x2,inBetween] = error_area(mean,std,x)
curve1 = mean + std;
curve2 = mean - std;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
end


end