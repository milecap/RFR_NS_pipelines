% Author: Milena Capiglioni, University of Bern
% Contact: milena.capiglioni@extern.insel.ch
% Last update: Nov.2023

% If you use/modify this code for your future publication, please cite the
% corresponding article:  "Stimulus-Induced Rotary Saturation imaging of
%visually evoked neuroelectric response: preliminary results and data
%analysis" (currently under review)


function RFR_statistical(Subjects_folder,RFR_folder,ROI)

% RFR statistical analysis
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

%% Define patient folders

cd(Subjects_folder);
files = dir(Subjects_folder);
dirFlag = contains({files.name},{'Sub_'})&[files.isdir];
subjects = files(dirFlag);

output_dir = [Subjects_folder RFR_folder filesep]; % filesep 'S_circular_insula_ant'

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

%% ------------ Run statistical analysis for each patient
for sub = 1:length(subjects)

    work_dir = strcat(Subjects_folder,subjects(sub).name,filesep);
    current_subject = subjects(sub).name;

    folders = dir(work_dir);
    dirFlag = contains({folders.name},{RFR_folder})&[folders.isdir];
    func_folders = folders(dirFlag);
    if isempty(func_folders)
        disp(['Subject', current_subject, ' does not have any sequence of interest']);
        continue;
    end

    current_dir = [work_dir func_folders.name filesep];

    for sl = 1:2
        if sl == 1
            % Read full data matrix
            noStim = load([current_dir 'RFR_outputSLoff_noStim.mat']);
            VisStim = load([current_dir 'RFR_outputSLoff_VisStim.mat']);
            txt_sl{sl} = 'SLoff';
        else
            noStim = load([current_dir 'RFR_outputSLon_noStim.mat']);
            VisStim = load([current_dir 'RFR_outputSLon_VisStim.mat']);
            txt_sl{sl} = 'SLon';
        end

        % Define mask in ROI for noStim and VisStim data
        mask_noStim = zeros(size(noStim.seg_data));
        mask_VisStim = zeros(size(VisStim.seg_data));
        switch ROI % numbers defined according to Destrieux atlas
            case 'V1'
                mask_noStim(noStim.seg_data == 11143) = 1;% V1 left side
                mask_noStim(noStim.seg_data == 12143) = 1;% V1 rigth side
                mask_VisStim(VisStim.seg_data == 11143) = 1;% V1 left side
                mask_VisStim(VisStim.seg_data == 12143) = 1;% V1 rigth side
            case 'G_subcallosal'
                mask_noStim(noStim.seg_data == 11132) = 1;% G_subcallosal left side
                mask_noStim(noStim.seg_data == 12132) = 1;% G_subcallosal rigth side
                mask_VisStim(VisStim.seg_data == 11132) = 1;% G_subcallosal left side
                mask_VisStim(VisStim.seg_data == 12132) = 1;% G_subcallosal rigth side
            case 'S_circular_insula_ant'
                mask_noStim(noStim.seg_data == 11148) = 1;% S_circular_insula_ant left side
                mask_noStim(noStim.seg_data == 12148) = 1;% S_circular_insula_ant rigth side
                mask_VisStim(VisStim.seg_data == 11148) = 1;% S_circular_insula_ant left side
                mask_VisStim(VisStim.seg_data == 12148) = 1;% S_circular_insula_ant rigth side
        end
        mask_noStim(mask_noStim == 0) = NaN;
        mask_VisStim(mask_VisStim == 0) = NaN;

        % (Optional) Count number of voxels in ROI
        for slice = 1:size(mask_noStim,3)
            voxel_per_roi_noStim(sub,slice) = nnz(~isnan(mask_noStim(:,:,slice)));
        end

        % (Optional) show ROI mask over low resolution anatomical data
        %         if show_im
        %             mask_noStim_mos = slices2mosaic(mask_noStim);
        %             figure(10)
        %             imshow(slices2mosaic(squeeze(noStim.func_data(:,:,:,1))),[]); hold on;
        %             green = cat(3, zeros(size(mask_noStim_mos)), ones(size(mask_noStim_mos)), zeros(size(mask_noStim_mos)));
        %             auxim = imshow(green);
        %             hold off;
        %             set(auxim,'AlphaData', 0.75*mask_noStim_mos)
        %         end

        if sl == 1
            anat = noStim.func_data(:,:,:,1);
        end

        std_array = std(noStim.func_data,[],4,'omitnan');   % no stim data before pp
        mean_std_0(sub,sl,1) = mean(std_array(~isnan(mask_noStim)));
        std_std_0(sub,sl,1) = std(std_array(~isnan(mask_noStim)));
        std_array = std(VisStim.func_data,[],4,'omitnan');  % Vis stim data before pp
        mean_std_0(sub,sl,2) = mean(std_array(~isnan(mask_VisStim)));
        std_std_0(sub,sl,2) = std(std_array(~isnan(mask_VisStim)));

        std_array = mean(noStim.m_hp_rect,4,'omitnan');     % no stim data after pp
        mean_std_RFR(sub,sl,1) = mean(std_array(~isnan(mask_noStim)));
        std_std_RFR(sub,sl,1) = std(std_array(~isnan(mask_noStim)));
        std_array = mean(VisStim.m_hp_rect,4,'omitnan');    % vis stim data after pp
        mean_std_RFR(sub,sl,2) = mean(std_array(~isnan(mask_VisStim)));
        std_std_RFR(sub,sl,2) = std(std_array(~isnan(mask_VisStim)));

        for t = 1:size(noStim.func_data,4)
            data_t = noStim.func_data(:,:,:,t);
            V1_mean_noStim_0(sub,sl,t) = mean(data_t(~isnan(mask_noStim)));
            V1_std_noStim_0(sub,sl,t) = std(data_t(~isnan(mask_noStim)));

            data_t = VisStim.func_data(:,:,:,t);
            V1_mean_VisStim_0(sub,sl,t) = mean(data_t(~isnan(mask_VisStim)));
            V1_std_VisStim_0(sub,sl,t) = std(data_t(~isnan(mask_VisStim)));

            data_t = noStim.m_hp_rect(:,:,:,t);
            V1_mean_noStim_RFR(sub,sl,t) = mean(data_t(~isnan(mask_noStim)));
            V1_std_noStim_RFR(sub,sl,t) = std(data_t(~isnan(mask_noStim)));

            data_t = VisStim.m_hp_rect(:,:,:,t);
            V1_mean_VisStim_RFR(sub,sl,t) = mean(data_t(~isnan(mask_VisStim)));
            V1_std_VisStim_RFR(sub,sl,t) = std(data_t(~isnan(mask_VisStim)));
        end
        clear data_t

        for step = 1:2

            if step == 1 % before procesing
                noStim_map = noStim.ft_amp_0;
                VisStim_map = VisStim.ft_amp_0;
                txt_step{1} = '0';
            else %after processing
                noStim_map = noStim.ft_amp_RFR;
                VisStim_map = VisStim.ft_amp_RFR;
                txt_step{2} = 'RFR';
            end

            % NoStim in ROI
            act_map_offres = squeeze(noStim_map(:,:,:));
            values_offres = act_map_offres(mask_noStim == 1);
            mean_offres(sub,sl,step) = mean(values_offres);
            std_offres(sub,sl,step) = std(values_offres);

            % VisStim in ROI
            act_map_onres = squeeze(VisStim_map(:,:,:));
            values_onres = act_map_onres(mask_VisStim == 1);
            mean_onres(sub,sl,step) = mean(values_onres);
            std_onres(sub,sl,step) = std(values_onres);

            % Perform t-test in ROI
            [h(sub,sl,step),p(sub,sl,step),ci{sub,sl,step},statsf{sub,sl,step}] = ttest(values_onres,values_offres,'Tail','right');

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%% t-test per slice
            for slice = 1:6
                % NoStim in ROI
                act_map_offres = squeeze(noStim_map(:,:,slice));
                values_offres = act_map_offres(mask_noStim(:,:,slice) == 1);
                m_values_offres(sub,sl,step,slice) = mean(values_offres);
                v_offres{sub,sl,step,slice} = values_offres;

                % VisStim in ROI
                act_map_onres = squeeze(VisStim_map(:,:,slice));
                values_onres = act_map_onres(mask_VisStim(:,:,slice) == 1);
                m_values_onres(sub,sl,step,slice) = mean(values_onres);
                v_onres{sub,sl,step,slice} = values_onres;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end

    end
    %     waitforbuttonpress;
    %     close all


    %%%%%%%%%%%%%%%%%%%%%% Figure 2 %%%%%%%%%%%%%%%%%%%%%%
    if show_im
        col = lines(4);
        f = figure;
%         f.Position = [2090 469 938 263];

        ax1 = subplot(2,3,[1 4]);
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

        ax2 = subplot(2,3,[2 5]);
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

        l = legend();
        l.ItemTokenSize = [15 16];
        grid on; box on;
        xline(samples,'k--','HandleVisibility','off');
        xlim([0 2*samples])
        xlabel('Sample number');
        ylabel('Amplitude (a.u.)');
        title('Before RFR')

        ax3 = subplot(2,3,3);
        hold all;

        [x2,inBetween] = error_area(squeeze(V1_mean_noStim_RFR(sub,1,:))',squeeze(V1_std_noStim_RFR(sub,1,:))',1:samples);
        fill(x2, inBetween, col(1,:),'EdgeColor',col(1,:),'EdgeAlpha',0.3,'FaceAlpha',0.2,'HandleVisibility','off');
        plot(1:samples,squeeze(V1_mean_noStim_RFR(sub,1,:)),'color',col(1,:),'DisplayName','SL off');

        [x2,inBetween] = error_area(squeeze(V1_mean_VisStim_RFR(sub,1,:))',squeeze(V1_std_VisStim_RFR(sub,1,:))',samples + 1: 2*samples);
        fill(x2, inBetween, col(1,:),'EdgeColor',col(1,:),'EdgeAlpha',0.3,'FaceAlpha',0.2,'HandleVisibility','off');
        plot(samples + 1: 2*samples,squeeze(V1_mean_VisStim_RFR(sub,1,:)),'color',col(1,:),'HandleVisibility','off');

        l = legend();
        l.ItemTokenSize = [15 16];
        grid on; box on;
        xline(samples,'k--','HandleVisibility','off');
        xlim([0 2*samples])
        xticks([]);
        ylabel('Amplitude (a.u.)');
        title('After RFR')

        ax4 = subplot(2,3,6);
        hold all
        [x2,inBetween] = error_area(squeeze(V1_mean_noStim_RFR(sub,2,:))',squeeze(V1_std_noStim_RFR(sub,2,:))',1:samples);
        fill(x2, inBetween, col(2,:),'EdgeColor',col(2,:),'EdgeAlpha',0.3,'FaceAlpha',0.2,'HandleVisibility','off');
        plot(1:samples,squeeze(V1_mean_noStim_RFR(sub,2,:)),'color',col(2,:),'DisplayName','SL on');

        [x2,inBetween] = error_area(squeeze(V1_mean_VisStim_RFR(sub,2,:))',squeeze(V1_std_VisStim_RFR(sub,2,:))',samples + 1: 2*samples);
        fill(x2, inBetween, col(2,:),'EdgeColor',col(2,:),'EdgeAlpha',0.3,'FaceAlpha',0.2,'HandleVisibility','off');
        plot(samples + 1: 2*samples,squeeze(V1_mean_VisStim_RFR(sub,2,:)),'color',col(2,:),'HandleVisibility','off');

        l = legend();
        l.ItemTokenSize = [15 16];
        grid on; box on;
        xline(samples,'k--','HandleVisibility','off');
        xlim([0 2*samples])
        xlabel('Sample number');
        ylabel('Amplitude (a.u.)');

        linkaxes([ax3,ax4],'y')

        sgtitle([current_subject ' - ' RFR_folder],'interpreter','none');

        ax4.Position(4) = ax4.Position(4) + 0.03;
        ax3.Position(4) = ax4.Position(4);

        aux = ax1.Position;
        ax2.Position(2) = ax4.Position(2);
        ax2.Position(4) = ax3.Position(4)*2 + ax3.Position(2)-(ax4.Position(2)+ax4.Position(4)); %equal to (ax3.Position(2)+ ax3.Position(4)) - ax4.Position(2)

        ax1.Position = aux;
        ax1.Position(2) = ax2.Position(2);
        ax2.FontSize = ax4.FontSize;

%         waitforbuttonpress;

        if save_im
            saveas(f,[output_dir 'timecourse_' current_subject '.svg'])
        end

        close(f);
    end
    %%%%%%%%%%%%%%%%% end figure
end


%%

cd(Subjects_folder)

SequenceName  = [txt_sl{sl} '_' txt_step{step}];
SubjectName = {subjects.name}';
SLoff_0_t = squeeze(h(:,1,1));
SLoff_RFR_t = squeeze(h(:,1,2));
SLon_0_t = squeeze(h(:,2,1));
SLon_RFR_t = squeeze(h(:,2,2));
SLoff_0_p = squeeze(p(:,1,1));
SLoff_RFR_p = squeeze(p(:,1,2));
SLon_0_p = squeeze(p(:,2,1));
SLon_RFR_p = squeeze(p(:,2,2));

T2 = table(SubjectName,SLoff_0_t,SLoff_RFR_t,SLon_0_t,SLon_RFR_t,SLoff_0_p,SLoff_RFR_p,SLon_0_p,SLon_RFR_p);
writetable(T2,[output_dir RFR_folder ROI '_allslices_allsubjects_all.xls'],'WriteVariableNames',1);


%%
if show_im
    %%%%%%%%%%%%%%%%%%%% Subjects before and after RFR %%%%%%%%%%%%%%%%%
   
    f = figure(101);
    subplot(2,1,1)
    data = [mean_std_RFR(:,1,1) mean_std_RFR(:,1,2)];
    b = bar(data);
    xlim([0 13]);grid on;
    xlabel('Subject number'); ylabel('Amplitude (a.u.)');title('SL off')
    legend('noStim','VisStim')
    subplot(2,1,2)
    data = [mean_std_RFR(:,2,1) mean_std_RFR(:,2,2)];
    bar(data);
    xlim([0 13]); grid on;
    xlabel('Subject number'); ylabel('Amplitude (a.u.)'); title('SL on')
    legend('noStim','VisStim')

    sgtitle(RFR_folder,'interpreter','none')

    if save_im
        saveas(f, [output_dir RFR_folder 'Barplot_allsubjects_on_off.svg'])
    end

    %%%%%%%%%%%%%%%%%%%% Slice analysis for positive subjects %%%%%%%%%%%%%%%%%
    data_idx = [1 2 3 4];

    % Swarmchart
    clear x2
    for slice = 1:6
        y1{slice} = cat(1,v_offres{data_idx,1,2,slice})';
        x1{slice} = (slice - 0.1).*ones(1,length(y1{slice}));

        y2{slice} = cat(1,v_offres{data_idx,2,2,slice})';
        x2{slice} = (slice - 0.1).*ones(1,length(cat(1,v_offres{data_idx,2,2,slice})));

        y3{slice} = cat(1,v_onres{data_idx,1,2,slice})';
        x3{slice} = (slice + 0.1).*ones(1,length(cat(1,v_onres{data_idx,1,2,slice})));

        y4{slice} = cat(1,v_onres{data_idx,2,2,slice})';
        x4{slice} = (slice + 0.1).*ones(1,length(cat(1,v_onres{data_idx,2,2,slice})));
    end

    f = figure();
    s(1) = subplot(2,1,1); title('SL off'); hold all
    for as = 1:6
        x = [x1{as}];
        y = [y1{as}];
        ss_noStim = swarmchart(x,y,4,'MarkerEdgeColor','none','MarkerFaceColor','b','XJitterWidth',0.5,'HandleVisibility','off');

        x = [x3{as}];
        y = [y3{as}];
        ss_VisStim = swarmchart(x,y,4,'MarkerEdgeColor','none','MarkerFaceColor','r','XJitterWidth',0.5,'HandleVisibility','off');
    end

    ss_noStim.HandleVisibility = 'on'; ss_noStim.DisplayName = 'noStim';
    ss_VisStim.HandleVisibility = 'on'; ss_VisStim.DisplayName = 'VisStim';

    %     legend()

    s(2) = subplot(2,1,2); title('SL on'); hold all
    for as = 1:6
        x = [x2{as}];
        y = [y2{as}];
        ss = swarmchart(x,y,4,'MarkerEdgeColor','none','MarkerFaceColor','b','XJitterWidth',0.5);

        x = [x4{as}];
        y = [y4{as}];
        ss = swarmchart(x,y,4,'MarkerEdgeColor','none','MarkerFaceColor','r','XJitterWidth',0.5);
    end

    grid(s,'on'); box(s,'on')
    xlabel(s,'Slice number')
    ylabel(s,'Contrast amplitude (a.u.)')
    linkaxes(s,'y')
    sgtitle(RFR_folder,'interpreter','none')

%     if save_im
%         saveas(f, [output_dir RFR_folder 'Swarmchart_slices.svg'])
%     end

    % overlay violin plot

    % figure() uncomment to make a separate figure
    x = [x1{1}(1) x1{2}(1) x1{3}(1) x1{4}(1) x1{5}(1) x1{6}(1) x3{1}(1) x3{2}(1) x3{3}(1) x3{4}(1) x3{5}(1) x3{6}(1)];
    y = {y1{1}' y1{2}' y1{3}' y1{4}' y1{5}' y1{6}' y3{1}' y3{2}' y3{3}' y3{4}' y3{5}' y3{6}'};
    subplot(2,1,1)
    hold all
    violin_2(y,'x',x,'mc',[],'medc',[],'facecolor',['b';'b';'b';'b';'b';'b';'r';'r';'r';'r';'r';'r'],'facealpha',0.3);

    x = [x2{1}(1) x2{2}(1) x2{3}(1) x2{4}(1) x2{5}(1) x2{6}(1) x4{1}(1) x4{2}(1) x4{3}(1) x4{4}(1) x4{5}(1) x4{6}(1)];
    y = {y2{1}' y2{2}' y2{3}' y2{4}' y2{5}' y2{6}' y4{1}' y4{2}' y4{3}' y4{4}' y4{5}' y4{6}'};
    subplot(2,1,2)
    hold all
    violin_2(y,'x',x,'mc',[],'medc',[],'facecolor',['b';'b';'b';'b';'b';'b';'r';'r';'r';'r';'r';'r'],'facealpha',0.3);

    set(s,'FontSize',9)
    set(s,'XLim',[0.1 6.8])
    set(s,'YLim',[-20 100])

    if save_im
        saveas(f, [output_dir RFR_folder 'Violinplot_slices.svg'])
    end
end

%% code functions

    function [x2,inBetween] = error_area(mean,std,x)
    curve1 = mean + std;
    curve2 = mean - std;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    end

end