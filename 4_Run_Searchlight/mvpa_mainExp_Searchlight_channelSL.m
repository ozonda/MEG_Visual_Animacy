% Run RSA Searchlight
% Define the neightborhood for each feature (channel), calculate pairwise
% neural dissimilarity in this neighborhood (usind decoding accuracy), 
% and model it as a combination of Category and Visual target matrices. 
% Visualize the results using topoplot

clc
clear all
close all

basedir='C:\Users\daria.proklova\Documents\MATLAB\MEG_CategoryShape\MEG_DATA_1back\';

%Load VisualSearch target matrix
load('MDSresults_AllSubs_CHECK.mat', 'RTMatrix_AllSubs');
VisSearchMatrix=RTMatrix_AllSubs;

load('MDSresults_18Subs_Patches.mat', 'RTMatrix_AllSubs');
PatchesMatrix=RTMatrix_AllSubs;

load('MDSresults_18Subs_Lines.mat', 'RTMatrix_AllSubs');
LinesMatrix=RTMatrix_AllSubs;

%Load Category target matrix
load('DM_Category.mat');
CategoryMatrix=target_dsm_category; 

ss=1:29;

d_channels=2; %1: all; 2: magnetometers; 3: gradiometers

classi=1; % 1:lda classification, 2:correlation

average_samples=1; %1:average single trials for more robust signal
nch=4; %Number of chunks
average_nr=0.25; %Which percentage to average
average_out=50; %How many new samples should created?

for s=ss
            
    display(['Subject ',num2str(s),' - Processing...'])   

    %Load data
    datapath=[basedir,'ERP\ERP_MainExp_oneback_nofilter_s',num2str(s)]; %not filtered
    data2=load(datapath);
    data1=data2.timelock;
    
    %Crop
    cfg=[];
    cfg.toilim=[-0.1,0.6];
    data1=ft_redefinetrial(cfg,data1);
    
    %Baseline Correct
    cfg=[];
    cfg.baseline=[-0.1,0];
    cfg.parameter={'trial'};
    data1=ft_timelockbaseline(cfg,data1); 

    % convert to cosmomvpa struct
    ds_tl=cosmo_meeg_dataset(data1); 
    
    % get rid of features with at least one NaN value across samples
    fa_nan_mask=sum(isnan(ds_tl.samples),1)>0;
    fprintf('%d / %d features have NaN\n', ...
                sum(fa_nan_mask), numel(fa_nan_mask));
    ds_tl=cosmo_slice(ds_tl, ~fa_nan_mask, 2);

    %Assign all possible target types
    ds_tl.sa.targets=data1.trialinfo(:,1);

    %Select the conditions of interest
    ds_tl=cosmo_slice(ds_tl,ismember(ds_tl.sa.targets,[1:16])); %Potentially recode them to be 1 or 2 if needed
    
    %Assign chunks
    ds_tl.sa.chunks=ones(size(ds_tl.sa.targets),1);

    ds_new=ds_tl; %a copy of this dataset with all the targets

    %Select channels (the old way)
    if d_channels>1
        if d_channels==2
            chanx='meg_axial';
        elseif d_channels==3
            chanx='meg_planar';
        end
        ax=cosmo_meeg_find_layout(ds_new,'chantype',chanx);
        ax_label2=ax.label(1:end-2,:);
        feature_msk=cosmo_dim_match(ds_new,'chan',ax_label2);
        ds_new=cosmo_slice(ds_new,feature_msk,2);
    end

    %Average some samples
    ds_new=cosmo_average_samples(ds_new,'ratio',1)
    
    % Time x channel Searchlight code
    % set measure
    measure=@cosmo_target_dsm_corr_measure;
    measure_args=struct();
%     measure_args.glm_dsm={target_dsm_category, target_dsm_Patches, target_dsm_Lines, target_dsm_VisSearch}; %This is the regression searchlight part
    measure_args.glm_dsm={VisSearchMatrix, CategoryMatrix};
    
    % New: channel-time neighborhood!
    % define the neighborhood for each dimension
    chan_count=20;
    
    chan_nbrhood=cosmo_meeg_chan_neighborhood(ds_new,'count', chan_count,'chantype', chanx);
    time_nbrhood=cosmo_interval_neighborhood(ds_new,'time','radius',0);
    
    % cross neighborhoods for chan-time searchlight
    nbrhood=cosmo_cross_neighborhood(ds_new,{chan_nbrhood,time_nbrhood});
    
    % print some info
    nbrhood_nfeatures=cellfun(@numel,nbrhood.neighbors);
    fprintf('Features have on average %.1f +/- %.1f neighbors\n', ...
        mean(nbrhood_nfeatures), std(nbrhood_nfeatures));
    
    % only keep features with at least 10 neighbors
    center_ids=find(nbrhood_nfeatures>10);
    
    % Running searchlight analysis with the set paramenters
    res_ds=cosmo_searchlight(ds_new,nbrhood,measure,measure_args,'center_ids',center_ids);

    %split into 2 maps (for glm_dsm only)
    ds_map1=cosmo_slice(res_ds,1);
    filename1=([output_path,'22.09_MegSL_sub',num2str(s),'.nii']);
    % deduce layout from output
    layout=cosmo_meeg_find_layout(ds_map1);
    fprintf('The output uses layout %s\n', layout.name);
    res_map1=cosmo_map2meeg(ds_map1);
   

    ds_map2=cosmo_slice(res_ds,2);
    filename2=([output_path,'22.09_MegSL_sub',num2str(s),'.nii']);
    % deduce layout from output
    layout=cosmo_meeg_find_layout(ds_map2);
    fprintf('The output uses layout %s\n', layout.name);
    res_map2=cosmo_map2meeg(ds_map2);

    %Save data  
    save([basedir,'class_data\GLM_VS_Cat_s',num2str(s)],'res_ds'); 

end   

%% visualize timeseries results

% Map 1
figure();
cfg = [];
cfg.interactive = 'yes';
cfg.zlim=[-1 1];
cfg.layout       = layout;

% show figure with plots for each sensor
ft_multiplotER(cfg, res_map1);

% Map 2
figure();
cfg = [];
cfg.interactive = 'yes';
cfg.zlim=[-1 1];
cfg.layout       = layout;

% show figure with plots for each sensor
ft_multiplotER(cfg, res_map2);

%% visualize topology results
% show figure with topology for 0 to 600ms after stimulus onset in bins of
% 100 ms
figure();
cfg.xlim=-0.1:0.1:0.5;
ft_topoplotER(cfg, res_map1);

figure();
cfg.xlim=-0.1:0.1:0.5;
ft_topoplotER(cfg, res_map2);

