%%%%Load the MEG RDMs and correlate them with the target RDM (Category or Visual)

clc
clear all
close all

%% Load the matrices
basedir='D:\MEG_CategoryShape\MEG_DATA_1back\'

%Load VisualSearch target matrix
load('MDSresults_AllSubs_CHECK.mat', 'RTMatrix_AllSubs');
VisSearchMatrix=tril(RTMatrix_AllSubs); %taking the lower triangle

%Load Category target matrix
load('DM_Category.mat');
CategoryMatrix=tril(target_dsm_category); %taking the lower triangle

%Load Shape target matrix
load('DM_Shape.mat');
ShapeMatrix=tril(target_dsm_shape); %taking the lower triangle

%Load Texture 
load('MDSresults_18Subs_Patches.mat', 'RTMatrix_AllSubs');
PatchesMatrix=tril(RTMatrix_AllSubs);

%Load Outline
load('MDSresults_18Subs_Lines.mat', 'RTMatrix_AllSubs');
LinesMatrix=tril(RTMatrix_AllSubs);

%% Which corr do we run
GLM_style=2; %1: Category+VisSearch, 2:Category + Outline + Texture
ss=1:29; %1:15 - 1-back task; 16:29 - oddball task; 1:29 - all

num_timepoints=71; %or 90

%FINAL NOV 2017: 800 ms, correct files
filename=['Final_Nov17_DecAcc_800ms_grad_3GLM_BothExp']
           
%% My code, using glmfit function. Load individual matrices!
for t=1:num_timepoints
    
    for s=ss
        
        display(['Subject ',num2str(s),' - Processing...'])
        
        % DECODING ACCURACY, INDIV SUBS
        datapath=[basedir,'class_data\Nov2017_FINAL_800ms_DM_MAG_NoFilter_s',num2str(s)];
        data2=load(datapath);
        DM_timepoint=data2.DM_allsubs(:,:,t); %Load the MEG RDM for one timepoint

%         % MAHALANOBIS
%         datapath=[basedir,'class_data\Final_Nov2017_800ms_SMOOTH_AvSamp_Mahalanobis_MAG_s',num2str(s)]; %mahalanobis mag
%         data2=load(datapath);
%         DM_tmpt=data2.DM_SmoothMahal_alltime(:,:,t);
%         DM_timepoint=[DM_tmpt(:,1:15) zeros(1,16)']; % add a column of zeros because it lacked it

        % normalize data
        DM_timepoint_vector=(squareform(DM_timepoint))';
        DM_timepoint_zscore=cosmo_normalize(DM_timepoint_vector,'zscore');
        
        % GLM regression
        if GLM_style==1 %Category + VisualSearch
            Predictors(:,1)=(squareform(CategoryMatrix))';
            Predictors(:,2)=(squareform(VisSearchMatrix))';
            % normalize matrices
            Predictors_zscore(:,1)=cosmo_normalize(Predictors(:,1),'zscore');
            Predictors_zscore(:,2)=cosmo_normalize(Predictors(:,2),'zscore');
            
            Predictors_intercept=[ones(size(Predictors_zscore,1),1) Predictors_zscore];
            
            % get betas
            betas = Predictors_zscore \ DM_timepoint_zscore;
            betas_intercept = Predictors_intercept \ DM_timepoint_zscore;
            Betas_Cat(s,t)=betas_intercept(2); %first element is a constant, starting from 2 are predictor betas
            Betas_VS(s,t)=betas_intercept(3);
            

        elseif GLM_style==2 %Category + Outline + Texture
            Predictors(:,1)=(squareform(CategoryMatrix))';
            Predictors(:,2)=(squareform(LinesMatrix))';
            Predictors(:,3)=(squareform(PatchesMatrix))';
            % normalize matrices
            Predictors_zscore(:,1)=cosmo_normalize(Predictors(:,1),'zscore');
            Predictors_zscore(:,2)=cosmo_normalize(Predictors(:,2),'zscore');
            Predictors_zscore(:,3)=cosmo_normalize(Predictors(:,3),'zscore');
            
            %Run the regression
            betas = Predictors_zscore \ DM_timepoint_zscore;
            
            Betas_Cat(s,t)=betas(1); 
            Betas_Outline(s,t)=betas(2);
            Betas_Texture(s,t)=betas(3);
            
           
        end

        clear('Predictors', 'Predictors_zscore', 'Betas','DM_timepoint_vector','DM_timepoint','DM_timepoint_zscore');
    end
    
end
% 
% % For subs 16-29: remove first 15 lines to exclude zeros
% if ss==16:29
%     if GLM_style==1
%         Betas_VS(1:15,:)=[];
%         Betas_Cat(1:15,:)=[];
%     elseif GLM_style==2
%         Betas_Cat(1:15,:)=[];
%         Betas_Outline(1:15,:)=[];
%         Betas_Texture(1:15,:)=[];
%     end
% end



% Now we have 3 matrices 9x71, each row is a subject and each column is a
% timepoint
%% We can average across subjects and get one timecourse.

if GLM_style==1
    Betas_Category_mean=mean(Betas_Cat,1);
    Betas_VisSearch_mean=mean(Betas_VS,1);
    %Find the peaks
    max_Cat=max(Betas_Category_mean)
    max_Cat_time=find(Betas_Category_mean==max_Cat)
    
    max_VS=max(Betas_VisSearch_mean)
    max_VS_time=find(Betas_VisSearch_mean==max_VS)
    
elseif GLM_style==2
    Betas_Category_mean=mean(Betas_Cat,1);
    Betas_Outline_mean=mean(Betas_Outline,1);
    Betas_Texture_mean=mean(Betas_Texture,1);
    %Find the peaks
    max_Cat=max(Betas_Category_mean)
    max_Cat_time=find(Betas_Category_mean==max_Cat)
    
    max_Outline=max(Betas_Outline_mean)
    max_Outline_time=find(Betas_Outline_mean==max_Outline)
    
    max_Texture=max(Betas_Texture_mean)
    max_Texture_time=find(Betas_Texture_mean==max_Texture)
    
end


%%Get Betas from 100-300 ms for the follow-up analysis
if GLM_style==1
    %Select 100-300ms window
    Betas_Cat_100to300=mean((Betas_Cat(:,20:40)),2);
    Betas_VS_100to300=mean((Betas_VS(:,20:40)),2);
    save('Betas_2GLM','Betas_Cat_100to300','Betas_VS_100to300')
    
elseif GLM_style==2
    Betas_Cat_100to300=mean((Betas_Cat(:,20:40)),2);
    Betas_Outline_100to300=mean((Betas_Outline(:,20:40)),2);
    Betas_Texture_100to300=mean((Betas_Texture(:,20:40)),2);
    save('Betas_3GLM','Betas_Cat_100to300','Betas_Outline_100to300','Betas_Texture_100to300')
    
end


%% Category and VS together
set(0,'DefaultAxesFontSize', 20)

if GLM_style==1
    
    % Significance testing with TFCE 
    
    % Put matrices in datasets
    ds_cat=struct(); ds_cat.samples=Betas_Cat;
    ds_vs=struct(); ds_vs.samples=Betas_VS;

    % insert time dimension
    time_axis=-.1:.01:.6;
    ds_cat=cosmo_dim_insert(ds_cat,2,0,{'time'},{time_axis},{1:num_timepoints});
    ds_vs=cosmo_dim_insert(ds_vs,2,0,{'time'},{time_axis},{1:num_timepoints});
    
    % Make a temporal neighborhood for clustering
    allow_clustering_over_time=false; % true or false
    
    % define the neighborhood
    nh_cl_cat=cosmo_cluster_neighborhood(ds_cat,'time',allow_clustering_over_time);
    nh_cl_vs=cosmo_cluster_neighborhood(ds_vs,'time',allow_clustering_over_time);
    
    % set subject information
    n_samples=size(ds_cat.samples,1);
    ds_cat.sa.chunks=(1:n_samples)';     % all subjects are independent
    ds_cat.sa.targets=ones(n_samples,1); % one-sample t-test
    
    ds_vs.sa.chunks=(1:n_samples)';
    ds_vs.sa.targets=ones(n_samples,1);
    
    % set clustering
    opt=struct();
    opt.h0_mean=0; % expected mean against which t-test is run
    opt.niter=10000; % 10,000 is recommdended for publication-quality analyses
    
    % run TFCE Monte Carlo multiple comparison correction
    % in the output map, z-scores above 1.65 (for one-tailed) or 1.96 (for
    % two-tailed) tests are significant
    ds_cat_tfce=cosmo_montecarlo_cluster_stat(ds_cat,nh_cl_cat,opt);
    ds_vs_tfce=cosmo_montecarlo_cluster_stat(ds_vs,nh_cl_vs,opt);
    
    % get significant time points for plotting. 1-tailed: |z|>1.6449, 2-tailed:
    % |z|>1.96
    signif_cat=zeros(1,num_timepoints); signif_vs=zeros(1,num_timepoints); 
    signif_cat(find(ds_cat_tfce.samples>1.6449))=1; %signif_cat(find(ds_cat_tfce.samples<-1.6449))=1;
    signif_vs(find(ds_vs_tfce.samples>1.6449))=1; %signif_vs(find(ds_vs_tfce.samples<-1.6449))=1;
%     
    %% Plot the figure
    e_cat=std(Betas_Cat,0,1)/sqrt(length(ss));
    plot_data_cat=Betas_Category_mean;
    e_vs=std(Betas_VS,0,1)/sqrt(length(ss));
    plot_data_vs=Betas_VisSearch_mean;

    figure
    set(gcf,'color','w')
    x_scale=linspace(-0.1,0.6,71);
    h_cat=boundedline(x_scale,plot_data_cat,e_cat,'cmap',[0,0,1],'transparency', 0.5,'alpha');
    set(h_cat,'linewidth',2);
    hold on
    h_vs=boundedline(x_scale,plot_data_vs,e_vs,'cmap',[1,0,0],'transparency', 0.5,'alpha');
    set(h_vs,'linewidth',2);
    hold on
    xlim([-0.1 0.6]);
    set(gca,'XTick',[-0.1:0.1:0.6]);
    ylim([-0.3 0.7]);
    yPos = 0;
    plot(get(gca,'xlim'), [yPos yPos],':','color','b','linewidth',1.5); %baseline
    hold on
    xPos=0;
    plot([xPos xPos],get(gca,'ylim'),'color','b','linewidth',1); %onset
    hold on
    Signif_cat=zeros(71,1)-5; Signif_VS=zeros(71,1)-5;
    Signif_cat(find(signif_cat))=0.57; Signif_VS(find(signif_vs))=0.6;
    plot(linspace(-0.1,0.6,71),Signif_cat,'o','color','b','linewidth',1);% signif for categ
    hold on
    plot(linspace(-0.1,0.6,71),Signif_VS,'o','color','r','linewidth',1); % signif for VS
    xlabel('time')
    ylabel('beta estimate')
    title('GLM Cat and VS, N=',num2str(length(ss)))
    grid off
    axis square

end



%% Cat, Text and Outline together
if GLM_style==2
    
    e_cat=std(Betas_Cat,0,1)/sqrt(length(ss));
    plot_data_cat=Betas_Category_mean;
    e_out=std(Betas_Outline,0,1)/sqrt(length(ss));
    plot_data_out=Betas_Outline_mean;
    e_text=std(Betas_Texture,0,1)/sqrt(length(ss));
    plot_data_text=Betas_Texture_mean;
    
 
    % Significance testing with TFCE 
    
    % Put matrices in datasets
    ds_cat=struct(); ds_cat.samples=Betas_Cat;
    ds_out=struct(); ds_out.samples=Betas_Outline;
    ds_text=struct(); ds_text.samples=Betas_Texture;
    
    % insert time dimension
    time_axis=-.1:.01:.6;
    ds_cat=cosmo_dim_insert(ds_cat,2,0,{'time'},{time_axis},{1:num_timepoints});
    ds_out=cosmo_dim_insert(ds_out,2,0,{'time'},{time_axis},{1:num_timepoints});
    ds_text=cosmo_dim_insert(ds_text,2,0,{'time'},{time_axis},{1:num_timepoints});
    
    % Make a temporal neighborhood for clustering
    allow_clustering_over_time=false; % true or false
    
    % define the neighborhood
    nh_cl_cat=cosmo_cluster_neighborhood(ds_cat,'time',allow_clustering_over_time);
    nh_cl_out=cosmo_cluster_neighborhood(ds_out,'time',allow_clustering_over_time);
    nh_cl_text=cosmo_cluster_neighborhood(ds_text,'time',allow_clustering_over_time);
    
    % set subject information
    n_samples=size(ds_cat.samples,1);
    ds_cat.sa.chunks=(1:n_samples)';     % all subjects are independent
    ds_cat.sa.targets=ones(n_samples,1); % one-sample t-test
    
    ds_out.sa.chunks=(1:n_samples)';
    ds_out.sa.targets=ones(n_samples,1);
    
    ds_text.sa.chunks=(1:n_samples)';
    ds_text.sa.targets=ones(n_samples,1);
    
    % set clustering
    opt=struct();
    opt.h0_mean=0; % expected mean against which t-test is run
    opt.niter=10000; % 10,000 is recommdended for publication-quality analyses
    
    % run TFCE Monte Carlo multiple comparison correction
    % in the output map, z-scores above 1.65 (for one-tailed) or 1.96 (for
    % two-tailed) tests are significant
    ds_cat_tfce=cosmo_montecarlo_cluster_stat(ds_cat,nh_cl_cat,opt);
    ds_out_tfce=cosmo_montecarlo_cluster_stat(ds_out,nh_cl_out,opt);
    ds_text_tfce=cosmo_montecarlo_cluster_stat(ds_text,nh_cl_text,opt);
    
    % get significant time points for plotting. 1-tailed: |z|>1.6449, 2-tailed:
    % |z|>1.96
    signif_cat=zeros(1,num_timepoints); signif_out=zeros(1,num_timepoints); signif_text=zeros(1,num_timepoints);
    signif_cat(find(ds_cat_tfce.samples>1.6449))=1; %signif_cat(find(ds_cat_tfce.samples<-1.6449))=1;
    signif_out(find(ds_out_tfce.samples>1.6449))=1; %signif_out(find(ds_out_tfce.samples<-1.6449))=1;
    signif_text(find(ds_text_tfce.samples>1.6449))=1; %signif_text(find(ds_text_tfce.samples<-1.6449))=1;
    
    
    
%     %% Plot Category, Outline, Texture together
%     figure
%     set(gcf,'color','w')
%     x_scale=linspace(-0.1,0.6,71);
%     h_cat=boundedline(x_scale,plot_data_cat,e_cat,'cmap',[0,0,1],'transparency', 0.5,'alpha');
%     set(h_cat,'linewidth',2);
%     hold on
%     h_out=boundedline(x_scale,plot_data_out,e_out,'cmap',[0,1,0],'transparency', 0.5,'alpha');
%     set(h_out,'linewidth',2);
%     hold on
%     h_text=boundedline(x_scale,plot_data_text,e_text,'cmap',[1,1,0],'transparency', 0.5,'alpha');
%     set(h_text,'linewidth',2);
%     xlim([-0.1 0.6]);
%     set(gca,'XTick',[-0.1:0.1:0.6]);
%     ylim([-0.3 0.7]);
%     yPos = 0;
%     plot(get(gca,'xlim'), [yPos yPos],':','color','b','linewidth',1.5); %baseline
%     hold on
%     xPos=0;
%     plot([xPos xPos],get(gca,'ylim'),'color','b','linewidth',1); %onset
%     hold on
%     Signif_cat=zeros(71,1)-5; Signif_Out=zeros(71,1)-5; Signif_Text=zeros(71,1)-5;
%     Signif_cat(find(signif_cat))=0.57; Signif_Out(find(signif_out))=0.6;  Signif_Text(find(signif_text))=0.55;
%     plot(linspace(-0.1,0.6,71),Signif_cat,'o','color','b','linewidth',1);% signif for categ
%     hold on
%     plot(linspace(-0.1,0.6,71),Signif_Out,'o','color','g','linewidth',1); % signif for outline
%     hold on
%     plot(linspace(-0.1,0.6,71),Signif_Text,'o','color',[1,1,0],'linewidth',1); % signif for VS
%     xlabel('time')
%     ylabel('beta estimate')
%     title('GLM Cat+Out+Text,N=',num2str(length(ss)))
%     grid off
%     axis square

end

% Save files
if GLM_style==1
    eval(['save ',filename,'  Betas_Cat Betas_VS signif_cat signif_vs']);
elseif GLM_style==2
    eval(['save ',filename,'  Betas_Cat Betas_Outline Betas_Texture signif_cat signif_out signif_text']);
end

