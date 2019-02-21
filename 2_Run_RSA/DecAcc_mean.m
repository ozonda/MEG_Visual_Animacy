clear all
close all

ss=1:15;
basedir='D:\FilesRovereto\MATLAB\MEG_CategoryShape\MEG_DATA_1back\';

filename=['Final_Nov2017_MeanDecAcc_grad_Exp2Oneback']; %filename to save

% Mean decoding accuracy
for t=1:71
    
    for s=ss
        
        display(['Subject ',num2str(s),' - Processing...'])
        
        % DECODING ACCURACY, INDIV SUBS
        datapath=[basedir,'class_data\Nov2017_FINAL_800ms_DM_GRAD_NoFilter_s',num2str(s)];
        data2=load(datapath);
        DM_timepoint=data2.DM_allsubs(:,:,t); %Load the MEG RDM for one timepoint

       %Turn into vector
        DM_timepoint_vector=(squareform(DM_timepoint))';
        
        %calculate the mean off diagonal
        mean_DecAcc_subs(s,t)=mean(DM_timepoint_vector);
        
        clear('DM_timepoint_vector','DM_timepoint');

    end
    
end

% % For oddball task only: clear first 15 rows of matrix
% mean_DecAcc_subs(1:15,:)=[];

%Find the peaks
mean_DecAcc=mean(mean_DecAcc_subs,1)
max_dec_acc=max(mean_DecAcc)
max_time=find(mean_DecAcc==max_dec_acc)



%% Significance testing TFCE

% Put matrices in datasets
ds_dec_acc=struct(); ds_dec_acc.samples=mean_DecAcc_subs;

% insert time dimension
time_axis=-.1:.01:.6;
ds_dec_acc=cosmo_dim_insert(ds_dec_acc,2,0,{'time'},{time_axis},{1:71});

% Make a temporal neighborhood for clustering
allow_clustering_over_time=false; % true or false

% define the neighborhood
nh_cl_dec_acc=cosmo_cluster_neighborhood(ds_dec_acc,'time',allow_clustering_over_time);

% set subject information
n_samples=size(ds_dec_acc.samples,1);
ds_dec_acc.sa.chunks=(1:n_samples)';     % all subjects are independent
ds_dec_acc.sa.targets=ones(n_samples,1); % one-sample t-test

% set clustering
opt=struct();
opt.h0_mean=0.5; % expected mean against which t-test is run
opt.niter=10000; % 10,000 is recommdended for publication-quality analyses

% run TFCE Monte Carlo multiple comparison correction
% in the output map, z-scores above 1.65 (for one-tailed) or 1.96 (for
% two-tailed) tests are significant
ds_dec_acc_tfce=cosmo_montecarlo_cluster_stat(ds_dec_acc,nh_cl_dec_acc,opt);

% get significant time points for plotting. 1-tailed: |z|>1.6449, 2-tailed:
% |z|>1.96
signif_dec_acc=zeros(1,71); 
signif_dec_acc(find(ds_dec_acc_tfce.samples>1.6449))=1; %signif_cat(find(ds_cat_tfce.samples<-1.6449))=1;
start_signif=find(signif_dec_acc)

% Save the file!
eval(['save ',filename,'  mean_DecAcc_subs mean_DecAcc max_dec_acc max_time signif_dec_acc start_signif']);




% %% Calculate the mean over subs
% mean_DecAcc=mean(mean_DecAcc_subs,1);
% e_DecAcc=std(mean_DecAcc_subs,0,1)/sqrt(length(ss));
% 
% %Plot the time course of the mean decoding accuracy
% figure
% set(0,'DefaultAxesFontSize', 20)
% set(gcf,'color','w')
% x_scale=linspace(-0.1,0.6,71);
% h_cat=boundedline(x_scale,mean_DecAcc,e_DecAcc,'cmap',[0,0,1],'transparency', 0.5,'alpha');
% set(h_cat,'linewidth',2);
% hold on
% xlim([-0.1 0.6]) ;
% ylim([0.4 1]);
% set(gca,'XTick',[-0.1:0.1:0.6]);
% set(gca,'YTick',[0.5:0.1:1]);
% set(gca,'Linewidth',2);
% yPos = 0.5;
% plot(get(gca,'xlim'), [yPos yPos],':','color','b','linewidth',1.5); %baseline
% hold on
% xPos=0;
% plot([xPos xPos],get(gca,'ylim'),'color','b','linewidth',1); %onset
% hold on
% xPos_peak=1.3; %the peak
% plot([xPos_peak xPos_peak],get(gca,'ylim'),'color','r','linewidth',1); %onset
% % Signif_dec_acc=zeros(71,1)-5; 
% % Signif_dec_acc(find(signif_dec_acc))=0.9; 
% % plot(linspace(-0.1,0.6,71),Signif_dec_acc,'o','color','b','linewidth',1);% signif for categ %'MarkerFaceColor','b'
% xlabel('time')
% ylabel('decoding accuracy')
% % title('Mean pairwise decoding accuracy')
% grid off
% axis square

