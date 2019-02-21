clc
clear all
close all

basedir='D:\FilesRovereto\MATLAB\MEG_CategoryShape\MEG_DATA_1back\';

ss=1:29;

d_channels=2; %1: all; 2: magnetometers; 3: gradiometers

classi=1; % 1:lda classification, 2:correlation

average_samples=1; %1:average single trials for more robust signal
nch=4; %Number of chunks
average_nr=0.25; %Which percentage to average
average_out=500; %How many new samples should created?


for s=1:29;
    
    display(['Subject ',num2str(s),' - Processing...'])
    DM_allsubs=zeros(16,16,90);
    
    for j=1:16
        for k=1:16
            if j>k
                
                %Load data
                datapath=[basedir,'ERP\ERP_MainExp_OneBack_Nov2017_NoFilter_s',num2str(s)];
                data2=load(datapath);
                data1=data2.timelock;
                
                %Crop
                cfg=[];
                cfg.toilim=[-0.1,0.8];
                data1=ft_redefinetrial(cfg,data1);
                
                %Baseline Correct
                cfg=[];
                cfg.baseline=[-0.1,0];
                cfg.parameter={'trial'};
                data1=ft_timelockbaseline(cfg,data1);
                
                % convert to cosmomvpa struct
                ds_tl=cosmo_meeg_dataset(data1);
                
                %Assign all possible target types
                ds_tl.sa.targets=data1.trialinfo(:,1)
                
                %Select the conditions of interest
                ds_tl=cosmo_slice(ds_tl,ismember(ds_tl.sa.targets,[j,k])); %select a pair of conditions
                
                %Balance the number of trials per condition
                if sum(ds_tl.sa.targets==j)>sum(ds_tl.sa.targets==k)
                    i=ds_tl.sa.targets==j;
                    takeouts=find(i==1);
                    takeouts=takeouts(randperm(length(takeouts)));
                    takeouts=takeouts(1:(sum(ds_tl.sa.targets==j)-sum(ds_tl.sa.targets==k)));
                    ds_tl.sa.targets(takeouts)=0;
                elseif sum(ds_tl.sa.targets==j)<sum(ds_tl.sa.targets==k)
                    i=ds_tl.sa.targets==k;
                    takeouts= find(i==1);
                    takeouts = takeouts(randperm(length(takeouts)));
                    takeouts=takeouts(1:(sum(ds_tl.sa.targets==k)-sum(ds_tl.sa.targets==j)));
                    ds_tl.sa.targets(takeouts)=0;
                end
                ds_tl=cosmo_slice(ds_tl,ds_tl.sa.targets>0);
                
                %Create Chunks
                if average_samples==1
                    ds_tl.sa.chunks=[1:length(ds_tl.sa.targets)]';
                    ds_tl.sa.chunks=cosmo_chunkize(ds_tl,nch);
                else
                    ds_tl.sa.chunks=zeros(size(ds_tl.sa.targets));
                    ds_tl.sa.chunks(ds_tl.sa.targets==1)=Shuffle(1:length(ds_tl.sa.chunks(ds_tl.sa.targets==j)));
                    ds_tl.sa.chunks(ds_tl.sa.targets==2)=Shuffle(1:length(ds_tl.sa.chunks(ds_tl.sa.targets==k)));
                end
                
                %Select channels
                if d_channels>1 %% and if it is 1?
                    if d_channels==2
                        chanx='meg_axial';
                    elseif d_channels==3
                        chanx='meg_planar';
                    end
                    ax=cosmo_meeg_find_layout(ds_tl,'chantype',chanx);
                    ax_label2=ax.label(1:end-2,:);
                    feature_msk=cosmo_dim_match(ds_tl,'chan',ax_label2);
                    ds_tl=cosmo_slice(ds_tl,feature_msk,2);
                end
                
                %Average some samples
                if average_samples==1
                    ds_tl=cosmo_average_samples(ds_tl,'ratio',average_nr,'nrep',average_out);
                end
                
                %Select classifier
                measure_args=struct();
                if classi==1
                    measure=@cosmo_crossvalidation_measure;
                    measure_args.classifier=@cosmo_classify_lda;
                    measure_args.partitions=cosmo_nchoosek_partitioner(ds_tl,1);
                else
                    measure=@cosmo_correlation_measure;
                    measure_args.partitions=cosmo_nchoosek_partitioner(ds_tl,nch/2);
                end
                
                % Computing neighbors for classification
                nbrhood=cosmo_interval_neighborhood(ds_tl,'time','radius',0);
                
                % Running searchlight analysis with the set paramenters
                res_ds=cosmo_searchlight(ds_tl,nbrhood,measure,measure_args);
                
                tp2avg = 3; % smooth results over n timepoints
                smoothed_acc=conv([res_ds.samples(1:2),res_ds.samples,res_ds.samples(end-1:end)],ones(1,tp2avg).*(1/tp2avg),'same');
                smoothed_acc=smoothed_acc(3:end-2);
                res_ds.samples=smoothed_acc;
                
                for ii=1:90
                    DM_allsubs(j,k,ii)=res_ds.samples(ii);
                end
                
            end
            
        end
    end
    
    save([basedir,'class_data\Nov2017_800ms_DM_MAG_NoFilter_s',num2str(s)],'DM_allsubs'); 
    
end

%%Same without filtering
    
% for s=1:2;
%     
%     display(['Subject ',num2str(s),' - Processing...'])
%     DM_allsubs=zeros(16,16,71); % changed from 71 to 91 timepoints (corresponding to 0.6 to 0.8 s)
%     
%     for j=1:16
%         for k=1:16
%             if j>k
%                 
%                 %Load data
%                 datapath=[basedir,'ERP\ERP_MainExp_oneback_nofilter_s',num2str(s)];
%                 data2=load(datapath);
%                 data1=data2.timelock;
%                 
%                 %Crop
%                 cfg=[];
%                 cfg.toilim=[-0.1,0.6]; %NEW 2017: changed from 0.6 to 0.8
%                 data1=ft_redefinetrial(cfg,data1);
%                 
%                 %Baseline Correct
%                 cfg=[];
%                 cfg.baseline=[-0.1,0];
%                 cfg.parameter={'trial'};
%                 data1=ft_timelockbaseline(cfg,data1);
%                 
%                 % convert to cosmomvpa struct
%                 ds_tl=cosmo_meeg_dataset(data1);
%                 
%                 %Assign all possible target types
%                 ds_tl.sa.targets=data1.trialinfo(:,1)
%                 
%                 %Select the conditions of interest
%                 ds_tl=cosmo_slice(ds_tl,ismember(ds_tl.sa.targets,[j,k])); %select a pair of conditions
%                 
%                 %Balance the number of trials per condition
%                 if sum(ds_tl.sa.targets==j)>sum(ds_tl.sa.targets==k)
%                     i=ds_tl.sa.targets==j;
%                     takeouts=find(i==1);
%                     takeouts=takeouts(randperm(length(takeouts)));
%                     takeouts=takeouts(1:(sum(ds_tl.sa.targets==j)-sum(ds_tl.sa.targets==k)));
%                     ds_tl.sa.targets(takeouts)=0;
%                 elseif sum(ds_tl.sa.targets==j)<sum(ds_tl.sa.targets==k)
%                     i=ds_tl.sa.targets==k;
%                     takeouts= find(i==1);
%                     takeouts = takeouts(randperm(length(takeouts)));
%                     takeouts=takeouts(1:(sum(ds_tl.sa.targets==k)-sum(ds_tl.sa.targets==j)));
%                     ds_tl.sa.targets(takeouts)=0;
%                 end
%                 ds_tl=cosmo_slice(ds_tl,ds_tl.sa.targets>0);
%                 
%                 %Create Chunks
%                 if average_samples==1
%                     ds_tl.sa.chunks=[1:length(ds_tl.sa.targets)]';
%                     ds_tl.sa.chunks=cosmo_chunkize(ds_tl,nch);
%                 else
%                     ds_tl.sa.chunks=zeros(size(ds_tl.sa.targets));
%                     ds_tl.sa.chunks(ds_tl.sa.targets==1)=Shuffle(1:length(ds_tl.sa.chunks(ds_tl.sa.targets==j))); %IS THIS CORRECT
%                     ds_tl.sa.chunks(ds_tl.sa.targets==2)=Shuffle(1:length(ds_tl.sa.chunks(ds_tl.sa.targets==k)));
%                 end
%                 
%                 %Select channels
%                 if d_channels>1
%                     if d_channels==2
%                         chanx='meg_axial';
%                     elseif d_channels==3
%                         chanx='meg_combined_from_planar';
%                     end
%                     ax=cosmo_meeg_find_layout(ds_tl,'chantype',chanx);
%                     ax_label2=ax.label(1:end-2,:);
%                     feature_msk=cosmo_dim_match(ds_tl,'chan',ax_label2);
%                     ds_tl=cosmo_slice(ds_tl,feature_msk,2);
%                 end
%                 
%                 %Average some samples
%                 if average_samples==1
%                     ds_tl=cosmo_average_samples(ds_tl,'ratio',average_nr,'repeats',average_out);
%                 end
%                 
%                 %Select classifier
%                 measure_args=struct();
%                 if classi==1
%                     measure=@cosmo_crossvalidation_measure;
%                     measure_args.classifier=@cosmo_classify_lda;
%                     measure_args.partitions=cosmo_nchoosek_partitioner(ds_tl,1);
%                 else
%                     measure=@cosmo_correlation_measure;
%                     measure_args.partitions=cosmo_nchoosek_partitioner(ds_tl,nch/2);
%                 end
%                 
%                 % Computing neighbors for classification
%                 nbrhood=cosmo_interval_neighborhood(ds_tl,'time','radius',0);
%                 
%                 % Running searchlight analysis with the set paramenters
%                 res_ds=cosmo_searchlight(ds_tl,nbrhood,measure,measure_args);
%                 
%                 tp2avg = 3; % smooth results over n timepoints
%                 smoothed_acc=conv([res_ds.samples(1:2),res_ds.samples,res_ds.samples(end-1:end)],ones(1,tp2avg).*(1/tp2avg),'same');
%                 smoothed_acc=smoothed_acc(3:end-2);
%                 res_ds.samples=smoothed_acc;
%                 
%                 for ii=1:71 %changed to 91
%                     DM_allsubs(j,k,ii)=res_ds.samples(ii);
%                 end
%                 
%             end
%             
%         end
%     end
%     
%     %         save([basedir,'class_data\Caos_OneBack_RSA_DM_grad_NoFilter_s',num2str(s),'rep',num2str(rep)],'DM_allsubs');
%     save([basedir,'class_data\Sept2017_RealityCheck_600ms_RSA_DM_ALLSENS_NoFilter_s',num2str(s)],'DM_allsubs');
%     
% end
% 
% 
