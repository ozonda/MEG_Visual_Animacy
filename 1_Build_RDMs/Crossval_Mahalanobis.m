%% Creating RDMs using Mahalanobis distance

clc
clear all
close all

basedir='D:\FilesRovereto\MATLAB\MEG_CategoryShape\MEG_DATA_1back\';

d_channels=3; %1: all; 2: magnetometers; 3: gradiometers

% NEW
average_samples=1; %1:average single trials for more robust signal
nch=2; %Number of chunks
average_nr=0.25; %Which percentage to average
average_out=500; %How many new samples should created?


%% Use the RSA toolbox to get LDAt values for training set A and testing B

for t=1:90 %We run the entire script for each timepoint, to get the average RDM across subs for each time pont!
    
    for s=1:29
        
        % Load the MEG dataset (trials X channels)
        display(['Subject ',num2str(s),' - Processing...'])
        
        datapath=[basedir,'ERP\ERP_MainExp_OneBack_Nov2017_NoFilter_s',num2str(s)]; %not filtered
        data2=load(datapath);
        data1=data2.timelock;
        
        %Crop
        cfg=[];
        cfg.toilim=[-0.1,0.8]; %NEW 2017: changed from 0.6 to 0.8
        data1=ft_redefinetrial(cfg,data1);
        
        %Baseline Correct
        cfg=[];
        cfg.baseline=[-0.1,0];
        cfg.parameter={'trial'};
        data1=ft_timelockbaseline(cfg,data1);
        
        % convert to cosmomvpa struct
        ds_tl=cosmo_meeg_dataset(data1);
        
        %Assign all possible target types
        ds_tl.sa.targets=data1.trialinfo(:,1);
        
        %Select the conditions of interest
        ds_tl=cosmo_slice(ds_tl,ismember(ds_tl.sa.targets,[1:16])); %Potentially recode them to be 1 or 2 if needed
        
        ds_tl=cosmo_slice(ds_tl,ds_tl.sa.targets>0);
        
        % Assign chunks: 1 to odd trials, 2 to even (for the splitting of
        % the dataset in two halves later)
        for trial=1:size(ds_tl.samples,1)
            if mod(trial,2)==1 %odd trial
                ds_tl.sa.chunks(trial,:)=1;
            else %even trial
                ds_tl.sa.chunks(trial,:)=2
            end
        end
        
        % just to check everything is ok
        cosmo_check_dataset(ds_tl);
        
        % for subs 7 and 13 - remove the features that are zero
        indices=ones(1,size(ds_tl.samples,2)); %building a mask - bad features will have 0, the rest 1
        
        for i=1:size(ds_tl.samples,2)
            if ds_tl.samples(:,i)==zeros(size(ds_tl.samples,1),1)
               indices(i)=0;
            end
        end
        
        % Remove these features
        feature_msk=find(indices);
        ds_tl=cosmo_slice(ds_tl,feature_msk,2);
        
        %Select channels
        if d_channels>1
            if d_channels==2
                chanx='meg_axial';
            elseif d_channels==3
                chanx='meg_combined_from_planar';
            end
            ax=cosmo_meeg_find_layout(ds_tl,'chantype',chanx);
            ax_label2=ax.label(1:end-2,:);
            feature_msk=cosmo_dim_match(ds_tl,'chan',ax_label2);
            ds_tl=cosmo_slice(ds_tl,feature_msk,2);
        end
        
        %Average some samples (NEW)
        if average_samples==1
            ds_tl=cosmo_average_samples(ds_tl,'ratio',average_nr,'nrep',average_out);
        end
        
        % Split it in two (even/odd, or first/second half)
        splits=cosmo_split(ds_tl,'chunks');
        Dataset_train=splits{1};
        Dataset_test=splits{2};
        
        % Check if the number of trials is the same, in case remove one trial
        if size(Dataset_train.samples,1)>size(Dataset_test.samples,1)
            Dataset_train.sa.targets(end)=0; %remove the last target
            Dataset_train=cosmo_slice(Dataset_train,Dataset_train.sa.targets>0); %slice the dataset to remove last trial
        end
        
        %Slice both train and test over time point
        sliced_ds_train=cosmo_slice(Dataset_train,find(Dataset_train.fa.time==t),2);
        sliced_ds_test=cosmo_slice(Dataset_test,find(Dataset_test.fa.time==t),2);
        
        %Create design matrices for train and test set.
        %Trials x Conditions. If trial N belonged to condition K, matrix (n,k)=1, otherwise 0
        DM_train=zeros(size(sliced_ds_train.samples,1),16);
        DM_test=zeros(size(sliced_ds_test.samples,1),16);
        
        for n=1:size(sliced_ds_train.samples,1) %trials
            for k=1:16 %conditions
                indices_train=find(sliced_ds_train.sa.targets==k);
                DM_train(indices_train,k)=1;
                
                indices_test=find(sliced_ds_test.sa.targets==k);
                DM_test(indices_test,k)=1;
            end
        end
        
        %Run the toolbox function to get LDAt RDM for this timepoint only
        [RDM_fdtFolded_ltv, cv2RDM_fdt_sq] = fisherDiscrTRDM(DM_train,sliced_ds_train.samples,DM_test,sliced_ds_test.samples);
        RDM_lda = squareform(RDM_fdtFolded_ltv); % diagonals will contain zeros
        RDMs(s).RDM = RDM_lda;
        RDMs(s).name = ['LDAtRDM | subject',num2str(s)];
        RDMs(s).color = [1 0 0];
        
        % Save the indiv matrix again...
        Indiv_mahal_DM=RDM_lda;
        
        %Clean up
        clear DM_train
        clear DM_test
        clear Dataset_train
        clear Dataser_test
        clear sliced_ds_train
        clear sliced_ds_test
        clear ds_tl
        
        %Save indiv subject's matrix
        save(['Nov2017_800ms_AvSamp_Mahalanobis_GRAD_s',num2str(s),'_t',num2str(t)],'Indiv_mahal_DM')
    end
   
end

% New: load and save average rdms
for t=1:90
    for s=1:29
%         load(['February_AvSamp_Mahalanobis_MAG_s',num2str(s),'_t',num2str(t)],'Indiv_mahal_DM');
        load(['Final_Nov2017_800ms_AvSamp_Mahalanobis_GRAD_s',num2str(s),'_t',num2str(t)],'Indiv_mahal_DM');
        RDMs(s).RDM = Indiv_mahal_DM;
        RDMs(s).name = ['LDAtRDM | subject',num2str(s)];
        RDMs(s).color = [1 0 0];
    end
        
    averageRDMs_LDt(t) = averageRDMs_subjectSession(RDMs, 'subject');
end
    
% Save the averaged matrix
save(['Final_Nov2017_800ms_AvSamp_Mahalanobis_GRAD_ALLSUBS'],'averageRDMs_LDt')


