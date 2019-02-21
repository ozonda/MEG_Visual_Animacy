
% Loading and averaging Mahalanobis RDMs across 3 timepoints (smoothing).
clc
clear all
close all

basedir='D:\FilesRovereto\MATLAB\MEG_CategoryShape\MEG_DATA_1back\';

for s=1:29
    display(['Subject ',num2str(s),' - Processing...'])
    for i=1:16
        for j=1:16
            if i>j
                for t=1:90
                    %Load existing non-smoothed matrix for this t, this sub.
                    load(['Nov2017_800ms_AvSamp_Mahalanobis_MAG_s',num2str(s),'_t',num2str(t)],'Indiv_mahal_DM'); %mag
                    NewVector(t)=Indiv_mahal_DM(i,j);
                end
                
                %Now we average the vector across 3 timepoints
                tp2avg = 3; % smooth over n timepoints
                smoothed_Newvector=conv([NewVector(1:2),NewVector,NewVector(end-1:end)],ones(1,tp2avg).*(1/tp2avg),'same');
                smoothed_Newvector=smoothed_Newvector(3:end-2);
                Newvector=smoothed_Newvector;
                
                %Now we put it back into matrix for each t
                for ii=1:90
                    DM_SmoothMahal_alltime(i,j,ii)=Newvector(ii);
                end
            end
        end
    end
    %After we went through all pairs, save the matrix for this subject. 
    %It contains all time points
    save([basedir,'class_data\Nov2017_800ms_SMOOTH_AvSamp_Mahalanobis_MAG_s',num2str(s)],'DM_SmoothMahal_alltime'); %mag
end                   