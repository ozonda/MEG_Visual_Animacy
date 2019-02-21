close all
clear all

%Specify what type of GLM (2 or 3 predictors) and number of participants
GLM_style=2;
ss=1:29; 

xlimit=0.8
num_timepoints=90; %analyzing till 600 ms (71 tp) or 800 ms (90 tp)

% Load the files for Magnetometers and Gradiometers here. pay attention to
% GLM type (2 or 3 predictors)
load('Final_Nov17_MAHALANOBIS_800ms_grad_3GLM_BothExp');

Betas_Cat_grad=Betas_Cat; Betas_Outline_grad=Betas_Outline; Betas_Texture_grad=Betas_Texture;
signif_cat_grad=signif_cat; signif_out_grad=signif_out; signif_text_grad=signif_text;

load('Final_Nov17_MAHALANOBIS_800ms_mag_3GLM_BothExp');

%%For 3GLM
Betas_Cat_mag=Betas_Cat; Betas_Outline_mag=Betas_Outline; Betas_Texture_mag=Betas_Texture;
signif_cat_mag=signif_cat; signif_out_mag=signif_out; signif_text_mag=signif_text;

Betas_Category_grad_mean=mean(Betas_Cat_grad,1); Betas_Category_mag_mean=mean(Betas_Cat_mag,1);
Betas_Outline_grad_mean=mean(Betas_Outline_grad,1); Betas_Outline_mag_mean=mean(Betas_Outline_mag,1);
Betas_Texture_grad_mean=mean(Betas_Texture_grad,1); Betas_Texture_mag_mean=mean(Betas_Texture_mag,1);

% % %% Load Magnetometers and Gradiometers 2 GLM
% load('Final_Nov17_DecAcc_800ms_grad_2GLM_BothExp');
% Betas_Cat_grad=Betas_Cat; Betas_VS_grad=Betas_VS; 
% signif_cat_grad=signif_cat; signif_vs_grad=signif_vs; 

% load('Final_Nov17_DecAcc_800ms_mag_2GLM_BothExp');
% 
% %%For 2GLM
% 
% Betas_Cat_mag=Betas_Cat; Betas_VS_mag=Betas_VS; 
% signif_cat_mag=signif_cat; signif_vs_mag=signif_vs; 
% 
% Betas_Category_grad_mean=mean(Betas_Cat_grad,1); Betas_Category_mag_mean=mean(Betas_Cat_mag,1);
% Betas_VS_grad_mean=mean(Betas_VS_grad,1); Betas_VS_mag_mean=mean(Betas_VS_mag,1);





%% Plot Category and VS together
set(0,'DefaultAxesFontSize', 20)

if GLM_style==1


    
    %% Plot the figure
    e_cat_grad=std(Betas_Cat_grad,0,1)/sqrt(length(ss));  e_cat_mag=std(Betas_Cat_mag,0,1)/sqrt(length(ss));
    plot_data_cat_grad=Betas_Category_grad_mean; plot_data_cat_mag=Betas_Category_mag_mean;
    
    e_vs_grad=std(Betas_VS_grad,0,1)/sqrt(length(ss)); e_vs_mag=std(Betas_VS_mag,0,1)/sqrt(length(ss));
    plot_data_vs_grad=Betas_VS_grad_mean; plot_data_vs_mag=Betas_VS_mag_mean;
    
%     %Use this to remove error bars
%     error=zeros(1,num_timepoints);
    

    figure
    set(gcf,'color','w')
    x_scale=linspace(-0.1,xlimit,num_timepoints);
    set(gca,'linewidth',2);
    
    h_cat_grad=boundedline(x_scale,plot_data_cat_grad(1:num_timepoints),e_cat_grad,'cmap',[0,0,1],'transparency', 0.4,'alpha'); %plot category grad
    set(h_cat_grad,'linewidth',2);
    hold on
    h_cat_mag=boundedline(x_scale,plot_data_cat_mag(1:num_timepoints),e_cat_mag,'cmap',[0,0,1],'transparency', 0.2,'alpha','--'); %plot category mag
    set(h_cat_mag,'linewidth',2);
    hold on 
    
    h_vs_grad=boundedline(x_scale,plot_data_vs_grad(1:num_timepoints),e_vs_grad,'cmap',[1,0,0],'transparency', 0.4,'alpha'); %plot vs grad
    set(h_vs_grad,'linewidth',2);
    hold on
    h_vs_mag=boundedline(x_scale,plot_data_vs_mag(1:num_timepoints),e_vs_mag,'cmap',[1,0,0],'transparency', 0.2,'alpha','--'); % plot vs mag
    set(h_vs_mag,'linewidth',2);
    
    hold on
    xlim([-0.1 xlimit]); %or 0.8
%     set(gca,'XTick',[-0.1:0.1:xlimit]);
    set(gca,'XTick',[0:0.2:xlimit]);
    ylim([-0.2 0.7]);
    yPos = 0;
    plot(get(gca,'xlim'), [yPos yPos],':','color',[0 0 0],'linewidth',1.5); %baseline
    hold on
    xPos=0;
    plot([xPos xPos],get(gca,'ylim'),'color',[0 0 0],'linewidth',1); %onset
    hold on
    
    %Significance grad
    Signif_cat_grad=zeros(num_timepoints,1)-5; Signif_VS_grad=zeros(num_timepoints,1)-5;
    Signif_cat_grad(find(signif_cat_grad))=0.54; Signif_VS_grad(find(signif_vs_grad))=0.60;
    
    %Significance mag
    Signif_cat_mag=zeros(num_timepoints,1)-5; Signif_VS_mag=zeros(num_timepoints,1)-5;
    Signif_cat_mag(find(signif_cat_mag))=0.51; Signif_VS_mag(find(signif_vs_mag))=0.57;
    
    plot(linspace(-0.1,xlimit,num_timepoints),Signif_cat_grad(1:num_timepoints),'o','color','b','MarkerFaceColor','b','linewidth',1);% signif for categ grad
    hold on
    plot(linspace(-0.1,xlimit,num_timepoints),Signif_VS_grad(1:num_timepoints),'o','color','r','MarkerFaceColor','r','linewidth',1); % signif for VS grad
    hold on
    
    plot(linspace(-0.1,xlimit,num_timepoints),Signif_cat_mag(1:num_timepoints),'o','color','b','linewidth',1);% signif for categ mag
    hold on
    plot(linspace(-0.1,xlimit,num_timepoints),Signif_VS_mag(1:num_timepoints),'o','color','r','linewidth',1); % signif for VS mag
    
    xlabel('time')
    ylabel('beta estimate')
    grid off
    axis square

end



%% Cat, Text and Outline together
if GLM_style==2
    
    %Gradiometers
    e_cat_grad=std(Betas_Cat_grad,0,1)/sqrt(length(ss));
    plot_data_cat_grad=Betas_Category_grad_mean;
    e_out_grad=std(Betas_Outline_grad,0,1)/sqrt(length(ss));
    plot_data_out_grad=Betas_Outline_grad_mean;
    e_text_grad=std(Betas_Texture_grad,0,1)/sqrt(length(ss));
    plot_data_text_grad=Betas_Texture_grad_mean;
    
    %Magnetometers
    e_cat_mag=std(Betas_Cat_mag,0,1)/sqrt(length(ss));
    plot_data_cat_mag=Betas_Category_mag_mean;
    e_out_mag=std(Betas_Outline_mag,0,1)/sqrt(length(ss));
    plot_data_out_mag=Betas_Outline_mag_mean;
    e_text_mag=std(Betas_Texture_mag,0,1)/sqrt(length(ss));
    plot_data_text_mag=Betas_Texture_mag_mean;
    
%     %Use this to remove error bars
%     error=zeros(1,num_timepoints);
    
    
    
    %% Plot Category, Outline, Texture together
    figure
    set(gcf,'color','w')
    x_scale=linspace(-0.1,xlimit,num_timepoints); %was 0.6, 71
    set(gca,'linewidth',2);
    
    %Category
    h_cat_grad=boundedline(x_scale,plot_data_cat_grad,e_cat_grad,'cmap',[0,0,1],'transparency', 0.5,'alpha');
    set(h_cat_grad,'linewidth',2);
    hold on
    h_cat_mag=boundedline(x_scale,plot_data_cat_mag,e_cat_mag,'cmap',[0,0,1],'transparency', 0.2,'alpha','--');
    set(h_cat_mag,'linewidth',2);
    hold on
    
    %Outline
    h_out_grad=boundedline(x_scale,plot_data_out_grad,e_out_grad,'cmap','g','transparency', 0.4,'alpha');
    set(h_out_grad,'linewidth',2);
    hold on
    h_out_mag=boundedline(x_scale,plot_data_out_mag,e_out_mag,'cmap','g','transparency', 0.2,'alpha','--');
    set(h_out_mag,'linewidth',2);
    hold on
    
    %Texture
    h_text_grad=boundedline(x_scale,plot_data_text_grad,e_text_grad,'cmap',[1,0.5,0],'transparency', 0.6,'alpha');
    set(h_text_grad,'linewidth',2);
    hold on
    h_text_mag=boundedline(x_scale,plot_data_text_mag,e_text_mag,'cmap',[1,0.5,0],'transparency', 0.2,'alpha','--');
    set(h_text_mag,'linewidth',2);
    hold on
    
    xlim([-0.1 xlimit]);
%     set(gca,'XTick',[-0.1:0.1:xlimit]);
    set(gca,'XTick',[0:0.2:xlimit]);
    ylim([-0.2 0.7]);
    yPos = 0;
    plot(get(gca,'xlim'), [yPos yPos],':','color',[0 0 0],'linewidth',1.5); %baseline
    hold on
    xPos=0;
    plot([xPos xPos],get(gca,'ylim'),'color',[0 0 0],'linewidth',1); %onset
    hold on
    
    %signif for grad
    Signif_cat_Grad=zeros(num_timepoints,1)-5; Signif_Out_Grad=zeros(num_timepoints,1)-5; Signif_Text_Grad=zeros(num_timepoints,1)-5;
    Signif_cat_Grad(find(signif_cat_grad))=0.48; Signif_Out_Grad(find(signif_out_grad))=0.60;  Signif_Text_Grad(find(signif_text_grad))=0.54;
    
    %same for mag
    Signif_cat_Mag=zeros(num_timepoints,1)-5; Signif_Out_Mag=zeros(num_timepoints,1)-5; Signif_Text_Mag=zeros(num_timepoints,1)-5;
    Signif_cat_Mag(find(signif_cat_mag))=0.45; Signif_Out_Mag(find(signif_out_mag))=0.57;  Signif_Text_Mag(find(signif_text_mag))=0.51;
    
    plot(linspace(-0.1,xlimit,num_timepoints),Signif_cat_Grad,'o','color','b','MarkerFaceColor','b','linewidth',1);% signif for categ grad
    hold on
    plot(linspace(-0.1,xlimit,num_timepoints),Signif_Out_Grad,'o','color','g','MarkerFaceColor','g','linewidth',1); % signif for outline grad
    hold on
    plot(linspace(-0.1,xlimit,num_timepoints),Signif_Text_Grad,'o','color',[1,0.5,0],'MarkerFaceColor',[1,0.5,0],'linewidth',1); % signif for texture grad
    hold on
    
    plot(linspace(-0.1,xlimit,num_timepoints),Signif_cat_Mag,'o','color','b','linewidth',1);% signif for categ mag
    hold on
    plot(linspace(-0.1,xlimit,num_timepoints),Signif_Out_Mag,'o','color','g','linewidth',1); % signif for outline mag
    hold on
    plot(linspace(-0.1,xlimit,num_timepoints),Signif_Text_Mag,'o','color',[1,0.5,0],'linewidth',1); % signif for texture mag
    
    xlabel('time')
    ylabel('beta estimate')
%     title('GLM Cat+Out+Text')
    grid off
    axis square

end