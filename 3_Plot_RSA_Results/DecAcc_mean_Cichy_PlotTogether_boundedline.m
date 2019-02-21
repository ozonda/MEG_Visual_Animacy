 %% Load the files 
close all
clear all

% %% Load Magnetometers and Gradiometers for one experiment
load('Final_Nov2017_MeanDecAcc_grad_Exp2Oneback')
mean_DecAcc_grad=mean_DecAcc; 
signif_grad=signif_dec_acc; 
e_DecAcc_grad=std(mean_DecAcc_subs,0,1) %/sqrt(14);

load('Final_Nov2017_MeanDecAcc_mag_Exp2Oneback') 
mean_DecAcc_mag=mean_DecAcc; 
signif_mag=signif_dec_acc; 
e_DecAcc_mag=std(mean_DecAcc_subs,0,1) %/sqrt(15);

% %we set error to zero to avoid plotting it with boundedline
% error=zeros(1,71);

%Plot the time course of the mean decoding accuracy
figure
set(0,'DefaultAxesFontSize', 20)
set(gcf,'color','w')
x_scale=linspace(-0.1,0.6,71);

h_cat=boundedline(x_scale,mean_DecAcc_grad,e_DecAcc_grad,'cmap',[0,0,1],'transparency', 0.4,'alpha');
set(h_cat,'linewidth',2);
hold on
h_cat=boundedline(x_scale,mean_DecAcc_mag,e_DecAcc_mag,'cmap',[0,0,1],'transparency', 0.3,'alpha','--');
set(h_cat,'linewidth',2);
hold on
xlim([-0.1 0.5]) ;
ylim([0.4 1]);
set(gca,'XTick',[-0.1:0.1:0.6]);
set(gca,'YTick',[0.5:0.1:1]);
set(gca,'Linewidth',2);

yPos = 0.5;
plot(get(gca,'xlim'), [yPos yPos],':','color',[0 0 0],'linewidth',1.5); %baseline
hold on
xPos=0;
plot([xPos xPos],get(gca,'ylim'),'color',[0 0 0],'linewidth',1); %onset
hold on
xPos_peak=1.3; %the peak
plot([xPos_peak xPos_peak],get(gca,'ylim'),'color','r','linewidth',1); %onset

Signif_grad=zeros(71,1)-5; 
Signif_grad(find(signif_grad))=0.9; 
plot(linspace(-0.1,0.6,71),Signif_grad,'o','color','b','MarkerFaceColor','b','linewidth',1);% signif for gradiometers
hold on

Signif_mag=zeros(71,1)-5; 
Signif_mag(find(signif_mag))=0.87; 
plot(linspace(-0.1,0.6,71),Signif_mag,'o','color','b','linewidth',1);% signif for magnetometers

xlabel('time')
ylabel('decoding accuracy')
% title('Mean pairwise decoding accuracy')
grid off
axis square