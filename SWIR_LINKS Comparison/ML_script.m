%Refresh screen (clear all information and graphs displayed)
close all
clear variables
clc
% **********************************************************
%load data
%M7H16 is hyperspectral time series and WL is  weightloss data
%cd('C:\...') %M%addpath(genpath('D:\Ronan_paper\Final version'))
load M7H16
Wts=load('WL.txt');

%define wavelength & time variables
wvgood=(880:7:1720);
Time=[0,0,0,1,1,1,3,3,3,5,5,5,7,7,7,24,24,24,30,30,30,17*24,17*24,17*24,18*24,18*24,18*24]';

%Plot weight loss over time to 30 hours
%sample replicate 1 (Figure 5 in paper)
figure,plot(Time(1:3:21),Wts(1:3:21),'o')
title('Weight loss during first 30 hours','FontSize',16);
xlabel('Time (hours)','FontSize',16);
ylabel('Weight loss (%)','FontSize',16);
hold on
%sample replicate 2
plot(Time(2:3:21),Wts(2:3:21),'o')
%sample replicate 3
plot(Time(3:3:21),Wts(3:3:21),'o')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Removing unwanted data - section 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%4.1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Generate image of acquisition and spectra (Figure 6 in paper)
for i=1:5
figure(1),subplot(1,5,i),imshow(M7H16(i*3-2).file(:,:,50),[])
title(sprintf('%s, %dnm',M7H16(i*3-2).M7H16,wvgood(1,50)),'FontSize',16);
end
%Figure 8 in paper
for i=1:5
%plot spectra of sample region (red) and background regions (blue)
figure(2),subplot(1,5,i),plot(squeeze(wvgood(1,:)),squeeze(unfold(M7H16(i*3-2).file(71:94,59:73,:))),'r')
xlim([850 1750])
title(M7H16(i*3-2).M7H16,'FontSize',16);
hold on
%bg1
plot(squeeze(wvgood(1,:)),squeeze(unfold(M7H16(i*3-2).file(64:91,5:23,:))),'b')
%bg2
plot(squeeze(wvgood(1,:)),squeeze(unfold(M7H16(i*3-2).file(160:170,56:81,:))),'b')
%bg3
plot(squeeze(wvgood(1,:)),squeeze(unfold(M7H16(i*3-2).file(1:11,49:56,:))),'b')
end
%create spectral subset
for i=1:size(M7H16,2)
M7H16(i).file=M7H16(i).file(:,:,13:115);
end
wvgood=wvgood(13:115);
%plot cut spectra of sample (red) and background (blue) regions  (Fig.8b)
for i=1:5
figure(3),subplot(1,5,i),plot(squeeze(wvgood(1,:)),squeeze(unfold(M7H16(i*3-2).file(71:94,59:73,:))),'r')
xlim([850 1750])
title(M7H16(i*3-2).M7H16,'FontSize',16);
hold on
%bg1
plot(squeeze(wvgood(1,:)),squeeze(unfold(M7H16(i*3-2).file(64:91,5:23,:))),'b')
%bg2
plot(squeeze(wvgood(1,:)),squeeze(unfold(M7H16(i*3-2).file(160:170,56:81,:))),'b')
%bg3
plot(squeeze(wvgood(1,:)),squeeze(unfold(M7H16(i*3-2).file(1:11,49:56,:))),'b')
end
%%

%4.2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Selecting threshold
%selecting pixel spectra from foreground(FG)and background (BG)
FG_all=[];
BG_all=[];

for i=1:size(M7H16,2)
FG=squeeze(unfold(M7H16(i).file(71:94,59:73,:)));

BG=[squeeze(unfold(M7H16(i).file(64:91,5:23,:)));squeeze(unfold(M7H16(i).file(169:179,56:81,:)));squeeze(unfold(M7H16(i).file(49:56,1:11,:)))];

FG_all=[FG_all;FG];
BG_all=[BG_all;BG];
end

for i=1:103
    for j=1:103

    Ratio_FG=FG_all(:,i)./FG_all(:,j);
    Ratio_BG=BG_all(:,i)./BG_all(:,j);

    %calculate the distance between FG & BG ratios

    Ratio_FG_mean=mean(Ratio_FG);
    Ratio_BG_mean=mean(Ratio_BG);

    Ratio_FG_sd=std(Ratio_FG);
    Ratio_BG_sd=std(Ratio_BG);

        if Ratio_FG_mean>Ratio_BG_mean

        Dist(i,j)=(Ratio_FG_mean-2*Ratio_FG_sd)-(Ratio_BG_mean+2*Ratio_BG_sd);
        Thresh(i,j)=0.5*((Ratio_FG_mean-2*Ratio_FG_sd)+(Ratio_BG_mean+2*Ratio_BG_sd));
        else

        Dist(i,j)=(Ratio_BG_mean-2*Ratio_BG_sd)-(Ratio_FG_mean+2*Ratio_FG_sd);
        Thresh(i,j)=0.5*((Ratio_BG_mean-2*Ratio_BG_sd)+(Ratio_FG_mean+2*Ratio_FG_sd));
        end
    end
end

%Figure 10 in paper
figure,imagesc(Dist);
title(sprintf('Difference between Sample and Background Wavelength Ratio Ranges ','FontSize',16));
xlabel('Wavelength Slice');
ylabel('Wavelength Slice');

Max_D=(max(Dist(:)));

[a,b]=find(Dist>floor(1000*Max_D)/1000);

FG_optima=FG_all(:,a)./FG_all(:,b);
BG_optima=BG_all(:,a)./BG_all(:,b);

%histogram (Figure 11 in paper)
hist(FG_optima,100)
set(gca,'FontSize',16)
hold on
hist(BG_optima,100)
bar(mean(FG_optima),2000,0.001,'EdgeColor','r')
bar(mean(BG_optima),2000,0.001,'EdgeColor','r')
bar(mean(FG_optima)-2*std(FG_optima),2000,0.001,'EdgeColor','g')
bar(mean(BG_optima)+2*std(BG_optima),2000,0.001,'EdgeColor','g')
bar(0.5*(mean(FG_optima)-2*std(FG_optima)+mean(BG_optima)+2*std(BG_optima)),2000,0.001,'EdgeColor','b')

%Apply mask to each hypercube
for(i=1:size(M7H16,2))
    M7H16(i).ratio=[];
    M7H16(i).ratio=M7H16(i).file(:,:,a)./M7H16(i).file(:,:,b)>Thresh(a,b);
    [x,y,z]=size(M7H16(i).file);
    A=(zeros(x*y,z));
    M_u=reshape(M7H16(i).ratio,x*y,1);
    [Ind]=find(M_u>0);
    Hcube=reshape(M7H16(i).file,x*y,z);
    A(Ind,:)=Hcube(Ind,:);
    M7H16(i).MSKd=reshape(A,x,y,z);
end
%Figure 12
for i=1:5
   subplot(1,5,i),imshow(M7H16(i*3-2).MSKd(:,:,38));
   title(sprintf('%s,%dnm',M7H16(i*3-2).M7H16,wvgood(1,38)),'FontSize',16);
end
    
%4.3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fig 13: normal, dead, spike pixels
%remove background
[x,y,z]=size(M7H16(15).MSKd);
X_nz=nonzeros(M7H16(15).MSKd);
X_nz_re=reshape(X_nz,size(X_nz,1)/z,z);

%select a spectrum
test_spec=X_nz_re(500,:);

%make a dead pixel at wavel 100
test_spec_0=test_spec;
test_spec_0(100)=0;

%make a spike pixel at wavel 100
test_spec_1=test_spec;
test_spec_1(100)=1;

%Normal, dead pixel and spike figures
figure(15),subplot(1,3,1),plot(wvgood,test_spec')
title('Normal Spectrum','FontSize',16);
xlabel('Wavelength (nm)','FontSize',16);
ylabel('Reflectance','FontSize',16);
ylim([0 1]);
figure(15),subplot(1,3,2),plot(wvgood,test_spec_0')
title('Spectrum with Dead Pixel','FontSize',16);
xlabel('Wavelength (nm)','FontSize',16);
ylabel('Reflectance','FontSize',16);
ylim([0 1]);
figure(15),subplot(1,3,3),plot(wvgood,test_spec_1')
title('Spectrum with Spike','FontSize',16);
xlabel('Wavelength (nm)','FontSize',16);
ylabel('Reflectance','FontSize',16);
ylim([0 1]);

%calculate difference to identify spikes
diff=test_spec(2:end)-test_spec(1:end-1);
diff0=test_spec_0(2:end)-test_spec_0(1:end-1);
diff1=test_spec_1(2:end)-test_spec_1(1:end-1);

%figure
figure(16),subplot(1,3,1),plot(wvgood(2:end),diff),ylim([-0.6 0.6])
xlabel('Wavelength (nm)','FontSize',16);
ylabel('Difference','FontSize',16);
figure(16),subplot(1,3,2),plot(wvgood(2:end),diff0),ylim([-0.6 0.6])
xlabel('Wavelength (nm)','FontSize',16);
ylabel('Difference','FontSize',16);
figure(16),subplot(1,3,3),plot(wvgood(2:end),diff1),ylim([-0.6 0.6])
xlabel('Wavelength (nm)','FontSize',16);
ylabel('Difference','FontSize',16);

%Identify spikes in hypercube
Ad=M7H16(1).MSKd(:,:,1:z-1);

for(i=1:size(M7H16,2))
[x,y,z]=size(M7H16(i).MSKd);
Ad = M7H16(i).MSKd(:,:,1:z-1)-M7H16(i).MSKd(:,:,2:z);
Ad1 = reshape(Ad,x*y,z-1);
Ad2 = (Ad1)';
Ad3 = std(Ad2);
M7H16(i).Ad_SD=reshape((Ad3),x,y);
end

%Histogram for selecting threshold level (Figure 14 in paper)
for i=1:5
   figure(15),subplot(1,5,i),hist(M7H16(i*3-2).Ad_SD,100,'b')
   set(gca,'xlim',[0 0.020],'ylim',[0 55])
   title(M7H16(i*3-2).M7H16,'FontSize',16);
   hold on   
   xlabel('Standard Deviation (SD)');
   ylabel('Frequency');
end

%add threshold line
for i=1:5
   figure(15),subplot(1,5,i),bar(0.0135,40,0.0001,'EdgeColor','r')
end

%Apply threshold
for i=1:size(M7H16,2)
    M7H16(i).MS=(M7H16(i).Ad_SD>0.0135);
end

%Generate images of standard deviation of difference spectra (Figure 15)
for i=1:5
figure (16)
subplot(1, 5, i);
imshow(M7H16(i*3-2).Ad_SD,[0 0.02]);
title(sprintf('(a) SD of difference\n spectrum %s',M7H16(i*3-2).M7H16),'FontSize',16);
end

%View spike mask from standard deviation difference spectra threshold
for i=1:5
figure (17)
subplot(1, 5, i);
imshow(M7H16(i*3-2).MS,[0 0.02]);
title(sprintf('(b) Thresholded\n (a) %s',M7H16(i*3-2).M7H16),'FontSize',16);
end

for i=1:size(M7H16,2)
M7H16(i).unfMSKd=unfold(M7H16(i).MSKd);
end

%plot spike pixel spectra 
for i=1:5
    figure(18)
    subplot(1,5,i);
    plot(wvgood,M7H16(i*3-2).unfMSKd(find(M7H16(i*3-2).MS),:))
    title(sprintf('(c) Spike and dead\npixel spectra, %s',M7H16(i*3-2).M7H16),'FontSize',16);
    xlabel('Wavelength (nm)');
    ylabel('Reflectance');
end

%to remove the spikes from the original hypercube, use the complement of
%the of the spike mask, i.e. mask the pixels with difference matrix pixels 
%with high standard deviation 
for(i=1:size(M7H16,2))
    for(j=1:z)
    M7H16(i).A_MS(:,:,j)=M7H16(i).MSKd(:,:,j).*(M7H16(i).MS<1);
    end
end

for i=1:5
    figure(19)
    subplot(1,5,i),imshow(M7H16(i*3-2).A_MS(:,:,38),[])
    title(sprintf('(d) %s\nSpike Removed',M7H16(i*3-2).M7H16),'FontSize',16)
end

%4.4 spectral pretreatments

for i=1:size(M7H16,2)

    [x,y,z]=size(M7H16(i).A_MS);
    
    %unfold masked image
    M7H16(i).uf_A_MS=reshape(M7H16(i).A_MS,x*y,z);
    
    %find nonzero indices
    M7H16(i).row=find(M7H16(i).uf_A_MS(:,1)>0); 
    M7H16(i).row_dim=size(M7H16(i).row,1);
end
%Figure 16
figure,subplot(1,2,1)
plot(wvgood,mean(M7H16(5).uf_A_MS(M7H16(5).row,:)));
title(sprintf('No Pre-treatment',i),'FontSize',16);
xlabel('Wavelength (nm)','FontSize',16);
ylabel('Reflectance','FontSize',16);
hold on
plot(wvgood,mean(M7H16(5).uf_A_MS(M7H16(5).row,:))+std(M7H16(5).uf_A_MS(M7H16(5).row,:)),'--r');
plot(wvgood,mean(M7H16(5).uf_A_MS(M7H16(5).row,:))-std(M7H16(5).uf_A_MS(M7H16(5).row,:)),'--r');

subplot(1,2,2)
plot(wvgood,mean(MSC(M7H16(5).uf_A_MS(M7H16(5).row,:),mean(M7H16(5).uf_A_MS(M7H16(5).row,:),1)')));
title(sprintf('MSC',i),'FontSize',16);
xlabel('Wavelength (nm)','FontSize',16);
ylabel('Reflectance','FontSize',16);
hold on
plot(wvgood,mean(MSC(M7H16(5).uf_A_MS(M7H16(5).row,:),mean(M7H16(5).uf_A_MS(M7H16(5).row,:),1)'))+std(MSC(M7H16(5).uf_A_MS(M7H16(5).row,:),mean(M7H16(5).uf_A_MS(M7H16(5).row,:),1)')),'--r');
plot(wvgood,mean(MSC(M7H16(5).uf_A_MS(M7H16(5).row,:),mean(M7H16(5).uf_A_MS(M7H16(5).row,:),1)'))-std(MSC(M7H16(5).uf_A_MS(M7H16(5).row,:),mean(M7H16(5).uf_A_MS(M7H16(5).row,:),1)')),'--r');
ylim([0.2,0.8])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PCA section 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%5.1.1 PCA on each hypercube individually



%find non-zero pixel indices
%project the hypercubes along each of their PC loadings
for i=1:size(M7H16,2)
    
    %find mean of nonzeros 
    M7H16(i).Am=mean(M7H16(i).uf_A_MS(M7H16(i).row,:));
    
    %do PCA on nonzero spectra
    [M7H16(i).Loadings,M7H16(i).A_Scores]=princomp(M7H16(i).uf_A_MS(M7H16(i).row,:));
    
    %mean centre each cube (only nonzero parts)
    M7H16(i).A_MS_MC=zeros(x*y,z);
    M7H16(i).A_MS_MC_NZ=M7H16(i).uf_A_MS(M7H16(i).row,:)-repmat(M7H16(i).Am,M7H16(i).row_dim,1); 
    M7H16(i).A_MS_MC(M7H16(i).row,:)= M7H16(i).A_MS_MC_NZ;
    
    %project mean centred along PC eigenvectors
    M7H16(i).A_Scores=reshape((M7H16(i).A_MS_MC)*M7H16(i).Loadings(:,1:10),x,y,10);
end

%histogram plot (Figure 18)
k=0;
for j=1:7
    for i=1:5
    figure(24)
    subplot(7, 5, i+k),hist(M7H16(i*3-2).A_Scores(:,:,j))
    set(gca,'xlim',[min(min(M7H16(i*3-2).A_Scores(:,:,j)-0.1)) max(max(M7H16(i*3-2).A_Scores(:,:,j)+0.1))],'ylim',[0 55])
    hold on
    title(sprintf('%s PC%d',M7H16(i*3-2).M7H16,j),'FontSize',5);
    bar(max(max(M7H16(i*3-2).A_Scores(:,:,j))),50,0.001,'EdgeColor','r')
    bar(min(min(M7H16(i*3-2).A_Scores(:,:,j))),50,0.001,'EdgeColor','r')
    end
    k=k+5;
end
%Figure 19
%Images of PC1
for i=1:5    
    figure(25)
    subplot(7, 5, i),imshow(M7H16(i*3-2).A_Scores(:,:,1),[-1.7 1.5])
    title(sprintf('%s PC%d',M7H16(i*3-2).M7H16,1),'FontSize',5);   
end

%Images of PC2
for i=1:5    
    
    subplot(7, 5, i+5),imshow(M7H16(i*3-2).A_Scores(:,:,2),[-0.9 1.25])
    title(sprintf('%s PC%d',M7H16(i*3-2).M7H16,2),'FontSize',5);   
end

%Images of PC3
for i=1:5    
  
    subplot(7, 5, i+10),imshow(M7H16(i*3-2).A_Scores(:,:,3),[-0.17 0.20])
    title(sprintf('%s PC%d',M7H16(i*3-2).M7H16,3),'FontSize',5);   
end

%Images of PC4
for i=1:5    
    
    subplot(7, 5, i+15),imshow(M7H16(i*3-2).A_Scores(:,:,4),[-0.15 0.12])
    title(sprintf('%s PC%d',M7H16(i*3-2).M7H16,4),'FontSize',5);   
end

%Images of PC5
for i=1:5    
    
    subplot(7, 5, i+20),imshow(M7H16(i*3-2).A_Scores(:,:,5),[-0.06 0.05])
    title(sprintf('%s PC%d',M7H16(i*3-2).M7H16,5),'FontSize',5);   
end

%Images of PC6
for i=1:5    
    
    subplot(7, 5, i+25),imshow(M7H16(i*3-2).A_Scores(:,:,6),[-0.06 0.05])
    title(sprintf('%s PC%d',M7H16(i*3-2).M7H16,6),'FontSize',5);   
end

%Images of PC7
for i=1:5    
   
    subplot(7, 5, i+30),imshow(M7H16(i*3-2).A_Scores(:,:,7),[-0.06 0.07])
    title(sprintf('%s PC%d',M7H16(i*3-2).M7H16,7),'FontSize',5);   
end

%Check correlation between each PC and weightloss (WL) 


%extract mean Scores for each h cube
Mean_scores_Ind=zeros(27,10);
for (i=1:size(M7H16,2))
    S_uf=unfold(M7H16(i).A_Scores);
    Mean_scores_Ind(i,:)=mean(S_uf(M7H16(i).row,:));
end

%calculate Pearson's correlation coeffecient between each mean
%score & WL
for j=1:10
    Corr_Ind_W(j)=corr(Wts,Mean_scores_Ind(:,j));
end

%5.1.2 CONCATENATED PCA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%concatenate unfolded, background removed hypercubes 
Masterimage=[];
for i=1:size(M7H16,2)
    Masterimage = [Masterimage; M7H16(i).uf_A_MS(M7H16(i).row,:)];
end

%apply PCA
[Master_Loadings,Master_Scores]=princomp(Masterimage);
Mean_Master=mean(Masterimage);

%project the individual hypercubes along the PC loadings defined by
%Masterimage 

for (i=1:size(M7H16,2))
    [x,y,z]=size(M7H16(i).A_MS);
    %mean centre each cube (only nonzero parts)
    M7H16(i).A_MS_Master_MC=zeros(x*y,z);
    M7H16(i).A_MS_Master_MC_NZ=M7H16(i).uf_A_MS(M7H16(i).row,:)-repmat(Mean_Master,M7H16(i).row_dim,1); 
    M7H16(i).A_MS_Master_MC(M7H16(i).row,:)= M7H16(i).A_MS_Master_MC_NZ;
    %project mean centred along PC eigenvectors
    M7H16(i).Master_Scores=reshape((M7H16(i).A_MS_Master_MC)*Master_Loadings(:,1:10),x,y,10);
end

%Check correlation between each PC and WL for non-concatenated and
%concatenated PCA

%extract mean Scores for each h cube
Mean_scores_Master=zeros(27,10);
for (i=1:size(M7H16,2))
    S_uf=unfold(M7H16(i).Master_Scores);
    Mean_scores_Master(i,:)=mean(S_uf(M7H16(i).row,:));
end

%use corr to calculate Pearson's correlation coeffecient between each mean
%score &  WL

for j=1:10
    Corr_Master_W(j)=corr(Wts,Mean_scores_Master(:,j));
end
  
%Concatenating mean spectra from each hcube

Mean_spec=[];
for (i=1:size(M7H16,2)) 
    %find mean of nonzeros 
    Mean_spec(i,:)= M7H16(i).Am; 
end

[L_mean,S_mean]=princomp(Mean_spec);

for j=1:10
    Corr_Mean_W(j)=corr(Wts,S_mean(:,j));
end


%Applying MSC to Master & mean PCA as an example of how to apply a pretreatment

[Master_MSC_Loadings,Master_MSC_Scores]=princomp(MSC(Masterimage,mean(Masterimage)'));
Mean_Master_MSC=mean(MSC(Masterimage,mean(Masterimage)'));

%project the individual hypercubes along the PC loadings defined by
%Masterimage 
for (i=1:size(M7H16,2))
    [x,y,z]=size(M7H16(i).A_MS);
    %mean centre each cube (only nonzero parts)
    M7H16(i).A_MS_Master_MSC_MC=zeros(x*y,z);
    M7H16(i).A_MS_Master_MSC_MC_NZ=MSC(M7H16(i).uf_A_MS(M7H16(i).row,:),mean(Masterimage)')-repmat(Mean_Master_MSC,M7H16(i).row_dim,1); 
    M7H16(i).A_MS_Master_MSC_MC(M7H16(i).row,:)= M7H16(i).A_MS_Master_MSC_MC_NZ;
    %project mean centred along PC eigenvectors
    M7H16(i).Master_MSC_Scores=reshape((M7H16(i).A_MS_Master_MSC_MC)*Master_MSC_Loadings(:,1:10),x,y,10);
end

%extract mean Scores for each hcube
Mean_scores_Master_MSC=zeros(27,10);
for (i=1:size(M7H16,2))
    S_uf=unfold(M7H16(i).Master_MSC_Scores);
    Mean_scores_Master_MSC(i,:)=mean(S_uf(M7H16(i).row,:));
end

%use corr to calculate Pearson's correlation coeffecient between each mean
%score & WL

for j=1:10
    Corr_Master_MSC_W(j)=corr(Wts,Mean_scores_Master_MSC(:,j));
end

%Using MSC pretreatment for concatenated mean spectra before PCA

[L_mean_MSC,S_mean_MSC]=princomp(MSC(Mean_spec,mean(Mean_spec)'));

for j=1:10
    Corr_Mean_MSC_W(j)=corr(Wts,S_mean_MSC(:,j));
end

%examine effect of MSC on PCA scores (Figure 21)
%Look at PCA concatenated score images for no pretreatment and for MSC

figure,
for i = 1  
  for j = 3:3:27
      subplot(2,9,j/3) 
      imshow((M7H16(j).Master_Scores(:,:,i)),[-0.5,2])
      title(sprintf('%s\nUntreated', M7H16(j).M7H16),'FontSize',16);
  end
end
for i = 1  
  for j = 3:3:27
      subplot(2,9,9+j/3)
     imshow((M7H16(j).Master_MSC_Scores(:,:,i)),[-0.5,0.25])
     title(sprintf('%s\nMSC',M7H16(j).M7H16),'FontSize',16);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PLS section 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%make PLS function to do cal, train & test 
%apply to mean spectra

%firstly split the data into calibration, training and test sets
Mean_spec_cal=Mean_spec(1:3:27,:);
Mean_spec_train=Mean_spec(2:3:27,:);
Mean_spec_test=Mean_spec(3:3:27,:);

Wt_cal=Wts(1:3:27);
Wt_train=Wts(2:3:27);           
Wt_test=Wts(3:3:27);

PLS_Raw_w=PLS_output(Mean_spec_cal,Wt_cal,Mean_spec_train,Wt_train,Mean_spec_test,Wt_test,8);
PLS_MSC_w=PLS_output(MSC(Mean_spec_cal,mean(Mean_spec_cal)'),Wt_cal,MSC(Mean_spec_train,mean(Mean_spec_cal)'),Wt_train,MSC(Mean_spec_test,mean(Mean_spec_cal)'),Wt_test,8);

%project to cubes and calculate RMSEPIm for raw data

RMSEPIM=zeros(27,7);
MeanPRed=zeros(27,7);
for i = 1:27 
    for j = 1:7
  M7H16(i).PLS_Pred(:,j)= M7H16(i).uf_A_MS*PLS_Raw_w.B(:,j);
  RMSEPIM(i,j)=(sum((M7H16(i).PLS_Pred(M7H16(i).row,j)-Wts(i)).^2)/M7H16(i).row_dim)^0.5;
    end
end

%Plot predictions on raw data for 1st 7 LVs (Figure 24)
figure,
for i = 1:7  
  for j = 3:3:27
      subplot(7,9,(i-1)*9+j/3)
    [x,y,z]=size(M7H16(j).A_MS);
    imshow(reshape(M7H16(j).PLS_Pred(:,i),x,y),[-2,6])
  end
end

%Make prediction maps for MSC pretreated data
for i = 1:27 
    for j = 1:7
  M7H16(i).PLS_Pred_MSC(:,j)= MSC(M7H16(i).uf_A_MS,mean(Mean_spec_cal)')*PLS_MSC_w.B(:,j);
    end
end

%Plot MSC 3 LV & Raw 3LV predictions together (Figure 22)
figure,
 i = 3;  
  for j = 3:3:27
    subplot(2,9,j/3)
    imshow(reshape(M7H16(j).PLS_Pred(:,i),x,y),[-2,6])
    title(sprintf('%s\nUntreated',M7H16(j).M7H16),'FontSize',16);
    subplot(2,9,9+j/3)
    imshow(reshape(M7H16(j).PLS_Pred_MSC(:,i),x,y),[-2,6])
    title(sprintf('%s\nMSC',M7H16(j).M7H16),'FontSize',16);
  end

%calculate pooled RMSEPIM for 1st 18 images (i.e. cal & train sets)(Fig 25)
for j = 1:7
  RMSEPIM_ALL(j)=(sum(RMSEPIM(1:18,j).^2)/18).^0.5;
end

%Plot pooled RMSEPIM 
figure,plot(RMSEPIM_ALL,'-o')
title('RMSEPIm','FontSize',16);
xlabel('Number of PLS LV','FontSize',16)
ylabel('Pooled RMSEPIM','FontSize',16)



%generating table 6
Table6(1,:)=PLS_Raw_w.stats;
Table6(2,:)=PLS_MSC_w.stats;

%generating table 5
Table5(:,1)=Corr_Ind_W;
Table5(:,2)=Corr_Master_W;
Table5(:,3)=Corr_Master_MSC_W;
Table5(:,4)=Corr_Mean_W;
Table5(:,5)=Corr_Mean_MSC_W;