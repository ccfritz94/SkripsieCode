%TITLE: SKRIPSIE (C.C Fritz) 
%VERSION: 01
%This program contains a main menu for analysis of a single 3D array of
%   pulp sample. It fistly displays information about the array.
%  *********************************************************************
%% Clear command window and all variables stored
close all
clear variables
clc
%define wavelength variable
Y=(950:5.4:2500);
%*********************** load hyperspectral data into workspace *****************************
load('S1_1.mat');
load('S2_1.mat');
load('S3_1.mat');
load('S4_1.mat');
load('S5_1.mat');
load('S6_1.mat');
load('S7_1.mat');
load('S8_1.mat');
load('S9_1.mat');
load('S10_1.mat');

load('S11_1.mat');
load('S12_1.mat');
load('S13_1.mat');
load('S14_1.mat');
load('S15_1.mat');
load('S16_1.mat');
load('S17_1.mat');
load('S18_1.mat');
load('S19_1.mat');
load('S20_1.mat');

load('S21_1.mat');
load('S22_1.mat');
load('S23_1.mat');
load('S24_1.mat');
load('S25_1.mat');
load('S26_1.mat');
load('S27_1.mat');
load('S28_1.mat');
load('S29_1.mat');
load('S30_1.mat');

load('S31_1.mat');
load('S32_1.mat');
load('S33_1.mat');

% Change S5_hypo to S23
S23_1=Sample_5_Hypo_number_no_17_66_1;
%%
%***********************create structure for hyperspectral data*****************************
Data(1).name='sample 1'; Data(1).cube=S1_1;Data(1).hypo=15.95;
Data(2).name='sample 2'; Data(2).cube=S2_1;Data(2).hypo=16.80;
Data(3).name='sample 3'; Data(3).cube=S3_1;Data(3).hypo=13.18;
Data(4).name='sample 4'; Data(4).cube=S4_1;Data(4).hypo=16.08;
Data(5).name='sample 5'; Data(5).cube=S5_1;Data(5).hypo=18.06;
Data(6).name='sample 6'; Data(6).cube=S6_1;Data(6).hypo=12.35;
Data(7).name='sample 7'; Data(7).cube=S7_1;Data(7).hypo=13.47;
Data(8).name='sample 8'; Data(8).cube=S8_1;Data(8).hypo=16.17;
Data(9).name='sample 9'; Data(9).cube=S9_1;Data(9).hypo=13.55;
Data(10).name='sample 10'; Data(10).cube=S10_1;Data(10).hypo=15.96;

Data(11).name='sample 11'; Data(11).cube=S11_1;Data(11).hypo=13.87;
Data(12).name='sample 12'; Data(12).cube=S12_1;Data(12).hypo=13.20;
Data(13).name='sample 13'; Data(13).cube=S13_1;Data(13).hypo=17.55;
Data(14).name='sample 14'; Data(14).cube=S14_1;Data(14).hypo=14.14;
Data(15).name='sample 15'; Data(15).cube=S15_1;Data(15).hypo=17.68;
Data(16).name='sample 16'; Data(16).cube=S16_1;Data(16).hypo=14.92;
Data(17).name='sample 17'; Data(17).cube=S17_1;Data(17).hypo=16.26;
Data(18).name='sample 18'; Data(18).cube=S18_1;Data(18).hypo=13.67;
Data(19).name='sample 19'; Data(19).cube=S19_1;Data(19).hypo=12.66;
Data(20).name='sample 20'; Data(20).cube=S20_1;Data(20).hypo=13.36;

Data(21).name='sample 21'; Data(21).cube=S21_1;Data(21).hypo=11.92;
Data(22).name='sample 22'; Data(22).cube=S22_1;Data(22).hypo=16.61;
Data(23).name='sample 23'; Data(23).cube=S23_1;Data(23).hypo=16.64;
Data(24).name='sample 24'; Data(24).cube=S24_1;Data(24).hypo=16.13;
Data(25).name='sample 25'; Data(25).cube=S25_1;Data(25).hypo=13.15;
Data(26).name='sample 26'; Data(26).cube=S26_1;Data(26).hypo=12.75;
Data(27).name='sample 27'; Data(27).cube=S27_1;Data(27).hypo=12.78;
Data(28).name='sample 28'; Data(28).cube=S28_1;Data(28).hypo=14.13;
Data(29).name='sample 29'; Data(29).cube=S29_1;Data(29).hypo=13.74;
Data(30).name='sample 30'; Data(30).cube=S30_1;Data(30).hypo=14.82;

Data(31).name='sample 31'; Data(31).cube=S31_1;Data(31).hypo=16.65;
Data(32).name='sample 32'; Data(32).cube=S32_1;Data(32).hypo=18.90;
Data(33).name='sample 33'; Data(33).cube=S33_1;Data(33).hypo=19.80;

%% *****************Preprocessing sample 1***************************
% Creat filter based on reference reflection of 0
imtool(Data(1).cube(:,:,50)),
figure (1),
level=0;     % Determined threshold based on information in figure (1). The threshold is normalized to the range [0,1]
BW=imbinarize(Data(1).cube(:,:,50),level);  %Convert the image into a binary image using the threshold
B=imshowpair(Data(1).cube(:,:,50),BW,'montage');%Display orginal image next to binary image
img=im2double(BW); % Actual filter
%% Print orginal image next to binary image
print(gcf,'figure 1.bmp','-dbmp','-r300');
%% Filter on rows                
BG=squeeze(unfold (Data(1).cube(:,:,:).*img));          %Create background matrix
FG=squeeze(unfold (Data(1).cube(:,:,:)))-BG;          %Create sample matrix

%Remove zero rows
BG( ~any(BG,2), : ) = [];  %rows
FG( ~any(FG,2), : ) = [];  %rows

%Determining mean spectra
FG_shift = FG +1;
FG_mean_shift = mean(FG_shift);
FG_mean = FG_mean_shift-1;

%Determing std dev of spectra
FG_std=std(FG_shift);
%%
%Plot mean spectra
%Displaying reflectance spectra for sample 1
figure(2), subplot (1,2,1),plot (squeeze(Y(1,:)),FG,'b','linewidth',2),
xlim ([950 2500]),set(gca,'box','off');set(gca, 'Fontsize',24),
xlabel('Wavelength (nm)');ylabel('Reflectance spectra');

subplot (1,2,2),plot (squeeze(Y(1,:)),FG_mean,'b','linewidth',2),
hold on, 
plot(squeeze(Y(1,:)),FG_mean + FG_std, '--r','linewidth',2),
plot(squeeze(Y(1,:)),FG_mean - FG_std, '--r','linewidth',2),
xlim ([950 2500]),set(gca,'box','off');set(gca, 'Fontsize',24),
xlabel('Wavelength (nm)'),ylabel('Reflectance spectra');
legend({'average' 'standard deviation'}, 'location', 'NW','box','off');
hold off
%%
print(gcf,'figure 2.bmp','-dbmp','-r300');
%% *****************Preprocessing all samples***************************
for i=1:33
    % Create a filter
    Data(i).filter =im2double(imbinarize(Data(i).cube(:,:,50),level));
    %Apply filter
    Data(i).XB = squeeze(unfold(Data(i).cube(:,:,:).*Data(i).filter));
    Data(i).XS =squeeze(unfold(Data(i).cube(:,:,:)))-Data(i).XB;
    
    %Remove zero rows
    Data(i).XB( ~any(Data(i).XB,2), : ) = [];  %rows
    Data(i).XS( ~any(Data(i).XS,2), : ) = [];  %row
    
    %Determining mean spectra
    FG_shift = Data(i).XS +1;
    FG_mean_shift = mean(FG_shift);
    Data(i).XS_mean = FG_mean_shift-1;
    
    %Ploting mean spectra for all samples
    
    figure(3), plot (squeeze(Y(1,:)),Data(i).XS_mean,'b')
    xlim ([950 2500]),set(gca, 'Fontsize',24),
    xlabel('Wavelength (nm)','Fontsize',28);set(gca,'box','off');
    ylabel('Reflectance','Fontsize',28),
    hold on, 
end

%********************* Developing Regression models***********************
%*************************************************************************
%Formating data for regression model development
for i=1:33
    NIR (i,:)=Data(i).XS_mean(1,:);
end
%%
print(gcf,'figure 3.bmp','-dbmp','-r300');
%%

for i= 1:33
    Hypo(i,1)=Data(i).hypo;
end

X = NIR; % Mean spectra
y = Hypo;% Hypo numbers
[n,p] = size(X);
%Plot 3D Diagram of hypo number x wavelength x spectra
[dummy,h] = sort(y);
figure (4),
oldorder = get(gcf,'DefaultAxesColorOrder');
set(gcf,'DefaultAxesColorOrder',jet(33));
plot3(repmat(950:5.4:2500,33,1)',repmat(y(h),1,288)',NIR(h,:)','linewidth',2);set(gca, 'Fontsize',24);
set(gcf,'DefaultAxesColorOrder',oldorder);
xlabel('Wavelength'); ylabel('Hypo number'); zlabel ('Reflectance spectra');axis('tight');
grid on
%%
print(gcf,'figure 4.bmp','-dbmp','-r300');
