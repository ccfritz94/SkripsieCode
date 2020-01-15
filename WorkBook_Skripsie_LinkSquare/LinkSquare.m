%TITLE: SKRIPSIE (C.C Fritz) 
%VERSION: 02
%This program contains a main menu for analysis of LinkSquare data of
%   pulp sample. It fistly displays information about the array.
%  *********************************************************************
%% *********************************************************************
%Refresh screen (clear all information and graphs displayed)
close all
clear variables
clc
%Define wavelength variable
Y=(401:1:1000);

%***********************Preparing data for preprossing********************
%%load data into matlab workspace
S1 = readtable('S1.csv');
S2 = readtable('S2.csv');
S3 = readtable('S3.csv');
S4 = readtable('S4.csv');
S5 = readtable('S5.csv');
S6 = readtable('S6.csv');
S7 = readtable('S7.csv');
S8 = readtable('S8.csv');
S9 = readtable('S9.csv');
S10 = readtable('S10.csv');

S11 = readtable('S11.csv');
S12 = readtable('S12.csv');
S13 = readtable('S13.csv');
S14 = readtable('S14.csv');
S15 = readtable('S15.csv');
S16 = readtable('S16.csv');
S17 = readtable('S17.csv');
S18 = readtable('S18.csv');
S19 = readtable('S19.csv');
S20 = readtable('S20.csv');

S21 = readtable('S21.csv');
S22 = readtable('S22.csv');
S23 = readtable('S23.csv');
S24 = readtable('S24.csv');
S25 = readtable('S25.csv');
S26 = readtable('S26.csv');
S27 = readtable('S27.csv');
S28 = readtable('S28.csv');
S29 = readtable('S29.csv');
S30 = readtable('S30.csv');

S31 = readtable('S31.csv');
S32 = readtable('S32.csv');
S33 = readtable('S33.csv');


%***********************create structure for linksquare data***************
%covert data from table to matlab data structure
Links(1).name='S1'; Links(1).data=S1;Links(1).hypo=15.95;
Links(2).name='S2'; Links(2).data=S2;Links(2).hypo=16.80;
Links(3).name='S3'; Links(3).data=S3;Links(3).hypo=13.18;
Links(4).name='S4'; Links(4).data=S4;Links(4).hypo=16.08;
Links(5).name='S5'; Links(5).data=S5;Links(5).hypo=18.06;
Links(6).name='S6'; Links(6).data=S6;Links(6).hypo=12.35;
Links(7).name='S7'; Links(7).data=S7;Links(7).hypo=13.47;
Links(8).name='S8'; Links(8).data=S8;Links(8).hypo=16.17;
Links(9).name='S9'; Links(9).data=S9;Links(9).hypo=13.55;
Links(10).name='S10'; Links(10).data=S10;Links(10).hypo=15.96;

Links(11).name='S11'; Links(11).data=S11;Links(11).hypo=13.87;
Links(12).name='S12'; Links(12).data=S12;Links(12).hypo=13.20;
Links(13).name='S13'; Links(13).data=S13;Links(13).hypo=17.55;
Links(14).name='S14'; Links(14).data=S14;Links(14).hypo=14.14;
Links(15).name='S15'; Links(15).data=S15;Links(15).hypo=17.68;
Links(16).name='S16'; Links(16).data=S16;Links(16).hypo=14.92;
Links(17).name='S17'; Links(17).data=S17;Links(17).hypo=16.26;
Links(18).name='S18'; Links(18).data=S18;Links(18).hypo=13.67;
Links(19).name='S19'; Links(19).data=S19;Links(19).hypo=12.66;
Links(20).name='S20'; Links(20).data=S20;Links(20).hypo=13.36;

Links(21).name='S21'; Links(21).data=S21;Links(21).hypo=11.92;
Links(22).name='S22'; Links(22).data=S22;Links(22).hypo=16.61;
Links(23).name='S23'; Links(23).data=S23;Links(23).hypo=16.64; % Known hypo number=17.66 P-sample
Links(24).name='S24'; Links(24).data=S24;Links(24).hypo=16.13;
Links(24).name='S24'; Links(24).data=S24;Links(24).hypo=16.13;
Links(25).name='S25'; Links(25).data=S25;Links(25).hypo=13.15;
Links(26).name='S26'; Links(26).data=S26;Links(26).hypo=12.75;
Links(27).name='S27'; Links(27).data=S27;Links(27).hypo=12.78;
Links(28).name='S28'; Links(28).data=S28;Links(28).hypo=14.13;
Links(29).name='S29'; Links(29).data=S29;Links(29).hypo=13.74;
Links(30).name='S30'; Links(30).data=S30;Links(30).hypo=14.82;

%D-samples
Links(31).name='S31'; Links(31).data=S31;Links(31).hypo=16.65;
Links(32).name='S32'; Links(32).data=S32;Links(32).hypo=18.90;
Links(33).name='S33'; Links(33).data=S33;Links(33).hypo=19.80;
%% Data formating
for counter=1:33
    Links(counter).data(:,1)=[];                % Remove first two columns (col 1-> timestamp)
    
    %Convert table to matrix
    Links(counter).actdata=table2array(Links(counter).data)
   
end
%% *****************Preprocessing sample 1 *********************************
%************************Removing noise **********************************
%A=1:2:9; % LED rows
%B=2:2:10;% Bulb rows

%Smoothing of data- removing noise using a savitzky-golay digital filter 
%for data using LED data
Links(1).smooth_led = sgolayfilt(Links(1).actdata (1,:),4,51); 
%
%for data using Bulb data
Links(1).smooth_bulb= sgolayfilt(Links(1).actdata (2,:),2,35);

% Compare smoothing technique to actual data
%For LED data
figure(1),subplot(1,2,2), plot (Y,Links(1).smooth_led,'b','linewidth',2),
xlim ([400 1000]),ylim([ 5 40]),set(gca, 'Fontsize',24),
xlabel('Wavelength (nm)','Fontsize',28);ylabel('Intensity (A.U)','Fontsize',28);set(gca,'box','off');
subplot(1,2,1),plot (Y,Links(1).actdata (1,:),'b','linewidth',2),
xlim ([400 1000]),ylim([ 5 40]),set(gca, 'Fontsize',24),
xlabel('Wavelength (nm)','Fontsize',28);ylabel('Intensity (A.U)','Fontsize',28);set(gca,'box','off');
%title('Light intensity of first scan of sample 1 using LED source','FontSize',16)
hold off
%% Print picture of removal of noise of LED scoure data- figure 1-Light intensity of first scan of sample 1 using LED source
print(gcf,'figure 1.bmp','-dbmp','-r300')

%% For Bulb data
figure (2), subplot(1,2,1),plot (Y,Links(1).actdata (2,:),'b','linewidth',2), 
xlim ([400 1000]),set(gca, 'Fontsize',24),
xlabel('Wavelength (nm)','Fontsize',28);ylabel('Intensity (A.U)','Fontsize',28);set(gca,'box','off');
subplot(1,2,2),plot (Y,Links(1).smooth_bulb,'b','linewidth',2),
xlim ([400 1000]),set(gca, 'Fontsize',24),
xlabel('Wavelength (nm)','Fontsize',28);ylabel('Intensity (A.U)','Fontsize',28);set(gca,'box','off');
hold off
%% Print picture of removal of noise of Bulb scoure data- figure 2-Light intensity of first scan of sample 1 using Bulb source
print(gcf,'figure 2.bmp','-dbmp','-r300')
%% Removal of noise for all LED and Bulb scans
%for data collected using LED Data
for i=1:2:9
    Links(1).smooth_led2(i,:) = sgolayfilt(Links(1).actdata (i,:),4,51);
    
end 

%for data collected using bulb
for j=2:2:10
    Links(1).smooth_bulb2(j,:)= sgolayfilt(Links(1).actdata (j,:),2,35);
end 
%% Need to removed empty entries in matix
Links(1).smooth_led2( ~any(Links(1).smooth_led2,2), : ) = [];  %from LED data 
Links(1).smooth_bulb2( ~any(Links(1).smooth_bulb2,2), : ) = [];%from bulb data

%% Displaying reflectance spectra for sample 1- using LED source
figure(3), plot (Y,Links(1).smooth_led2 ,'b','linewidth',2),
xlim ([400 1000]),set(gca, 'Fontsize',24),
xlabel('Wavelength (nm)','Fontsize',28);ylabel('Intensity (A.U)','Fontsize',28);set(gca,'box','off');
%% Print picture of all scanse of sample 1- figure 3 -using the Led source 
print(gcf,'figure 3.bmp','-dbmp','-r300')
%% Displaying reflectance spectra for sample 1-using Bulb source
figure(4), plot (Y,Links(1).smooth_bulb2 ,'b','linewidth',2),
xlim ([400 1000]),set(gca, 'Fontsize',24),
xlabel('Wavelength (nm)','Fontsize',28);ylabel('Intensity (A.U)','Fontsize',28);set(gca,'box','off');

%% Print picture of all scanse of sample 1- figure 3 -using the Buld source 
print(gcf,'figure 4.bmp','-dbmp','-r300')

%% ***************************************************************************
%find the avergae of the NIR data 
% For white LED light
Links(1).mean_led = mean(Links(1).smooth_led2);
Links(1).std_led=std(Links(1).smooth_led2);

% For bulb
Links(1).mean_bulb = mean(Links(1).smooth_bulb2);
Links(1).std_bulb=std(Links(1).smooth_bulb2);
%% Displaying mean reflectance spectra for sample 1-using White LED light
% Plot smooth data
%Displaying reflectance spectra for sample 1-using White LED light
figure(5), subplot (1,2,1),plot (Y,Links(1).smooth_led2,'b','linewidth',2),
xlim ([400 1000]),ylim([0 65]),set(gca,'box','off');set(gca, 'Fontsize',24),
ylabel('Intensity (A.U)','Fontsize',28);xlabel('Wavelength (nm)','Fontsize',28);
hold off
subplot (1,2,2),plot(Y,Links(1).mean_led,'b','linewidth',2),
hold on,
plot (Y,Links(1).mean_led + Links(1).std_led,'--r', 'linewidth',2),
plot (Y,Links(1).mean_led - Links(1).std_led,'--r','linewidth',2),
xlim ([400 1000]),ylim ([0 65]),set(gca,'box','off');set(gca, 'Fontsize',24),
ylabel('Intensity (A.U)','Fontsize',28);xlabel('Wavelength (nm)','Fontsize',28);
legend({'average' 'standard deviation'}, 'location', 'NW', 'box','off');
hold off
%% Print picture of  the mean spectra of sample 1- figure 3 -using the LED source 
print(gcf,'figure 5.bmp','-dbmp','-r300');
%% Displaying mean reflectance spectra for sample 1-using Bulb
figure(6), subplot(1,2,1),plot (Y,Links(1).smooth_bulb2,'b','linewidth',2),
xlim ([400 1000]),ylim([0 140]),set(gca,'box','off');set(gca, 'Fontsize',24),
ylabel('Intensity (A.U)','Fontsize',28);xlabel('Wavelength (nm)','Fontsize',28);
hold off
subplot (1,2,2),plot(Y,Links(1).mean_bulb,'b','linewidth',2), 
hold on
plot (Y,Links(1).mean_bulb + Links(1).std_bulb ,'--r','linewidth',2),
plot (Y,Links(1).mean_bulb - Links(1).std_bulb,'--r','linewidth',2),
xlim ([400 1000]),ylim([0 140]),set(gca,'box','off');set(gca, 'Fontsize',24),
ylabel('Intensity (A.U)','Fontsize',28);xlabel('Wavelength (nm)','Fontsize',28);
hold off
%% Print picture of  the mean spectra of sample 1- figure 3 -using the Bulb source 
print(gcf,'figure 6.bmp','-dbmp','-r300');
%% *****************Preprocessing all samples ******************************
%************************Removing noise **********************************
%A=1:2:9; % LED rows
%B=2:2:10;% Bulb rowss

%Smoothing of data 
% smothing for all the rows 
%LED Data
for counter =1:33 % number of samples
    for i=1:2:9
        Links(counter).smooth_led(i,:) = sgolayfilt(Links(counter).actdata (i,:),4,51);
    end 

%for data using bulb
    for j=2:2:10
        Links(counter).smooth_bulb(j,:)= sgolayfilt(Links(counter).actdata (j,:),2,35);
	end 

%Detele zero rows
Links(counter).smooth_led( ~any(Links(counter).smooth_led,2), : ) = [];  %from LED data 
Links(counter).smooth_bulb( ~any(Links(counter).smooth_bulb,2), : ) = [];%from bulb data

%plot mean spectra for LED
Links(counter).mean_led= mean(Links(counter).smooth_led);
Links(counter).std_led= std(Links(counter).smooth_led);         % Standard Deviation
figure(7),plot (squeeze(Y(1,:)),Links(counter).mean_led,'b')
xlim ([400 1000])
title('Mean light intensity for all pulp samples using LED source','FontSize',16),xlabel('Wavelength (nm)');
ylabel('Intensity (A.U)');
hold on

%plot mean spectra for LED
Links(counter).mean_bulb= mean(Links(counter).smooth_bulb);
Links(counter).std_bulb= std(Links(counter).smooth_bulb);       % Standard Deviation
figure(8), plot (squeeze(Y(1,:)),Links(counter).mean_bulb,'b')
xlim ([400 1000])
title('Mean light intensity for all pulp samples using bulb source','FontSize',16),xlabel('Wavelength (nm)');
ylabel('Intensity (A.U)','FontSize',16);
hold on

XS_m_led (counter,:)=Links(counter).mean_led(1,:);
XS_m_bulb (counter,:)=Links(counter).mean_bulb(1,:);
end 
%% ********************* Developing Regression models***********************
%*************************************************************************
%Formating data for regression model development
for i= 1:33
    Hypo(i,1)=Links(i).hypo;
end
LINKS =XS_m_led;
X = LINKS; % Mean spectra
y = Hypo;% Hypo numbers
[n,p] = size(X);
%Plot 3D Diagram of hypo number x wavelength x spectra
[dummy,h] = sort(y);
figure (9),
oldorder = get(gcf,'DefaultAxesColorOrder');
set(gcf,'DefaultAxesColorOrder',jet(33));
plot3(repmat(401:1000,33,1)',repmat(y(h),1,600)',LINKS(h,:),'linewidth',2);
set(gcf,'DefaultAxesColorOrder',oldorder);set(gca, 'Fontsize',24);set(gca,'box','off'),
xlabel('Wavelength'); ylabel('Hypo number'); zlabel ('Reflectance spectra');axis('tight');
grid on
%% Print picture of  the mean spectra of sample 1- figure 3 -using the Bulb source 
print(gcf,'figure 9.bmp','-dbmp','-r300');
