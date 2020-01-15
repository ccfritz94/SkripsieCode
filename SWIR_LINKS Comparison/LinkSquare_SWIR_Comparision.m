%%Comparing LinkSquare to the SWIR
%%
%Grettin SWIR Data in the correct format
%*********************** load NIR, LinkSquare and Experimental data into workspace*****************************
load('SWIR_spectra.mat');
load('LINKS_spectra.mat');
load('HYPO.mat');
X=NIR;
Z=LINKS;
y=Hypo;

[n,p]=size (X);
[q,r]=size (Z);
%%
%Plot 3D Diagram of hypo number x wavelength x spectra
[dummy,h] = sort(y);
figure (1),
%%oldorder = get(gcf,'DefaultAxesColorOrder');
set(gcf,'DefaultAxesColorOrder',jet(33));
subplot (1,2,1), plot3(repmat(950:5.4:2500,33,1)',repmat(y(h),1,288)',NIR(h,:)','linewidth',2);set(gca, 'Fontsize',28),
%%set(gcf,'DefaultAxesColorOrder',oldorder);
xlabel('Wavelength','FontWeight','bold'); ylabel('Hypo number','FontWeight','bold'); zlabel ('Mean reflectance spectra','FontWeight','bold');axis('tight');
title ('SWIR','FontWeight','bold'),
grid on
ax = gca
ax.GridAlpha = 0.35;  % Make grid lines less transparent.
ax.GridColor = [0,0,0]; % Dark Green.
oldorder = get(gcf,'DefaultAxesColorOrder');
set(gcf,'DefaultAxesColorOrder',jet(33));
subplot (1,2,2), plot3(repmat(401:1000,33,1)',repmat(y(h),1,600)',LINKS(h,:)','linewidth',2);set(gca, 'Fontsize',28),
xlabel('Wavelength','FontWeight','bold'); ylabel('Hypo number','FontWeight','bold'); zlabel ('Mean reflectance spectra','FontWeight','bold');axis('tight');
title ('LINKS','FontWeight','bold');
grid on
ax = gca
ax.GridAlpha = 0.35;  % Make grid lines less transparent.
ax.GridColor = [0,0,0]; % Dark Green.
%%
print(gcf,'figure 1.bmp','-dbmp','-r600');
%% Score plots
rng(0)
%PLS SWIR
c1 = cvpartition(y,'KFold',3,'Stratify',false) 
[PLS_Xloadings_SWIR,PLS_Yloadings_SWIR,PLS_Xscore_SWIR,PLS_Yscore_SWIR,PLS_beta_SWIR,PLS_Var_SWIR,PLS_msep_SWIR] = plsregress(X,y,32,'CV',c1);
%PLS LINKS
[PLS_Xloadings_LINKS,PLS_Yloadings_LINKS,PLS_Xscore_LINKS,PLS_Yscore_LINKS,PLS_beta_LINKS,PLS_Var_LINKS,PLS_msep_LINKS] = plsregress(Z,y,32,'CV',c1);
%PCA SWIR
[PCA_Loadings_SWIR,PCA_Scores_SWIR,PCA_Var_SWIR] = pca(X,'Economy',false);
PCR_msep_SWIR= sum(crossval(@pcrsse,X,y,'kfold',10),1)/n;
%PCA LINKS
[PCA_Loadings_LINKS,PCA_Scores_LINKS,PCA_Var_LINKS] = pca(Z,'Economy',false);
PCR_msep_LINKS= sum(crossval(@pcrsse,Z,y,'kfold',10),1)/n;

PCRmsep = sum(crossval(@pcrsse,X,y,'KFold',10),1) / n;
%Display MSEP on graph
figure (2), subplot(2,2,1), plot(0:32, PLS_msep_SWIR(2,:),'b-o','MarkerFaceColor','b','Markersize',15,'linewidth',3),set(gca,'box','off');set(gca, 'Fontsize',30),
title('PLS'),xlabel('Components','FontWeight','bold'),ylabel({'SWIR' 'MSEP'},'FontWeight','bold'), xlim([0 21]),ylim([0 8]),
subplot(2,2,2),plot(0:10,PCR_msep_SWIR(1,:),'r-o','MarkerFaceColor','r','Markersize',15,'linewidth',3),set(gca,'box','off');set(gca, 'Fontsize',30),
title('PCA'),xlabel('Components','FontWeight','bold'),ylabel('MSEP','FontWeight','bold'),ylim([0 8]),
subplot(2,2,3), plot(0:32, PLS_msep_LINKS(2,:),'b-^','MarkerFaceColor','b','Markersize',15,'linewidth',3),set(gca,'box','off');set(gca, 'Fontsize',30),
xlabel('Components','FontWeight','bold'),ylabel({'LINKS' 'MSEP'},'FontWeight','bold'), xlim([0 21]),ylim([0 6]),
subplot(2,2,4), plot(0:10, PCR_msep_LINKS(1,:),'r-^','MarkerFaceColor','r','Markersize',15,'linewidth',3),set(gca,'box','off');set(gca, 'Fontsize',30),
xlabel('Components','FontWeight','bold'),ylabel('MSEP','FontWeight','bold'),ylim([0 6]),
%sgtitle('Cross validation'),
%%
print(gcf,'figure 2.bmp','-dbmp','-r600');
%% Optimum components
Nopt_SWIR_PLS=5,MSEP_SWIR_PLS=min(PLS_msep_SWIR(2,:)),
Nopt_LINKS_PLS=3,MSEP_LINKS_PLS=min(PLS_msep_LINKS(2,:)),
Nopt_SWIR_PCA=4,MSEP_SWIR_PCA=min(PCR_msep_SWIR(1,:)),
Nopt_LINKS_PCA=3,MSEP_LINKS_PCA=min(PCR_msep_LINKS(1,:)),
%% Display variance explained by each principal component for SWIR and LINKS
figure(4), subplot(1,2,1),bar(1:288, PCA_Var_SWIR(:,1));set(gca,'box','off');set(gca, 'Fontsize',24),
title('SWIR'),xlim([0 10]),xlabel('Components'),ylabel('Variance'),
subplot(1,2,2),bar(1:600, PCA_Var_LINKS(:,1));set(gca,'box','off');set(gca, 'Fontsize',24),
title('LINKS'),xlim([0 10]),xlabel('Components'),ylabel('Variance'),
%%
print(gcf,'figure 4.bmp','-dbmp','-r600');
%%
%MAPE Plots for SWIR
figure (3),subplot(2,2,1), plot(1:32,PLS_mape_SWIR(1,:),'b-o'),            %Top left figure
xlim([1 32]), set(gca,'box','off');
legend({'PLS'},'location','NW');
xlabel('Components'),ylabel({'SWIR' 'MAPE (%)'});
subplot(2,2,2), plot(1:32,PCR_mape_SWIR(1,:),'r-^'),                       %Top right figure
xlim([1 32]), set(gca,'box','off');
legend({'PCR'},'location','NW');
%hold on, plot(0:10,PCRmsep_SWIR2(1,:),'b-^'),
xlabel('Components'),ylabel('MAPE (%)'),
%MAPE Plots for SWIR LinkSquare
subplot(2,2,3), plot(1:32,PLS_mape_LINKS(1,:),'b-o'),                      %Bottom left figure
xlim([1 32]), set(gca,'box','off');
legend({'PLS'},'location','NW');
xlabel('Components'),ylabel({'LINKS' 'MAPE (%)'}),
subplot(2,2,4), plot(1:32,PCR_mape_LINKS(1,:),'r-^'),                       %Bottom right figure
xlim([1 32]), set(gca,'box','off');
legend({'PCR'},'location','NW','box off');
xlabel('Components'),ylabel('MAPE (%)'),
print(gcf,'Regression model performance based on MAPE','-dbmp','-r300');
%%
%Show regression PLS and PCA regression model comparing LINKS square and SWIR 

%PLS fit response for 2 compnent model.
%SWIR
c1 = cvpartition(y,'KFold',3,'Stratify',false) 
[PLS_Xloadings_SWIR_opt,PLS_Yloadings_SWIR_opt,PLS_Xs_SWIR_opt,PLS_Ys_SWIR_opt,PLS_beta_SWIR_opt] = plsregress(X,y,5,'CV',c1);
PLS_yfit_SWIR_opt = [ones(n,1) X]*PLS_beta_SWIR_opt;
%PLS_TSS_SWIR = PLSmsep_LINKS (2,6)*n;
%PLS_RSS_PLS_SWIR = sum((y-PLS_yfit_SWIR_opt).^2);
%PLS_rsquared_SWIR = 1 - PLS_RSS_PLS_SWIR./PLS_TSS_SWIR


%LINKS
[PLS_Xl_LINKS_opt,PLS_Yl_LINKS_opt,PLS_Xs_LINKS_opt,PLS_Ys_LINKS_opt,PLS_beta_LINKS_opt] = plsregress(Z,y,2,'CV',c1);
PLS_yfit_LINKS_opt = [ones(n,1) Z]*PLS_beta_LINKS_opt;

%PLS_TSS_LINKS = PLSmsep_LINKS (2,3)*n;
%PLS_RSS_PLS_LINKS = sum((y-PLS_yfit_LINKS_opt).^2);
%PLS_rsquared_LINKS = 1 - PLS_RSS_PLS_LINKS./PLS_TSS_LINKS
% 10% lines
x = linspace(0,20);
figure (5), subplot(1,2,1),plot(y,PLS_yfit_SWIR_opt,'bo','MarkerFaceColor','b','Markersize',15,'linewidth',3),
hold on
plot(y,PLS_yfit_LINKS_opt,'r^','MarkerFaceColor','r','Markersize',15,'linewidth',3),
plot(x , x, 'k', x, 1.1*x, 'k--','linewidth',2),
plot(x ,x*0.9 ,'k--','linewidth',2),set(gca,'box','off');
%title('PLS regression model');
set(gca, 'Fontsize',28),
xlabel('Observed Response','FontWeight','bold');xlim([10 20]);set(gca,'XTick',(10:2:20));
ylabel('Fitted Response','FontWeight','bold');ylim([10 20]);set(gca,'YTick',(10:2:20));
legend({'SWIR n_o_p_t = 5' 'LinkSquare n_o_p_t = 3' 'ideal' '10% interval'},'location','NW','box','off');
pbaspect([1 1 1]),

%PCA fit in response to n=XXX components
%SWIR
[PCA_Loadings_SWIR,PCA_Scores_SWIR,PCA_Var_SWIR] = pca(X,'Economy',false);
PCR_beta_SWIR = regress(y-mean(y), PCA_Scores_SWIR(:,1:6));

PCR_beta_SWIR = PCA_Loadings_SWIR(:,1:6)*PCR_beta_SWIR;
PCR_beta_SWIR = [mean(y) - mean(X)*PCR_beta_SWIR; PCR_beta_SWIR];
PCR_yfit_SWIR= [ones(n,1) X]*PCR_beta_SWIR;

%PCR_TSS_SWIR = PCRmsep_SWIR (1,3).*n;
%PCR_RSS_PLS_SWIR = sum((y-PCR_yfit_SWIR).^2);
%PCR_rsquared_SWIR = 1 - PCR_RSS_PLS_SWIR./PCR_TSS_SWIR 

%LINKS
[PCA_Loadings_LINKS,PCA_Scores_LINKS,PCA_Var_LINKS] = pca(Z,'Economy',false);
PCR_beta_LINKS = regress(y-mean(y), PCA_Scores_LINKS(:,1:2));

PCR_beta_LINKS = PCA_Loadings_LINKS(:,1:2)*PCR_beta_LINKS;
PCR_beta_LINKS = [mean(y) - mean(Z)*PCR_beta_LINKS; PCR_beta_LINKS];
PCR_yfit_LINKS= [ones(n,1) Z]*PCR_beta_LINKS;

%PCR_TSS_LINKS = PCRmsep_LINKS (1,3)*n;
%PCR_RSS_PLS_LINKS = sum((y-PCR_yfit_LINKS).^2);
%PCR_rsquared_LINKS = 1 - PCR_RSS_PLS_LINKS./PCR_TSS_LINKS 
% 10% lines
x = linspace(0,20);

subplot(1,2,2), plot(y,PCR_yfit_SWIR,'bo','MarkerFaceColor','b','Markersize',15,'linewidth',3);
hold on
plot(y,PCR_yfit_LINKS,'r^','MarkerFaceColor','r','Markersize',15,'linewidth',3);
plot(x , x, 'k','linewidth',2),
plot(x, 1.1*x, 'k--','linewidth',2),
plot(x ,x*0.9 ,'k--','linewidth',2);set(gca,'box','off');
%title('PCR regression model'),
set(gca, 'Fontsize',28),
xlabel('Observed Response','FontWeight','bold');xlim([10 20]);set(gca,'XTick',(10:2:20));
ylabel('Fitted Response','FontWeight','bold');ylim([10 20]);set(gca,'YTick',(10:2:20));
legend({'SWIR n_o_p_t = 6' 'LinkSquare n_o_p_t = 2' 'ideal' '10% interval'},'location','NW','box','off');
hold off
pbaspect([1 1 1]),


%%
print(gcf,'figure 5.bmp','-dbmp','-r300');
%%
look at the amout of variation (in x) shown each regression model
%SWIR
figure (6), plot(1:10,100*cumsum(pctVar_SWIR(1,:)),'b-o',1:10,100*cumsum(PCA_Var_SWIR(1:10))/sum(PCA_Var_SWIR(1:10)),'r-^');
title ('SWIR')
xlabel('Number of Principal Components');
ylabel('Percent Variance Explained in X');
legend({'PLSR' 'PCR'},'location','SE');

%LINKS
figure (7), plot(1:10,100*cumsum(pctVar_LINKS(1,:)),'b-o',1:10,100*cumsum(PCA_Var_LINKS(1:10))/sum(PCA_Var_LINKS(1:10)),'r-^');
title ('LINKS')
xlabel('Number of Principal Components');
ylabel('Percent Variance Explained in X');
legend({'PLSR' 'PCR'},'location','SE');


%%
[Xl_S,Yl_S,Xs_S,Ys_S,beta_S,pctVar_S,mse_S,stats_S] = plsregress(X,y,5);
figure (18),subplot (2,2,1),plot(950:5.4:2500,stats_S.W,'-','linewidth',2);set(gca,'box','off');set(gca, 'Fontsize',19),
xlabel('Wavelength (nm)', 'Fontsize', 24);
ylabel({'SWIR' 'PLS weight'}, 'Fontsize', 24);
legend({'1^s^t' '2^n^d' '3^r^d' '4^t^h' '5^t^h'},'location','NW','box','off');
%LINKS
[Xl_L,Yl_L,Xs_L,Ys_L,beta_L,pctVar_L,mse_L,stats_L] = plsregress(Z,y,3);
subplot (2,2,3),plot(401:1000,stats_L.W,'-','linewidth',2);set(gca,'box','off');set(gca, 'Fontsize',19),
xlabel('Wavelength (nm)', 'Fontsize', 24);
ylabel({'LINKS' 'PLS weight'}, 'Fontsize', 24);
legend({'1^s^t ' '2^n^d' '3^r^d'},'location','NW','box','off');
%Wavelegth contributions to PCA 
%SWIR
subplot (2,2,2),plot(950:5.4:2500,PCA_Loadings_SWIR(:,1:4),'-','linewidth',2);set(gca,'box','off');set(gca, 'Fontsize',19),
xlabel('Wavelength (nm)', 'Fontsize', 24);
ylabel('PCA loadings', 'Fontsize', 24);
legend({'1^s^t' '2^n^d' '3^r^d' '4^t^h'},'location','NW','box','off');
%LINKS
subplot (2,2,4),plot(401:1000,PCA_Loadings_LINKS(:,1:3),'-','linewidth',2);set(gca,'box','off');set(gca, 'Fontsize',19),
xlabel('Wavelength (nm)', 'Fontsize', 24);
ylabel('PCA loadings', 'Fontsize', 24);
legend({'1^s^t ' '2^n^d' '3^r^d'},'location','NW','box','off');
%%
print(gcf,'figure 18.bmp','-dbmp','-r300');
%%
%***********************Model Parsimony***********************************
%Wavelegth contributions to PLS 
%SWIR
[Xl_S,Yl_S,Xs_S,Ys_S,beta_S,pctVar_S,mse_S,stats_S] = plsregress(X,y,10);
figure (18),subplot (1,2,1),plot(950:5.4:2500,stats_S.W,'-','linewidth',2);set(gca,'box','off');set(gca, 'Fontsize',19),
xlabel('Wavelength (nm)', 'Fontsize', 24);
ylabel('PLS weight', 'Fontsize', 24);
legend({'1st Component' '2nd Component' '3rd Component' '4th Component' '5td Component' '6th Component' '7th Component' '8th Component' '9th Component' '10 th Component'},  ...
	'location','NW','box','off');
%LINKS
[Xl_L,Yl_L,Xs_L,Ys_L,beta_L,pctVar_L,mse_L,stats_L] = plsregress(Z,y,5);
subplot (1,2,2),plot(401:1000,stats_L.W,'-','linewidth',2);set(gca,'box','off');set(gca, 'Fontsize',19),
xlabel('Wavelength (nm)', 'Fontsize', 24);
ylabel('PLS weight', 'Fontsize', 24);
legend({'1st Component' '2nd Component' '3rd Component' '4th Component' '5th Component'},  ...
	'location','NW','box','off');
%%
print(gcf,'figure 18.bmp','-dbmp','-r300');
%%
%Wavelegth contributions to PCA 
%SWIR
figure (19),subplot (1,2,1),plot(950:5.4:2500,PCA_Loadings_SWIR(:,1:6),'-','linewidth',2);set(gca,'box','off');set(gca, 'Fontsize',19),
xlabel('Wavelength (nm)', 'Fontsize', 24);
ylabel('PCA loadings', 'Fontsize', 24);
legend({'1st Component' '2nd Component' '3rd Component' '4th Component' '5th Component' '6th component'},  ...
	'location','NW','box','off');
%LINKS
subplot (1,2,2),plot(401:1000,PCA_Loadings_LINKS(:,1:3),'-','linewidth',2);set(gca,'box','off');set(gca, 'Fontsize',19),
xlabel('Wavelength (nm)', 'Fontsize', 24);
ylabel('PCA loadings', 'Fontsize', 24);
legend({'1st Component' '2nd Component' '3rd Component'},'location','NW','box','off');
%%
print(gcf,'PCR wavelength weightage for SWIR and LINKS.bmp','-dbmp','-r300');
%%
