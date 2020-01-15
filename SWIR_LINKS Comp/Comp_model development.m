%Load data into workspace
load('SWIR_spectra.mat');
load('LINKS_spectra.mat');
load('HYPO.mat');

% Assign variables to each matrix
X=SWIR;
Z=LINKS;
Y=Hypo;

% Assign variables to each row an column of spectra matrices
[n,p]=size (X);
[q,r]=size (Z);
%% ******Display Data ***************

%Plot 3D Diagram of hypo number x wavelength x spectra
[dummy,h] = sort(Y);
figure (1),
%%oldorder = get(gcf,'DefaultAxesColorOrder');
set(gcf,'DefaultAxesColorOrder',jet(33));
subplot (1,2,1), plot3(repmat(950:5.4:2500,33,1)',repmat(Y(h),1,288)',SWIR(h,:)','linewidth',2);set(gca, 'Fontsize',28),
%%set(gcf,'DefaultAxesColorOrder',oldorder);
xlabel('Wavelength','FontWeight','bold'); ylabel('Hypo number','FontWeight','bold'); zlabel ('Mean reflectance spectra','FontWeight','bold');axis('tight');
title ('SWIR','FontWeight','bold'),
grid on
ax = gca
ax.GridAlpha = 0.35;  % Make grid lines less transparent.
ax.GridColor = [0,0,0]; % Dark Green.
oldorder = get(gcf,'DefaultAxesColorOrder');
set(gcf,'DefaultAxesColorOrder',jet(33));
subplot (1,2,2), plot3(repmat(401:1000,33,1)',repmat(Y(h),1,600)',LINKS(h,:)','linewidth',2);set(gca, 'Fontsize',28),
xlabel('Wavelength','FontWeight','bold'); ylabel('Hypo number','FontWeight','bold'); zlabel ('Mean reflectance spectra','FontWeight','bold');axis('tight');
title ('LINKS','FontWeight','bold');
grid on
ax = gca
ax.GridAlpha = 0.35;  % Make grid lines less transparent.
ax.GridColor = [0,0,0]; % Dark Green.
%%  Principal Component Analysis (SWIR data X)
% Method 1: EVD
%Taking the eigen decomposition of X transpose times X
C_SWIR=(X'*X);
% Assign the output of the eig decomposition function to matrices W and
% Lambda, which are equal in size
[W_SWIR, Lambda_SWIR] = eig (C_SWIR);
%W represents the loadings 
% Lambda represents the principal components 
% Reorder the  Lambda matrix so that the largest value appears first instead
% of last

W_SWIR= W_SWIR (:,end:-1:1);
lv_SWIR = flip (diag (Lambda_SWIR));

W_r_SWIR = W_SWIR(:, 1:2);
b_SWIR = X*W_r_SWIR;
size (b_SWIR);
figure (2); plot (b_SWIR(:,1), b_SWIR(:,2), 'o');
xlabel('PC 1'); ylabel ('PC2');

%Method 2 : SVD
[U_SWIR, Sigma_SWIR, V_SWIR] = svd(X); 
% Show that V and W are identical or very very close
V_r_SWIR = V_SWIR(:, 1:2);
scores_SWIR = X*V_r_SWIR;
size (a_SWIR);
figure (3); plot (scores_SWIR(:,1), scores_SWIR(:,2), 'o');
xlabel('PC 1'); ylabel ('PC2');


size (Sigma_SWIR);
sv_SWIR = diag (Sigma_SWIR);

figure (4); stairs (cumsum(sv_SWIR)/sum (sv_SWIR));
xlabel('number of components'); ylabel ('% Variance');
figure (5); bar (sv_SWIR);
xlabel('number of components'); ylabel ('% sigma values');
xlim([1 10]);
figure (6); bar(lv_SWIR);
xlabel('number of components'); ylabel ('% Lambda values');
xlim([1 10]);
%% Prinicipal componet analysis (LinkSquare data, Z)
% Method 1: EVD
%Taking the eigen decomposition of X transpose times X
C_LINKS=(Z'*Z);
% Assign the output of the eig decomposition function to matrices W and
% Lambda, which are equal in size
[W_LINKS, Lambda_LINKS] = eig (C_LINKS);
%W represents the loadings 
% Lambda represents the principal components 
% Reorder the  Lambda matrix so that the largest value appears first instead
% of last

W_LINKS= W_LINKS (:,end:-1:1);
lv_LINKS = flip (diag (Lambda_LINKS));

%Method 2 : SVD
[U_LINKS, Sigma_LINKS, V_LINKS] = svd(Z); 
% Show that V and W are identical or very very close
V_r_LINKS = V_LINKS(:, 1:2);
a_LINKS = X*V_r_LINKS;
size (a_LINKS);
figure (2); plot (a_LINKS(:,1), a_LINKS(:,2), 'o');
xlabel('PC 1'); ylabel ('PC2');

W_r_LINKS = W_LINKS(:, 1:2);
b_LINKS = X*W_r_LINKS;
size (b_SWIR);
figure (3); plot (b_LINKS(:,1), b_LINKS(:,2), 'o');
xlabel('PC 1'); ylabel ('PC2');

size (Sigma_LINKS);
sv_LINKS = diag (Sigma_LINKS);

figure (7); stairs (cumsum(sv_LINKS)/sum (sv_LINKS));
xlabel('number of components'); ylabel ('% Variance');
figure (8); bar (sv_LINKS);
xlabel('number of components'); ylabel ('% sigma values');
xlim([1 10]);
figure (9); bar(lv_LINKS);
xlabel('number of components'); ylabel ('% Lambda values');
xlim([1 10]);
%% Principal component regeression
%SWIR
PCA_loadings = V (:, 1:33);
PCA_scores = X * PCA_loadings;
% Principal component regression - (test) for the calibration model 
PCR_beta = regress (Y - mean (Y), PCA_scores (:, 1:2));
PCR_beta = PCA_loadings (:, 1:20)* PCR_beta;
PCR_beta = [mean(Y) - mean(X)*PCR_beta; PCR_beta];
PCR_yfit = [ones(n,1) X]*PCR_beta;

plot(Y,PCR_yfit,'r^');
xlabel('Observed Response');ylabel('Fitted Response');
xlim ([11 20]); ylim ([11 20]);
('Fitted Response');
legend({'PCR with X Components'},'location','NW');

%% PCA scores and loadings analysis
figure, plot(PCA_scores(:,1),PCA_scores(:,2),'o'),
xlabel('PC1');ylabel('PC2');

%% ******* Partical least square regression ******** - (test) for the calibration 
%input = Z
%output = Hypo number 
[PLS_beta,XSCORES, XLOADINGS,YSCORES, YLOADINGS, WEIGHTS, XPCTVAR, YPCTVAR] = pls (X, Y, 5);
PLS_yfit = [X]* PLS_beta;

plot(Y,PCR_yfit,'r^',Y,PLS_yfit,'bo');
xlabel('Observed Response');ylabel('Fitted Response');
xlim ([11 20]); ylim ([11 20]);
('Fitted Response');
legend({'PCR with X Components' 'PLS with X Components'},'location','NW');
%% PCR - validation of model predictions (k-fold cross validation)
%for simplication purpose K=11 -> 
% Training set = 3 Observations
% Test Set = 30 Observations
%PCR SWIR
c =cvpartition (Y, 'KFold',10,'Stratify',false);


