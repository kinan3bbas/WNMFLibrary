%% Data simulation
clear
clc
rng(1)

r = 4; % rank of the data matrix X
m = 500; n = 500; % size of the data matrix X (m x n)

Gtheo = rand(m,r); % simulation of theoretical matrix Gtheo
Ftheo = rand(r,n); % simulation of theoretical matrix Ftheo
Xtheo = Gtheo*Ftheo; % simulation of theoretical data matrix Xtheo

N = zeros(m,n); % simulating noise matrice N
X = Xtheo+N; % simulating data matrix X


%% Matrix initialisation
Ginit = rand(m,r);
Finit = rand(r,n);


%% NMF parameters
Iter_max = 1.e4;
rho_G = .01;
rho_F = .01;

%% Simulating matrices for weighted methods (missing entry example)
W = ones(m,n);
prop_missing = .1; % proportion of missing entries in X
idx_missing = randperm(m*n);
W(idx_missing(1:round(prop_missing*m*n))) = 0;
X = X.*W;

%% Running EM_WMU_NMF
Iter_max_E_step=10;
% WMU_NMF
fprintf('\nRunning NMF Lee and Seung multiplicative update...');
% [ G_WMU_NMF , F_WMU_NMF ] = WMU_NMF( W , X , Ginit , Finit , Iter_max );

% EM_WMU_NMF
fprintf('\nRunning EM_WMU_NMF...');
[ G_EM_WMU_NMF , F_EM_WMU_NMF ] = WNMFLibrary.EM_WMU_NMF( W , X , Ginit , Finit , Iter_max,Iter_max_E_step );

% EM_WNE_NMF
fprintf('\nRunning EM_WNE_NMF...');
[ G_EM_WNE_NMF , F_EM_WNE_NMF ] = WNMFLibrary.EM_WNE_NMF( W , X , Ginit , Finit ,Iter_max_E_step );


