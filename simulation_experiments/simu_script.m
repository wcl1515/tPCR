%% please go to the folder 'data' where stores the data you need

tframes=1:50;% it can change from 1:2 to 1:1360

%% genearation kspace data and decomposition
simu_kspace_generation_svd(tframes,48,1);

%% SR
traj=48;%trajectory number, i.e 'traj2d_48.grad'
off=1;% field map scaling factor
tol=5e-4;%tolerance
iter_ex=[0:-10:-100];

L='l2';%regularization types, 'l2' or 'l1'
lambda=0.05;%regularization parameter
sr_recon_simu(tframes,L,lambda,traj,off,tol) ;

L='l1';
lambda=5e-6;
sr_recon_simu(tframes,L,lambda,traj,off,tol) ;% SR

%% tPCR

% step 1: decomposition
% already done

% step 2: reconstruction
L='l2';%regularization types, 'l2' or 'l1'
lambda=0.05;%regularization parameter
tpcr_recon_simu(tframes,L,lambda,traj,off,tol,iter_ex);

L='l1';
lambda=6e-5;
tpcr_recon_simu(tframes,L,lambda,traj,off,tol,iter_ex);

% step 3: composition
tpcr_combination(tframes,'sr',1,'l2_0.01','simu/simu_1/simu_48')
tpcr_combination(tframes,'tpcr',1:10,'l2_0.01','simu/simu_1/simu_48')
tpcr_combination(tframes,'sr',1,'l1_6e-05','simu/simu_1/simu_48')
tpcr_combination(tframes,'tpcr',1:10,'l1_6e-05','simu/simu_1/simu_48')

%% comparing error
error_comparison(length(tframes))