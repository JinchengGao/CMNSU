clear,clc,close all
addpath("./data")

% Load data
abundance = load('dc1_75.mat').Xtrue;
%abundance = abundance(2:10,:);
E = load('dc1_75.mat').A;
%E = E(:,2:10);
y_real = E*abundance;
y = load("DC1_case3.mat").y;
% y = addSNRNoise(y_real,[25,35]); %对每个波段添加高斯噪声
% y = addSpNoise(y,0.05);

iter = 300;
tol = 1e-4;
bands = 224;
R = 240;
n_row = 75;
n_col = 75;

lambda_sunsal = 1;
tic
[X_sunsal,res_p_sunsal,res_d_sunsal] = sunsal(E,y,'LAMBDA',lambda_sunsal,'AL_ITERS',iter, ...
    'POSITIVITY','yes','TOL',1e-4);
toc
SRE_sunsal = SRE(abundance,X_sunsal);
fprintf('Sunsal SRE:%f\n',SRE_sunsal)

param_sunsal_tv.lambda1 = 1e-1;
param_sunsal_tv.lambda2 = 1e-1;
tic
[X_sunsal_tv] = sunsal_tv(E, y, 'lambda_1',param_sunsal_tv.lambda1 , ...
        'lambda_tv', param_sunsal_tv.lambda2, 'IM_SIZE', [75, 75], 'ADDONE', 'no', 'POSITIVITY', 'yes', ...
        'AL_iters', iter, 'verbose', 'no', 'MU',  0.1);
toc
SRE_sunsal_tv = SRE(abundance, X_sunsal_tv);
fprintf('Sunsal TV SRE:%f\n',SRE_sunsal_tv)

param_CMN.mu = 5;
param_CMN.gamma = 5;
tic
[X_CMN] = CMN(E, y, 'MU', param_CMN.mu,'GAMMA',param_CMN.gamma,'AL_ITERS', iter,'STEP',1.1);
toc
SRE_CMN = SRE(abundance, X_CMN);
fprintf('CMN Sunsal SRE:%f\n',SRE_CMN)

param_CMN_TV.mu1 = 5e-1;
param_CMN_TV.mu2 = 1e-1;
param_CMN_TV.gamma = 10;
tic
%If you have GPU,you can use CMN_TV_GPU,It is faster than CMN_TV
% [X_CMN_TV] = CMN_TV_GPU (E, y, 'MU_1', param_CMN_TV.mu1,'MU_2', param_CMN_TV.mu2,'GAMMA',param_CMN_TV.gamma, ...
%     'AL_ITERS', iter,'STEP',1.1,'IM_SIZE',[75,75]);
[X_CMN_TV] = CMN_TV(E, y, 'MU_1', param_CMN_TV.mu1,'MU_2', param_CMN_TV.mu2,'GAMMA',param_CMN_TV.gamma, ...
    'AL_ITERS', iter,'STEP',1.1,'IM_SIZE',[75,75]);
toc
SRE_CMN_TV = SRE(abundance, X_CMN_TV);
fprintf('CMN Sunsal TV SRE:%f\n',SRE_CMN_TV)


figure('Position', [100, 100, 1500, 400]); % 调整窗口的大小
[ha, pos] = tight_subplot(1, 5, [.025 .015], [.05 .05], [.05 .15]);

axes(ha(1));
abun_3D = permute(reshape(abundance,[R,n_row,n_col]),[2,3,1]);
imagesc(abun_3D(:,:,2), [0 1])
axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])
title('GT')


axes(ha(2));
X_sunsal_3D = permute(reshape(X_sunsal,[R,n_row,n_col]),[2,3,1]);
imagesc(X_sunsal_3D(:,:,2), [0 1])
axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])
title('Sunsal')

axes(ha(3));
X_sunsal_tv_3D = permute(reshape(X_sunsal_tv,[R,n_row,n_col]),[2,3,1]);
imagesc(X_sunsal_tv_3D(:,:,2), [0 1])
axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])
title('Sunsal TV')

axes(ha(4));
X_CMN_3D = permute(reshape(X_CMN,[R,n_row,n_col]),[2,3,1]);
imagesc(X_CMN_3D(:,:,2), [0 1])
axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])
title('CMN')


axes(ha(5));
X_CMN_TV_3D = permute(reshape(X_CMN_TV,[R,n_row,n_col]),[2,3,1]);
imagesc(X_CMN_TV_3D(:,:,2), [0 1])
axis square, set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])
title('CMN TV')


colormap jet


