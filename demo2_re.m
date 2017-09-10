%restore teapot

clear; close all;

In = im2double(imread('./teapot/teapot_vis.png'));
G = im2double(imread('./teapot/teapot_nir.png'));
T = im2double(imread('./teapot/teapot_res.png'));

maxIter = 6;

% for eps = [0.002, 0.004, 0.006] 
% for eta_sqr = [0.01, 0.001, 0.0001]
% for phi_alpha = [0.7, 0.8, 0.9];
% for phi_eps = [0.0005, 0.001, 0.005]
% for lambda = [1, 3, 5]
% for beta = [0.6, 0.8, 1]
    
eps = 0.006;
eta_sqr = 0.001;
phi_alpha = 0.9;
phi_eps = 0.001;
lambda = 1;
beta = 2;

tic
S = zeros(size(In));
for c = 1:size(S, 3)
S(:,:,c) = cross_field_re(In(:, :, c),G,lambda,beta,eps,eta_sqr,phi_alpha,phi_eps,maxIter);
end
toc

m = psnr(T, S);

filename = sprintf('eps_%0.3f_eta_%0.4f_alpha_%0.1f_eps_%0.4f_lambda_%d_beta_%0.1f_psnr_%.4f.png', ...
    eps, eta_sqr, phi_alpha, phi_eps, lambda, beta, m);
disp(filename);
% imwrite(S, filename);

figure,imshow(In)
figure,imshow(G)
figure,imshow(S)
figure,imshow(T)

% end
% end
% end
% end
% end
% end
