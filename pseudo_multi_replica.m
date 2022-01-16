%%
twix = mapVBVD('meas_MID161_ep2d_pF7_8_Grp2_mtx74_TE19_TR1000_ref24_FID51863.dat');
recon_images = recon(twix);
%%
noise = squeeze(twix.noise());
noise_cov = cov(noise);
%%
I = eye(148);
noise_cov_bar = kron(noise_cov,I);