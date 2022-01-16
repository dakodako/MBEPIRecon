%   examples.m
%   mchiew@fmrib.ox.ac.uk
%%
%[mask,twix] = mapVBVD('ep2d_74/meas_MID161_ep2d_pF7_8_Grp2_mtx74_TE19_TR1000_ref24_FID51863.dat');
[mask,twix] = mapVBVD('/Volumes/Samsung_T5/visual/rawdata/3218/os2/98/meas_MID36_Didi_ep2d_bold_os2_mtx98_grp_24_pF_7_8_FID5627.dat');
twix.image.flagRemoveOS = true;
twix.image.flagIgnoreSeg = true;
twix.refscan.flagRemoveOS = true;
twix.refscan.flagIgnoreSeg = true;

twix.image.flagRampSampRegrid = true;
twix.refscan.flagRampSampRegrid = true;
refscan = squeeze(twix.refscan());
raw_img = squeeze(twix.image());
%%
data_temp = new_image_data(:,:,:,1,1);
calib_temp = new_refscan_data(:,:,:,1);
% data = zeros(32,98,86);
% calib = zeros(32,98,32);
% data = zeros(32,74,65);
% calib = zeros(32,74,24);
data = permute(data_temp, [2,1,3]);
calib = permute(calib_temp, [2,1,3]);
% data = zeros(size(data));
% calib_temp = zeros(size(calib_temp));
% for c = 1:size(data,1)
%     for n = 1:size(data,3)
%         if(mod(n,2) ==1)
%             data(c,:,n) = flip(squeeze(data_temp(c,:,n)));
%         else
%             data(c,:,n) = data_temp(c,:,n);
%         end
%     end
% end
% for i = 1:32
%     data(i,:,:) =squeeze(data_temp(:,i,:));
%     calib(i,:,:) = squeeze(calib_temp(:,i,:));
% end

%data = permute(data_temp,[2 1 3]);
%calib = permute(calib_temp,[2 1 3]);
%%
%   Load example data
input   =   matfile('data/data.mat');
truth   =   input.truth;
calib   =   input.calib;
%%
%   ============================================================================
%   The R=2 problem
%   ============================================================================

R       =   [1,2];
kernel  =   [3,4];
%%
mask    =   false(32,96,96);
mask(:,:,1:2:end)   =   true;

data    =   truth.*mask;

recon   =   grappa(data, calib, R, kernel);

show_quad(data, recon, 'R=2');

%%
%   ============================================================================
%   The R=3 problem
%   ============================================================================

R       =   [1,3];
kernel  =   [3,4];

mask    =   false(32,96,96);
mask(:,:,1:3:end)   =   true;

data    =   truth.*mask;

recon   =   grappa(data, calib, R, kernel);

show_quad(data, recon, 'R=3');


%   ============================================================================
%   The R=6 problem
%   ============================================================================

R       =   [1,6];
kernel  =   [3,2];

mask    =   false(32,96,96);
mask(:,:,1:6:end)   =   true;

data    =   truth.*mask;

recon   =   grappa(data, calib, R, kernel);

show_quad(data, recon, 'R=6');


%   ============================================================================
%   The noisy R=6 problem
%   ============================================================================

R       =   [1,6];
kernel  =   [3,2];

mask    =   false(32,96,96);
mask(:,:,1:6:end)   =   true;

noise   =   1E-6*(randn(size(mask)) + 1j*randn(size(mask)));
data    =   (truth + noise).*mask;

recon   =   grappa(data, calib, R, kernel);

show_quad(data, recon, 'R=6 with noise');
