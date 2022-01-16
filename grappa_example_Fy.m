load sens.mat
sens1 = squeeze(sens(:,:,1,:));
img_full = img.*sens1;
img_full_Fy = fftshift(fft(img_full,[],2),2);
img_ds_Fy = complex(zeros(size(img_full_Fy)));
img_ds_Fy(:,1:2:end) = img_full_Fy(:,1:2:end);
calib = img_full_Fy(:,33-4:33+5,:);
img_ds_Fy_permute = permute(img_ds_Fy,[3,1,2]);
calib = img_full_Fy(:,33-4:33+5,:);
calib_permute = permute(calib,[3,1,2]);
recon_Fy   =   grappa(img_ds_Fy_permute, calib_permute, R, [3,4]);