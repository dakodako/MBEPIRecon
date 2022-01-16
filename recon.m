function recon_images = recon(twix)
twix.image.flagIgnoreSeg = true;
twix.image.flagRemoveOS = true;
twix.refscan.flagRemoveOS = true;
twix.refscan.flagIgnoreSeg = true;
%%
raw_data = squeeze(twix.image());
out_raw = zeros(74,32,74,12);
for s = 1:12
raw_data_meas1 = raw_data(:,:,:,s,1);
% raw_data_meas1_padded = zeros([size(raw_data_meas1,1), size(raw_data_meas1,2), size(raw_data_meas1,3) + 9, size(raw_data_meas1,4)]);
% %%
% raw_data_meas1_padded(:,:,1:65,:) = raw_data_meas1;
% raw_data_meas1_padded(:,:,66:end,:) = raw_data_meas1(:,:,1:9,:);
%%

refscan = squeeze(twix.refscan());
refscan_meas1 = refscan(:,:,:,s);

refscan_meas1_rotated = zeros([size(refscan_meas1,3),size(refscan_meas1,2),size(refscan_meas1,1)]);
raw_data_meas1_rotated = zeros([size(raw_data_meas1,3),size(raw_data_meas1,2),size(raw_data_meas1,1)]);

for i = 1:32
    refscan_meas1_rotated(:,i,:) = imrotate(squeeze(refscan_meas1(:,i,:)),90);
    raw_data_meas1_rotated(:,i,:) = imrotate(squeeze(raw_data_meas1(:,i,:)),90);
end

%%
raw_data_meas1_temp = zeros((size(refscan_meas1_rotated,1)-2)*(size(refscan_meas1_rotated,3)-2),9,32);%zeros(2592, 32,9);
kNo = 1; % kernel/patch number
for ny = 1:(size(refscan_meas1_rotated,1)-2)
    for nx = 1:(size(refscan_meas1_rotated,3)-2)
        for nc = 1:32
            temp_data = squeeze(refscan_meas1_rotated(ny:ny+2,nc,nx:nx+2));
            raw_data_meas1_temp(kNo,:,nc) = reshape(temp_data',[1,9]);
        end
        %raw_data_meas1_temp(kNo,:,:) = ...
            %reshape(refscan_meas1_rotated(ny:ny+2,:,nx:nx+2),[1,32,9]); % ch1: each patch 3x3 is reshaped into vector and put into matrix one line after another
%         S_ch_2_temp(kNo,:) = ...
%             reshape(phantom_ch_2_k_acl(ny:ny+2,nx:nx+2)',[1,9]); 
%         S_ch_3_temp(kNo,:) = ...
%             reshape(phantom_ch_3_k_acl(ny:ny+2,nx:nx+2)',[1,9]); 
%         S_ch_4_temp(kNo,:) = ...
%             reshape(phantom_ch_4_k_acl(ny:ny+2,nx:nx+2)',[1,9]); 
%         S_ch_5_temp(kNo,:) = ...
%             reshape(phantom_ch_5_k_acl(ny:ny+2,nx:nx+2)',[1,9]); 
%         S_ch_6_temp(kNo,:) = ...
%             reshape(phantom_ch_6_k_acl(ny:ny+2,nx:nx+2)',[1,9]); 
        kNo = kNo + 1; % to move through all patches
    end
end

%%
S_data_meas1 = squeeze(raw_data_meas1_temp(:,[1:3,7:9],:));
T_data_meas1 = squeeze(raw_data_meas1_temp(:,5,:));
%S = reshape(S_data_meas1, [size(S_data_meas1,1),size(S_data_meas1,2)*size(S_data_meas1,3)]);
S = [squeeze(S_data_meas1(:,:,1)) squeeze(S_data_meas1(:,:,2)) squeeze(S_data_meas1(:,:,3)) squeeze(S_data_meas1(:,:,4)) ...
    squeeze(S_data_meas1(:,:,5)) squeeze(S_data_meas1(:,:,6)) squeeze(S_data_meas1(:,:,7)) squeeze(S_data_meas1(:,:,8)) ...
    squeeze(S_data_meas1(:,:,9)) squeeze(S_data_meas1(:,:,10)) squeeze(S_data_meas1(:,:,11)) squeeze(S_data_meas1(:,:,12)) ...
    squeeze(S_data_meas1(:,:,13)) squeeze(S_data_meas1(:,:,14)) squeeze(S_data_meas1(:,:,15)) squeeze(S_data_meas1(:,:,16)) ...
    squeeze(S_data_meas1(:,:,17)) squeeze(S_data_meas1(:,:,18)) squeeze(S_data_meas1(:,:,19)) squeeze(S_data_meas1(:,:,20)) ...
    squeeze(S_data_meas1(:,:,21)) squeeze(S_data_meas1(:,:,22)) squeeze(S_data_meas1(:,:,23)) squeeze(S_data_meas1(:,:,24)) ...
    squeeze(S_data_meas1(:,:,25)) squeeze(S_data_meas1(:,:,26)) squeeze(S_data_meas1(:,:,27)) squeeze(S_data_meas1(:,:,28)) ...
    squeeze(S_data_meas1(:,:,29)) squeeze(S_data_meas1(:,:,30)) squeeze(S_data_meas1(:,:,31)) squeeze(S_data_meas1(:,:,32))];

W = pinv(S) * T_data_meas1;
%%
clear raw_data_meas1_temp
%%
KNo = 1;
S_new_temp = zeros(72*32, 9,32);
for ny = 1:2:(65-2)
    for nx = 1:(74-2)
        for nc = 1:32
            temp_data = squeeze(raw_data_meas1_rotated(ny:ny + 2, nc,nx:nx + 2));
            S_new_temp(KNo,:,nc) = reshape(temp_data',[1,9]);
        end
       % S_new_temp(KNo,:, :) = reshape(raw_data_meas1_rotated(ny:ny+2,:,nx:nx+2),[1,32,9]);
        KNo = KNo + 1;
    end
end

S_new = S_new_temp(:,[1:3,7:9],:);
%%
clear S_new_temp
S_new = reshape(S_new, [size(S_new,1),size(S_new,2)*size(S_new,3)]);
%%
T_new = S_new*W;
T_new_M = reshape(T_new, [72 32 32]);
for nc = 1:32
    T_new_M_transpose(:,:,nc) = T_new_M(:,:,nc)';
end
P_f_u_new = raw_data_meas1_rotated;
for nc = 1:32
    P_f_u_new(2:2:end-1,nc,2:end-1) = T_new_M_transpose(:,:,nc);
end
%imgPF = squeeze(P_f_u_new(:,1,:));
phs = zeros(size(P_f_u_new,3),32,size(P_f_u_new,3));
phs(1:65,:,:) = P_f_u_new;
%img_conj = ifftdim(phs, 2);
for nc = 1:32
    phs(66:end,nc,:) = flip(conj(squeeze(P_f_u_new(2:10,nc,:))),1);
end
%phs(66:end,:,:) = flip(conj(imgPF(2:10,:)),1);
%img_conj = abs(ifftdim(phs,2));
out_raw(:,:,:,s) = phs;
end
%%
recon_images = zeros(74,74,12);
for nS = 1:12
    temp = zeros(74,32,74);
    for nC = 1:32
        sli_chn_raw = squeeze(out_raw(:,nC,:,nS));
        temp(:,nC,:) = fftshift(ifft2(sli_chn_raw));
    end
    sli_img = sqrt(sum(temp.^2,2));
    %sli_img = sum(sqrt(temp),2);
    recon_images(:,:,nS) = sli_img;
end
end
%%
function out = nrmse(a,b)
    out = norm(a(:)-abs(b(:)))/norm(abs(b(:)));
end

function show_pair(dataL, cscaleL, dataR, cscaleR)
    clf();
    subplot('position',[0   0   0.5 1]);imshow(dataL,cscaleL)
    subplot('position',[0.5 0   0.5 1]);imshow(dataR,cscaleR)
end
function plt(data,label)
    clf();
    plot(data,'linewidth',2);
    grid on;
    title(label);
end

function c = img_diff(a,b)
    c = abs(a - abs(b));
end

function x = fftdim(x,dims)
    for i = dims
        x = fftshift(fft(ifftshift(x,i),[],i),i);
    end
end
function x = ifftdim(x,dims)
    for i = dims
        x = fftshift(ifft(ifftshift(x,i),[],i),i);
    end
end

function w = hann(N)
    w = 0.5*(1-cos(2*pi*(0:N-1)'/(N-1)));
end

function b = real_nonneg(a)
    b = max(real(a),0);
end