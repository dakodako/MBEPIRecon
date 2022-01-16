filename = '3751/meas_MID00095_FID14279_ep2d_bold_p2_sms4_1p8iso_20slc_10meas_matchAjuVol_offiso.dat';
twix_obj = mapVBVD(filename);
%%
NSlc = 84;
MB = 4;
NSlc_mb = NSlc/MB;
numCha = 32;
NMeas = 116;
nADC = 116;
twix_obj.image.flagRemoveOS = true;
twix_obj.image.flagRampSampRegrid = true;
image_raws = squeeze(twix_obj.image.unsorted());
image_meas_raws = image_raws(:,:,57*NSlc + 1:end);
image_meas_raws_reshaped = reshape(image_meas_raws,[nADC, numCha, 57,NSlc_mb, NMeas]);
clear image_meas_raws
%%
for m = 1:NMeas
    image_meas_raw_reshaped = image_meas_raws_reshaped(:,:,:,:,m);
    filename = strcat('3751/19_ep2d_bold_p2_sms4_1p8iso_20slc_10meas_matchAjuVol_offiso/ImageRawData/image_meas_',num2str(m));
    save(filename,'image_meas_raw_reshaped');
end
sg_refscan_raws = image_raws(:,:,1:57*NSlc); % 57 NLin
save('3751/19_ep2d_bold_p2_sms4_1p8iso_20slc_10meas_matchAjuVol_offiso/ImageRawData/sg_refscan_raws','sg_refscan_raws');
clear image_raws
%%
twix_obj.phasecor.flagRemoveOS = true;
twix_obj.phasecor.flagRampSampRegrid = true;
phscor_raws = squeeze(twix_obj.phasecor.unsorted());
phscor_raws_reshaped = reshape(phscor_raws,[size(phscor_raws,1),size(phscor_raws,2),...
    3,size(phscor_raws,3)/3]);
%%

phscor_meas_raws = phscor_raws_reshaped(:,:,:,NSlc + 1:end);
phscor_meas_raws_mb = reshape(phscor_meas_raws,[size(phscor_meas_raws,1),numCha,...
    size(phscor_meas_raws,3),NSlc_mb,size(phscor_meas_raws,4)/NSlc_mb]);
phasecor_raws = phscor_raws;
phscor_raws = phscor_meas_raws_mb;
%%
% clear phscor_meas_raws phscor_raws_reshaped
% phscor_imgs = mean(squeeze(phscor_raws(:,:,:,:,1:4)),5);
phscor_even1 = squeeze(phscor_raws(:,:,1,:,:));
phscor_even2 = squeeze(phscor_raws(:,:,3,:,:));
phscor_odd = squeeze(phscor_raws(:,:,2,:,:));
phscor_even = (phscor_even1 + phscor_even2)./2;
%%
% PCs_all = estimatePC(phscor_even, phscor_odd);
% %%

%%
% load 3751/19_ep2d_bold_p2_sms4_1p8iso_20slc_10meas_matchAjuVol_offiso/ImageRawData/image_meas_1.mat
% PC_all_meas1 = estimatePC(phscor_even(:,:,:,1),phscor_odd(:,:,:,1));
% new_image_raw_meas = PCCorrection(image_meas_raw_reshaped,PC_all_meas1,isReflected_imgMeas);
%%
% image_meas_raws_reshaped1 = image_meas_raws_reshaped(:,:,:,:,1);
%%
isReflected_img = twix_obj.image.IsReflected();
isReflected_img_meas = isReflected_img(57 * NSlc + 1:end);
isReflected_phscor = twix_obj.phasecor.IsReflected();
isReflected_phscor_img_meas = isReflected_phscor(3*NSlc+1:end);
%%
isReflected_phscor_sgRef = isReflected_phscor(1:NSlc*3);
isReflected_sgRef = isReflected_img(1:57 * NSlc);
%%
isReflected_imgMeas = reshape(isReflected_img_meas,[57,NSlc_mb,NMeas]);
isReflected_phscor_imgMeas = reshape(isReflected_phscor_img_meas,[3,NSlc_mb,NMeas]);
isReflected_phscor_sgRef = reshape(isReflected_phscor_sgRef,[3,NSlc]);
isReflected_sgRef = reshape(isReflected_sgRef,[57,NSlc]);
%% Phase correction on the image measurement data
new_image_raw_meas = zeros(size(image_meas_raws_reshaped1));
for m = 1:size(image_meas_raws_reshaped1,5)
    image_raw_meas = image_meas_raws_reshaped1(:,:,:,:,m);
    for slc = 1:size(image_raw_meas,4)
        image_raw_meas_slc = image_raw_meas(:,:,:,slc);
        new_image_raw = zeros(size(image_raw_meas_slc));
        for ch = 1:32
            nADC = size(image_raw_meas,1);
            pc_even = phscor_even(:,ch,slc);
            pc_odd = phscor_odd(:,ch,slc);
            pc_odd_img = fftshift(ifft(fftshift(pc_odd,1),[],1),1);
            pc_even_img = fftshift(ifft(fftshift(pc_even,1),[],1),1);
            cmplx_diff = pc_even_img.*conj(pc_odd_img);
            cmplx_slope = cmplx_diff(2:end).*conj(cmplx_diff(1:end-1));
            mslope_phs = angle(sum(cmplx_slope(:)));
%             mslope_phs = PCs_all(ch,slc);
            mslope_phs2 = mslope_phs/2;
            beta_linear = mslope_phs2.*(linspace(0,size(pc_even,1)-1,size(pc_odd,1)) - size(pc_odd,1)/2);
            image_raw_1 = squeeze(image_raw_meas_slc(:,ch,:));
            for n = 1:size(image_raw_1,2)
                S = image_raw_1(:,n);
                s = fftshift(ifft(fftshift(S,1),[],1),1);
                if(isReflected_imgMeas(n,slc,m) == 0)
                    new_image_raw(:,ch,n) = fftshift(fft(fftshift(s.*exp(1i * beta_linear'),1),[],1),1);
                else
                    new_image_raw(:,ch,n) = fftshift(fft(fftshift(s.*exp(1i * -beta_linear'),1),[],1),1);
                end
            end
        end
        new_image_raw_meas(:,:,:,slc,m) = new_image_raw;
    end
end
%%
% sg_refscan_raws = image_raws(:,:,1:57*NSlc);
load 3751/19_ep2d_bold_p2_sms4_1p8iso_20slc_10meas_matchAjuVol_offiso/ImageRawData/sg_refscan_raws.mat
sg_refscan_raws_reshaped = reshape(sg_refscan_raws,[nADC, numCha, 57,NSlc]);
%%
sg_phscor_raws = phasecor_raws(:,:,1:3*NSlc);
sg_phscor_raws = reshape(sg_phscor_raws,[nADC,numCha,3,NSlc]);

%% phase correction on slice-grappa reference data
PCs_gp_all_new = estimatePC(squeeze(sg_phscor_raws(:,:,1,:) + sg_phscor_raws(:,:,3,:))./2,squeeze(sg_phscor_raws(:,:,2,:)));
new_sgRef_raw = PCCorrection(sg_refscan_raws_reshaped,PCs_gp_all_new,isReflected_sgRef);

%%
twix_obj.refscan.flagRemoveOS = true;
twix_obj.refscan.flagRampSampRegrid = true;
gp_refscan_raws = squeeze(twix_obj.refscan.unsorted());
%%
gp_refscan_raws = reshape(gp_refscan_raws,[116,numCha,50,NSlc]);
twix_obj.refscanPC.flagRemoveOS = true;
twix_obj.refscanPC.flagRampSampRegrid = true;
gp_phs_raws = squeeze(twix_obj.refscanPC.unsorted());
gp_phscor_raws = reshape(gp_phs_raws,[116,numCha,3,NSlc]);
isReflected_gp = twix_obj.refscan.IsReflected();
isReflected_gp = reshape(isReflected_gp,[50 NSlc]);
%% phase correction on grappa reference data

PCs_gp_all_new = estimatePC(squeeze(gp_phscor_raws(:,:,1,:) + gp_phscor_raws(:,:,3,:))./2,squeeze(gp_phscor_raws(:,:,2,:)));
new_gpRef_raw = PCCorrection(gp_refscan_raws,PCs_gp_all_new,isReflected_gp);

%% constant phase shift in gpRef
new_gpRef_raw_zp = zeros(116,32,116,NSlc);
new_gpRef_raw_zp(:,:,57-25:57+24,:) = new_gpRef_raw;
new_gpRef_raw_zp_slc1 = new_gpRef_raw_zp(:,:,:,1);
[xf,yf] = meshgrid( ((-116/2):1:(116/2-1 ))./(116) , ((-116/2):1:(116/2 - 1))./(116) );
shift = [-10 0]; % slc 1 [-15 0] for slc 2
new_gpRef_raw_zp_slc1_ps = new_gpRef_raw_zp_slc1.*exp(-1i*(2*pi.*(xf*shift(1) + yf*shift(2))));
%%
new_sgRef_raw_zp = zeros([116,32,114,NSlc]);
new_sgRef_raw_zp(:,:,1:2:end,:) = new_sgRef_raw;
new_sgRef_raw_zp_per = permute(new_sgRef_raw_zp,[2,1,3,4]);
new_gpRef_raw_per = permute(new_gpRef_raw,[2,1,3,4]);
R = [1,2];
kernel = [3,4];
%%
recon_sgRef = zeros([32 116 114 NSlc]);
for s = 1:size(new_sgRef_raw_zp,4)
    recon_sgRef(:,:,:,s) = grappa(new_sgRef_raw_zp_per(:,:,:,s),new_gpRef_raw_per(:,:,:,s),R,kernel);
end
%%
mapping = [2:2:NSlc,1:2:NSlc];
mapping_back = zeros(size(mapping));
mapping_back(2:2:end) = 1:NSlc/2;
mapping_back(1:2:end) = NSlc/2+1:NSlc;
MB = 4;
NLins = 57;
slc_order = twix_obj.image.Sli(NLins*NSlc +1:NLins:NLins*NSlc + NSlc_mb*NLins);
%%
curr_slc = mapping(slc_order(1));
mb_slcs = curr_slc:(NSlc/MB):NSlc;
curr_mb_slcs = mapping_back(mb_slcs);

%%
% load 3750/outAllzp_meas4.mat
for m = 1:NMeas
    filename = strcat('3751/19_ep2d_bold_p2_sms4_1p8iso_20slc_10meas_matchAjuVol_offiso/ImageRawData/image_meas_',num2str(m)); % load img raw data
    load(filename);
    PCs_all = estimatePC(phscor_even(:,:,:,m),phscor_odd(:,:,:,m));
    new_image_raw_meas = PCCorrection(image_meas_raw_reshaped,PCs_all,isReflected_imgMeas);
    for i_s = 1:length(slc_order)
        curr_slc = mapping(slc_order(i_s));
        mb_slcs = curr_slc:(NSlc/MB):NSlc;
        curr_mb_slcs = mapping_back(mb_slcs);
        %
        new_sgRef_raw_per = permute(new_sgRef_raw,[2,1,3,4]);
        new_sgRef_raw_per_slcs = new_sgRef_raw_per(:,:,:,curr_mb_slcs);
        new_img_raw_meas_slc = new_image_raw_meas(:,:,:,i_s,1);
        new_img_raw_meas_slc_per = permute(new_img_raw_meas_slc,[2,1,3]);
        %
        phsShift = zeros(size(new_img_raw_meas_slc_per));
        phsShift(:,:,1:2:end) = exp(1i*pi/2);
        phsShift(:,:,2:2:end) = exp(-1i*pi/2);

        %
        new_sgRef_raw_per_slcs_ps = new_sgRef_raw_per_slcs;
        new_sgRef_raw_per_slcs_ps(:,:,:,2) = new_sgRef_raw_per_slcs_ps(:,:,:,2).*phsShift;
        new_sgRef_raw_per_slcs_ps(:,:,:,4) = new_sgRef_raw_per_slcs_ps(:,:,:,4).*phsShift;
        % new_sgRef_raw_per_slcs_ps(:,:,:,4) = new_sgRef_raw_per_slcs_ps(:,:,:,4).*conj(phsShift);
        out_sg = sg(new_img_raw_meas_slc_per,new_sgRef_raw_per_slcs_ps,[3,3]);
        % out_sg = sg(new_img_raw_meas_slc_per,new_sgRef_raw_per_slcs,[5,5]);
        %
        out_sg_psBack = out_sg;
        out_sg_psBack(:,:,:,2) = out_sg_psBack(:,:,:,2).*phsShift;
        out_sg_psBack(:,:,:,4) = out_sg_psBack(:,:,:,4).*phsShift;
        %
        out_sg_zp = zeros([32 116 114 4]);
        out_sg_zp(:,:,1:2:end,:) = out_sg_psBack;
        %
        for s = 1:4
            calib_1 = new_gpRef_raw_per(:,:,:,curr_mb_slcs(s));
            input = out_sg_zp(:,:,:,s);
            outAll(:,:,:,mb_slcs(s)) = grappa(input, calib_1, R, kernel);
        end
    end
    outAllzp = zeros(32,116,116,84);
    outAllzp(:,:,2:end-1,:) = outAll;
%     if(m == 4)
%         outAllzp4 = outAllzp;
%         save('3750/outAllzp_meas4','outAllzp4');
%     end
%     if(m > 4)
%         phsDiff = outAllzp./outAllzp4;
%         phsDiffAll = squeeze(sum(phsDiff,1));
%         x = ((-nADC/2):1:(nADC/2-1 ))./(nADC);
%         x = x(2:end-1);
%         X = zeros(length(x),2);
%         X(:,1) = 1;
%         X(:,2) = x;
%         for s = 1:NSlc
%             phsDiff_unwrap = unwrap_phase(angle(phsDiffAll(:,:,s)));
%             phsDiff_unwrap_sum = mean(phsDiff_unwrap,1);
%             Y = phsDiff_unwrap_sum(2:end-1)';
%             Bs = X\Y;
%             pShift(s) = Bs(2);
%         end
%         mslp = -mean(pShift);
%         [xf,yf] = meshgrid( ((-nADC/2):1:(nADC/2-1 ))./(nADC) , ((-nADC/2):1:(nADC/2 - 1))./(nADC) );
%         shiftMap = exp(1i.*(xf*mslp + yf*0));
%         outAllzp_per = permute(outAllzp,[2,3,1,4]);
%         outAllzp_shifted = outAllzp_per.*shiftMap;
%     end
    outAllimg(:,:,:,m) = combineCoilSOS(ifft2d(outAllzp,2,3),1);
end