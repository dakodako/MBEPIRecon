filename = 'meas_MID00170_FID14366_cmrr_mbep2d_bold_1mm_iso_mb6_120slc_TE43TR2500_PA-039.dat';
twix_obj = mapVBVD(filename);
%%
NSlc = twix_obj.image.sqzSize(4);
MB = twix_obj.hdr.Meas.NMBSliceBands;
NSlc_mb = NSlc/MB;
numCha = twix_obj.image.sqzSize(2);
RefLinesPE = twix_obj.hdr.Meas.NRefLin;
RawLin = twix_obj.hdr.Meas.RawLin;
NMeas = twix_obj.image.sqzSize(5);
nADC = twix_obj.hdr.Meas.NColMeas/2;
nLin = RawLin;
twix_obj.image.flagRemoveOS = true;
twix_obj.image.flagRampSampRegrid = true;
%%
image_raws = squeeze(twix_obj.image.unsorted());
%%
image_meas_raws = image_raws(:,:,nLin*NSlc + 1:end);
image_meas_raws_reshaped = reshape(image_meas_raws,[nADC, numCha, nLin,NSlc_mb, NMeas]);
clear image_meas_raws
%%
for m = 1:NMeas
    image_meas_raw_reshaped = image_meas_raws_reshaped(:,:,:,:,m);
    filename = strcat('raw/image_meas_',num2str(m));
    save(filename,'image_meas_raw_reshaped');
end
%%
sg_refscan_raws = image_raws(:,:,1:nLin*NSlc); % 57 NLin % calibration scan for slice-grappa
save('raw/sg_refscan_raws','sg_refscan_raws');
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
isReflected_img = twix_obj.image.IsReflected();
isReflected_img_meas = isReflected_img(nLin * NSlc + 1:end);
isReflected_phscor = twix_obj.phasecor.IsReflected();
isReflected_phscor_img_meas = isReflected_phscor(3*NSlc+1:end);
%%
isReflected_phscor_sgRef = isReflected_phscor(1:NSlc*3);
isReflected_sgRef = isReflected_img(1:nLin * NSlc);
%%
isReflected_imgMeas = reshape(isReflected_img_meas,[nLin,NSlc_mb,NMeas]);
isReflected_phscor_imgMeas = reshape(isReflected_phscor_img_meas,[3,NSlc_mb,NMeas]);
isReflected_phscor_sgRef = reshape(isReflected_phscor_sgRef,[3,NSlc]);
isReflected_sgRef = reshape(isReflected_sgRef,[nLin,NSlc]);
%%
sg_refscan_raws_reshaped = reshape(sg_refscan_raws,[nADC, numCha, nLin,NSlc]);
%%
sg_phscor_raws = phasecor_raws(:,:,1:3*NSlc);
sg_phscor_raws = reshape(sg_phscor_raws,[nADC,numCha,3,NSlc]);
%% phase correction on slice-grappa reference data
PCs_gp_all_new = estimatePC(squeeze(sg_phscor_raws(:,:,1,:) + sg_phscor_raws(:,:,3,:))./2,squeeze(sg_phscor_raws(:,:,2,:)));
new_sgRef_raw = PCCorrection(sg_refscan_raws_reshaped,PCs_gp_all_new,isReflected_sgRef);
%% read calibration scan for GRAPPA
twix_obj.refscan.flagRemoveOS = true;
twix_obj.refscan.flagRampSampRegrid = true;
gp_refscan_raws = squeeze(twix_obj.refscan());
%%
% gp_refscan_raws = reshape(gp_refscan_raws,[nADC,numCha,RefLinesPE,NSlc]);
% twix_obj.refscanPC.flagRemoveOS = true;
% twix_obj.refscanPC.flagRampSampRegrid = true;
% gp_phs_raws = squeeze(twix_obj.refscanPC.unsorted());
gp_phs_raws = phscor_raws_reshaped(:,:,:,1:NSlc);
gp_phscor_raws = reshape(gp_phs_raws,[nADC,numCha,3,NSlc]);
isReflected_gp = twix_obj.refscan.IsReflected();
isReflected_gp = reshape(isReflected_gp,[RefLinesPE NSlc]);
%%
%% phase correction on grappa reference data

PCs_gp_all_new = estimatePC(squeeze(gp_phscor_raws(:,:,1,:) + gp_phscor_raws(:,:,3,:))./2,squeeze(gp_phscor_raws(:,:,2,:)));
new_gpRef_raw = PCCorrection(gp_refscan_raws,PCs_gp_all_new,isReflected_gp);
%%
phsShift = zeros(size(new_sgRef_raw(:,:,:,1)));
phsShift(:,:,1:2:end) = exp(1i*pi/2);
phsShift(:,:,2:2:end) = exp(-1i*pi/2);
new_new_sgRef_raw = new_sgRef_raw;
new_new_sgRef_raw(:,:,:,21:40) = new_sgRef_raw(:,:,:,21:40).*phsShift;
new_new_sgRef_raw(:,:,:,61:80) = new_sgRef_raw(:,:,:,61:80).*phsShift;
new_new_sgRef_raw(:,:,:,101:120) = new_sgRef_raw(:,:,:,101:120).*phsShift;
%%
new_sgRef_raw_zp = zeros([nADC,numCha,nLin * 2,NSlc]);
new_sgRef_raw_zp(:,:,1:2:end,:) = new_new_sgRef_raw;
new_sgRef_raw_zp_per = permute(new_sgRef_raw_zp,[2,1,3,4]);
%%
% new_gpRef_raw_per = permute(new_gpRef_raw,[2,1,3,4]);
mapping = [2:2:NSlc,1:2:NSlc];
gp_refscan_raws_mapped = zeros(size(gp_refscan_raws));
gp_refscan_raws_mapped(:,:,:,mapping) = gp_refscan_raws;
%%

new_gpRef_raw_per = permute(gp_refscan_raws_mapped,[2,1,3,4]);
R = [1,2];
kernel = [3,4];
%% run GRAPPA to get the single band reference images
recon_sgRef = zeros([numCha nADC nLin*2 NSlc]);
for s = 1:size(new_sgRef_raw_zp,4)
    recon_sgRef(:,:,:,s) = grappa(new_sgRef_raw_zp_per(:,:,:,s),new_gpRef_raw_per(:,:,:,s),R,kernel);
end
%%
figure;montage(imrotate(combineCoilSOS(ifft2d(recon_sgRef,2,3),1),-90),'DisplayRange',[0 0.5e-5])
%%
% mapping = [2:2:NSlc,1:2:NSlc];
mapping_back = zeros(size(mapping));
mapping_back(2:2:end) = 1:NSlc/2;
mapping_back(1:2:end) = NSlc/2+1:NSlc;
% MB = 4;
% NLins = 57;
slc_order = twix_obj.image.Sli(nLin*NSlc +1:nLin:nLin*NSlc + NSlc_mb*nLin);
%%
curr_slc = mapping(slc_order(1));
mb_slcs = curr_slc:(NSlc/MB):NSlc;
curr_mb_slcs = mapping_back(mb_slcs);
%%
new_new_new_sgRef_raw(:,:,:,mapping_back) = new_sgRef_raw;
%%
new_new_gpRef_raw_per(:,:,:,mapping_back) = new_gpRef_raw_per;
%%
for m = 1:1
    filename = strcat('raw/image_meas_',num2str(m)); % load img raw data
    load(filename);
    PCs_all = estimatePC(phscor_even(:,:,:,m),phscor_odd(:,:,:,m));
    % phase correction on the image measurement
    new_image_raw_meas = PCCorrection(image_meas_raw_reshaped,PCs_all,isReflected_imgMeas);
    for i_s = 1:length(slc_order)
        curr_slc = mapping(slc_order(i_s));
        mb_slcs = curr_slc:(NSlc/MB):NSlc;
        curr_mb_slcs = mapping_back(mb_slcs);
        %
%         new_sgRef_raw_per = permute(new_sgRef_raw,[2,1,3,4]);
        new_sgRef_raw_per = permute(new_new_new_sgRef_raw,[2,1,3,4]);
        new_sgRef_raw_per_slcs = new_sgRef_raw_per(:,:,:,curr_mb_slcs);
        new_img_raw_meas_slc = new_image_raw_meas(:,:,:,i_s,1);
        new_img_raw_meas_slc_per = permute(new_img_raw_meas_slc,[2,1,3]);
        %
        phsShift = zeros(size(new_img_raw_meas_slc_per));
        phsShift(:,:,1:2:end) = exp(1i*pi/2);
        phsShift(:,:,2:2:end) = exp(-1i*pi/2);

       
        new_sgRef_raw_per_slcs_ps = new_sgRef_raw_per_slcs;
%         new_sgRef_raw_per_slcs_ps(:,:,:,2) = new_sgRef_raw_per_slcs_ps(:,:,:,2).*phsShift;
%         new_sgRef_raw_per_slcs_ps(:,:,:,4) = new_sgRef_raw_per_slcs_ps(:,:,:,4).*phsShift;
%         new_sgRef_raw_per_slcs_ps(:,:,:,6) = new_sgRef_raw_per_slcs_ps(:,:,:,6).*phsShift;
        
        % new_sgRef_raw_per_slcs_ps(:,:,:,4) = new_sgRef_raw_per_slcs_ps(:,:,:,4).*conj(phsShift);
        % run slice-GRAPPA
        out_sg = sg(new_img_raw_meas_slc_per,new_sgRef_raw_per_slcs_ps,[3,3]);
        % out_sg = sg(new_img_raw_meas_slc_per,new_sgRef_raw_per_slcs,[5,5]);
        %
        out_sg_psBack = out_sg;
        out_sg_psBack(:,:,:,2) = out_sg_psBack(:,:,:,2).*phsShift;
        out_sg_psBack(:,:,:,4) = out_sg_psBack(:,:,:,4).*phsShift;
        out_sg_psBack(:,:,:,6) = out_sg_psBack(:,:,:,6).*phsShift;
        %
        out_sg_zp = zeros([numCha nADC nLin*2 MB]);
        out_sg_zp(:,:,1:2:end,:) = out_sg_psBack;
        % run GRAPPA
        for s = 1:MB
            calib_1 = new_new_gpRef_raw_per(:,:,:,curr_mb_slcs(s));
            input = out_sg_zp(:,:,:,s);
            outAll(:,:,:,mb_slcs(s)) = grappa(input, calib_1, R, kernel);
        end
    end
%     outAllzp = zeros(numCha,nADC,nADC,NSlc);
%     outAllzp(:,:,2:end-1,:) = outAll;

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
    outAllimg(:,:,:,m) = combineCoilSOS(ifft2d(outAll,2,3),1);
end