% ================ Created on 10/01/2022 by D.Chi ================
function new_image_raw_meas = PCCorrection(image_meas_raws_reshaped1,PCs_all,isReflected_imgMeas)
new_image_raw_meas = zeros(size(image_meas_raws_reshaped1));
for m = 1:size(image_meas_raws_reshaped1,5)
    image_raw_meas = image_meas_raws_reshaped1(:,:,:,:,m);
    for slc = 1:size(image_raw_meas,4)
        image_raw_meas_slc = image_raw_meas(:,:,:,slc);
        new_image_raw = zeros(size(image_raw_meas_slc));
        for ch = 1:32
            nADC = size(image_raw_meas,1);
%             pc_even = phscor_even(:,ch,slc);
%             pc_odd = phscor_odd(:,ch,slc);
%             pc_odd_img = fftshift(ifft(fftshift(pc_odd,1),[],1),1);
%             pc_even_img = fftshift(ifft(fftshift(pc_even,1),[],1),1);
%             cmplx_diff = pc_even_img.*conj(pc_odd_img);
%             cmplx_slope = cmplx_diff(2:end).*conj(cmplx_diff(1:end-1));
%             mslope_phs = angle(sum(cmplx_slope(:)));
            mslope_phs = PCs_all(ch,slc);
            mslope_phs2 = mslope_phs/2;
            beta_linear = mslope_phs2.*(linspace(0,nADC-1,nADC) - nADC/2);
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
end