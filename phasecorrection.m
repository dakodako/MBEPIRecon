function corrected_data = phasecorrection(phasecor_data,input_data)
[NCol, NCha, NRow, NSlc, NMeas] = size(input_data);
corrected_data = zeros([NCol, NCha, NRow, NSlc, NMeas]);
pc_range = 10;
for m = 1:NMeas
    for sl = 1:NSlc
        image_data_meas_slc = squeeze(input_data(:,:,:,sl,m));
        phasecor_data_meas_slc = squeeze(phasecor_data(:,:,:,sl,m));
%         new_image_data_meas_slc = zeros(size(image_data_meas_slc));
        for ch = 1:NCha
            S2 = phasecor_data_meas_slc(:,ch,2);
            S1 = phasecor_data_meas_slc(:,ch,1);
            S3 = phasecor_data_meas_slc(:,ch,3);
            S2_s = (S1+S3)./2;
            s2_s = ifft(S2_s);
            s2 = ifft(S2);
            Ra = angle(fftshift(s2)./fftshift(s2_s));
            R_slope = (Ra(NCol/2+pc_range) - Ra(NCol/2-(pc_range-1)))./(NCol/2+pc_range - (NCol/2-(pc_range-1)));
            beta = R_slope/2;
            beta_s = -R_slope/2;
%             [s2_autocorr,lags] = xcorr(conj(s2),circshift(conj(s2),-1));
%             beta = angle(s2_autocorr(lags == 0));
%             [s2_s_autocorr,lags] = xcorr(conj(s2_s),circshift(conj(s2_s),-1));
%             beta_s = angle(s2_s_autocorr(lags == 0));
%             s2_autocorr = conj(s2).*circshift(s2,-1);
%             s2_s_autocorr = conj(s2_s).*circshift(s2_s,-1);
%             beta = angle(sum(s2_autocorr(2:end)));
%             beta = -0.0057;
%             beta_s = 0.0057;
%             beta_s = angle(sum(s2_s_autocorr(2:end)));
            beta_linear = beta.*(linspace(0,size(input_data,1)-1,size(input_data,1))-size(input_data,1)/2);
            beta_s_linear = beta_s.*(linspace(0,size(input_data,1)-1,size(input_data,1))-size(input_data,1)/2);
            s2_corrected = fftshift(s2).*exp(-1i.*beta_linear');
            s2_s_corrected = fftshift(s2_s).*exp(-1i.*beta_s_linear');
            s2_xcro = sum(conj(s2_corrected).*s2_s_corrected);
%             phi = angle(s2_xcro);
            phi = 0;
            for l = 1:NRow
                
                S = image_data_meas_slc(:,ch,l);
                s = fftshift(ifft(S));
                if(mod(l,2) == 1)
%                         corrected = fftshift(fft(fftshift(s.*exp(-1i*beta_linear').*exp(1i*phi))));
                    corrected_data(:,ch,l,sl,m) = (fft(fftshift(s.*exp(-1i*beta_linear').*exp(1i*phi))));
%                     new_image_data(:,ch,l,sl,m) = s;
                else
%                             corrected = fftshift(fft(fftshift(s.*exp(-1i*beta_s_linear').*exp(1i*phi))));
                    corrected_data(:,ch,l,sl,m) = (fft(fftshift(s.*exp(-1i*beta_s_linear').*exp(1i*phi))));
                end
%                     new_image_data(:,ch,l,sl,m) = corrected;
                
            end
%             new_image_data(:,ch,1:4:end,sl,m) = new_image_data_meas_slc(:,ch,1:4:end,1); % 168
%                 new_image_data(:,ch,2:4:end,sl,m) = new_image_data_meas_slc(:,ch,2:4:end,1); % 98
%                 new_image_data(:,ch,4:4:end,sl,m) = new_image_data_meas_slc(:,ch,4:4:end,2); % 98
%             new_image_data(:,ch,3:4:end,sl,m) = new_image_data_meas_slc(:,ch,3:4:end,2);   %168
        end
    end
end
end