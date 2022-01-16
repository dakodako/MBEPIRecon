function padded_kspace = partialFourier(kspace,endLin,Niter)
    [M,I] = max(kspace,[],'all','linear');
    [m,n] = ind2sub(size(kspace),I);
    
    if(n - (endLin-n-1) < 1|| (endLin - n < 2))
        centerLin = round(endLin/2);
        gapLin = round((endLin - centerLin)/2);
%         padded_kspace = kspace;
%         padded_kspace(:,endLin:end) = 0;
%         return
    else
        centerLin = n;
        gapLin = endLin- centerLin;
    end
%     centerLin = n;
    % endLin = 117;
%     Niter = 10;

    phase_part = zeros(size(kspace));
    phase_part(:,centerLin - (gapLin-1):centerLin + (gapLin)) = kspace(:,centerLin - (gapLin-1):centerLin + (gapLin));
    p = angle(fftshift(fftshift(ifft(fftshift(ifft(fftshift(phase_part,1),[],1),2),[],2),1),2));
    temp_raw_rotate = kspace;
    for i = 1:Niter
        m = abs(fftshift(fftshift(ifft(fftshift(ifft(fftshift(temp_raw_rotate,1),[],1),2),[],2),1),2));
        new_temp_raw_rotate = m.*exp(1j.*p);
        new_temp_raw_rotate_kspace = fftshift(fftshift(fft(fftshift(fft(fftshift(new_temp_raw_rotate,1),[],1),2),[],2),1),2);
        temp_raw_rotate(:,endLin+1:end) = new_temp_raw_rotate_kspace(:,endLin+1:end);
    end
    padded_kspace = temp_raw_rotate;
end