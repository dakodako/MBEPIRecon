% ================ Created on 10/01/2022 by D.Chi ================
function PCs_all = estimatePC(phscor_even, phscor_odd)
phscor_even_ft = fftshift(ifft(fftshift(phscor_even,1),[],1),1);
phscor_odd_ft = fftshift(ifft(fftshift(phscor_odd,1),[],1),1);
nADC = size(phscor_even,1);
x = linspace(0,nADC-1,nADC) - nADC/2;
ctr = nADC/2-1;
x = x((ctr-22):(ctr+21))';%39:64 for human
X = zeros(length(x),2);
X(:,1) = 1;
X(:,2) = x;
phsDiff_PC_all = phscor_even_ft./phscor_odd_ft;
for s = 1:size(phscor_even,3)
    for ch = 1:32

        Y = angle(phsDiff_PC_all(ctr-22:ctr+21,ch,s));
        Bs = X\Y;
        mslp = Bs(2);
        PCs_all(ch,s) = mslp;
        
    end
end
end