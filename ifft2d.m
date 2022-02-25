% ================ Created on 14/09/2021 by D.Chi ================
function imgV = ifft2d(imgF,dim1, dim2)
imgV = fftshift(fftshift(ifft(fftshift(ifft(fftshift(imgF,dim1),[],dim1),dim2),[],dim2),dim1),dim2);
