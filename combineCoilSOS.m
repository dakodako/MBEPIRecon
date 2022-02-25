% ================ Created on 14/09/2021 by D.Chi ================
function imgCmb = combineCoilSOS(imgChs, dimCh)
imgCmb = squeeze(sqrt(sum(abs(imgChs).^2,dimCh)));