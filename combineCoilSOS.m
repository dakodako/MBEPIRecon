function imgCmb = combineCoilSOS(imgChs, dimCh)
imgCmb = squeeze(sqrt(sum(abs(imgChs).^2,dimCh)));