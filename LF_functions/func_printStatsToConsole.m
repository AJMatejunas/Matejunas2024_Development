function outData = func_printStatsToConsole(varName,varUnit,inData,QSVal)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 6/9/2017

outData.mean = mean(inData);
outData.std = std(inData);
outData.covPc = outData.std/outData.mean*100;
outData.pcDiffToQS = (outData.mean-QSVal)/QSVal*100;

fprintf('\n-----------------------------------------------------------------\n')
fprintf('STATS for %s\n',varName)
fprintf('%s mean: \t %.2f %s\n',varName,outData.mean,varUnit)
fprintf('%s STD: \t %.2f %s\n',varName,outData.std,varUnit)
fprintf('%s COV: \t %.2f %% \n',varName,outData.covPc)
fprintf('%s Diff to QS: \t %.2f %% \n',varName,outData.pcDiffToQS)
fprintf('-----------------------------------------------------------------\n\n')

end

