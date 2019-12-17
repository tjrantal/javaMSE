function res = getBoutAnalysisJava(resultantIn,stepDurationLims)
    res = struct();
%     %Entropy analyses
    rcmeo = javaObject('edu.deakin.timo.RefinedCompositeMultiscaleEntropy',resultantIn,4,0.3,20);
    res.rcme = rcmeo.getRCME()';
	rmpeo = javaObject('edu.deakin.timo.RefinedMultiscalePermutationEntropy',resultantIn,4,20);
    res.rmpe = rmpeo.getRMPE()';

    %Get f0    
    temp = struct();
    temp.samplingFreq = 100;
    temp.stepDurationLims = stepDurationLims;  %Stride rate/2 to get step rate
    fFreq = estimateFundamentalFreq(temp,resultantIn);
    res.f0 = fFreq.f0;
    res.duration = length(resultantIn)/(temp.samplingFreq *60);
