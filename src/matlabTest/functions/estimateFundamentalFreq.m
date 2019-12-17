%Function to estimate fundamental frequency
%	@param:		constants, a structure with field samplingFreq containing sampling frequency.
%	@param:		signalIn, the signal to be analysed for fundamental frequency
%	@returns:	fundamentalFreq, a result struct. Field f0 has the exctracted fundamental frequency
%Written by Timo Rantalainen 2014 tjrantal at gmail dot com

function fundamentalFreq = estimateFundamentalFreq(constants,signalIn)
	%remove DC offset
	signalIn = signalIn-mean(signalIn);
	%Apply hann windowing
	hannWindowed = signalIn.*hanning(length(signalIn));
	if size(hannWindowed,1) > size(hannWindowed,2)
		hannWindowed = hannWindowed';
	end
	%Double the signal length by appending zeroes
	appended = [hannWindowed zeros(1,length(hannWindowed))];
	%fftSignal = fft(signalIn);						%Take the fft
    fftSignal = fft(appended);						%Take the fft
	
	fftSignal= fftSignal./(length(fftSignal)/2+1);	%Normalize the fft
	fftSignal = fftSignal(1:length(fftSignal)/2+1);	%Take only the real part
    fftSignal(1) = fftSignal(1)/2;					%First coefficient is double to what it should be
	fftSignal(length(fftSignal)) = fftSignal(length(fftSignal))/2;					%Nyquist coefficient is double to what it 
	power = abs(fftSignal).^2;
	amp = abs(fftSignal);
	freq = [0:(length(fftSignal)-1)]*((constants.samplingFreq/2)/(length(fftSignal)-1));
	
%     figure,plot(freq,power);
%     keyboard;
    %Create fundamental frequency candidates, and figure out corresponding frequency bins
	%Idea taken from Klapuri 2006 ISMIR congress publication "Multiple Fundamental Frequency Estimation by Summing Harmonic Amplitudes."
	f0cands = (1/constants.stepDurationLims(2)):(1/(5*constants.samplingFreq)):(1/constants.stepDurationLims(1));	%Allow stride duration to be between 0.5 to 2.0 s
	%keyboard;
    %Pre-calculate frequency bins to include for a specific f0 candidate
	harmonics = 20;
	halfBinWidth = ((constants.samplingFreq/2)/(length(fftSignal)-1))/2;%;((constants.samplingFreq/2)/(length(power)))/2;
	for i = 1:length(f0cands)
		binIndices = [];
		for h = 1:harmonics
			hIndices = find(freq >= (f0cands(i)*h-halfBinWidth) & freq < (f0cands(i)*h+halfBinWidth));
			binIndices = [binIndices hIndices];
		end
		f0candsFreqBins(i).binIndices = binIndices;
	end
	
	alpha = 0.5;	%Parameter to weigh harmonics, arbitrarily set at 0.5 Hz
	beta = 2;		%Parameter to weigh harmonics, arbitrarily set at 2 Hz
	salience = zeros(1,length(f0cands));
	for i= 1:length(f0cands)
		findices = f0candsFreqBins(i).binIndices;
		j=1:length(findices);
		salience(i) = sum((constants.samplingFreq*freq(findices(j))+alpha)./(j.*constants.samplingFreq.*freq(findices(j))+beta).*power(findices(j)));
    end

%     figure,plot(f0cands,salience);
%     keyboard;
    
    %Fundamental frequency is selected as the highest peak (=step
    %frequency)
	%Select the fundamental frequency
	%salPeaks = detectPeaks(salience,0.3*max(salience),0);
	%Select the correct peak, compare against autocorrelation
	%selectPeak = 1;
	%[ignore, salTempP] = max(salience(salPeaks.inits(selectPeak):salPeaks.ends(selectPeak)));
    [ignore, salInd] = max(salience);
	fundamentalFreq = struct();
	%fundamentalFreq.f0 = f0cands(salPeaks.inits(selectPeak)+salTempP-1);
    fundamentalFreq.f0 = f0cands(salInd); %Multiply by two to have stride frequency
	fundamentalFreq.f0cands = f0cands;
	fundamentalFreq.f0candsFreqBins = f0candsFreqBins;
	fundamentalFreq.power = power;
	fundamentalFreq.amp = amp;
	fundamentalFreq.freq = freq;
	fundamentalFreq.salience = salience;
	fundamentalFreq.halfBinWidth = halfBinWidth;
	fundamentalFreq.harmonics = harmonics;
	if isfield(constants,'debug')
		keyboard;
	end
end