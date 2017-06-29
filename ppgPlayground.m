%% ppgPlayground
%
%  Dave Haas · 18 May 2017
%       Last modified: 2 June 2017:
%       + added PSD plots
%       + changed time series plots to time in seconds, from raw samples
%       + clean-up a lot of maintenance / housekeeping routines

%% Do some housekeeping

clear;

% Define some important variables

figureCount       = 0;       % start a count for auto-labeling figures
ver               = version; % get Matlab version number
testCase          = '9.1';   % is Matlab at least a version 9.1 copy?
newerMatlab       = strncmpi(ver,testCase,3);

%% Set fullPathName that has the data of interest...

% Use Ashley's finger RPOX signal to visualize human PPG patterns
% fullPathName = '/Users/dave/Dropbox/TrialsMay17/CondoRun2.xlsx';
% fullPathName = '/Users/dave/Dropbox/TrialsMay17/June2AshFinger.xlsx';

% Use Dave's locally captured PPG pattern
% fullPathName = '/Users/dave/Documents/MATLAB/_ppgTools/DQ2017/human/20170604_Dave4min.xls';

% Use a T. truncatus sample
% fullPathName = '/Users/dave/Dropbox/TrialsMay17/May26Pos4DQ4Min_S1_M.xlsx';
% fullPathName = '/Users/dave/Dropbox/TrialsMay17/June2DaveDQnew_S1_M.xlsx';
% fullPathName = '/Users/dave/Dropbox/TrialsMay17/June3Pos2DQ4Minnew_S1_M.xlsx';
fullPathName = '/Users/dave/Documents/MATLAB/_ppgTools/DQ2017/June5BlowholeDQNew4Min_S1_M.xlsx';
% fullPathName = '/Users/dave/Documents/MATLAB/_ppgTools/DQ2017/June6BlowHoleDQNew4Min_S1_M.xlsx';


%% Load up some PPG data

% if Excel and mult-sheet, define which sheet of data to use

[fileType, sheetNames, xlsFormat] = xlsfinfo(fullPathName);
numSheets = length(sheetNames);

if ( strcmp(fileType, 'Microsoft Excel Spreadsheet') && numSheets >= 2 )
    
    % >>>>>>>>> DEFINE sheetNumber MANUALLY HERE <<<<<<<<<<<<<<
    
    % Multiple sheets in the Excel file!
    sheetNumber = char(sheetNames(2));
    fprintf('There are %d Excel sheets in this file. Make sure to analyze both!\n', ...
    numSheets);
    fprintf('Loading data from %s now...\n', ...
        sheetNumber) ;
    tempSpace = xlsread(fullPathName, char(sheetNumber));

elseif ( strcmp(fileType, 'Microsoft Excel Spreadsheet') && numSheets == 1 )

    % The typical / default case, with a single Sheet1 in Excel file
    fprintf('Data file is a legit Excel file. Loading Sheet1 data...\n');
    sheetNumber = 'Sheet1';    % default is Sheet1 of .xlsx file
    tempSpace = xlsread(fullPathName, char(sheetNumber));
    
else
    
    % Bail if this isn't a usable Excel file
    fprintf('That`s not the filetype you`re looking for... bailing!\n\n');

end


%% Separate vars into useful names

time_ms     = tempSpace(:,1);   %  time of sample, in milliseconds
% green       = tempSpace(:,2);   %  raw PPG waveform for green wavelength
red         = tempSpace(:,2);   %  raw PPG waveform for red wavelength
nir         = tempSpace(:,3);   %  raw PPG waveform for nearIR wavelength
%instaBPM    = tempSpace(:,5);   %  instant BPM estimate (high variance)
%avgBPM      = tempSpace(:,6);   %  time-averaged BPM estimate (more stable)
%lastBPM_ms  = tempSpace(:,7);   %  time in ms for last stable avgBPM calc
%temp        = tempSpace(:,8);   %  temperature

% define time in seconds
time_s = time_ms / 1000;

%% Do NaN cleanup, if needed, setting equal to zero

redNaNs = isnan(red);
nirNaNs = isnan(nir);

red(find(redNaNs == 1)) = 0;
nir(find(nirNaNs == 1)) = 0;

%% aggregate red and nir into a matrix

[X] = [red nir];

clear tempSpace;                 % to free up memory

%% Evaluate PPG signals for zero values; skip zero'd matrix (sensor = off)

anyGreen = exist('green');

if anyGreen == 1
    if (green(:,1) == 0 )
        fprintf('Green PPGs absent; skipping green PPG.\n');
        useGreen = 0;

    else
        fprintf('Green PPGs present; using green PPG.\n');
        useGreen = 1;
    end
end


if (red(:,1)) == 0
    fprintf('Red PPGs absent; skipping red PPG.\n');
    useRed = 0;
else
    fprintf('Red PPGs present; using red PPG.\n');
    useRed = 1;
end


if (nir(:,1)) == 0
    fprintf('NIR PPGs absent; skipping NIR PPG.\n');
    useNIR = 0;
else
    fprintf('NIR PPGs present; using NIR PPG.\n');
    useNIR = 1;
end

% TO-DO: until green, red, NIR detection works, use green baseline = 0
green(1:length(red),1) = zeros;

%% do a quick diagnostic plot of the entire data set


figureCount = 1;

figure;
plot(time_s, nir,'b.');
hold on;
plot(time_s,red,'r.');
plot(time_s,green,'g-');
hold off;
xlabel('Time, seconds');
ylabel('Raw ADC values');
tempTitle = sprintf('Figure %d: Raw PPG Time Series', ...
    figureCount);
title(tempTitle);

%% compute and display masimoDST

if newerMatlab == 1

    fprintf('This Matlab version is too new to run MasimoDST.\n');
    fprintf('Kludging together Masimo ppg and noise using NIR ppg.\n');
    ppg = nir; noise = nir;

else
    
    fprintf('\n-------------------------------------------------\n');
    fprintf('Hang on a sec... running Masimo DST stuff...\n');
    [r ppg noise] = masimoDST(X);
    fprintf('\nOkay... done!\n');
    fprintf('-------------------------------------------------\n\n');

end

%% plot raw ppg vs masimo ppg

figureCount = figureCount + 1;

figure;

fig2a = subplot(311);
plot(time_s, red, 'r-');
xlabel('Time, s');
ylabel('Raw red ADC ppg');
tempTitle = sprintf('Figure %d: Raw PPG, Masimo PPG, and Masimo Noise', ...
    figureCount);
title(tempTitle);

fig2b = subplot(312);
plot(time_s, ppg,'b-');
xlabel('Time, s');
ylabel('Masimo PPG');

fig2c = subplot(313);
plot(time_s, noise,'k-');
xlabel('Time, s');
ylabel('Masimo Noise');

linkaxes([fig2a, fig2b, fig2c], 'x');


%% use all data (= 1) , or choose a sub-sample time range (= 0)

% Choose whether to grab some or all of the raw PPG signal

all = 1;

if (all == 1)
    % newSig = nir(1:length(nir),1);      % capture the whole sig
    fprintf('Proceeding with analysis of the entire time series\n');
    newSig = red(1:length(red),1);      % capture the whole sig
else
    % NOTE: THIS IS CURRENTLY NOT FULLY FUNCTIONAL
    % TO-DO: Define lower and upper bounds in time_s
    %        then grab PPG / sensor data of interest in that range
    timeLower = 60;     % e.g.: from 60 seconds
    timeUpper = 120;    %       to 120 seconds
    newSig = nir(2300:3000,1);          % capture a part of the sig
    fprintf('Proceeding with analysis of one time series segment:\n');
    time
end

%% Divide AC and DC components of all signals... may want to do this before

%  Do a quick sanity check to make sure all signal matrix lengths are equal

if (length(red) == length(nir) && length(ppg) == length(nir) )
    fprintf('All signals are equal length matricies. Continuing!\n');
    
    sigLength = length(ppg);
    
    Ts     = mean(diff(time_s));          % average time sample
    Fs     = 1 / Ts;                      % effective sampling rate
    nFFT   = 2 ^ nextpow2(sigLength);     % next power of 2 from sig length
    nirFFT = fft(nir, nFFT) / sigLength;  % NIR FFT
    redFFT = fft(red, nFFT) / sigLength;  % Red FFT
    ppgFFT = fft(ppg, nFFT) / sigLength;  % Masimo PPG FFT
    
    f      = Fs/2 * linspace(0,1,nFFT/2+1);
    freq   = -Fs/2 + Fs/nFFT:Fs/nFFT:Fs/2;
    
    redCenter = fftshift(redFFT);
    nirCenter = fftshift(nirFFT);
    ppgCenter = fftshift(ppgFFT);
    
    redAC     = fftshift(redFFT(2:end));
    nirAC     = fftshift(nirFFT(2:end));
    ppgAC     = fftshift(ppgFFT(2:end));
    
    redDC     = redFFT(1);
    nirDC     = nirFFT(1);
    ppgDC     = ppgFFT(1);
    
    figureCount = figureCount + 1;
    
    figure;
    
    figFFT1a = subplot(311);
    plot(freq(2:end), abs(redAC),'r-');
    xlim([ 0.166666667 4.3 ]);
    ylabel('abs(red(f))');
    tempTitle = sprintf('Figure %d.1: Red AC Single-Sided Amplitude Spectrum', ...
        figureCount);
    title(tempTitle);
    
    figFFT1b = subplot(312);
    plot(freq(2:end), abs(nirAC),'b-');
    xlim([ 0.166666667 4.3 ]);
    ylabel('abs(nir(f))');
    tempTitle = sprintf('Figure %d.2: NIR AC Single-Sided Amplitude Spectrum', ...
        figureCount);
    title(tempTitle);
    
    figFFT1c = subplot(313);
    plot(freq(2:end), abs(ppgAC),'m-');
    xlim([ 0.166666667 4.3 ]);
    xlabel('Frequency (Hz)');
    ylabel('abs(ppg(f))');
    tempTitle = sprintf('Figure %d.3: Masimo AC PPG Single-Sided Amplitude Spectrum', ...
        figureCount);
    title(tempTitle);
    
    linkaxes( [figFFT1a, figFFT1b, figFFT1c], 'x');
    
else
    
    fprintf('Signal length mismatch. Check the data for problems!\n');
    
end
%% apply 2nd order butterworth filter

sampleInterval    = Fs;                     
samplingRate      = 1 / Fs;
% At ~50 Hz, cutoffValues tranlate to cutoffFreq and BPM, as follows:
% cutoffValue =  250 -> cutOffFreq = 0.2005 Hz -> BPM = 12.030
% cutoffValue =  333 -> cutOffFreq = 0.1505 Hz -> BPM = 9.0300
% cutoffValue =  500 -> cutOffFreq = 0.1003 Hz -> BPM = 6.0180
% cutoffValue = 1000 -> cutOffFreq = 0.0501 Hz -> BPM = 3.0060
% cutoffValue = 2000 -> cutOffFreq = 0.0251 Hz -> BPM = 1.5060
% cutoffValue = 3000 -> cutOffFreq = 0.0167 Hz -> BPM = 1.0020                             
% cutoffValue = 4000 -> cutOffFreq = 0.0125 Hz -> BPM = 0.7500                               
cutoffValue       = 4000;           
cutoffFrequency   = 1 / (samplingRate * cutoffValue) ;    
butterOrder       = 6;
butterType        = 'high';
% butterNyquistNormalizedCutoff = cutoffFrequency / (samplingRate / 2);
butterNyquistNormalizedCutoff = cutoffFrequency;


[b, a] = butter(butterOrder, butterNyquistNormalizedCutoff, butterType);
    
filteredPPGgreen  = filtfilt( b, a, green(:,1) );
filteredPPGred    = filtfilt( b, a, red(:,1) );
filteredPPGnir    = filtfilt( b, a, nir(:,1) );
filteredPPG       = filtfilt( b, a, ppg(:,1) ); % use the Masimo PPG too

%% Plot overlaid filtered ppgWaveforms + standalone Masimo in two panels

figureCount = figureCount + 1;

figure;

figFilt1a = subplot(211);
plot(time_s, filteredPPGred,'r-');
xlabel('Time, s');
ylabel('G/R/NIR Reflectance');
tempTitle = sprintf('Figure %d.1: Filtered G/R/NIR PPG (%d° BW w %0.3g Hz Cutoff)', ...
    figureCount, butterOrder, cutoffFrequency);
title(tempTitle);
hold on;
plot(time_s, filteredPPGnir, 'b-');
plot(time_s, filteredPPGgreen,'g-');
hold off;

figFilt1b = subplot(212);
plot(time_s, filteredPPG, 'm-');
xlabel('Time, s');
ylabel('Masimo PPG');
tempTitle = sprintf('Figure %d.1: Filtered Masimo PPG ((%d° BW w %0.3g Hz Cutoff)', ...
    figureCount, butterOrder, cutoffFrequency);
title(tempTitle);
linkaxes([figFilt1a, figFilt1b], 'x');

%% Plot the filtered ppgWaveforms in individual panels

figureCount = figureCount + 1;

figure;
fig3a = subplot(411);
plot(time_s, filteredPPGgreen,'g-');
ylabel('Green Reflectance');
ylim([-200 200]);
tempTitle = sprintf('Figure %d: Filtered PPG ((%d° BW w %0.3g Hz Cutoff)', ...
    figureCount, butterOrder, cutoffFrequency);
title(tempTitle);
fig3b = subplot(412);
plot(time_s, filteredPPGred,'r-');
ylabel('Red Reflectance');
ylim([-200 200]);

fig3c = subplot(413);
plot(time_s, filteredPPGnir, 'b-');
xlabel('Time, s');
ylabel('NIR Reflectance');
ylim([-200 200]);

fig3d = subplot(414);
plot(time_s, filteredPPG, 'm-');
xlabel('Time, s');
ylabel('Masimo PPG');
ylim([-2000 2000]);
linkaxes( [fig3a, fig3b, fig3c, fig3d], 'x');


%% Make some spectral analysis plots

welchNFFT     = 2^7;        % try segment lengths of 2^8 = 256
winSize  = hanning(welchNFFT);   % set hanning window shape
nOverlap = welchNFFT / 2;        % set 50% overlap between segments

[Pred, Fred] = periodogram(red, [], welchNFFT, Fs, 'power');
figure; 
plot(Fred,10*log10(Pred),'r-');
xlabel('Frequency, Hz'); ylabel('Power spectrum (dBW)');

[PredPower,FredPower] = pwelch(red,ones(welchNFFT,1),0,welchNFFT,Fs,'power');
figure; 
plot(FredPower,10*log10(PredPower),'r-');
xlabel('Frequency, Hz'); ylabel('Power spectrum (dBW)');

[psd_green, f_green] = pwelch(filteredPPGgreen(:,1), winSize, nOverlap, welchNFFT);
[psd_red, f_red] = pwelch(filteredPPGred(:,1), winSize, nOverlap, welchNFFT);
[psd_nir, f_nir] = pwelch(filteredPPGnir(:,1), winSize, nOverlap, welchNFFT);
[psd_ppg, f_ppg] = pwelch(ppg(:,1), winSize, nOverlap, welchNFFT);

%% Try to make some spectrograms and periodograms of pulse ox data

[S, F, T, P] = spectrogram(ppg,winSize,nOverlap,nFFT,Fs);

% make surface3d plot
figure;
colormap(parula(5));
% surface with edgecolor
% surfc(T,F,10*log10(abs(P)));
% surface with no edgecolor
surfc(T,F,10*log10(abs(P)), 'EdgeColor','none');
xlabel('Time,s'); ylabel('Frequency, Hz'); zlabel('Magnitude, dB');
axis tight;
view(-45,45);

%% make contour3d plot

figure;
colormap(parula(5));
% surfc(T,F,10*log10(abs(P)), 'EdgeColor','none');
contour3(T,F,10*log10(abs(P)));
xlabel('Time,s'); ylabel('Frequency, Hz'); 
c = colorbar;
c.Label.String = 'Intensity, dB';
axis tight;


%% Try a static spectrogram plot

colormap(parula(5));
colorLimit = [40 90] ;  % color axis limits in dB for specgram
figure; 
imagesc(T, F, 10*log10(abs(P)));
% contourf(10*log10(abs(P)));
xlabel('Time,s'); ylabel('Frequency, Hz'); 
c = colorbar;
c.Label.String = 'Intensity, dB';
axis tight;
axis xy;

%% Make some PSDs and plots of PSDs

doPSDPlots = 1;

if (doPSDPlots == 1)
    
    figureCount = figureCount + 1;
        
    figure;

    figPSD1 = subplot(411);
    loglog(f_green, psd_green, '-', 'Color', 'Green');
    xlabel ('Frequency (Hz)');
    ylabel ('PSD');
    figureCount = figureCount + 1;
    tempTitle = sprintf('Figure %d.1: Filtered Green PPG PSD', ...
        figureCount);
    title(tempTitle);


    figPSD2 = subplot(412);
    loglog(f_red, psd_red, 'r-');
    xlabel ('Frequency (Hz)');
    ylabel ('PSD');
    % draw the cutoff line on filter plot
    hold on;
    loglog( [-log10(cutoffFrequency) -log10(cutoffFrequency)], [10^2 10^8], 'k:');
    hold off;
    tempTitle = sprintf('Figure %d.2: Filtered Red PPG PSD', ...
        figureCount);
    title(tempTitle);


    figPSD3 = subplot(413);
    loglog(f_nir, psd_nir, '-', 'Color', 'blue');
    xlabel ('Frequency (Hz)');
    ylabel ('PSD');
    hold on;
    loglog( [-log10(cutoffFrequency) -log10(cutoffFrequency)], [10^2 10^8], 'k:');
    hold off;
    tempTitle = sprintf('Figure %d.3: Filtered NIR PPG PSD', ...
        figureCount);
    title(tempTitle);

    
    figPSD4 = subplot(414);
    loglog(f_ppg, psd_ppg, '-', 'Color', 'magenta');
    xlabel ('Frequency (Hz)');
    ylabel ('PSD');
    hold on;
    loglog( [-log10(cutoffFrequency) -log10(cutoffFrequency)], [10^2 10^8], 'k:');
    hold off;    
    tempTitle = sprintf('Figure %d.4: Filtered Masimo PPG PSD', ...
        figureCount);
    title(tempTitle);
    
    linkaxes([ figPSD1, figPSD2, figPSD3, figPSD4], 'x');
    
end


%% Get 1st and 2nd derivatives of the filtered PPG waveform

d1red = diff(filteredPPGred) / (samplingRate * 0.1) ; d1red(end+1,1) = 0;
d1nir = diff(filteredPPGnir) / (samplingRate * 0.1) ; d1nir(end+1,1) = 0;
% d1ppg = diff(filteredPPG)    / (samplingRate * 0.1) ; d1ppg(end+1,1) = 0;
d1ppg = diff(filteredPPG, 1);
d1ppg(end + 1, 1) = 0;


d2red = diff(d1red)/ (samplingRate * 0.1) ; d2red(end+1,1) = 0;
d2nir = diff(d1nir)/ (samplingRate * 0.1) ; d2nir(end+1,1) = 0;
% d2ppg = diff(d1ppg)/ (samplingRate * 0.1) ; d2ppg(end+1,1) = 0;
d2ppg = diff(filteredPPG, 2);
d2ppg(end + 1, 1) = 0; d2ppg(end + 1, 1) = 0; 




%% plot 1st and 2nd ppg derivatives for red

figureCount = figureCount + 1;

figure; 

fig4a = subplot(311);

plot(time_s, filteredPPGred, 'r-');
% [maxtab, mintab] = peakdet(filteredPPGred, 0.5);
%hold on; plot(mintab(:,1), mintab(:,2), 'g*');
%plot(maxtab(:,1), maxtab(:,2), 'r*');
ylabel('Red Intensity');
ylim([-200, 200]);
tempTitle = sprintf('Figure %d: Filtered Red PPG w/ 1° and 2nd Derivatives', ...
    figureCount);
title(tempTitle);

fig4b = subplot(312);
plot(time_s, d1red, 'k-');
ylabel('1° Red Intensity');

fig4c = subplot(313);
plot(time_s, d2red,'k-');
xlabel('Time, s');
ylabel('2° Red Intensity');

linkaxes( ([fig4a, fig4b, fig4c]), 'x');

%% plot 1st and 2nd ppg derivatives for NIR

figureCount = figureCount + 1;

figure; 

fig4a = subplot(311);
plot(time_s, filteredPPGnir, 'b-');
% [maxtab, mintab] = peakdet(filteredPPGnir, 0.5);
% hold on; plot(mintab(:,1), mintab(:,2), 'g*');
% plot(maxtab(:,1), maxtab(:,2), 'r*');
ylabel('NIR Intensity');
ylim([-200, 200]);

tempTitle = sprintf('Figure %d: Filtered PPG w/ 1° and 2nd Derivatives', ...
    figureCount);
title(tempTitle);

fig4b = subplot(312);
plot(time_s, d1nir, 'k-');
ylabel('1° NIR Intensity');

fig4c = subplot(313);
plot(time_s, d2nir,'k-');
xlabel('Time, ms');
ylabel('2° NIR Intensity');

linkaxes( ([fig4a, fig4b, fig4c]), 'x');

%% plot 1st and 2nd ppg derivatives for Masimo PPG

if newerMatlab == 0
    
    figureCount = figureCount + 1;

    figure; 

    fig4a = subplot(311);
    plot(time_s, filteredPPG, 'm-');
    % [maxtab, mintab] = peakdet(filteredPPG, 0.5);
    % hold on; plot(mintab(:,1), mintab(:,2), 'g*');
    % plot(maxtab(:,1), maxtab(:,2), 'r*');
    ylabel('Masimo PPG Intensity');
    ylim([-2000, 2000]);
    tempTitle = sprintf('Figure %d: Filtered PPG w/ 1° and 2nd Derivatives', ...
        figureCount);
    title(tempTitle);

    fig4b = subplot(312);
    plot(time_s, d1nir, 'k-');
    ylabel('1° NIR Intensity');

    fig4c = subplot(313);
    plot(time_s, d2nir,'k-');
    xlabel('Time, ms');
    ylabel('2° NIR Intensity');

    linkaxes( ([fig4a, fig4b, fig4c]), 'x');

end
%% do some experimental frequency domain analysis using Masimo PPG

if newerMatlab == 0

    nFFT = length(filteredPPG(:,1));
    Y    = fft(filteredPPG(:,1), nFFT); 
    F    = ((0:1/nFFT:1-1/nFFT)*Fs);

    magnitudePPGFFT = abs(Y);        % Magnitude of the FFT
    phasePPGFFT     = unwrap(angle(Y));  % Phase of the FFT

    figureCount = figureCount + 1;

    figure;
    fig9a = subplot(211);
    plot(F, magnitudePPGFFT,'m-');
    xlabel('Frequency, Hz');
    ylabel('Magnitude, dB');
    tempTitle = sprintf('Figure %d.1: Masimo PPG FFT', ...
        figureCount);
    title(tempTitle);
    fig9b = subplot(212);
    plot(F, phasePPGFFT, 'm-');
    xlabel('Frequency, Hz');
    ylabel('Phase, radians');
    tempTitle = sprintf('Figure %d.2: Masimo PPG Phase FFT', ...
        figureCount);
    title(tempTitle);

end

%% do some experimental frequency domain analysis on filtered red PPG


nFFTred = length(filteredPPGred(:,1));
Y    = fft(filteredPPGred(:,1), nFFTred); 
F    = ((0:1/nFFTred:1-1/nFFTred)*Fs);

magnitudeRedFFT = abs(Y);        % Magnitude of the FFT
phaseRedFFT     = unwrap(angle(Y));  % Phase of the FFT

figureCount = figureCount + 1;

figure;
fig9a = subplot(211);
plot(F, magnitudeRedFFT,'r-');
xlabel('Frequency, Hz');
ylabel('Magnitude, dB');
tempTitle = sprintf('Figure %d.1: Filtered Red PPG FFT', ...
    figureCount);
title(tempTitle);
fig9b = subplot(212);
plot(F, phaseRedFFT, 'r-');
xlabel('Frequency, Hz');
ylabel('Phase, radians');
tempTitle = sprintf('Figure %d.2: Filtered Red PPG Phase FFT', ...
    figureCount);
title(tempTitle);

%% do some experimental frequency domain analysis on filtered NIR PPG


nFFTnir = length(filteredPPGnir(:,1));
Y    = fft(filteredPPGnir(:,1), nFFTnir); 
F    = ((0:1/nFFTnir:1-1/nFFTnir)*Fs);

magnitudeNIRFFT = abs(Y);        % Magnitude of the FFT
phaseNIRFFT     = unwrap(angle(Y));  % Phase of the FFT

figureCount = figureCount + 1;

figure;
fig9a = subplot(211);
plot(F, magnitudeNIRFFT,'b-');
xlabel('Frequency, Hz');
ylabel('Magnitude, dB');
tempTitle = sprintf('Figure %d.1: Filtered NIR PPG FFT', ...
    figureCount);
title(tempTitle);

fig9b = subplot(212);
plot(F, phaseNIRFFT, 'b-');
xlabel('Frequency, Hz');
ylabel('Phase, radians');
tempTitle = sprintf('Figure %d.2: Filtered NIR PPG Phase FFT', ...
    figureCount);
title(tempTitle);

%% overplot red and nir FFT 

figureCount = figureCount + 1;

figure;
plot(F, magnitudeNIRFFT,'b-');
hold on;
plot(F, magnitudeRedFFT,'r-');
xlabel('Frequency, Hz'); ylabel('Magnitude, dB');
xlim([0 3]);
tempTitle = sprintf('Figure %d: Filtered Red and NIR PPG FFT', ...
    figureCount);
title(tempTitle);

%% check if new version of Matlab and do some wavelets

if (newerMatlab == 1)
    
    fprintf('Matlab wavelets supported! Wait while these are generated...\n');
    
    doCWT = 1;
    
    if doCWT == 1
        
        fprintf('\tWorking on cwt...\n');
        % do cwt for raw red and nir

        figure;
        subplot(211);
        cwt(red,Fs);
        caxis([0 200]);
        title('cwt(red)');
        subplot(212);
        caxis([0 500]);
        cwt(nir,Fs);
        
        title('cwt(nir)');

        % do cwt for filtered red and nir

        figure;
        subplot(211);
        cwt(filteredPPGred,Fs);
        caxis([0 200]);
        title('cwt(filteredRed)');
        subplot(212);
        cwt(filteredPPGnir,Fs);
        caxis([0 500]);
        title('cwt(filteredNIR)');

    end
    
    % do wsst for for raw red and nir
    
    doWSST = 1;
    
    if doWSST == 1

        fprintf('\tWorking on wsst...\n');
        
        % do wsst for for raw red and nir
        
        figure;
        subplot(211);
        wsst(red,Fs);
        ylim([0 5]);
        caxis([0 10]);
        title('wsst(red)');
        subplot(212);
        wsst(nir,Fs);
        ylim([0 5]);
        caxis([0 10]);
        title('wsst(nir)');

        % do wsst for for filtered red and nir

        figure;
        figWSSTred = subplot(211);
        wsst(filteredPPGred,Fs);
        ylim([0 5]);
        caxis([0 5]);
        title('wsst(filteredRed)');
        figWSSTnir = subplot(212);
        wsst(filteredPPGnir,Fs);
        ylim([0 5]);
        caxis([0 10]);
        title('wsst(filteredNIR)');

        % do some manual wsst stuff for plotting contours of least
        % penalized WSST results
        
        doManualWSST = 1;
        
        if doManualWSST == 1
            
            % build wsst for signals of interest
            
            [wsstRed, wsstRedF]         = wsst(red,Fs);
            [wsstNIR, wsstNIRF]         = wsst(nir,Fs);
            [wsstFiltRed, wsstFiltRedF] = wsst(filteredPPGred,Fs);
            [wsstFiltNIR, wsstFiltNIRF] = wsst(filteredPPGnir,Fs);
        
            % build wsst ridges for signals of interest, using 'NumRidges'
            % set to 2, to get primary signal and 1st harmonic
            % and run a loop (for now) to explore the effect of 'penalty'
            
            doMultiPenaltyTest = 0;
            
            if doMultiPenaltyTest == 0
                
                penalty = 10;
                
                wsstRedRidge = wsstridge(wsstRed, penalty, wsstRedF, ...
                    'NumRidges',2);
                wsstFiltRedRidge = wsstridge(wsstFiltRed, penalty, ...
                    wsstFiltRedF, 'NumRidges',2);
                wsstNIRRidge = wsstridge(wsstNIR, penalty, wsstNIRF, ...
                    'NumRidges',2);
                wsstFiltNIRRidge = wsstridge(wsstFiltNIR, penalty, ...
                    wsstFiltNIRF, 'NumRidges',2);

                redPulseF = wsstFiltRedRidge(:,1);
                nirPulseF = wsstFiltNIRRidge(:,1);
                redPulse1stH = wsstFiltRedRidge(:,2);
                nirPulse1stH = wsstFiltNIRRidge(:,2);

                redPulseBPM   = 60 * redPulseF;
                nirPulseBPM   = 60 * nirPulseF;
                redPulseBPMh1 = 60 * redPulse1stH;
                nirPulseBPMh1 = 60 * nirPulse1stH;

                % do contour plots with ridge overlays

                figure;
                contour(time_s,wsstRedF,abs(wsstRed));
                hold on;
                plot(time_s, wsstRedRidge, 'r.');
                hold off;
                xlabel('Time, s'); ylabel('Frequency, Hz');
                grid on;
                ylim([0 4]);
                caxis([0 5]);
                title('Synchrosqueezed Transform - Raw Red PPG');

                figure;
                contour(time_s,wsstFiltRedF,abs(wsstFiltRed));
                hold on;
                plot(time_s, wsstFiltRedRidge, 'r.');
                hold off;
                xlabel('Time, s'); ylabel('Frequency, Hz');
                grid on;
                ylim([0 4]);
                caxis([0 5]);
                title('Synchrosqueezed Transform for Filtered Red PPG');            

                figure;
                contour(time_s,wsstNIRF,abs(wsstNIR));
                hold on;
                plot(time_s, wsstNIRRidge, 'b.');
                hold off;
                xlabel('Time, s'); ylabel('Frequency, Hz');
                grid on;
                ylim([0 4]);
                caxis([0 10]);
                title('Synchrosqueezed Transform for Raw NIR PPG');  

                figure;
                contour(time_s,wsstFiltNIRF,abs(wsstFiltNIR));
                hold on;
                plot(time_s, wsstFiltNIRRidge, 'b.');
                hold off;
                xlabel('Time, s'); ylabel('Frequency, Hz');
                grid on;
                ylim([0 4]);
                caxis([0 10]);
                title('Synchrosqueezed Transform for Filtered NIR PPG'); 

                % plot 1st wsstridge for filtered red and nir as -
                % plot 2nd wsstridge for filtered red and nir as -.

                figure;
                plot(time_s, wsstFiltRedRidge(:,1), 'r-');
                hold on; 
                plot(time_s, wsstFiltRedRidge(:,2), 'r--');
                plot(time_s, wsstFiltNIRRidge(:,1), 'b-');
                plot(time_s, wsstFiltNIRRidge(:,2), 'b--');
                hold off;
                xlabel('Time, s'); ylabel('Frequency, Hz');
                ylim([0 2.5]);
                title('Synchrosqueezed Ridges - Raw and Filtered PPGs');
                legend('Red Pulse','NIR Pulse','Red 1st Harmonic', ...
                    'NIR 1st Harmonic');

                figure;
                plot(time_s,nirPulseBPM,'b-');
                hold on;
                plot(time_s,redPulseBPM,'r-');
                hold off;
                xlabel('Time, s'); ylabel('Estimated Heart Rate, BPM');
                title('Synchrosqueeze Ridges - Estimated Heart Rate');
                grid on;

                figure;
                plot(time_s,nirPulseBPMh1,'b-');
                hold on;
                plot(time_s,redPulseBPMh1,'r-');
                hold off;
                xlabel('Time, s'); ylabel('Estimated 1H Heart Rate, BPM');
                title('Synchrosqueezed Ridges - Estimated HR 1st Harmonic');
                grid on;
                
            else
                
                for penalty = 6:4:11

                    wsstRedRidge = wsstridge(wsstRed, penalty, wsstRedF, ...
                        'NumRidges',2);
                    wsstFiltRedRidge = wsstridge(wsstFiltRed, penalty, ...
                        wsstFiltRedF, 'NumRidges',2);
                    wsstNIRRidge = wsstridge(wsstNIR, penalty, wsstNIRF, ...
                        'NumRidges',2);
                    wsstFiltNIRRidge = wsstridge(wsstFiltNIR, penalty, ...
                        wsstFiltNIRF, 'NumRidges',2);

                    redPulseF = wsstFiltRedRidge(:,1);
                    nirPulseF = wsstFiltNIRRidge(:,1);
                    redPulse1stH = wsstFiltRedRidge(:,2);
                    nirPulse1stH = wsstFiltNIRRidge(:,2);

                    redPulseBPM   = 60 * redPulseF;
                    nirPulseBPM   = 60 * nirPulseF;
                    redPulseBPMh1 = 60 * redPulse1stH;
                    nirPulseBPMh1 = 60 * nirPulse1stH;

                    % do contour plots with ridge overlays

                    figure;
                    contour(time_s,wsstRedF,abs(wsstRed));
                    hold on;
                    plot(time_s, wsstRedRidge, 'r.');
                    hold off;
                    xlabel('Time, s'); ylabel('Frequency, Hz');
                    grid on;
                    ylim([0 4]);
                    caxis([0 5]);
                    title('Synchrosqueezed Transform - Raw Red PPG');

                    figure;
                    contour(time_s,wsstFiltRedF,abs(wsstFiltRed));
                    hold on;
                    plot(time_s, wsstFiltRedRidge, 'r.');
                    hold off;
                    xlabel('Time, s'); ylabel('Frequency, Hz');
                    grid on;
                    ylim([0 4]);
                    caxis([0 5]);
                    title('Synchrosqueezed Transform for Filtered Red PPG');            

                    figure;
                    contour(time_s,wsstNIRF,abs(wsstNIR));
                    hold on;
                    plot(time_s, wsstNIRRidge, 'b.');
                    hold off;
                    xlabel('Time, s'); ylabel('Frequency, Hz');
                    grid on;
                    ylim([0 4]);
                    caxis([0 10]);
                    title('Synchrosqueezed Transform for Raw NIR PPG');  

                    figure;
                    contour(time_s,wsstFiltNIRF,abs(wsstFiltNIR));
                    hold on;
                    plot(time_s, wsstFiltNIRRidge, 'b.');
                    hold off;
                    xlabel('Time, s'); ylabel('Frequency, Hz');
                    grid on;
                    ylim([0 4]);
                    caxis([0 10]);
                    title('Synchrosqueezed Transform for Filtered NIR PPG'); 

                    % plot 1st wsstridge for filtered red and nir as -
                    % plot 2nd wsstridge for filtered red and nir as -.

                    figure;
                    plot(time_s, wsstFiltRedRidge(:,1), 'r-');
                    hold on; 
                    plot(time_s, wsstFiltRedRidge(:,2), 'r--');
                    plot(time_s, wsstFiltNIRRidge(:,1), 'b-');
                    plot(time_s, wsstFiltNIRRidge(:,2), 'b--');
                    hold off;
                    xlabel('Time, s'); ylabel('Frequency, Hz');
                    ylim([0 2.5]);
                    title('Synchrosqueezed Ridges - Raw and Filtered PPGs');
                    legend('Red Pulse','NIR Pulse','Red 1st Harmonic', ...
                        'NIR 1st Harmonic');

                    figure;
                    plot(time_s,nirPulseBPM,'b-');
                    hold on;
                    plot(time_s,redPulseBPM,'r-');
                    hold off;
                    xlabel('Time, s'); ylabel('Estimated Heart Rate, BPM');
                    title('Synchrosqueeze Ridges - Estimated Heart Rate');
                    grid on;

                    figure;
                    plot(time_s,nirPulseBPMh1,'b-');
                    hold on;
                    plot(time_s,redPulseBPMh1,'r-');
                    hold off;
                    xlabel('Time, s'); ylabel('Estimated 1H Heart Rate, BPM');
                    title('Synchrosqueezed Ridges - Estimated HR 1st Harmonic');
                    grid on;
            end

            end
            
        end
        
    end
else
    
    fprintf('This version of Matlab is too old for wavelets!\n');

end


%% try dddtree cleanup on a segment of unfiltered red ppg


% next to be an even multiple of 2 ^ J, e.g.: 8192, 16384
newRed  = red(150:150+2047);         
dt      = 1/fix(Fs);
t       = 0:dt:(length(newRed)*dt)-dt;
J       = 6;
dtcplx1 = dddtree('cplxdt', newRed, J, 'dtf3');
QQQ     = dtcplx1.cfs{2}(:,1,1) + 1i * dtcplx1.cfs{2}(:,1,2);
QQQQ    = dtcplx1.cfs{3}(:,1,1) + 1i * dtcplx1.cfs{3}(:,1,2);

figure;
figD3Tree1b = subplot(211);
stem(real(QQQ),'r-');
figD3Tree1c = subplot(212);
stem(real(QQQQ),'r-');
% figD3Tree1a = subplot(311);
% plot(t, newRed,'r-');

% linkaxes([ figD3Tree1a, figD3Tree1b, figD3Tree1c], 'x');
linkaxes([ figD3Tree1b, figD3Tree1c], 'x');