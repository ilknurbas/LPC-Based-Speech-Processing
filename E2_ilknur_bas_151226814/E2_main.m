clear all

curdir = fileparts(which('E2_main.m'));
addpath([curdir '/aux_scripts/']);


%% Step 0) Waveform reading and resampling
% Read waveform 251-136532-0016.flac
% and resample to 16 kHz (if not already; always ensure correct sample rate!)
%
% Required functions:   audioread()
%                       resample() 

[x,fs] = audioread([curdir '/data/251-136532-0016.flac']);

% Read textGrid annotation of the signal

TG = readTextGrid([curdir '/data/251-136532-0016.TextGrid']);

% Resample to 16 kHz (if not 16 kHz originally)

if(fs ~= 16000)
   x = resample(x,16000,fs);
   fs = 16000;
end

%% Step 1) Spectrogram, energy, and zero-crossing rate
%--------------------------------------------------------------------------
% Calculate logarithmic FFT magnitude spectrogram, log-energy, 
% and zero-crossing rate using a sliding 20-ms window with 10-ms steps. 
% Use Hamming windowing for FFT and rectangular
% windowing for ZCR and energy. 
% 
%--------------------------------------------------------------------------
%
% Helpful functions:    fft(), hamming(), log10(), sum(), diff(), sign()
%                       figure(), plot(), subplot()
% 
%
% Definitions: 
%       D1a:    Zero-crossing rate (ZCR) is the number of times that the
%               signal amplitude crosses zero value (i.e., changes sign) 
%               per time unit. 
%               Usually: zcr = N_crossings / (number of samples-1).
%
% Hints:             
%       H1a: for zero-crossing rate, combining diff(), sign(), and
%       sum() functions in a correct way will do the trick!
%
%       H1b: Voiced speech spectrum typically has approx. 60 dB dynamic
%       range. 
%
%       H1c: set(gca,'YDir','normal'); can be used to flip y-axis of an
%       image (spectrogram) so that frequencies increase from bottom to
%       top.
%
%       H1d: One-sided magnitude spectrum calculated from an even number of 
%       N samples has N/2+1 unique values from 0 to pi, and the rest are
%       symmetric to those. 
%
%       H1e: Frequency-axis of the spectrum runs from 0 to fs/2 Hz (Nyquist
%       frequency), consisting of the N/2+1 values mentioned above.


% Define window length and window step size
 
% wl = ?    % window length in samples
% ws = ?    % step size in samples
wl = 0.02 * fs;
ws = 0.01 * fs;

% ww = ?    % windowing function?
ww = hamming(wl);

% Variables for storing FFT, energy and ZCR for each window position.
FFT_y = zeros(round(length(x)/ws)-2,wl/2+1);
E = zeros(round(length(x)/ws)-2,1);
zcr = zeros(round(length(x)/ws)-2,1);
% variables for comparison
E_ham = zeros(round(length(x)/ws)-2,1);
zcr_ham = zeros(round(length(x)/ws)-2,1);
%%
% Loop over signal with pre-defined steps
c = 1;
for winpos = 1:ws:length(x)-wl+1
    
    % Get signal window
    % y =  ?
    y = x(winpos:winpos+wl-1); % frame
    y_rect = y.* rectwin(wl) ; % ones(wl,1) or rectwin(wl) same thing
    % y_win for hamming 
    y_win = y.* ww;

    % Calculate FFT for the Hamming-windowed frame.
    % Remember to trim the mirrored (redundant) part from the spectrum 
    % (hint H3d).
    % FFT_y(c,:) = ?
    temp = 20*log10(abs(fft(y_win)));
    FFT_y(c,:) = temp(1:wl/2+1);

    % Calculate energy of rectangular window
    % E(c) = ?
    E(c) = 10 * log10(sum(y_rect .^2));
    % for comparison
    E_ham(c) = 10 * log10(sum(y_win .^2));
        
    % Calculate zero-crossing of rectangular window
    % zcr(c) = ?
    % sign finds the sign of each sample
    % diff calculates differences between adjacent elements
    n_crossings = sum(abs(diff(sign(y_rect)))~=0);
    zcr(c) = n_crossings/(length(y_rect)-1);
    % for comparison
    n_crossings_ham = sum(abs(diff(sign(y_win)))~=0);
    zcr_ham(c) = n_crossings_ham/(length(y_win)-1);
    
    c = c+1;  
end

if zcr_ham == zcr
    display("zcr_ham == zcr")
end
if E_ham == E
    display("E_ham == E")
end
%%
figure(2);clf;
% plot magnitude spectrogram and appropriate axes
% t = ?   % time axis
num_frames = floor((length(x) - wl) / ws) + 1;
t = (0:num_frames-1)*ws/fs;

% f = ?   % frequency axis
f = linspace(0, fs/2, length(FFT_y)); % trimmed version
 
subplot(3,1,1); 
axis xy; imagesc(t, f, FFT_y'); 
set(gca,'YDir','normal'); % flip y-axis to have frequencies in ascending order
xlabel('time (s)'); ylabel('frequency (Hz)');title('Logarithmic Magnitude Spectrogram');

subplot(3,1,2);
% plot energy as a function of time with appropriate axes
plot(t, E); xlabel('time (s)'); ylabel('energy (dB)'); title('Energy vs. Time');
% can be commented
xlim([0, ceil(length(x)/ws)*ws/fs]); % To set the last value of x as the end of the x-axis

subplot(3,1,3);
% plot ZCR as a functon of time with appropriate axes
plot(t, zcr); xlabel('time(s)'); ylabel('zcr'); title('Zero-Crossing Rate');
% can be commented
xlim([0, ceil(length(x)/ws)*ws/fs]); % To set the last value of x as the end of the x-axis 

%%
% for comparison
figure(3);clf;
subplot(2,1,1);
plot(t, E); xlabel('time (s)'); ylabel('energy (dB)'); title('Energy vs. Time using rectangular window');
xlim([0, ceil(length(x)/ws)*ws/fs]); % To set the last value of x as the end of the x-axis 
subplot(2,1,2);
plot(t, E_ham); xlabel('time (s)'); ylabel('energy (dB)'); title('Energy vs. Time using Hamming window');
xlim([0, ceil(length(x)/ws)*ws/fs]); % To set the last value of x as the end of the x-axis 

%% Step 2) Linear prediction coefficients
%--------------------------------------------------------------------------
% Implement LPC estimation using autocorrelation method. 
% Calculate LPC coefficients using 20 ms hamming window and 10 ms step from
% pre-emphasized speech signal, using model order 20.
% Create a plot that has three panels:
% the top panel shows original magnitude spectrogram of the signal
% the middle pane shows magnitude spectrogram of the LPC filters
% the bottom one shows magnitude spectrogram of the LPC residual
%--------------------------------------------------------------------------
%
% Helpful functions:    fft(), hamming(), log10(), sum(), diff(), sign()
%                       filter()
%
% Hints:             
%       H2a: For LPC, you will always need 2 poles to decsribe one formant. 
%           As a rule of thumb: A good LPC order is usually fs/1000 + N, 
%           where N is a small constant (e.g., 2 or 4). With 8 kHz signals,
%           that would give, e.g., order of 10, allowing to describe 4
%           formants + general level and tilt of the spectrum.
%       H2b: You can sanity check your getlpc() implementation aganst lpc()
%           built-in to MATLAB. Note that lpc() uses Levinson Durbin
%           recursion for solving the LPC normal equation, whereas your
%           implementation does that through matrix inverse. Therefore yo
%           cannot use lpc() implementation as it is, but the outputs
%           should be identical (to some numerical precision).%
%       H2c: LPC residual is the difference ("prediction error") between 
%           the original speech signal and the signal captured by the LPC 
%           filter. For each window, you can simply filter the signal with 
%           the estimated LPC analysis filter (FIR filter) to obtain the
%           residual.
%       H2d: Since LPC filter focuses on modeling the energetic resonances 
%           of the vocal tract, the residual should contain what is left of
%           the signal after the resonances have been removed. In contrast,
%           LPC filter should capture the key aspects of the spectral
%           envelope but cannot describe spectral details unless LPC order
%           is very high.
%            


% Define window length and window step size
 
% wl = ?    % window length in samples
% ws = ?    % step size in samples
wl = 0.02 * fs;
ws = 0.01 * fs;

% ww = ?    % windowing function?
ww = hamming(wl);

%lpc_order = ?; % Choose an LPC order suitable for 16 kHz data.
lpc_order = 20;

% Apply 1st order pre-emphasis to the signal

% x_emph = ? 
h = [1 -0.95];
x_emph = filter(h, 1, x);

% Create matrices for storing calculated values

LPC_y = zeros(round(length(x)/ws)-2,lpc_order+1); % LPC coefficients 
g_y = zeros(round(length(x)/ws)-2,1); % gain coefficients
RES_y = zeros(round(length(x)/ws)-2,wl); % residuals
FFT_LPC_y = zeros(round(length(x)/ws)-2,wl/2+1); % spectra of LPC
FFT_RES_y = zeros(round(length(x)/ws)-2,wl/2+1); % spectra of residual

c = 1;
for winpos = 1:ws:length(x_emph)-wl+1
     
    % Extract signal window
    % y = ?
    y = x_emph(winpos:winpos+wl-1); % frame
    y = y .* ww;
    
    % Calculate LPC for the window using your implementation of getlpc()  
    [LPC_y(c,:),g_y(c),~] = getlpc(y,lpc_order);
    
    % LPC_y values are corrrect
    % display(LPC_y(c,:))
    % real_a = lpc(y, lpc_order)

    % Calculate prediction error signal aka. residual for the window using 
    % the LPC analysis filter (Eq. 4.5 in instructions). 
    
    % RES_y(c,:) = ? % for each frame and for each sample in that frame 
    RES_y(c,:) = filter(LPC_y(c,:),1, y);
 
    % Calculate LPC filter magnitude spectrum for the window using FFT (Eq.
    % 4.16)
    % FFT_LPC_y(c,:) = ?
    a_k = abs(fft(LPC_y(c,1:end), wl));
    temp = 20*log10(g_y(c)) - 20*log10(a_k);
    FFT_LPC_y(c,:) = temp(1:wl/2+1);
        
    % Calculate residual magnitude spectrum for the window using FFT
    % FFT_RES_y(c,:) = ?
    temp = -20*log10(abs(fft(RES_y(c,:))));
    FFT_RES_y(c,:) = temp(1:wl/2+1);
    

    c = c+1;
end



%%
num_frames = round(length(x)/ws)-2; 
t = (0:num_frames-1) * ws/fs;
f = linspace(0, fs/2, length(FFT_RES_y));

figure(4);clf;
subplot(3,1,1);
% Plot original speech magnitude spectrogram (calculated in step 3)
axis xy; imagesc(t, f, FFT_y'); 
set(gca,'YDir','normal'); % flip y-axis to have frequencies in ascending order
xlabel('time (s)'); ylabel('frequency (Hz)');title('Original Logarithmic Magnitude Spectrogram');
colorbar

subplot(3,1,2);
% Plot LPC synthesis filter magnitude spectrogram
imagesc(t, f, FFT_LPC_y'); 
set(gca,'YDir','normal'); % flip y-axis to have frequencies in ascending order
xlabel('time (s)'); ylabel('frequency (Hz)');title('LPC synthesis filter Logarithmic magnitude spectrogram');
colorbar

subplot(3,1,3);
% Plot LPC residual magnitude spectrogram 
imagesc(t, f, FFT_RES_y'); 
set(gca,'YDir','normal'); % flip y-axis to have frequencies in ascending order
xlabel('time (s)'); ylabel('frequency (Hz)');title('LPC residual Logarithmic magnitude spectrogram');
colorbar

%% Step 3) Signal re-synthesis with LPC
%--------------------------------------------------------------------------
% Synthesize the original utterance back from the LPC coefficients using
% A) original residual as the excitation, and
% B) impulse train as the excitation. For this, use a meaningful
%    fundamental frequency (F0) of your choice. 
%
%       (optional): you may combine impulse-based and noise-based
%       excitation to simulate voiced and unvoiced frames. You can use zero
%       crossing rate (ZCR) extracted in step 1) to roughly determine 
%       voiced and unvoiced segments.
%
% Save the synthesized waveforms as synth_impulse.wav and 
% synth_residual.wav
%--------------------------------------------------------------------------
% Helpful functions:    filter(), audiowrite()
%
% Hints:             
%       H3a: LPC synthesis requires an all-pole model, which is an IIR
%       filter (i.e., not a FIR). 
%       H3b: Remember to remove signal pre-emphasis afterwards!
%       H3c: Since LPC analysis separates the signal into two components
%           without loss of information, a proper reconstruction of the original
%           signal from the same two components should sound similar to the
%           original signal.
 

% f0 = ?    % set f0 of the signal
f0 = 100 ;

for excitation_type = {'impulse','residual'}
 
    % Call signal re-synthesis function
    [x_synth,~] = lpcResynthesis(x,LPC_y,RES_y,g_y,ws,wl,fs,f0,excitation_type);
    
    % Save to file
    audiowrite(sprintf([curdir '/synth_%s.wav'],char(excitation_type)),x_synth,fs);
    
end






