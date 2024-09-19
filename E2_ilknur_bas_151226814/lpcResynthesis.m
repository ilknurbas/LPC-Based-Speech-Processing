function [x_synth,impulse_train] = lpcResynthesis(orig_signal,LPC_y,RES_y,g_y,ws,wl,fs,f0,excitation_type,zcr)
% function [x_synth,impulse_train] = lpcResynthesis(orig_signal,LPC_y,RES_y,g_y,ws,wl,fs,f0,excitation_type,zcr)
%
% Synthesizes the original utterance back from the LPC coefficients using
%     A) original residual as the excitation, or
%     B) impulse train as the excitation.
%
% Inputs:
%       orig_signal:       the original signal-to-be-synthesized
%       LPC_y:             LPC inverse filter coefficients (from getlpc.m)
%       RES_y:             LPC residual
%       g_y:               gain (from getlpc.m)
%       ws:                step size in samples
%       wl:                window length in samples
%       fs:                sampling rate
%       f0:                fundamental frequency
%       excitation_type:   the type of the excitation ('impulse' / 'residual')
%       zcr:               (Optional) zero-crossing rate
%
% Outputs:
%       x_synth:           synthesized signal
%       impulse_train:     impulse train excitation signal

if nargin == 10
    use_voicing = 1;
else
    use_voicing = 0;
end


% Create impulse train excitation that corresponds to the original input
% signal in duration.

% impulse_train = ?
T0 = 1/f0; 
% impulse_train = zeros(size(orig_signal));
impulse_train = zeros(size(orig_signal,1), 1);
%impulse_train(1:round(fs*T0):end) = 1;
impulse_train(1:(fs*T0):end) = 1;

wpos = 1;    
x_synth = zeros(size(orig_signal)); % empty signal to synthesize

% (optional) Extract voicing from zero-crossing rates
if use_voicing
    % voicing = ?
end

% Go through each extracted frame and create signal window by window
for c = 1:size(LPC_y,1)

    a = LPC_y(c,:); % LPC coefficients for the current frame

    % Create corresponding excitation signal
    if(strcmp(char(excitation_type),'residual'))
        excitation = RES_y(c,:);
        %excitation = reshape( excitation, [320,1]);
        excitation = excitation';
        g = 1;

    elseif(strcmp(char(excitation_type),'impulse'))            
        g = g_y(c);             % gain for window
        % excitation = ?        % excitation for window
        excitation = impulse_train(wpos:wpos+wl-1,:);
         

        % (optional) create a mixture of impulses and noise for voiced
        % and voiceless sounds based on "voicing" estimated from ZCR.
        if use_voicing
            % excitation = ?
        end
    end

    % Synthesize waveform for the frame using the LPC synthesis filter and 
    % excitation signal. Hint: use filter() with appropriate use of
    % gain, filter coefficients, and excitation. 

    % y_synth = filter(?,?,?)
    y_synth = zeros(wl, 1);
    y_synth = filter(g, a, excitation);

    % Overlap+add add the current signal window to the full signal to
    % synthesize

    size(x_synth(wpos:wpos+wl-1)); %320     1
    size(y_synth); %     1 320
    % x_synth(?) = x_synth(?)+y_synth
    x_synth(wpos:wpos+wl-1) = x_synth(wpos:wpos+wl-1)+y_synth;    


    wpos = wpos+ws; 
end

% Remove pre-emphasis from x_synth    
% x_synth = filter(?,?,x_synth)
% x_emph = filter(h, 1, x);
h = [1 -0.95];
x_synth = filter(1, h, x_synth);

% rescale to avoid clipping in .wav format ([-1,1])
x_synth = x_synth./max(abs(x_synth)); 
    
    
