function Hd = getBPFilter
%GETFILTER Returns a discrete-time filter System object.

% MATLAB Code
% Generated by MATLAB(R) 9.9 and DSP System Toolbox 9.11.
% Generated on: 11-Mar-2021 13:39:21

Fstop1 = 6;   % First Stopband Frequency
Fpass1 = 8;   % First Passband Frequency
Fpass2 = 30;  % Second Passband Frequency
Fstop2 = 35;  % Second Stopband Frequency
Astop1 = 60;   % First Stopband Attenuation (dB)
Apass  = 1;    % Passband Ripple (dB)
Astop2 = 80;   % Second Stopband Attenuation (dB)
Fs     = 800;  % Sampling Frequency

h = fdesign.bandpass('fst1,fp1,fp2,fst2,ast1,ap,ast2', Fstop1, Fpass1, ...
    Fpass2, Fstop2, Astop1, Apass, Astop2, Fs);

Hd = design(h, 'butter', ...
    'MatchExactly', 'stopband', ...
    'SystemObject', true);


