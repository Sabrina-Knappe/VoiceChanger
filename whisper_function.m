function [outputArg1,outputArg2] = whisper_function(wavfile, increment)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% VX_whisper.m [DAFXbook, 2nd ed., chapter 7]
%===== This program makes the whisperization of a sound,
%===== by randomizing the phase
%----- user data -----
s_win = 512; % analysis window length [samples]
n1 = s_win/8; % analysis step [samples]
n2 = n1; % synthesis step [samples]
[DAFx_in,FS] = audioread(wavfile);
%----- initialize windows, arrays, etc -----
w1 = hann(s_win, 'periodic'); % analysis window
w2 = w1; % synthesis window
L = length(DAFx_in);
thing = zeros(s_win, 1);
DAFx_in = [zeros(s_win, 2); DAFx_in; zeros(s_win-mod(L,n1),2)] / max(abs(DAFx_in));
DAFx_out = zeros(length(DAFx_in),1);
tic
%UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
pin = 0;
pout = 0;
pend = length(DAFx_in) - s_win;
while pin<pend
    grain = DAFx_in(pin+1:pin+s_win).* w1;
    %===========================================
    f = fft(fftshift(grain));
    r = abs(f);
    phi = increment*pi*rand(s_win,1); %add a 2*
    ft = (r.* exp(i*phi));
    grain = fftshift(real(ifft(ft))).*w2;
    % ===========================================
    DAFx_out(pout+1:pout+s_win) = DAFx_out(pout+1:pout+s_win) + grain;
    pin = pin + n1;
    pout = pout + n2;
end
%UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
toc
%----- listening and saving the output -----
% DAFx_in = DAFx_in(s_win+1:s_win+L);
DAFx_out = DAFx_out(s_win+1:s_win+L) / max(abs(DAFx_out));
%soundsc(DAFx_out, FS);
audiowrite('whisper.wav', DAFx_out, FS);
outputArg1 = DAFx_out;
outputArg2 = FS;
end

