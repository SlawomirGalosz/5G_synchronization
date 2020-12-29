%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tomasz Zieliński (CR)
% „Starting Digital Signal Processing in Telecommunication Engineering.
% A Laboratory-Based Course.”
% Springer 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Xtf, t, f] = spectrogram(x, win,  Nover, f   , fs, freqaxis )
%                       spectrogram(x, win,  Nover, Nfft, fs, 'xaxis')
%                       spectrogram(x, Nwin, Nover, f   , fs, 'yaxis') 
%                       spectrogram(x, Nwin, Nover, Nfft, fs, 'yaxis') 
% INPUT:
% x     - signal to analyze
% win   - window function weights
% Nwin  - window length
% Nover - window overlap
% Nfft  - FFT length
% f - frequencies of interest
% fs - sampling frequency
% freqaxis - 'xaxis','yaxis': position of the frequency axis
% OUTPUT:
% X1 = STFT(t,f) [V (dB)]

% One should specify at least 5 first input arguments
if( nargin < 6 ) freqaxis = 'xaxis'; end
    
N = length(x); x = reshape(x,1,N);                    % signal length, to horizontal
if( length(win)==1 ) Nwin=win; win=hamming(Nwin)';    % default window 
else  Nwin = length(win); win = reshape(win,1,Nwin);  % window length, to horizontal
end

if( length(f)==1 )
    Nfft = f;
    if( isreal(x) ) what=1; k=1:Nfft/2+1;             % one-side  for real
    else            what=2; k=1:Nfft;                 % two-sides for complex
    end
    f = fs/Nfft*(k-1);
else
    Nfft = length(f); k=1:Nfft;                       % always two-sides 
    if( f(1)==0 )   what=2; f=fs/Nfft*(k-1);          % not-centered
    else            what=3; f=fs/Nfft*(k-Nfft/2-1);   % centered
    end    
end
Nstep = Nwin-Nover;                                   % window shift
Nspec = floor((N-Nwin)/Nstep)+1;                      % number of STFT spectra 
Xtf = zeros(Nspec,length(k));                         % STFT initialization 

for m = 1 : Nspec                                     % processing loop
    bx = x( 1+(m-1)*Nstep : Nwin+(m-1)*Nstep );       % signal fragment 
    bx = bx .* win;                                   % windowing
    X = fft( bx, Nfft )/sum(win);                     % FFT with scaling
    if(what==1 || what==2) Xtf(m,k) = X(k);           % FFT: whole or half
    else                   Xtf(m,k) = fftshift(X(k)); % FFT: whole centered 
    end                                               %
end                                                   % loop end
Xtf = 10*log10( real(Xtf.*conj(Xtf)) );               % to decibels

dt=1/fs; t=(Nwin/2+1/2)*dt + Nstep*dt*(0:Nspec-1);    % time
if( strcmp( freqaxis,'xaxis') )
    h=imagesc(f,t,Xtf);        % amplitude spectrum matrix as an image
    xlabel('Frequency (Hz)'); ylabel('Time (s)'); set( gca,'YDir','normal');
else
    imagesc(t,f,Xtf.');      % amplitude spectrum matrix as an image
    xlabel('Time (s)'); ylabel('Frequency (Hz)'); set( gca,'YDir','normal');

end    
cb=colorbar('location','EastOutside'); set( get(cb,'Ylabel'),'String','V (dB)');

