%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tomasz Zieliński (CR)
% „Starting Digital Signal Processing in Telecommunication Engineering.
% A Laboratory-Based Course.”
% Springer 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Xf, f] = pwelch(x, win,  Nover, f   , fs, freqrange )
%                   pwelch(x, win,  Nover, Nfft, fs, 'onesided')
%                   pwelch(x, Nwin, Nover, f   , fs, 'twosided') 
%                   pwelch(x, Nwin, Nover, Nfft, fs, 'centered') 
% INPUT:
% x     - signal to analyze
% win   - window function used
% Nwin  - window length
% Nover - window overlap
% Nfft  - FFT length
% f     - frequencies of interest
% fs    - sampling frequency
% OUTPUT:
% Xf = WelchPowerSpectralDensity(f) [V^2 (dB/Hz)]

% One should specify at least 5 first input arguments
if    ( nargin<6 && ~isreal(x) ) freqrange = 'centered';
elseif( nargin<6 &&  isreal(x) ) freqrange = 'onesided';
else
end

N = length(x); x = reshape(x,1,N);                     % signal length, to horizontal
if( length(win)==1 ) Nwin=win; win=hamming(Nwin)';     % default window 
else  Nwin = length(win); win = reshape(win,1,Nwin);   % window length, to horizontal
end

    if( strcmp( freqrange,'centered') )
    else                         
  end
  
if( length(f)==1 )
    Nfft = f;
    if( strcmp( freqrange,'onesided') ) what=1; k=1:Nfft/2+1; offs=0;  end % onesided  for real
    if( strcmp( freqrange,'twosided') ) what=2; k=1:Nfft; offs=0;      end % twosided  for complex    else
    if( strcmp( freqrange,'centered') ) what=3; k=1:Nfft; offs=Nfft/2; end
    f = fs/Nfft*(k-1-offs);
else
    Nfft = length(f); k=1:Nfft;                       % always two-sides 
    if( f(1)==0 )   what=2; f=fs/Nfft*(k-1);          % not-centered
    else            what=3; f=fs/Nfft*(k-Nfft/2-1);   % centered
    end    
end
Nstep = Nwin-Nover;                                   % window shift
Nspec = floor((N-Nwin)/Nstep)+1;                      % number of FFT spectra 
Xf = zeros(1,length(k)); work = zeros(1,length(k));   % PSD initialization 

for m = 1 : Nspec                                     % processing loop
    bx = x( 1+(m-1)*Nstep : Nwin+(m-1)*Nstep );       % signal fragment 
    bx = bx .* win;                                   % windowing
    X = fft( bx, Nfft )/sum(win);                     % FFT with scaling
    if(what==1 || what==2) work(k) = X(k);            % PSD: whole or half
    else                   work(k) = fftshift(X(k));  % PSD: whole centered 
    end                                               %
    Xf = Xf + abs(work).^2;                           % power accumulation
end                                                   % loop end
Xf = (1/Nspec)*Xf/fs;                                 % PSD normalization

semilogy(f,Xf); grid; title('Welch estimation of PSD'); xlabel('f [Hz]'); ylabel('V^2 / Hz');
