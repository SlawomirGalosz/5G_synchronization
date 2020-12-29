%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tomasz Zieliński (CR)
% „Starting Digital Signal Processing in Telecommunication Engineering.
% A Laboratory-Based Course.”
% Springer 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_freq(x, fs)
% One long DFT spectrum calculation and plotting (in logarithmic scale (dB)).
  N = length(x);                           % signal length
  X = abs( fftshift( fft(x) ) ) / N;       % FFT calculation, centering, magnitude
  f = (fs/N * (-N/2 : N/2-1) ) / 1e+6;     % frequency axis in MHz  
  figure;
  plot(f, 20 * log10(X), 'b'); xlim( [f(1), f(end)] ); grid;
  xlabel('Frequency [MHz]'); ylabel('Magnitude [dB]');
end