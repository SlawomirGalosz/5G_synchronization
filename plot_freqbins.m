function plot_freqbins(x, Nfft)
% Mean Nfft-point DFT spectrum calculation and plotting (in linear scale).
% Signal length should be equal to a multiple of Nfft. 
  N = length( x );
  k = -Nfft/2 : Nfft/2-1;
  x = reshape( x, Nfft, N/Nfft );
  X = fftshift( fft(x) / Nfft );
  X = mean( abs(X.') );
  figure; stem( k, X, 'b' ); grid; xlim([k(1),k(end)]);
  %figure; plot( k, X, 'b.-' ); grid; xlim([k(1),k(end)]);
  grid;
  xlabel('Magnitude');
  xlabel('FFT frequency bin');
end