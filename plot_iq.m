%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_iq(waveform, sample)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_iq(waveform, sample)
    %plot real and imag part of first n signal samples
    figure;
    subplot(211);
    n = 1:sample;
    plot(n,real(waveform(n)), 'bo-'); title('Signal real part');
    subplot(212);
    plot(n,imag(waveform(n)), 'rx-'); title('Signal imaginary part');    
end