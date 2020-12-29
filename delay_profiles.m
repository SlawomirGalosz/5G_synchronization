%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [delay, gain] = delay_profiles(name)
% Input: 
%     name - name of selected delay profile relate to TS 38.104 G.2.1.1 or TR 125.943 5
%     sample_rate - sampling rate of signal
% Output:
%     delay - additional delay [sample]
%     gain - normalized relative power [0-1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [delay, gain] = delay_profiles(name, sample_rate)

    if (strcmp(name, 'TDLA30'))
        delay = [0 10 15 20 25 50 65 75 105 135 150 290];
        gain = [-15.5 0 -5.1 -5.1 -9.6 -8.2 -13.1 -11.5 -11.0 -16.2 -16.6 -26.2];
    elseif (strcmp(name, 'TDLB100'))
        delay = [0 10 20 30 35 45 55 120 170 245 330 480];
        gain = [0 -2.2 -0.6 -0.6 -0.3 -1.2 -5.9 -2.2 -0.8 -6.3 -7.5 -7.1];
    elseif (strcmp(name, 'TDLC300'))
        delay = [0 65 70 190 195 200 240 325 520 1045 1510 2595];
        gain = [-6.9 0 -7.7 -2.5 -2.4 -9.9 -8.0 -6.6 -7.1 -13.0 -14.2 -16.0];
    elseif (strcmp(name, 'no delay'))
        delay = [0];
        gain = [0];
    %these below have stronger impact over processing, especially 
    %'typical urban', this one need lower threshold for pss detection
    elseif (strcmp(name, 'typical urban'))
        delay = [0 217 512 514 517 674 882 1230 1287 1311 1349 1533 1535 1622 1818 1836 1884 1943 2048 2140];
        gain  = [-5.7 -7.6 -10.1 -10.2 -10.2 -11.5 -13.4 -16.3 -16.9 -17.1 -17.4 -19.0 -19.0 -19.8 -21.5 -21.6 -22.1 -22.6 -23.5 -24.3];
    elseif (strcmp(name, 'rural area'))
        delay = [0 42 101 129 149 245 312 410 469 528];
        gain  = [-5.2 -6.4 -8.4 -9.3 -10.0 -13.1 -15.3 -18.5 -20.4 -22.4];
    elseif (strcmp(name, 'hilly terrain'))
        delay = [0 356 441 528 546 609 625 842 916 941 15000 16172 16492 16876 16882 16978 17615 17827 17849 18016];
        gain  = [-3.6 -8.9 -10.2 -11.5 -11.8 -12.7 -13.0 -16.2 -17.3 -17.7 -17.6 -22.7 -24.1 -25.8 -25.8 -26.2 -29.0 -29.9 -30.0 -30.7];
    else
       error('Model not defined') 
    end
    
    delay = round(delay * 1e-9 * sample_rate);
    gain = 10 .^ (gain/10);
    
    if (length(delay) ~= length(gain))
        error('wrong matrix size')
    end

end