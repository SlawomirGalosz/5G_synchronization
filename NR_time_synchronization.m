%downlink synchronization script
clear all; close all;
signal_choice = 1;    %different signals
do_disturbance = 0;   %add disturbance to signal
do_repeat = 0;        %repeat signal 5 times
do_CFO = 1;           %do signal corrections

if (signal_choice == 1)
    load('Signals/NR_DL_2159.91MHz_10MHz.mat'); %ncellid = 440
    nfft = 1024;
    sample_rate = 15360000;
    scs = 15*1e3;
    nb_rb = 52;
    CPlen1 = 80;
    CPlen2 = 72;
    sym_len1 = nfft + CPlen1;
    sym_len2 = nfft + CPlen2;
    threshold = 60; %
    ssb_offset = -48; %SSB position related to center subcarrier
%For generated signals we have one frame lasting 10 ms
%Since CP autocorrelation is calculated for 4 frames, option do_repeat should be set to 1
elseif (signal_choice == 2)
    load('Signals/5g_downlink_10MHzBW_15kHzSCS.mat'); %ncellid = 23
    waveform = waveStruct1.waveform(:, 1).';
    nfft = waveStruct1.config.Nfft;
    sample_rate = waveStruct1.Fs;
    scs = waveStruct1.config.SubcarrierSpacing*1e3;
    nb_rb = 52;
    nr_sample_rate = sample_rate;
    CPlen1 = 80;
    CPlen2 = 72;
    sym_len1 = nfft + CPlen1;
    sym_len2 = nfft + CPlen2;
    threshold = 20;
    ssb_offset = 0;
elseif (signal_choice == 3)
    load('Signals/5g_downlink_10MHzBW_30kHzSCS.mat'); %ncellid = 23
    waveform = waveStruct.waveform.';
    nfft = waveStruct.config.waveform.Nfft;
    sample_rate = waveStruct.Fs;
    scs = waveStruct.config.waveform.SubcarrierSpacing*1e3;
    nb_rb = 24;
    nr_sample_rate = sample_rate;
    CPlen1 = 40;
    CPlen2 = 36;
    sym_len1 = nfft + CPlen1;
    sym_len2 = nfft + CPlen2;
    threshold = 20;
    ssb_offset = 0;
elseif (signal_choice == 4)
    load('Signals/5g_downlink_100MHzBW_30kHzSCS.mat'); %ncellid = 512
    waveform = waveStruct.waveform(:, 1).';
    nfft = 4096;
    sample_rate = waveStruct.Fs;
    scs = 30*1e3;
    nb_rb = 273;
    nr_sample_rate = sample_rate;
    CPlen1 = 320;
    CPlen2 = 288;
    sym_len1 = nfft + CPlen1;
    sym_len2 = nfft + CPlen2;
    threshold = 20;
    ssb_offset = 0;
end

N = length(waveform);

plot_iq(waveform, 100);
plot_freq(waveform, sample_rate); title('Captured signal spectrum');
figure; pwelch(waveform,2048,2048-1024,2048,sample_rate,'centered');
figure; f = sample_rate/2048*(-1024:1023); spectrogram(waveform,2048,1024,f,sample_rate); pause

%signal repeating - preferred for shorter signals
if (do_repeat)
    waveform = [waveform waveform waveform waveform waveform];
    N = length(waveform);
end


if (do_disturbance)
    % Parameter values 
    frc_cfo  = 2000;           %fractional offset [Hz]
    int_cfo  = 0;              %integer offset    [sample]
    noise_pw = -25;            %sigal to noise    [dB]
    fading_model = 'TDLC300';  %for available models refer to delay_profiles.m
    
    %add delay spread
    [delay, gain] = delay_profiles(fading_model, sample_rate);
    fading = zeros(12, N);
    for n = 1:length(delay)
        fading(n, :) = gain(n) * [waveform(end-delay(n)+1:end) waveform(1:end-delay(n))];
    end
    waveform = sum(fading) / sum(gain);
    
    %add cfo 
    waveform = waveform .* exp(1j*2*pi * (int_cfo*scs + frc_cfo) / sample_rate *(0:N-1));
    %add white noise
    if (noise_pw ~= 0)
        %be aware that different signals could have different powers
        noise = sqrt(1 / (10^(-noise_pw/10)*2)) * (randn(1, N) + 1j*randn(1, N));
        waveform = waveform + noise;
    end

    figure; pwelch(waveform,2048,2048-1024,2048,sample_rate,'centered'); 
    title('Welch estimation of PSD wth all interferences added');
end


%Signal autocorelation
ccCP = [];
for n = 1:sample_rate/25 %4 ofdm frames
    n1 = n:n+CPlen2-1; 
    n2 = n1 + nfft;
    a1 = abs(sum( waveform(n1) .* conj(waveform(n2))));
    b1 = abs(sum( waveform(n1) .* conj(waveform(n1)) + sum(waveform(n2) .* conj(waveform(n2)))));
    ccCP(n) = (a1*a1) / (b1*b1);
end
figure; plot(ccCP(1:length(ccCP)/4));
title('Signal autocorrelation using CP'); xlabel('sample index n');
pause

%franctional CFO estimation and correction
if (do_CFO)
    npeak = 50;
    
    %find a few max peaks
    cp_pos = zeros(1, npeak);
    for n = 1:npeak
        [aa, cp_pos(n)] = max(ccCP);
        %don't take indexes from beyond array length
        start = max(1, cp_pos(n)-sym_len2);
        stop = min(length(ccCP), cp_pos(n)+sym_len2);
        %clean this symbol to ensure that we don't take it again
        ccCP(start:stop) = zeros(1, stop-start+1);
    end
    
    cp_pos = sort(cp_pos);
    %calculate CFO for choosen symbols
    for n = 1:npeak 
        n1 = cp_pos(n):cp_pos(n)+CPlen2-1;
        n2 = n1 + nfft;
        indv_CFO(n) = sample_rate/nfft * mean(angle(conj(waveform(n1)) .* waveform(n2)) / (2*pi));
    end

    figure; plot(cp_pos(1:npeak), indv_CFO(1:npeak), 'o-');
    title('Estimated fractional CFO'); xlabel('sample index n'); ylabel('[Hz]');
    
    %correct signal with mean of different CFO values
    frc_CFO = mean(indv_CFO)
    waveform = waveform .* exp(-1j*2*pi * frc_CFO/sample_rate *(0:N-1));
    pause
end


%integer CFO estimation and correction
if (do_CFO)
    search_space = 10;
    for n = -search_space:search_space
        [test_indxs, test_nid] = pss_decoding(waveform(1:end), nfft, threshold, n+ssb_offset);
        if (~isempty(test_indxs))
            %get NID2 and PSS position indexes
            int_CFO = n
            NID2 = test_nid
            pss_indxs = test_indxs;
            %signal correction
            waveform = waveform .* exp(-1j*2*pi * int_CFO/nfft *(0:N-1));
            break
        end
    end
else
    [pss_indxs, NID2] = pss_decoding(waveform(1:end), nfft, threshold, 0+ssb_offset);
    NID2
end
pss_indxs = pss_indxs - sym_len2 + 1;
pause

for pss = 1:length(pss_indxs)
    pss_pos = pss_indxs(pss);
    %second CFO estimation using PSS
    if (do_CFO)
        n1 = pss_pos:pss_pos+CPlen2-1;
        n2 = n1 + nfft;
        pss_CFO = sample_rate/nfft * mean(angle(conj(waveform(n1)) .* waveform(n2)) / (2*pi));

        %correction only for block we want to extract
        block_size = (pss_pos:pss_pos+4*sym_len2-1);
        waveform(block_size) = waveform(block_size) .* exp(-1j*2*pi * pss_CFO/sample_rate *(0:length(block_size)-1));
    end

    SS_block = zeros(240, 4);
    for n = 1:4
        pos = pss_pos + sym_len2*(n-1);
        samples = waveform(pos+CPlen2:pos+CPlen2+nfft-1);
        time_samples = fftshift(fft(samples)) / sqrt(nfft);
        SS_block(:, n) = time_samples(nfft/2-119+ssb_offset:nfft/2+120+ssb_offset);
    end

    %searching for NID1
    NID1 = sss_decoding(SS_block(57:183, 3).', NID2);

    %PCI calculation
    ncellid = 3 * NID1 + NID2
    pause

    %find correct issb for extracted SS block
    issb = pbchdmrs_decode(ncellid, SS_block);
    pause

    %find indicates for PBCH
    pbch_pos = PBCHposition(ncellid);

    %channel estimation 
    ref_dmrs = generate_dmrs(ncellid, issb);
    rx_dmrs = SS_block(PBCHDMRSposition(ncellid)).';
    hest = rx_dmrs .* conj(ref_dmrs);
    %extend channel estimation to PBCH samples
    hest = hest + [0;0;0];
    pbch_hest = hest(:);

    %PBCH correction
    pbch_sym = SS_block(pbch_pos);
    if(do_CFO)
        pbch_eq = pbch_sym ./ pbch_hest;
    else
        pbch_eq = pbch_sym;
    end
    %show PBCH constellation
    figure; plot(pbch_eq, 'o');
    title('PBCH constellation'); xlabel('In-Phase'); ylabel('Quadrature');
    pause
    
    %MIB decoding using 5G toolbox functions
    v = mod(issb, 4);
    pbchBits = nrPBCHDecode(pbch_eq, ncellid, v, 1e-2);
    
    polarListLength = 8;
    [~, crcBCH, trblk, sfn4lsb, nHalfFrame, msbidxoffset] = ...
        nrBCHDecode(pbchBits, polarListLength, 4, ncellid);
    
    % Display the BCH CRC and ssb index
    disp([' BCH CRC: ' num2str(crcBCH)]);
    disp([' SSB index: ' num2str(v)]);

    k_SSB = msbidxoffset * 16;
    commonSCSs = [15 30];
    
    % Create a structure of MIB fields from the decoded MIB bits. The BCH
    % transport block 'trblk' is the RRC message BCCH-BCH-Message, consisting
    % of a leading 0 bit then 23 bits corresponding to the MIB
    mib.NFrame = bi2de([trblk(2:7); sfn4lsb] .','left-msb');
    mib.SubcarrierSpacingCommon = commonSCSs(trblk(8) + 1);
    mib.k_SSB = k_SSB + bi2de(trblk(9:12).','left-msb');
    mib.DMRSTypeAPosition = 2 + trblk(13);
    mib.PDCCHConfigSIB1 = bi2de(trblk(14:21).','left-msb');
    mib.CellBarred = trblk(22);
    mib.IntraFreqReselection = trblk(23);

    % Display the MIB structure
    disp(' BCH/MIB Content:')
    disp(mib);
    pause
end