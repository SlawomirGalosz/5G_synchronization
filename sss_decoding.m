%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sss_decoding(sequence)
% Input:
%   sequence - sss sequence from waveform  
%   nid2 - nid2 from pss sequence
% Output:
%   NID1 - signal correct NID1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NID1 = sss_decoding(sequence, NID2)
    %correlation between extracted sss sequence and every possible 
    for n = 0:335
        gen_seq = sss_generate(n, NID2);
        all_seq(n+1) = sqrt(abs(sum(sequence .* conj(gen_seq(1:end))) .^2));
    end
    figure; stem(0:335, all_seq, 'o');
    title('Detected SSS sequence'); xlabel('Sequence number NID1');
    
    [aa, NID1] = max(all_seq);
    NID1 = NID1-1
end