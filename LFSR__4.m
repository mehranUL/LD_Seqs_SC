% Please visit: https://ieeexplore.ieee.org/document/8049760
% and https://en.wikipedia.org/wiki/Linear-feedback_shift_register

function seed_scalar = LFSR__4(seed, N)

%--------------------------------------------------------------------------
% N = 256
% initial seed: [true false true true true false false false]
% taps: 8 1 5 6
if N == 256
    
    for i = 1:1:256 % (2^3-1) Cycle
        seed_scalar(i) = b2d(flip(seed));
        
        insertion = xor(xor(xor(seed(8), seed(1)), seed(5)), seed(6));
        
        seed(8) = seed(7); %after the first cycle each LFSR value is obtained; I still keep to name it as "seed"
        seed(7) = seed(6);
        seed(6) = seed(5);
        seed(5) = seed(4);
        seed(4) = seed(3);
        seed(3) = seed(2);
        seed(2) = seed(1);
        seed(1) = insertion;
        
    end
end

end



