% Generating bit-streams related to Niederreiter sequence.
[ nd1, ~ ] = niederreiter2_generate ( 20, 1024, 2, 27 ); % First argument: Max # of dimensions. Second argument: # of desired points. Third argument: integer key. Fourth argument: # of bits the default was 31.
nd1 = nd1';
PP_nd1 = ones(D,20);
for i = 1:20    
        for z = 1:D
            if 8 <= nd1(z,i)
                PP_nd1(z,i) = -1;
            end            
        end    
end
PP_nd1 = PP_nd1';