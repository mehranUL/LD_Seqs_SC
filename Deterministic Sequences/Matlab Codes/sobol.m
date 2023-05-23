% Creating MATLAB built-in Sobol sequence 
%D is the number of Sobol points in each dimension.
%maximum number of dimensions are 1111.

D = 1024;
sobol_seq = net(sobolset(1111), D);

%Generating the related bit-stream sequence based on the Sobol values
% Here we consider 144 bit-streams of length 1024 or 1K.
PP_sobol = ones(D,144);

for i = 1:144   
        for z = 1:D
            if 0.5 <= sobol_seq(z,i)
                PP_sobol(z,i) = -1;
            end            
        end    
end