% Generating bit-stream related to Hammersley sequence.
HH = Hammersley(1024,144);
HH = HH';
PP_hammersley = ones(1024,144);
for i = 1:144 
    for z = 1:1024
         if 0.43 <= HH(z,i) % 0.43 is a good threshold according to several trials.
              PP_hammersley(z,i) = -1;
         end            
    end    
end 
PP_hammersley = PP_hammersley';


%Reference: mansour torabi (2023). MATLAB Hammersley Sampling for Design of Experiments DOE (https://github.com/Mansourt/MATLAB_Hammersley_Sampling_for_Design_of_Experiments_DOE/releases/tag/v1.0), GitHub. Retrieved April 23, 2023.

function H = Hammersley(P,N)
%% 
% Author: Mansour Torabi
% Email:  smtoraabi@ymail.com
%%
% P: Number of samples in each dimensions
% N: Number of dimensions
% Hammersley Sequence:
% Ref: Wong, Tien-Tsin, Wai-Shing Luk, and Pheng-Ann Heng. "Sampling with 
% Hammersley and Halton points." Journal of graphics tools - 1997
H(1,:)=(0:P-1)/P;
PP = primes(1000);
%PP = PP(randperm(numel(PP)));
for i = 2:N  
  for j = 0:P-1
      J = j;
      H(i,j+1)=0;
      C=1/PP(i-1);
      while J > 0
          R0 = rem(J,PP(i-1));
          H(i,j+1) = H(i,j+1) + R0*C;
          
          C = C / PP(i-1);
          J = floor(J/PP(i-1));
      end
  end
end