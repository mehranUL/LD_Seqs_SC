% How to generate bit-stream based on faure sequence
%Maximum dimension: 40  ----> https://www.extremeoptimization.com/Documentation/Mathematics/Random-Numbers/Quasi-Random-Sequences.aspx
ff = faure(D,image_row_size*image_column_size,10); %last parameter is base. base>2 are better
ff = ff(:,2:end);
ff = ff';
PP_ff = ones(D,image_row_size*image_column_size);
for i = 1:image_row_size*image_column_size    
        for z = 1:D
            if 0.5 <= ff(z,i)
                PP_ff(z,i) = -1;
            end            
        end    
end
PP_ff = PP_ff';


%Reference: Dimitri Shvorob (2023). Generate a Faure sequence (https://www.mathworks.com/matlabcentral/fileexchange/15373-generate-a-faure-sequence), MATLAB Central File Exchange. Retrieved April 25, 2023.

function[s] = faure(k,d,b)
% FAURE     Faure sequence, elements 0,..,k
% INPUTS  : k - maximum sequence index, non-negative integer
%           d - sequence dimension, positive integer
%           b - sequence base, integer exceeding 1
% OUTPUTS : s - d*(k+1) array, with s(:,i) storing element (i+1)
%               of base-b d-dimensional Faure sequence
% EXAMPLE : faure(0,1,2) = 0
%           faure(1,1,2) = [0 0.5]'
%           faure(2,1,2) = [0 0.5 0.25]'
% AUTHOR  : Dimitri Shvorob, dimitri.shvorob@vanderbilt.edu, 6/20/07
if ~(isint(k) && k >= 0)
   error('Input argument "k" must be a non-negative integer')
end
if ~(isint(d) && d > 0)
   error('Input argument "d" must be a positive integer')
end
if ~(isint(b) && b > 1)
   error('Input argument "b" must be an integer greater than 1')
end
s = zeros(d,k+1);
K = k;
D = d;
parfor k = 0:K                   
    a = basexpflip(k,b);
    J = length(a);
    L = J - 1;                
    y = zeros(J,1);
    g = b.^(1:J)';
    for d = 1:D          
        for j = 1:J      
            S = 0;
            for l = 0:L 
                c = comb(l,j-1);
                h = (d-1)^(l-j+1); 
                if isinf(h)
                   h = 0; 
                end   
                S = S + c*h*a(l+1);
            end
            y(j) = mod(S,b);
        end
        s(d,k+1) = sum(y./g);
    end
end   
function[a] = basexpflip(k,b) % reversed base-b expansion of non-negative integer k
if k 
   j = fix(log(k)/log(b)) + 1;
   a = zeros(1,j);
   q = b^(j-1);
   for i = 1:j
       a(i) = floor(k/q);
       k = k - q*a(i);
       q = q/b;
   end
   a = fliplr(a);
else
   a = 0;
end
function[c] = comb(n,k)       % number of combinations, C(n,k)
if n < k
   c = 0;
else
   c = nchoosek(n,k);
end   
function[i] = isint(x)        % check if integer 
i = (x == floor(x));