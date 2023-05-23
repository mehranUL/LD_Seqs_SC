% This sequence is generated based on irrational numbers.

% Define the irrational number
alpha = (sqrt(5) - 1) / 2; % golden ratio
%alpha = pi;
%alpha = sqrt(2) - 1; % Silver ratio
%alpha = exp(1);

weyl = mod((1:1024)*alpha, 1);
ww = zeros(1024,144);
for i = 1:144
    ww(:,i) = weyl(randperm(1024));
end

PP_weyl = ones(1024,144);

for i = 1:144   
        for z = 1:1024
            if 0.5 <= ww(z,i)
                PP_weyl(z,i) = -1;
            end            
        end    
end
PP_weyl = PP_weyl';




