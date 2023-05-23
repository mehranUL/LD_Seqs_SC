%Max dimension is 40 ---->https://www.extremeoptimization.com/Documentation/Mathematics/Random-Numbers/Quasi-Random-Sequences.aspx
p = haltonset(1024,'Skip',1e3,'Leap',1e2);
p = scramble(p,'RR2');
X = net(p,144);
X = X';

%Generating bit-stream related to Halton sequence.
PP_halton = ones(1024,144);
for i = 1:144    
        for z = 1:1024
            if 0.499 <= X(z,i)
                PP_halton(z,i) = -1;
            end            
        end    
end
PP_halton = PP_halton';