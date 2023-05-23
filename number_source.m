function X_BS = number_source(X_quant, N, sequences)

X_BS = zeros(1, N);

for a = 1:1:N
    if sequences(a) <= X_quant
        X_BS(a) = 1;
    end
end

end