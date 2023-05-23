clear all
close all

N = 256;
test_run = 100000;

%sequences = transpose(net(sobolset(7),(N)));
sequences = VonDerCorput_6(N);
%sequences = transpose(rand(N,7));
%sequences = niederreiter2_generate(20, N, 2, 31);

%sequences = LFSR_Bulk(N)/N;

for iter = 1:1:test_run
    %i = randi(N+1,1)-1;
    %j = randi(N+1,1)-1;
    %m = randi(N+1,1)-1;
    %n = randi(N+1,1)-1;

    i = rand(1,1);
    j = rand(1,1);
    m = rand(1,1);
    n = rand(1,1);

    expected = min(1, ((i)*(j)) + ((m)*(n)))

%     X1_Sobol_BS = number_source(i, N, sequences(1, :));
%     X2_Sobol_BS = number_source(j, N, sequences(7, :));
%     X3_Sobol_BS = number_source(m, N, sequences(2, :));
%     X4_Sobol_BS = number_source(n, N, sequences(5, :));

    X1_Sobol_BS = number_source(i, N, sequences(1, :));
    X2_Sobol_BS = number_source(j, N, sequences(2, :));
    X3_Sobol_BS = number_source(m, N, sequences(3, :));
    X4_Sobol_BS = number_source(n, N, sequences(4, :));

    ANDControl1 = ((X1_Sobol_BS)&(X2_Sobol_BS));
    ANDControl2 = ((X3_Sobol_BS)&(X4_Sobol_BS));
    
    AND3_minfinder = or(ANDControl1, ANDControl2);

    obtained = sum(AND3_minfinder)/N;

    err(iter) = abs(expected - obtained);

end

sum(err)/test_run

%sequences = (Hammersley(N,6))
