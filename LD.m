image_row_size = 12;
image_column_size = 12;
D = 1024;

R2_ws = load("R2_1.mat");
R2 = R2_ws.z;
%PP_R2 = ones(D,image_row_size*image_column_size);
PP_R2 = ones(D,2);

for i = 1:2   
        for z = 1:D
            if 0.5 <= R2(z,i)
                PP_R2(z,i) = -1;
            end            
        end    
end

%Weyl Sequence Start
alpha = (sqrt(5) - 1) / 2;
weyl = mod((1:D)*alpha, 1);
ww = zeros(D,image_row_size*image_column_size);
for i = 1:image_row_size*image_column_size
    ww(:,i) = weyl(randperm(D));
end

PP_weyl = ones(D,image_row_size*image_column_size);

for i = 1:image_row_size*image_column_size   
        for z = 1:D
            if 0.5 <= ww(z,i)
                PP_weyl(z,i) = -1;
            end            
        end    
end
PP_weyl = PP_weyl';
sum_weyl = zeros(1,image_row_size*image_column_size);
for t = 1:image_row_size*image_column_size
    sum_weyl(t) = sum(PP_weyl(t,:));
end
sim_weyl = zeros(image_row_size*image_column_size,image_row_size*image_column_size);
for m = 1:image_row_size*image_column_size
    for n = 1:image_row_size*image_column_size
        sim_weyl(m,n) = cosine_sim(PP_weyl(m,:),PP_weyl(n,:));
    end
end
%Weyl Sequence End
%--------------------------------------------------------------------------
%Kasami Sequence Start
K_Kasami = kasami(log2(D));
K_Kasami = K_Kasami(:,1:image_row_size*image_column_size);
K_Kasami = K_Kasami';
sum_kasami = zeros(1,image_row_size*image_column_size);
for t = 1:image_row_size*image_column_size
    sum_kasami(t) = sum(K_Kasami(t,:));
end
sim_kasami = zeros(image_row_size*image_column_size,image_row_size*image_column_size);
for m = 1:image_row_size*image_column_size
    for n = 1:image_row_size*image_column_size
        sim_kasami(m,n) = cosine_sim(K_Kasami(m,:),K_Kasami(n,:));
    end
end
%Kasami Sequence Start
%--------------------------------------------------------------------------
%Sobol Sequence Start

% ws8k1 = load('sobol_bl_optimized_8k_92.mat','aa8k1');
% aa8k = ws8k1.aa8k1;

% N_sobol = 1:1111;   %Vector of sobol sequence indices
% ws = load('sobol_pairs_mul_xnor1k.mat','x2_1k');% Loads the matrix of MAE using xor operator for sobol sequences
% x2_1k = ws.x2_1k;
% a1k = find(x2_1k(1,:) ~= 0); %Find worst case sobol sequence indices
% dd = setdiff(N_sobol,a1k);
% aa = [1,dd];

ws16k1 = load('sobol_bl_optimized_16k_91.mat','aa8k1');
aa16k = ws16k1.aa8k1;
sobol_seq1 = net(sobolset(1111), D);
sobol_seq_new = sobol_seq1(:,aa16k);
PP_sobol = ones(D,image_row_size*image_column_size);

for i = 1:image_row_size*image_column_size   
        for z = 1:D
            if 0.5 <= sobol_seq_new(z,i)
                PP_sobol(z,i) = -1;
            end            
        end    
end
PP_sobol = PP_sobol';
sum_sobol = zeros(1,image_row_size*image_column_size);
for t = 1:image_row_size*image_column_size
    sum_sobol(t) = sum(PP_sobol(t,:));
end
sim_sobol = zeros(image_row_size*image_column_size,image_row_size*image_column_size);
for m = 1:image_row_size*image_column_size
    for n = 1:image_row_size*image_column_size
        sim_sobol(m,n) = cosine_sim(PP_sobol(m,:),PP_sobol(n,:));
    end
end
%Sobol Sequence End
%--------------------------------------------------------------------------
%Latin Hypercube Sequence Start
PP_Y = ones(D,2*image_row_size*image_column_size);
Y = lhsdesign(D,2*image_row_size*image_column_size);
for i = 1:2*image_row_size*image_column_size   
        for z = 1:D
            if 0.5 <= Y(z,i)
                PP_Y(z,i) = -1;
            end            
        end    
end
PP_Y = PP_Y';
sum_latin = zeros(1,2*image_row_size*image_column_size);
for t = 1:2*image_row_size*image_column_size
    sum_latin(t) = sum(PP_Y(t,:));
end
sim_latin = zeros(image_row_size*image_column_size,image_row_size*image_column_size);
for m = 1:image_row_size*image_column_size
    for n = 1:image_row_size*image_column_size
        sim_latin(m,n) = cosine_sim(PP_Y(m,:),PP_Y(n,:));
    end
end
%Latin Hypercube Sequence End
%--------------------------------------------------------------------------
%Half_LFSR Sequence Start
N = log2(D);
gg = logical(randi(2, [1 N]) - 1);
%
rnd1 = zeros(D,image_row_size*image_column_size);

for i = 1:image_row_size*image_column_size  
    rnd1(1:D/2,i) = lfsr_value512(D/2,i);
    rnd1((D/2+1):D,i) = 1 - rnd1(1:D/2,i);
end

TT = ones(D,image_row_size*image_column_size);

for i = 1:image_row_size*image_column_size
    for z = 1:D
        if 0.5 <= rnd1(z,i)
            TT(z,i) = -1;
        end
    end
end
TT = TT';
sum_lfsr = zeros(1,image_row_size*image_column_size);
for t = 1:image_row_size*image_column_size
    sum_lfsr(t) = sum(TT(t,:));
end
sim_halflfsr = zeros(image_row_size*image_column_size,image_row_size*image_column_size);
for m = 1:image_row_size*image_column_size
    for n = 1:image_row_size*image_column_size
        sim_halflfsr(m,n) = cosine_sim(TT(m,:),TT(n,:));
    end
end
%Half_LFSR Sequence End
%--------------------------------------------------------------------------
%Hadamard Sequence Start
h = hadamard(2048);
h = h(2:1025,2:145);
h = h';
sum_hadamard = zeros(1,image_row_size*image_column_size);
for t = 1:image_row_size*image_column_size
    sum_hadamard(t) = sum(h(t,:));
end
sim_hadamard = zeros(image_row_size*image_column_size,image_row_size*image_column_size);
for m = 1:image_row_size*image_column_size
    for n = 1:image_row_size*image_column_size
        sim_hadamard(m,n) = cosine_sim(h(m,:),h(n,:));
    end
end
%Hadamard Sequence End
%--------------------------------------------------------------------------
%Gold Code Sequence Start
G = gold_code();
G = G(:,1:image_row_size*image_column_size);
G = G';
sum_gold = zeros(1,image_row_size*image_column_size);
for t = 1:image_row_size*image_column_size
    sum_gold(t) = sum(G(t,:));
end
sim_gold = zeros(image_row_size*image_column_size,image_row_size*image_column_size);
for m = 1:image_row_size*image_column_size
    for n = 1:image_row_size*image_column_size
        sim_gold(m,n) = cosine_sim(G(m,:),G(n,:));
    end
end
%Gold Code Sequence End
%--------------------------------------------------------------------------
%Faure Sequence Start
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
sum_faure = zeros(1,image_row_size*image_column_size);
for t = 1:image_row_size*image_column_size
    sum_faure(t) = sum(PP_ff(t,:));
end
sim_faure = zeros(image_row_size*image_column_size,image_row_size*image_column_size);
for m = 1:image_row_size*image_column_size
    for n = 1:image_row_size*image_column_size
        sim_faure(m,n) = cosine_sim(PP_ff(m,:),PP_ff(n,:));
    end
end
%Faure Sequence End
%--------------------------------------------------------------------------
%Halton Sequence Start
p = haltonset(1024,'Skip',1e3,'Leap',1e2);
p = scramble(p,'RR2');
X = net(p,image_row_size*image_column_size);
X = X';
PP_halton = ones(D,image_row_size*image_column_size);
for i = 1:image_row_size*image_column_size    
        for z = 1:D
            if 0.499 <= X(z,i)
                PP_halton(z,i) = -1;
            end            
        end    
end
PP_halton = PP_halton';
sum_halton = zeros(1,image_row_size*image_column_size);
for t = 1:image_row_size*image_column_size
    sum_halton(t) = sum(PP_halton(t,:));
end
sim_halton = zeros(image_row_size*image_column_size,image_row_size*image_column_size);
for m = 1:image_row_size*image_column_size
    for n = 1:image_row_size*image_column_size
        sim_halton(m,n) = cosine_sim(PP_halton(m,:),PP_halton(n,:));
    end
end
%Halton Sequence End
%--------------------------------------------------------------------------
%Hammersley Sequence Start
HH = Hammersley(1024,144);
HH = HH';
PP_hammersley = ones(D,image_row_size*image_column_size);
for i = 1:image_row_size*image_column_size 
    for z = 1:D
         if 0.43 <= HH(z,i)
              PP_hammersley(z,i) = -1;
         end            
    end    
end 
PP_hammersley = PP_hammersley';
sum_hammersley = zeros(1,image_row_size*image_column_size);
for t = 1:image_row_size*image_column_size
    sum_hammersley(t) = sum(PP_hammersley(t,:));
end
sim_hammersley = zeros(image_row_size*image_column_size,image_row_size*image_column_size);
for m = 1:image_row_size*image_column_size
    for n = 1:image_row_size*image_column_size
        sim_hammersley(m,n) = cosine_sim(PP_hammersley(m,:),PP_hammersley(n,:));
    end
end
%Hammersley Sequence End

%deneme = [sum_weyl, sum_faure];

%--------------------------------------------------------------------------
%ZadoffChu Seq Start
a3 = zeros(1,16383);
for t = 2:16383
    if (mod(t,3) ~= 0) && (mod(t,43) ~= 0) && (mod(t,127) ~= 0)
        if mod(16383,t) ~= 0
            a3(t) = t;
        end
    end
end
b3 = find(a3 ~= 0);

re = zeros(16383,numel(b3));
for zd = 1:numel(b3)
    re(:,zd) = real(zadoffChuSeq(b3(zd),16383));
end
re(16384,:) = -0.01;
re = re(3000:4023,:); %for MNIST 7000-8023  2000-11999
%re = abs(re);

PP_zadoff = ones(D,image_row_size*image_column_size);
for i = 1:image_row_size*image_column_size 
    for z = 1:D
         if -0.01 <= re(z,i)
              PP_zadoff(z,i) = -1;
         end            
    end    
end 
PP_zadoff = PP_zadoff';
sum_zadoff = zeros(1,image_row_size*image_column_size);
for t = 1:image_row_size*image_column_size
    sum_zadoff(t) = sum(PP_zadoff(t,:));
end
sim_zadoff = zeros(image_row_size*image_column_size,image_row_size*image_column_size);
for m = 1:image_row_size*image_column_size
    for n = 1:image_row_size*image_column_size
        sim_zadoff(m,n) = cosine_sim(PP_zadoff(m,:),PP_zadoff(n,:));
    end
end
%ZadoffChu Seq End
%--------------------------------------------------------------------------
%VanderCorupt Seq Start
vd = VanDerCorputSequence(D);
PP_vd = ones(D,image_row_size*image_column_size);
for i = 1:image_row_size*image_column_size    
        for z = 1:D
            if 0.5 <= vd(z,i)
                PP_vd(z,i) = -1;
            end            
        end    
end
PP_vd = PP_vd';
sum_vd = zeros(1,image_row_size*image_column_size);
for t = 1:image_row_size*image_column_size
    sum_vd(t) = sum(PP_vd(t,:));
end
sim_vd = zeros(image_row_size*image_column_size,image_row_size*image_column_size);
for m = 1:image_row_size*image_column_size
    for n = 1:image_row_size*image_column_size
        sim_vd(m,n) = cosine_sim(PP_vd(m,:),PP_vd(n,:));
    end
end
%VanderCorupt Seq End
%--------------------------------------------------------------------------
%Niederreiter Seq Start
[ nd1, ~ ] = niederreiter2_generate ( 20, 1024, 2, 27 );
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

% sim_nd1 = zeros(20,20);
% for m = 1:20
%     for n = 1:20
%         sim_nd1(m,n) = cosine_sim(PP_nd1(m,:),PP_nd1(n,:));
%     end
% end
% sum_nd1 = zeros(1,20);
% for t = 1:20
%     sum_nd1(t) = sum(PP_nd1(t,:));
% end

[ nd2, ~ ] = niederreiter2_generate ( 20, 1024, 2, 28 );
nd2 = nd2';
PP_nd2 = ones(D,20);
for i = 1:20    
        for z = 1:D
            if 4 <= nd2(z,i)
                PP_nd2(z,i) = -1;
            end            
        end    
end
PP_nd2 = PP_nd2';
% sum_nd2 = zeros(1,20);
% for t = 1:20
%     sum_nd2(t) = sum(PP_nd2(t,:));
% end

[ nd3, ~ ] = niederreiter2_generate ( 20, 1024, 2, 29 );
nd3 = nd3';
PP_nd3 = ones(D,20);
for i = 1:20    
        for z = 1:D
            if 2 <= nd3(z,i)
                PP_nd3(z,i) = -1;
            end            
        end    
end
PP_nd3 = PP_nd3';
% sum_nd3 = zeros(1,20);
% for t = 1:20
%     sum_nd3(t) = sum(PP_nd3(t,:));
% end

[ nd4, ~ ] = niederreiter2_generate ( 20, 1024, 2, 30 );
nd4 = nd4';
PP_nd4 = ones(D,20);
for i = 1:20    
        for z = 1:D
            if 1 <= nd4(z,i)
                PP_nd4(z,i) = -1;
            end            
        end    
end
PP_nd4 = PP_nd4';
% sum_nd4 = zeros(1,20);
% for t = 1:20
%     sum_nd4(t) = sum(PP_nd4(t,:));
% end

[ nd5, ~ ] = niederreiter2_generate ( 20, 1024, 2, 31 );
nd5 = nd5';
PP_nd5 = ones(D,20);
for i = 1:20    
        for z = 1:D
            if 0.5 <= nd5(z,i)
                PP_nd5(z,i) = -1;
            end            
        end    
end
PP_nd5 = PP_nd5';
% sum_nd5 = zeros(1,20);
% for t = 1:20
%     sum_nd5(t) = sum(PP_nd5(t,:));
% end

[ nd6, ~ ] = niederreiter2_generate ( 20, 1024, 2, 32 );
nd6 = nd6';
PP_nd6 = ones(D,20);
for i = 1:20    
        for z = 1:D
            if 0.25 <= nd6(z,i)
                PP_nd6(z,i) = -1;
            end            
        end    
end
PP_nd6 = PP_nd6';
% sum_nd6 = zeros(1,20);
% for t = 1:20
%     sum_nd6(t) = sum(PP_nd6(t,:));
% end

[ nd7, ~ ] = niederreiter2_generate ( 20, 1024, 2, 33 );
nd7 = nd7';
PP_nd7 = ones(D,20);
for i = 1:20    
        for z = 1:D
            if 0.125 <= nd7(z,i)
                PP_nd7(z,i) = -1;
            end            
        end    
end
PP_nd7 = PP_nd7';
% sum_nd7 = zeros(1,20);
% for t = 1:20
%     sum_nd7(t) = sum(PP_nd7(t,:));
% end

[ nd8, ~ ] = niederreiter2_generate ( 20, 1024, 2, 34 );
nd8 = nd8';
PP_nd8 = ones(D,20);
for i = 1:20    
        for z = 1:D
            if 0.0625 <= nd8(z,i)
                PP_nd8(z,i) = -1;
            end            
        end    
end
PP_nd8 = PP_nd8';
% sum_nd8 = zeros(1,20);
% for t = 1:20
%     sum_nd8(t) = sum(PP_nd8(t,:));
% end


PP_nd = [PP_nd1', PP_nd2', PP_nd3', PP_nd4', PP_nd5', PP_nd6', PP_nd7', PP_nd8'];
PP_nd = PP_nd';
sum_nd = zeros(1,160);
for t = 1:160
    sum_nd(t) = sum(PP_nd(t,:));
end
sim_nd = zeros(image_row_size*image_column_size,image_row_size*image_column_size);
for m = 1:image_row_size*image_column_size
    for n = 1:image_row_size*image_column_size
        sim_nd(m,n) = cosine_sim(PP_nd(m,:),PP_nd(n,:));
    end
end
%Niederreiter Seq End
%--------------------------------------------------------------------------