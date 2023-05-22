D = 1024;
image_row_size = 12;
image_column_size = 12;

%Weyl Sequence Start
% alpha = (sqrt(5) - 1) / 2;
% %alpha = (sqrt(2) - 1); % Silver ratio
% weyl = mod((1:D)*alpha, 1);
% ww = zeros(D,10);
% for i = 1:10
%     ww(:,i) = weyl(randperm(D));
% end
alpha = (sqrt(2) - 1); % Silver ratio
weyl(:,1) = mod((1:D)*alpha, 1);
beta = pi;
weyl(:,2) = mod((1:D)*beta, 1);

PP_weyl = ones(D,2);
shift_weyl = ones(D,2);

for i = 1:2   
        for z = 1:D
            if 0.5 <= weyl(z,i)
                PP_weyl(z,i) = -1;
            end            
        end 
        shift_weyl(:,i) = circshift(PP_weyl(:,i),i-1);
end
%PP_weyl = PP_weyl';
%Weyl Sequence End

%Kasami Sequence Start
K_Kasami = kasami(log2(D));
K_Kasami = K_Kasami(:,1:2);
shift_kasami = ones(D,2);
for i =1:2
    shift_kasami(:,i) = circshift(K_Kasami(:,i),i-1);
end
%K_Kasami = K_Kasami';
%Kasami Sequence End

%Sobol Start
sobol_seq_new = net(sobolset(2), D);
PP_sobol = ones(D,2);
shift_sobol = ones(D,2);

for i = 1:2   
        for z = 1:D
            if 0.5 <= sobol_seq_new(z,i)
                PP_sobol(z,i) = -1;
            end            
        end
        shift_sobol(:,i) = circshift(PP_sobol(:,i),i-1);
end
%Sobol End
PP_Y = ones(D,2);
shift_latin = ones(D,2);
Y = lhsdesign(D,2);
for i = 1:2   
        for z = 1:D
            if 0.5 <= Y(z,i)
                PP_Y(z,i) = -1;
            end            
        end 
        shift_latin(:,i) = circshift(PP_Y(:,i),i-1);
end
%PP_Y = PP_Y';
%Latin Hypercube Sequence End
%Half_LFSR Sequence Start
%N = log2(D);
%gg = logical(randi(2, [1 N]) - 1);
%

% rnd1 = zeros(D,2);
% shift_half = ones(D,2);
% 
% for i = 1:2  
%     rnd1(1:D/2,i) = lfsr_value512(D/2,i);
%     rnd1((D/2+1):D,i) = 1 - rnd1(1:D/2,i);
% end
% 
% TT = ones(D,10);
% 
% for i = 1:10
%     for z = 1:D
%         if 0.5 <= rnd1(z,i)
%             TT(z,i) = -1;
%         end
%     end
%     shift_half(:,i) = circshift(TT(:,i),i-1);
% end
%TT = TT';
%Half_LFSR Sequence End

%Hadamard Start
h = hadamard(4096);
h1 = h(2:1025,2:3);
shift_hadamard = ones(D,2);
for i =1:2
    shift_hadamard(:,i) = circshift(h1(:,i),i-1);
end
%Hadamard End
G = gold_code();
G = G(:,1:2);
shift_gold = ones(D,2);
for i =1:2
    shift_gold(:,i) = circshift(G(:,i),i-1);
end
%G = G';
%Gold Code Sequence End
%Faure Sequence Start
ff = faure(D,2,2); %last parameter is base. base>2 are better
%ff = ff(:,2:end);
ff = ff';
PP_ff = ones(D,2);
shift_faure = ones(D,2);
for i = 1:2    
        for z = 1:D
            if 0.5 <= ff(z,i)
                PP_ff(z,i) = -1;
            end            
        end 
        shift_faure(:,i) = circshift(PP_ff(:,i),i-1);
end
%PP_ff = PP_ff';
%Faure Sequence End
p = haltonset(D,'Skip',1e3,'Leap',1e2);
p = scramble(p,'RR2');
X = net(p,2);
X = X';
PP_halton = ones(D,2);
shift_halton = ones(D,2);
for i = 1:2    
        for z = 1:D
            if 0.499 <= X(z,i)
                PP_halton(z,i) = -1;
            end            
        end 
        shift_halton(:,i) = circshift(PP_halton(:,i),i-1);
end
%PP_halton = PP_halton';
%Halton Sequence End
%Niederreiter Start
[ nd5, ~ ] = niederreiter2_generate ( 20, 1024, 2, 31 );
nd5 = nd5';
nd5 = nd5(:,1:2);
PP_nd5 = ones(D,2);
shift_nd = ones(D,2);
for i = 1:2    
        for z = 1:D
            if 0.5 <= nd5(z,i)
                PP_nd5(z,i) = -1;
            end            
        end
        shift_nd(:,i) = circshift(PP_nd5(:,i),i-1);
end
%Niederreiter End
%VanDerCorput Start
% vd(:,1) = vdcorput(1023,2);
% vd(:,2) = vdcorput(1023,3);
% vd(:,3) = vdcorput(1023,4);
% vd(:,4) = vdcorput(1023,8);
% vd(:,5) = vdcorput(1023,16);
% vd(:,6) = vdcorput(1023,5);
% vd(:,7) = vdcorput(1023,6);
% vd(:,8) = vdcorput(1023,7);
% vd(:,9) = vdcorput(1023,9);
% vd(:,10) = vdcorput(1023,10);

vd(:,1) = vdcorput(D-1,4);% 2,16
vd(:,2) = vdcorput(D-1,8);%for addition 4,8 is good

PP_vd = ones(D,2);
shift_vd = ones(D,2);
for i = 1:2    
        for z = 1:D
            if 0.5 <= vd(z,i)
                PP_vd(z,i) = -1;
            end            
        end 
        shift_vd(:,i) = circshift(PP_vd(:,i),i-1);
end
%VanDerCorput End
%R2 Starts
R2_ws = load("R2_1024_50.mat");
R2 = R2_ws.z;
%R2 = R2(1:D,1:2);
PP_R2 = ones(D,2);
shift_R2 = ones(D,2);
for i = 1:2    
        for z = 1:D
            if 0.5 <= R2(z,i)
                PP_R2(z,i) = -1;
            end            
        end
        shift_R2(:,i) = circshift(PP_R2(:,i),i-1);
end
%R2 Ends
%Hammersley Start
HH = Hammersley(D,2);
HH = HH';
PP_hammersley = ones(D,2);
shift_hammersley = ones(D,2);
for i = 1:2 
    for z = 1:D
         if 0.43 <= HH(z,i)
              PP_hammersley(z,i) = -1;
         end            
    end
    shift_hammersley(:,i) = circshift(PP_hammersley(:,i),i-1);
end 
%Hammersley End


sim_weyl = zeros(2, 2);
sim_weyl_shift = zeros(2, 2);
%sim_kasami = zeros(10, 10);
sim_sobol = zeros(2, 2);
sim_sobol_shift = zeros(2, 2);
%sim_latin = zeros(10, 10);
%sim_half = zeros(10, 10);
sim_hadamard = zeros(2, 2);
sim_hadamard_shift = zeros(2, 2);
sim_gold = zeros(2, 2);
sim_gold_shift = zeros(2, 2);
sim_faure = zeros(2, 2);
sim_faure_shift = zeros(2, 2);
sim_halton = zeros(2, 2);
sim_halton_shift = zeros(2, 2);
sim_nd = zeros(2, 2);
sim_nd_shift = zeros(2, 2);
sim_vd = zeros(2, 2);
sim_vd_shift = zeros(2, 2);
sim_R2 = zeros(2, 2);
sim_R2_shift = zeros(2, 2);
sim_hammersley = zeros(2, 2);
sim_hammersley_shift = zeros(2, 2);

for i = 1:2
    for j = 1:2
        sim_weyl(i,j) = cosine_sim(PP_weyl(:,i),PP_weyl(:,j));
        sim_weyl_shift(i,j) = cosine_sim(shift_weyl(:,i),PP_weyl(:,j));
        %sim_kasami(i,j) = cosine_sim(K_Kasami(:,i),K_Kasami(:,j));
        sim_sobol(i,j) = cosine_sim(PP_sobol(:,i),PP_sobol(:,j));
        sim_sobol_shift(i,j) = cosine_sim(shift_sobol(:,i),PP_sobol(:,j));
        %sim_latin(i,j) = cosine_sim(PP_Y(:,i),PP_Y(:,j));
        %sim_half(i,j) = cosine_sim(TT(:,i),TT(:,j));
        sim_hadamard(i,j) = cosine_sim(h1(:,i),h1(:,j));
        sim_hadamard_shift(i,j) = cosine_sim(shift_hadamard(:,i),h1(:,j));
        sim_gold(i,j) = cosine_sim(G(:,i),G(:,j));
        sim_gold_shift(i,j) = cosine_sim(shift_gold(:,i),G(:,j));
        sim_faure(i,j) = cosine_sim(PP_ff(:,i),PP_ff(:,j));
        sim_faure_shift(i,j) = cosine_sim(shift_faure(:,i),PP_ff(:,j));
        sim_halton(i,j) = cosine_sim(PP_halton(:,i),PP_halton(:,j));
        sim_halton_shift(i,j) = cosine_sim(shift_halton(:,i),PP_halton(:,j));
        sim_nd(i,j) = cosine_sim(PP_nd5(:,i),PP_nd5(:,j));
        sim_nd_shift(i,j) = cosine_sim(shift_nd(:,i),PP_nd5(:,j));
        sim_vd(i,j) = cosine_sim(PP_vd(:,i),PP_vd(:,j));
        sim_vd_shift(i,j) = cosine_sim(shift_vd(:,i),PP_vd(:,j));
        sim_R2(i,j) = cosine_sim(PP_R2(:,i),PP_R2(:,j));
        sim_R2_shift(i,j) = cosine_sim(shift_R2(:,i),PP_R2(:,j));
        sim_hammersley(i,j) = cosine_sim(PP_hammersley(:,i),PP_hammersley(:,j));
        sim_hammersley_shift(i,j) = cosine_sim(shift_hammersley(:,i),PP_hammersley(:,j));
    end
end
