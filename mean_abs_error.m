
clear all
close all

N = (2^8); %stream size
i = 1; %for iteration and vector index
number_of_test_iteration = 1000; % total random tests

gold = (sqrt(5) - 1) / 2;
alpha = (sqrt(2) - 1); % Silver ratio
beta = pi;
alpha2 = exp(1);
weyl(:,1) = mod((1:N)*alpha, 1);
weyl(:,2) = mod((1:N)*beta, 1);
weyl(:,3) = mod((1:N)*gold, 1);
weyl(:,4) = mod((1:N)*alpha2, 1);

% ws16k1 = load('sobol_bl_optimized_16k_91.mat','aa8k1');
% aa16k = ws16k1.aa8k1;
% sobol_seq_new = net(sobolset(1111), N);
%sobol_seq_new = sobol_seq1(:,aa16k);
sobol_seq_new = net(sobolset(4), N);

% R2 Sequence
R2_ws = load("R2_65536.mat");
%R2_ws = load("R2_2048.mat");
%R2_ws = load("R2_4096.mat");
%R2_ws = load("R2_8192.mat");
R2 = R2_ws.z;
R2 = R2(1:N,1:4);

%Kasami Sequence
%K_Kasami = kasami(log2(N));

%Latin Hypercube
L = lhsdesign(N,20);
L = L(:,13:16);

%Faure sequence
%ff = faure(N,2,2);
F = faure(N,4,4);
F = F';
%ff = ff(2:end,:);

%Halton Sequence
p = haltonset(N,'Skip',1e3,'Leap',1e2);
p = scramble(p,'RR2');
HL = net(p,4);
HL = HL';

%Hammersley sequence
%HH = Hammersley(N,3);
HM = Hammersley(N,5);
HM = HM';
HM = HM(:,2:5);

%ZadoffChu Sequence
% a = zeros(1,257);
% %f = factor(1023);
% for i =2:257
%     if mod(i,257) ~= 0 
%         if mod(257,i) ~= 0
%             a(i) = i;
%         end
%     end
% end
% b = find(a ~= 0);
% 
% % f = factor(1023);
% % b = 1:1023;
% % a = setdiff(b,f);
% % a = a(2:numel(a));
% 
% re = zeros(257,numel(b));
% for zd = 1:numel(b)
%     re(:,zd) = real(zadoffChuSeq(b(zd),257));
% end
% re = abs(re(2:end,1:2));

a = zeros(1,1023);
%f = factor(1023);
for i =2:1023
    if mod(i,3) ~= 0 && mod(i,11) ~= 0 && mod(i,31) ~= 0
        if mod(1023,i) ~= 0
            a(i) = i;
        end
    end
end
b = find(a ~= 0);

% f = factor(1023);
% b = 1:1023;
% a = setdiff(b,f);
% a = a(2:numel(a));

re = zeros(1023,numel(b));
for zd = 1:numel(b)
    re(:,zd) = real(zadoffChuSeq(b(zd),1023));
end
re(1024,:) = -0.01;
%re = abs(re(2:end,1:2));

% ZadoffChu sequence of length 2K.
% a2 = zeros(1,2047);
% %f = factor(1023);
% for i =2:2047
%     if mod(i,23) ~= 0 && mod(i,89) ~= 0
%         if mod(2047,i) ~= 0
%             a2(i) = i;
%         end
%     end
% end
% b2 = find(a2 ~= 0);
% 
% re = zeros(2047,numel(b2));
% for zd = 1:numel(b2)
%     re(:,zd) = real(zadoffChuSeq(b2(zd),2047));
% end
% re(2048,:) = -0.01;

% a2 = zeros(1,4095);
% %f = factor(1023);
% for i =2:4095
%     if mod(i,3) ~= 0 && mod(i,5) ~= 0 && mod(i,7) ~= 0 && mod(i,13) ~= 0
%         if mod(4095,i) ~= 0
%             a2(i) = i;
%         end
%     end
% end
% b2 = find(a2 ~= 0);
% 
% re = zeros(4095,numel(b2));
% for zd = 1:numel(b2)
%     re(:,zd) = real(zadoffChuSeq(b2(zd),4095));
% end
% re(4096,:) = -0.01;

% re = zeros(8191,2);
% for zd = 1:2
%     re(:,zd) = real(zadoffChuSeq(zd,8191));
% end
% re(8192,:) = -0.01;


%VanDerCorput sequence
vd(:,1) = vdcorput(N-1,2);
vd(:,2) = vdcorput(N-1,8);%for addition 4,8 is good
vd(:,3) = vdcorput(N-1,256);
vd(:,4) = vdcorput(N-1,4);

%Niederreiter sequence
[ nd1, ~ ] = niederreiter2_generate ( 20, N, 2, 31 );
nd1 = nd1';
%nd = nd1(:,1:4);
nd(:,1) = nd1(:,1);
nd(:,2) = nd1(:,2);
nd(:,3) = nd1(:,3);
nd(:,4) = nd1(:,4);

%Poisson Disk sequence
%pt = poissonDisc([100,100,100],8,N); %6 for 2k-- 5 for 4k --- 4 for 8k
pt = poissonDisc([100,100,100,100],8,N);
pt(:,1) = pt(:,1)/100;
pt(:,2) = pt(:,2)/100;
pt(:,3) = pt(:,3)/100;
pt(:,4) = pt(:,4)/100;

%Half-LFSR method
%rnd1 = zeros(N,1);
%for i = 1:image_row_size*image_column_size  
%     rnd1(1:N/2) = lfsr_value512(N/2,1);
%     rnd1((N/2+1):N) = 1 - rnd1(1:N/2);
%end

%TT = zeros(1,N);
%rnd1 = rand(1,N);
%for i = 1:image_row_size*image_column_size
%     for z = 1:N
%         if 0.5 <= rnd1(z)
%             TT(z) = 1;
%         end
%     end
%end

% TT_sobol = zeros(1,N);
% TT_sobol2 = zeros(1,N);
% TT_vd = zeros(1,N);
% TT_vd2 = zeros(1,N);
% TT_R2 = zeros(1,N);
% TT_R22 = zeros(1,N);
% TT_weyl = zeros(1,N);
% TT_weyl2 = zeros(1,N);
% TT_latin = zeros(1,N);
% TT_latin2 = zeros(1,N);
% TT_halton = zeros(1,N);
% TT_halton2 = zeros(1,N);
% TT_hammersley = zeros(1,N);
% TT_hammersley2 = zeros(1,N);
% TT_faure = zeros(1,N);
% TT_faure2 = zeros(1,N);
% TT_nd = zeros(1,N);
% TT_nd2 = zeros(1,N);
% TT_pt = zeros(1,N);
% TT_pt2 = zeros(1,N);

%  for z = 1:N
%      if 0.5 <= sobol_seq_new(z,3)
%              TT_sobol(z) = 1;
%      end
%     if 0.5 <= sobol_seq_new(z,4)
%             TT_sobol2(z) = 1;
%     end
%      if 0.5 <= vd(z,3)
%               TT_vd(z) = 1;
%      end
%     if 0.5 <= vd(z,4)
%              TT_vd2(z) = 1;
%     end
%     if 0.5 <= R2(z,3)
%               TT_R2(z) = 1;
%      end
%     if 0.5 <= R2(z,4)
%              TT_R22(z) = 1;
%     end
%     if 0.5 <= weyl(z,3)
%               TT_weyl(z) = 1;
%      end
%     if 0.5 <= weyl(z,4)
%              TT_weyl2(z) = 1;
%     end
%     if 0.5 <= L(z,3)
%               TT_latin(z) = 1;
%      end
%     if 0.5 <= L(z,4)
%              TT_latin2(z) = 1;
%     end
%     if 0.5 <= F(z,3)
%               TT_faure(z) = 1;
%      end
%     if 0.5 <= F(z,4)
%              TT_faure2(z) = 1;
%     end
%     if 0.5 <= HL(z,3)
%               TT_halton(z) = 1;
%      end
%     if 0.5 <= HL(z,4)
%              TT_halton2(z) = 1;
%     end
%     if 0.5 <= HM(z,3)
%               TT_hammersley(z) = 1;
%      end
%     if 0.5 <= HM(z,4)
%              TT_hammersley2(z) = 1;
%     end
%     if 0.5 <= nd(z,3)
%               TT_nd(z) = 1;
%      end
%     if 0.5 <= nd(z,4)
%              TT_nd2(z) = 1;
%     end
%     if 0.5 <= pt(z,3)
%               TT_pt(z) = 1;
%      end
%     if 0.5 <= pt(z,4)
%              TT_pt2(z) = 1;
%     end
%  end


%load("sobolset.mat"); %for random only change here
%-load("randomset.mat");
%SOURCE
%seq = zeros(256,2);
%seq(:,1) = round((255-0).*rand(1,256) + 0)/255;
%seq(:,2) = round((255-0).*rand(1,256) + 0)/255;


for k = 1:number_of_test_iteration %total test iteration, above, you may change the final value

    %----------------Random X1 X2 X3 scalar selection----------------------
    % https://www.mathworks.com/help/matlab/ref/randi.html
    X1 = randi((N+1),1)-1; %pseudorandom integer from a uniform discrete distribution
    % range -> [0,N]

    X2 = randi((N+1),1)-1; %pseudorandom integer from a uniform discrete distribution
    % range -> [0,N]
    X3 = randi((N+1),1)-1;
    X4 = randi((N+1),1)-1;
    %-----------------------Expected Value---------------------------------
    %expected_value(i) = (X1/N)*(X2/N);
    expected_value(i) = (X1/N)*(X2/N)*(X3/N);
    %addition_value(i) = ((X1/N) + (X2/N))/2;
    %addition_value(i) = ((X1/N) + (X2/N) + (X3/N) + (X4/N))/4;
    %-----------------------Expected Value---------------------------------

    X1_stream_sobol = zeros(1, N);
    X2_stream_sobol = zeros(1, N);
    X3_stream_sobol = zeros(1, N);
    X4_stream_sobol = zeros(1, N);
    X1_stream_weyl = zeros(1, N);
    X2_stream_weyl = zeros(1, N);
    X3_stream_weyl = zeros(1, N);
    X4_stream_weyl = zeros(1, N);
    X1_stream_R2 = zeros(1, N);
    X2_stream_R2 = zeros(1, N);
    X3_stream_R2 = zeros(1, N);
    X4_stream_R2 = zeros(1, N);
%     X1_stream_kasami = zeros(1, N);
%     X2_stream_kasami = zeros(1, N);
    X1_stream_latin = zeros(1, N);
    X2_stream_latin = zeros(1, N);
    X3_stream_latin = zeros(1, N);
    X4_stream_latin = zeros(1, N);
    X1_stream_faure = zeros(1, N);
    X2_stream_faure = zeros(1, N);
    X3_stream_faure = zeros(1, N);
    X4_stream_faure = zeros(1, N);
    X1_stream_halton = zeros(1, N);
    X2_stream_halton = zeros(1, N);
    X3_stream_halton = zeros(1, N);
    X4_stream_halton = zeros(1, N);
    X1_stream_hammersley = zeros(1, N);
    X2_stream_hammersley = zeros(1, N);
    X3_stream_hammersley = zeros(1, N);
    X4_stream_hammersley = zeros(1, N);
%     X1_stream_zadoff = zeros(1, N);
%     X2_stream_zadoff = zeros(1, N);
    X1_stream_vandercorput = zeros(1, N);
    X2_stream_vandercorput = zeros(1, N);
    X3_stream_vandercorput = zeros(1, N);
    X4_stream_vandercorput = zeros(1, N);
    X1_stream_nd = zeros(1, N);
    X2_stream_nd = zeros(1, N);
    X3_stream_nd = zeros(1, N);
    X4_stream_nd = zeros(1, N);
    X1_stream_ps = zeros(1, N);
    X2_stream_ps = zeros(1, N);
    X3_stream_ps = zeros(1, N);
    X4_stream_ps = zeros(1, N);


    for z=1:N
        if ((X1/N) > sobol_seq_new(z,1))
            X1_stream_sobol(1,z) = 1;
        end
        if ((X2/N) > sobol_seq_new(z,2))
            X2_stream_sobol(1,z) = 1;
        end
        if ((X3/N) > sobol_seq_new(z,3))
            X3_stream_sobol(1,z) = 1;
        end
        if ((X4/N) > sobol_seq_new(z,4))
            X4_stream_sobol(1,z) = 1;
        end
        if ((X1/N) > weyl(z,1))
            X1_stream_weyl(1,z) = 1;
        end
        if ((X2/N) > weyl(z,2))
            X2_stream_weyl(1,z) = 1;
        end
        if ((X3/N) > weyl(z,3))
            X3_stream_weyl(1,z) = 1;
        end
        if ((X4/N) > weyl(z,4))
            X4_stream_weyl(1,z) = 1;
        end
        if ((X1/N) > R2(z,1))
            X1_stream_R2(1,z) = 1;
        end
        if ((X2/N) > R2(z,2))
            X2_stream_R2(1,z) = 1;
        end
        if ((X3/N) > R2(z,3))
            X3_stream_R2(1,z) = 1;
        end
        if ((X4/N) > R2(z,4))
            X4_stream_R2(1,z) = 1;
        end
%         if ((X1/N) > K_Kasami(z,1))
%             X1_stream_kasami(1,z) = 1;
%         end
        if ((X1/N) > L(z,1))
            X1_stream_latin(1,z) = 1;
        end
        if ((X2/N) > L(z,2))
            X2_stream_latin(1,z) = 1;
        end
        if ((X3/N) > L(z,3))
            X3_stream_latin(1,z) = 1;
        end
        if ((X4/N) > L(z,4))
            X4_stream_latin(1,z) = 1;
        end
        if ((X1/N) > F(z,1))
            X1_stream_faure(1,z) = 1;
        end
        if ((X2/N) > F(z,2))
            X2_stream_faure(1,z) = 1;
        end
        if ((X3/N) > F(z,3))
            X3_stream_faure(1,z) = 1;
        end
        if ((X4/N) > F(z,4))
            X4_stream_faure(1,z) = 1;
        end
        if ((X1/N) > HL(z,1))
            X1_stream_halton(1,z) = 1;
        end
        if ((X2/N) > HL(z,2))
            X2_stream_halton(1,z) = 1;
        end
        if ((X3/N) > HL(z,3))
            X3_stream_halton(1,z) = 1;
        end
        if ((X4/N) > HL(z,4))
            X4_stream_halton(1,z) = 1;
        end
        if ((X1/N) > HM(z,1))
            X1_stream_hammersley(1,z) = 1;
        end
        if ((X2/N) > HM(z,2))
            X2_stream_hammersley(1,z) = 1;
        end
        if ((X3/N) > HM(z,3))
            X3_stream_hammersley(1,z) = 1;
        end
        if ((X4/N) > HM(z,4))
            X4_stream_hammersley(1,z) = 1;
        end
%         if ((X1/N) > re(z,1))
%             X1_stream_zadoff(1,z) = 1;
%         end
        if ((X1/N) > vd(z,1))
            X1_stream_vandercorput(1,z) = 1;
        end
        if ((X2/N) > vd(z,2))
            X2_stream_vandercorput(1,z) = 1;
        end
        if ((X3/N) > vd(z,3))
            X3_stream_vandercorput(1,z) = 1;
        end
        if ((X4/N) > vd(z,4))
            X4_stream_vandercorput(1,z) = 1;
        end
        if ((X1/N) > nd(z,1))
            X1_stream_nd(1,z) = 1;
        end
        if ((X2/N) > nd(z,2))
            X2_stream_nd(1,z) = 1;
        end
        if ((X3/N) > nd(z,3))
            X3_stream_nd(1,z) = 1;
        end
        if ((X4/N) > nd(z,4))
            X4_stream_nd(1,z) = 1;
        end
        if ((X1/N) > pt(z,1))
            X1_stream_ps(1,z) = 1;
        end
        if ((X2/N) > pt(z,2))
            X2_stream_ps(1,z) = 1;
        end
        if ((X3/N) > pt(z,3))
            X3_stream_ps(1,z) = 1;
        end
        if ((X4/N) > pt(z,4))
            X4_stream_ps(1,z) = 1;
        end
        

%         if ((X1/N) <= seq(z,1))
%             X1_stream(1,z) = 0;
%         end
    end

%     for z=1:N
%         if ((X2/N) > sobol_seq_new(z,2))
%             X2_stream_sobol(1,z) = 1;
%         end
%         if ((X2/N) > weyl(z,2))
%             X2_stream_weyl(1,z) = 1;
%         end
%         if ((X2/N) > R2(z,2))
%             X2_stream_R2(1,z) = 1;
%         end
%         if ((X2/N) > K_Kasami(z,2))
%             X2_stream_kasami(1,z) = 1;
%         end
%         if ((X2/N) > Y(z,2))
%             X2_stream_latin(1,z) = 1;
%         end
%         if ((X2/N) > ff(z,2))
%             X2_stream_faure(1,z) = 1;
%         end
%         if ((X2/N) > X(z,2))
%             X2_stream_halton(1,z) = 1;
%         end
%         if ((X2/N) > HH(z,2))
%             X2_stream_hammersley(1,z) = 1;
%         end
%         if ((X2/N) > re(z,2))
%             X2_stream_zadoff(1,z) = 1;
%         end
%         if ((X2/N) > vd(z,2))
%             X2_stream_vandercorput(1,z) = 1;
%         end
%         if ((X2/N) > nd(z,2))
%             X2_stream_nd(1,z) = 1;
%         end
%         if ((X2/N) > pt(z,2))
%             X2_stream_ps(1,z) = 1;
%         end
% 
% %         if ((X2/N) <= seq(z,2))
% %             X2_stream(1,z) = 0;
% %         end
%     end

%               2-input Multiplication

     %X_output_sobol = and(X1_stream_sobol, X2_stream_sobol);
%     X_output_weyl = and(X1_stream_weyl, X2_stream_weyl);
%     X_output_R2 = and(X1_stream_R2, X2_stream_R2);
%     X_output_kasami = and(X1_stream_kasami, X2_stream_kasami);
%     X_output_latin = and(X1_stream_latin, X2_stream_latin);
%     X_output_faure = and(X1_stream_faure, X2_stream_faure);
%     X_output_halton = and(X1_stream_halton, X2_stream_halton);
%     X_output_hammersley = and(X1_stream_hammersley, X2_stream_hammersley);
%     X_output_zadoff = and(X1_stream_zadoff, X2_stream_zadoff);
     %X_output_vandercorput = and(X1_stream_vandercorput, X2_stream_vandercorput);
%     X_output_nd = and(X1_stream_nd, X2_stream_nd);
%     X_output_ps = and(X1_stream_ps, X2_stream_ps);


%           3-input Multiplication
    X_output_sobol = and(X1_stream_sobol, and(X2_stream_sobol,X3_stream_sobol));
    X_output_weyl = and(X1_stream_weyl, and(X2_stream_weyl,X3_stream_weyl));
    X_output_R2 = and(X1_stream_R2, and(X2_stream_R2,X3_stream_R2));
    %X_output_kasami = and(X1_stream_kasami, and(X2_stream_kasami,X3_st));
    X_output_latin = and(X1_stream_latin, and(X2_stream_latin,X3_stream_latin));
    X_output_faure = and(X1_stream_faure, and(X2_stream_faure,X3_stream_faure));
    X_output_halton = and(X1_stream_halton, and(X2_stream_halton,X3_stream_halton));
    X_output_hammersley = and(X1_stream_hammersley, and(X2_stream_hammersley,X3_stream_hammersley));
    %X_output_zadoff = and(X1_stream_zadoff, X2_stream_zadoff);
    X_output_vandercorput = and(X1_stream_vandercorput, and(X2_stream_vandercorput,X3_stream_vandercorput));
    X_output_nd = and(X1_stream_nd, and(X2_stream_nd,X3_stream_nd));
    X_output_ps = and(X1_stream_ps, and(X2_stream_ps,X3_stream_ps));

    %and_tmp1 = and(not(TT),X1_stream_sobol);
    %and_tmp2 = and(TT,X2_stream_sobol);

%               2-input Addition
     %X_add_sobol = or(and(not(TT_sobol),X1_stream_sobol),and(TT_sobol,X2_stream_sobol));
%     X_add_weyl = or(and(not(TT),X1_stream_weyl),and(TT,X2_stream_weyl));
%     X_add_R2 = or(and(not(TT),X1_stream_R2),and(TT,X2_stream_R2));
%     X_add_kasami = or(and(not(TT),X1_stream_kasami),and(TT,X2_stream_kasami));
%     X_add_latin = or(and(not(TT),X1_stream_latin),and(TT,X2_stream_latin));
%     X_add_faure = or(and(not(TT),X1_stream_faure),and(TT,X2_stream_faure));
%     X_add_halton = or(and(not(TT),X1_stream_halton),and(TT,X2_stream_halton));
%     X_add_hammersley = or(and(not(TT),X1_stream_hammersley),and(TT,X2_stream_hammersley));
%     X_add_zadoff = or(and(not(TT),X1_stream_zadoff),and(TT,X2_stream_zadoff));
     %X_add_vandercorput = or(and(not(TT_vd),X1_stream_vandercorput),and(TT_vd,X2_stream_vandercorput));
%     X_add_nd = or(and(not(TT),X1_stream_nd),and(TT,X2_stream_nd));
%     X_add_ps = or(and(not(TT),X1_stream_ps),and(TT,X2_stream_ps));


%                       4-input Addition
% temp1s = and(and(not(TT_sobol),not(TT_sobol2)),X1_stream_sobol);
% temp2s = and(and(TT_sobol,not(TT_sobol2)),X2_stream_sobol);
% temp3s = and(and(not(TT_sobol),TT_sobol2),X3_stream_sobol);
% temp4s = and(and(TT_sobol,TT_sobol2),X4_stream_sobol);
% temp5s = or(temp1s,temp2s);
% temp6s = or(temp3s,temp4s);
% X_add_sobol = or(temp5s,temp6s) ;
% 
% temp1v = and(and(not(TT_vd),not(TT_vd2)),X1_stream_vandercorput);
% temp2v = and(and(TT_vd,not(TT_vd2)),X2_stream_vandercorput);
% temp3v = and(and(not(TT_vd),TT_vd2),X3_stream_vandercorput);
% temp4v = and(and(TT_vd,TT_vd2),X4_stream_vandercorput);
% temp5v = or(temp1v,temp2v);
% temp6v = or(temp3v,temp4v);
% X_add_vandercorput = or(temp5v,temp6v) ;
% 
% temp1r = and(and(not(TT_R2),not(TT_R22)),X1_stream_R2);
% temp2r = and(and(TT_R2,not(TT_R22)),X2_stream_R2);
% temp3r = and(and(not(TT_R2),TT_R22),X3_stream_R2);
% temp4r = and(and(TT_R2,TT_R22),X4_stream_R2);
% temp5r = or(temp1r,temp2r);
% temp6r = or(temp3r,temp4r);
% X_add_R2 = or(temp5r,temp6r);
% 
% temp1w = and(and(not(TT_weyl),not(TT_weyl2)),X1_stream_weyl);
% temp2w = and(and(TT_weyl,not(TT_weyl2)),X2_stream_weyl);
% temp3w = and(and(not(TT_weyl),TT_weyl2),X3_stream_weyl);
% temp4w = and(and(TT_weyl,TT_weyl2),X4_stream_weyl);
% temp5w = or(temp1w,temp2w);
% temp6w = or(temp3w,temp4w);
% X_add_weyl = or(temp5w,temp6w) ;
% 
% temp1l = and(and(not(TT_latin),not(TT_latin2)),X1_stream_latin);
% temp2l = and(and(TT_latin,not(TT_latin2)),X2_stream_latin);
% temp3l = and(and(not(TT_latin),TT_latin2),X3_stream_latin);
% temp4l = and(and(TT_latin,TT_latin2),X4_stream_latin);
% temp5l = or(temp1l,temp2l);
% temp6l = or(temp3l,temp4l);
% X_add_latin = or(temp5l,temp6l) ;
% 
% temp1f = and(and(not(TT_faure),not(TT_faure2)),X1_stream_faure);
% temp2f = and(and(TT_faure,not(TT_faure2)),X2_stream_faure);
% temp3f = and(and(not(TT_faure),TT_faure2),X3_stream_faure);
% temp4f = and(and(TT_faure,TT_faure2),X4_stream_faure);
% temp5f = or(temp1f,temp2f);
% temp6f = or(temp3f,temp4f);
% X_add_faure = or(temp5f,temp6f) ;
% 
% temp1hl = and(and(not(TT_halton),not(TT_halton2)),X1_stream_halton);
% temp2hl = and(and(TT_halton,not(TT_halton2)),X2_stream_halton);
% temp3hl = and(and(not(TT_halton),TT_halton2),X3_stream_halton);
% temp4hl = and(and(TT_halton,TT_halton2),X4_stream_halton);
% temp5hl = or(temp1hl,temp2hl);
% temp6hl = or(temp3hl,temp4hl);
% X_add_halton = or(temp5hl,temp6hl) ;
% 
% temp1hm = and(and(not(TT_hammersley),not(TT_hammersley2)),X1_stream_hammersley);
% temp2hm = and(and(TT_hammersley,not(TT_hammersley2)),X2_stream_hammersley);
% temp3hm = and(and(not(TT_hammersley),TT_hammersley2),X3_stream_hammersley);
% temp4hm = and(and(TT_hammersley,TT_hammersley2),X4_stream_hammersley);
% temp5hm = or(temp1hm,temp2hm);
% temp6hm = or(temp3hm,temp4hm);
% X_add_hammersley = or(temp5hm,temp6hm) ;
% 
% temp1n = and(and(not(TT_nd),not(TT_nd2)),X1_stream_nd);
% temp2n = and(and(TT_nd,not(TT_nd2)),X2_stream_nd);
% temp3n = and(and(not(TT_nd),TT_nd2),X3_stream_nd);
% temp4n = and(and(TT_nd,TT_nd2),X4_stream_nd);
% temp5n = or(temp1n,temp2n);
% temp6n = or(temp3n,temp4n);
% X_add_nd = or(temp5n,temp6n) ;
% 
% temp1p = and(and(not(TT_pt),not(TT_pt2)),X1_stream_ps);
% temp2p = and(and(TT_pt,not(TT_pt2)),X2_stream_ps);
% temp3p = and(and(not(TT_pt),TT_pt2),X3_stream_ps);
% temp4p = and(and(TT_pt,TT_pt2),X4_stream_ps);
% temp5p = or(temp1p,temp2p);
% temp6p = or(temp3p,temp4p);
% X_add_ps = or(temp5p,temp6p) ;

    %AbsErr_CM_vs_UM_no_flip(i) = abs((sum(X_output)/N) - (expected_value(i)));
    AbsErr_sobol(i) = abs((sum(X_output_sobol)/N) - (expected_value(i)));
    AbsErr_weyl(i) = abs((sum(X_output_weyl)/N) - (expected_value(i)));
    AbsErr_R2(i) = abs((sum(X_output_R2)/N) - (expected_value(i)));
    %AbsErr_kasami(i) = abs((sum(X_output_kasami)/N) - (expected_value(i)));
    AbsErr_latin(i) = abs((sum(X_output_latin)/N) - (expected_value(i)));
    AbsErr_faure(i) = abs((sum(X_output_faure)/N) - (expected_value(i)));
    AbsErr_halton(i) = abs((sum(X_output_halton)/N) - (expected_value(i)));
    AbsErr_hammersley(i) = abs((sum(X_output_hammersley)/N) - (expected_value(i)));
    %AbsErr_zadoff(i) = abs((sum(X_output_zadoff)/N) - (expected_value(i)));
    AbsErr_vandercorput(i) = abs((sum(X_output_vandercorput)/N) - (expected_value(i)));
    AbsErr_nd(i) = abs((sum(X_output_nd)/N) - (expected_value(i)));
    AbsErr_ps(i) = abs((sum(X_output_ps)/N) - (expected_value(i)));

%     AbsErr_sobol_add(i) = abs((sum(X_add_sobol)/N) - (addition_value(i)));
%     AbsErr_weyl_add(i) = abs((sum(X_add_weyl)/N) - (addition_value(i)));
%     AbsErr_R2_add(i) = abs((sum(X_add_R2)/N) - (addition_value(i)));
%     %AbsErr_kasami_add(i) = abs((sum(X_add_kasami)/N) - (addition_value(i)));
%     AbsErr_latin_add(i) = abs((sum(X_add_latin)/N) - (addition_value(i)));
%     AbsErr_faure_add(i) = abs((sum(X_add_faure)/N) - (addition_value(i)));
%     AbsErr_halton_add(i) = abs((sum(X_add_halton)/N) - (addition_value(i)));
%     AbsErr_hammersley_add(i) = abs((sum(X_add_hammersley)/N) - (addition_value(i)));
%     %AbsErr_zadoff_add(i) = abs((sum(X_add_zadoff)/N) - (addition_value(i)));
%      AbsErr_vandercorput_add(i) = abs((sum(X_add_vandercorput)/N) - (addition_value(i)));
%     AbsErr_nd_add(i) = abs((sum(X_add_nd)/N) - (addition_value(i)));
%     AbsErr_ps_add(i) = abs((sum(X_add_ps)/N) - (addition_value(i)));

      i = i + 1; %for counting iteration

end

%sum(AbsErr_CM_vs_UM_no_flip)/number_of_test_iteration
MAE_sobol = sum(AbsErr_sobol)/number_of_test_iteration
MAE_weyl = sum(AbsErr_weyl)/number_of_test_iteration
MAE_R2 = sum(AbsErr_R2)/number_of_test_iteration
%MAE_kasami = sum(AbsErr_kasami)/number_of_test_iteration
MAE_latin = sum(AbsErr_latin)/number_of_test_iteration
MAE_faure = sum(AbsErr_faure)/number_of_test_iteration
MAE_halton = sum(AbsErr_halton)/number_of_test_iteration
MAE_hammersley = sum(AbsErr_hammersley)/number_of_test_iteration
%MAE_zadoff = sum(AbsErr_zadoff)/number_of_test_iteration
MAE_vandercorput = sum(AbsErr_vandercorput)/number_of_test_iteration
MAE_nd = sum(AbsErr_nd)/number_of_test_iteration
MAE_ps = sum(AbsErr_ps)/number_of_test_iteration

% MAE_sobol_add = sum(AbsErr_sobol_add)/number_of_test_iteration
% MAE_weyl_add = sum(AbsErr_weyl_add)/number_of_test_iteration
% MAE_R2_add = sum(AbsErr_R2_add)/number_of_test_iteration
% %MAE_kasami_add = sum(AbsErr_kasami_add)/number_of_test_iteration
% MAE_latin_add = sum(AbsErr_latin_add)/number_of_test_iteration
% MAE_faure_add = sum(AbsErr_faure_add)/number_of_test_iteration
% MAE_halton_add = sum(AbsErr_halton_add)/number_of_test_iteration
% MAE_hammersley_add = sum(AbsErr_hammersley_add)/number_of_test_iteration
% % MAE_zadoff_add = sum(AbsErr_zadoff_add)/number_of_test_iteration
% MAE_vandercorput_add = sum(AbsErr_vandercorput_add)/number_of_test_iteration
% MAE_nd_add = sum(AbsErr_nd_add)/number_of_test_iteration
% MAE_ps_add = sum(AbsErr_ps_add)/number_of_test_iteration


%XX_output_sobol = and(XX1, XX2);

