clear;
bp = 8; % Bit-precision
%D = 2^(2*bp); % for Multiplication

%D = 2^(bp)*2;   % for Addition D = 512 & for Subtraction use D = 1024
D = 2^(bp);
%D = 2^(bp)/2;   %Subtraction

  %sobol_seq_new = net(sobolset(2), D);
% % 
R2_ws = load("R2_64k.mat");
R2 = R2_ws.z;
%R2 = R2(1:D,1:2);
% 
% alpha = (sqrt(2) - 1); % Silver ratio
% weyl(:,1) = mod((1:D)*alpha, 1);
% beta = pi;
% weyl(:,2) = mod((1:D)*beta, 1);
% 
% Y = lhsdesign(D,20);
% Y = Y(:,13:14);
% 
% ff = faure(D,2,7);
% ff = ff';
% ff = ff(1:D,:);
% 
% p = haltonset(D,'Skip',1e3,'Leap',1e2);
% p = scramble(p,'RR2');
% HT = net(p,2);
% HT = HT';
% HT = net(haltonset(6),D);
% HT = HT(:,5:6);

% % 
% HH = Hammersley(D,3);
% HH = HH';
% HH = HH(:,2:3);

% ws_mohsen = load("mohsen.mat");
% ws_mohsen2 = load("mohsen2.mat");
% temp1 = ws_mohsen2.B;
% HH(:,2) = temp1(1:D);
% HH(:,2) = HH(:,2)/D;
% temp2 = ws_mohsen.A;
% HH(:,3) = temp2(1:D);
% HH(:,3) = HH(:,3)/D;
% HH = HH(:,2:3);


% vd(:,1) = vdcorput(D-1,2);% 2,16
% vd(:,2) = vdcorput(D-1,D);%for addition 8,16 is good--->6.7139e-04
% vd(:,3) = vdcorput(D-1,8);
% vd(:,4) = vdcorput(D-1,16);
% 
% [ nd1, ~ ] = niederreiter2_generate ( 20, D, 2, 31 );
% nd1 = nd1';
% %nd = nd1(:,1:4);
% nd(:,1) = nd1(:,1);
% nd(:,2) = nd1(:,2);

% pt = poissonDisc([100,100,100],7,D); %6 for 2k-- 5 for 4k --- 4 for 8k
% pt(:,1) = pt(:,1)/100;
% pt(:,2) = pt(:,2)/100;

% rnd1(1:D/2) = lfsr_value512(D/2,1);
% rnd1((D/2+1):D) = 1 - rnd1(1:D/2);
%end

% TT_sobol = zeros(1,D);
% % TT_sobol2 = zeros(1,D);
% TT_R2 = zeros(1,D);
% TT_weyl = zeros(1,D);
% TT_latin = zeros(1,D);
% TT_halton = zeros(1,D);
% TT_hammersley = zeros(1,D);
% TT_faure = zeros(1,D);
% TT_nd = zeros(1,D);
% TT_vd = zeros(1,D);
% % TT_vd2 = zeros(1,D);
% TT_pt = zeros(1,D);

%rnd1 = rand(1,N);
%for i = 1:image_row_size*image_column_size
%  for z = 1:D
% %         %if 0.5 <= rnd1(z)
% %         if 0.5 <= sobol_seq_new(z,2)
% %             TT_sobol(z) = 1;
% %         end
%         if 0.5 <= R2(z,2)
%             TT_R2(z) = 1;
%         end
% %         if 0.5 <= weyl(z,2)
% %             TT_weyl(z) = 1;
% %         end
% %         if 0.5 <= Y(z,2)
% %             TT_latin(z) = 1;
% %         end
% %         if 0.5 <= HT(z,2)
% %             TT_halton(z) = 1;
% %         end
% %         if 0.5 <= HH(z,2)
% %             TT_hammersley(z) = 1;
% %         end
% %         if 0.5 <= ff(z,2)
% %             TT_faure(z) = 1;
% %         end
% %         if 0.5 <= nd(z,2)
% %             TT_nd(z) = 1;
% %         end
% %          if 0.5 <= vd(z,2)
% %              TT_vd(z) = 1;
% %          end
% %         if 0.5 <= pt(z,2)
% %             TT_pt(z) = 1;
% %         end
%  end

EE = zeros(2^bp,2^bp);
%AA = zeros(2^bp,2^bp);
%AA_OR = zeros(2^bp,2^bp);
% SS = zeros(2^bp,2^bp);
%SS_mux = zeros(2^bp,2^bp);

%               Multiplication
%  MAE_sobol = zeros(2^bp,2^bp);
MAE_R2 = zeros(2^bp,2^bp);
 %MAE_weyl = zeros(2^bp,2^bp);
 %MAE_vandercorput = zeros(2^bp,2^bp);
 %MAE_latin = zeros(2^bp,2^bp);
%MAE_halton = zeros(2^bp,2^bp);
 %MAE_hammersley = zeros(2^bp,2^bp);
%MAE_ps = zeros(2^bp,2^bp);
 %MAE_nd = zeros(2^bp,2^bp);
 %MAE_faure = zeros(2^bp,2^bp);


%               Addition
%  MAE_sobol_add = zeros(2^bp,2^bp);
  %MAE_vandercorput_add = zeros(2^bp,2^bp);
% MAE_R2_add = zeros(2^bp,2^bp);
% MAE_weyl_add = zeros(2^bp,2^bp);
% MAE_latin_add = zeros(2^bp,2^bp);
% MAE_faure_add = zeros(2^bp,2^bp);
% MAE_nd_add = zeros(2^bp,2^bp);
% MAE_halton_add = zeros(2^bp,2^bp);
% MAE_hammersley_add = zeros(2^bp,2^bp);
% MAE_ps_add = zeros(2^bp,2^bp);

%               Subtraction using XOR gate
% MAE_sobol_sub = zeros(2^bp,2^bp);
% MAE_R2_sub = zeros(2^bp,2^bp);
% MAE_weyl_sub = zeros(2^bp,2^bp);
% MAE_latin_sub = zeros(2^bp,2^bp);
% MAE_faure_sub = zeros(2^bp,2^bp);
% MAE_halton_sub = zeros(2^bp,2^bp);
% MAE_hammersley_sub = zeros(2^bp,2^bp);
% MAE_nd_sub = zeros(2^bp,2^bp);
% MAE_ps_sub = zeros(2^bp,2^bp);
% MAE_vandercorput_sub = zeros(2^bp,2^bp);

%               Subtraction using MUX
%MAE_sobol_sub_mux = zeros(2^bp,2^bp);
% MAE_R2_sub_mux = zeros(2^bp,2^bp);
% MAE_weyl_sub_mux = zeros(2^bp,2^bp);
% MAE_latin_sub_mux = zeros(2^bp,2^bp);
% MAE_faure_sub_mux = zeros(2^bp,2^bp);
% MAE_halton_sub_mux = zeros(2^bp,2^bp);
% MAE_hammersley_sub_mux = zeros(2^bp,2^bp);
% MAE_nd_sub_mux = zeros(2^bp,2^bp);
% MAE_ps_sub_mux = zeros(2^bp,2^bp);
% MAE_vandercorput_sub_mux = zeros(2^bp,2^bp);

% XX1_sobol = zeros(2^bp,D);
% XX2_sobol = zeros(2^bp,D);
XX1_R2 = zeros(2^bp,D);
XX2_R2 = zeros(2^bp,D);
% XX1_weyl = zeros(2^bp,D);
% XX2_weyl = zeros(2^bp,D);
% XX1_vandercorput = zeros(2^bp,D);
% XX2_vandercorput = zeros(2^bp,D);
% XX1_latin = zeros(2^bp,D);
% XX2_latin = zeros(2^bp,D);
% XX1_faure = zeros(2^bp,D);
% XX2_faure = zeros(2^bp,D);
% XX1_halton = zeros(2^bp,D);
% XX2_halton = zeros(2^bp,D);
% XX1_hammersley = zeros(2^bp,D);
% XX2_hammersley = zeros(2^bp,D);
% XX1_nd = zeros(2^bp,D);
% XX2_nd = zeros(2^bp,D);
% XX1_ps = zeros(2^bp,D);
% XX2_ps = zeros(2^bp,D);

%               Multiplication
%  XX_output_sobol = zeros(2^bp,2^bp,D);
XX_output_R2 = zeros(2^bp,2^bp,D);
 %XX_output_weyl = zeros(2^bp,2^bp,D);
%XX_output_vandercorput = zeros(2^bp,2^bp,D);
 %XX_output_latin = zeros(2^bp,2^bp,D);
 %XX_output_halton = zeros(2^bp,2^bp,D);
 %XX_output_hammersley = zeros(2^bp,2^bp,D);
% XX_output_ps = zeros(2^bp,2^bp,D);
 %XX_output_nd = zeros(2^bp,2^bp,D);
 %XX_output_faure = zeros(2^bp,2^bp,D);

%               Addition
% X_add_sobol = zeros(2^bp,2^bp,D);
% X_add_vandercorput = zeros(2^bp,2^bp,D);
% X_add_R2 = zeros(2^bp,2^bp,D);
% X_add_weyl = zeros(2^bp,2^bp,D);
% X_add_latin = zeros(2^bp,2^bp,D);
% X_add_faure = zeros(2^bp,2^bp,D);
% X_add_nd = zeros(2^bp,2^bp,D);
% X_add_halton = zeros(2^bp,2^bp,D);
% X_add_hammersley = zeros(2^bp,2^bp,D);
% X_add_ps = zeros(2^bp,2^bp,D);

%               Subtraction using XOR gate
% X_sub_sobol = zeros(2^bp,2^bp,D);
% X_sub_R2 = zeros(2^bp,2^bp,D);
% X_sub_weyl = zeros(2^bp,2^bp,D);
% X_sub_latin = zeros(2^bp,2^bp,D);
% X_sub_faure = zeros(2^bp,2^bp,D);
% X_sub_halton = zeros(2^bp,2^bp,D);
% X_sub_hammersley = zeros(2^bp,2^bp,D);
% X_sub_nd = zeros(2^bp,2^bp,D);
% X_sub_ps = zeros(2^bp,2^bp,D);
% X_sub_vandercorput = zeros(2^bp,2^bp,D);

%               Subtraction_MUX
% X_sub_mux_sobol = zeros(2^bp,2^bp,D);
% X_sub_mux_R2 = zeros(2^bp,2^bp,D);
% X_sub_mux_weyl = zeros(2^bp,2^bp,D);
% X_sub_mux_latin = zeros(2^bp,2^bp,D);
% X_sub_mux_faure = zeros(2^bp,2^bp,D);
% X_sub_mux_halton = zeros(2^bp,2^bp,D);
% X_sub_mux_hammersley = zeros(2^bp,2^bp,D);
% X_sub_mux_nd = zeros(2^bp,2^bp,D);
% X_sub_mux_ps = zeros(2^bp,2^bp,D);
% X_sub_mux_vandercorput = zeros(2^bp,2^bp,D);



tic

for i = 1:2^bp

    for j = 1:2^bp
        EE(i,j) = ((i-1)/(2^bp))*((j-1)/(2^bp)); % Expected multiplication results
        %AA(i,j) = (((i-1)/(2^bp)) + ((j-1)/(2^bp)))/2;  %Expected Scaled Addition results
        %AA_OR(i,j) = (((i-1)/(2^bp)) + ((j-1)/(2^bp)));
        %SS(i,j) = abs(((i-1)/(2^bp)) - ((j-1)/(2^bp))); % Expected Subtraction results using XOR
        %if (j-1)/(2^bp) < (i-1)/(2^bp) % Checking X2 < X1
         %   SS_mux(i,j) = -100;
        %else    % X2 >= X1
        %SS_mux(i,j) = (j-1)/(2^bp) - (i-1)/(2^bp);
        %SS_mux(i,j) = (((i-1)/(2^bp)) - ((j-1)/(2^bp)))/2;  %Expected Scaled Subtraction
        %end

         for k = 1:D
%              if ((i-1)/(2^bp)) > sobol_seq_new(k,1) %Unipolar
%              %if (1+(i-1)/(2^bp))/2 > sobol_seq_new(k,1)  %Bipolar Encoding ----> BE = (x+1)/2 ----- -1<=x<=1
%                  XX1_sobol(i,k) = 1;
%              end
%              if ((j-1)/(2^bp)) > sobol_seq_new(k,1)
%              %if (1+(j-1)/(2^bp))/2 > sobol_seq_new(k,1)
%                  XX2_sobol(j,k) = 1;
%              end
            if ((i-1)/(2^bp)) > R2(k,1)
            %if (1+(i-1)/(2^bp))/2 > R2(k,1)
                    XX1_R2(i,k) = 1;
            end
            if ((j-1)/(2^bp)) > R2(k,2)
            %if (1+(j-1)/(2^bp))/2 > R2(k,1)
                    XX2_R2(j,k) = 1;
            end
%             if ((i-1)/(2^bp)) > weyl(k,1)
%             %if (1+(i-1)/(2^bp))/2 > weyl(k,1)
%                 XX1_weyl(i,k) = 1;
%             end
%             if ((j-1)/(2^bp)) > weyl(k,1)
%             %if (1+(j-1)/(2^bp))/2 > weyl(k,1)
%                 XX2_weyl(j,k) = 1;
%             end
%              if ((i-1)/(2^bp)) > vd(k,1)
%              %if (1+(i-1)/(2^bp))/2 > vd(k,1)
%                  XX1_vandercorput(i,k) = 1;
%              end
%              if ((j-1)/(2^bp)) > vd(k,1)
%              %if (1+(j-1)/(2^bp))/2 > vd(k,1)
%                  XX2_vandercorput(j,k) = 1;
%              end
%             if ((i-1)/(2^bp)) > Y(k,1)
%             %if (1+(i-1)/(2^bp))/2 > Y(k,1)
%                 XX1_latin(i,k) = 1;
%             end
%             if ((j-1)/(2^bp)) > Y(k,1)
%             %if (1+(j-1)/(2^bp))/2 > Y(k,1)
%                 XX2_latin(j,k) = 1;
%             end
%             if ((i-1)/(2^bp)) > ff(k,1)
%             %if (1+(i-1)/(2^bp))/2 > ff(k,1)
%                 XX1_faure(i,k) = 1;
%             end
%             if ((j-1)/(2^bp)) > ff(k,1)
%             %if (1+(j-1)/(2^bp))/2 > ff(k,1)
%                 XX2_faure(j,k) = 1;
%             end
%             if ((i-1)/(2^bp)) > HT(k,1)
%             %if (1+(i-1)/(2^bp))/2 > HT(k,1)
%                 XX1_halton(i,k) = 1;
%             end
%             if ((j-1)/(2^bp)) > HT(k,1)
%             %if (1+(j-1)/(2^bp))/2 > HT(k,1)
%                 XX2_halton(j,k) = 1; 
%             end
%             if ((i-1)/(2^bp)) > HH(k,1)
%             %if (1+(i-1)/(2^bp))/2 > HH(k,1)
%                 XX1_hammersley(i,k) = 1;
%             end
%             if ((j-1)/(2^bp)) > HH(k,1)
%             %if (1+(j-1)/(2^bp))/2 > HH(k,1)
%                 XX2_hammersley(j,k) = 1;
%             end
             
%             if ((i-1)/(2^bp)) > nd(k,1)
%             %if (1+(i-1)/(2^bp))/2 > nd(k,1)
%                 XX1_nd(i,k) = 1;
%             end
%             if ((j-1)/(2^bp)) > nd(k,2)
%             %if (1+(j-1)/(2^bp))/2 > nd(k,1)
%                 XX2_nd(j,k) = 1;
%             end
%             if ((i-1)/(2^bp)) > pt(k,1)
%             %if (1+(i-1)/(2^bp))/2 > pt(k,1)
%                 XX1_ps(i,k) = 1;
%             end
%             if ((j-1)/(2^bp)) > pt(k,1)
%             %if (1+(j-1)/(2^bp))/2 > pt(k,1)
%                 XX2_ps(j,k) = 1;
%             end
% 
         end
        %               Multiplication
%          XX_output_sobol(i,j,:) = and(XX1_sobol(i,:),XX2_sobol(j,:));    %Bit-stream multiplication (Stochastic AND operation)
        XX_output_R2(i,j,:) = and(XX1_R2(i,:),XX2_R2(j,:));
         %XX_output_weyl(i,j,:) = and(XX1_weyl(i,:),XX2_weyl(j,:));
         %XX_output_vandercorput(i,j,:) = and(XX1_vandercorput(i,:),XX2_vandercorput(j,:));
         %XX_output_latin(i,j,:) = and(XX1_latin(i,:),XX2_latin(j,:));
         %XX_output_halton(i,j,:) = and(XX1_halton(i,:),XX2_halton(j,:));
         %XX_output_hammersley(i,j,:) = and(XX1_hammersley(i,:),XX2_hammersley(j,:));
         %XX_output_ps(i,j,:) = and(XX1_ps(i,:),XX2_ps(j,:));
         %XX_output_nd(i,j,:) = and(XX1_nd(i,:),XX2_nd(j,:));
         %XX_output_faure(i,j,:) = and(XX1_faure(i,:),XX2_faure(j,:));

        %               Addition
        
%           X_add_sobol(i,j,:) = or(and(not(TT_sobol),XX1_sobol(i,:)),and(TT_sobol,XX2_sobol(j,:)));
           %X_add_vandercorput(i,j,:) = or(and(not(TT_vd),XX1_vandercorput(i,:)),and(TT_vd,XX2_vandercorput(j,:)));
         %X_add_R2(i,j,:) = or(and(not(TT_R2),XX1_R2(i,:)),and(TT_R2,XX2_R2(j,:)));
%         X_add_weyl(i,j,:) = or(and(not(TT_weyl),XX1_weyl(i,:)),and(TT_weyl,XX2_weyl(j,:)));
%         X_add_latin(i,j,:) = or(and(not(TT_latin),XX1_latin(i,:)),and(TT_latin,XX2_latin(j,:)));
 %        X_add_faure(i,j,:) = or(and(not(TT_faure),XX1_faure(i,:)),and(TT_faure,XX2_faure(j,:)));
%         X_add_nd(i,j,:) = or(and(not(TT_nd),XX1_nd(i,:)),and(TT_nd,XX2_nd(j,:)));
        % X_add_halton(i,j,:) = or(and(not(TT_halton),XX1_halton(i,:)),and(TT_halton,XX2_halton(j,:)));
%          X_add_hammersley(i,j,:) = or(and(not(TT_hammersley),XX1_hammersley(i,:)),and(TT_hammersley,XX2_hammersley(j,:)));
%         X_add_ps(i,j,:) = or(and(not(TT_pt),XX1_ps(i,:)),and(TT_pt,XX2_ps(j,:)));

%               Addition using OR gate
%         X_add_sobol(i,j,:) = or(XX1_sobol(i,:),XX2_sobol(j,:));
         %X_add_vandercorput(i,j,:) = or(XX1_vandercorput(i,:),XX2_vandercorput(j,:));
%         X_add_R2(i,j,:) = or(XX1_R2(i,:),XX2_R2(j,:));
%         X_add_weyl(i,j,:) = or(XX1_weyl(i,:),XX2_weyl(j,:));
%         X_add_latin(i,j,:) = or(XX1_latin(i,:),XX2_latin(j,:));
%         X_add_faure(i,j,:) = or(XX1_faure(i,:),XX2_faure(j,:));
%         X_add_nd(i,j,:) = or(XX1_nd(i,:),XX2_nd(j,:));
%         X_add_halton(i,j,:) = or(XX1_halton(i,:),XX2_halton(j,:));
%         X_add_hammersley(i,j,:) = or(XX1_hammersley(i,:),XX2_hammersley(j,:));
%         X_add_ps(i,j,:) = or(XX1_ps(i,:),XX2_ps(j,:));


        %               Subtraction using XOR gate
%         X_sub_sobol(i,j,:) = xor(XX1_sobol(i,:),XX2_sobol(j,:));
%         X_sub_R2(i,j,:) = xor(XX1_R2(i,:),XX2_R2(j,:));
%         X_sub_weyl(i,j,:) = xor(XX1_weyl(i,:),XX2_weyl(j,:));
%         X_sub_latin(i,j,:) = xor(XX1_latin(i,:),XX2_latin(j,:));
%         X_sub_faure(i,j,:) = xor(XX1_faure(i,:),XX2_faure(j,:));
%         X_sub_halton(i,j,:) = xor(XX1_halton(i,:),XX2_halton(j,:));
%         X_sub_hammersley(i,j,:) = xor(XX1_hammersley(i,:),XX2_hammersley(j,:));
%         X_sub_nd(i,j,:) = xor(XX1_nd(i,:),XX2_nd(j,:));
%         X_sub_ps(i,j,:) = xor(XX1_ps(i,:),XX2_ps(j,:));
%         X_sub_vandercorput(i,j,:) = xor(XX1_vandercorput(i,:),XX2_vandercorput(j,:));

        
        %               Subtraction using MUX
%         X_sub_mux_sobol(i,j,:) = or(and(not(TT_sobol),not(XX1_sobol(i,:))),and(TT_sobol,XX2_sobol(j,:)));
%         X_sub_mux_vandercorput(i,j,:) = or(and(not(TT_vd),not(XX1_vandercorput(i,:))),and(TT_vd,XX2_vandercorput(j,:)));
%         X_sub_mux_R2(i,j,:) = or(and(not(TT_R2),not(XX1_R2(i,:))),and(TT_R2,XX2_R2(j,:)));
%         X_sub_mux_weyl(i,j,:) = or(and(not(TT_weyl),not(XX1_weyl(i,:))),and(TT_weyl,XX2_weyl(j,:)));
%         X_sub_mux_latin(i,j,:) = or(and(not(TT_latin),not(XX1_latin(i,:))),and(TT_latin,XX2_latin(j,:)));
%         X_sub_mux_faure(i,j,:) = or(and(not(TT_faure),not(XX1_faure(i,:))),and(TT_faure,XX2_faure(j,:)));
%         X_sub_mux_nd(i,j,:) = or(and(not(TT_nd),not(XX1_nd(i,:))),and(TT_nd,XX2_nd(j,:)));
%         X_sub_mux_halton(i,j,:) = or(and(not(TT_halton),not(XX1_halton(i,:))),and(TT_halton,XX2_halton(j,:)));
%         X_sub_mux_hammersley(i,j,:) = or(and(not(TT_hammersley),not(XX1_hammersley(i,:))),and(TT_hammersley,XX2_hammersley(j,:)));
%         X_sub_mux_ps(i,j,:) = or(and(not(TT_pt),not(XX1_ps(i,:))),and(TT_pt,XX2_ps(j,:)));



        %               Multiplication
%          MAE_sobol(i,j) = abs(sum(XX_output_sobol(i,j,:)/D) - EE(i,j));    % Mean Absolute Error (MAE)
        MAE_R2(i,j) = abs(sum(XX_output_R2(i,j,:)/D) - EE(i,j));
         %MAE_weyl(i,j) = abs(sum(XX_output_weyl(i,j,:)/D) - EE(i,j));
         %MAE_vandercorput(i,j) = abs(sum(XX_output_vandercorput(i,j,:)/D) - EE(i,j));
         %MAE_latin(i,j) = abs(sum(XX_output_latin(i,j,:)/D) - EE(i,j));
         %MAE_halton(i,j) = abs(sum(XX_output_halton(i,j,:)/D) - EE(i,j));
         %MAE_hammersley(i,j) = abs(sum(XX_output_hammersley(i,j,:)/D) - EE(i,j));
         %MAE_ps(i,j) = abs(sum(XX_output_ps(i,j,:)/D) - EE(i,j));
         %MAE_nd(i,j) = abs(sum(XX_output_nd(i,j,:)/D) - EE(i,j));
         %MAE_faure(i,j) = abs(sum(XX_output_faure(i,j,:)/D) - EE(i,j));


        %               MUX Addition
%           MAE_sobol_add(i,j) = abs(sum(X_add_sobol(i,j,:)/D) - AA(i,j));
%           MAE_vandercorput_add(i,j) = abs(sum(X_add_vandercorput(i,j,:)/D) - AA(i,j));
         %MAE_R2_add(i,j) = abs(sum(X_add_R2(i,j,:)/D) - AA(i,j));
%         MAE_weyl_add(i,j) = abs(sum(X_add_weyl(i,j,:)/D) - AA(i,j));
%         MAE_latin_add(i,j) = abs(sum(X_add_latin(i,j,:)/D) - AA(i,j));
%         MAE_faure_add(i,j) = abs(sum(X_add_faure(i,j,:)/D) - AA(i,j));
%         MAE_nd_add(i,j) = abs(sum(X_add_nd(i,j,:)/D) - AA(i,j));
        % MAE_halton_add(i,j) = abs(sum(X_add_halton(i,j,:)/D) - AA(i,j));
%          MAE_hammersley_add(i,j) = abs(sum(X_add_hammersley(i,j,:)/D) - AA(i,j));
%         MAE_ps_add(i,j) = abs(sum(X_add_ps(i,j,:)/D) - AA(i,j));


%               Addition using OR gate
%         MAE_sobol_add(i,j) = abs(sum(X_add_sobol(i,j,:)/D) - AA_OR(i,j));
         %MAE_vandercorput_add(i,j) = abs(sum(X_add_vandercorput(i,j,:)/D) - AA_OR(i,j));
%         MAE_R2_add(i,j) = abs(sum(X_add_R2(i,j,:)/D) - AA_OR(i,j));
%         MAE_weyl_add(i,j) = abs(sum(X_add_weyl(i,j,:)/D) - AA_OR(i,j));
%         MAE_latin_add(i,j) = abs(sum(X_add_latin(i,j,:)/D) - AA_OR(i,j));
%         MAE_faure_add(i,j) = abs(sum(X_add_faure(i,j,:)/D) - AA_OR(i,j));
%         MAE_nd_add(i,j) = abs(sum(X_add_nd(i,j,:)/D) - AA_OR(i,j));
%         MAE_halton_add(i,j) = abs(sum(X_add_halton(i,j,:)/D) - AA_OR(i,j));
%         MAE_hammersley_add(i,j) = abs(sum(X_add_hammersley(i,j,:)/D) - AA_OR(i,j));
%         MAE_ps_add(i,j) = abs(sum(X_add_ps(i,j,:)/D) - AA_OR(i,j));


        %               Subtraction using XOR gate
%         MAE_sobol_sub(i,j) = abs(sum(X_sub_sobol(i,j,:)/D) - SS(i,j));
%         MAE_R2_sub(i,j) = abs(sum(X_sub_R2(i,j,:)/D) - SS(i,j));
%         MAE_weyl_sub(i,j) = abs(sum(X_sub_weyl(i,j,:)/D) - SS(i,j));
%         MAE_latin_sub(i,j) = abs(sum(X_sub_latin(i,j,:)/D) - SS(i,j));
%         MAE_faure_sub(i,j) = abs(sum(X_sub_faure(i,j,:)/D) - SS(i,j));
%         MAE_halton_sub(i,j) = abs(sum(X_sub_halton(i,j,:)/D) - SS(i,j));
%         MAE_hammersley_sub(i,j) = abs(sum(X_sub_hammersley(i,j,:)/D) - SS(i,j));
%         MAE_nd_sub(i,j) = abs(sum(X_sub_nd(i,j,:)/D) - SS(i,j));
%         MAE_ps_sub(i,j) = abs(sum(X_sub_ps(i,j,:)/D) - SS(i,j));
%         MAE_vandercorput_sub(i,j) = abs(sum(X_sub_vandercorput(i,j,:)/D) - SS(i,j));


        %               Subtraction with MUX
%         if SS_mux(i,j) == -100
%             MAE_ps_sub_mux(i,j) = 0;
%             MAE_sobol_sub_mux(i,j) = 0;
%             MAE_vandercorput_sub_mux(i,j) = 0;
%             MAE_R2_sub_mux(i,j) = 0;
%             MAE_weyl_sub_mux(i,j) = 0;
%             MAE_latin_sub_mux(i,j) = 0;
%             MAE_faure_sub_mux(i,j) = 0;
%             MAE_halton_sub_mux(i,j) = 0;
%             MAE_hammersley_sub_mux(i,j) = 0;
%             MAE_nd_sub_mux(i,j) = 0;
%         else

          if i<j
             %MAE_sobol_sub_mux(i,j) = abs((2*sum(X_sub_mux_sobol(i,j,:))/D - 1) - (-SS_mux(i,j)));
%             MAE_R2_sub_mux(i,j) = abs((2*sum(X_sub_mux_R2(i,j,:))/D - 1) - (-SS_mux(i,j)));
%             MAE_weyl_sub_mux(i,j) = abs((2*sum(X_sub_mux_weyl(i,j,:))/D - 1) - (-SS_mux(i,j)));
%             MAE_latin_sub_mux(i,j) = abs((2*sum(X_sub_mux_latin(i,j,:))/D - 1) - (-SS_mux(i,j)));
%             MAE_faure_sub_mux(i,j) = abs((2*sum(X_sub_mux_faure(i,j,:))/D - 1) - (-SS_mux(i,j)));
%             MAE_halton_sub_mux(i,j) = abs((2*sum(X_sub_mux_halton(i,j,:))/D - 1) - (-SS_mux(i,j)));
%             MAE_hammersley_sub_mux(i,j) = abs((2*sum(X_sub_mux_hammersley(i,j,:))/D - 1) - (-SS_mux(i,j)));
%             MAE_nd_sub_mux(i,j) = abs((2*sum(X_sub_mux_nd(i,j,:))/D - 1) - (-SS_mux(i,j)));
%             MAE_ps_sub_mux(i,j) = abs((2*sum(X_sub_mux_ps(i,j,:))/D - 1) - (-SS_mux(i,j)));
             %MAE_vandercorput_sub_mux(i,j) = abs((2*sum(X_sub_mux_vandercorput(i,j,:))/D - 1) - (-SS_mux(i,j)));
          else
             %MAE_sobol_sub_mux(i,j) = abs(abs(2*sum(X_sub_mux_sobol(i,j,:))/D - 1) - SS_mux(i,j));        %Bipolar: (sum1's -sum0's)/N = 2*P(1) -1
%             MAE_R2_sub_mux(i,j) = abs(abs(2*sum(X_sub_mux_R2(i,j,:))/D - 1) - SS_mux(i,j));
%             MAE_weyl_sub_mux(i,j) = abs(abs(2*sum(X_sub_mux_weyl(i,j,:))/D - 1) - SS_mux(i,j));
%             MAE_latin_sub_mux(i,j) = abs(abs(2*sum(X_sub_mux_latin(i,j,:))/D - 1) - SS_mux(i,j));
%             MAE_faure_sub_mux(i,j) = abs(abs(2*sum(X_sub_mux_faure(i,j,:))/D - 1) - SS_mux(i,j));
%             MAE_halton_sub_mux(i,j) = abs(abs(2*sum(X_sub_mux_halton(i,j,:))/D - 1) - SS_mux(i,j));
%             MAE_hammersley_sub_mux(i,j) = abs(abs(2*sum(X_sub_mux_hammersley(i,j,:))/D - 1) - SS_mux(i,j));
%             MAE_nd_sub_mux(i,j) = abs(abs(2*sum(X_sub_mux_nd(i,j,:))/D - 1) - SS_mux(i,j));
%             MAE_ps_sub_mux(i,j) = abs(abs(2*sum(X_sub_mux_ps(i,j,:))/D - 1) - SS_mux(i,j));
             %MAE_vandercorput_sub_mux(i,j) = abs(abs(2*sum(X_sub_mux_vandercorput(i,j,:))/D - 1) - SS_mux(i,j));
          end

%             MAE_sobol_sub_mux(i,j) = abs((sum(X_sub_mux_sobol(i,j,:))/D ) - SS_mux(i,j));        %Unipolar: sum1's/N = P(1)
%             MAE_R2_sub_mux(i,j) = abs((sum(X_sub_mux_R2(i,j,:))/D ) - SS_mux(i,j));
%             MAE_weyl_sub_mux(i,j) = abs((sum(X_sub_mux_weyl(i,j,:))/D ) - SS_mux(i,j));
%             MAE_latin_sub_mux(i,j) = abs((sum(X_sub_mux_latin(i,j,:))/D ) - SS_mux(i,j));
%             MAE_faure_sub_mux(i,j) = abs((sum(X_sub_mux_faure(i,j,:))/D ) - SS_mux(i,j));
%             MAE_halton_sub_mux(i,j) = abs((sum(X_sub_mux_halton(i,j,:))/D ) - SS_mux(i,j));
%             MAE_hammersley_sub_mux(i,j) = abs((sum(X_sub_mux_hammersley(i,j,:))/D ) - SS_mux(i,j));
%             MAE_nd_sub_mux(i,j) = abs((sum(X_sub_mux_nd(i,j,:))/D ) - SS_mux(i,j));
%             MAE_ps_sub_mux(i,j) = abs((sum(X_sub_mux_ps(i,j,:))/D ) - SS_mux(i,j));
%             MAE_vandercorput_sub_mux(i,j) = abs((sum(X_sub_mux_vandercorput(i,j,:))/D ) - SS_mux(i,j));

%         end

       

    end
end
% final_MAE_sobol = sum(MAE_sobol_sub_mux,"all")/(D*D)
% final_MAE_R2 = sum(MAE_R2_sub_mux,"all")/(D*D)
% final_MAE_weyl = sum(MAE_weyl_sub_mux,"all")/(D*D)
% final_MAE_latin = sum(MAE_latin_sub_mux,"all")/(D*D)
% final_MAE_faure = sum(MAE_faure_sub_mux,"all")/(D*D)
% final_MAE_halton = sum(MAE_halton_sub_mux,"all")/(D*D)
% final_MAE_hammersley = sum(MAE_hammersley_sub_mux,"all")/(D*D)
% final_MAE_nd = sum(MAE_nd_sub_mux,"all")/(D*D)
% final_MAE_vd = sum(MAE_vandercorput_sub_mux,"all")/(D*D)
% final_MAE_ps = sum(MAE_ps_sub_mux,"all")/(D*D)

%  MAE_S = sum(MAE_sobol_add,"all")/((2^bp)*(2^bp))
 MAE_R = sum(MAE_R2,"all")/((2^bp)*(2^bp))
%  MAE_W = sum(MAE_weyl_add,"all")/((2^bp)*(2^bp))
% MAE_V = sum(MAE_vandercorput_add,"all")/((2^bp)*(2^bp))
% MAE_L = sum(MAE_latin_add,"all")/((2^bp)*(2^bp))
% MAE_F = sum(MAE_faure_add,"all")/((2^bp)*(2^bp))
% MAE_HL = sum(MAE_halton_add,"all")/((2^bp)*(2^bp))
%  MAE_HM = sum(MAE_hammersley_add,"all")/((2^bp)*(2^bp))
 %MAE_N = sum(MAE_nd,"all")/((2^bp)*(2^bp))
%  MAE_P = sum(MAE_ps_add,"all")/((2^bp)*(2^bp))

toc