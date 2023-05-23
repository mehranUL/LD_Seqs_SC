clear all
close all

N = 256;
test_run = 100000;

%sequences = transpose(net(sobolset(6),(N)));
sequences = VonDerCorput_6(N);
%sequences = transpose(rand(N,6));
%sequences = (Hammersley(N,6))
%sequences = niederreiter2_generate(20, N, 2, 31);

%sequences = LFSR_Bulk(N)/N;

for iter = 1:1:test_run
    i = randi(N+1,1)-1;
    j = randi(N+1,1)-1;
    m = randi(N+1,1)-1;
    n = randi(N+1,1)-1;

    X1_Sobol_BS = number_source(i/N, N, sequences(1, :));
    X2_Sobol_BS = number_source(j/N, N, sequences(2, :));
    X3_Sobol_BS = number_source(m/N, N, sequences(4, :));
    X4_Sobol_BS = number_source(n/N, N, sequences(, :));

    Sel = number_source(0.5, N, sequences(6, :));

    ANDControl1 = double((X1_Sobol_BS)&(X2_Sobol_BS));
    ANDControl2 = double(((X3_Sobol_BS)&(X4_Sobol_BS)));
    AND_AND(iter) = correlation(ANDControl1, ANDControl2, N);

    MUXControl1 = double((X1_Sobol_BS)&(~Sel));
    MUXControl2 = double((X2_Sobol_BS)&(Sel));
    MUX(iter) = correlation(MUXControl1, MUXControl2, N);

    ORControl1 = double((X1_Sobol_BS)|(X2_Sobol_BS));
    ORControl2 = double((X3_Sobol_BS)|(X4_Sobol_BS));
    OR_OR(iter) = correlation(ORControl1, ORControl2, N);

    XORControl1 = xor((X1_Sobol_BS),(X2_Sobol_BS));
    XORControl2 = xor((X3_Sobol_BS),(X4_Sobol_BS));
    XOR_XOR(iter) = correlation(XORControl1, XORControl2, N);

    ANDControl1 = double((X1_Sobol_BS)&(X2_Sobol_BS));
    ORControl2 = double((X3_Sobol_BS)|(X4_Sobol_BS));
    AND_OR(iter) = correlation(ANDControl1, ORControl2, N);

    ANDControl1 = double((X1_Sobol_BS)&(X2_Sobol_BS));
    XORControl2 = xor((X3_Sobol_BS), (X4_Sobol_BS));
    AND_XOR(iter) = correlation(ANDControl1, XORControl2, N);

    ORControl1 = double((X1_Sobol_BS)|(X2_Sobol_BS));
    XORControl2 = xor((X3_Sobol_BS), (X4_Sobol_BS));
    OR_XOR(iter) = correlation(ORControl1, XORControl2, N);

end

% sum(abs(AND))/test_run
% sum(abs(OR))/test_run
% sum(abs(XOR))/test_run
% 
% sign_AND = sign(AND);
% sign_OR = sign(OR);
% sign_XOR = sign(XOR);
% 
% SCC_POSITIVE_AND = sum(sign_AND(:)==1)
% SCC_NEGATIVE_AND=sum(sign_AND(:)==-1)
% SCC_ZERO_AND=sum(sign_AND(:)==0)
% 
% SCC_POSITIVE_OR=sum(sign_OR(:)==1)
% SCC_NEGATIVE_OR=sum(sign_OR(:)==-1)
% SCC_ZERO_OR=sum(sign_OR(:)==0)
% 
% SCC_POSITIVE_XOR=sum(sign_XOR(:)==1)
% SCC_NEGATIVE_XOR=sum(sign_XOR(:)==-1)
% SCC_ZERO_XOR=sum(sign_XOR(:)==0)




subplot(2,3,1)
hist(AND_AND);
xlabel("SCC");
ylabel("Cartesian Couples");
title("AND-AND");
ax = gca;
ax.LineWidth = 2;
set(gca,'FontWeight','bold')


subplot(2,3,2)
hist(OR_OR)
xlabel("SCC");
ylabel("Cartesian Couples");
title("OR-OR");
ax = gca;
ax.LineWidth = 2;
set(gca,'FontWeight','bold')

subplot(2,3,3)
hist(XOR_XOR)
xlabel("SCC");
ylabel("Cartesian Couples");
title("XOR-XOR");
ax = gca;
ax.LineWidth = 2;
set(gca,'FontWeight','bold')

subplot(2,3,4)
hist(AND_OR)
xlabel("SCC");
ylabel("Cartesian Couples");
title("AND-OR");
ax = gca;
ax.LineWidth = 2;
set(gca,'FontWeight','bold')


subplot(2,3,5)
hist(AND_XOR)
xlabel("SCC");
ylabel("Cartesian Couples");
title("AND-XOR");
ax = gca;
ax.LineWidth = 2;
set(gca,'FontWeight','bold')

subplot(2,3,6)
hist(OR_XOR)
xlabel("SCC");
ylabel("Cartesian Couples");
title("OR-XOR");
ax = gca;
ax.LineWidth = 2;
set(gca,'FontWeight','bold')

hold on

sgtitle('LFSR-based Generator') 

set(findall(gcf,'-property','FontSize'),'FontSize', 30, 'FontName', 'consolas')


% subplot(1,7,7)
% hist(MUX)
% xlabel("SCC");
% ylabel("Cartesian Couples");
% title("MUX");




