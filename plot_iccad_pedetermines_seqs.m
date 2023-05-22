D = 1024;

%load("LD_seqs");
load("zadoff");
%load("vd2");
sobolseq2 = net(sobolset(2),(D));

R2_ws = load("R2_test.mat");
R2 = R2_ws.z;
%R2 = R2(1:D,1:2);

ff = faure(D,2,7);
ff = ff';

alpha = (sqrt(2) - 1); % Silver ratio
weyl(:,1) = mod((1:D)*alpha, 1);
beta = pi;
weyl(:,2) = mod((1:D)*beta, 1);
% 
Y = lhsdesign(D,20);
Y = Y(:,13:14);

HT = net(haltonset(6),D);
HT = HT(:,5:6);

HH = Hammersley(D,3);
HH = HH';
HH = HH(:,2:3);

vd(:,1) = vdcorput(D-1,2);% 2,16
vd(:,2) = vdcorput(D-1,D);

[ nd1, ~ ] = niederreiter2_generate ( 20, D, 2, 31 );
nd1 = nd1';

pt = poissonDisc([100,100,100],7,D); %6 for 2k-- 5 for 4k --- 4 for 8k
pt(:,1) = pt(:,1)/100;
pt(:,2) = pt(:,2)/100;

% seq = zadoffChuSeq(2,1023);
% zadoffchu = real(seq);


figure
subplot(3,4,1)
plot(weyl(1:500, 1), weyl(1:500, 2), 'o', 'linewidth', 1)
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', xt) 
set(findall(gcf,'-property','FontSize'),'FontSize', 28, 'FontName', 'consolas')
title("Weyl (W)")


subplot(3,4,2)
plot(R2(1:500, 1), R2(1:500, 2), 'o', 'linewidth', 1)
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', xt) 
set(findall(gcf,'-property','FontSize'),'FontSize', 28, 'FontName', 'consolas')
title("R2 (R)")


subplot(3,4,3)
plot(Y(1:500, 1), Y(1:500, 2), 'o', 'linewidth', 1)
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', xt) 
set(findall(gcf,'-property','FontSize'),'FontSize', 28, 'FontName', 'consolas')
title("Latin Hypercube (L)")

subplot(3,4,4)
plot(ff(1:500, 1), ff(1:500, 2), 'o', 'linewidth', 1)
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', xt) 
set(findall(gcf,'-property','FontSize'),'FontSize', 28, 'FontName', 'consolas')
title("Faure (F)")

subplot(3,4,5)
plot(sobolseq2(1:500, 1), sobolseq2(1:500, 2), 'o', 'linewidth', 1)
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', xt) 
set(findall(gcf,'-property','FontSize'),'FontSize', 28, 'FontName', 'consolas')
title("Sobol (S)")



subplot(3,4,6)
plot(HT(1:500, 1), HT(1:500, 2), 'o', 'linewidth', 1)
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', xt) 
set(findall(gcf,'-property','FontSize'),'FontSize', 28, 'FontName', 'consolas')
title("Halton (HL)")


subplot(3,4,7)
plot(HH(1:500, 1), HH(1:500, 2), 'o', 'linewidth', 1)
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', xt) 
set(findall(gcf,'-property','FontSize'),'FontSize', 28, 'FontName', 'consolas')
title("Hemmersly (HH)")

subplot(3,4,8)
plot(zadoffchu(1:500, 1), zadoffchu(1:500, 2), 'o', 'linewidth', 1)
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', xt) 
set(findall(gcf,'-property','FontSize'),'FontSize', 28, 'FontName', 'consolas')
title("Zadoff–Chu (Z)")

%sequences = transpose(sequences);

subplot(3,4,9)
plot(vd(1:500, 1), vd(1:500, 2), 'o', 'linewidth', 1)
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', xt) 
set(findall(gcf,'-property','FontSize'),'FontSize', 28, 'FontName', 'consolas')
title("Van der Corput (VDC)")



subplot(3,4,10)
plot(nd1(1:500, 1), nd1(1:500, 2), 'o', 'linewidth', 1)
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', xt) 
set(findall(gcf,'-property','FontSize'),'FontSize', 28, 'FontName', 'consolas')
title("Niederreitter (N)")


subplot(3,4,11)
plot(pt(1:500, 1), pt(1:500, 2), 'o', 'linewidth', 1)
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', xt) 
set(findall(gcf,'-property','FontSize'),'FontSize', 28, 'FontName', 'consolas')
title("Poisson Disk (P)")


subplot(3,4,12)
r1 = rand(500,1);
r2 = rand(500,1);
plot(r1(:, 1), r2(:, 1), 'ro', 'linewidth', 1)
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', xt) 
set(findall(gcf,'-property','FontSize'),'FontSize', 28, 'FontName', 'consolas')
title("Random")

