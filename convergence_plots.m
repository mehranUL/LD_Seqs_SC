figure
S = [0.92 0.45 0.19 0.092 0.041 0.019 0.009 0.0035 0.0013 0.0003 0.0000];
R = [1.14 1.07 0.48 0.220 0.130 0.055 0.037 0.0164 0.0099 0.0078 0.0024];
W = [1.46 1.40 0.83 0.530 0.400 0.220 0.190 0.1800 0.1300 0.0090 0.0065];
L = [3.19 0.93 0.85 0.380 0.390 0.250 0.240 0.1200 0.0795 0.0508 0.0424];
F = [2.60 1.40 0.88 0.480 0.210 0.110 0.077 0.0360 0.0136 0.0113 0.0040];
HL = [3.31 1.42 1.14 0.780 0.380 0.150 0.093 0.0570 0.0330 0.0150 0.0083];
HM = [1.31 0.85 0.37 0.200 0.120 0.061 0.030 0.0170 0.0098 0.0043 0.0019];
N = [0.95 0.51 0.34 0.130 0.072 0.032 0.019 0.0067 0.0039 0.0015 0.0011];
V = [1.76 0.88 0.39 0.170 0.073 0.030 0.012 0.0045 0.0015 0.0003 0.0000];
P = [4.16 2.61 1.06 0.960 0.710 0.340 0.440 0.3800 0.3700 0.6000 0.5200];

X = [2^6 2^7 2^8 2^9 2^10 2^11 2^12 2^13 2^14 2^15 2^16];
%X = [6 7 8 9 10 11 12 13 14 15 16];

semilogx(X,S, '-o', 'MarkerSize', 30, 'LineWidth', 8, 'DisplayName', "S")
hold on
xlabel("N = 2^n (log scale)")
ylabel("MAE(%)")
title("MAE(%) of 2-input AND Multiplier")


grid on
box on
ax = gca;
ax.LineWidth = 4;
legend



semilogx(X,R, '-->', 'MarkerSize', 30, 'LineWidth', 8, 'DisplayName', "R")
hold on

semilogx(X,W, '-', 'MarkerSize', 30, 'LineWidth', 8, 'DisplayName', "W")
hold on

semilogx(X,L, '--^', 'MarkerSize', 24, 'LineWidth', 8, 'DisplayName', "L")
hold on

semilogx(X,F, '-*', 'MarkerSize', 15, 'LineWidth', 8, 'DisplayName', "F")
hold on

semilogx(X,V, '-x', 'MarkerSize', 30, 'LineWidth', 8, 'DisplayName', "V")
hold on

semilogx(X,HL, ':*', 'MarkerSize', 30, 'LineWidth', 8, 'DisplayName', "HL")
hold on

semilogx(X,HM, ':', 'MarkerSize', 20, 'LineWidth', 8, 'DisplayName', "HM")
hold on

semilogx(X,N, '-*', 'MarkerSize', 30, 'LineWidth', 8, 'DisplayName', "N")
hold on

semilogx(X,P, '--o', 'MarkerSize', 30, 'LineWidth', 8, 'DisplayName', "P")
hold on

set(findall(gcf,'-property','FontSize'),'FontSize', 36, 'FontName', 'consolas')
set(gca,'FontWeight','bold')










figure

S = [14.98 7.43 3.66 1.77 0.83 0.37 0.150 0.00];
R = [12.69 13.37 8.68 5.04 1.91 1.66 0.600 0.25];
W = [12.94 6.74 4.52 3.07 2.19 1.17 0.920 0.79];
L = [8.54 5.30 3.92 2.95 1.81 0.97 1.020 0.86];
F = [31.35 13.39 11.48 6.52 2.96 2.59 1.620 0.78];
HL = [37.94 21.66 14.86 9.74 3.95 2.09 1.580 0.75];
HM = [14.98 7.49 5.55 2.49 1.19 0.85 0.290 0.18];
N = [10.43 5.74 3.18 1.65 0.81 0.36 0.150 0.00];
V = [13.40 6.63 3.24 1.55 0.71 0.29 0.097 0.00];
P = [26.46 16.31 5.66 6.57 8.16 2.69 1.530 1.65];

X = [2^2 2^3 2^4 2^5 2^6 2^7 2^8 2^9];

semilogx(X,S, '-o', 'MarkerSize', 30, 'LineWidth', 8, 'DisplayName', "S")
hold on
xlabel("N = 2^n (log scale)")
ylabel("MAE(%)")
title("MAE(%) of 2-input MUX Adder")


grid on
box on
ax = gca;
ax.LineWidth = 4;
legend

semilogx(X,R, '-->', 'MarkerSize', 30, 'LineWidth', 8, 'DisplayName', "R")
hold on

semilogx(X,W, '-', 'MarkerSize', 30, 'LineWidth', 8, 'DisplayName', "W")
hold on

semilogx(X,L, '--^', 'MarkerSize', 24, 'LineWidth', 8, 'DisplayName', "L")
hold on

semilogx(X,F, '-*', 'MarkerSize', 15, 'LineWidth', 8, 'DisplayName', "F")
hold on

semilogx(X,V, '-x', 'MarkerSize', 30, 'LineWidth', 8, 'DisplayName', "V")
hold on

semilogx(X,HL, ':*', 'MarkerSize', 30, 'LineWidth', 8, 'DisplayName', "HL")
hold on

semilogx(X,HM, ':', 'MarkerSize', 20, 'LineWidth', 8, 'DisplayName', "HM")
hold on

semilogx(X,N, '-*', 'MarkerSize', 30, 'LineWidth', 8, 'DisplayName', "N")
hold on

semilogx(X,P, '--o', 'MarkerSize', 30, 'LineWidth', 8, 'DisplayName', "P")
hold on

set(findall(gcf,'-property','FontSize'),'FontSize', 36, 'FontName', 'consolas')
set(gca,'FontWeight','bold')
