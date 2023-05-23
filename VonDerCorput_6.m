function seq = VonDerCorput_6(N)

    V1 = transpose(vdcorput(N, 2));
    V2 = transpose(vdcorput(N, 4));
    V3 = transpose(vdcorput(N, 8));
    V4 = transpose(vdcorput(N, 16));
    V5 = transpose(vdcorput(N, 64));
    V6 = transpose(vdcorput(N, 128));
    V7 = transpose(vdcorput(N, 256));

    seq = [V1; V2; V3; V4; V5; V6; V7];

end

