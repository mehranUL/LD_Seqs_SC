function seq = LFSR_Bulk(N)

    L1 = LFSR__1([true false true true true false false false], N);
    L2 = LFSR__2([false false false true true true false false], N);
    L3 = LFSR__3([true false true true false false false true], N);
    L4 = LFSR__4([false false true true true false true false], N);
    L5 = LFSR__5([true false true true true false false false], N);
    L6 = LFSR__6([true false true false true false false true], N);

    seq = [L1; L2; L3; L4; L5; L6];

end
