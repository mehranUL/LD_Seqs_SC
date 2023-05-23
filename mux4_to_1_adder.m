function [result] = mux4_to_1_adder(x1, x2, x3, x4, s1, s2)

%Bitwise operations
result = (((~s1)&(~s2)& x1) | ((s1)&(~s2)& x2) | ...
 ((~s1)&(s2)& x3) | ((s1)&(s2)& x4));

end