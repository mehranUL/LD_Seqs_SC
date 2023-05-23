%Generating bit-stream related to Kasami sequence. 
K_Kasami = kasami(log2(1024));
K_Kasami = K_Kasami(:,1:144);
K_Kasami = K_Kasami';


%Reference: Travis Wiens (2023). Kasami Sequences, m-sequences, Linear Feedback Shift Registers (https://www.mathworks.com/matlabcentral/fileexchange/22716-kasami-sequences-m-sequences-linear-feedback-shift-registers), MATLAB Central File Exchange. Retrieved April 25, 2023.

function [ K ] = kasami(m, feedback)
%Generates the small set of Kasami sequences. Each row in K is a sequence.
%Inputs:
% m - order of sequence (must be even)
% feedback - feedback term for lfsr to generate maximal length sequence.
% See http://www.ece.cmu.edu/~koopman/lfsr/index.html for some values
% (but note these values must be converted from hex to decimal).
%For more details, see Pingzhui Fan and Michael Darnell, "Sequence Design
% for communications applications," Research  Studies Press, Tauton, 1996.

%
%Copyright (c) 2009, Travis Wiens
%All rights reserved.
%
%Redistribution and use in source and binary forms, with or without 
%modification, are permitted provided that the following conditions are 
%met:
%
%    * Redistributions of source code must retain the above copyright 
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright 
%      notice, this list of conditions and the following disclaimer in 
%      the documentation and/or other materials provided with the distribution
%      
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%POSSIBILITY OF SUCH DAMAGE.
%
% If you would like to request that this software be licensed under a less
% restrictive license (i.e. for commercial closed-source use) please
% contact Travis at travis.mlfx@nutaksas.com

use_mex=false;%set to true if you have lfsr_mex working

if nargin<1
    m=8;
end

if mod(m,2)~=0
    error('m must be even')
end

N=2^m-1;%period of sequence

%N_calc=2^(m/2);%number of sequences to calculate

%******Modified By Meh
N_calc=2^m;

%feedbacks=[nan hex2dec('3') hex2dec('5') hex2dec('9') hex2dec('12')...
%    hex2dec('21') hex2dec('41') hex2dec('8e') hex2dec('108') hex2dec('204')...
%    hex2dec('402') hex2dec('829') hex2dec('100D') hex2dec('2015')...
%    hex2dec('4001') hex2dec('8016')];
%maximal length lfsr feedback values from http://www.ece.cmu.edu/~koopman/lfsr/index.html
feedbacks=[nan 3 5 9 18 33 65 142 264 516 1026 2089 4109 8213 16385 32790];
%same as above for those without hex2dec

if nargin<2
    feedback=feedbacks(m);%select appropriate mls feedback
end

f=1+2^(m/2);%sampler
if use_mex
    a=lfsr_mex(feedback,1,N,1)*2-1;%generate base MLS
    b=lfsr_mex(feedback,1,N,f)*2-1;
else
    a=lfsr(feedback,1,N,1)*2-1;%generate base MLS
    b=lfsr(feedback,1,N,f)*2-1;
end

K=zeros(N_calc,N);%allocate memory

K(1,:)=a;%first sequence is the first mls

for i=0:(N_calc-2)
    K(i+2,:)=xor(a==1,circshift(b,[0 i])==1)*2-1;%kasami sequences
end