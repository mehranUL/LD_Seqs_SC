% The gold code is ready to use because it consists of -1s and +1s 
G = gold_code();
G = G(:,1:144);
G = G';

% Reference: Sanjeet Kumar (2023). Gold Code Generation (https://www.mathworks.com/matlabcentral/fileexchange/24474-gold-code-generation), MATLAB Central File Exchange. Retrieved April 25, 2023.

function Gold_Seq = gold_code()

G=3067;  % Code length 931 ---> 3*1022 + 1
x=[];
%...............Generation of first perferred PN sequence............
sd1 =[0 1 0 1 1];      % Initial State of Register.
PS1=[];                       
for j=1:G        
    PS1=[PS1 sd1(5)];
    if sd1(1)==sd1(4)
        temp1=0;
    else temp1=1;
    end
    sd1(1)=sd1(2);
    sd1(2)=sd1(3);
    sd1(3)=sd1(4);
    sd1(4)=sd1(5);
    sd1(5)=temp1;
end
x=[x PS1];
%.................Generation of Second Preferred sequnces..............
PS2=[];
PS2(1)=x(1);
for i=1:1021 %309
    j=(3*i)+1;
    PS2(i+1)=x(j);
end
PS2=[PS2];
 
%.................Shifting and Storing of PS1 in Matrix 'y'............
for k=1:1022
    for j=1:1022
        y(k,j)=x(j+k-1);
    end
end
%..................Generation of Gold Sequences........................
for i=1:1022 %310
Gold_Seq(1,:)=[PS1(1,(1:1022))]; %310
Gold_Seq(2,:)=[PS2];
Gold_Seq(i+2,:)=xor(PS2,y(i,(1:1022))); %310
end
% for j=1:33
% subplot(11,3,j)
% stem(Gold_Seq(j,:))
% axis([1 32 0 1.5])
% end
[r,c] = size(Gold_Seq);
for t = 1:c
    temp1 = find(Gold_Seq(:,t)==0);
    Gold_Seq(temp1,t)=-1;
end

end
