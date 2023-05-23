%ZadoffChu Sequence

% ZadoffChu sequence of length 16K.
% f = factor(16383);    % finds prime factors of the given number. Here for 16383 we have three prime factors: 3,43,127.
a3 = zeros(1,16383);
for t = 2:16383
    if (mod(t,3) ~= 0) && (mod(t,43) ~= 0) && (mod(t,127) ~= 0)
        if mod(16383,t) ~= 0
            a3(t) = t;
        end
    end
end
b3 = find(a3 ~= 0);

re = zeros(16383,numel(b3));
for zd = 1:numel(b3)
    re(:,zd) = real(zadoffChuSeq(b3(zd),16383));
end
re(16384,:) = -0.01;

% ZadoffChu sequence of length 1K.
% a = zeros(1,1023);
% %f = factor(1023);
% for i =2:1023
%     if mod(i,3) ~= 0 && mod(i,11) ~= 0 && mod(i,31) ~= 0
%         if mod(1023,i) ~= 0
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
% re = zeros(1023,numel(b));
% for zd = 1:numel(b)
%     re(:,zd) = real(zadoffChuSeq(b(zd),1023));
% end
% re(1024,:) = -0.01;

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