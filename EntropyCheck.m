clc
clear

X = randn(1000,1);
X=[4 14 16 24 36];
enp=entropy_single(X);


X1=[2 7 8 12 18];
enp1=entropy_single(X1);

H_default = Entropy_Array(X);
H_fd = Entropy_Array(X,'fd');



X=[51 120 200 500];
H_default1 = Entropy_Array(X);
H_fd1 = Entropy_Array(X,'fd');

t=0;

function [en]=entropy_single(Z)
%p = hist(Z,sqrt(length(Z)));
%p = hist(Z,unique(Z));
p = hist(Z,length(Z));
% remove zero entries in p
p(p==0) = [];
% normalize p so that sum(p) is one.
p = p ./ numel(Z);
en = -sum(p.*log(p));
end

