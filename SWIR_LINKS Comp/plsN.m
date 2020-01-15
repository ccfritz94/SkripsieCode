function [B,T,TT,P,W,Q] = plsN(X,Y,A)
%NIPALS PLS

for a = 1:A,
v = X'*Y;
W(:,a) = v/sqrt(v'*v);
T(:,a) = X*W(:,a);
TT = T(:,a)'*T(:,a);
P(:,a) = X'*T(:,a)/TT;
X = X-T(:,a)*P(:,a)';
Q(a,1) = T(:,a)'*Y/TT;
B(:,a) = W*inv(P'*W)*Q;
end

    
