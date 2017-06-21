%% simple DMD

function [Atilde]= DMD_AK(x,y)

X = x;
Y = y;

[Phi, Sig,W] = svd(X);
diagSig = diag(Sig);

Atilde = Phi'*Y*W/Sig;

%[V,Lam] = eig(Atilde);

