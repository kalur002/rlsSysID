%% 
function [Tout,Yout,Yn,u,Xout] = SimpleHarmonicOsc(M,K,F_ext,tf,sigma,Graphic)
% Graphic is a binary input to plot the output obtained
% m        mass [units KG]
% K        spring constant
% F_ext    external force

narginchk(5,6)

%% SHM system definition
A= [0 1;-K/M 0];
B = [0;F_ext/M];
C = [1 0;0 1];
D = [0];

t0 = 0;
dt = 1;
%tf = 100;
tspan = [t0:dt:tf]';
x0 = [1;.5];
u = [0;zeros(length(tspan(2:end)),1)];
sys = ss(A,B,C,D);
sysd = c2d(sys,dt,'zoh');
%% Simulate SHM with noise
[Yout,Tout,Xout]=lsim(sysd,u,tspan,x0);
%[Youtd,Toutd] = dlsim(A,B,C,D,u,tspan,x0);
[row,col] = size(Yout);
noise = sigma*rand(row,col);
Yn = Yout+noise;

if Graphic==1
figure()
plot(Tout,Yout)
legend('True Posn','True Vel')
elseif Graphic ==2
figure()
plot(Tout,Yout)
legend('True Posn','True Vel')
hold on
plot(Tout,Yn)
legend('True Posn','True Vel','Noisy Posn','Noisy Vel')
end
title('True vs Measurement')
Ad = expm(A*dt)
eval = eig(Ad)
end