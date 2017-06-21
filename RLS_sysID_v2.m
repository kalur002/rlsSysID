%% RunSHM
clear all
close all 
clc

%% Generate data for Simple Harmoinc Oscillator
M = 5;         % mass in KG
K = 0.2;       % spring constant
Fext = 1;     % external force magnitude
tf = 1000;      % time span
gMode = 2;     % graphic mode
sigma = 0.04;    % standard deviation
[t,y,ytilde,u,X] = SimpleHarmonicOsc(M,K,Fext,tf,sigma,gMode);

%% 
% y(t) = -a1*y(t-Ts)-a2*y(t-2*Ts); here Ts is sampling time
% y(t) = -a1*y(t-1)-a2*y(t-2); where Ts = 1
% y = H*x; where H = [y(t-1) y(t-2)] and x = [-a1 -a2]'
% Now use RLS to estimate x which is xhat(k+1) = xhat(k) + K(k+1)*(ytilde(k+1)-H(k+1)*xhat(k)) 

start = tic;
sysorder = 2; % number of output coeff in a 

% take only 50 points for training
N = 50; 
lamda = 0.5;
del = 1e2;
P = del*eye(sysorder);
w = zeros(sysorder,1);

for ii = 3:N
    outp = ytilde(ii-1:-1:ii-sysorder);
    H = [outp]';
    phi = H'*P;
    Kg = phi'/(lamda+phi*H);
    est(ii) = w'*H;
    e(ii) = ytilde(ii,1)-est(ii);
    w = w+Kg*e(ii);
    P = (P-Kg*phi)/lamda;
    Recordedw(1:sysorder,ii)=w;
end

for ii = (N+1):length(ytilde)
    outp = ytilde(ii-1:-1:ii-sysorder);
    H = [outp]';
    est(ii) = w'*H;
    e(ii) = ytilde(ii,1) - est(ii);
end


figure()
plot(ytilde(:,1),'r--','linewidth',2)
hold on
plot(est,'k','linewidth',1)
title('Original  VS Estimated')
legend('Original Signal','Estimated Signal')

figure()
plot(t,e)
title('Estimation Error Vs Time')



%% Eigen Values of reconstructed Matrix

A_est = [0 1;w(2) w(1)]
eval_est = eig(A_est)

A= [0 1;-K/M 0];
B = [0;Fext/M];
C = [1 0;0 1];
D = [0];

sys_est = ss(A_est,B,C,D)
u = zeros(length(t),1);
x0 = [1;.5];
lsim(sys_est,u,t,x0)

eval_act = eig(expm(A))

figure()
plot(eval_est,'kx','linewidth',4)
hold on
ezplot('x^2+Y^2=1')
plot(eval_act,'ro','linewidth',4)
title('Eigen Values of Identified Systems')
legend('Identified System','Unit circle','Actual System')
axis tight

stop = toc;

lamda
total_time = start-stop


% clear all
% close all
% hold off
% % Number of system points
% N=2000;
% inp = randn(N,1);
% n = randn(N,1);
% [b,a] = butter(2,0.25);
% Gz = tf(b,a,-1);
% %y(n)=b1*u(n)+b2*u(n-1)+b3*u(n-2)-a1*y(n-1)-a2*y(n-2)
% h=[b -a(2:length(a))];

