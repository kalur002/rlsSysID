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
lamda = 1;
del = 1e2;
P = del*eye(sysorder);
w = zeros(sysorder,1);

for ii = 3:length(ytilde)
    outp = ytilde(ii-1:-1:ii-sysorder,:);
    H = [outp];
    phi = H'*P;
    Kg = phi'/(lamda+phi*H);
    est(ii,:) = w'*H;
    e(ii,:) = ytilde(ii,:)-est(ii,:);
    w = w+Kg*e(ii,:)';
    P = (P-Kg*phi)/lamda;
    Recordedw(1:sysorder,ii)=w;
end

% for ii = (N+1):length(ytilde)
%     outp = ytilde(ii-1:-1:ii-sysorder,:);
%     H = [outp];
%     est(ii,:) = w'*H;
%     e(ii,:) = ytilde(ii,:) - est(ii,:);
% end


figure()
subplot(2,1,1)
plot(ytilde(:,1),'r--','linewidth',2)
hold on
plot(est(:,1),'k','linewidth',1)
title('Original Posn  VS Estimated Posn')
legend('Original Signal','Estimated Signal')

subplot(2,1,2)
plot(ytilde(:,2),'b--','linewidth',2)
hold on
plot(est(:,2),'r','linewidth',1)
title('Original Vel  VS Estimated Vel')
legend('Original Signal','Estimated Signal')


figure()
subplot(2,1,1)
plot(t,e(:,1))
title('Estimation Posn Error Vs Time')
subplot(2,1,2)
plot(t,e(:,2))
title('Estimation Vel Error Vs Time')

%% Estimated A matrix

A_est = [0 1;w(2) w(1)]
eval_est = eig(A_est)

