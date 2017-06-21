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

% figure()
% plot(ytilde(:,1),ytilde(:,2))


%% RLS for Parameter estimation
% y(t) = -a1*y(t-Ts)-a2*y(t-2*Ts); here Ts is sampling time
% y(t) = -a1*y(t-1)-a2*y(t-2); where Ts = 1
% y = H*x; where H = [y(t-1) y(t-2)] and x = [-a1 -a2]'
% Now use RLS to estimate x which is xhat(k+1) = xhat(k) + K(k+1)*(ytilde(k+1)-H(k+1)*xhat(k)) 

% p = 2; % order of the system
% lambda = 1;
% invlambda = 1/lambda;
% xhat(:,:,1) = [0.25;1]
% %xhat(:,:,2) = zeros(p,1);
% H = [ytilde(2,1) ytilde(1,1)];
% %H(:,:,1) = zeros(p);
% %H(:,:,2) = zeros(p);
% %H = [y(1,:)' y(2,:)'];
% %H = [u(2) ytilde(2,:) ytilde(1,:)]
% P(:,:,1) = eye(p,p)*100;
% %P(:,:,2) = eye(p,p)*100;%[10^6 1.038*10^3;1.038*10^3 1.077*10^0];%eye(3,3)%;%zeros(p,p);
% Kg(:,:,1)  = P(:,:,1)*H(:,:,1)'*pinv(lambda+H(:,:,1)*P(:,:,1)*H(:,:,1)');
% %Kg(:,:,2)  = P(:,:,1)*H(:,:,2)'*pinv(lambda+H(:,:,2)*P(:,:,1)*H(:,:,2)');
% for i = 1:(length(ytilde)-2)
% H(:,:,i+1) = [ytilde(i+2,1) ytilde(i+1,1)];
% %H(:,:,i+1) = [ytilde(i,:)' ytilde(i-1,:)'];
% Kg(:,:,i+1) = P(:,:,i)*H(:,:,i+1)'*pinv(lambda+H(:,:,i+1)*P(:,:,i)*H(:,:,i+1)');
% P(:,:,i+1) = invlambda*(P(:,:,i)-Kg(:,:,i+1)*H(:,:,i+1)*P(:,:,i));
% xhat(:,:,i+1) = xhat(:,:,i)+ Kg(:,:,i+1)*(ytilde(i+1,1)-H(:,:,i+1)*xhat(:,:,i));
% end
% 
% A_est = [0 1;xhat(2,1,end) xhat(1,1,end)]
% eval_est = eig(A_est)

%% using matlab RLS filter

% P0 = 10*eye(2);
% lam = 1;
% d = ytilde;
% rlsf1 = dsp.RLSFilter(2,'ForgettingFactor',1)
% w = rlsf1.Coefficients;



%% Parameter estimation using function call
% lambda = 1; % forgetting factor
% delta = 1;
% %ytilde = horzcat(ytilde,u);
% 
% [H,e,xhat] = RLSfilta(lambda,ytilde);
% 
% % figure()
% % %for i = 1:(length(ytilde))
% % plot(t,xhat(1,:))
% % hold on
% % %end
% % title('Xhat(1) value over time')
% 
% 
% figure()
% plot(t,e(:,1))
% title('Error over time')
% 
% 
% ynew = H*xhat;
% 
% 
% figure()
% plot(t,ynew(1,:))
% hold on
% plot(t,ytilde(:,1))
% legend('est posn','true posn')
% title('Estimated signal vs Measured Signal')
% 
% % figure
% % for i=1:1001
% %     plot(ynew(:,1001))
% %     hold on
% % end

%% recursive least sqaures
% what = recursiveParamLS(ytilde)

%% DMD

Xdmd = ytilde(1:end-1,:)';
Ydmd = ytilde(2:end,:)';
Atilde = DMD_AK(Xdmd,Ydmd);
[V,Lam] = eig(Atilde);

%% actual

Ad = [0.9801 0.9933;-0.0397 0.9801];
eval_Act = eig(Ad)

%% least squares

A_est1 = Ydmd*pinv(Xdmd)
eval_pls = eig(A_est1)

B = [0;1/M];
C = eye(2,2);
D = 0;
t0 = 0;
dt = 1;
%tf = 100;
tspan = [t0:dt:100]';
x0 = [1;.5];
u = [0;zeros(length(tspan(2:end)),1)];
sys_lse = ss(A_est1,B,C,D);

%% least squares

H = [ytilde(1:end-1,1) ytilde(1:end-1,2)];
xhat1 = pinv(H'*H)*H'*ytilde(2:end,:)

eval_bls = eig(xhat1)