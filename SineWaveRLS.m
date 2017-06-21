clear all
clc
close all

%% generate data

n = 100;                   % num of points
p = 3;                      % filter order
t = linspace(0,1000,n);     % time    
omega  = 1;                % freq
s = sin(2*pi*omega*t)       % sine wave
SNR = -20;
noise = sqrt(0.5)*rand(n,1)
noise = filter(100,1,noise);
%s = s/sqrt(var(s)/(10^(SNR/10)*var(noise)));
%    disp(['Calculated SNR = ' num2str(10*log10(var(s)/var(noise)))])
signal = s+noise';

%% recursive estimation filtering

p = 2 ;                      % filter order
lambda = 1.0  ;              % forgetting factor
del = 1.0
invlambda = 1/lambda;

w = zeros(p,1);
x = signal;
e = x*0;
est_sig = [];
P  = del*eye(p);

% algorithm taken from wikipedia

for ii  = p:n
    y = noise(ii:-1:ii-p+1);
    est_sig(ii) = y'*w;
    e(ii) = x(ii) - y'*w;
    K = P*y*pinv(lambda+y'*P*y);
    P = invlambda*(P-K*y'*P);
    w = w+K*e(ii);
end
    
figure()
plot(t,x)
hold on
plot(t,est_sig)
legend('original','recreated signal')