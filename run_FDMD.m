clear all
close all 
clc


%% Generate data for Simple Harmoinc Oscillator
M = 5;         % mass in KG
K = 0.2;       % spring constant
Fext = 1;     % external force magnitude
tf = 1000;      % time span
gMode = 1;     % graphic mode
sigma = 0.04;    % standard deviation
[t,y,ytilde,u,X] = SimpleHarmonicOsc(M,K,Fext,tf,sigma,gMode);

%% Forgetting DMD
q = 10;
Fdmd = ForgettingDMD(0,q,1)
Fdmd = Fdmd.Initialize(ytilde)

for ii = q+1:length(ytilde)-1
Fdmd = Fdmd.update(ytilde(ii,:),ytilde(ii+1,:))
[eval,evecs] = Fdmd.computeEval;
end