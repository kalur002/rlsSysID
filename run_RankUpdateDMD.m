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

%% run rank update DMD

rdmd = RankUpdateDMd

rdmd = rdmd.Initialize(ytilde)

rdmd = rdmd.update(ytilde)

[eval,evecs] = rdmd.computeEval