% RankUpdateDMD incrementally perform rank one update on dynamic mode
% decomposition
%
% usage: rdmd = RankUpdateDMD
%
%
% RankUpdateDMD has one method
%
%
% Author : Aniketh Kalur
% Date   : 29/6/2017

classdef RankUpdateDMd
    
    properties
        streamIndex     % index to start streaming
        initIndex = 2   % index to initialize the process from
        max_rank = 0    % max rank
        xk              % initial measurements in x
        yk              % initial measurements in y
        xnew            % incoming stream of x measurements
        ynew            % incoming stream of y measurements
        M               % Matrix to store yx'
        G               % Matrix to store the sherman-morrison update
        Phi             % Matrix to store the second term of sherman morrison formula
        A               % System Matrix
        K               % Matrix to store the update of x*x'
    end
    
    methods
        function obj = RankUpdateDMD(max_rank,posn2initialize)
            % constructor created to take input as max rank with the follwoing
            % properties
            if nargin==2
                obj.max_Rank = max_rank;
                obj.initIndex = posn2initialize;
            end
        end
        
        function obj = Initialize(obj,ytilde)
            % Intialize the K matrix to avoid rank deficiency to compute
            % the inverse using sherman morrison formula. This step
            % computes the A,K,M and G matrix using standard DMD algortihm
            % upto the measurement snapshot number mentioned in iniIndex
            obj.xk = ytilde(1:obj.initIndex,:)';
            obj.yk = ytilde(2:obj.initIndex+1,:)';
            obj.K = obj.xk*obj.xk';
            obj.G = inv(obj.K);
            obj.M = obj.yk*obj.xk';
            obj.A = obj.M*obj.G;
        end
        
        function obj = update(obj,ytilde)
            % Updates the DMD algorithm on the fly as and when new
            % measurements are available. The streaming starts from the
            % measurement number mentioned in streamIndex
            obj.streamIndex = obj.initIndex+2;
            obj.ynew = ytilde(obj.streamIndex,:)';
            
            for ii = obj.streamIndex:length(ytilde)-1
                obj.xnew = obj.ynew;
                obj.ynew = ytilde(ii,:)';
                obj.M = obj.M+obj.ynew*obj.xnew';
                obj.Phi = obj.G*obj.xnew*obj.xnew'*obj.G/(1+obj.xnew'*obj.G*obj.xnew);
                obj.A = obj.A + obj.ynew*obj.xnew'*obj.G-obj.M*obj.Phi;
                obj.G = obj.G - obj.Phi;
            end
        end
        
        function [eval,evecs] = computeEval(obj)
            [evecs,eval] = eig(obj.A);
        end
    end
end

