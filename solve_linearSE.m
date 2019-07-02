function [LinearStateEstimate] = solve_linearSE(Yi_lin,Yj_lin,Yk_lin,CKT_Node_Map,RTUClass,PMUClass,n_elem)
% FUNCTION DESCRIPTION:
% Function that stamps the linear unbounded PMU and relaxed (linear) RTU
% models. This function provides the solution to the linear SE problem and is
% further used to initialize the nonlinear SE problem.
%% References:
% More information on linear SE formulation can be found:
% [1] A. Jovicic, M. Jereminov, L. Pileggi, and G. Hug, ?A Linear Formulation for Power System State Estimation 
%     including RTU and PMU Measurements, in Proc. of ISGT Europe 2019, Sept. 2019, Bucharest, Romania?
% [2] Martin R. Wagner, Marko Jereminov, Amritanshu Pandey, Larry Pileggi, "A Probabilistic Approach to Power System 
%    State Estimation using a Linear Algorithm," in Proc. of IEEE/EEEIC Genoa 2019, June 2019, Genoa, Italy
%___________________________________________________________________________________________________    
% INPUT:
    % Yi_lin: vector of concatinated i indices of the NR Sensitivity matrix that correspond to the linear network elements
    % Yj_lin: vector of concatinated j indices of the NR Sensitivity matrix that correspond to the linear network elements
    % Yk_lin: vector of concatinated k values of the NR Sensitivity matrix that correspond to the linear network elements
    % CKT_Node_Map: a node map of the ECP circuit
    % RTUClass: RTU measurement device class
    % PMUClass: RTU measurement device Class
    % n_elem: structure that defines number of element of the respective system
%___________________________________________________________________________________________________    
% OUTPUT:
    % Yi_lin: vector of concatinated i indices of the NR Sensitivity matrix that correspond to the linear network elements
    % Yj_lin: vector of concatinated j indices of the NR Sensitivity matrix that correspond to the linear network elements
    % Yk_lin: vector of concatinated k values of the NR Sensitivity matrix that correspond to the linear network elements
%___________________________________________________________________________________________________    
% AUTHOR: Marko Jereminov
%         m.jereminov92@gmail.com
%         Carnegie Mellon University
%         Department of Electrical and Computer Engineering
%         Pittsburgh, PA
%         United States
%___________________________________________________________________________________________________   
%% LICENSE:
%   This file is part of open source version of ECP based Static State Estimator.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%___________________________________________________________________________________________________  
%% Stamping linear RTU model [1]-[2]:
if ~isempty(RTUClass)
    [Yi_lin,Yj_lin,Yk_lin] = RTUClass.linearSE_stamp(Yi_lin,Yj_lin,Yk_lin,CKT_Node_Map);
end
%% Stamping unbounded PMU model:
[Yi_lin,Yj_lin,Yk_lin,Ji_lin,Jk_lin]=PMUClass.linearSE_stamp(Yi_lin,Yj_lin,Yk_lin,CKT_Node_Map);
%% Building the sparse system equations that represent the neccessary and sufficient optimality conditions of the linear SE problem
% The linearity of the problem guarantees its global optimality
H_LinSE = sparse(Yi_lin,Yj_lin,Yk_lin); % Hessian matrix
J_LinSE = sparse([Ji_lin;4*n_elem.Bus],1,[Jk_lin;0]); % Constant vector
%% Solving for the linear State Estimate:
LinearStateEstimate = H_LinSE\J_LinSE;
end