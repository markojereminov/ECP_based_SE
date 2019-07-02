function [Y_linear,J_linear] = build_linear(Yi_lin,Yj_lin,Yk_lin,PMUClass,varCount,CKT_Node_Map)
% FUNCTION DESCRIPTION:
% Function that builds the linear sparse contribution to the Hessian of necessary optimality conditions
%___________________________________________________________________________________________________    
% INPUT:
    % Yi_lin: vector of concatinated i indices of the NR Sensitivity matrix that correspond to the linear network elements
    % Yj_lin: vector of concatinated j indices of the NR Sensitivity matrix that correspond to the linear network elements
    % Yk_lin: vector of concatinated k values of the NR Sensitivity matrix that correspond to the linear network elements
    % PMUClass: RTU measurement device Class    
    % varCount: variable count
    % CKT_Node_Map: a node map of the ECP circuit
%___________________________________________________________________________________________________    
% OUTPUT:
    % Y_linear: Linear Y contributon of the Hessian natrix
    % J_linear: Linear J vector contribution
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
%% Stamp the linear part of the PMU stamps: (bounds will be stamped later in NR function since they introduce nonlinearities)
[Yi_lin,Yj_lin,Yk_lin,Ji_lin,Jk_lin]=PMUClass.stamp_linear(Yi_lin,Yj_lin,Yk_lin,CKT_Node_Map);
%% Building the Linear linear sparse contribution to the Hessian
Y_linear = sparse([Yi_lin;varCount],[Yj_lin;varCount],[Yk_lin;0]);
J_linear = sparse([Ji_lin;varCount],1,[Jk_lin;0]);
end

