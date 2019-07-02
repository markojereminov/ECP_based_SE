function [i_NR,j_NR,k_NR,J_NR_i,J_NR_k] = StampLinearized(CKT_Node_Map,RTUClass,PMUClass,epsilon,Sol)
% FUNCTION DESCRIPTION:
% Function that stamps the linearized measurement device equvalent circuit elements which include PMU bounds and RTU devices:
%___________________________________________________________________________________________________    
% INPUT: 
    % CKT_Node_Map: a node map of the ECP circuit
    % RTUClass: RTU device class
    % PMUClass: PMU device class
    % epsilon: diode parameter
    % Sol: Solution Vector
%___________________________________________________________________________________________________    
% OUTPUT:
    % i_NR: vector of concatinated i indices of the NR Sensitivity matrix that correspond to the linearized network elements
    % j_NR: vector of concatinated j indices of the NR Sensitivity matrix that correspond to the linearized network elements
    % k_NR: vector of concatinated k values of the NR Sensitivity matrix that correspond to the linearized network elements
    % J_NR_i: vector of concatinated i indices of the NR constant vector
    % J_NR_k: vector of concatinated respective values of the NR constant vector
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
i_NR = []; %Initial empty i linear vector corresponding to the linearized elements
j_NR = []; %Initial empty j linear vector corresponding to the linearized elements
k_NR = []; %Initial empty k linear vector corresponding to the linearized elements

J_NR_i = [];
J_NR_k = [];
%% Stamping the linearized RTU device:
if ~isempty(RTUClass)
    [i_NR, j_NR,k_NR,J_NR_i,J_NR_k] = RTUClass.stamp_linearized(i_NR,j_NR,k_NR,J_NR_i,J_NR_k,Sol,epsilon,CKT_Node_Map); 
end
%% Stamping the PMU boubds:
[i_NR, j_NR,k_NR,J_NR_i,J_NR_k] = PMUClass.stamp_linearized(i_NR,j_NR,k_NR,J_NR_i,J_NR_k,Sol,epsilon,CKT_Node_Map);   
end