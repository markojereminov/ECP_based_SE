function [Yi_lin,Yj_lin,Yk_lin] = StampLinearNetwork(CKT_Node_Map,TxLineClass,XformsClass,ShuntClass,PShiftClass)
% FUNCTION DESCRIPTION:
% Function that stamps the linear network equvalent circuit elements which include Tx Line, fixed shunts and xfmrs:
%___________________________________________________________________________________________________    
% INPUT: 
    % CKT_Node_Map: a node map of the ECP circuit
    % TxLineClass: Tx Line device class
    % XformsClass: Xform devicenClass
    % ShuntClass: Shunt devicenClass
    % PShiftClass: Phase Shifting xfmr class
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
Yi_lin = []; %Initial empty i linear vector corresponding to the network elements
Yj_lin = []; %Initial empty j linear vector corresponding to the network elements
Yk_lin = []; %Initial empty k linear vector corresponding to the network elements
%% Stamping Transmission Lines:
[Yi_lin, Yj_lin,Yk_lin] = TxLineClass.stamp_linear(Yi_lin, Yj_lin,Yk_lin,CKT_Node_Map); 
%% Stamping the transformeres:
if ~isempty(XformsClass)
    [Yi_lin, Yj_lin,Yk_lin] = XformsClass.stamp_linear(Yi_lin, Yj_lin,Yk_lin,CKT_Node_Map); 
end
%% Stamping the shunt elements:
if ~isempty(ShuntClass)
    [Yi_lin, Yj_lin,Yk_lin] = ShuntClass.stamp_linear(Yi_lin, Yj_lin,Yk_lin,CKT_Node_Map);   
end
%% Stamping the phase shifting transformeres:
if ~isempty(PShiftClass)
    [Yi_lin, Yj_lin,Yk_lin] = PShiftClass.stamp_linear(Yi_lin, Yj_lin,Yk_lin,CKT_Node_Map); 
end
end