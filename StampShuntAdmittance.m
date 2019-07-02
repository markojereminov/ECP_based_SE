function [Yi_vect,Yj_vect,Yk_vect] = StampShuntAdmittance(Yi_vect,Yj_vect,Yk_vect,G_shunt,B_shunt,Node,CKT_Node_Map)
% FUNCTION DESCRIPTION:
%A function that stamps a linear shunt admittnces for a given vector of Node indices.
%___________________________________________________________________________________________________    
% INPUT: 
    % Yi_vect: vector of prior i indices of the NR Sensitivity matrix
    % Yj_vect: vector of prior j indices of the NR Sensitivity matrix
    % Yk_vect: vector of prior k values of the NR Sensitivity matrix
    % G_shunt: vector of series conductances related to the stamped model
    % B_shunt: vector of series susceptances related to the stamped model
    % Node: Node indices
    % CKT_Node_Map: a node map of the ECP circuit
%___________________________________________________________________________________________________    
% OUTPUT:
    % Yi_vect: vector of concatinated i indices of the NR Sensitivity matrix 
    % Yj_vect: vector of concatinated j indices of the NR Sensitivity matrix
    % Yk_vect: vector of concatinated k values of the NR Sensitivity matrix
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
%
% Equations stamped for the the series admittance element:
% IR_sh = G_shunt*VRF - B_shunt*VIF
% II_sh = G_shunt*VIF + B_shunt*VRF
% adjIR_sh = G_shunt*LRF + B_shunt*LIF
% adjII_sh = G_shunt*LIF - B_shunt*LRF

% Obtaining the real and imaginary circuit nodes from the map related to FromNode
NodeReal = CKT_Node_Map.Bus.NR(Node);
NodeImag = CKT_Node_Map.Bus.NI(Node);
NodeRealAdj = CKT_Node_Map.Bus.LNR(Node);
NodeImagAdj = CKT_Node_Map.Bus.LNI(Node);

%Appending the G_series stamp to the triplet vector definition of Sensitivity matrix:
if ~all(G_shunt==0)
    % Real ckt:
    Yi_vect = [Yi_vect;NodeReal];
    Yj_vect = [Yj_vect;NodeReal];
    Yk_vect = [Yk_vect;G_shunt];
    % Imag ckt:
    Yi_vect = [Yi_vect;NodeImag];
    Yj_vect = [Yj_vect;NodeImag];
    Yk_vect = [Yk_vect;G_shunt];    
    % AdjReal ckt:
    Yi_vect = [Yi_vect;NodeRealAdj];
    Yj_vect = [Yj_vect;NodeRealAdj];
    Yk_vect = [Yk_vect;G_shunt];
    % AdjImag ckt:
    Yi_vect = [Yi_vect;NodeImagAdj];
    Yj_vect = [Yj_vect;NodeImagAdj];
    Yk_vect = [Yk_vect;G_shunt];       
end

%Appending the B_series stamp to the triplet vector definition of Sensitivity matrix:
if ~all(B_shunt==0)
    % Real ckt:
    Yi_vect = [Yi_vect;NodeReal];
    Yj_vect = [Yj_vect;NodeImag];
    Yk_vect = [Yk_vect;-B_shunt];
    % Imag ckt:
    Yi_vect = [Yi_vect;NodeImag];
    Yj_vect = [Yj_vect;NodeReal];
    Yk_vect = [Yk_vect;B_shunt];    
    % AdjReal ckt:
    Yi_vect = [Yi_vect;NodeRealAdj];
    Yj_vect = [Yj_vect;NodeImagAdj];
    Yk_vect = [Yk_vect;B_shunt];
    % AdjImag ckt:
    Yi_vect = [Yi_vect;NodeImagAdj];
    Yj_vect = [Yj_vect;NodeRealAdj];
    Yk_vect = [Yk_vect;-B_shunt];   
end
end