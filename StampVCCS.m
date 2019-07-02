function [Yi_vect,Yj_vect,Yk_vect] = StampVCCS(Yi_vect,Yj_vect,Yk_vect,G_vccs,B_vccs,Node,ControlNode,CKT_Node_Map)
% FUNCTION DESCRIPTION:
%A function that stamps a voltage controlled current source from node to ground for a given vector of Node indices.
%___________________________________________________________________________________________________    
% INPUT: 
    % Yi_vect: vector of prior i indices of the NR Sensitivity matrix
    % Yj_vect: vector of prior j indices of the NR Sensitivity matrix
    % Yk_vect: vector of prior k values of the NR Sensitivity matrix
    % G_vccs: vector of real VCCS paramteres related to the stamped model
    % B_vccs: vector of imag VCCS paramteres related to the stamped model
    % Node: Node indices
    % ControlNode: Node of the controlled voltage
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
% Equations stamped for the the series admittance element:
% IR = G_vccs*VR_control - B_vccs*VI_control;
% II = G_vccs*VI_control + B_vccs*VR_control;
% adjIRcontr = G_vccs*VR + B_vccs*VI;
% adjIIcontr = G_vccs*VI - B_vccs*VR;

% Obtaining the real and imaginary circuit nodes from the map:
NodeReal = CKT_Node_Map.Bus.NR(Node);
NodeImag = CKT_Node_Map.Bus.NI(Node);
NodeRealAdj = CKT_Node_Map.Bus.LNR(Node);
NodeImagAdj = CKT_Node_Map.Bus.LNI(Node);
NodeRealContr = CKT_Node_Map.Bus.NR(ControlNode);
NodeImagContr = CKT_Node_Map.Bus.NI(ControlNode);
NodeRealContrAdj = CKT_Node_Map.Bus.LNR(ControlNode);
NodeImagContrAdj = CKT_Node_Map.Bus.LNI(ControlNode);
%Appending the G_series stamp to the triplet vector definition of Sensitivity matrix:
if ~all(G_vccs==0)
    % Real ckt:
    Yi_vect = [Yi_vect;NodeReal];
    Yj_vect = [Yj_vect;NodeRealContr];
    Yk_vect = [Yk_vect;G_vccs];
    % Imag ckt:
    Yi_vect = [Yi_vect;NodeImag];
    Yj_vect = [Yj_vect;NodeImagContr];
    Yk_vect = [Yk_vect;G_vccs];    
    % AdjReal ckt:
    Yi_vect = [Yi_vect;NodeRealContrAdj];
    Yj_vect = [Yj_vect;NodeRealAdj];
    Yk_vect = [Yk_vect;G_vccs];
    % AdjImag ckt:
    Yi_vect = [Yi_vect;NodeImagContrAdj];
    Yj_vect = [Yj_vect;NodeImagAdj];
    Yk_vect = [Yk_vect;G_vccs];       
end
%Appending the B_series stamp to the triplet vector definition of Sensitivity matrix:
if ~all(B_vccs==0)
    % Real ckt:
    Yi_vect = [Yi_vect;NodeReal];
    Yj_vect = [Yj_vect;NodeImagContr];
    Yk_vect = [Yk_vect;-B_vccs];
    % Imag ckt:
    Yi_vect = [Yi_vect;NodeImag];
    Yj_vect = [Yj_vect;NodeRealContr];
    Yk_vect = [Yk_vect;B_vccs];    
    % AdjReal ckt:
    Yi_vect = [Yi_vect;NodeRealContrAdj];
    Yj_vect = [Yj_vect;NodeImagAdj];
    Yk_vect = [Yk_vect;B_vccs];
    % AdjImag ckt:
    Yi_vect = [Yi_vect;NodeImagContrAdj];
    Yj_vect = [Yj_vect;NodeRealAdj];
    Yk_vect = [Yk_vect;-B_vccs];   
end
end