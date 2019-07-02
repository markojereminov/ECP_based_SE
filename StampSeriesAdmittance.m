function [Yi_vect,Yj_vect,Yk_vect] = StampSeriesAdmittance(Yi_vect,Yj_vect,Yk_vect,G_series,B_series,FromNode,ToNode,CKT_Node_Map)
% FUNCTION DESCRIPTION:
%A function that stamps a linear series admittnces between nodes FROM and TO.
%___________________________________________________________________________________________________    
% INPUT: 
    % Yi_vect: vector of prior i indices of the NR Sensitivity matrix
    % Yj_vect: vector of prior j indices of the NR Sensitivity matrix
    % Yk_vect: vector of prior k values of the NR Sensitivity matrix
    % G_series: vector of series conductances related to the stamped model
    % B_series: vector of series susceptances related to the stamped model
    % FromNode: From Node indices
    % ToNode: To Node indices
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
% IR_from = G_series*(VRF-VRT) - B_series*(VIF-VIT)
% IR_to = G_series*(VRT-VRF) - B_series*(VIT-VIF)
% II_from = G_series*(VIF-VIT) + B_series*(VRF-VRT)
% II_to = G_series*(VIT-VIF) + B_series*(VRT-VRF)
% adjIR_from = G_series*(LRF-LRT) + B_series*(LIF-LIT)
% adjIR_to = G_series*(LRT-LRF) + B_series*(LIT-LIF)
% adjII_from = G_series*(LIF-VLT) - B_series*(LRF-LRT)
% adjII_to = G_series*(LIT-LIF) - B_series*(LRT-LRF)

% Obtaining the real and imaginary circuit nodes from the map related to FromNode
FromReal = CKT_Node_Map.Bus.NR(FromNode);
FromImag = CKT_Node_Map.Bus.NI(FromNode);
FromRealAdj = CKT_Node_Map.Bus.LNR(FromNode);
FromImagAdj = CKT_Node_Map.Bus.LNI(FromNode);
% Obtaining the real and imaginary circuit nodes from the map related to ToNode
ToReal = CKT_Node_Map.Bus.NR(ToNode);
ToImag = CKT_Node_Map.Bus.NI(ToNode);
ToRealAdj = CKT_Node_Map.Bus.LNR(ToNode);
ToImagAdj = CKT_Node_Map.Bus.LNI(ToNode);

%Appending the G_series stamp to the triplet vector definition of Sensitivity matrix:
if ~all(G_series==0)
    % Real ckt:
    Yi_vect = [Yi_vect;FromReal;FromReal;ToReal;ToReal];
    Yj_vect = [Yj_vect;FromReal;ToReal;ToReal;FromReal];
    Yk_vect = [Yk_vect;G_series;-G_series;G_series;-G_series];
    % Imag ckt:
    Yi_vect = [Yi_vect;FromImag;FromImag;ToImag;ToImag];
    Yj_vect = [Yj_vect;FromImag;ToImag;ToImag;FromImag];
    Yk_vect = [Yk_vect;G_series;-G_series;G_series;-G_series];    
    % AdjReal ckt:
    Yi_vect = [Yi_vect;FromRealAdj;FromRealAdj;ToRealAdj;ToRealAdj];
    Yj_vect = [Yj_vect;FromRealAdj;ToRealAdj;ToRealAdj;FromRealAdj];
    Yk_vect = [Yk_vect;G_series;-G_series;G_series;-G_series];
    % AdjImag ckt:
    Yi_vect = [Yi_vect;FromImagAdj;FromImagAdj;ToImagAdj;ToImagAdj];
    Yj_vect = [Yj_vect;FromImagAdj;ToImagAdj;ToImagAdj;FromImagAdj];
    Yk_vect = [Yk_vect;G_series;-G_series;G_series;-G_series];       
end

%Appending the B_series stamp to the triplet vector definition of Sensitivity matrix:
if ~all(B_series==0)
    % Real ckt:
    Yi_vect = [Yi_vect;FromReal;FromReal;ToReal;ToReal];
    Yj_vect = [Yj_vect;FromImag;ToImag;ToImag;FromImag];
    Yk_vect = [Yk_vect;-B_series;B_series;-B_series;B_series];
    % Imag ckt:
    Yi_vect = [Yi_vect;FromImag;FromImag;ToImag;ToImag];
    Yj_vect = [Yj_vect;FromReal;ToReal;ToReal;FromReal];
    Yk_vect = [Yk_vect;B_series;-B_series;B_series;-B_series];    
    % AdjReal ckt:
    Yi_vect = [Yi_vect;FromRealAdj;FromRealAdj;ToRealAdj;ToRealAdj];
    Yj_vect = [Yj_vect;FromImagAdj;ToImagAdj;ToImagAdj;FromImagAdj];
    Yk_vect = [Yk_vect;B_series;-B_series;B_series;-B_series];
    % AdjImag ckt:
    Yi_vect = [Yi_vect;FromImagAdj;FromImagAdj;ToImagAdj;ToImagAdj];
    Yj_vect = [Yj_vect;FromRealAdj;ToRealAdj;ToRealAdj;FromRealAdj];
    Yk_vect = [Yk_vect;-B_series;B_series;-B_series;B_series];   
end
end