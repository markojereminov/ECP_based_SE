function [CKT_Node_Map,varCount] = BuildCktNodeMap(n_elem)
% DESCRIPTION: 
% A funtion that creates a node map of the ECP circuit
%__________________________________________________________________________
% INPUT: 
    % n_elem: a structure that defines number of devices within the considered systsem
%___________________________________________________________________________________________________    
% OUTPUT:
    % CKT_Node_Map: Circuit Node Map
    % varCount: variable count
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
% Obtaining the total numbers of devices:
nRTU = n_elem.RTU; %total number of RTUs
nPMU = n_elem.PMU; %total number of PMUs
nBus = n_elem.Bus; % total number of buses

CKT_Node_Map.Bus.NR = (1:1:nBus)'; % Real Node Voltages of the primal circuit
CKT_Node_Map.Bus.NI = ((nBus+1):1:2*nBus)'; % Imag Node Voltages of the primal circuit
CKT_Node_Map.Bus.LNR = ((2*nBus+1):1:3*nBus)'; % Real Node Voltages of the adjoint circuit
CKT_Node_Map.Bus.LNI = ((3*nBus+1):1:4*nBus)';  % Imag Node Voltages of the adjoint circuit

varCount = 4*nBus;

% Generate RTU contribution to the ckt map:
CKT_Node_Map.RTU.deltaG = ((varCount+1):1:(nRTU+varCount))';
CKT_Node_Map.RTU.deltaB = ((varCount+1+nRTU):1:(2*nRTU+varCount))';
CKT_Node_Map.RTU.MGmax = ((varCount+1+2*nRTU):1:(3*nRTU+varCount))';
CKT_Node_Map.RTU.MGmin = ((varCount+1+3*nRTU):1:(4*nRTU+varCount))'; 
CKT_Node_Map.RTU.MBmax = ((varCount+1+4*nRTU):1:(5*nRTU+varCount))'; 
CKT_Node_Map.RTU.MBmin = ((varCount+1+5*nRTU):1:(6*nRTU+varCount))'; 

varCount = varCount + 6*nRTU;

% Generate PMU contribution to the ckt map:
CKT_Node_Map.PMU.VR = ((varCount+1):1:(nPMU+varCount))';
CKT_Node_Map.PMU.VI = ((varCount+1+nPMU):1:(2*nPMU+varCount))'; 
CKT_Node_Map.PMU.IsR = ((varCount+1+2*nPMU):1:(3*nPMU+varCount))';
CKT_Node_Map.PMU.IsI = ((varCount+1+3*nPMU):1:(4*nPMU+varCount))';
CKT_Node_Map.PMU.IR = ((varCount+1+4*nPMU):1:(5*nPMU+varCount))';
CKT_Node_Map.PMU.II = ((varCount+1+5*nPMU):1:(6*nPMU+varCount))'; 
CKT_Node_Map.PMU.LR = ((varCount+1+6*nPMU):1:(7*nPMU+varCount))';
CKT_Node_Map.PMU.LI = ((varCount+1+7*nPMU):1:(8*nPMU+varCount))'; 
CKT_Node_Map.PMU.adjIsR = ((varCount+1+8*nPMU):1:(9*nPMU+varCount))';
CKT_Node_Map.PMU.adjIsI = ((varCount+1+9*nPMU):1:(10*nPMU+varCount))';   
CKT_Node_Map.PMU.MIRmax = ((varCount+1+10*nPMU):1:(11*nPMU+varCount))';
CKT_Node_Map.PMU.MIRmin = ((varCount+1+11*nPMU):1:(12*nPMU+varCount))';
CKT_Node_Map.PMU.MIImax = ((varCount+1+12*nPMU):1:(13*nPMU+varCount))';
CKT_Node_Map.PMU.MIImin = ((varCount+1+13*nPMU):1:(14*nPMU+varCount))';
CKT_Node_Map.PMU.MVRmax = ((varCount+1+14*nPMU):1:(15*nPMU+varCount))';
CKT_Node_Map.PMU.MVRmin = ((varCount+1+15*nPMU):1:(16*nPMU+varCount))';
CKT_Node_Map.PMU.MVImax = ((varCount+1+16*nPMU):1:(17*nPMU+varCount))';
CKT_Node_Map.PMU.MVImin = ((varCount+1+17*nPMU):1:(18*nPMU+varCount))';

varCount = varCount + 18*nPMU;
end