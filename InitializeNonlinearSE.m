function [Sol] = InitializeNonlinearSE(varCount,LinearStateEstimate,RTUClass,PMUClass,CKT_Node_Map,n_elem,epsilon)
% FUNCTION DESCRIPTION:
% Function that initializes the nonlinear SE solution vector using the previously obtained linear State Estimate 
%% References:
% [1] M. Jereminov, A. Jovicic, M. Wagner, G. Hug, L. Pileggi, ?Equivalent Circuit Programming for Estimating
%     the State of a Power System,? in Proc. IEEE PowerTech Milan, June 2019.
%___________________________________________________________________________________________________    
% INPUT:
    % varCount: variable count
    % LinearStateEstimate: Linear State Estimate
    % RTUClass: RTU measurement device class
    % PMUClass: RTU measurement device Class
    % CKT_Node_Map: a node map of the ECP circuit
    % n_elem: structure that defines number of element of the respective system
    % epsilon: complementary slackness approximation (diode coefficient) see [1]
%___________________________________________________________________________________________________    
% OUTPUT:
    % Sol: initialized solution vector of primal and adjoint (dual) circuit variables
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

% Allocate the memory space for the solution vector:
Sol = zeros(varCount,1);
Sol(1:4*n_elem.Bus) = LinearStateEstimate; % Initialize the part of the Sol from Linear State Estimate

% Bounds on PMU currents and voltages:
IR_PMU_MIN = PMUClass.IR_min;
IR_PMU_MAX = PMUClass.IR_max;
II_PMU_MIN = PMUClass.II_min;
II_PMU_MAX = PMUClass.II_max;
VR_PMU_MIN = PMUClass.VR_min;
VR_PMU_MAX = PMUClass.VR_max;
VI_PMU_MIN = PMUClass.VI_min;
VI_PMU_MAX = PMUClass.VI_max;

% Computing mean values of PMU currents and voltages:
IR_PMU_mean =  0.5*(IR_PMU_MIN+IR_PMU_MAX);
II_PMU_mean =  0.5*(II_PMU_MIN+II_PMU_MAX);
VR_PMU_mean =  0.5*(VR_PMU_MIN+VR_PMU_MAX);
VI_PMU_mean =  0.5*(VI_PMU_MIN+VI_PMU_MAX);
% Initialization:
Sol(CKT_Node_Map.PMU.IR) = IR_PMU_mean; %initializing IR
Sol(CKT_Node_Map.PMU.II) = II_PMU_mean; %initializing II
Sol(CKT_Node_Map.PMU.VR) = VR_PMU_mean; %initializing IR
Sol(CKT_Node_Map.PMU.VI) = VI_PMU_mean; %initializing II
%Initializing adjoint diode current (Mu) variables that correspond to the bounded PMU model [1]:
Sol(CKT_Node_Map.PMU.MIRmax) = -epsilon./(IR_PMU_mean - IR_PMU_MAX);
Sol(CKT_Node_Map.PMU.MIRmin) = -epsilon./(IR_PMU_MIN - IR_PMU_mean);
Sol(CKT_Node_Map.PMU.MIImax) = -epsilon./(II_PMU_mean - II_PMU_MAX);
Sol(CKT_Node_Map.PMU.MIImin) = -epsilon./(II_PMU_MIN - II_PMU_mean); 
Sol(CKT_Node_Map.PMU.MVRmax) = -epsilon./(VR_PMU_mean - VR_PMU_MAX);
Sol(CKT_Node_Map.PMU.MVRmin) = -epsilon./(VR_PMU_MIN - VR_PMU_mean);
Sol(CKT_Node_Map.PMU.MVImax) = -epsilon./(VI_PMU_mean - VI_PMU_MAX);
Sol(CKT_Node_Map.PMU.MVImin) = -epsilon./(VI_PMU_MIN - VI_PMU_mean); 
%Initializing adjoint diode current (Mu) variables that correspond to the bounded PMU model [1] to zero if bounds are the same:
Sol(CKT_Node_Map.PMU.MIRmax((IR_PMU_MIN-IR_PMU_MAX)==0)) = 0;
Sol(CKT_Node_Map.PMU.MIRmin((IR_PMU_MIN-IR_PMU_MAX)==0)) = 0;
Sol(CKT_Node_Map.PMU.MIImax((II_PMU_MIN-II_PMU_MAX)==0)) = 0;
Sol(CKT_Node_Map.PMU.MIImin((II_PMU_MIN-II_PMU_MAX)==0)) = 0; 
Sol(CKT_Node_Map.PMU.MVRmax((VR_PMU_MIN-VR_PMU_MAX)==0)) = 0;
Sol(CKT_Node_Map.PMU.MVRmin((VR_PMU_MIN-VR_PMU_MAX)==0)) = 0;
Sol(CKT_Node_Map.PMU.MVImax((VI_PMU_MIN-VI_PMU_MAX)==0)) = 0;
Sol(CKT_Node_Map.PMU.MVImin((VI_PMU_MIN-VI_PMU_MAX)==0)) = 0; 

% Obtaining the mean values of RTU GB values:
G_RTU_MIN = RTUClass.G_min;
G_RTU_MAX = RTUClass.G_max;
B_RTU_MIN = RTUClass.B_min;
B_RTU_MAX = RTUClass.B_max;

G_M = RTUClass.G_mean;% mean RTU condictance
B_M = RTUClass.B_mean;% mean RTU susceptance

Sol(CKT_Node_Map.RTU.deltaG) = G_M; % initialize G RTU
Sol(CKT_Node_Map.RTU.deltaB) = B_M; % initialize B RTU
%Initializing adjoint diode current (Mu) variables that correspond to the RTU model from [1]:
Sol(CKT_Node_Map.RTU.MGmax) = -epsilon./(G_M - G_RTU_MAX);
Sol(CKT_Node_Map.RTU.MGmin) = -epsilon./(G_RTU_MIN - G_M);
Sol(CKT_Node_Map.RTU.MBmax) = -epsilon./(B_M - B_RTU_MAX);
Sol(CKT_Node_Map.RTU.MBmin) = -epsilon./(B_RTU_MIN - B_M);
end