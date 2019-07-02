function [Sol] = circuit_limiting(Sol,Sol_old,CKT_Node_Map,RTUClass,PMUClass,epsilon,damp_coeff)
% FUNCTION DESCRIPTION:
%A function where the circuit limiting techniques are applied; In the case
%of this version of the open-sourced code, we have implemented only the
%basic version of diode limiting in order to handle the inequality
%constraints that can be introduced by the bounds on PMU and RTU measurements
%___________________________________________________________________________________________________    
% INPUT: 
    % Sol: new solution vector
    % Sol_old: solution vector of the previous iteration
    % CKT_Node_Map: a node map of the ECP circuit
    % RTUClass: RTU device class needed for the information about the bounds
    % PMUClass: PMU device class needed for the information about the bounds
    % epsilon: diode coefficient (used to approximate CS conditions)
    % damp_coeff: diode damping coeff that should be (0,1)
%___________________________________________________________________________________________________    
% OUTPUT:
    % Sol: new solution vector (after application of step limiting techniques)
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
%% RTU devices:
% Obtaining respective node indices:
Gnode = CKT_Node_Map.RTU.deltaG;
Bnode = CKT_Node_Map.RTU.deltaB;
G_MINnode = CKT_Node_Map.RTU.MGmin;
G_MAXnode = CKT_Node_Map.RTU.MGmax;
B_MINnode = CKT_Node_Map.RTU.MBmin;
B_MAXnode = CKT_Node_Map.RTU.MBmax;
% Obtaining the limiting parameters for RTU variables (note that each variable step size is limited separately):
[deltaG,deltaGmax,deltaGmin] = diode_limiting(Sol,Sol_old,Gnode,G_MAXnode,G_MINnode,RTUClass.G_max,RTUClass.G_min,epsilon,damp_coeff);
[deltaB,deltaBmax,deltaBmin] = diode_limiting(Sol,Sol_old,Bnode,B_MAXnode,B_MINnode,RTUClass.B_max,RTUClass.B_min,epsilon,damp_coeff);
% Limiting the NR step of each RTU variable:  
Sol(Gnode) = Sol_old(Gnode) + deltaG;
Sol(Bnode) = Sol_old(Bnode) + deltaB;
Sol(G_MAXnode) = Sol_old(G_MAXnode) + deltaGmax;
Sol(G_MINnode) = Sol_old(G_MINnode) + deltaGmin;
Sol(B_MAXnode) = Sol_old(B_MAXnode) + deltaBmax;
Sol(B_MINnode) = Sol_old(B_MINnode) + deltaBmin;
%% PMU devices:
% Obtaining respective node indices:
IRnode = CKT_Node_Map.PMU.IR;
IInode = CKT_Node_Map.PMU.II;
VRnode = CKT_Node_Map.PMU.VR;
VInode = CKT_Node_Map.PMU.VI;
IR_MINnode = CKT_Node_Map.PMU.MIRmin;
IR_MAXnode = CKT_Node_Map.PMU.MIRmax;
II_MINnode = CKT_Node_Map.PMU.MIImin;
II_MAXnode = CKT_Node_Map.PMU.MIImax; 
VR_MINnode = CKT_Node_Map.PMU.MVRmin;
VR_MAXnode = CKT_Node_Map.PMU.MVRmax;
VI_MINnode = CKT_Node_Map.PMU.MVImin;
VI_MAXnode = CKT_Node_Map.PMU.MVImax; 
 % Obtaining the limiting parameters for PMU variables (note that each variable step size is limited separately):
[deltaIR,deltaIRmax,deltaIRmin] = diode_limiting(Sol,Sol_old,IRnode,IR_MAXnode,IR_MINnode,PMUClass.IR_max,PMUClass.IR_min,epsilon,damp_coeff);
[deltaII,deltaIImax,deltaIImin] = diode_limiting(Sol,Sol_old,IInode,II_MAXnode,II_MINnode,PMUClass.II_max,PMUClass.II_min,epsilon,damp_coeff);
[deltaVR,deltaVRmax,deltaVRmin] = diode_limiting(Sol,Sol_old,VRnode,VR_MAXnode,VR_MINnode,PMUClass.VR_max,PMUClass.VR_min,epsilon,damp_coeff);
[deltaVI,deltaVImax,deltaVImin] = diode_limiting(Sol,Sol_old,VInode,VI_MAXnode,VI_MINnode,PMUClass.VI_max,PMUClass.VI_min,epsilon,damp_coeff);
% Limiting the NR-step of each PMU variable:  
Sol(IRnode) = Sol_old(IRnode) + deltaIR;
Sol(IInode) = Sol_old(IInode) + deltaII;
Sol(VRnode) = Sol_old(VRnode) + deltaVR;
Sol(VInode) = Sol_old(VInode) + deltaVI;
Sol(IR_MAXnode) = Sol_old(IR_MAXnode) + deltaIRmax;
Sol(IR_MINnode) = Sol_old(IR_MINnode) + deltaIRmin;
Sol(II_MAXnode) = Sol_old(II_MAXnode) + deltaIImax;
Sol(II_MINnode) = Sol_old(II_MINnode) + deltaIImin;
Sol(VR_MAXnode) = Sol_old(VR_MAXnode) + deltaVRmax;
Sol(VR_MINnode) = Sol_old(VR_MINnode) + deltaVRmin;
Sol(VI_MAXnode) = Sol_old(VI_MAXnode) + deltaVImax;
Sol(VI_MINnode) = Sol_old(VI_MINnode) + deltaVImin;
end