 classdef PhaseShiftXfmr < handle
%% Class that defines the Phase Shifting xfmr device
%% This class will contain the properties required for the respective equivalent circuit stamps as well as the following functions:
% PhaseShiftXfmr(.)=>the function that initializes the properties of the class
%       INPUT: ShiftXfmrData phase shifting xfmr data from the parsed input file
% stamp_linear(.)=> the function that stamps an equvalent circuit of a respective device
%       INPUT:
%           Yi_vect: vector of prior i indices of the NR Sensitivity matrix
%           Yj_vect: vector of prior j indices of the NR Sensitivity matrix
%           Yk_vect: vector of prior k values of the NR Sensitivity matrix
%           CKT_Node_Map: a node map of the ECP circuit
%       OUTPUT:
%           Yi_vect: vector of concatinated i indices of the NR Sensitivity matrix 
%           Yj_vect: vector of concatinated j indices of the NR Sensitivity matrix
%           Yk_vect: vector of concatinated k values of the NR Sensitivity matrix
%% AUTHOR: 
%         Marko Jereminov
%         m.jereminov92@gmail.com
%         Carnegie Mellon University
%         Department of Electrical and Computer Engineering
%         Pittsburgh, PA
%         United States
%___________________________________________________________________________________________________ 
%% LICENSE:
%   This file is part of open source version of ECP based Static State Estimator.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%% Properties of the PhaseShiftXfmr class:
    properties
        from_bus % from bus index
        to_bus % to bus index
        G_from % conductance at the from side of the phase shifting xfmr
        B_from % susceptance at the from side of the phase shifting xfmr
        MR_from % real gain of the VCCS at the from side defined by the phase shifting xfmr model
        MI_from % imag gain of the VCCS at the from side defined by the phase shifting xfmr model
        MR_to % real gain of the VCCS at the to side defined by the phase shifting xfmr model
        MI_to % imag gain of the VCCS at the to side defined by the phase shifting xfmr model
        G_to % conductance at the to side of the phase shifting xfmr
        B_to % susceptance at the to side of the phase shifting xfmr
    end
%% Methods of the PhaseShiftXfmr class:
    methods
        function this = PhaseShiftXfmr(ShiftXfmrData)
           this.from_bus=ShiftXfmrData(:,1);
           this.to_bus=ShiftXfmrData(:,2);
           this.G_from=ShiftXfmrData(:,3);
           this.B_from=ShiftXfmrData(:,4);
           this.MR_from=ShiftXfmrData(:,5);
           this.MI_from=ShiftXfmrData(:,6); 
           this.MR_to=ShiftXfmrData(:,7); 
           this.MI_to=ShiftXfmrData(:,8); 
           this.G_to=ShiftXfmrData(:,9);
           this.B_to=ShiftXfmrData(:,10);                 
        end

        function [Yi_vect,Yj_vect,Yk_vect]=stamp_linear(this,Yi_vect,Yj_vect,Yk_vect,CKT_Node_Map)
            % Stamping a shunt admittance of the phase shifting xfmr model at the circuit nodes that correspond to bus from:
            [Yi_vect,Yj_vect,Yk_vect] = StampShuntAdmittance(Yi_vect,Yj_vect,Yk_vect,this.G_from,this.B_from,this.from_bus,CKT_Node_Map);
            % Stamping a shunt admittance of the phase shifting xfmr model at the circuit nodes that correspond to bus to:
            [Yi_vect,Yj_vect,Yk_vect] = StampShuntAdmittance(Yi_vect,Yj_vect,Yk_vect,this.G_to,this.B_to,this.to_bus,CKT_Node_Map);
            % Stamping a voltage controlled current source (VCCS) of the phase shifting xfmr model at the circuit nodes that correspond to bus from:
            [Yi_vect,Yj_vect,Yk_vect] = StampVCCS(Yi_vect,Yj_vect,Yk_vect,this.MR_from,this.MI_from,this.from_bus,this.to_bus,CKT_Node_Map);
            % Stamping a voltage controlled current source (VCCS) of the phase shifting xfmr model at the circuit nodes that correspond to bus to:
            [Yi_vect,Yj_vect,Yk_vect] = StampVCCS(Yi_vect,Yj_vect,Yk_vect,this.MR_to,this.MI_to,this.to_bus,this.from_bus,CKT_Node_Map);
        end
    end
end