classdef InPhaseXfmr < handle
%% Class that defines the InPhase Transformer device
%% This class will contain the properties required for the resoective equvalent circuit stamps as well as the following functions:
% InPhaseXfmr(.)=>the function that initializes the properties of the class
%       INPUT: XfmrData inPhase xfmr data from the parsed input file
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
%% Properties of the InPhaseXfmr class:
    properties
        from_bus % from bus index
        to_bus % to bus index
        Gseries % series conductance of transformer branch pi-model
        Bseries % series susceptance of transformer branch pi-model 
        Gfrom % conductance from the branch pi-model at from bus
        Bfrom % susceptance from the branch pi-model at from bus
        Gto % conductance from the branch pi-model at to bus
        Bto % susceptance from the branch pi-model at to bus
    end
%% Methods of the InPhaseXfmr class:
    methods
        function this = InPhaseXfmr(XfmrData) 
           this.from_bus=XfmrData(:,1); 
           this.to_bus=XfmrData(:,2);
           this.Gseries=XfmrData(:,3); 
           this.Bseries=XfmrData(:,4);
           this.Gfrom=XfmrData(:,5);
           this.Bfrom=XfmrData(:,6);
           this.Gto=XfmrData(:,7);
           this.Bto=XfmrData(:,8);
        end
        
        function [Yi_vect,Yj_vect,Yk_vect]=stamp_linear(this,Yi_vect,Yj_vect,Yk_vect,CKT_Node_Map)
            % Stamping a series admittance connecting the respective nodes that correpond to the from and to buses of the xfmr branch model:
           [Yi_vect,Yj_vect,Yk_vect] = StampSeriesAdmittance(Yi_vect,Yj_vect,Yk_vect,this.Gseries,this.Bseries,this.from_bus,this.to_bus,CKT_Node_Map);
            % Stamping a shunt admittance of the xfmr branch model at the circuit nodes that correspond to bus from:
           [Yi_vect,Yj_vect,Yk_vect] = StampShuntAdmittance(Yi_vect,Yj_vect,Yk_vect,this.Gfrom,this.Bfrom,this.from_bus,CKT_Node_Map);
            % Stamping a shunt admittance of the xfmr branch model at the circuit nodes that correspond to bus to:
           [Yi_vect,Yj_vect,Yk_vect] = StampShuntAdmittance(Yi_vect,Yj_vect,Yk_vect,this.Gto,this.Bto,this.to_bus,CKT_Node_Map);
        end
    end
end