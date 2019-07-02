classdef TxLine < handle
%% Class that defines the Tx Line device
%% This class will contain the properties required for the respective equivalent circuit stamps as well as the following functions:
% TxLine(.)=>the function that initializes the properties of the class
%       INPUT: TxLineData Tx Line data from the parsed input file
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
%% Properties of the TxLine class:
    properties
        from_bus % from bus index
        to_bus % to bus index
        Gs % series conductance of the Tx Line pi-model
        Bs % series susceptance of the Tx Line pi-model
        Bsh % shunt susceptance of the Tx Line pi-model
    end
%% Methods of the TxLine class:
    methods
        function this = TxLine(TxLineData)
           this.from_bus=TxLineData(:,1);
           this.to_bus=TxLineData(:,2);
           this.Gs=TxLineData(:,3);
           this.Bs=TxLineData(:,4);
           this.Bsh=TxLineData(:,5);
        end
        
        function [Yi_vect,Yj_vect,Yk_vect]=stamp_linear(this,Yi_vect,Yj_vect,Yk_vect,CKT_Node_Map)
            % Stamping a series admittance connecting the respective nodes that correpond to the from and to buses of the tx line branch model:
            [Yi_vect,Yj_vect,Yk_vect] = StampSeriesAdmittance(Yi_vect,Yj_vect,Yk_vect,this.Gs,this.Bs,this.from_bus,this.to_bus,CKT_Node_Map);
            % Stamping a shunt admittance of the tx line branch model at the circuit nodes that correspond to bus from:
            [Yi_vect,Yj_vect,Yk_vect] = StampShuntAdmittance(Yi_vect,Yj_vect,Yk_vect,0,0.5*this.Bsh,this.from_bus,CKT_Node_Map);
            % Stamping a shunt admittance of the tx line branch model at the circuit nodes that correspond to bus to:
            [Yi_vect,Yj_vect,Yk_vect] = StampShuntAdmittance(Yi_vect,Yj_vect,Yk_vect,0,0.5*this.Bsh,this.to_bus,CKT_Node_Map);
        end   
    end    
end