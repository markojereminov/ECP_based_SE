classdef FixedShunt < handle
%% Class that defines the fixed shunt device
%% This class will contain properties required for the resoective equvalent circuit stamps as well as the following functions:
% FixedShunt(.)=>the function that initializes the properties of the class
%       INPUT: bus_data: shunt bus number
%              G_data: shunt conductance 
%              B_data: shunt susceptance 
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
%% Properties of the FixedShunt class:
    properties
        bus_num % shunt bus index
        Gsh % shunt conductance
        Bsh % shunt susceptance
    end
%% Methods of the Fixed Shunt class:
    methods
        function this = FixedShunt(bus_data,G_data,B_data)
           this.bus_num=bus_data;
           this.Gsh=G_data;
           this.Bsh=B_data;
        end
        
        function [Yi_vect,Yj_vect,Yk_vect]=stamp_linear(this,Yi_vect,Yj_vect,Yk_vect,CKT_Node_Map)
            % Stamping a shunt admittance of the shunt model at the circuit nodes that correspond to respective bus:
            [Yi_vect,Yj_vect,Yk_vect] = StampShuntAdmittance(Yi_vect,Yj_vect,Yk_vect,this.Gsh,this.Bsh,this.bus_num,CKT_Node_Map);
        end
    end
end