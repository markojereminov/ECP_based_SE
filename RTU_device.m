classdef RTU_device < handle
%% Class that defines the RTU measurement device
%% This class will contain the properties required for the respective equivalent circuit stamps as well as the following functions:
% RTU_device(.)=>the function that initializes the properties of the class
%       INPUT: RTU_weight: vector of weights for the RTU measurements (see [1]-[2])
%              bus_data: RTU device bus number
%              data_G: bounds on G-rtu variable
%              data_B: bounds on B-rtu variable
%              data_Gm: mean value of G-rtu variable
%              data_Bm: mean value of B-rtu variable
% linearSE_stamp(.)=> the function that stamps an equvalent circuit of a respective device for the linear SE formulation (relaxed linear RTU [2]-[3])
%        INPUT: Yi_lin: vector of prior i indices of the NR Sensitivity matrix
%               Yj_lin: vector of prior j indices of the NR Sensitivity matrix
%               Yk_lin: vector of prior k values of the NR Sensitivity matrix
%               CKT_Node_Map: a node map of the ECP circuit
%        OUTPUT:Yi_lin: vector of concatinated i indices of the NR Sensitivity matrix 
%               Yj_lin: vector of concatinated j indices of the NR Sensitivity matrix
%               Yk_lin: vector of concatinated k values of the NR Sensitivity matrix
% stamp_linearized(.) => the function that stamps an equvalent circuit of respective RTU devices
%        INPUT: Yi_lz: vector of prior i indices of the NR Sensitivity matrix
%               Yj_lz: vector of prior j indices of the NR Sensitivity matrix
%               Yk_lz: vector of prior k values of the NR Sensitivity matrix
%               Ji_lz: vector of concatinated i indices of the NR constant vector
%               Jk_lz: vector of concatinated respective values of the NR constant vector
%               Sol: Solution Vector
%               epsilon: diode parameter for modeling the CS equations
%               CKT_Node_Map: a node map of the ECP circuit
%        OUTPUT:Yi_lz: vector of concatinated i indices of the NR Sensitivity matrix 
%               Yj_lz: vector of concatinated j indices of the NR Sensitivity matrix
%               Yk_lz: vector of concatinated k values of the NR Sensitivity matrix
%               Ji_lz: vector of concatinated i indices of the NR constant vector
%               Jk_lz: vector of concatinated respective values of the NR constant vector
%% AUTHOR: 
%         Marko Jereminov
%         m.jereminov92@gmail.com
%         Carnegie Mellon University
%         Department of Electrical and Computer Engineering
%         Pittsburgh, PA
%         United States
%% References:
% [1] M. Jereminov, A. Jovicic, M. Wagner, G. Hug, L. Pileggi, Equivalent Circuit Programming for Estimating
%     the State of a Power System, in Proc. IEEE PowerTech Milan, June 2019.
% [2] A. Jovicic, M. Jereminov, L. Pileggi, and G. Hug, A Linear Formulation for Power System State Estimation 
%     including RTU and PMU Measurements, in Proc. of ISGT Europe 2019, Sept. 2019, Bucharest, Romania?
% [3] Martin R. Wagner, Marko Jereminov, Amritanshu Pandey, Larry Pileggi, "A Probabilistic Approach to Power System 
%    State Estimation using a Linear Algorithm," in Proc. of IEEE/EEEIC Genoa 2019, June 2019, Genoa, Italy
%___________________________________________________________________________________________________ 
%% LICENSE:
%   This file is part of open source version of ECP based Static State Estimator.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%% Properties of the RTU measurement class:
    properties
        RTU_weight
        busNum
        G_max % MAXIMUM G RTU measurement
        G_min % MINIMUM G RTU measurement
        B_max % MAXIMUM B RTU measurement 
        B_min % MINIMUM B RTU measurement
        G_mean %mean G value
        B_mean %mean B value
    end
%% Methods of the RTU measurement class:
    methods
        function this = RTU_device(weight,bus_data,data_G, data_B, data_Gm, data_Bm) %initialize
           this.RTU_weight = weight;
           this.busNum = bus_data;
           this.G_max=data_G(:,1);  
           this.G_min=data_G(:,2);  
           this.B_max=data_B(:,1);  
           this.B_min=data_B(:,2); 
           this.G_mean = data_Gm(:,1);
           this.B_mean = data_Bm(:,1);
        end        
        function [Yi_lin,Yj_lin,Yk_lin]=linearSE_stamp(this,Yi_lin,Yj_lin,Yk_lin,CKT_Node_Map)              
 
            RTUbus = this.busNum;
            OneVect = ones(length(RTUbus),1);
            Gain = 1./(2*this.RTU_weight);
            
            RealRTUnode = CKT_Node_Map.Bus.NR(RTUbus);
            ImagRTUnode = CKT_Node_Map.Bus.NI(RTUbus);
            adjRealRTUnode = CKT_Node_Map.Bus.LNR(RTUbus);
            adjImagRTUnode = CKT_Node_Map.Bus.LNI(RTUbus);            

            [Yi_lin,Yj_lin,Yk_lin] = StampShuntAdmittance(Yi_lin,Yj_lin,Yk_lin,this.G_mean,-(this.B_mean),RTUbus,CKT_Node_Map);
            
            % Adding the RTU error sources (Feasibility sources):
            Yi_lin = [Yi_lin;RealRTUnode;ImagRTUnode];
            Yj_lin = [Yj_lin;adjRealRTUnode;adjImagRTUnode];
            Yk_lin = [Yk_lin;-Gain.*OneVect;-Gain.*OneVect];       
        end   

        function [Yi_lz,Yj_lz,Yk_lz,Ji_lz,Jk_lz]=stamp_linearized(this,Yi_lz,Yj_lz,Yk_lz,Ji_lz,Jk_lz,Sol,epsilon,CKT_Node_Map)               
            RTUbus = this.busNum;
            Weight = this.RTU_weight;            
            G_Node = CKT_Node_Map.RTU.deltaG;
            B_Node = CKT_Node_Map.RTU.deltaB;
            Gmax_Node = CKT_Node_Map.RTU.MGmax;
            Gmin_Node = CKT_Node_Map.RTU.MGmin;
            Bmax_Node = CKT_Node_Map.RTU.MBmax;
            Bmin_Node = CKT_Node_Map.RTU.MBmin;      
            RealNode = CKT_Node_Map.Bus.NR(RTUbus);
            ImagNode = CKT_Node_Map.Bus.NI(RTUbus);
            adjRealNode = CKT_Node_Map.Bus.LNR(RTUbus);
            adjImagNode = CKT_Node_Map.Bus.LNI(RTUbus);
            
            % k^th Iterate Values:    
            G_k = Sol(G_Node); % 
            B_k = Sol(B_Node); %            
            VR_k = Sol(RealNode);
            VI_k = Sol(ImagNode);
            LR_k = Sol(adjRealNode);
            LI_k = Sol(adjImagNode);      
            OneVect = ones(length(RTUbus),1);

            % stamp linearized Y_RTU
            [Yi_lz,Yj_lz,Yk_lz] = StampShuntAdmittance(Yi_lz,Yj_lz,Yk_lz,G_k,-B_k,RTUbus,CKT_Node_Map);
            % Node indices corresponding to voltages of RTU
            VR_VI_adjVR_adjVI = [RealNode;ImagNode;adjRealNode;adjImagNode]; 
            % Real and Imaginary RTU currents at k^th iteration:
            IRrtu = G_k.*VR_k+B_k.*VI_k;
            IIrtu = G_k.*VI_k-B_k.*VR_k;
            adjIRrtu = G_k.*LR_k-B_k.*LI_k;
            adjIIrtu = G_k.*LI_k+B_k.*LR_k;            
            % Stamp G part:
            Yi_lz = [Yi_lz;VR_VI_adjVR_adjVI];
            Yj_lz = [Yj_lz;repmat(G_Node,[4,1])];
            Yk_lz = [Yk_lz;VR_k;VI_k;LR_k;LI_k];
            %Stamp B part:
            Yi_lz = [Yi_lz;VR_VI_adjVR_adjVI];
            Yj_lz = [Yj_lz;repmat(B_Node,[4,1])];
            Yk_lz = [Yk_lz;VI_k;-VR_k;-LI_k;LR_k];
            
            Ji_lz = [Ji_lz;VR_VI_adjVR_adjVI];
            Jk_lz = [Jk_lz;IRrtu;IIrtu;adjIRrtu;adjIIrtu];
 
            % G equation:
            Yi_lz = [Yi_lz;repmat(G_Node,[7,1])];
            Yj_lz = [Yj_lz;VR_VI_adjVR_adjVI;G_Node;Gmax_Node;Gmin_Node];
            Yk_lz = [Yk_lz;LR_k;LI_k;VR_k;VI_k;2*Weight.*OneVect;OneVect;-OneVect];
            Ji_lz = [Ji_lz;G_Node];
            Jk_lz = [Jk_lz;(2*Weight.*this.G_mean+LR_k.*VR_k+VI_k.*LI_k)];
            % B equation:
            Yi_lz = [Yi_lz;repmat(B_Node,[7,1])];
            Yj_lz = [Yj_lz;VR_VI_adjVR_adjVI;B_Node;Bmax_Node;Bmin_Node];
            Yk_lz = [Yk_lz;-LI_k;LR_k;VI_k;-VR_k;2*Weight.*OneVect;OneVect;-OneVect];
            Ji_lz = [Ji_lz;B_Node];
            Jk_lz = [Jk_lz;(2*Weight.*this.B_mean+VI_k.*LR_k-VR_k.*LI_k)];
 % Stamping Control Ckts (Complementary Slackness conditions):  
            % Bounds on G variable:
            [Yi_lz,Yj_lz,Yk_lz,Ji_lz,Jk_lz] = StampUpLowBounds(Sol,this.G_min,this.G_max,epsilon,Ji_lz,Jk_lz,Yi_lz,Yj_lz,Yk_lz,Gmax_Node,Gmin_Node,G_Node);
            % Bounds on B variable
            [Yi_lz,Yj_lz,Yk_lz,Ji_lz,Jk_lz] = StampUpLowBounds(Sol,this.B_min,this.B_max,epsilon,Ji_lz,Jk_lz,Yi_lz,Yj_lz,Yk_lz,Bmax_Node,Bmin_Node,B_Node);                       
        end  
    end
end