classdef PMU_device < handle
%% Class that defines the PMU measurement device
%% This class will contain the properties required for the respective equivalent circuit stamps as well as the following functions:
% PMU_device(.)=>the function that initializes the properties of the class
%       INPUT: bus_data: PMU device bus number
%              data_Gerr: PMU error conductance (see [1]-[2])
%              data_IR: bounds on IR measurement of PMU
%              data_II: bounds on II measurement of PMU
%              data_VR: bounds on VR measurement of PMU
%              data_VI: bounds on VI measurement of PMU
% linearSE_stamp(.)=> the function that stamps an equvalent circuit of a respective device for the linear SE formulation (unbounded PMU [2]-[3])
%        INPUT: Yi_lin: vector of prior i indices of the NR Sensitivity matrix
%               Yj_lin: vector of prior j indices of the NR Sensitivity matrix
%               Yk_lin: vector of prior k values of the NR Sensitivity matrix
%               CKT_Node_Map: a node map of the ECP circuit
%        OUTPUT:Yi_lin: vector of concatinated i indices of the NR Sensitivity matrix 
%               Yj_lin: vector of concatinated j indices of the NR Sensitivity matrix
%               Yk_lin: vector of concatinated k values of the NR Sensitivity matrix
%               Ji_lin: vector of concatinated i indices of the NR Constant vector
%               Jk_lin: vector of concatinated k values of the NR Constant vector
% stamp_linear(.)=> the function that stamps an equvalent circuit of a respective device
%         INPUT:Yi_vect: vector of prior i indices of the NR Sensitivity matrix
%               Yj_vect: vector of prior j indices of the NR Sensitivity matrix
%               Yk_vect: vector of prior k values of the NR Sensitivity matrix
%               CKT_Node_Map: a node map of the ECP circuit
%         OUTPUT:Yi_vect: vector of concatinated i indices of the NR Sensitivity matrix 
%                Yj_vect: vector of concatinated j indices of the NR Sensitivity matrix
%                Yk_vect: vector of concatinated k values of the NR Sensitivity matrix
%                Ji_lin: vector of concatinated i indices of the NR Constant vector
%                Jk_lin: vector of concatinated k values of the NR Constant vector
% stamp_linearized(.) =>
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
%% LICENSE:
%   This file is part of open source version of ECP based Static State Estimator.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%% References:
% [1] M. Jereminov, A. Jovicic, M. Wagner, G. Hug, L. Pileggi, Equivalent Circuit Programming for Estimating
%     the State of a Power System, in Proc. IEEE PowerTech Milan, June 2019.
% [2] A. Jovicic, M. Jereminov, L. Pileggi, and G. Hug, A Linear Formulation for Power System State Estimation 
%     including RTU and PMU Measurements, in Proc. of ISGT Europe 2019, Sept. 2019, Bucharest, Romania.
% [3] Martin R. Wagner, Marko Jereminov, Amritanshu Pandey, Larry Pileggi, "A Probabilistic Approach to Power System 
%    State Estimation using a Linear Algorithm," in Proc. of IEEE/EEEIC Genoa 2019, June 2019, Genoa, Italy
%___________________________________________________________________________________________________ 
%% Properties of the PMU measurement class:
    properties
        busNum
        G_PMU % PMU error conductance
        IR_max % real voltage PMU measurement
        IR_min % real voltage PMU measurement
        II_max % imag voltage PMU measurement  
        II_min % imag voltage PMU measurement
        VR_max % real voltage PMU measurement
        VR_min % real voltage PMU measurement
        VI_max % imag voltage PMU measurement  
        VI_min % imag voltage PMU measurement
    end
%% Methods:
    methods
        
        function this = PMU_device(bus_data,data_Gerr,data_IR,data_II,data_VR,data_VI) %initialize variables
            this.busNum = bus_data;
            this.G_PMU=data_Gerr; %  PMU error conductance
            this.IR_max=data_IR(:,1); % real PMU voltage
            this.IR_min=data_IR(:,2); % real PMU voltage
            this.II_max=data_II(:,1); % imag PMU voltage
            this.II_min=data_II(:,2); % imag PMU voltage
            this.VR_max=data_VR(:,1); % imag PMU voltage
            this.VR_min=data_VR(:,2); % imag PMU voltage
            this.VI_max=data_VI(:,1); % imag PMU voltage
            this.VI_min=data_VI(:,2); % imag PMU voltage
        end
        
        function [Yi_lin,Yj_lin,Yk_lin,Ji_lin,Jk_lin]=linearSE_stamp(this,Yi_lin,Yj_lin,Yk_lin,CKT_Node_Map)

            PMU_BusID = this.busNum;
            G_err = this.G_PMU; % error conductance
            VRpmu = 0.5*(this.VR_max + this.VR_min); % real mean PMU voltage
            VIpmu = 0.5*(this.VI_max + this.VI_min); % imaginary mean PMU voltage
            IRpmu = 0.5*(this.IR_max + this.IR_min); % real mean PMU voltage
            IIpmu = 0.5*(this.II_max + this.II_min); % imaginary mean PMU voltage               

            % Finding the corresponding nodes in the circuit:
            Node_Real = CKT_Node_Map.Bus.NR(PMU_BusID); % real powerflow node index related to the PMU bus
            Node_Imag = CKT_Node_Map.Bus.NI(PMU_BusID); % imag powerflow node index related to the PMU bus
            adjNode_Real = CKT_Node_Map.Bus.LNR(PMU_BusID); % real adjoint node index related to the PMU bus
            adjNode_Imag = CKT_Node_Map.Bus.LNI(PMU_BusID); % imag adjoint node index related to the PMU bus  
            
            % KCL at the node related to real PMU voltage: 
            [Yi_lin,Yj_lin,Yk_lin] = StampShuntAdmittance(Yi_lin,Yj_lin,Yk_lin,G_err,0,PMU_BusID,CKT_Node_Map);
            Ji_lin = [Node_Real;Node_Imag];
            Jk_lin = [(IRpmu+G_err.*VRpmu);(IIpmu+G_err.*VIpmu)]; 
            
            % Adding the gradient of objective function (source to the adjoint ckt):
            Yi_lin = [Yi_lin;adjNode_Real;adjNode_Imag];
            Yj_lin = [Yj_lin;Node_Real;Node_Imag];
            Yk_lin = [Yk_lin;2*G_err.^2;2*G_err.^2];
            
            Ji_lin = [Ji_lin;adjNode_Real;adjNode_Imag];
            Jk_lin = [Jk_lin;2*VRpmu.*G_err.^2;2*VIpmu.*G_err.^2];                
        end

        function [Yi_lin,Yj_lin,Yk_lin,Ji_lin,Jk_lin]=stamp_linear(this,Yi_lin,Yj_lin,Yk_lin,CKT_Node_Map)

            PMU_BusID = this.busNum;
            OneVect = ones(length(PMU_BusID),1);
            G_err = this.G_PMU; % error conductance
            IR_MAX = this.IR_max;
            IR_MIN = this.IR_min;
            II_MAX = this.II_max;
            II_MIN = this.II_min;            
            VR_MAX = this.VR_max;
            VR_MIN = this.VR_min;
            VI_MAX = this.VI_max;
            VI_MIN = this.VI_min;
            
            % Finding the corresponding nodes in the ECP circuit:
            Node_Real = CKT_Node_Map.Bus.NR(PMU_BusID); % real powerflow node index related to the PMU bus
            Node_Imag = CKT_Node_Map.Bus.NI(PMU_BusID); % imag powerflow node index related to the PMU bus
            adjNode_Real = CKT_Node_Map.Bus.LNR(PMU_BusID); % real adjoint node index related to the PMU bus
            adjNode_Imag = CKT_Node_Map.Bus.LNI(PMU_BusID); % imag adjoint node index related to the PMU bus  
            VR_PMU_Node = CKT_Node_Map.PMU.VR;
            VI_PMU_Node = CKT_Node_Map.PMU.VI;
            IRs_Node = CKT_Node_Map.PMU.IsR;
            IIs_Node = CKT_Node_Map.PMU.IsI;
            IR_PMU_Node = CKT_Node_Map.PMU.IR;
            II_PMU_Node = CKT_Node_Map.PMU.II;             
            adjIRs_Node = CKT_Node_Map.PMU.adjIsR;
            adjIIs_Node = CKT_Node_Map.PMU.adjIsI; 
            adjVR_PMU_Node = CKT_Node_Map.PMU.LR;
            adjVI_PMU_Node = CKT_Node_Map.PMU.LI;            
            IRmx_Node = CKT_Node_Map.PMU.MIRmax;
            IRmn_Node = CKT_Node_Map.PMU.MIRmin;
            IImx_Node = CKT_Node_Map.PMU.MIImax;
            IImn_Node = CKT_Node_Map.PMU.MIImin;
            VRmx_Node = CKT_Node_Map.PMU.MVRmax;
            VRmn_Node = CKT_Node_Map.PMU.MVRmin;
            VImx_Node = CKT_Node_Map.PMU.MVImax;
            VImn_Node = CKT_Node_Map.PMU.MVImin;            
            % KCL at the node related to real PMU voltage source: 
            % in Modified Nodal Analysis (MNA) for every voltage source we can add a current variable
            % for which we add the voltage constraint, i.e. add I_source in
            % KCL equations and for those variables add voltage equations,
            % e.g. VR_PMU = VR_measured.
            Yi_lin = [Yi_lin;VR_PMU_Node;VR_PMU_Node;VI_PMU_Node;VI_PMU_Node;Node_Real;Node_Imag];
            Yj_lin = [Yj_lin;IRs_Node;IR_PMU_Node;IIs_Node;II_PMU_Node;IR_PMU_Node;II_PMU_Node];
            Yk_lin = [Yk_lin;repmat(OneVect,[4,1]);-repmat(OneVect,[2,1])]; 
            
            % Stamping Error conductance at PMU node to KCL equations of
            % the real and imaginary ckts:
            % @node VR_PMU: G_err*(VR_PMU-VR_term)
            % @node VI_PMU: G_err*(VI_PMU-VI_term)
            Yi_lin = [Yi_lin;VR_PMU_Node;VR_PMU_Node;VI_PMU_Node;VI_PMU_Node];
            Yj_lin = [Yj_lin;VR_PMU_Node;Node_Real;VI_PMU_Node;Node_Imag];
            Yk_lin = [Yk_lin;G_err;-G_err;G_err;-G_err];
            % Stamping Error conductance at terminal node to KCL equations of
            % the real and imaginary ckts:
            % @node VR_PMU: G_err*(-VR_PMU+VR_term)
            % @node VI_PMU: G_err*(-VI_PMU+VI_term)
            Yi_lin = [Yi_lin;Node_Real;Node_Real;Node_Imag;Node_Imag];
            Yj_lin = [Yj_lin;Node_Real;VR_PMU_Node;Node_Imag;VI_PMU_Node];
            Yk_lin = [Yk_lin;G_err;-G_err;G_err;-G_err]; 
          
            %@ adjoint real PMU node: G_err*(adjVRpmu - adjVR)+adjIsR
            %@ adjoint real terminal node: G_err*(adjVR - adjVRpmu)
            Yi_lin = [Yi_lin;repmat(adjVR_PMU_Node,[3,1]);repmat(adjNode_Real,[2,1])];
            Yj_lin = [Yj_lin;adjVR_PMU_Node;adjNode_Real;adjIRs_Node;adjNode_Real;adjVR_PMU_Node];
            Yk_lin = [Yk_lin;G_err;-G_err;OneVect;G_err;-G_err];            
            %@ adjoint imag PMU node: G_err*(adjVIpmu - adjVI)+adjIsI
            %@ adjoint imag terminal node: G_err*(adjVI - adjVIpmu)
            Yi_lin = [Yi_lin;repmat(adjVI_PMU_Node,[3,1]);repmat(adjNode_Imag,[2,1])];
            Yj_lin = [Yj_lin;adjVI_PMU_Node;adjNode_Imag;adjIIs_Node;adjNode_Imag;adjVI_PMU_Node];
            Yk_lin = [Yk_lin;G_err;-G_err;OneVect;G_err;-G_err];             
            
            % Now add the "source" - Gradient of objective function to real adjoint ckt:
            Yi_lin = [Yi_lin;repmat(adjVR_PMU_Node,[2,1]);repmat(adjNode_Real,[2,1])];
            Yj_lin = [Yj_lin;Node_Real;VR_PMU_Node;Node_Real;VR_PMU_Node];
            Yk_lin = [Yk_lin;-2*G_err.^2; 2*G_err.^2;2*G_err.^2; -2*G_err.^2];     
            % Now add the "source" - Gradient of objective function to real adjoint ckt:
            Yi_lin = [Yi_lin;repmat(adjVI_PMU_Node,[2,1]);repmat(adjNode_Imag,[2,1])];
            Yj_lin = [Yj_lin;Node_Imag;VI_PMU_Node;Node_Imag;VI_PMU_Node];
            Yk_lin = [Yk_lin;-2*G_err.^2; 2*G_err.^2;2*G_err.^2; -2*G_err.^2];             
           
            % Set the adjoint PMU voltages to zero:
            Yi_lin = [Yi_lin;IRs_Node;IIs_Node];
            Yj_lin = [Yj_lin;adjVR_PMU_Node;adjVI_PMU_Node];
            Yk_lin = [Yk_lin;repmat(OneVect,[2,1])];
            
            % for same bounds (exact PMU measurement) 
            SameVR = OneVect;
            SameVI = OneVect;
            SameIR = OneVect;
            SameII = OneVect;       
            SameVRadj = 0*OneVect;
            SameVIadj = 0*OneVect;
            SameIRadj = 0*OneVect;
            SameIIadj = 0*OneVect; 
            VRset = 0*OneVect;
            VIset = 0*OneVect;
            IRset = 0*OneVect;
            IIset = 0*OneVect;
            
            VRset((VR_MAX-VR_MIN)==0) = VR_MAX((VR_MAX-VR_MIN)==0);
            VIset((VI_MAX-VI_MIN)==0) = VI_MAX((VI_MAX-VI_MIN)==0);
            IRset((IR_MAX-IR_MIN)==0) = IR_MAX((IR_MAX-IR_MIN)==0);
            IIset((II_MAX-II_MIN)==0) = II_MAX((II_MAX-II_MIN)==0);             
            
            SameVR((VR_MAX-VR_MIN)~=0) = 0;
            SameVI((VI_MAX-VI_MIN)~=0) = 0;
            SameVRadj((VR_MAX-VR_MIN)~=0) = 1;
            SameVIadj((VI_MAX-VI_MIN)~=0) = 1;
           
            % adjIsR and adjIsI equations:
            Yi_lin = [Yi_lin;repmat(adjIRs_Node,[4,1]);repmat(adjIIs_Node,[4,1])];            
            Yj_lin = [Yj_lin;VRmx_Node;VRmn_Node;VR_PMU_Node;adjIRs_Node;VImx_Node;VImn_Node;VI_PMU_Node;adjIIs_Node];
            Yk_lin = [Yk_lin;-OneVect;OneVect;SameVR;SameVRadj;-OneVect;OneVect;SameVI;SameVIadj]; 
            
            Ji_lin = [adjIRs_Node;adjIIs_Node];
            Jk_lin = [VRset;VIset];

            % IRpmu and IIpmu equations:   
            % LR_PMU - LR + MU_IR_max - MU_IR_min = 0 if IR is bounded or IR = IR_set if IR is fixed
            % LI_PMU - LI + MU_II_max - MU_II_min = 0 if II is bounded or II = II_set if II is fixed
            SameIR((IR_MAX-IR_MIN)~=0) = 0;
            SameII((II_MAX-II_MIN)~=0) = 0;             
            SameIRadj((IR_MAX-IR_MIN)~=0) = 1;
            SameIIadj((II_MAX-II_MIN)~=0) = 1; 
            
            Yi_lin = [Yi_lin;repmat(IR_PMU_Node,[5,1]);repmat(II_PMU_Node,[5,1])];
            Yj_lin = [Yj_lin;IRmx_Node;IRmn_Node;adjVR_PMU_Node;IR_PMU_Node;adjNode_Real;IImx_Node;IImn_Node;adjVI_PMU_Node;II_PMU_Node;adjNode_Imag];
            Yk_lin = [Yk_lin;OneVect;-OneVect;OneVect;SameIR;-SameIRadj;OneVect;-OneVect;OneVect;SameII;-SameIIadj;];
            
            Ji_lin = [Ji_lin;IR_PMU_Node;II_PMU_Node];
            Jk_lin = [Jk_lin;IRset;IIset];            
        end
        
        function [Yi_vect,Yj_vect,Yk_vect,Ji_vect,Jk_vect]=stamp_linearized(this,Yi_vect,Yj_vect,Yk_vect,Ji_vect,Jk_vect,Sol,epsilon,CKT_Node_Map)    
            IR_Node = CKT_Node_Map.PMU.IR;
            II_Node = CKT_Node_Map.PMU.II;
            VR_Node = CKT_Node_Map.PMU.VR;
            VI_Node = CKT_Node_Map.PMU.VI;            
            IRmx_Node = CKT_Node_Map.PMU.MIRmax;
            IRmn_Node = CKT_Node_Map.PMU.MIRmin;
            IImx_Node = CKT_Node_Map.PMU.MIImax;
            IImn_Node = CKT_Node_Map.PMU.MIImin;
            VRmx_Node = CKT_Node_Map.PMU.MVRmax;
            VRmn_Node = CKT_Node_Map.PMU.MVRmin;
            VImx_Node = CKT_Node_Map.PMU.MVImax;
            VImn_Node = CKT_Node_Map.PMU.MVImin;                           
            %% IR control ckts:
            [Yi_vect,Yj_vect,Yk_vect,Ji_vect,Jk_vect] = StampUpLowBounds(Sol,this.IR_min,this.IR_max,epsilon,Ji_vect,Jk_vect,Yi_vect,Yj_vect,Yk_vect,IRmx_Node,IRmn_Node,IR_Node);
            %% II control ckts:
            [Yi_vect,Yj_vect,Yk_vect,Ji_vect,Jk_vect] = StampUpLowBounds(Sol,this.II_min,this.II_max,epsilon,Ji_vect,Jk_vect,Yi_vect,Yj_vect,Yk_vect,IImx_Node,IImn_Node,II_Node);        
            %% VR control ckts:
            [Yi_vect,Yj_vect,Yk_vect,Ji_vect,Jk_vect] = StampUpLowBounds(Sol,this.VR_min,this.VR_max,epsilon,Ji_vect,Jk_vect,Yi_vect,Yj_vect,Yk_vect,VRmx_Node,VRmn_Node,VR_Node);
            %% VI control ckts
            [Yi_vect,Yj_vect,Yk_vect,Ji_vect,Jk_vect] = StampUpLowBounds(Sol,this.VI_min,this.VI_max,epsilon,Ji_vect,Jk_vect,Yi_vect,Yj_vect,Yk_vect,VImx_Node,VImn_Node,VI_Node);               
        end  
    end
end