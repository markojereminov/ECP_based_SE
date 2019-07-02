function [EstimatedState] = runSE(casedata,NonlinearSEswitch,PMU_Gerr,RTU_weight,damp_coeff,epsilon,tol,maxIter)
% FUNCTION DESCRIPTION: Main function that runs the SE problem
%___________________________________________________________________________________________________    
% INPUT: 
    % casedata: input augmented mpc file with PMU and RTU data
    % NonlinearSEswitch: if 1: Solving a nonlinear SE from [1]-[2], if 0 solving a linear version presented in [3]-[4]
    % PMU_Gerr: PMU error condutance values
    % RTU_weight
    % damp_coeff: diode damping coefficient
    % epsilon: diode epsilon pramter - complementary slackness approximation (see [1])
    % tol: NR tolerance for convergence
    % maxIter: maximum NR iteration count 
%___________________________________________________________________________________________________    
% OUTPUT:
    % EstimatedState: Structures that give estimated states
%___________________________________________________________________________________________________ 
% References:
% [1] M. Jereminov, A. Jovicic, M. Wagner, G. Hug, L. Pileggi, Equivalent Circuit Programming for Estimating
%     the State of a Power System, in Proc. IEEE PowerTech Milan, June 2019.
% [2] A. Jovicic, M. Jereminov, L. Pileggi, and G. Hug An equivalent circuit formulation for power system state
%     estimation including PMUs, 2018 North American Power Symposium (NAPS).
% [3] A. Jovicic, M. Jereminov, L. Pileggi, and G. Hug, A Linear Formulation for Power System State Estimation 
%     including RTU and PMU Measurements, in Proc. of ISGT Europe 2019, Sept. 2019, Bucharest, Romania.
% [4] Martin R. Wagner, Marko Jereminov, Amritanshu Pandey, Larry Pileggi, "A Probabilistic Approach to Power System 
%    State Estimation using a Linear Algorithm," in Proc. of IEEE/EEEIC Genoa 2019, June 2019, Genoa, Italy
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
%% Setting the default input parameters:
if nargin < 1
   casedata = case14_SE; %default case if not specified in input
end
if nargin < 2
    NonlinearSEswitch = 0; % Solve linear SE
end
if nargin < 3
    PMU_Gerr = 10*ones(length(casedata.Vpmu(:,1)),1); % Defining PMU error conductance:
end
if nargin < 4
    RTU_weight = 0.1*ones(length(casedata.RTU(:,1)),1); % RTU objective function weight
end
if nargin < 5
    damp_coeff = 0.9; % diode limiting coeff
end
if nargin < 6
    epsilon = 1e-6; % diode epsilon pramter
end
if nargin < 7
    tol= 1e-4; % Setting the tolerance for convergence
end
if nargin < 8
    maxIter = 50; 
end
%% Parsing the data 
tic;[RTUClass,PMUClass,TxLineClass,XformsClass,ShuntClass,PShiftClass,n_elem,CKT_Node_Map,varCount] = mpc_parser(casedata,PMU_Gerr,RTU_weight); 
fprintf('Parsing Completed. ParsingTime = %d seconds\n', toc);
%% Stamp Linear network Elements: 
tic;[Yi_lin,Yj_lin,Yk_lin] = StampLinearNetwork(CKT_Node_Map,TxLineClass,XformsClass,ShuntClass,PShiftClass);
fprintf('Linear Network Elements Stamped. Time = %d seconds\n\n', toc);
%% Solve the Linear SE [xx]:
tic;[LinearStateEstimate] = solve_linearSE(Yi_lin,Yj_lin,Yk_lin,CKT_Node_Map,RTUClass,PMUClass,n_elem);
%% if NonlinearSEswitch is ON use Linear SE to initialize the nonlinear one, else use the linear SE results:
if NonlinearSEswitch == 1
    disp('Solving nonlinear SE')
    % Initializing the solution vector:
    [Sol] = InitializeNonlinearSE(varCount,LinearStateEstimate,RTUClass,PMUClass,CKT_Node_Map,n_elem,epsilon);
    % Linear sparse contribution to the Hessian of necessary optimality conditions:
    [Y_linear,J_linear] = build_linear(Yi_lin,Yj_lin,Yk_lin,PMUClass,varCount,CKT_Node_Map);
    % Solve the linearized ECP representation of nonlinear SE (Newton Raphson iterations):
    [Sol,sucess] = runNR(Sol,varCount,Y_linear,J_linear,CKT_Node_Map,RTUClass,PMUClass,epsilon,tol,maxIter,damp_coeff);
else
    disp('Solving Linear SE')
    Sol = LinearStateEstimate;
    sucess = 1;
end
% if solver converged display results:
if sucess == 1
    fprintf('\nEstimated State Found. Simulation Time = %d seconds\n\n', toc);
    
    Bus = casedata.bus(:,1);
    RealEstVoltage = Sol(CKT_Node_Map.Bus.NR);
    ImagEstVoltage = Sol(CKT_Node_Map.Bus.NI);
    VoltMag = abs(RealEstVoltage + 1i*ImagEstVoltage);
    VoltAng = angle(RealEstVoltage + 1i*ImagEstVoltage)*180/pi;
    
    EstimatedState.BusNum = Bus;
    EstimatedState.RectVolt.VR = RealEstVoltage;
    EstimatedState.RectVolt.VI = ImagEstVoltage;
    EstimatedState.PolarVolt.Vmag = VoltMag;
    EstimatedState.PolarVolt.Vang = VoltAng;
    RefVoltMag = casedata.bus(:,8);
    RefVoltAng = casedata.bus(:,9);
    
    if n_elem.Bus < 1000
        fprintf('Estimated Voltage States vs "True" State from input file:\n\n');
        disp(table(Bus,VoltMag,RefVoltMag,VoltAng,RefVoltAng))
    end
        
    % Compute L1, L2 ans LINF norms of Estimated States (with respect to the "true" one form the input file)
    VR_ref = RefVoltMag.*cosd(RefVoltAng);
    VI_ref = RefVoltMag.*sind(RefVoltAng);

    L1_norm = norm([VR_ref-RealEstVoltage;VI_ref-ImagEstVoltage],1);
    L2_norm = norm([VR_ref-RealEstVoltage;VI_ref-ImagEstVoltage],2);
    Linf_norm = norm([VR_ref-RealEstVoltage;VI_ref-ImagEstVoltage],inf);
    
    fprintf('L1, L2 and L_inf norms of estimated state vector in p.u. :\n\n');
    disp(table(L1_norm,L2_norm,Linf_norm))
else
    EstimatedState = [];
end
end