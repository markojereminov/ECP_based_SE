function [Sol, sucess] = runNR(Sol,varCount,Y_linear,J_linear,CKT_Node_Map,RTUClass,PMUClass,epsilon,tol,maxIter,damp_coeff)
% FUNCTION DESCRIPTION:
%A function that iteratively solves the linearized circuits, which corresponds to iteratively solving a NR problem (see [1])
%___________________________________________________________________________________________________    
% INPUT: 
    % Sol: initial solution vector
    % varCount: total variable count of the problem
    % Y_linear: linear contribution to Y matrix
    % J_linear: linear contribution to J vector
    % CKT_Node_Map: circuit node map
    % RTUClass
    % PMUClass
    % epsilon: diode coefficient that represents CS approximation
    % tol: NR tolerance (absolute tolerance on NR step)
    % maxIter: maximum iteration count
    % damp_coeff:  diode damping coeff.
%___________________________________________________________________________________________________    
% OUTPUT:
    % Sol: Solution vector 
    % sucess: indicator for convergence of NR
%___________________________________________________________________________________________________ 
% References:
% [1] L. Pileggi, R. Rohrer, and C. Visweswariah, Electronic Circuit & System Simulation Methods. New York, NY, USA: McGraw-Hill, 1995
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
% Initial NR parameters:
ERROR = inf; %initial error
iterCount = 0; %initialize iteration counter

% NR iterations:
while (ERROR > tol) && (iterCount <=maxIter)
    %store the solution from previous iteration
    Sol_old = Sol; 
    % Re-Stamp the linearized terms:
    [i_NR,j_NR,k_NR,J_NR_i,J_NR_k] = StampLinearized(CKT_Node_Map,RTUClass,PMUClass,epsilon,Sol);
    % Define sparse Y and J:
    Y_matrix = Y_linear+sparse([i_NR;varCount],[j_NR;varCount],[k_NR;0]);
    J_vector = J_linear+sparse([J_NR_i;varCount],1,[J_NR_k;0]);
    %Solving for the next iteration:
    Sol = full(Y_matrix\J_vector);
    [Sol] = circuit_limiting(Sol,Sol_old,CKT_Node_Map,RTUClass,PMUClass,epsilon,damp_coeff);
    % Calculate Error and increase iter counter:                                    
    Mismatch = abs(Sol-Sol_old); 
    ERROR = max(Mismatch); 
    iterCount = iterCount+1; %increasing the counter
    % display max mismatch:
    fprintf('iteration #%d, MaxRealVoltageError = %E, MaxImagVoltageError = %E\n',iterCount, max(Mismatch(CKT_Node_Map.Bus.NR)),max(Mismatch(CKT_Node_Map.Bus.NI)));
end
% Setting the sucess indicator:
if (iterCount >= maxIter) || (ERROR > tol)
    sucess = 0;
    fprintf('\n NR didnt converge. Try reducing diode_coeff\n')
else
    sucess = 1;
    fprintf('\n NR converged\n')
end
end