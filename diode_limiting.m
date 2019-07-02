function [deltaX,deltaXmax,deltaXmin] = diode_limiting(Sol,Sol_old,Xnode,XmaxNode,XminNode,X_MAX,X_MIN,epsilon,damp_coeff)
% FUNCTION DESCRIPTION:
%A function that applies basic diode limiting technique (see [1]) to the
%bounded variable
%___________________________________________________________________________________________________    
% INPUT: 
    % Sol: new solution vector
    % Sol_old: solution vector of the previous iteration
    % Xnode: Node indices corresponding to the bounded primal circuit variable
    % XmaxNode: Node indices corresponding to the bounded max_adjoint circuit variable 
    % XminNode: Node indices corresponding to the bounded min_adjoint circuit variable
    % X_MAX: upper bounds on the primal ckt variable
    % X_MIN: lower bounds on the primal ckt variable
    % epsilon: diode coefficient (used to approximate CS conditions)
    % damp_coeff: diode damping coeff that should be (0,1)
%___________________________________________________________________________________________________    
% OUTPUT:
    % deltaX: limited NR step of respective primal ckt variables:
    % deltaXmax: limited NR step of respective adjoint ckt variables:
    % deltaXmin: limited NR step of respective adjoint ckt variables:
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
% Initialize damping vectors:
X_damp = ones(length(Xnode),1);
Xmax_damp = ones(length(Xnode),1);
Xmin_damp = ones(length(Xnode),1);

% Obtaining the lower bound on MU (note that due to the CS approximation MU
% is not > 0, but some small number defined by epsilon and upper and lower
% bounds)
MU_bound = epsilon./(X_MAX-X_MIN);
MU_bound((X_MAX-X_MIN)==0) = 0;

% Finding X_kp1 and X_k of NR iterations
X_kp1 = Sol(Xnode);
X_k = Sol_old(Xnode);
Xmax_kp1 = Sol(XmaxNode);
Xmax_k = Sol_old(XmaxNode);
Xmin_kp1 = Sol(XminNode);
Xmin_k = Sol_old(XminNode);

% not limitited NR steps
deltaX = X_kp1-X_k;
deltaXmax = Xmax_kp1-Xmax_k;
deltaXmin = Xmin_kp1-Xmin_k;

% Ensuring the step feasibility and damping it if approaches towards the bound:
Xmin_Ratio = min(1,damp_coeff*(X_MIN-X_k)./deltaX);
Xmax_Ratio = min(1,damp_coeff*(X_MAX-X_k)./deltaX);
MXmin_Ratio = min(1,damp_coeff*(MU_bound-Xmin_k)./deltaXmin);
MXmax_Ratio = min(1,damp_coeff*(MU_bound-Xmax_k)./deltaXmax);
% Removing the parameters if upp and low bounds are the same:
Xmin_Ratio((X_MAX-X_MIN)==0) = 1;
Xmax_Ratio((X_MAX-X_MIN)==0) = 1;
MXmin_Ratio((X_MAX-X_MIN)==0) = 1;
MXmax_Ratio((X_MAX-X_MIN)==0) = 1;

% Finding the respective damping parameters:
X_damp(deltaX<0) = Xmin_Ratio(deltaX<0);
X_damp(deltaX>0) = Xmax_Ratio(deltaX>0);
Xmax_damp(deltaXmax<0) = MXmax_Ratio(deltaXmax<0); % NR_step approaches MU_boundary only if < 0
Xmin_damp(deltaXmin<0) = MXmin_Ratio(deltaXmin<0); % NR_step approaches MU_boundary only if < 0

% Limiting respective NR steps:
deltaX = deltaX.*X_damp;
deltaXmax = deltaXmax.*Xmax_damp;
deltaXmin = deltaXmin.*Xmin_damp;
end