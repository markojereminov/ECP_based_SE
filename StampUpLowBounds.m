function [Yi_vect,Yj_vect,Yk_vect,Ji_vect,Jk_vect] = StampUpLowBounds(Sol,X_MIN,X_MAX,epsilon,Ji_vect,Jk_vect,Yi_vect,Yj_vect,Yk_vect,Xmax_Node,Xmin_Node,X_Node)
% FUNCTION DESCRIPTION:
%A function that stamps the upper and lower bounds of a device.
%___________________________________________________________________________________________________    
% INPUT: 
    % Sol: Solution Vector
    % X_MIN: vector of lower bounds on respective circuit variables
    % X_MAX: vector of upper bounds on respective circuit variables
    % epsilon: diode (steepness) coefficient, i.e. complementary slackness approximation
    % Ji_vect: vector of prior i indices of the NR Constant Vector
    % Jk_vect: vector of prior k values of the NR Constant Vector
    % Yi_vect: vector of prior i indices of the NR Sensitivity matrix
    % Yj_vect: vector of prior j indices of the NR Sensitivity matrix
    % Yk_vect: vector of prior k values of the NR Sensitivity matrix
    % Xmax_Node: Node indices corresponding to the adjoint variables related to the upper limit
    % Xmin_Node: Node indices corresponding to the adjoint variables related to the lower limit
    % X_Node: Node indices corresponding to the adjoint variables related to the bounded variable
%___________________________________________________________________________________________________    
% OUTPUT:
    % Yi_vect: vector of concatinated i indices of the NR Sensitivity matrix 
    % Yj_vect: vector of concatinated j indices of the NR Sensitivity matrix
    % Yk_vect: vector of concatinated k values of the NR Sensitivity matrix
    % Ji_vect: vector of concatinated i indices of the NR Constant Vector
    % Jk_vect: vector of concatinated k values of the NR Constant Vector
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
SameBoundID = X_MIN==X_MAX; % check for the fixed variables
X_k = Sol(X_Node); % vector of bounded circuit variables at k^th Newton Raphson iteration
MUmin_k = Sol(Xmin_Node); % vector of adjoint variables that correspond to the lower bounds of circuit variables
MUmax_k = Sol(Xmax_Node); % vector of adjoint variables that correspond to the upper bounds of circuit variables 

% Linearized Upper bounds at k^th NR iteration:
% Upper Bound: (X_k-X_MAX).*MUmax_kp1 + MUmax_k.*X_kp1 = MUmax_k.*X_k - epsilon
% Lower Bound: (X_MIN-X_k).*MUmin_kp1 - MUmin_k.*X_kp1 = -MUmin_k.*X_k - epsilon
%Appending the stamps to the triplet vector definition of Sensitivity matrix and Constant vector:
Yi_vect = [Yi_vect;Xmax_Node(~SameBoundID);Xmax_Node(~SameBoundID);Xmin_Node(~SameBoundID);Xmin_Node(~SameBoundID)];
Yj_vect = [Yj_vect;Xmax_Node(~SameBoundID);X_Node(~SameBoundID);Xmin_Node(~SameBoundID);X_Node(~SameBoundID)];
Yk_vect = [Yk_vect;(X_k(~SameBoundID)-X_MAX(~SameBoundID));MUmax_k(~SameBoundID);(X_MIN(~SameBoundID)-X_k(~SameBoundID));-MUmin_k(~SameBoundID)];

Ji_vect = [Ji_vect;Xmax_Node(~SameBoundID);Xmin_Node(~SameBoundID)];
Jk_vect = [Jk_vect;(MUmax_k(~SameBoundID).*X_k(~SameBoundID) - epsilon);(-MUmin_k(~SameBoundID).*X_k(~SameBoundID)-epsilon)];     

%Appending the stamps to the triplet vector definition of Sensitivity matrix and Constant vector for the same bounds (fixed variable):
% MU_max = 0
% MU_min = 0
% equation of X = X_set is added in the other part of stamping function
Yi_vect = [Yi_vect;Xmax_Node(SameBoundID);Xmin_Node(SameBoundID)];
Yj_vect = [Yj_vect;Xmax_Node(SameBoundID);Xmin_Node(SameBoundID)];
Yk_vect = [Yk_vect;ones(2*sum(SameBoundID),1)];
end