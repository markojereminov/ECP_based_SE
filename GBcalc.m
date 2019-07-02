function [ Gmax, Gmin, Bmax, Bmin] = GBcalc(Pmin,Pmax,Qmin,Qmax,Vmin,Vmax)
% FUNCTION DESCRIPTION: Function that computes upper and lower bound on GB RTU variables, 
% for more info see [1]:
% INPUT: 
    % Pmin: minimum P bound
    % Pmax: maximum P bound
    % Qmin: minimum Q bound
    % Qmax: maximum Q bound
    % Vmin: minimum voltage bound
    % Vmax: maximum voltage bound

%___________________________________________________________________________________________________    
% OUTPUT:
    % Gmax: upper bound on RTU G var  
    % Gmin: lower bound on RTU G var
    % Bmax: upper bound on RTU B var
    % Bmin: lower bound on RTU B var
% Reference:
% [1] A. Jovicic, M. Jereminov, L. Pileggi, and G. Hug An equivalent circuit formulation for power system state
%     estimation including PMUs, 2018 North American Power Symposium (NAPS).
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
Tmax=max([Vmax.*Vmax; Vmin.*Vmax; Vmin.*Vmin]);
Tmin=min([Vmax.*Vmax; Vmin.*Vmax; Vmin.*Vmin]);

% negative P means generation of real power, while negative Q means
% generated Q (capacitor). 
Gmin=-max([Pmax./Tmax; Pmax./Tmin; Pmin./Tmax; Pmin./Tmin])';
Gmax=-min([Pmax./Tmax; Pmax./Tmin; Pmin./Tmax; Pmin./Tmin])';
Bmin=-max([Qmax./Tmax; Qmax./Tmin; Qmin./Tmax; Qmin./Tmin])';
Bmax=-min([Qmax./Tmax; Qmax./Tmin; Qmin./Tmax; Qmin./Tmin])';
end

