function [RTUClass,PMUClass,TxLineClass,XformsClass,ShuntClass,PShiftClass,n_elem,CKT_Node_Map,varCount] = mpc_parser(mpc,data_Gerr,RTU_weight)
% FUNCTION DESCRIPTION: Function that parses in the mpc file and builds the 
% structures for each of the elements of the ECP representation of a SE problem:
%___________________________________________________________________________________________________    
% INPUT: 
    % mpc: augmented mpc input file
    % data_Gerr: vector of PMU error conductances
    % RTU_weight: vector of RTU weights
%___________________________________________________________________________________________________    
% OUTPUT:
    % RTUClass: defined class of RTU devices
    % PMUClass: defined class of PMU devices
    % TxLineClass: defined class of Tx Lines
    % XformsClass: defined class of in phase xfrms
    % ShuntClass: defined class of fixed shunt devices
    % PShiftClass: defined class of phase shifting xfrms
    % n_elem: structure defining number of respective elements 
    % CKT_Node_Map: node map (bus->ckt nodes)
    % varCount: total variable count
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
%% convert to internal indexing
OriginalBus = mpc.bus(:,1);
BusNum = length(OriginalBus);
Correction_Bus_Map(OriginalBus,1) = (1:BusNum);
% Converting Bus numbers to internal indexing:
mpc.bus(:,1) = Correction_Bus_Map(OriginalBus);
% Converting RTU devices to internal indexing:
mpc.RTU(:,1) = Correction_Bus_Map(mpc.RTU(:,1));
% Converting PMU devices to internal indexing:
mpc.Vpmu(:,1) = Correction_Bus_Map(mpc.Vpmu(:,1));
mpc.Ipmu(:,1) = Correction_Bus_Map(mpc.Ipmu(:,1));
% Converting Network data to internal indexing:
mpc.branch(:,1) = Correction_Bus_Map(mpc.branch(:,1));
mpc.branch(:,2) =Correction_Bus_Map(mpc.branch(:,2));

[baseMVA,bus,branch,RTU,Vpmu,Ipmu] = deal(mpc.baseMVA, mpc.bus, mpc.branch,mpc.RTU,mpc.Vpmu,mpc.Ipmu);
%% Processing fixed shunt data: 
YshuntID = all([bus(:,5),bus(:,6)]  ==0,2); % finding zero shunt indices
Ysh_bus_N = bus(:,1);
Ysh_bus_N(YshuntID,:)=[]; % deleting zero shunts
% Filling the Shunt Class
if isempty(Ysh_bus_N)
    ShuntClass = [];
else
    ShuntClass = FixedShunt(Ysh_bus_N,bus(Ysh_bus_N,5)/baseMVA,bus(Ysh_bus_N,6)/baseMVA);   
end
%% Processing TxLine data: 
from_bus = branch(:,1);
to_bus = branch(:,2);
T_line_tap_zero_ID = all([branch(:,9) branch(:,10)] ==0,2); %transmission line is t = 0;
branch(T_line_tap_zero_ID,9)=1;% Set zero turns ratio to 1

T_line_tap_ID = all(branch(:,9) ==1 & branch(:,10) ==0,2); %transmission line is t = 1;
TL_from_bus = from_bus(T_line_tap_ID); %from bus number
TL_to_bus = to_bus(T_line_tap_ID); %to bus number

branch_num = size(branch);
Series_index = (1:1:branch_num(1))'; 
TL_index = Series_index(T_line_tap_ID); %index of transmission line elements

%Transmission Line Pi_model admittances:
TxLineData(:,1) = TL_from_bus; %From bus
TxLineData(:,2) = TL_to_bus; %To bus
TxLineData(:,3) = branch(TL_index,3)./(branch(TL_index,3).^2+branch(TL_index,4).^2); %Gs
TxLineData(:,4) = -branch(TL_index,4)./(branch(TL_index,3).^2+branch(TL_index,4).^2); %Bs
TxLineData(:,5) = branch(TL_index,5); %Bsh  
% Filling the Transmission Line class
TxLineClass= TxLine(TxLineData);
%% Processing in-Phase xfmr data: 
Trafo_ID = setdiff(Series_index,TL_index); %index of all trafo elements
Trafo_index = Trafo_ID((branch(Trafo_ID,10)) == 0);

%Calculate PI-model parameters and store them in the new field:
trafo_tap = branch(Trafo_index,9);
Y_km = 1./(branch(Trafo_index,3)+1i*branch(Trafo_index,4));
Gkm = real(Y_km); %G_trafo
Bkm = imag(Y_km); %B_trafo
Bsh_trafo = branch(Trafo_index,5)/2;

XfmrData(:,1) = from_bus(Trafo_index); % from bus
XfmrData(:,2) = to_bus(Trafo_index); % to bus
XfmrData(:,3) = Gkm./trafo_tap; %G_series
XfmrData(:,4) = Bkm./trafo_tap; %B_series
XfmrData(:,5) = ((1-trafo_tap)./(trafo_tap.^2)).*Gkm; %G_from
XfmrData(:,6) = ((1-trafo_tap)./(trafo_tap.^2)).*Bkm + Bsh_trafo./(trafo_tap.^2); %B_from
XfmrData(:,7) = Gkm.*(1-1./trafo_tap); %G_to
XfmrData(:,8) = Bkm.*(1-1./trafo_tap) + Bsh_trafo; %B_to

if isempty(Trafo_index)
    XformsClass = [];
else
    XformsClass = InPhaseXfmr(XfmrData);  
end
%% Processing phase shifting xfmr data: 
PShift_index = Trafo_ID((branch(Trafo_ID,10)) ~= 0);
PShift_tap = branch(PShift_index,9);
PShift_tap(PShift_tap ==0) = 1;
 
Shift_angle = branch(PShift_index,10);
Y_km_shift = 1./(branch(PShift_index,3)+1i*branch(PShift_index,4));
alpha_shift = real(Y_km_shift); %real admittance
beta_shift = imag(Y_km_shift)+branch(PShift_index,5)/2; %imaginary admittance

PShiftData(:,1) = from_bus(PShift_index); % from bus
PShiftData(:,2) = to_bus(PShift_index); % to bus
PShiftData(:,3) = alpha_shift./PShift_tap.^2; %
PShiftData(:,4) = beta_shift./PShift_tap.^2; %
PShiftData(:,5) = -(cosd(Shift_angle).*real(Y_km_shift) - sind(Shift_angle).*imag(Y_km_shift))./PShift_tap; %
PShiftData(:,6) = -(sind(Shift_angle).*real(Y_km_shift) + cosd(Shift_angle).*imag(Y_km_shift))./PShift_tap;%
PShiftData(:,7) = -(cosd(Shift_angle).*real(Y_km_shift) + sind(Shift_angle).*imag(Y_km_shift))./PShift_tap;%
PShiftData(:,8) = -(cosd(Shift_angle).*imag(Y_km_shift) - sind(Shift_angle).*real(Y_km_shift))./PShift_tap;%
PShiftData(:,9) = alpha_shift;%
PShiftData(:,10) = beta_shift;%

if isempty(PShift_index)
    PShiftClass = [];
else
    PShiftClass = PhaseShiftXfmr(PShiftData);  
end
%% Processing RTU data:
if ~isempty(mpc.RTU)
    RTU_data.BusNum = RTU(:,1); % RTU bus number
    
    % Generate confidence intervals on RTU GB variables:
    [Gmax, Gmin, Bmax, Bmin] = GBcalc(RTU(:,2)',RTU(:,3)',RTU(:,4)',RTU(:,5)',RTU(:,6)',RTU(:,7)');
    
    RTU_data.G_Bounds = [Gmax,Gmin]; % bounds on G variable
    RTU_data.G_mean = 0.5*(Gmax+Gmin);
    RTU_data.B_Bounds = [Bmax,Bmin];
    RTU_data.B_mean = 0.5*(Bmax+Bmin);
    
    if any((RTU_data.G_Bounds(:,1) - RTU_data.G_Bounds(:,2)) == 0)
        error('Current Version of the code doesnt support the exact RTU measurements')
    end
    if any((RTU_data.B_Bounds(:,1) - RTU_data.B_Bounds(:,2)) == 0)
        error('Current Version of the code doesnt support the exact RTU measurements')
    end
    if any((RTU_data.G_Bounds(:,1) - RTU_data.G_Bounds(:,2)) <= 0)
       error('RTU G bounds are not consistent')
    end
    if any((RTU_data.B_Bounds(:,1) - RTU_data.B_Bounds(:,2)) <= 0)
       error('RTU B bounds are not consistent')
    end
    
    RTUClass = RTU_device(RTU_weight,RTU_data.BusNum,RTU_data.G_Bounds, RTU_data.B_Bounds, RTU_data.G_mean, RTU_data.B_mean); 
else
    RTUClass = [];
end
%% Processing PMU data:
PMU_data.BusNum = Vpmu(:,1); %reading the PMU bus number

PMU_data.VR_Bounds = [Vpmu(:,2),Vpmu(:,3)]; %reading the PMU VR bounds
PMU_data.VI_Bounds = [Vpmu(:,4),Vpmu(:,5)]; %reading the PMU VI bounds
        
for pmuCount = 1:length(PMU_data.BusNum)  
    I_Temp = Ipmu(PMU_data.BusNum(pmuCount) == Ipmu(:,1),2:5); % summed current bounds related to i^th bus
    PMU_data.IR_Bounds(pmuCount,1:2) = sum(I_Temp(:,1:2),1);
    PMU_data.II_Bounds(pmuCount,1:2) = sum(I_Temp(:,3:4),1);  
end
PMU_data.IR_Bounds = [PMU_data.IR_Bounds(:,1), PMU_data.IR_Bounds(:,2)];
PMU_data.II_Bounds = [PMU_data.II_Bounds(:,1), PMU_data.II_Bounds(:,2)];

if any((PMU_data.IR_Bounds(:,1) - PMU_data.IR_Bounds(:,2)) < 0)
    error('PMU IR bounds are not consistent')
end  
if any((PMU_data.II_Bounds(:,1) - PMU_data.II_Bounds(:,2)) < 0)
    error('PMU II bounds are not consistent')
end
if any((PMU_data.VR_Bounds(:,1) - PMU_data.VR_Bounds(:,2)) < 0)
    error('PMU VR bounds are not consistent')
end  
if any((PMU_data.VI_Bounds(:,1) - PMU_data.VI_Bounds(:,2)) < 0)
    error('PMU VI bounds are not consistent')
end    
PMUClass = PMU_device(PMU_data.BusNum,data_Gerr,PMU_data.IR_Bounds, PMU_data.II_Bounds,PMU_data.VR_Bounds, PMU_data.VI_Bounds);
%% Generate the structure that contains the info about number of devices within the analyzed network:
n_elem.Shunt = length(Ysh_bus_N); %number of shunt elements
n_elem.Xforms = length(XfmrData(:,1)); %number of transformers
n_elem.PShift = length(PShiftData(:,1)); %number of phase shifters
n_elem.RTU = length(RTU_data.BusNum); %total number of RTUs
n_elem.Bus = length(bus(:,1)); 
n_elem.PMU = length(PMU_data.BusNum); %total number of PMUs
%% Generating Map Classes:
[CKT_Node_Map,varCount] = BuildCktNodeMap(n_elem);
end