function [sstTemp,ns,altitude,panels,rho,v,radiusOfEarth,MeanMotion,mu,satelliteMass,panelSurface,...
  sstDesiredFunction,windOn,sunOn,deltaAngle,timetemp,totalTime,wakeAerodynamics,masterSatellite]...
  =IRSRendezvousInitial()
% initial conditions for IRS's Discoverers mission:
% r0=6778.137         %% km
% inclinition=10;     %% degree
% RAAN=45;            %% degree
% argumentOfPerigee=130;  %% degree
% truAnomaly=45;      %% degree
% mass=10             %% kg
% panelsurface=2*1.1  %% m
%% position of deputy
% z=82.5;             %% is x in publication
% x=-930.46;          %% is y in publication
% y=55.27;            %% is z in publication
% w=-0.17;            %% is u in publication
% u=-0.04;            %% is v in publication
% v=0.29;             %% is w in publication
%% in-plane maneuver is successful if deputy is within  10m in z-direction (x in publication) and 15m in x direction (y in publication).
%% out-of-plane maneuver is success if y-direction (z in publication) distance is 1m and velocity is 1cm/s (I dont know what velocity)
%% maneuver duration for the case described is: 4.2+14.2

  totalTime         =4*90*60;   %% simulation period (approximate multiples of orbit periods),[s]
  compStep          =9;          %% computational step size in [s]
  lengthControlLoop =90;         %% in [s]
  timetemp  =0:compStep:lengthControlLoop; %% duration and interpolation timestep for each control loop
  wakeAerodynamics=0;             %% use of wake aerodynamics
  masterSatellite=2;              %% if 0 no master, 1 active master, 2 passive master,

  
  sstDesiredFunction=@IRSRendezvousDesired;
  windOn      =1;
  sunOn       =0;
  deltaAngle  =45;              %% roll,pitch,yaw angle resolution
  ns          =1;               %% number of satellites
  satelliteMass=10;             %% [kg]
  argumentOfPerigeeAtTe0=0;     %% not used yet
  trueAnomalyAtTe0=0;           %% not used yet
  %% other constants
  [~,~,radiusOfEarth,~,~]=orbitalproperties(500000); %% fake value is passed to get the physics value only
  altitude=6778137-radiusOfEarth;
  [rho,v,radiusOfEarth,mu,MeanMotion]=orbitalproperties(altitude);
  
  %% satelliteshapeproperties, number of 10cmx10cm faces to x,y,z(body coordinates, normally aligned with):
  panels=[0 0 2]; 
  panelSurface=1.1;                %% [m^2]  
  %% initial condition
  sstTemp=zeros(9,ns,size(timetemp,2));
end
