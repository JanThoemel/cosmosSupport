function [sstTemp,ns,altitude,panels,rho,v,radiusOfEarth,meanMotion,mu,satelliteMass,panelSurface,...
          sstDesiredFunction,windOn,sunOn,deltaAngle,timetemp,totalTime,wakeAerodynamics,masterSatellite]...
          =IvanovFormationFlightInitial()
%% initial conditions for Ivanov

  totalTime         =50*90*60;   %% simulation period (approximate multiples of orbit periods),[s]
  compStep          =9;          %% computational step size in [s]
  lengthControlLoop =90;         %% in [s]
  timetemp  =0:compStep:lengthControlLoop; %% duration and interpolation timestep for each control loop
  wakeAerodynamics=0;             %% use of wake aerodynamics
  masterSatellite=0;              %% if 0 no master, 1 active master, 2 passive master,

  sstDesiredFunction=@IvanovFormationFlightDesired;
  windOn=1;
  sunOn=0;
  deltaAngle        =45;                   %% roll,pitch,yaw angle resolution
  ns=4;
  satelliteMass=1;
  altitude=340000;                %% in m
  argumentOfPerigeeAtTe0=0;       %% not used yet
  trueAnomalyAtTe0=0;             %% not used yet
  ejectionVelocity=0.5;           %% m/s
  timeBetweenEjections=10;        %% in s
  %startSecondPhase=8*2*pi/MeanMotion;  %% in s
  panelSurface=0.01;              %% m^2  

  %% other constants
  [rho,v,radiusOfEarth,mu,meanMotion]=orbitalproperties(altitude);
  r0=radiusOfEarth+altitude;    %% in m
  panels=[0 0 2]; 
  sstTemp=zeros(9,ns,size(timetemp,2));
  for i=2:ns
      sstTemp(1,i,1)=(i-2)*ejectionVelocity*timeBetweenEjections; %% x position
      sstTemp(4,i,1)=0;            %% u velocity
      sstTemp(7,i,1)=0;           %% alpha
      sstTemp(8,i,1)=0;           %% beta
      sstTemp(9,i,1)=0;           %% gamma
  end
end