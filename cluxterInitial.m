function [sstTemp,ns,altitude,panels,rho,v,radiusOfEarth,MeanMotion,mu,satelliteMass,panelSurface,...
  sstDesiredFunction,windOn,sunOn,deltaAngle,timetemp,totalTime,wakeAerodynamics,masterSatellite]=cluxterInitial()
%% initial conditions for Cluxter Mission


%% simulation time constants
  totalTime         =50*90*60;   %% simulation period (approximate multiples of orbit periods),[s]
  compStep          =9;          %% computational step size in [s]
  lengthControlLoop =90;         %% in [s]
  timetemp  =0:compStep:lengthControlLoop; %% duration and interpolation timestep for each control loop
  %fprintf('\nnumber of float variables for each control loop: %d, size %f kbyte',size(timetemp,2)*(1+6+6+3)+6+1,(size(timetemp,2)*(1+6+6+3)+6+1)*4/1024); %% size of time vector, state vector, desired state vector and anglevector, C and MeanMotion of analytical solution
  wakeAerodynamics=0;             %% use of wake aerodynamics
  masterSatellite=0;              %% if 0 no master, 1 active master, 2 passive master,


  sstDesiredFunction=@cluxterDesired;
  windOn       =1;
  sunOn        =1;
  deltaAngle   =45;                   %% roll,pitch,yaw angle resolution

  ns=3;
  satelliteMass=1;
  altitude=300000;                %% in m
  argumentOfPerigeeAtTe0=0;       %% not used yet
  trueAnomalyAtTe0=0;             %% not used yet
  ejectionVelocity=1;           %% m/s
  timeBetweenEjections=4;        %% in s
  %startSecondPhase=8*2*pi/MeanMotion;  %% in s
  panelSurface=0.01;              %% m^2  

  %% other constants
  [rho,v,radiusOfEarth,mu,MeanMotion]=orbitalproperties(altitude);
  r0=radiusOfEarth+altitude;    %% in m
  panels=[0 0 2]; 
  sstTemp=zeros(9,ns,size(timetemp,2));
  for i=1:ns
      sstTemp(1,i,1)=(i-1)*ejectionVelocity*timeBetweenEjections; %% x position
      sstTemp(4,i,1)=0;           %% u velocity
      sstTemp(7,i,1)=0;           %% alpha
      sstTemp(8,i,1)=0;           %% beta
      sstTemp(9,i,1)=0;           %% gamma
  end
end