function [sstTemp,ns,altitude,panels,rho,Tatmos,v,radiusOfEarth,meanMotion,mu,satelliteMass,panelSurface,...
          sstDesiredFunction,windOn,sunOn,deltaAngle,timetemp,totalTime,wakeAerodynamics,masterSatellite,...
          SSCoeff,inclination,SSParameters,meanAnomalyOffSet]=IvanovFormationFlightInitial()
%% initial conditions for Ivanov
%% input variables
%%  none
%% output variables
%%  sstTemp
%%  ns
%%  altitude
%%  panels
%%  rho
%%  Tatmos
%%  v
%%  radiusOfEarth
%%  meanMotion    = mean motion, [rad/s], double
%%  mu
%%  satelliteMass
%%  panelSurface
%%  sstDesiredFunction
%%  windOn
%%  sunOn
%%  deltaAngle
%%  timetemp
%%  totalTime
%%  wakeAerodynamics
%%  masterSatellite
%%  SSCoeff
%%  inclination
%%  SSParameters
%%  meanAnomalyOffSet      

  totalTime         =50*90*60; %% simulation period (approximate multiples of orbit periods),[s], it is not relevant for cosmosFS
  compStep          =9;         %% computational step size in [s]
  lengthControlLoop =90;        %% in [s]
  timetemp  =0:compStep:lengthControlLoop; %% duration and interpolation timestep for each control loop
  wakeAerodynamics  =0;         %% use of wake aerodynamics
  masterSatellite   =0;         %% if 0 no master, 1 active master, 2 passive master,

  sstDesiredFunction=@IvanovFormationFlightDesired;
  windOn            =1;
  sunOn             =0;
  deltaAngle        =30;        %% roll,pitch,yaw angle resolution
  ns                =4;
  satelliteMass     =1;
  altitude          =340e3 ;    %% [m]
  argumentOfPerigeeAtTe0=0;     %% not used yet
  trueAnomalyAtTe0  =0;         %% not used yet
  ejectionVelocity  =0.5;       %% m/s
  timeBetweenEjections=10;      %% in s
  %startSecondPhase=8*2*pi/MeanMotion;  %% in s
  panelSurface      =0.01;      %% m^2  

  %% other constants
  %[rho,Tatmos,v,radiusOfEarth,mu,meanMotion,~]=orbitalproperties(altitude);
  [rho,Tatmos,v,radiusOfEarth,mu,meanMotion,SSOinclination,~]=orbitalproperties(altitude);
  %%% meanMotion is in [rad/s]
  fprintf('\noriginal atmospheric density (rho) is overwritten with %1.3e Ivanovs (rho) %1.3e \n',rho,1e-11);
  rho=10e-11;

  r0=radiusOfEarth+altitude;    %% [m]
  inclination=SSOinclination;   %%
  %inclination=0                %% manualinclination
  
  SSCoeff=1;%sqrt(1+3*J2*radiusOfEarth^2/8/r0^2*(1+3*cosd(2*inclination))) ;
    
  panels=[0 0 2]; 
  sstTemp=zeros(9,ns,size(timetemp,2));
  for i=1:ns
    sstTemp(1,i,1)=(i-1)*ejectionVelocity*timeBetweenEjections; %% x position
    sstTemp(4,i,1)=0;          %% u velocity
    sstTemp(7,i,1)=0;          %% alpha
    sstTemp(8,i,1)=0;          %% beta
    sstTemp(9,i,1)=0;          %% gamma
  end
  
  %% parameters for desired conditions  
  meanAnomalyOffSet =pi/4; %% for pi/2 satellite cross on the poles; for 0 satellite crosses at equator 
  numberOfModes     =10;
  SSParameters      =zeros(6,ns,numberOfModes);
end