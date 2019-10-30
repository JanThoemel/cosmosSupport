function  [sstTemp,ns,altitude,panels,rho,Tatmos,v,radiusOfEarth,MeanMotion,mu,satelliteMass,panelSurface,...
          sstDesiredFunction,windOn,sunOn,deltaAngle,timetemp,totalTime,wakeAerodynamics,masterSatellite,...
          SSCoeff,SSParameters,meanAnomalyOffSet]=cluxterInitial()
%% initial conditions for Cluxter Mission


%% simulation time constants
  totalTime         =30*90*60;   %% simulation period (approximate multiples of orbit periods),[s]
  compStep          =3;          %% computational step size in [s]
  lengthControlLoop =90;         %% in [s]
  timetemp  =0:compStep:lengthControlLoop; %% duration and interpolation timestep for each control loop
  %fprintf('\nnumber of float variables for each control loop: %d, size %f kbyte',size(timetemp,2)*(1+6+6+3)+6+1,(size(timetemp,2)*(1+6+6+3)+6+1)*4/1024); %% size of time vector, state vector, desired state vector and anglevector, C and MeanMotion of analytical solution
  wakeAerodynamics=0;             %% use of wake aerodynamics
  masterSatellite=0;              %% if 0 no master, 1 active master, 2 passive master,


  sstDesiredFunction=@cluxterDesired;
  windOn       =1;
  sunOn        =1;
  J2On         =0; %% not yet implemented
  if sunOn
    fprintf('\n only dusk-dawn orbits for the time being! Arbitrary orbits require rotation of sunlight vector and accounting for eclipse\n');
  end
  deltaAngle   =30;                   %% roll,pitch,yaw angle resolution

  ns=3;
  satelliteMass=1.5;
  altitude=500000;                %% [m]
  argumentOfPerigeeAtTe0=0;       %% not used yet
  trueAnomalyAtTe0=0;             %% not used yet
  ejectionVelocity=0;%0.01;          %% [m/s]
  timeBetweenEjections=0.01;      %% [s]
  %startSecondPhase=8*2*pi/MeanMotion;  %% in s
  panelSurface=0.01;              %% [m^2]
  panels=[0 0 2]; 

  %% other constants
  [rho,Tatmos,v,radiusOfEarth,mu,MeanMotion,SSOinclination,J2]=orbitalproperties(altitude);
  r0=radiusOfEarth+altitude;    %% in m

  inclination=SSOinclination;   %%
  %inclination=0                %% manualinclination
  
  SSCoeff=sqrt(1+3*J2*radiusOfEarth^2/8/r0^2*(1+3*cosd(2*inclination)))^J2On; 
   
  %% initial conditions
  sstTemp=zeros(9,ns,size(timetemp,2));
  for i=1:ns
      sstTemp(1,i,1)=(i-1)*ejectionVelocity*timeBetweenEjections; %% x position
      sstTemp(4,i,1)=0;           %% u velocity
      sstTemp(7,i,1)=0;           %% alpha
      sstTemp(8,i,1)=0;           %% beta
      sstTemp(9,i,1)=0;           %% gamma
  end
  
  %% parameters for desired conditions  
  meanAnomalyOffSet=0;% %% for 0 satellite cross on the poles; for pi/2 satellite crosses at equator
  numberOfModes=10;
  SSParameters=zeros(6,ns,numberOfModes);
  Aold=7; Dold=15;

  %% generalized analytical solution, from Traub with variable renaming and accounting for different coordinate system, accelerations do not belong here
  SSParameters(1,1,1)=0;     %%
  SSParameters(2,1,1)=-Dold;     %%
  SSParameters(3,1,1)=0;     %%
  SSParameters(4,1,1)=0;     %%
  SSParameters(5,1,1)=0;     %%
  SSParameters(6,1,1)=0;     %%

  SSParameters(1,2,1)=2*Aold;     %% xmax
  SSParameters(2,2,1)=0;              %% xpermanent.only here, could assume any value
  SSParameters(3,2,1)=Aold*sqrt(3);   %% ymax
  SSParameters(4,2,1)=0;              %% ymaxdot0
  SSParameters(5,2,1)=0;%Aold;           %% zmax
  SSParameters(6,2,1)=0;              %% zpermanent, it is a permanent nadir-zenit movement should always be zero because it implies a permant ram-antiram movement. should this ever be non-zero, all has to be verified

  SSParameters(1,3,1)=0;     %%
  SSParameters(2,3,1)=Dold;     %%
  SSParameters(3,3,1)=0;     %%
  SSParameters(4,3,1)=0;     %%
  SSParameters(5,3,1)=0;     %%
  SSParameters(6,3,1)=0;     %%

  
end



%{
%%%%%%%%THIS WORKS PERFECTLY WITH A THETARANGE OF 90 DEG
%% simulation time constants
  totalTime         =150*90*60;   %% simulation period (approximate multiples of orbit periods),[s]
  compStep          =9;          %% computational step size in [s]
  lengthControlLoop =90;         %% in [s]
  timetemp  =0:compStep:lengthControlLoop; %% duration and interpolation timestep for each control loop
  %fprintf('\nnumber of float variables for each control loop: %d, size %f kbyte',size(timetemp,2)*(1+6+6+3)+6+1,(size(timetemp,2)*(1+6+6+3)+6+1)*4/1024); %% size of time vector, state vector, desired state vector and anglevector, C and MeanMotion of analytical solution
  wakeAerodynamics=0;             %% use of wake aerodynamics
  masterSatellite=0;              %% if 0 no master, 1 active master, 2 passive master,


  sstDesiredFunction=@cluxterDesired;
  windOn       =1;
  sunOn        =0;
  if sunOn
    fprintf('\n only dusk-dawn orbits for the time being! Arbitrary orbits require rotation of sunlight vector and accounting for eclipse\n');
  end
  deltaAngle   =45;                   %% roll,pitch,yaw angle resolution

  ns=3;
  satelliteMass=1;
  altitude=500000;                %% in m
  argumentOfPerigeeAtTe0=0;       %% not used yet
  trueAnomalyAtTe0=0;             %% not used yet
  ejectionVelocity=0.01;          %% m/s
  timeBetweenEjections=0.01;      %% in s
  %startSecondPhase=8*2*pi/MeanMotion;  %% in s
  panelSurface=0.01;              %% m^2  
  panels=[0 0 4]; 

  %% other constants
  [rho,v,radiusOfEarth,mu,MeanMotion]=orbitalproperties(altitude);
  r0=radiusOfEarth+altitude;    %% in m
  
  sstTemp=zeros(9,ns,size(timetemp,2));
  for i=1:ns
      sstTemp(1,i,1)=(i-1)*ejectionVelocity*timeBetweenEjections; %% x position
      sstTemp(4,i,1)=0;           %% u velocity
      sstTemp(7,i,1)=0;           %% alpha
      sstTemp(8,i,1)=0;           %% beta
      sstTemp(9,i,1)=0;           %% gamma
  end


%}





%{
THIS WORK VERY WELL WITH A 
-thetaRange OF 45 DEGREES
-simplified GSI and LI
-it seems not to work with 700km altitude

%% simulation time constants
  totalTime         =20*90*60;   %% simulation period (approximate multiples of orbit periods),[s]
  compStep          =3;          %% computational step size in [s]
  lengthControlLoop =90;         %% in [s]
  timetemp  =0:compStep:lengthControlLoop; %% duration and interpolation timestep for each control loop
  %fprintf('\nnumber of float variables for each control loop: %d, size %f kbyte',size(timetemp,2)*(1+6+6+3)+6+1,(size(timetemp,2)*(1+6+6+3)+6+1)*4/1024); %% size of time vector, state vector, desired state vector and anglevector, C and MeanMotion of analytical solution
  wakeAerodynamics=0;             %% use of wake aerodynamics
  masterSatellite=0;              %% if 0 no master, 1 active master, 2 passive master,


  sstDesiredFunction=@cluxterDesired;
  windOn       =1;
  sunOn        =1;
  J2On         =0; %% not yet implemented
  if sunOn
    fprintf('\n only dusk-dawn orbits for the time being! Arbitrary orbits require rotation of sunlight vector and accounting for eclipse\n');
  end
  deltaAngle   =30;                   %% roll,pitch,yaw angle resolution

  ns=3;
  satelliteMass=1;
  altitude=500000;                %% [m]
  argumentOfPerigeeAtTe0=0;       %% not used yet
  trueAnomalyAtTe0=0;             %% not used yet
  ejectionVelocity=0;%0.01;          %% [m/s]
  timeBetweenEjections=0.01;      %% [s]
  %startSecondPhase=8*2*pi/MeanMotion;  %% in s
  panelSurface=0.01;              %% [m^2]
  panels=[0 0 2]; 

  %% other constants
  [rho,v,radiusOfEarth,mu,MeanMotion]=orbitalproperties(altitude);
  r0=radiusOfEarth+altitude;    %% in m
  
  sstTemp=zeros(9,ns,size(timetemp,2));
  for i=1:ns
      sstTemp(1,i,1)=(i-1)*ejectionVelocity*timeBetweenEjections; %% x position
      sstTemp(4,i,1)=0;           %% u velocity
      sstTemp(7,i,1)=0;           %% alpha
      sstTemp(8,i,1)=0;           %% beta
      sstTemp(9,i,1)=0;           %% gamma
  end
%}




%{
28/10/2019, this works with simple GSI and thetarange 45

%% simulation time constants
  totalTime         =100*90*60;   %% simulation period (approximate multiples of orbit periods),[s]
  compStep          =3;          %% computational step size in [s]
  lengthControlLoop =90;         %% in [s]
  timetemp  =0:compStep:lengthControlLoop; %% duration and interpolation timestep for each control loop
  %fprintf('\nnumber of float variables for each control loop: %d, size %f kbyte',size(timetemp,2)*(1+6+6+3)+6+1,(size(timetemp,2)*(1+6+6+3)+6+1)*4/1024); %% size of time vector, state vector, desired state vector and anglevector, C and MeanMotion of analytical solution
  wakeAerodynamics=0;             %% use of wake aerodynamics
  masterSatellite=0;              %% if 0 no master, 1 active master, 2 passive master,


  sstDesiredFunction=@cluxterDesired;
  windOn       =1;
  sunOn        =1;
  J2On         =0; %% not yet implemented
  if sunOn
    fprintf('\n only dusk-dawn orbits for the time being! Arbitrary orbits require rotation of sunlight vector and accounting for eclipse\n');
  end
  deltaAngle   =30;                   %% roll,pitch,yaw angle resolution

  ns=3;
  satelliteMass=1;
  altitude=550000;                %% [m]
  argumentOfPerigeeAtTe0=0;       %% not used yet
  trueAnomalyAtTe0=0;             %% not used yet
  ejectionVelocity=0;%0.01;          %% [m/s]
  timeBetweenEjections=0.01;      %% [s]
  %startSecondPhase=8*2*pi/MeanMotion;  %% in s
  panelSurface=0.01;              %% [m^2]
  panels=[0 0 2]; 

  %% other constants
  [rho,v,radiusOfEarth,mu,MeanMotion]=orbitalproperties(altitude);
  r0=radiusOfEarth+altitude;    %% in m
  
  sstTemp=zeros(9,ns,size(timetemp,2));
  for i=1:ns
      sstTemp(1,i,1)=(i-1)*ejectionVelocity*timeBetweenEjections; %% x position
      sstTemp(4,i,1)=0;           %% u velocity
      sstTemp(7,i,1)=0;           %% alpha
      sstTemp(8,i,1)=0;           %% beta
      sstTemp(9,i,1)=0;           %% gamma
  end


%}


%{
THIS SEEMS TO WORK
WITH simple SSI, GSI
BUT NOT WITH SUNOFF 
 
%% simulation time constants
  totalTime         =150*90*60;   %% simulation period (approximate multiples of orbit periods),[s]
  compStep          =3;          %% computational step size in [s]
  lengthControlLoop =90;         %% in [s]
  timetemp  =0:compStep:lengthControlLoop; %% duration and interpolation timestep for each control loop
  %fprintf('\nnumber of float variables for each control loop: %d, size %f kbyte',size(timetemp,2)*(1+6+6+3)+6+1,(size(timetemp,2)*(1+6+6+3)+6+1)*4/1024); %% size of time vector, state vector, desired state vector and anglevector, C and MeanMotion of analytical solution
  wakeAerodynamics=0;             %% use of wake aerodynamics
  masterSatellite=0;              %% if 0 no master, 1 active master, 2 passive master,


  sstDesiredFunction=@cluxterDesired;
  windOn       =1;
  sunOn        =1;
  J2On         =0; %% not yet implemented
  if sunOn
    fprintf('\n only dusk-dawn orbits for the time being! Arbitrary orbits require rotation of sunlight vector and accounting for eclipse\n');
  end
  deltaAngle   =30;                   %% roll,pitch,yaw angle resolution

  ns=3;
  satelliteMass=1.5;
  altitude=600000;                %% [m]
  argumentOfPerigeeAtTe0=0;       %% not used yet
  trueAnomalyAtTe0=0;             %% not used yet
  ejectionVelocity=0;%0.01;          %% [m/s]
  timeBetweenEjections=0.01;      %% [s]
  %startSecondPhase=8*2*pi/MeanMotion;  %% in s
  panelSurface=0.01;              %% [m^2]
  panels=[0 0 2]; 

  %% other constants
  [rho,v,radiusOfEarth,mu,MeanMotion]=orbitalproperties(altitude);
  r0=radiusOfEarth+altitude;    %% in m
  
  sstTemp=zeros(9,ns,size(timetemp,2));
  for i=1:ns
      sstTemp(1,i,1)=(i-1)*ejectionVelocity*timeBetweenEjections; %% x position
      sstTemp(4,i,1)=0;           %% u velocity
      sstTemp(7,i,1)=0;           %% alpha
      sstTemp(8,i,1)=0;           %% beta
      sstTemp(9,i,1)=0;           %% gamma
  end
end





%}