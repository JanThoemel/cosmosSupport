function sstDesired=cluxterDesired(time,meanMotion,i,goFoFli)
%% cluxterDesired does ....
%% input variables:
%% time             = time variable, used only for formation reconfiguration [s], double
%% meanMotion       = mean motion [rad/s], double %% this should be [deg/s] in the future
%% i                = which satellite number in formation, int
%% goFoFli          = switch, used only in X
%% output variables:
%% sstDesired       = desired state vector, array of double, 9 state variables, length dependend on length of time

%% desired solution for Cluxter mission
  sstDesired=zeros(9,1,size(time,2));
  %% analytical solution according to Ivanov
  if time<100000*90*60
    A=7;    D=15;
  else
    A=14;    D=30;    
  end
  switch i
    case 1
      sstDesired(1,1,:)=-D;
    case 2
      sstDesired(1,1,:)=2*A*        cos(meanMotion*(time)-acos(1/3)+pi/2);  
      sstDesired(2,1,:)=  A*sqrt(3)*sin(meanMotion*(time)+pi/2);
      sstDesired(3,1,:)=  A*        sin(meanMotion*(time)-acos(1/3)+pi/2);
      sstDesired(4,1,:)=2*A*       -sin(meanMotion*(time)-acos(1/3)+pi/2)*meanMotion;  
      sstDesired(5,1,:)=  A*sqrt(3)*cos(meanMotion*(time)+pi/2)*meanMotion;
      sstDesired(6,1,:)=  A*        cos(meanMotion*(time)-acos(1/3)+pi/2)*meanMotion;
    case 3
      sstDesired(1,1,:)=D;
  end
end