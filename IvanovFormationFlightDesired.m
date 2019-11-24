function [sstDesired]=IvanovFormationFlightDesired(time,meanMotion,i,goFoFli,SSCoeff,SSParameters,meanAnomalyOffSet)
%% desired solution for Ivanov
%% input variables
%%  time
%%  meanMotion    = mean motion, [rad/s], double
%%  i
%%  goFoFli
%%  SSCoeff
%%  SSParameters
%%  meanAnomalyOffSet
%% output variables
%%  sstDesired
%% to do
%! change this to use the generalized formulation as in cluxterDesired

  sstDesired=zeros(9,size(time,2));
  %% analytical solution according to Ivanov
  switch goFoFli
    case 1
      A=100;    D=115;
    case 2
      factor=1000;
      A=factor*100;    D=factor*115;
  end
  switch i
    case 1
      sstDesired(1,:)=-D;
    case 2
      sstDesired(1,:)=D;
    case 3
      sstDesired(1,:)=2*A*        cos(meanMotion*(time)-acos(1/3));  
      sstDesired(2,:)=  A*sqrt(3)*sin(meanMotion*(time));
      sstDesired(3,:)=  A*        sin(meanMotion*(time)-acos(1/3));
      sstDesired(4,:)=2*A*       -sin(meanMotion*(time)-acos(1/3))*meanMotion;  
      sstDesired(5,:)=  A*sqrt(3)*cos(meanMotion*(time))*meanMotion;
      sstDesired(6,:)=  A*        cos(meanMotion*(time)-acos(1/3))*meanMotion;
    case 4
      sstDesired(1,:)=2*A*        cos(meanMotion*(time));
      sstDesired(2,:)=  A*sqrt(3)*sin(meanMotion*(time)+acos(1/3));
      sstDesired(3,:)=  A*        sin(meanMotion*(time));
      sstDesired(4,:)=2*A*       -sin(meanMotion*(time))*meanMotion;
      sstDesired(5,:)=  A*sqrt(3)*cos(meanMotion*(time)+acos(1/3))*meanMotion;
      sstDesired(6,:)=  A*        cos(meanMotion*(time))*meanMotion;
  end  
end