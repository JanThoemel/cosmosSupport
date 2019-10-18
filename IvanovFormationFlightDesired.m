function [sstDesired]=IvanovFormationFlightDesired(time,meanMotion,i,mode)
%% desired solution for Ivanov
  sstDesired=zeros(9,1,size(time,2));
  %% analytical solution according to Ivanov
  switch mode
    case 1
      A=100;    D=115;
    case 2
      factor=2;
      A=factor*100;    D=factor*115;
  end
  switch i
    case 1
      sstDesired(1,1,:)=-D;
    case 2
      sstDesired(1,1,:)=D;
    case 3
      sstDesired(1,1,:)=2*A*        cos(meanMotion*(time)-acos(1/3));  
      sstDesired(2,1,:)=  A*sqrt(3)*sin(meanMotion*(time));
      sstDesired(3,1,:)=  A*        sin(meanMotion*(time)-acos(1/3));
      sstDesired(4,1,:)=2*A*       -sin(meanMotion*(time)-acos(1/3))*meanMotion;  
      sstDesired(5,1,:)=  A*sqrt(3)*cos(meanMotion*(time))*meanMotion;
      sstDesired(6,1,:)=  A*        cos(meanMotion*(time)-acos(1/3))*meanMotion;
    case 4
      sstDesired(1,1,:)=2*A*        cos(meanMotion*(time));
      sstDesired(2,1,:)=  A*sqrt(3)*sin(meanMotion*(time)+acos(1/3));
      sstDesired(3,1,:)=  A*        sin(meanMotion*(time));
      sstDesired(4,1,:)=2*A*       -sin(meanMotion*(time))*meanMotion;
      sstDesired(5,1,:)=  A*sqrt(3)*cos(meanMotion*(time)+acos(1/3))*meanMotion;
      sstDesired(6,1,:)=  A*        cos(meanMotion*(time))*meanMotion;
  end
end