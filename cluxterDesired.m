function sstDesired=cluxterDesired(timetemptemp,MeanMotion,i,goFoFli)
%% desired solution for Cluxter mission
  sstDesired=zeros(9,1,size(timetemptemp,2));
  %% analytical solution according to Ivanov
  if timetemptemp<100000*90*60
    AAA=7;    DDD=15;
  else
    AAA=14;    DDD=30;    
  end
  switch i
    case 1
      sstDesired(1,1,:)=-DDD;
    case 2
      sstDesired(1,1,:)=2*AAA*        cos(MeanMotion*(timetemptemp)-acos(1/3)+pi/2);  
      sstDesired(2,1,:)=  AAA*sqrt(3)*sin(MeanMotion*(timetemptemp)+pi/2);
      sstDesired(3,1,:)=  AAA*        sin(MeanMotion*(timetemptemp)-acos(1/3)+pi/2);
      sstDesired(4,1,:)=2*AAA*       -sin(MeanMotion*(timetemptemp)-acos(1/3)+pi/2)*MeanMotion;  
      sstDesired(5,1,:)=  AAA*sqrt(3)*cos(MeanMotion*(timetemptemp)+pi/2)*MeanMotion;
      sstDesired(6,1,:)=  AAA*        cos(MeanMotion*(timetemptemp)-acos(1/3)+pi/2)*MeanMotion;
    case 3
      sstDesired(1,1,:)=DDD;
  end
end