function sstDesired=cluxterDesired(time,meanMotion,i,goFoFli,SSCoeff,SSParameters,meanAnomalyOffSet)
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
    A=7;    D=15;

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
 % fprintf('\n old\n')
 % sstDesired(:,1,1)
  
    sstDesired=zeros(9,size(time,2));
  SSCoeff=1;
    otherOffSet=-acos(1/3);
    
    A=2*SSCoeff/(2-SSCoeff^2);
    B=(2-5*SSCoeff^2)/2/SSCoeff;
    D=sqrt(3*SSCoeff^2-2);
 %   SSParameters
    switch i
      case 1
        sstDesired(1,:)=  SSParameters(2,goFoFli) + B*meanMotion*SSParameters(6,goFoFli)*time + SSParameters(1,goFoFli)*cos(sqrt(2*SSCoeff/A)*meanMotion*time+meanAnomalyOffSet) - SSParameters(5,goFoFli)*sqrt(2*SSCoeff*A)*sin(sqrt(2*SSCoeff/A)*meanMotion*time+meanAnomalyOffSet);
        sstDesired(2,:)=  SSParameters(3,goFoFli)*cos(D*meanMotion*time+meanAnomalyOffSet) + SSParameters(4,goFoFli)/D/meanMotion*sin(D*meanMotion*time+meanAnomalyOffSet);
        sstDesired(3,:)=  SSParameters(6,goFoFli) + SSParameters(5,goFoFli)*cos(sqrt(2*SSCoeff/A)*meanMotion*time+meanAnomalyOffSet) + SSParameters(1,goFoFli)/sqrt(2*SSCoeff*A)*sin(sqrt(2*SSCoeff/A)*meanMotion*time+meanAnomalyOffSet);

        sstDesired(4,:)=                                                                       -SSParameters(1,goFoFli)*sin(sqrt(2*SSCoeff/A)*meanMotion*time+meanAnomalyOffSet)*sqrt(2*SSCoeff/A)*meanMotion - SSParameters(5,goFoFli)*sqrt(2*SSCoeff*A)*cos(sqrt(2*SSCoeff/A)*meanMotion*time+meanAnomalyOffSet)*sqrt(2*SSCoeff/A)*meanMotion;
        sstDesired(5,:)=  -SSParameters(3,goFoFli)*sin(D*meanMotion*time + meanAnomalyOffSet)*D*meanMotion + SSParameters(4,goFoFli)/D/meanMotion*cos(D*meanMotion*time + meanAnomalyOffSet)*D*meanMotion;
        sstDesired(6,:)=                            -SSParameters(5,goFoFli)*sin(sqrt(2*SSCoeff/A)*meanMotion*time+meanAnomalyOffSet)*sqrt(2*SSCoeff*A)*meanMotion + SSParameters(1,goFoFli)/sqrt(2*SSCoeff/A)*cos(sqrt(2*SSCoeff/A)*meanMotion*time+meanAnomalyOffSet)*sqrt(2*SSCoeff/A)*meanMotion;
      case 2   
        sstDesired(1,:)=  SSParameters(2,goFoFli) + B*meanMotion*SSParameters(6,goFoFli)*time + SSParameters(1,goFoFli)*cos(sqrt(2*SSCoeff/A)*meanMotion*time+meanAnomalyOffSet) - SSParameters(5,goFoFli)*sqrt(2*SSCoeff*A)*sin(sqrt(2*SSCoeff/A)*meanMotion*time+meanAnomalyOffSet);
        sstDesired(2,:)=  SSParameters(3,goFoFli)*cos(D*meanMotion*time+meanAnomalyOffSet) + SSParameters(4,goFoFli)/D/meanMotion*sin(D*meanMotion*time+meanAnomalyOffSet);
        sstDesired(3,:)=  SSParameters(6,goFoFli) + SSParameters(5,goFoFli)*cos(sqrt(2*SSCoeff/A)*meanMotion*time+meanAnomalyOffSet) + SSParameters(1,goFoFli)/sqrt(2*SSCoeff*A)*sin(sqrt(2*SSCoeff/A)*meanMotion*time+meanAnomalyOffSet);

        sstDesired(4,:)=                                                                       -SSParameters(1,goFoFli)*sin(sqrt(2*SSCoeff/A)*meanMotion*time+meanAnomalyOffSet)*sqrt(2*SSCoeff/A)*meanMotion - SSParameters(5,goFoFli)*sqrt(2*SSCoeff*A)*cos(sqrt(2*SSCoeff/A)*meanMotion*time+meanAnomalyOffSet)*sqrt(2*SSCoeff/A)*meanMotion;
        sstDesired(5,:)=  -SSParameters(3,goFoFli)*sin(D*meanMotion*time + meanAnomalyOffSet)*D*meanMotion + SSParameters(4,goFoFli)/D/meanMotion*cos(D*meanMotion*time + meanAnomalyOffSet)*D*meanMotion;
        sstDesired(6,:)=                            -SSParameters(5,goFoFli)*sin(sqrt(2*SSCoeff/A)*meanMotion*time+meanAnomalyOffSet)*sqrt(2*SSCoeff*A)*meanMotion + SSParameters(1,goFoFli)/sqrt(2*SSCoeff/A)*cos(sqrt(2*SSCoeff/A)*meanMotion*time+meanAnomalyOffSet)*sqrt(2*SSCoeff/A)*meanMotion;
      case 3
        sstDesired(1,:)=  SSParameters(2,goFoFli) + B*meanMotion*SSParameters(6,goFoFli)*time + SSParameters(1,goFoFli)*cos(sqrt(2*SSCoeff/A)*meanMotion*time+meanAnomalyOffSet) - SSParameters(5,goFoFli)*sqrt(2*SSCoeff*A)*sin(sqrt(2*SSCoeff/A)*meanMotion*time+meanAnomalyOffSet);
        sstDesired(2,:)=  SSParameters(3,goFoFli)*cos(D*meanMotion*time+meanAnomalyOffSet) + SSParameters(4,goFoFli)/D/meanMotion*sin(D*meanMotion*time+meanAnomalyOffSet);
        sstDesired(3,:)=  SSParameters(6,goFoFli) + SSParameters(5,goFoFli)*cos(sqrt(2*SSCoeff/A)*meanMotion*time+meanAnomalyOffSet) + SSParameters(1,goFoFli)/sqrt(2*SSCoeff*A)*sin(sqrt(2*SSCoeff/A)*meanMotion*time+meanAnomalyOffSet);

        sstDesired(4,:)=                                                                       -SSParameters(1,goFoFli)*sin(sqrt(2*SSCoeff/A)*meanMotion*time+meanAnomalyOffSet)*sqrt(2*SSCoeff/A)*meanMotion - SSParameters(5,goFoFli)*sqrt(2*SSCoeff*A)*cos(sqrt(2*SSCoeff/A)*meanMotion*time+meanAnomalyOffSet)*sqrt(2*SSCoeff/A)*meanMotion;
        sstDesired(5,:)=  -SSParameters(3,goFoFli)*sin(D*meanMotion*time + meanAnomalyOffSet)*D*meanMotion + SSParameters(4,goFoFli)/D/meanMotion*cos(D*meanMotion*time + meanAnomalyOffSet)*D*meanMotion;
        sstDesired(6,:)=                            -SSParameters(5,goFoFli)*sin(sqrt(2*SSCoeff/A)*meanMotion*time+meanAnomalyOffSet)*sqrt(2*SSCoeff*A)*meanMotion + SSParameters(1,goFoFli)/sqrt(2*SSCoeff/A)*cos(sqrt(2*SSCoeff/A)*meanMotion*time+meanAnomalyOffSet)*sqrt(2*SSCoeff/A)*meanMotion;
    end
%    fprintf('\n new\n')
%    sstDesired(:,1)
%    input('');
end