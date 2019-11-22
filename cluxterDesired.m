function sstDesired=cluxterDesired(time,MM,i,goFoFli,SSC,SSP,MAO)
%% cluxterDesired does ....
%% it follows the principles of 
%% C. Traub et al., “On the exploitation of differential aerodynamic lift and drag as a means to
%% control satellite formation flight,” CEAS Sp. J., 2019
%% but renames variables and uses Ivanov's coordinate system (x-ram,y-orbital plane vector,z-zenith)

%% input variables:
%% time             = time variable, used only for formation reconfiguration [s], double
%% MM               =  mean motion (shortened variable name for better readibility)[rad/s], double %% this should be [deg/s] in the future
%% i                = which satellite number in formation, int
%% goFoFli          = switch, used only in X
%% SSC              =SSCoefficient, shortened variable name for better readability
%% SSP              =SSParameters, shortened variable name for better readability
%% MAO              =meanAnomalyOffSet(shortened variable name for better readibility)
%% output variables:
%% sstDesired       = desired state vector, array of double, 9 state variables, length dependend on length of time

%% desired solution for Cluxter mission

    sstDesired=zeros(9,size(time,2));
    %! remove the following overwriting of the SSC asap
    SSC=1;
    
    A=2*SSC/(2-SSC^2);
    B=(2-5*SSC^2)/2/SSC;
    D=sqrt(3*SSC^2-2);
 
    sstDesired(1,:)= SSP(2)   + B*MM*SSP(6)*time + SSP(1)*cos(sqrt(2*SSC/A)*MM*time+SSP(7)+MAO)                  - SSP(5)*sqrt(2*SSC*A)*sin(sqrt(2*SSC/A)*MM*time+SSP(7)+MAO);
    sstDesired(2,:)=                               SSP(3)*cos(D*MM*time+SSP(8)+MAO)                              + SSP(4)/D/MM*         sin(D*MM*time+SSP(8)+MAO);
    sstDesired(3,:)= SSP(6)                      + SSP(5)*cos(sqrt(2*SSC/A)*MM*time+SSP(7)+MAO)                  + SSP(1)/sqrt(2*SSC*A)*sin(sqrt(2*SSC/A)*MM*time+SSP(7)+MAO);

    sstDesired(4,:)=            B*MM*SSP(6)       -SSP(1)*sin(sqrt(2*SSC/A)*MM*time+SSP(7)+MAO)*sqrt(2*SSC/A)*MM - SSP(5)*sqrt(2*SSC*A)*cos(sqrt(2*SSC/A)*MM*time+SSP(7)+MAO)*sqrt(2*SSC/A)*MM;
    sstDesired(5,:)=                              -SSP(3)*sin(D*MM*time+SSP(8)+MAO)*D*MM                         + SSP(4)/D/MM*         cos(D*MM*time+SSP(8)+MAO)*D*MM;
    sstDesired(6,:)=                              -SSP(5)*sin(sqrt(2*SSC/A)*MM*time+SSP(7)+MAO)*sqrt(2*SSC*A)*MM + SSP(1)/sqrt(2*SSC/A)*cos(sqrt(2*SSC/A)*MM*time+SSP(7)+MAO)*sqrt(2*SSC/A)*MM;
end