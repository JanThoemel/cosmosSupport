function  [ssttemptemp,controlVector]=HCWEquation(IR,P,A,B,deltat,sst0,e,windPressure,...
          solarPressure,alphas,betas,gammas,aeropressureforcevector,solarpressureforcevector,...
          oldAlphaOpt,oldBetaOpt,oldGammaOpt,refSurf,satelliteMass,wakeAerodynamics,masterSatellite,...
          currentTime0,radiusOfEarth,altitude,meanMotion)

%% HCWequation, Hill-Clohessy-Wiltshire equation
%% this function is doing X
%% input variables:
%% variable name,   explanation,  physical dimension,   data type
%% IR
%% P
%% A
%% B
%% deltat
%% sst0
%% e
%% windPressure
%% solarPressure
%% alphas
%% betas
%% gammas
%% aeropressureforcevector
%% solarpressureforcevector
%% oldAlphaOpt
%% oldBetaOpt
%% oldGammaOpt
%% refSurf
%% satelliteMass
%% wakeAerodynamics
%% masterSatellite
%% currentTime0
%% radiusOfEarth
%% altitude
%% meanMotion
%% output variables:
%% variable name,   explanation,  physical dimension,   data type
%% ssttemptemp
%% controlVector

  usedTotalForceVector=zeros(3,size(alphas,2),size(betas,2),size(gammas,2));
  rotatedSunForceVector=zeros(3,size(alphas,2),size(betas,2),size(gammas,2));

  %% compute control vector
  controlVector=-IR*B'*P*e;

  solarForceVectorOfMaster=[0 0 0]';
  maxSolarForce=[0 0 0]';

  rotatedSunForceVector=solarpressureforcevector;
%{
  %% rotate sunforcevector if necessary
  if solarPressure>0 && wrapTo2Pi(meanMotion*currentTime0)*180/pi<acos(radiusOfEarth/(radiusOfEarth+altitude))...
    && wrapTo2Pi(meanMotion*currentTime0)*180/pi<360-acos(radiusOfEarth/(radiusOfEarth+altitude))
      %% satellite is in sun light
    for k=1:size(gammas,2)
     for j=1:size(betas,2)
       for i=1:size(alphas,2)
         rotatedSunForceVector(:,i,j,k)=roty( wrapTo2Pi(meanMotion*currentTime0)*180/pi )*solarpressureforcevector(:,i,j,k);
         maxSolarForce                 =roty(wrapTo2Pi(meanMotion*currentTime0)*180/pi)*[0 0 -1]'*2.8*solarPressure*refSurf;
         solarForceVectorOfMaster      =roty(wrapTo2Pi(meanMotion*currentTime0)*180/pi)*[0 0 -1]'*1/2*2.8*solarPressure*refSurf;
       end
     end
    end
  else %% satellite is eclipsed
    rotatedSunForceVector=0*solarpressureforcevector;
  end  
  %}
  
  if masterSatellite==0 %% no master satellite
    %aeropressureforcevector
    %rotatedSunForceVector
    for k=1:size(gammas,2)
      for j=1:size(betas,2)
        for i=1:size(alphas,2)

          usedTotalForceVector(:,i,j,k) =aeropressureforcevector(:,i,j,k)+rotatedSunForceVector(:,i,j,k);
        end
      end
    end    
  elseif masterSatellite==1 %% active master satellite           
    maxforceVectorOfMaster=-2.8*windPressure*refSurf*[1 0 0]'-maxSolarForce;
    if wakeAerodynamics %% wake aerodynamics
      if abs(sst0(2))<1.5 && abs(sst0(3))< 1.5 %% is sat2 aligned with sat1?
        if sst0(1) <= 0 %% is sat2 before sat1?
          maxforceVectorOfMaster=maxforceVectorOfMaster/10;
        else            %% is sat2 before sat2?
          aeropressureforcevector(1,:,:,:)=aeropressureforcevector(1,:,:,:)/10;
        end
      end
    end
    for k=1:size(gammas,2)
      for j=1:size(betas,2)
        for i=1:size(alphas,2)
          usedTotalForceVector(:,i,j,k)=2*aeropressureforcevector(:,i,j,k)+rotatedSunForceVector(:,i,j,k)-maxforceVectorOfMaster(:);
        end
      end
    end
  %% passive master; wake aerodynamics not implemented
  elseif masterSatellite==2     %% passive master                  
    forceVectorOfMaster=-1/2*2.8*windPressure*refSurf*[1 0 0]'-solarForceVectorOfMaster;  
    for k=1:size(gammas,2)
      for j=1:size(betas,2)
        for i=1:size(alphas,2)
          usedTotalForceVector(:,i,j,k)=aeropressureforcevector(:,i,j,k)+rotatedSunForceVector(:,i,j,k)-forceVectorOfMaster(:);
        end
      end
    end
  else
    fprintf('\n error in defining master satellite \n');
    input('error');
  end  %% passive/active master switch

  [forceVector,alphaOpt,betaOpt,gammaOpt]=findBestAerodynamicAngles(usedTotalForceVector,controlVector,alphas,betas,gammas,oldAlphaOpt,oldBetaOpt,oldGammaOpt);
  %% solve ODE with backward Euler step
  ssttemptemp(1:6)=(A*sst0(1:6)+B*forceVector/satelliteMass)*deltat+sst0(1:6); 
 
  ssttemptemp(7:9)=[alphaOpt betaOpt gammaOpt]';
  %fprintf('f')
 % forceVector'
end