function [forceVector,alphaOpt,betaOpt,gammaOpt]=findBestAttitude(totalForceVector,controlVector,alphas,betas,...
                                                 gammas,oldAlphaOpt,oldBetaOpt,oldGammaOpt) 
%findBestAerodynamicAngles to find best angle
%%%%%%% input variables
%% totalforcevector
%% controlvector
%% alphas
%% betas
%% gammas
%% oldAlphaOpt
%% oldBetaOpt
%% oldGammaOpt
%%%%%%% output variables
%% forcevector
%% alphaOpt
%% betaOpt
%% gammaOpt

  thetaRange=90;
  %% is this a universally good thetaRange?
  %% What is the relation to angle granularity?
  %% Do we need a formula like 1.5xgranularity,
  %% or a check that thetaRange is larger than granularity?
  
  theta=zeros(size(gammas,2),size(betas,2),size(alphas,2));
  phi=zeros(size(gammas,2),size(betas,2),size(alphas,2));
  l=1;
  for k=1:size(gammas,2) %% yaw
    for j=1:size(betas,2) %% pitch
      for i=1:size(alphas,2) %% roll
        
        %totalForceVector(:,i,j,k)'/norm(totalForceVector(:,i,j,k))
        %controlVector'/norm(controlVector)
        [theta(i,j,k),phi(i,j,k)]=thetaphi(totalForceVector(:,i,j,k)/norm(totalForceVector(:,i,j,k)),controlVector/norm(controlVector));        
        %theta(i,j,k)
        %phi(i,j,k)
        if theta(i,j,k)<thetaRange %% find all force vectors that have theta<thetarange deg
          %fprintf('-x-\n')
          goodTheta(l)=theta(i,j,k);
          goodThetai(l)=i;
          goodThetaj(l)=j;
          goodThetak(l)=k;
          goodAlpha(l)=alphas(i);
          goodBeta(l)=betas(j);
          goodGamma(l)=gammas(k);
          %! some proper weighting factors are needed
          meritFactor(l)=1000/(  (1e-5+abs(theta(i,j,k)))  *    (1+wrapTo360(oldAlphaOpt-alphas(i)))  *(1+wrapTo360(oldBetaOpt-betas(j))) * (1+wrapTo360(oldGammaOpt-gammas(k))) );
          l=l+1;
        end
        %input('xxx')
      end
    end
  end
  if not(exist('meritFactor'))
    alphaOpt=oldAlphaOpt;
    betaOpt=oldBetaOpt;
    gammaOpt=oldGammaOpt;
    ii=find(alphas==oldAlphaOpt);
    jj=find(betas==oldBetaOpt);
    kk=find(gammas==oldGammaOpt);
    ii=1;jj=1;kk=1;
    forceVector=squeeze(totalForceVector(:,ii,jj,kk));
    %controlVector/norm(controlVector)
    %forceVector/norm(forceVector)
    %fprintf('\n findBestAerodynamicAngles: problem with theta\n')
    %input('')
  elseif oldAlphaOpt==0 && oldBetaOpt==0 && oldGammaOpt==0
    [~,optIndex]=min(goodTheta);
    alphaOpt=goodAlpha(optIndex);
    betaOpt=goodBeta(optIndex);
    gammaOpt=goodGamma(optIndex);
    forceVector=squeeze(totalForceVector(:,goodThetai(optIndex),goodThetaj(optIndex),goodThetak(optIndex)));
  else
    [~,optIndex]=max(meritFactor);
    alphaOpt=goodAlpha(optIndex);
    betaOpt=goodBeta(optIndex);
    gammaOpt=goodGamma(optIndex);
    forceVector=squeeze(totalForceVector(:,goodThetai(optIndex),goodThetaj(optIndex),goodThetak(optIndex)));
  end  

    %magcv=sqrt(controlvector(1)^2+controlvector(2)^2+controlvector(3)^2);
    %fprintf('\n %1.1e %1.3f %1.3f %1.3f',magcv,controlvector(1)/magcv,controlvector(2)/magcv,controlvector(3)/magcv);
    %magfv=sqrt(forcevector(1)^2+forcevector(2)^2+forcevector(3)^2);
    %fprintf('\n %1.1e %1.3f %1.3f %1.3f',magfv,forcevector(1)/magfv,forcevector(2)/magfv,forcevector(3)/magfv);
    %fprintf('-------')

%{
    figure
      histogram(meritfactor)
    figure
      plot(meritfactor)
    figure
      subplot(2,1,1)
        histogram(theta)
      subplot(2,1,2)
        histogram(phi)
    figure
      for n=1:l-1
        vectarrow([0 0 0],totalforcevector(:,squeeze(goodthetai(n)),squeeze(goodthetaj(n)),squeeze(goodThetak(n))));
        hold on;axis equal;
      end
      hold off;
    Iz=[0 0 1]';
    figure
      for n=1:l-1
        RzY =[cosd(goodgamma(n)) -sind(goodgamma(n)) 0; sind(goodgamma(n)) cosd(goodgamma(n)) 0; 0 0 1]; %%roll
        Ry =[cosd(goodbeta(n))  0 sind(goodbeta(n))  ; 0 1 0                          ; -sind(goodbeta(n)) 0 cosd(goodbeta(n))]; %% pitch
        RzR=[cosd(goodalpha(n)) -sind(goodalpha(n)) 0; sind(goodalpha(n)) cosd(goodalpha(n)) 0; 0 0 1]; %% yaw
        Ig=RzY*Ry*RzR*Iz;
        vectarrow([0 0 0],Ig);
        hold on;axis equal;
      end
    hold off;
    
%}
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theta,phi]=thetaphi(refvec, vec)
  theta = atan2d(norm(cross(refvec,vec)), dot(refvec,vec));
  phi=atand( (refvec(3)-vec(3)) / (refvec(2)-vec(2)) );
end



