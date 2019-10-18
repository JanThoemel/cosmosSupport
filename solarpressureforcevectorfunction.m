function solarpressureforcevector = solarpressureforcevectorfunction(sunlight,panelSurface,noxpanels,noypanels,nozpanels,alphas,betas,gammas)
  solarpressureforcevector=zeros(3,size(alphas,2),size(betas,2),size(gammas,2));
  if norm(sunlight)==0
    return
  end  
  rotspeed=30;draw=0;
  solarPressure=norm(sunlight);
  sunlight     =sunlight/solarPressure;
  fprintf('computing sunlightdynamics');  
  %% must be in dimensions of force, i.e. N
  axislength=1.1*solarPressure*panelSurface;
  Ry90=roty(90);
  %Rz90=rotz(90);
  Rx90=rotx(90);
  Iz = [0 0 0.66*axislength]';
  Ix=Ry90*Iz;
  Iy=Rx90*Iz;
  %% zpanel
  pz1  = [axislength*0.9,axislength*0.9,0];pz2 = [axislength*0.9,-axislength*0.9,0];pz3 = [-axislength*0.9,-axislength*0.9,0];pz4 = [-axislength*0.9,axislength*0.9,0];
  pz12 = 0.33*pz1; pz22 = 0.33*pz2; pz32 = 0.33*pz3; pz42 = 0.33*pz4;
  pz13 = 0.66*pz1; pz23 = 0.66*pz2; pz33 = 0.66*pz3; pz43 = 0.66*pz4;
  thetasun=zeros(size(gammas,2),size(betas,2),size(alphas,2));
  phisun=zeros(size(gammas,2),size(betas,2),size(alphas,2));
  sundragcoef=zeros(size(gammas,2),size(betas,2),size(alphas,2));
  sunliftcoef=zeros(size(gammas,2),size(betas,2),size(alphas,2));
  sunforcevectorz=[0 0 0]';
  sunforcevectorx=[0 0 0]';
  sunforcevectory=[0 0 0]';
  for k=1:size(gammas,2) %% yaw
    for j=1:size(betas,2) %% pitch
      for i=1:size(alphas,2) %% roll
        %% rotation matrizes
        RzY =[cosd(gammas(k)) -sind(gammas(k)) 0; sind(gammas(k)) cosd(gammas(k)) 0; 0 0 1]; %%yaw
        Ry =[cosd(betas(j))  0 sind(betas(j))  ; 0 1 0                          ; -sind(betas(j)) 0 cosd(betas(j))]; %% pitch
        RzR=[cosd(alphas(i)) -sind(alphas(i)) 0; sind(alphas(i)) cosd(alphas(i)) 0; 0 0 1]; %% roll
        if nozpanels %% zpanel
          Igz=RzY*Ry*RzR*Iz;
          [thetasun(i,j,k),phisun(i,j,k),Igz2]=thetaphi1(sunlight,Igz);
          [sundragcoef(i,j,k), sunliftcoef(i,j,k) ]=sundragliftcoef(thetasun(i,j,k));                                                               
          sunforcevectorz=-sunlight   * sundragcoef(i,j,k)*solarPressure*panelSurface;
          ax=cross(sunlight,Igz2) ;
          if norm(ax)~=0
            liftvector = rodrigues_rot(sunlight,ax,90/180*pi);
          else
            liftvector = [0 0 0]';
          end            
          sunforcevectorz=nozpanels*(-liftvector  *sunliftcoef(i,j,k)*solarPressure*panelSurface  +  sunforcevectorz);
        end
        if noxpanels %% xpanel
          Igx=RzY*Ry*RzR*Ix;
          [thetasun(i,j,k),phisun(i,j,k),Igx2]=thetaphi1(sunlight,Igx);
          [sundragcoef(i,j,k),sunliftcoef(i,j,k)]=sundragliftcoef(thetasun(i,j,k));
          sunforcevectorx=-sunlight           *  sundragcoef(i,j,k)*solarPressure*panelSurface;
          ax=cross(sunlight,Igx2);
          if norm(ax)~=0
            liftvector = rodrigues_rot(sunlight,ax,90/180*pi);
          else
            liftvector = [0 0 0]';
          end
          sunforcevectorx=noxpanels*(-liftvector         *  sunliftcoef(i,j,k)*solarPressure*panelSurface  +  sunforcevectorx);
        end
        if noypanels %% ypanel
          Igy=RzY*Ry*RzR*Iy;
          [thetasun(i,j,k),phisun(i,j,k),Igy2]=thetaphi1(sunlight,Igy);
          [sundragcoef(i,j,k),sunliftcoef(i,j,k)]=sundragliftcoef(thetasun(i,j,k));                                       
          sunforcevectory=-sunlight     *   sundragcoef(i,j,k)*solarPressure*panelSurface;
          ax=cross(sunlight,Igy2);
          if norm(ax)~=0
            liftvector = rodrigues_rot(sunlight,ax,90/180*pi);
          else
            liftvector = [0 0 0]';
          end                        
          sunforcevectory=noypanels*(-liftvector   *    sunliftcoef(i,j,k)*solarPressure*panelSurface+sunforcevectory);
        end
        %%draw
        if draw
            vectarrow([0 0 0],sunlight*solarPressure*panelSurface);hold on;text(sunlight(1)*solarPressure*panelSurface,sunlight(2)*solarPressure*panelSurface,sunlight(3)*solarPressure*panelSurface,"sunlight",'HorizontalAlignment','left','FontSize',6);
            pg = [(RzY*Ry*RzR*pz1')' ; (RzY*Ry*RzR*pz2')' ; (RzY*Ry*RzR*pz3')' ; (RzY*Ry*RzR*pz4')' ; (RzY*Ry*RzR*pz1')'];
            pg2 = [(RzY*Ry*RzR*pz12')' ; (RzY*Ry*RzR*pz22')' ; (RzY*Ry*RzR*pz32')' ; (RzY*Ry*RzR*pz42')' ; (RzY*Ry*RzR*pz12')'];
            pg3 = [(RzY*Ry*RzR*pz13')' ; (RzY*Ry*RzR*pz23')' ; (RzY*Ry*RzR*pz33')' ; (RzY*Ry*RzR*pz43')' ; (RzY*Ry*RzR*pz13')'];
            if nozpanels
              vectarrow([0 0 0],sunforcevectorz);hold on;text(sunforcevectorz(1),sunforcevectorz(2),sunforcevectorz(3),"sunforce",'HorizontalAlignment','left','FontSize',6);
              vectarrow([0 0 0],Igz);hold on;text(Igz(1),Igz(2),Igz(3),"normal",'HorizontalAlignment','left','FontSize',6);
              line(pg(:,1), pg(:,2), pg(:,3));line(pg2(:,1), pg2(:,2), pg2(:,3));line(pg3(:,1), pg3(:,2), pg3(:,3));hold on;
              axis([-axislength axislength -axislength axislength -axislength axislength]);
            end
            if noxpanels
              vectarrow([0 0 0],Igx);hold on;text(Igx(1),Igx(2),Igx(3),"normalx",'HorizontalAlignment','left','FontSize',6);hold on;
              vectarrow([0 0 0],sunforcevectorx);hold on;text(sunforcevectorx(1),sunforcevectorx(2),sunforcevectorx(3),"sunforcex",'HorizontalAlignment','left','FontSize',6);
              line(pgx(:,1), pgx(:,2), pgx(:,3));line(pgx2(:,1), pgx2(:,2), pgx2(:,3));line(pgx3(:,1), pgx3(:,2), pgx3(:,3));hold on;
              axis([-axislength axislength -axislength axislength -axislength axislength]);
            end
            if noypanels
              vectarrow([0 0 0],Igy);hold on;text(Igy(1),Igy(2),Igy(3),"normaly",'HorizontalAlignment','left','FontSize',6);
              vectarrow([0 0 0],sunforcevectory);hold on;text(sunforcevectory(1),sunforcevectory(2),sunforcevectory(3),"sunforcey",'HorizontalAlignment','left','FontSize',6); 
              line(pgy(:,1), pgy(:,2), pgy(:,3));line(pgy2(:,1), pgy2(:,2), pgy2(:,3));line(pgy3(:,1), pgy3(:,2), pgy3(:,3));hold on;
              axis([-axislength axislength -axislength axislength -axislength axislength]);
            end
          %text(0,0,0,strcat('scalingfactor: ',int2str(scalingfactor)),'HorizontalAlignment','left','FontSize',6);
          xlabel('fx [N]');ylabel('fy [N]');zlabel('fz [N]');
          axis([-axislength axislength -axislength axislength -axislength axislength]);
          hold off;
          pause(1/rotspeed);
        end
        solarpressureforcevector(:,i,j,k)=sunforcevectorz+sunforcevectorx+sunforcevectory;
      end
    end
  end
  fprintf(' - done\n');
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                       
function [theta,phi,Ig2]=thetaphi1(refvec, vec)
%% this function computes theta,phi, Ig2
  theta = atan2d(norm(cross(refvec,vec)), dot(refvec,vec));                                                                                                       
  if theta>90                                                                                                                                                      
    theta=180-theta;                                                                                                                                             
    Ig2=-vec;                                                                                                                                                     
  else                                                                                                                                                             
    Ig2=vec;                                                                                                                                                     
  end                                                                                                                                                              
  phi=atand( (refvec(3)-vec(3)) / (refvec(2)-vec(2)+1e-30) );  
end     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [drag,lift]=sundragliftcoef(theta)
  drag=-abs(sind(theta-90));      %% simplified formula
  lift=-abs(sind(2*(theta-90)));  %% simplified formula     
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%