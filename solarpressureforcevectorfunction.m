function  solarpressureforcevector = solarpressureforcevectorfunction(sunlight,panelSurface,...
          noxpanels,noypanels,nozpanels,alphas,betas,gammas)
%%  computes the experienced force caused by solar pressure for all satellite attitudes (pitch, yaw, roll)
%%  input variables
%%  -sunlight
%%  -panelSurface
%%  -noxpanels,noypanels,nozpanels
%%  -alphas,betas,gammas
%%  output variables
%%  -solarpressureforcevector  array     [F]
%%  local variables
%%  -Iz,Ix,Iy = panel normal vectors array [-] 
%%  remarks:
%%  -the vector sunlight is nondimensionalized, which seems not necessary
%%  -there is visualization possible, which should not be necessary anymore

  solarpressureforcevector=zeros(3,size(alphas,2),size(betas,2),size(gammas,2));

  if norm(sunlight)==0
    return
  end  
  rotSpeed=30;draw=0;
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
  sunforcevectorZ=[0 0 0]';
  sunforcevectorX=[0 0 0]';
  sunforcevectorY=[0 0 0]';
  for k=1:size(gammas,2) %% yaw
    for j=1:size(betas,2) %% pitch
      for i=1:size(alphas,2) %% roll
        %% rotation matrizes
        RzY=[cosd(gammas(k)) -sind(gammas(k)) 0; sind(gammas(k)) cosd(gammas(k)) 0; 0 0 1];        %%yaw
        Ry =[cosd(betas(j))  0 sind(betas(j))  ; 0 1 0                            ; -sind(betas(j)) 0 cosd(betas(j))]; %% pitch
        RzR=[cosd(alphas(i)) -sind(alphas(i)) 0; sind(alphas(i)) cosd(alphas(i)) 0; 0 0 1];        %% roll
        if nozpanels %% zpanel
          normalZ=RzY*Ry*RzR*Iz; %% determine the normal vector of a z panel for yaw,pitch, roll
          [thetasun(i,j,k),phisun(i,j,k),normalZ2]=thetaphi1(sunlight,normalZ);       %% determine theta and phi of the normal wrt sunlight vector
          %[sundragcoef(i,j,k), sunliftcoef(i,j,k) ]=sunDragLiftCoef(thetasun(i,j,k)); %% determine CD, CL                                                              
          %sunforcevectorZ=-sunlight   * sundragcoef(i,j,k)*solarPressure*panelSurface;%% determine drag vector
          %ax=cross(sunlight,normalZ2) ;
          %if norm(ax)~=0
          %  liftvector = rodrigues_rot(sunlight,ax,90/180*pi);
          %else
          %  liftvector = [0 0 0]';
          %end            
          %%% determine lift vector
          %sunforcevectorZ=nozpanels*(-liftvector  *sunliftcoef(i,j,k)*solarPressure*panelSurface  +  sunforcevectorZ)
          sunforcevectorZ=solarPressureForce(solarPressure, sunlight,normalZ/norm(normalZ), thetasun(i,j,k),panelSurface,nozpanels);
        end
        if noxpanels %% xpanel
          normalX=RzY*Ry*RzR*Ix;
          [thetasun(i,j,k),phisun(i,j,k),normalX2]=thetaphi1(sunlight,normalX);
          %[sundragcoef(i,j,k),sunliftcoef(i,j,k)]=sunDragLiftCoef(thetasun(i,j,k));
          %sunforcevectorX=-sunlight           *  sundragcoef(i,j,k)*solarPressure*panelSurface;
          %ax=cross(sunlight,normalX2);
          %if norm(ax)~=0
          %  liftvector = rodrigues_rot(sunlight,ax,90/180*pi);
          %else
          %  liftvector = [0 0 0]';
          %end
          %sunforcevectorX=noxpanels*(-liftvector         *  sunliftcoef(i,j,k)*solarPressure*panelSurface  +  sunforcevectorX);
          sunforcevectorX=solarPressureForce(solarPressure, sunlight,normalX/norm(normalX), thetasun(i,j,k),panelSurface,noxpanels);
        end
        if noypanels %% ypanel
          normalY=RzY*Ry*RzR*Iy;
          %[thetasun(i,j,k),phisun(i,j,k),normalY2]=thetaphi1(sunlight,normalY);
          %[sundragcoef(i,j,k),sunliftcoef(i,j,k)]=sunDragLiftCoef(thetasun(i,j,k));                                       
          %sunforcevectorY=-sunlight     *   sundragcoef(i,j,k)*solarPressure*panelSurface;
          %ax=cross(sunlight,normalY2);
          %if norm(ax)~=0
          %  liftvector = rodrigues_rot(sunlight,ax,90/180*pi);
          %else
          %  liftvector = [0 0 0]';
          %end                        
          %sunforcevectorY=noypanels*(-liftvector   *    sunliftcoef(i,j,k)*solarPressure*panelSurface+sunforcevectorY);
          sunforcevectorY=solarPressureForce(solarPressure, sunlight,normalY/norm(normalY), thetasun(i,j,k),panelSurface,noypanels);
        end
        %%draw
        if draw
          vectarrow([0 0 0],sunlight*solarPressure*panelSurface);hold on;text(sunlight(1)*solarPressure*panelSurface,sunlight(2)*solarPressure*panelSurface,sunlight(3)*solarPressure*panelSurface,"sunlight",'HorizontalAlignment','left','FontSize',6);
          pg = [(RzY*Ry*RzR*pz1')' ; (RzY*Ry*RzR*pz2')' ; (RzY*Ry*RzR*pz3')' ; (RzY*Ry*RzR*pz4')' ; (RzY*Ry*RzR*pz1')'];
          pg2 = [(RzY*Ry*RzR*pz12')' ; (RzY*Ry*RzR*pz22')' ; (RzY*Ry*RzR*pz32')' ; (RzY*Ry*RzR*pz42')' ; (RzY*Ry*RzR*pz12')'];
          pg3 = [(RzY*Ry*RzR*pz13')' ; (RzY*Ry*RzR*pz23')' ; (RzY*Ry*RzR*pz33')' ; (RzY*Ry*RzR*pz43')' ; (RzY*Ry*RzR*pz13')'];
          if nozpanels
            vectarrow([0 0 0],sunforcevectorZ);hold on;text(sunforcevectorZ(1),sunforcevectorZ(2),sunforcevectorZ(3),"sunforce",'HorizontalAlignment','left','FontSize',6);
            vectarrow([0 0 0],normalZ);hold on;text(normalZ(1),normalZ(2),normalZ(3),"normal",'HorizontalAlignment','left','FontSize',6);
            line(pg(:,1), pg(:,2), pg(:,3));line(pg2(:,1), pg2(:,2), pg2(:,3));line(pg3(:,1), pg3(:,2), pg3(:,3));hold on;
            axis([-axislength axislength -axislength axislength -axislength axislength]);
          end
          if noxpanels
            vectarrow([0 0 0],normalX);hold on;text(normalX(1),normalX(2),normalX(3),"normalx",'HorizontalAlignment','left','FontSize',6);hold on;
            vectarrow([0 0 0],sunforcevectorX);hold on;text(sunforcevectorX(1),sunforcevectorX(2),sunforcevectorX(3),"sunforcex",'HorizontalAlignment','left','FontSize',6);
            line(pgx(:,1), pgx(:,2), pgx(:,3));line(pgx2(:,1), pgx2(:,2), pgx2(:,3));line(pgx3(:,1), pgx3(:,2), pgx3(:,3));hold on;
            axis([-axislength axislength -axislength axislength -axislength axislength]);
          end
          if noypanels
            vectarrow([0 0 0],normalY);hold on;text(normalY(1),normalY(2),normalY(3),"normaly",'HorizontalAlignment','left','FontSize',6);
            vectarrow([0 0 0],sunforcevectorY);hold on;text(sunforcevectorY(1),sunforcevectorY(2),sunforcevectorY(3),"sunforcey",'HorizontalAlignment','left','FontSize',6); 
            line(pgy(:,1), pgy(:,2), pgy(:,3));line(pgy2(:,1), pgy2(:,2), pgy2(:,3));line(pgy3(:,1), pgy3(:,2), pgy3(:,3));hold on;
            axis([-axislength axislength -axislength axislength -axislength axislength]);
          end
          %text(0,0,0,strcat('scalingfactor: ',int2str(scalingfactor)),'HorizontalAlignment','left','FontSize',6);
          xlabel('fx [N]');ylabel('fy [N]');zlabel('fz [N]');
          axis([-axislength axislength -axislength axislength -axislength axislength]);
          hold off;
          pause(1/rotSpeed);
        end
        solarpressureforcevector(:,i,j,k)=sunforcevectorZ+sunforcevectorX+sunforcevectorY;
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
function [drag,lift]=sunDragLiftCoef(theta)
  drag=-abs(2.8*sind(theta-90));      %% simplified formula
  lift=-abs(0.42*sind(2*(theta-90)));  %% simplified formula     
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

