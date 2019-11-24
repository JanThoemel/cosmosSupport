function [CD,CL]=aeroDragLiftSentman(theta,Tatmos,v,rho)

  GSImodel=1;
  
  if GSImodel==0 %% simple GSI model
    CD=-abs(2.8*sind(theta-90));      %% simplified formula, adopted from Traub
    CL=-abs(0.42*sind(2*(theta-90))); %% simplified formula, adopted from Traub
  elseif GSImodel==1 %% Sentman GSI model
    %% currently, only atomic oxygen is considered, this needs to be generalized to account for Helium too.
    R=8.31446261815324;           %% ideal gas constant 	[J·K-1·mol-1]
    M=0.016;                      %% molar mass oxygen [kg/mol]
    no=rho/M*6.022140857e23;      %% number density oxygen [kg/m^3/kg*mol*number/mol]=number/m^3
    %% Dornboos
    cmpj=sqrt(2*R/M*Tatmos);
    s=v/cmpj; %% 

    %% Sentman, L.H.: Free molecule flow theory and its application
    %% to the determination of aerodynamic forces. Technical Report,
    %% Lockheed Aircraft Corporation (1961)
    gammaCHEM=cosd(theta);
    L=sind(theta);
    G=0.5/s^2;
    P=1/s*exp(-gammaCHEM^2*s^2);
    Q=1+G;
    Z=1+erf(gammaCHEM*s);
    %% Piliniski
    %%Pilinski, M.D., Argrow, B.M., Palo, S.E.: Semiempirical model
    %% for satellite energy-
    alphaCHEM=7.5e-17*no*Tatmos/(1+7.5e-17*no*Tatmos);
    %% Koppenwallner, G.: Energy accommodation coefficient and
    %% momentum transfer modeling. In: HTG-TN-08-11. HTG Hyperschall
    %% Technologie Göttingen (2009)
    vrvinf=sqrt(1/2*(1+alphaCHEM*(4*R/M*Tatmos/v^2-1)));

    CD=-(P/sqrt(pi)+gammaCHEM*Q*Z+gammaCHEM/2*vrvinf*(gammaCHEM*sqrt(pi)*Z+P));
    CL=-(L*G*Z+L/2*vrvinf*(gammaCHEM*sqrt(pi)*Z+P));
    if theta==90 
      CD=0;
      CL=0;
    end
    %CD=-abs(2.8*sind(theta-90))      %% simplified formula, adopted from Traub
    %CL=-abs(0.42*sind(2*(theta-90))) %% simplified formula, adopted from Traub
%input('ss' )
  elseif GSImodel==2 %% fast Sentman, to be implemented
    CD=0;
    CL=0;
  else
    fprintf('\n error aeroDragLiftSentman: no GSI model chosen \n');
  end
end
