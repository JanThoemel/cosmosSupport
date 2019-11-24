function sunforcevectorz2=solarPressureForce(solarPressure, sunlight,normal, theta,panelSurface,noPanels)
  %% model as per:
  %% 'Modelling of Solar Radiation Pressure Effects: Parameter Analysis for the MICROSCOPE Mission'
  %% by Meike List, Stefanie Bremer, Benny Rievers, and Hanns Selig
  %??,BOL  0.0727 
  %??,EOL  0.05 
  %??,BOL  0.007
  %??,EOL  0.030
    sunlight=-sunlight;
    gammaSunSpecular=0.072;
    gammaSunDiffusive=0.007;
    sunforcevectorz2=-2*solarPressure*((1-gammaSunSpecular)*sunlight+2*(gammaSunSpecular*cosd(theta)+1/3*gammaSunDiffusive)*normal)*cosd(theta)*panelSurface*noPanels;
  %% magnitude pretty  small
  %% direction in y wrong
end
