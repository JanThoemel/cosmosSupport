function [P,IR,A,B]=riccatiequation(meanMotion)
  %R=diag([1e-13 1e-14 1e-14]); %% this is R given in Ivanov's IAC paper. It seems to be wrong
  %R=diag([1e12 1e14 1e14]);    %% 
  %R=diag([1e6 1e6 1e6]);    %% 550, -2900
  %R=diag([1e9 1e9 1e9]);    %% 550, -2900
  %R=diag([1e10 1e10 1e10]);    %% 550, -2900
  %R=diag([1e11 1e11 1e11]);    %% 550, -2900
  %R=diag([1e12 1e12 1e12]);    %% This is an R of Traub %550, -2900 
  %R=diag([1e13 1e14 1e14]);    %% this is my R assuming Ivanov made a sign error (reported IR/R^-1 instead of R)
  %R=diag([1e13 1e13 1e13]);     %% 550, -2900 
  R=diag([1e14 1e14 1e14]);    %% 550, -2900
  %R=diag([1e15 1e15 1e15]);    %% 550, -2900
  %R=diag([1e16 1e16 1e16]);    %% 550, -2900
  %R=diag([1e17 1e17 1e17]);    %% 550, -2900
  %R=diag([1e18 1e18 1e18]);    %% 550, -2900
  %R=diag([1e19 1e19 1e19]);    %% 550, -2900
  %R=diag([1e20 1e20 1e20]);    %% 550, -2900
  %R=diag([1e22 1e22 1e22]);    %% 
  Q=eye(6);
  E=eye(3);
  Z=zeros(3,3);
  C=[0 0 0; 0 -meanMotion^2 0;0 0 3*meanMotion^2];
  D=[0 0 -2*meanMotion;0 0 0;2*meanMotion 0 0];

  A=[Z E; C D];
  B=[Z;E];
  IR=inv(R);
  %%https://nl.mathworks.com/help/control/ref/care.html
  %% the function care is replaced by icare in later matlab versions
  S=zeros(6,3);
  E2=eye(6);
 [P,~,~] = care(A,B,Q,R,S,E2);
end