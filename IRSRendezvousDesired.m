function [sstDesired]=IRSRendezvousDesired(timetemptemp,~,i,goFoFli)
%% desired IRS rendezvous solution
  startSecondPhase=2*90*60;             %% [s]
  sstDesired=zeros(9,1,size(timetemptemp,2));

  if timetemptemp<startSecondPhase
    sstDesired(1,1,:)=-930.46;          %% x for rendezvous
    %sstDesired(2,2,1)=55.27;           %% y for rendezvous
    %sstDesired(3,2,1)=82.5;            %% z for rendezvous
    %sstDesired(4,2,1)=-0.04;           %% u for rendezvous
    %sstDesired(5,2,1)=0.29;            %% v for rendezvous
    %sstDesired(6,2,1)=-0.17;           %% w for rendezvous
  else
    sstDesired(1,1,:)=0.1;
  end
 % sstDesired
end