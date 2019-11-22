function plotting(angles,sst,refPosChange,time,ns,meanMotion,u,e)

  if 1
    figure
      for i=1:3
        subplot(3,1,i)
        for j=1:ns
          plot(time/2/pi*meanMotion,squeeze(e(i,j,:)));hold on;
        end
        legend;
        title(strcat('error coordinate',num2str(i)))
      end
    hold off;
  end
 
  if 1 %% reference position change
    figure
      plot(time,refPosChange(1,:),time,refPosChange(2,:),time,refPosChange(3,:));
      legend('x','y','z');title('reference position change')
    hold off;  
  end
  
  if 1 %% plot general variables
    figure
     subplot(4,3,1)%% roll
       for i=1:ns
         plot(squeeze(time/2/pi*meanMotion),squeeze(angles(1,i,:)));hold on
         names(i)=[{strcat('sat',int2str(i))}];
       end
      ylabel('roll angle [deg]');xlabel('no. of orbits');grid on;hold off;legend(names);
      subplot(4,3,2)%% pitch
      for i=1:ns
         plot(squeeze(time/2/pi*meanMotion),squeeze(angles(2,i,:)));hold on
      end
      ylabel('pitch angle [deg]');xlabel('no. of orbits');grid on;hold off;
      subplot(4,3,3)%%yaw
      for i=1:ns
         plot(squeeze(time/2/pi*meanMotion),squeeze(angles(3,i,:)));hold on
      end
      ylabel('yaw angle [deg]');xlabel('no. of orbits');grid on;hold off
      subplot(4,3,4)%%x
      for i=1:ns
         plot(squeeze(time/2/pi*meanMotion),squeeze(sst(1,i,:)));hold on
      end
      ylabel('x [m]');xlabel('no. of orbits');grid on;hold off
      subplot(4,3,5)%%y
      for i=1:ns
         plot(squeeze(time/2/pi*meanMotion),squeeze(sst(2,i,:)));hold on
      end
      ylabel('y [m]');xlabel('no. of orbits');grid on;hold off
      subplot(4,3,6)%%z
      for i=1:ns
         plot(squeeze(time/2/pi*meanMotion),squeeze(sst(3,i,:)));hold on
      end
      ylabel('z [m]');xlabel('no. of orbits');grid on;hold off
      subplot(4,3,7)%%u1
      for i=1:ns
         plot(squeeze(time/2/pi*meanMotion),squeeze(u(1,i,:)));hold on
      end
      ylabel('u1');xlabel('no. of orbits');grid on;hold off
      subplot(4,3,8)%%u2
      for i=1:ns
         plot(squeeze(time/2/pi*meanMotion),squeeze(u(2,i,:)));hold on
      end
      ylabel('u2');xlabel('no. of orbits');grid on;hold off
      subplot(4,3,9)%%u3
      for i=1:ns
         plot(squeeze(time/2/pi*meanMotion),squeeze(u(3,i,:)));hold on
      end
      ylabel('u3');xlabel('no. of orbits');grid on;hold off  
      subplot(4,3,10)%%u
      for i=1:ns
         plot(squeeze(time/2/pi*meanMotion),squeeze(sst(4,i,:)));hold on
      end
      ylabel('u [m/s]');xlabel('no. of orbits');grid on;hold off
      subplot(4,3,11)%%v
      for i=1:ns
         plot(squeeze(time/2/pi*meanMotion),squeeze(sst(5,i,:)));hold on
      end
      ylabel('v [m/s]');xlabel('no. of orbits');grid on;hold off
      subplot(4,3,12)%%w
      for i=1:ns
         plot(squeeze(time/2/pi*meanMotion),squeeze(sst(6,i,:)));hold on
      end
      ylabel('w [m/s]');xlabel('no. of orbits');grid on;hold off

    %traub = csvread('TraubFig5.csv');
    %csvwrite('zxyplane.csv',squeeze(sst(1:3,2,:)))

    %xzplane=csvread('zxplane.csv');
    %xzyplane=csvread('zxyplane.csv');   
    % set(0, 'DefaultFigureRenderer', 'painters');
  end
  
  if 1 %% 3d printing
    figure
      fontSize=40;
      lineWidth=1;
      for i=1:ns
        plot3(squeeze(sst(1,i,:)),squeeze(sst(2,i,:)),squeeze(sst(3,i,:)),'-','LineWidth',lineWidth);hold on;
        names(i)=[{strcat('sat',int2str(i))}];
        %plot3(xzyplane(1,:),xzyplane(2,:),xzyplane(3,:),squeeze(sst(1,i,:)),squeeze(sst(2,i,:)),squeeze(sst(3,i,:)),'LineWidth',lineWidth);hold on
      end
      xlabel('X [m]','FontSize', fontSize);ylabel('Y [m]','FontSize', fontSize);zlabel('Z [m]','FontSize', fontSize);
      %legend('this research','FontSize', fontSize);  
      legend(names,'FontSize', fontSize,'units','centimeters','Position',[19 15.5 1 1])
      %legend('w/o solar pressure','w/ solar pressure','FontSize', fontSize,'units','centimeters','Position',[15 17 1 1]);
      grid on;axis equal;
      set(gcf,'units','centimeters','position',[0,0,27,20]);
      set(gca,'units','centimeters','position',[5.5 4 16 15]);
      set(gca,'FontSize',fontSize);
      %axis([-2300 150 -320 320 -320 320]);
      %xticks([-2000 -1000 0 ]);yticks([-200 200]);zticks([-200 200]);
      %view([1 -0.5 0.5]);
      %set(gcf, 'InvertHardCopy', 'off');
      %title('3d');
      hold off;

      %savefig('sst3d.fig');
      %print('-painters','-dmeta','sst3d.emf')

  end
    if 0 %% 2d printing
      figure
        plot(traub(:,1),traub(:,2),xzplane(1,:),xzplane(3,:),squeeze(sst(1,i,:)),squeeze(sst(3,i,:)))
        xlabel('X [m]','FontSize', fontSize);ylabel('Y [m]','FontSize', fontSize);
        grid on;axis equal;hold off;
        axis([-1400 100 -400 600]);
        box off;
        legend('Traub et al. [4] ','this research in-plane','this research in/out-of-plane','FontSize', fontSize,'Units','centimeters','Position',[16 16 1 1]);
        
        %set(gca,'FontSize',fontSize);
        %set(gcf,'units','centimeters','position',[0,0,25,20])
        %legend('w/o solar pressure','w/ solar pressure','FontSize', fontSize,'Units','centimeters','Position',[12 17 1 1]);
        set(gcf,'units','centimeters','position',[0,0,27,20]);
        set(gca,'units','centimeters','position',[5.5 3 20 15]);
        set(gca,'FontSize',fontSize);
        print('-painters','-dmeta','sst2d.emf')
    end
end