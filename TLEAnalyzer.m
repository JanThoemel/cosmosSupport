clear all;clc;close all;
%% location of the TLE files, they shall be zipped text files in twolines, i.e. with a name line
%% if you use the CelesTrak TLE retriever, then this setting will work
tleFiles=dir(strcat(getenv('USERPROFILE'),'\Documents\My TLEs\2019\Full Catalog\*'));
%% which satellites to be processed, use NORAD ID
catalogueID=[43765,43794,43799]; catalogueNames=["HAWK-A" "HAWK-B" "HAWK-C"];
%catalogueID=[43197,43196]; catalogueNames=["GOMX-4A" "GOMX-4B"]; 

%% initialize arrays
epochTime=zeros( size(catalogueID,2) ,1);
RAAN=zeros( size(catalogueID,2) ,1);
sma=zeros( size(catalogueID,2) ,1);
arOfPeri=zeros( size(catalogueID,2) ,1);
dn=zeros( size(catalogueID,2) ,1);
inclination=zeros( size(catalogueID,2) ,1);
meanAnomalyFromANAtMidnight=zeros( size(catalogueID,2) ,1);

epochTime2=zeros( size(catalogueID,2) ,1);
RAAN2=zeros( size(catalogueID,2) ,1);
sma2=zeros( size(catalogueID,2) ,1);
arOfPeri2=zeros( size(catalogueID,2) ,1);
dn2=zeros( size(catalogueID,2) ,1);
inclination2=zeros( size(catalogueID,2) ,1);
meanAnomalyFromANAtMidnight2=zeros( size(catalogueID,2) ,1);

for i=3:size(tleFiles,1) %% %% cycle overall all files, first two entries of tleFiles are '.' and '..'. They shall be skipped
  %% copy file here
  copyfile(strcat(tleFiles(i).folder,'\',tleFiles(i).name), 'temp\');
  %% unzip file
  unzippedfilename=unzip(strcat('temp\',tleFiles(i).name),'temp\');
  %% readfile into variable
  [epochTimeNext,RAANNext,smaNext,arOfPeriNext,inclinationNext,meanMotionNext,meanAnomalyNext,foundno]=readTLE( unzippedfilename{1},catalogueID);
  %% add TLE of catalogueID in arrary
  if size(epochTimeNext,2)==1 && size(epochTimeNext,1)==size(catalogueID,2) && sum(foundno)~=0
    dnNext=zeros( size(catalogueID,2) ,1);
    meanAnomalyFromANAtMidnightNext=zeros( size(catalogueID,2) ,1);
    for j=1:size(catalogueID,2)
        if epochTimeNext(j)~=0
          datestringNext=strcat('20',num2str(epochTimeNext(j)));
          dnNext(j)=datenum(str2num(datestringNext(1:4)),1,str2num(datestringNext(5:10)));
          %% compute second of day at epoch  [s]
          dayFraction=(str2num(datestringNext)-round(str2num(datestringNext),0))*24*60*60;
          %% compute traveled distance [deg]
          distanceAtDay=dayFraction*meanMotionNext(j);
          %% compute 
          meanAnomalyAtMidnight=wrapTo360(meanAnomalyNext(j)-distanceAtDay);          
          %% compute 
          meanAnomalyFromANAtMidnightNext(j)=wrapTo360(meanAnomalyAtMidnight-arOfPeriNext(j));
        else
          dnNext(j)=0;
          meanAnomalyFromANAtMidnightNext(j)=0;
        end
    end
    %% add new entries
    epochTime=[epochTime epochTimeNext];
    dn=[dn dnNext];
    RAAN=[RAAN RAANNext];
    sma=[sma smaNext];
    arOfPeri=[arOfPeri arOfPeriNext];
    inclination=[inclination inclinationNext];
    meanAnomalyFromANAtMidnight=[meanAnomalyFromANAtMidnight meanAnomalyFromANAtMidnightNext];
    
    %for j=1:size(catalogueID,2) %% add data already here in the array of structs
    %  sat2(j).epochTime=[sat2(j).epochTime epochTimeN(j)];
    %end
  elseif size(epochTimeNext,2)>1 && size(epochTimeNext,1)==size(catalogueID,2) && sum(foundno)~=0
    dnNext2=zeros(size(epochTimeNext));   
    meanAnomalyFromANAtMidnightNext2=zeros(size(epochTimeNext));   
    for j=1:size(epochTimeNext,2)
      for k=1:size(epochTimeNext,1)
        if epochTimeNext(k,j)~=0
          datestringNext=strcat('20',num2str(epochTimeNext(k,j)));
          dnNext2(k,j)=datenum(str2num(datestringNext(1:4)),1,str2num(datestringNext(5:end)));
          %% compute second of day at epoch  [s]
          dayFraction=(str2num(datestringNext)-round(str2num(datestringNext),0))*24*60*60;
          %% compute traveled distance [deg]
          distanceAtDay=dayFraction*meanMotionNext(k,j);
          %% compute 
          meanAnomalyAtMidnight=wrapTo360(meanAnomalyNext(k,j)-distanceAtDay);          
          %% compute 
          meanAnomalyFromANAtMidnightNext2(k,j)=wrapTo360(meanAnomalyAtMidnight-arOfPeriNext(k,j));
        else
          dnNext2(k,j)=0;
          meanAnomalyFromANAtMidnightNext2(k,j)=0;
        end
      end
    end
    %% add new entries
    epochTime2=[epochTime2 epochTimeNext];
    dn2=[dn2 dnNext2];
    RAAN2=[RAAN2 RAANNext];
    sma2=[sma2 smaNext];
    arOfPeri2=[arOfPeri2 arOfPeriNext];
    inclination2=[inclination2 inclinationNext];
    meanAnomalyFromANAtMidnight2=[meanAnomalyFromANAtMidnight2 meanAnomalyFromANAtMidnightNext2];
  elseif size(epochTimeNext,2)==1 && size(epochTimeNext,1)==size(catalogueID,2) && sum(foundno)~=1
    %fprintf('\n error 3000');
    %input('');
    ;
  else
    fprintf('\n error 4000');
    input('');
  end
  %% delete zip and unzipped file
  delete(strcat('temp\',tleFiles(i).name));
  delete(unzippedfilename{1});  
  
  if i>3
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
  end
  fprintf('number of loaded TLEs %4.0f/%4.0f (%4.0f%%)', i, size(tleFiles,1),double(i)/double(size(tleFiles,1))*100 )
  pause(0.001);
  
end

%% remove first column (the one that was created when creating the variables)
epochTime(:,1)=[];
RAAN(:,1)=[];
sma(:,1)=[];
arOfPeri(:,1)=[];
dn(:,1)=[];
inclination(:,1)=[];
meanAnomalyFromANAtMidnight(:,1)=[];

epochTime2(:,1)=[];
RAAN2(:,1)=[];
sma2(:,1)=[];
arOfPeri2(:,1)=[];
dn2(:,1)=[];
inclination2(:,1)=[];
meanAnomalyFromANAtMidnight2(:,1)=[];


for i=1:size(catalogueID,2)
  %% patch data together
  sat(i).catalogueID=catalogueID(i);
  sat(i).catalogueNames=catalogueNames(i);
  sat(i).epochTime=epochTime(i,:);
  sat(i).RAAN=RAAN(i,:);
  sat(i).sma=sma(i,:);
  sat(i).arOfPeri=arOfPeri(i,:);
  sat(i).dn=dn(i,:);
  sat(i).inclination=inclination(i,:);
  sat(i).meanAnomalyFromANAtMidnight=meanAnomalyFromANAtMidnight(i,:);
  
  sat(i).epochTime=[sat(i).epochTime epochTime2(i,epochTime2(i,:)~=0)];
  sat(i).RAAN=[sat(i).RAAN RAAN2(i,RAAN2(i,:)~=0)];
  sat(i).sma=[sat(i).sma sma2(i,sma2(i,:)~=0)];
  sat(i).arOfPeri=[sat(i).arOfPeri arOfPeri2(i,arOfPeri2(i,:)~=0)];
  sat(i).dn=[sat(i).dn dn2(i,dn2(i,:)~=0)];
  sat(i).inclination=[sat(i).inclination inclination2(i,inclination2(i,:)~=0)];
  sat(i).meanAnomalyFromANAtMidnight=[sat(i).meanAnomalyFromANAtMidnight meanAnomalyFromANAtMidnight2(i,meanAnomalyFromANAtMidnight2(i,:)~=0)] ;

  %% define time and order arrarys accordingly
  [sat(i).epochTime idx]=sort(sat(i).epochTime);
  sat(i).RAAN=sat(i).RAAN(idx);
  sat(i).sma=sat(i).sma(idx);
  sat(i).arOfPeri=sat(i).arOfPeri(idx);
  sat(i).dn=sat(i).dn(idx);
  sat(i).inclination=sat(i).inclination(idx);
  sat(i).meanAnomalyFromANAtMidnight=sat(i).meanAnomalyFromANAtMidnight(idx);
  
  %% interpolate on first satellite' time instances; this is required to plot relative quantities later
  if i~=1
    [TEMP, idx]                           =unique(sat(i).dn);
    sat(i).RAANInt                        =interp1(TEMP,sat(i).RAAN(idx),sat(1).dn);
    sat(i).smaInt                         =interp1(TEMP,sat(i).sma(idx),sat(1).dn);
    sat(i).arOfPeriInt                    =interp1(TEMP,sat(i).arOfPeri(idx),sat(1).dn);
    sat(i).inclinationInt                 =interp1(TEMP,sat(i).inclination(idx),sat(1).dn);
    sat(i).meanAnomalyFromANAtMidnightInt =interp1(TEMP,sat(i).meanAnomalyFromANAtMidnight(idx),sat(1).dn);
  end  
end %% loop over satellites


%% absolut plot argument of perigee, RAAN, excentrity

figure
set(gcf, 'Position',  [50, 50, 1500, 500]);
  subplot(2,5,1)
    for i=1:size(catalogueID,2)
      plot(sat(i).dn(:)-sat(i).dn(1),sat(i).RAAN);hold on;
      dataNameA(i)=sat(i).catalogueNames;
    end
    ylabel('RAAN [deg]'); xlabel(strcat('time from',{' '},datestr(datetime( sat(1).dn(1),'ConvertFrom','datenum') ),{' '},'[d]') );
    legend(dataNameA);
  subplot(2,5,2)
    for i=1:size(catalogueID,2)
      plot(sat(i).dn(:)-sat(i).dn(1),sat(i).inclination);hold on;
    end
    ylabel('inclination [deg]'); xlabel(strcat('time from',{' '},datestr(datetime( sat(1).dn(1),'ConvertFrom','datenum') ),{' '},'[d]') );
  subplot(2,5,3)
    for i=1:size(catalogueID,2)    
      plot(sat(i).dn(:)-sat(i).dn(1),sat(i).sma-6371000);hold on;
    end
    ylabel('altitude [m]'); xlabel(strcat('time from',{' '},datestr(datetime( sat(1).dn(1),'ConvertFrom','datenum') ),{' '},'[d]') );
    %axis([0 sat(1).dn(end)-sat(1).dn(1) mean(sat(1).sma-6371000)-0.001*std(sat(1).sma-6371000) mean(sat(1).sma-6371000)+0.06*std(sat(1).sma-6371000)])
  subplot(2,5,4)
    for i=1:size(catalogueID,2)  
      plot(sat(i).dn(:)-sat(i).dn(1),sat(i).arOfPeri);hold on;
    end
    ylabel('a o perigee [deg]'); xlabel(strcat('time from',{' '},datestr(datetime( sat(1).dn(1),'ConvertFrom','datenum') ),{' '},'[d]') );
  subplot(2,5,5)
    for i=1:size(catalogueID,2)  
      plot(sat(i).dn(:)-sat(i).dn(1),sat(i).meanAnomalyFromANAtMidnight);hold on;
    end
    ylabel('mean anomaly at 00:00 [deg]'); xlabel(strcat('time from',{' '},datestr(datetime( sat(1).dn(1),'ConvertFrom','datenum') ),{' '},'[d]') );

%% relative plot argument of perigee, RAAN, excentrity (relative to first satellite)
  subplot(2,5,6)
    for i=2:size(catalogueID,2)
      plot(sat(1).dn-sat(1).dn(1),sat(i).RAANInt-sat(1).RAAN);
      hold on;
      dataNameR(i-1)=sat(i).catalogueNames;
    end
    ylabel('rel RAAN [deg]'); xlabel(strcat('time from',{' '},datestr(datetime( sat(1).dn(1),'ConvertFrom','datenum') ),{' '},'[d]') );
    axis([0 sat(1).dn(end)-sat(1).dn(1) -0.001*std(sat(1).RAAN) 0.001*std(sat(1).RAAN)]);
    grid on; legend(dataNameR); 
  subplot(2,5,7)
    for i=2:size(catalogueID,2)
      plot(sat(1).dn-sat(1).dn(1),sat(i).inclinationInt-sat(1).inclination);hold on;
    end
    ylabel('rel inclincation [deg]'); xlabel(strcat('time from',{' '},datestr(datetime( sat(1).dn(1),'ConvertFrom','datenum') ),{' '},'[d]') );
    axis([0 sat(1).dn(end)-sat(1).dn(1) -0.001*std(sat(1).inclination) 0.001*std(sat(1).inclination)]);
    grid on;  
  subplot(2,5,8)
    for i=2:size(catalogueID,2)    
      plot(sat(1).dn-sat(1).dn(1),sat(i).smaInt-sat(1).sma);hold on;
    end
    ylabel('rel altitude [m]'); xlabel(strcat('time from',{' '},datestr(datetime( sat(1).dn(1),'ConvertFrom','datenum') ),{' '},'[d]') );
    %axis([0 sat(1).dn(end)-sat(1).dn(1) -0.002*std(sat(1).sma) 0.002*std(sat(1).sma)]);
    grid on;
  subplot(2,5,9)
    for i=2:size(catalogueID,2)  
      plot(sat(1).dn-sat(1).dn(1),sat(i).arOfPeriInt-sat(1).arOfPeri);hold on;
    end
    ylabel('rel a o perigee [deg]');xlabel(strcat('time from',{' '},datestr(datetime( sat(1).dn(1),'ConvertFrom','datenum') ),{' '},'[d]') );
    axis([0 sat(1).dn(end)-sat(1).dn(1) -0.01*std(sat(1).arOfPeri) 0.01*std(sat(1).arOfPeri)]);
    grid on; 
  subplot(2,5,10)
    for i=2:size(catalogueID,2)  
      plot(sat(1).dn-sat(1).dn(1),sat(i).meanAnomalyFromANAtMidnightInt-sat(1).meanAnomalyFromANAtMidnight);hold on;
    end
    ylabel('rel mean anomaly at 00:00 [deg]');xlabel(strcat('time from',{' '},datestr(datetime( sat(1).dn(1),'ConvertFrom','datenum') ),{' '},'[d]') );
    axis([0 sat(1).dn(end)-sat(1).dn(1) -0.1*std(sat(1).meanAnomalyFromANAtMidnight) 0.1*std(sat(1).meanAnomalyFromANAtMidnight)]);
    grid on; 

function [epochTime,RAAN,sma,arOfPeri,inclination,meanMotion,meanAnomaly,foundno]=readTLE(file, catalog)

% READTLE Read satellite ephemeris data from a NORAD two-line element (TLE) file.
%
% INPUTS:
%   file    - Path to any standard two-line element file.
%   catalog - Optional array of NORAD catalog numbers for the satellites of
%             interest. The default action is to display data from every
%             satellite in the file.
%
% Brett Pantalone
% North Carolina State University
% Department of Electrical and Computer Engineering
% Optical Sensing Laboratory
% mailto:bapantal@ncsu.edu
% http://research.ece.ncsu.edu/osl/
  if nargin < 2
    catalog = [];
  end
  fd = fopen(file,'r');
  if fd < 0, fd = fopen([file '.tle'],'r'); end
  assert(fd > 0,['Can''t open file ' file ' for reading.'])
  n = 0;
  foundno=zeros(size(catalog,2),1);
  
  epochTime=zeros(size(catalog,2),1);
  inclination=zeros(size(catalog,2),1);
  RAAN=zeros(size(catalog,2),1);
  ecc=zeros(size(catalog,2),1);
  arOfPeri=zeros(size(catalog,2),1);
  meanAnomaly=zeros(size(catalog,2),1);
  meanMotion=zeros(size(catalog,2),1);
  T=zeros(size(catalog,2),1);
  sma=zeros(size(catalog,2),1);
  b=zeros(size(catalog,2),1);
  
  %A0 = fgetl(fd);
  A1 = fgetl(fd);
  A2 = fgetl(fd);
  while ischar(A2)
    n = n + 1;
    satnum = str2num(A1(3:7));
    if isempty(catalog) || ismember(satnum, catalog)
      satcount=find(satnum==catalog);
      foundno(satcount)=foundno(satcount)+1;
      %fprintf('%s\n', repmat('-',1,50));
      %fprintf('Satellite: %s\n', A0)
      assert(chksum(A1), 'Checksum failure on line 1')
      assert(chksum(A2), 'Checksum failure on line 2')
      %fprintf('Catalog Number: %d\n', satnum)
      %fprintf('Epoch time: %s\n', A1(19:32)) % YYDDD.DDDDDDDD
      A1(19:32);
      epochTime(satcount,foundno(satcount))=str2double(A1(19:32));
      inclination(satcount,foundno(satcount)) = str2num(A2(9:16));
      %fprintf('inclinationination: %f deg\n', inclination)
      RAAN(satcount,foundno(satcount)) = str2num(A2(18:25));
      %fprintf('RA of ascending node: %f deg\n', RAAN)
      ecc(satcount,foundno(satcount)) = str2num(['.' A2(27:33)]);
      %fprintf('Eccentricity: %f\n', ecc)
      arOfPeri(satcount,foundno(satcount)) = str2num(A2(35:42));
      %fprintf('Arg of perigee: %f deg\n', arOfPeri)
      meanAnomaly(satcount,foundno(satcount)) = str2num(A2(44:51));
      %fprintf('Mean anomaly: %f deg\n', meanAnomaly)
      meanMotion(satcount,foundno(satcount)) = str2num(A2(53:63));
      %fprintf('Mean motion: %f rev/day\n', meanMotion)
      T(satcount,foundno(satcount)) = 86400/meanMotion(satcount,foundno(satcount));
      %fprintf('Period of rev: %.0f s/rev\n', T)
      sma(satcount,foundno(satcount)) = ((T(satcount,foundno(satcount))/(2*pi))^2*398.6e12)^(1/3);
      %fprintf('Semi-major axis: %.0f meters\n', a)
      b(satcount,foundno(satcount)) = sma(satcount,foundno(satcount))*sqrt(1-ecc(satcount,foundno(satcount))^2);
      %fprintf('Semi-minor axis: %.0f meters\n', b)

      %% change of dimenstions for meanMotion
      meanMotion(satcount,foundno(satcount)) = meanMotion(satcount,foundno(satcount))*360/(24*60*60);
      %fprintf('Mean motion: %f deg/second\n', meanMotion)

    
    end
    %A0 = fgetl(fd);
    A1 = fgetl(fd);
    A2 = fgetl(fd);
  end
  fclose(fd);
end
%%
% Checksum (Modulo 10)
% Letters, blanks, periods, plus signs = 0; minus signs = 1
function result = chksum(str)
  result = false; c = 0;
  
  for k = 1:68
    if str(k) > '0' && str(k) <= '9'
      c = c + str(k) - 48;
    elseif str(k) == '-'
      c = c + 1;
    end
  end
  if mod(c,10) == str(69) - 48
    result = true;
  end
  
end
