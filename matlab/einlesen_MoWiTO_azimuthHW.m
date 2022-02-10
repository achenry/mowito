
function []=einlesen_MoWiTO(data,Fs1,Fs2,assign,hwcalibfile,offsetfile,Dichte,temperature,numberhw)

% Kan?le zuweisen
for i = 0:15 
eval(['AI' num2str(i) ' = data(:,' num2str(i+1) ');']);
end
    
U_Mby_3 = AI2;
U_Mby_2 = AI1;
U_Mby_1 = AI0;
U_Welle = AI3;
U_Twr_FA = AI5;
U_Twr_SS = AI6;
U_Torque = AI7;
Prandtl = AI11;
HW1 = AI11;
HW2 = AI13;
HW3 = AI14;


%% Kalibrierungen

load(offsetfile);

off_Mby1 = offsets_env(1,1);
off_Mby2 = offsets_env(2,1);
off_Mby3 = offsets_env(3,1);
off_TwrFA = offsets_env(4,1);
off_TwrSS = offsets_env(5,1);
off_Prandtl = offsets_env(6,1);

% off_Mby1 = -6.465;
% off_Mby2 = 0.972;
% off_Mby3 = -2.875;
% off_TwrFA = 1.746;
% off_TwrSS = -4.084;
% off_Prandtl = 5.998;
load(hwcalibfile);

% Hotwire calibration already done.
% hw_calib_tempcorr;
% tempkorr1 = (215-mean(hw1(:,3)))/ (215 - temperature);
% k1_corr = k1.*tempkorr1;
k1 = mean(k,1);
k1_corr = k1;

p_HW1 = fliplr(k1_corr);
% p_HW2 = fliplr(k2_corr);
% p_HW3 = fliplr(k3_corr);




Torque = 2*U_Torque ;
Mby_1 = (U_Mby_1 - off_Mby1)./ (-0.2448) ;
Mby_2 = (U_Mby_2 - off_Mby2)./ (-0.2722);
Mby_3 = (U_Mby_3 - off_Mby3)./ (-0.2586);
Welle = U_Welle;
Twr_FA = (U_Twr_FA - off_TwrFA)./ (-0.0337);
Twr_SS = (U_Twr_SS - off_TwrSS)./ (-0.0263);
% prandtl_v = sqrt(abs((Prandtl-off_Prandtl).*474.6084.*2./(Dichte)));
hw_v1 = p_HW1(1,1) + p_HW1(1,2).*HW1 + p_HW1(1,3).*HW1.^2 + p_HW1(1,4).*HW1.^3 + p_HW1(1,5).*HW1.^4;
% hw_v2 = p_HW2(1,1) + p_HW2(1,2).*HW2 + p_HW2(1,3).*HW2.^2 + p_HW2(1,4).*HW2.^3 + p_HW2(1,5).*HW2.^4;
% hw_v3 = p_HW3(1,1) + p_HW3(1,2).*HW3 + p_HW3(1,3).*HW3.^2 + p_HW3(1,4).*HW3.^3 + p_HW3(1,5).*HW3.^4;
% hw_v = (hw_v1 +hw_v2 + hw_v3) ./3;
Wind = 6.0;

PowerWind = 0.5 * Dichte * Wind^3 * 3.14 * 0.9^2;
SchubRef = 0.5 * Dichte * Wind^2 * 3.14 * 0.9^2;

 count = data(:,21);
 time = data(:,23);

  for i=1:length(count);
     if  count(i) == 9999
         count(i) = NaN;
         time(i) = NaN;
     end
  end
  
 RPM = (count./time) .*(1000 /1440*60);
 hz=RPM./60;
 omega = 2*pi.*hz;
 countzero = data(:,22); 
 azimuth(1,1) = 0;
 omega2 = resample(nonan(omega),5000,100).*(360/(Fs1*2*pi));
 
 for i = 2: length(data)
    if countzero(i,1) == 1
        if countzero(i-1,1) == 0
            azimuth(i,1) = 0;
        else
            azimuth(i,1) = azimuth(i-1,1)+(omega2(i,1));    
        end
    else 
        azimuth(i,1) = azimuth(i-1,1)+(omega2(i,1));
    end  

 end
    azimuth(:,1) = azimuth(:,1) - 86;
   
    for j = 2: length(data)
        if azimuth(j,1) < 0
            azimuth(j,1)=azimuth(j,1)+360;
        end
    end
    

 PowerWind = 0.5 * Dichte .* hw_v1.^3 * 3.14 * 0.9^2;
 SchubRef = 0.5 * Dichte .* hw_v1.^2 * 3.14 * 0.9^2;

 Power = omega.* smooth(Torque,50);
 Cp =  Power ./ PowerWind;
 Ct =  Twr_FA ./ SchubRef;
%  Cpp = Power ./ PowerWindplus;
%  Ctp = Twr_FA ./ SchubRefplus;
%  Cpm = Power ./ PowerWindminus;
%  Ctm = Twr_FA ./ SchubRefminus;
t1=linspace(0,length(AI1)/Fs1,length(AI1));
t2=linspace(0,length(nonan(omega))/Fs2,length(nonan(omega)));

if assign==1
    assignin('base','omega',omega);
    assignin('base','RPM',RPM);
    assignin('base','Torque',Torque);
    assignin('base','Power',Power);
    assignin('base','t1',t1);
    assignin('base','t2',t2);
    assignin('base','Mby_1',Mby_1);
    assignin('base','Mby_2',Mby_2);
    assignin('base','Mby_3',Mby_3);
    assignin('base','Welle',Welle);
    assignin('base','Twr_FA',Twr_FA);
    assignin('base','Twr_SS',Twr_SS);
    assignin('base','Cp',Cp);
    assignin('base','Ct',Ct);
    assignin('base','azimuth',azimuth);
    assignin('base','p_HW1',p_HW1);
    assignin('base','hw_v1',hw_v1);
%     assignin('base','hw_v2',hw_v2);
%     assignin('base','hw_v3',hw_v3);
%     assignin('base','hw_v',hw_v);
%     assignin('base','Cpp',Cpp);
%     assignin('base','Ctp',Ctp);
%     assignin('base','Cpm',Cpm);
%     assignin('base','Ctm',Ctm);
end
end
