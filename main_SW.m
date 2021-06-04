function output = main_SW(theta,y0,tspan,ver,mutant)
% close all
%1/24/2021
% clear all
global T_e1
T_e1=35;
if isempty(theta)
load('T_6.mat');%load parameters
end
if isempty(y0)
%load y0
 load('y0.mat');%SW IC 
end

global p;
parameters(1,theta,ver,mutant);%all parameters



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Integration parameters
t0  =  0;       % Start time
tf  = tspan;%120;%150;     % End time  %120min Z-ring closed

options  =  odeset('Events',@podJ_event,'RelTol',1e-4,'AbsTol',1e-6);

tout = t0;
y0 = y0.';
yout = y0;

teout  =  [];
yeout  =  [];
ieout  =  [];

while t0<tf
%  [t,y,te,ye,ie] = ode15s(@ODE_CtrA,[t0 tf],y0,options);
[t,y,te,ye,ie] = ode15s(@ODE_DivK,[t0 tf],y0,options);
 nt = length(t);
%     
 tout = [tout;t(2:nt)]; 
 yout = [yout;y(2:nt,:)]; 
 teout  =  [teout;te];  
 yeout  =  [yeout;ye]; ieout  =  [ieout;ie];
    y0  =  y(nt,:);
     if isscalar(ie)  ==  0
        ie  =  0;
     end
    if ie  ==  1%DNA replication initiates
        T_e1 = te;
%         disp(sprintf('T_e1= %8.5f',T_ei))
    elseif ie == 2%Fork passes ctrA
y0(104)=1;%SctrA changed from 0 to 1
    elseif ie==3
y0(105)=1;%SpleC
    elseif ie==4
        y0(106)=1;%SperP
    elseif ie ==5
        y0(103)=1;%SpodJ
         elseif ie==6
       y0(103:106)=0;
    end
    t0 = t(nt);
    if t0 >= tf
        break;
    end
end


%%
    
 %% construct grid for plotting
 L=length(yout(:,1));
M(:,5)=yout(:,107);
M(:,6)=2*yout(:,107);
M(:,7)=3*yout(:,107);
M(:,4)=zeros(1,L);
M(:,1)=-3*yout(:,107);
M(:,2)=-2*yout(:,107);
M(:,3)=-yout(:,107);

    

PodJ(:,1:6) = yout(:,1:6);%PodJm
PodJL(:,1:6)=yout(:,7:12);%PodJL polymer
PodJS(:,1:6)=yout(:,13:18);


SpmXm(:,1:6)=yout(:,19:24);
SpmXp(:,1:6)=yout(:,25:30);

PopZm(:,1:6)=yout(:,31:36);
PopZp(:,1:6)=yout(:,37:42);

CtrA(:,1:6)=yout(:,43:48);
CtrAP(:,1:6)=yout(:,49:54);

PleCf(:,1:6)=yout(:,55:60);
PleCb(:,1:6)=yout(:,61:66);

DivJf(:,1:6)=yout(:,67:72);
DivJb(:,1:6)=yout(:,73:78);

DivK(:,1:6)=yout(:,79:84);
DivKPT(:,1:6)=yout(:,85:90)+yout(:,91:96);

PerP(:,1:6)=yout(:,97:102);
%%%%%


PodJ = fliplr(PodJ);
PodJL=fliplr(PodJL);
PodJS=fliplr(PodJS);
SpmXm = fliplr(SpmXm);
SpmXp = fliplr(SpmXp);
PopZm = fliplr(PopZm);
PopZp = fliplr(PopZp);
CtrA = fliplr(CtrA);
CtrAP = fliplr(CtrAP);
PleCf = fliplr(PleCf);
PleCb = fliplr(PleCb);
DivJf = fliplr(DivJf);
DivJb = fliplr(DivJb);
DivK = fliplr(DivK);
DivKPT = fliplr(DivKPT);
PerP=fliplr(PerP);



PodJ = PodJ.';
M = M.';
PodJL = PodJL.';
PodJS = PodJS.';
PerP=PerP.';
SpmXm=SpmXm.';
SpmXp=SpmXp.';
PopZm=PopZm.';
PopZp=PopZp.';
CtrA=CtrA.';
CtrAP=CtrAP.';
PleCf=PleCf.';
PleCb=PleCb.';
DivJf=DivJf.';
DivJb=DivJb.';
DivK=DivK.';
DivKPT=DivKPT.';
PerP=PerP.';

output.PodJ = PodJ;
output.PodJL = PodJL;
output.PodJS = PodJS;
output.SpmXm = SpmXm;
output.SpmXp = SpmXp;
output.PopZm = PopZm;
output.PopZp = PopZp;
output.CtrA= CtrA;
output.CtrAP= CtrAP;
output.PleCf = PleCf;
output.PleCb= PleCb;
output.DivJf = DivJf;
output.DivJb= DivJb;
output.DivK= DivK;
output.DivKPT= DivKPT;
output.PerP = PerP;
output.time = tout;
output.grid = M;
 output.yout = yout;


