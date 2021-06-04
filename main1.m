function [Y, time, y0_,TEOUT1,IEOUT1]=main1(T,y0,celltype,ver,mutant)
%main
%12/27/2020
% clear all
% close all
global T_e1
T_e1=35;
if isempty(T)
load('T_6.mat');%load parameters
end
if isempty(y0)
 load('y0.mat');%
end
if isempty(celltype)
    celltype='SW';
end
global p;
parameters(1,T,ver,mutant);%all parameters

 %% simulation
 if strcmp(celltype,'SW')
     tspan=125;%120;
t0  =  0;       % Start time
tf  = tspan;%     % End time  %Z-ring closed
options  =  odeset('Events',@podJ_event,'RelTol',1e-4,'AbsTol',1e-6);

tout = t0;
y0 = y0.';
yout = y0;

teout  =  [];
yeout  =  [];
ieout  =  [];

while t0<tf
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
 elseif strcmp(celltype,'ST')
     tspan=91;%125-34
     
 t0  =  0;       % Start time 
tf  = tspan;% 
options  =  odeset('Events',@podJ_eventST,'RelTol',1e-4,'AbsTol',1e-6);

tout = t0;
y0=y0.';
yout = y0;

teout  =  [];
yeout  =  [];
ieout  =  [];

while t0<tf
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
%    if ie  ==  1%DNA replication initiates
%         T_e1 = te;
% else
    if ie == 2%Fork passes ctrA
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
 end
 %% simulation results1
 Y_output1=yout';
 time1=tout;
 L=length(yout(:,1));
M1(:,5)=yout(:,107);
M1(:,6)=2*yout(:,107);
M1(:,7)=3*yout(:,107);
M1(:,4)=zeros(1,L);
M1(:,1)=-3*yout(:,107);
M1(:,2)=-2*yout(:,107);
M1(:,3)=-yout(:,107);


M1=M1';
TEOUT1=teout; IEOUT1=ieout;
%% simulation results after z-ring closed

y0 = Y_output1(:,end);
Output2 = main_DIV(T,y0,ver,mutant);
Y_output2 = Output2.yout;
Y_output2 = Y_output2';
M2=Output2.grid;
time2=Output2.time;


Y=[Y_output1 Y_output2];
M=[M1 M2];
time=[time1 ;time2+tspan];
y0_=Y(:,end);

