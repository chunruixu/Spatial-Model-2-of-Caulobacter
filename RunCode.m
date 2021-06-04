%%runcode
%%1/16/2021

clear all
close all
%% first cycle
load('T_6.mat');%load parameters
%% version of params
ver=7;  %4 is from the 4compartment folder 

 celltype='SW';
 mutant='WT';%'deltaPodJ';
%  mutant='deltaPodJ';
%initial values
y0=zeros(107,1);%SW IC - first cell cycle
    y0(7:12)=10e-6;%PodJp
y0(13:17)=0.001; y0(18)=0.1;%PodJS
y0(37:41)=0.5; y0(42)=2;%25;%PopZp
y0(43:48)=0.2; y0(49:54)=0.5;%CtrA and CtrAP
y0(55:60)=0.05; %PleCf
y0(66)=0.05;%PleCb
y0(79:84)=0.2;%DivK
y0(107)=0.02*100/6;%length of polar and central compartment

if strcmp(celltype,'ST')
     output=main_SW(T,y0,34,ver,mutant);
     yout=output.yout; yout=yout';
     y0=yout(:,end);
end

 CycleNum=5;%# of simulated cycles
 for i=1:CycleNum
   TITLE= [num2str(i) 'cellcyle'];
[Y, time, y0_,TE,IE]=main1(T,y0,celltype,ver,mutant);%simulation
graphcellcycle(Y,time,celltype,mutant,TITLE,1,0,1)%total graph
y0=IniValue(Y,celltype);%update y0 of next cycle
 end
% [Y, time, y0_,TE,IE]=main1(T,y0,celltype,ver,mutant);%simulation