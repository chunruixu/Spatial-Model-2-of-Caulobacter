function graphcellcycle4(yout,tout,celltype,mutant,TITLE,spatial,temporal,separate)
% close all
%spatial plots
%%simulation
PodJLp=yout(7:12,:);
PodJLp = flipud(PodJLp);
PodJm=yout(1:6,:);
PodJm = flipud(PodJm);
PodJL=PodJLp+PodJm;
PodJS=yout(13:18,:);
PodJS = flipud(PodJS);
PopZp = yout(37:42,:);
PopZp = flipud(PopZp);
PopZ = yout(31:36,:)+yout(37:42,:);
PopZ = flipud(PopZ);
SpmX=yout(19:24,:)+yout(25:30,:);
SpmX = flipud(SpmX);

DivJb=yout(73:78,:);
DivJ=yout(67:72,:)+DivJb;
DivJ = flipud(DivJ);
DivJb = flipud(DivJb);
PleC=yout(55:60,:)+yout(61:66,:);
PleC = flipud(PleC);
DivK=yout(79:84,:);
DivK = flipud(DivK);
DivKP=yout(85:90,:)+yout(91:96,:);
DivKP = flipud(DivKP);
CtrAP=yout(49:54,:);
CtrAP = flipud(CtrAP);
CtrA=yout(43:48,:);
CtrA = flipud(CtrA);
CtrAt = CtrA+CtrAP;
PerP=yout(97:102,:);
PerP=flipud(PerP);
%%Mgrid
Y=yout';
L=length(Y(:,1));
M(:,5)=Y(:,107);
M(:,6)=2*Y(:,107);
M(:,7)=3*Y(:,107);
M(:,4)=zeros(1,L);
M(:,1)=-3*Y(:,107);
M(:,2)=-2*Y(:,107);
M(:,3)=-Y(:,107);

M=M';
a=zeros(1,L);
time=tout;
%  if strcmp(celltype,'ST')
% time=time+34;
% end
%% plots
% ha = tight_subplot(2,4,[.06 .03],[.1 .04],[.02 .02]);
if spatial==1
figure()%%8 subplots heatmap1
set(gcf,'Name',TITLE);
ha = tight_subplot(4,2,[.1 .03],[.1 .04],[.1 .04]);

% set(gcf,'position',[400 100 800 500])
% set (gcf,'Position',[0,0,512,512])
axes(ha(1)); 
% plot(2,4,1)%PodJLp
PodJL(7,:)=a;
pcolor(time, M, PodJL)
shading flat
colorbar
caxis([0 0.6]);
% xlabel('time (min)')
title('(a) PodJL')

axes(ha(2));
% plot(2,4,2)%PodJS
PodJS(7,:)=a;
pcolor(time, M, PodJS)
shading flat
colorbar
caxis([0 0.6]);
% xlabel('time (min)')
title('(b) PodJS')
% subplot(2,4,3)%PopZp
axes(ha(3)); 
PopZ(7,:)=a;
pcolor(time, M, PopZ)
shading flat
colorbar
caxis([0 12]);
% xlabel('time (min)')
title('(c) PopZ')
% subplot(2,4,4)%SpmX_T
axes(ha(4)); 
SpmX(7,:)=a;
pcolor(time, M, SpmX)
shading flat
colorbar
% caxis([0 1.5]);
% xlabel('time (min)')
title('(d) SpmX')

% subplot(2,4,5)%PleC
axes(ha(5)); 
PleC(7,:)=a;
pcolor(time, M, PleC)
shading flat
colorbar
% caxis([0 20]);
% xlabel('time (min)')
title('(e) PleC')
% subplot(2,4,6)%DivJ
axes(ha(6)); 
DivJ(7,:)=a;
pcolor(time, M, DivJ)
shading flat
colorbar
% caxis([0 20]);
% xlabel('time (min)')
title('(f) DivJ')
% subplot(2,4,7)%DivKP
axes(ha(7)); 
DivKP(7,:)=a;
pcolor(time, M, DivKP)
shading flat
colorbar
% caxis([0 50]);
xlabel('time (min)')
title('(g) DivKP')

% subplot(2,4,8)%CtrAP
axes(ha(8)); 
CtrAP(7,:)=a;
pcolor(time, M, CtrAP)
shading flat
colorbar
% caxis([0 1]);
xlabel('time (min)')
title('(h) CtrAP')
set(findall(gcf,'-property','FontSize'),'FontSize',14)
end

if separate==1
    if strcmp(mutant,'deltaPodJ')
    figure()
    set(gcf,'position',[10 300 300 510])
subplot(3,1,1)
    PopZ(7,:)=a;
pcolor(time, M, PopZ)
shading flat
colorbar
caxis([0 12]);
% xlabel('time (min)')
title('PopZ')
subplot(3,1,2)
SpmX(7,:)=a;
pcolor(time, M, SpmX)
shading flat
colorbar
% caxis([0 1.5]);
% xlabel('time (min)')
title('SpmX')

subplot(3,1,3)
PleC(7,:)=a;
pcolor(time, M, PleC)
shading flat
colorbar
% caxis([0 20]);
xlabel('time (min)')
title('PleC')
    elseif strcmp(mutant,'PodJ+')
            figure()
    set(gcf,'position',[10 300 500 210])
    subplot(1,2,1)
    PopZ(7,:)=a;
pcolor(time, M, PopZ)
shading flat
colorbar
caxis([0 12]);
xlabel('time (min)')
title('PopZ')
subplot(1,2,2)
PleC(7,:)=a;
pcolor(time, M, PleC)
shading flat
colorbar
% caxis([0 20]);
xlabel('time (min)')
title('PleC')
    end
end

%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
tspan=[tout(1) tout(end)];
%% temperal DATA
SpmX = [0 2340.397; 20 11034.71; 40 13828.5; 60 14676.86; 80 13635.98; 100 9316.033; 120 8849.205; 140 12085.38; 160 13972.28];
SpmX(:,2) = SpmX(:,2)/max(SpmX(:,2));
dpSpmX=SpmX;%paper152
SpmX2 = [0 2857.246; 20 10014.882; 40 12674.083;60 13516.276;80 12773.64; 100 8493.619; 120 8413.962; 140 11470.74;160 11327.459];
SpmX2(:,2) = SpmX2(:,2)/max(SpmX2(:,2));
dpSpmX2=SpmX2;%paper118
SpmX3 =[0 1506.527;20 4945.134;40 8203.305;60 9740.79;80 10103.205;100 10964.376;120 10477.861;140 9446.033;160 4349.246];%paper126
SpmX3(:,2) = SpmX3(:,2)/max(SpmX3(:,2));
dpSpmX3=SpmX3;
PodJL=[0 0; 20 0; 40 1726.012; 60 12068.518; 80 23215.326; 100 22996.104; 120  9466.861; 140  3336.154; 160 5506.548];
PodJL(:,2) = PodJL(:,2)/max(PodJL(:,2));
dpPodJL=PodJL;%paper116
PodJL2=[0 790.719;20 1297.841;40 4612.861;60 9123.397;80 11505.154;100 13792.497;120 15360.276;140 8009.589;160 5775.439];
PodJL2(:,2) = PodJL2(:,2)/max(PodJL2(:,2));
dpPodJL2=PodJL2;%126
PodJS=[0  16869.841; 20  16960.447; 40  4387.841; 60  434.364;  80  0; 100  2627.527; 120  8780.79; 140  8984.033; 160  6333.598];
PodJS(:,2) = PodJS(:,2)/max(PodJS(:,2));
dpPodJS=PodJS;%paper116
PodJS2=[0 7656.953;20 11798.569;40 9196.033;60 7763.104;80 7449.276;100 8469.74; 120 11575.711; 140 11777.761; 160 10382.853];
PodJS2(:,2) = PodJS2(:,2)/max(PodJS2(:,2));
dpPodJS2=PodJS2;
dpDivKTOT =[0 0.67; 20 0.77; 40 0.88; 60 0.87; 80 0.72; 100 0.71; 120 0.89; 140 1];
DivKTOT2=[0 4395.619;20 5939.134;40 6636.012;60 6630.841;80 6025.305;100 5316.184;120 6691.891;140 7454.719;160 5877.569];
DivKTOT2(:,2) = DivKTOT2(:,2)/max(DivKTOT2(:,2));
dpDivKTOT2=DivKTOT2;
dpDivJ = [0 0.364; 20 0.574; 40 0.81; 60 0.997; 80 0.924; 100 0.9; 120 1; 140 0.94]; % Wheeler et al. 1999
dpDivJ2 = [0 0.26; 20 0.56; 40 0.86; 60 0.81; 80 0.97; 100 0.93; 120 0.86; 140 1]; %Sanselicio 2015
DivJ3 = [0 2638.933; 20 5578.497; 40 7750.447;60 10436.912;80 11503.497; 100 10182.79; 120 13502.083; 140 12420.326;160 13399.489];
DivJ3(:,2) = DivJ3(:,2)/max(DivJ3(:,2));%paper118
dpDivJ3=DivJ3;
dpPleCtot = [0 .85/0.85; 20 0.5/0.85; 40 .28/0.85; 60 .38/0.85; 80 .45/0.85; 100 .65/0.85; 120 0.75/0.85; 140 0.8/0.85]; %Viollier 2002 - ommitted last datapoint bc it was at time point 160 and is unlikely behavior
PleCtot2 =[0 7535.004;20 4253.184; 40 3452.841; 60 5912.426; 80 9557.225; 100 11047.447; 120 11043.326; 140 11704.104; 160 6926.953];
PleCtot2(:,2) = PleCtot2(:,2)/max(PleCtot2(:,2));
dpPleCtot2=PleCtot2;
dpCtrAP = [10*1.07 0.4; 80*1.07 0.5; 105*1.07 .8; 130*1.07 1]; %Jacobs 2003 (relative to max = 1) removed second point (t=40, P=0.15) as CtrAT suggests it is incorrect
Mcgrath = [0	12915.912; 20	9480.648; 40	2687.648; 60	515.698; 80	5769.719; 100	10582.426; 120	13602.205;140	13106.439];
Mcgrath(:,2)=Mcgrath(:,2)/max(Mcgrath(:,2));
dpCtrAT = Mcgrath;%SLOW
dpCtrAT2 = [0 0.8; 20 0; 40 0; 60 0; 80 0.5; 100 1; 120 .85; 140 0.7]; %Collier 2006 - DnaA couples DNA... (Normalized based on max expression)
%QUICK
CtrAT3=[0 16004.861;20 7879.962;40 1532.690;60 939.447;80 5083.154;100 11343.933; 120 12177.326; 140 13739.397; 160 13562.326];
CtrAT3(:,2) = CtrAT3(:,2)/max(CtrAT3(:,2));
dpCtrAT3=CtrAT3;
if temporal==1
%% temporal plots
%%
LegendSize=9; LabelSize=14;

figure()
set(gcf,'Name',TITLE);
set(gcf,'position',[10 50 500 760])
subplot(3,2,1)
SpmX=yout(19:24,:)+yout(25:30,:);
SpmX=sum(SpmX).*yout(107,:);
p5 = line(tout,SpmX , 'Color', 'k' , 'LineWidth', 2, 'Linestyle', '-');
% SpmX=(yout(13,:)+yout(16,:)+yout(17,:)+yout(20,:)).*yout(107,:)+(yout(14,:)+yout(15,:)+yout(18,:)+yout(19,:)).*yout(74,:);
% p5 = line(tout,SpmX, 'Color', 'k' , 'LineWidth', 2, 'Linestyle', '-');
% hold on
% [score,scalary] = findSquares(tout,sum(SpmX), dpSpmX(1:end-1,:));
% p7 = scatter(dpSpmX(:,1),dpSpmX(:,2).*scalary, 'k+', 'LineWidth', 1);
hold on
% [score,scalary] = findSquares(tout,sum(SpmX), dpSpmX2(1:end-1,:));
[score,scalary] = findSquares(tout,SpmX, dpSpmX2(1:end-1,:));
p7 = scatter(dpSpmX2(:,1),dpSpmX2(:,2).*scalary, 'k', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,SpmX, dpSpmX3(1:end-1,:));
p7 = scatter(dpSpmX3(:,1),dpSpmX3(:,2).*scalary, 'k^', 'LineWidth', 1);
hold off
axis([tspan 0 max(SpmX)*1.1]);
% xlabel('time (min)')
% ylabel('scaled concentration')
h = legend('SpmX_T','Radhakrishnan 07','Guo 04','Location', 'North','fontsize',LegendSize);
% h = legend('SpmX_T','Kaczmarczyk et al. ''20','Radhakrishnan et al.''07','Guo 04','Location', 'North','fontsize',LegendSize);
  set(h, 'Interpreter', 'tex', 'Box', 'off', 'Orientation', 'vertical');%horizontal
%    set(gca,'FontSize',12);

subplot(3,2,2)
PodJL=yout(1:6,:)+yout(7:12,:);   PodJS=yout(13:18,:);
PodJL_T=sum(PodJL).*yout(107,:);  PodJS_T=sum(PodJS).*yout(107,:);
% PodJL_T=(PodJL(1,:)+PodJL(4,:)).*yout(73,:)+(PodJL(2,:)+PodJL(3,:)).*yout(74,:);  PodJS_T=(PodJS(1,:)+PodJS(4,:)).*yout(73,:)+(PodJS(2,:)+PodJS(3,:)).*yout(74,:); 
PodJ_T=PodJL_T+PodJS_T;
Fraction_L=PodJL_T./PodJ_T;
Fraction_S=PodJS_T./PodJ_T;
%%data
Time=[0 20 40 60 80 100 120 140 160];
%%Calculate fraction of relative levels by getdata
L=[100 91 19 2 1 12 46 32 53];
S=[0 0 7 50 100 100 41 14 22];
Sum=L+S;

p5 = line(tout,PodJL_T , 'Color', 'k' , 'LineWidth', 2, 'Linestyle', '-');
p5 = line(tout,PodJS_T , 'Color', 'b' , 'LineWidth', 2, 'Linestyle', '-');
hold on
[score,scalary] = findSquares(tout,PodJL_T, dpPodJL(1:end-1,:));
p6 = scatter(dpPodJL(:,1),dpPodJL(:,2).*scalary, 'k+', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,PodJL_T, dpPodJL2(1:end-1,:));
p6 = scatter(dpPodJL2(:,1),dpPodJL2(:,2).*scalary, 'k^', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,PodJS_T, dpPodJS(1:end-1,:));
p7 = scatter(dpPodJS(:,1),dpPodJS(:,2).*scalary, 'b+', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,PodJS_T, dpPodJS2(1:end-1,:));
p7 = scatter(dpPodJS2(:,1),dpPodJS2(:,2).*scalary, 'b^', 'LineWidth', 1);
hold off
axis([tspan 0 max(PodJL_T+PodJS_T)*1.1]);
% xlabel('time (min)')
% ylabel('scaled concentration')
% h = legend('SpmX_T','Kaczmarczyk\newlineet al. ''20','Location', 'North');
h = legend('PodJL','PodJS','Chen 05','Guo 04','Chen 05','Guo 04','Location', 'North','fontsize',LegendSize);%,'NumColumns',3);
  set(h, 'Interpreter', 'tex', 'Box', 'off', 'Orientation', 'vertical');
%    set(gca,'FontSize',12);
   
subplot(3,2,3)

  DivJf=yout(67:72,:); DivJb=yout(73:78,:);
DivJf=sum(DivJf).*yout(107,:); DivJb=sum(DivJb).*yout(107,:);
%   DivJf=(yout(45,:)+yout(48,:)).*yout(73,:)+(yout(46,:)+yout(47,:)).*yout(74,:); DivJb=(yout(49,:)+yout(52,:)).*yout(73,:)+(yout(50,:)+yout(51,:)).*yout(74,:);
DivJT=DivJf+DivJb;

p5 = line(tout,DivJT , 'Color', 'k' , 'LineWidth', 2, 'Linestyle', '-');
p5 = line(tout,DivJf , 'Color', 'r' , 'LineWidth', 2, 'Linestyle', '--');
p5 = line(tout,DivJb , 'Color', 'r' , 'LineWidth', 2, 'Linestyle', ':');
hold on
[score,scalary] = findSquares(tout,DivJT, dpDivJ(1:end-1,:));
p7 = scatter(dpDivJ(:,1),dpDivJ(:,2).*scalary, 'k+', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,DivJT, dpDivJ2(1:end-1,:));
p7 = scatter(dpDivJ2(:,1),dpDivJ2(:,2).*scalary, 'k^', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,DivJT, dpDivJ3(1:end-1,:));
p7 = scatter(dpDivJ3(:,1),dpDivJ3(:,2).*scalary, 'k*', 'LineWidth', 1);
hold off
axis([tspan 0 max(DivJT)*1.1]);
% ylabel('scaled concentration')
h = legend('DivJ_T','DivJ_f','DivJ_b','Wheeler 1999','Sanselicio 2015','Radhakrishnan 07','Location', 'North','fontsize',LegendSize);%,'NumColumns',3);
 set(h, 'Interpreter', 'tex', 'Box', 'off', 'Orientation', 'vertical');
%    set(gca,'FontSize',12);

subplot(3,2,4) 
% subplot(4,1,3)%PleC
  PleCf=yout(55:60,:); PleCb=yout(61:66,:);
PleCf=sum(PleCf).*yout(107,:); PleCb=sum(PleCb).*yout(107,:);
%   PleCf=(yout(37,:)+yout(40,:)).*yout(73,:)+(yout(38,:)+yout(39,:)).*yout(74,:); PleCb=(yout(41,:)+yout(44,:)).*yout(73,:)+(yout(42,:)+yout(43,:)).*yout(74,:);
PleCT=PleCf+PleCb;
p5 = line(tout,PleCT , 'Color', 'k' , 'LineWidth', 2, 'Linestyle', '-');
p5 = line(tout,PleCf , 'Color', 'r' , 'LineWidth', 2, 'Linestyle', '--');
p5 = line(tout,PleCb , 'Color', 'r' , 'LineWidth', 2, 'Linestyle', ':');
hold on
[score,scalary] = findSquares(tout,PleCT, dpPleCtot(1:end-1,:));
p7 = scatter(dpPleCtot(:,1),dpPleCtot(:,2).*scalary, 'k+', 'LineWidth', 2);
hold on
[score,scalary] = findSquares(tout,PleCT, dpPleCtot2(1:end-1,:));
p7 = scatter(dpPleCtot2(:,1),dpPleCtot2(:,2).*scalary, 'k^', 'LineWidth', 2);
hold off
axis([tspan 0 max(PleCT)*1.1]);
% xlabel('time (min)','fontsize',LabelSize)
% ylabel('scaled concentration')
h = legend('PleC_T','PleC_f','PleC_b','Viollier 2002','Guo 04','Location', 'North','fontsize',LegendSize);
 set(h, 'Interpreter', 'tex', 'Box', 'off', 'Orientation', 'vertical');
%    set(gca,'FontSize',12);
   
subplot(3,2,5) 
%    subplot(4,1,4)%DivK
  DivK=yout(79:84,:); DivKP=yout(85:90,:)+yout(91:96,:);
DivK=sum(DivK).*yout(107,:);  DivKP=sum(DivKP).*yout(107,:);
%     DivK=(yout(53,:)+yout(56,:)).*yout(73,:)+(yout(54,:)+yout(55,:)).*yout(74,:); DivKP=(yout(57,:)+yout(60,:)+yout(61,:)+yout(64,:)).*yout(73,:)+(yout(58,:)+yout(59,:)+yout(62,:)+yout(63,:)).*yout(74,:);
DivKT=DivK+DivKP;
p5 = line(tout,DivKT , 'Color', 'k' , 'LineWidth', 2, 'Linestyle', '-');
p5 = line(tout,DivK , 'Color', 'r' , 'LineWidth', 2, 'Linestyle', '--');
p5 = line(tout,DivKP , 'Color', 'r' , 'LineWidth', 2, 'Linestyle', ':');
% hold on
% [score,scalary] = findSquares(tout,sum(DivKT), dpDivKTOT(1:end-1,:));
% p7 = scatter(dpDivKTOT(:,1),dpDivKTOT(:,2).*scalary, 'k+', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,DivKT, dpDivKTOT2(1:end-1,:));
p7 = scatter(dpDivKTOT2(:,1),dpDivKTOT2(:,2).*scalary, 'k^', 'LineWidth', 1);
hold off
axis([tspan 0 max(DivKT)*1.1]);
xlabel('time (min)','fontsize',LabelSize)
% ylabel('scaled concentration')
h = legend('DivK_T','DivK','DivK~P','Gao 04','Location', 'North','fontsize',LegendSize);
 set(h, 'Interpreter', 'tex', 'Box', 'off', 'Orientation', 'vertical');
%    set(gca,'FontSize',12);

subplot(3,2,6) 
% subplot(4,1,1)%temporal CtrA
CtrAP = yout(49:54,:); CtrAP=sum(CtrAP).*yout(107,:);
CtrA = yout(55:60,:); CtrA=sum(CtrA).*yout(107,:);
CtrAT=CtrA+CtrAP;
p5 = line(tout,CtrAT , 'Color', 'k' , 'LineWidth', 2, 'Linestyle', '-');
p5 = line(tout,CtrAP , 'Color', 'r' , 'LineWidth', 2, 'Linestyle', '-');
hold on
[score,scalary] = findSquares(tout,CtrAT, dpCtrAT(1:end-1,:));
p7 = scatter(dpCtrAT(:,1),dpCtrAT(:,2).*scalary, 'k+', 'LineWidth', 1);
% hold on
% [score,scalary] = findSquares(tout,sum(CtrAT), dpCtrAT2(1:end-1,:));
% p7 = scatter(dpCtrAT2(:,1),dpCtrAT2(:,2).*scalary, 'k^', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,CtrAT, dpCtrAT3(1:end-1,:));
p7 = scatter(dpCtrAT3(:,1),dpCtrAT3(:,2).*scalary, 'k*', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,CtrAP, dpCtrAP(1:end-1,:));
p7 = scatter(dpCtrAP(:,1),dpCtrAP(:,2).*scalary, 'r+', 'LineWidth', 1);
hold off
axis([tspan 0 max(CtrAT)*1.1]);
xlabel('time (min)','fontsize',LabelSize)
% ylabel('scaled concentration')
% h = legend('SpmX_T','Kaczmarczyk\newlineet al. ''20','Location', 'North');
% h = legend('CtrA_T','CtrA~P','Mcgraph','Collier et al. ''06','Jacobs et al. ''03','Location', 'North');
h = legend('CtrA_T','CtrA~P','Mcgraph','Radhakrishnan 07','Jacobs 03','Location', 'North','fontsize',LegendSize);
set(h, 'Interpreter', 'tex', 'Box', 'off', 'Orientation', 'vertical');
%    set(gca,'FontSize',12);
% set(findall(gcf,'-property','FontSize'),'FontSize',14)
set(gca,'ytick',[])
end

% figure()
% subplot(2,1,1)
% PodJL=yout(1:4,:)+yout(5:8,:);   PodJS=yout(9:12,:);
% PodJL_T=sum(PodJL);  PodJS_T=sum(PodJS);
% PodJ_T=PodJL_T+PodJS_T;
% Fraction_L=PodJL_T./PodJ_T;
% Fraction_S=PodJS_T./PodJ_T;
% %%data
% Time=[0 20 40 60 80 100 120 140 160];
% %%Calculate fraction of relative levels by getdata
% L=[100 91 19 2 1 12 46 32 53];
% S=[0 0 7 50 100 100 41 14 22];
% Sum=L+S;
% 
% p5 = line(tout,PodJL_T , 'Color', 'k' , 'LineWidth', 2, 'Linestyle', '-');
% p5 = line(tout,PodJS_T , 'Color', 'b' , 'LineWidth', 2, 'Linestyle', '-');
% hold on
% [score,scalary] = findSquares(tout,PodJL_T, dpPodJL(1:end-1,:));
% p6 = scatter(dpPodJL(:,1),dpPodJL(:,2).*scalary, 'k', 'LineWidth', 1.5);
% hold on
% [score,scalary] = findSquares(tout,PodJL_T, dpPodJL2(1:end-1,:));
% p6 = scatter(dpPodJL2(:,1),dpPodJL2(:,2).*scalary, 'k^', 'LineWidth', 1.5);
% hold on
% [score,scalary] = findSquares(tout,PodJS_T, dpPodJS(1:end-1,:));
% p7 = scatter(dpPodJS(:,1),dpPodJS(:,2).*scalary, 'b', 'LineWidth', 1.5);
% hold on
% [score,scalary] = findSquares(tout,PodJS_T, dpPodJS2(1:end-1,:));
% p7 = scatter(dpPodJS2(:,1),dpPodJS2(:,2).*scalary, 'b^', 'LineWidth', 1.5);
% hold off
% axis([tspan 0 max(PodJL_T+PodJS_T)*1.1]);
% % xlabel('time (min)')
% ylabel('scaled concentration')
% % h = legend('SpmX_T','Kaczmarczyk\newlineet al. ''20','Location', 'North');
% h = legend('PodJL','PodJS','Chen 05','Guo 04','Chen 05','Guo 04','Location', 'North','fontsize',14,'NumColumns',3);
%   set(h, 'Interpreter', 'tex', 'Box', 'on', 'Orientation', 'vertical');
%    set(gca,'FontSize',14);
%    
% subplot(2,1,2)
% PodJL=yout(1:4,:)+yout(5:8,:);   PodJS=yout(9:12,:);
% PodJL_T=sum(PodJL);  PodJS_T=sum(PodJS);
% PodJ_T=PodJL_T+PodJS_T;
% Fraction_L=PodJL_T./PodJ_T;
% Fraction_S=PodJS_T./PodJ_T;
% %%data
% Time=[0 20 40 60 80 100 120 140 160];
% %%Calculate fraction of relative levels by getdata
% L=[100 91 19 2 1 12 46 32 53];
% S=[0 0 7 50 100 100 41 14 22];
% Sum=L+S;
% L2=PodJL2(:,2);
% S2=PodJS2(:,2);
% Sum2=L2+S2;
% plot(tout,Fraction_L,'k','linewidth',2);
% hold on;
% plot(tout,Fraction_S,'b','linewidth',2);
% hold on;
% scatter(Time,S./Sum,'k','LineWidth',1.5);
% hold on;
% scatter(Time,L./Sum,'b','LineWidth',1.5);
% hold on;
% scatter(Time,S2./Sum2,'b','^','LineWidth',1.5);
% hold on;
% scatter(Time,L2./Sum2,'k','^','LineWidth',1.5);
% axis([tspan 0 1.1]);
% xlabel('time (min)')
% ylabel('fraction')
% % h=legend('Simulation: fraction of PodJL','Simulation: fraction of PodJS','Experimental fraction of PodJL','Experimental fraction of PodJS');
% h=legend('Simulated PodJ_L','Simulated PodJ_S','Chen et al.','Chen et al.','Guo et al.','Guo et al.','fontsize',14);
% set(h, 'Interpreter', 'tex', 'Box', 'on', 'Orientation', 'horizontal');
%    set(gca,'FontSize',14)
   


