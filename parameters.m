function parameters(cond,theta,ver,mutant)

% Global variables
global p;
%cond is for growing or ungrowing cell
% Growth conditions
if cond == 0              % cell size is fixed
    p.mu = 0;             % growth rate constant    
elseif cond == 1          % cell is growing 
    p.mu = 0.0053;%0.00526;%0.0055;        % growth rate constant = 0.0055 (units => 1/min)    
end

if ver==2
 % Diffusion parameters (units => um^2/min)
p.D_podJm =100;%for PopZ only & 100 in manuscript;  
p.D_podJL =0.0005;%40;%0.0005;  

p.syn_podJ = 0.01;%0.2*theta(1);   % monomer synthesis 0.2
p.syn_podJ2 = 0.006;
p.Ji_PodJCtrA = 0.5;
p.n_PodJCtrA = 4;
% p.deg_podJm =theta(2);% 0.05;     0.046 paper116
% p.deg_podJp =theta(2);% 0.05;       
p.deg_podJm =0.2;% 0.05;    0.01?- 0.046 paper116
p.deg_podJp =p.deg_podJm;% 0.05;    
p.deg_podJ1 =0.002;%0.007;%0.007;%paper116

p.dnv_podJ =0.01*theta(3);% 15              % denovo polymerization 12 in dissertation for table C.2; 1 for table C.1
p.aut1_podJ =50*theta(4);%10;               % autocatalytic polymerization (pole)
p.aut1_podJ1=0*0.1*p.aut1_podJ;%central compartment autcatalytic polymerization
% p.depol_basal=0.2;
p.depol_podJ =0.2*theta(5);%0.5;             % deplymerization 

p.podj=2;
p.deg_s=0.5*theta(6);%0.05  %0.01 paper116ref
   p.alpha_PodJSpmX =100*theta(17);%2; %SpmX on PodJ dnv
%% PopZ
p.syn_popz=0.1*theta(8);%1.5;% popZ monomer synthesis
    p.deg_popzm=0.01;%theta(9);%0.05 % popZ monomer degradation
    p.deg_popzp=p.deg_popzm;%theta(9);%0.05;
    p.dnv_popz=0.21*theta(10);%60;% denovo polymerization
    p.aut_popz=3;%3*theta(11);%3;% autocatalytic polymerization
    p.alpha_PopZPodJ=50*theta(12);%1;%PodJ on PopZ
    p.depol_popz=0.5*theta(13);%0.1; % deplymerization 
    % diffussion rates
    p.D_popzm =835;%750; % PopZ monomer dffusion  
p.D_popzp =0.0005;%40;%0.0005;  % PopZ polymer diffusion
   %% SpmX
   p.syn_spmx =6*theta(15);%0.04;
   p.deg_spmx = 0.1*theta(16);%0.1;
     p.dnv_spmx =1*theta(18);%0.1;%
   p.depol_spmx = 0.36*p.dnv_spmx;%theta(18);%
   p.aut_spmx=1*theta(18);
   p.Ja_SpmXCtrA=1;
   p.alpha_SpmXPopZ=50;
p.D_spmx =200;%um^2/min
p.D_spmxp=0.0005;
%% CtrA
p.syn_ctrA1=0.008;%0.0083;%Shenghua2008-exp; 0.026;%Murray-fitted
p.Ji_CtrACtrA=3;
p.syn_ctrA2=0.08;%0.073;%Shenghua2008-exp
p.Ja_CtrACtrA=3;
p.deg_ctrA1=0.0038;%exp-calculated
p.deg_ctrA2=0.12;%exp-estimated
p.Jd_CtrA=0.2;
p.dephoCtrA=100000*log(2)/5/2;%Murray-estimated~0.14
p.phoCtrA=100000*0.03;%7/93*p.dephoCtrA/0.5%~0.021;%Murray-estimated
p.D_CtrA=427;
p.D_CtrAP=p.D_CtrA;
%% PleC
p.syn_pleC=0.01;%0.053%Bronson
p.deg_pleC=0.02;%0.028;%Bronson
p.fb_PleC=10*0.1;
p.alpha_PleCPodJ=10;%0.5;
p.bf_PleC=10*0.05;
p.D_PleC=71;
%% DivJ
p.syn_divJ=0.008*2;%0.005;%Bronson
p.deg_divJ=0.035*2;%0.035;%Bronson
p.fb_DivJ=1;
p.alpha_DivJSpmX=5;
p.bf_DivJ=0.5;
p.D_DivJ=108;
%% DivK
p.syn_divK1=0.0004;%0.0004;%Bronson
p.syn_divK2=0.05*4;%0.125;%Bronson
p.Ja_DivKCtrA=0.5;%1.7/2;%Bronson
p.deg_divK=0.014*4;%Bronson; 0.002;%Shenghua2008-exp
p.phoDivK=0.01;
p.alpha_DivKDivJ=50;
p.dephoDivK=1;
p.D_DivK=1319;
p.D_DivKP=p.D_DivK;%p.D_DivK;

p.fb_DivKP=10000;
p.bf_DivKP=0.1;
%% PerP
p.syn_perP=0.5;%0.0428;%Bronson
p.Ja_PerPCtrA=1/2;%2.8/2;%Bronon
p.deg_perP=0.04;%0.04;%Bronson
p.D_PerP=853;
elseif ver==3
      % Diffusion parameters (units => um^2/min)
p.D_podJm =100;%for PopZ only & 100 in manuscript;  
p.D_podJL =0.0005;%40;%0.0005;  

p.syn_podJ = 0.01;%0.2*theta(1);   % monomer synthesis 0.2
p.syn_podJ2 = 0.006;
p.Ji_PodJCtrA = 0.5;
p.n_PodJCtrA = 4;
% p.deg_podJm =theta(2);% 0.05;     0.046 paper116
% p.deg_podJp =theta(2);% 0.05;       
p.deg_podJm =0.3;% 0.05;    0.01?- 0.046 paper116
p.deg_podJp =p.deg_podJm;% 0.05;    
p.deg_podJ1 =0.007;%0.007;%0.007;%paper116

p.dnv_podJ =0.01*theta(3);% 15              % denovo polymerization 12 in dissertation for table C.2; 1 for table C.1
p.aut1_podJ =50*theta(4);%10;               % autocatalytic polymerization (pole)
p.aut1_podJ1=0*0.1*p.aut1_podJ;%central compartment autcatalytic polymerization
% p.depol_basal=0.2;
p.depol_podJ =0.2*theta(5);%0.5;             % deplymerization 

p.podj=2;
p.deg_s=0.5*theta(6);%0.05  %0.01 paper116ref
   p.alpha_PodJSpmX =100*theta(17);%2; %SpmX on PodJ dnv
%% PopZ
p.syn_popz=0.1*theta(8);%1.5;% popZ monomer synthesis
    p.deg_popzm=0.01;%theta(9);%0.05 % popZ monomer degradation
    p.deg_popzp=p.deg_popzm;%theta(9);%0.05;
    p.dnv_popz=0.21*theta(10);%60;% denovo polymerization
    p.aut_popz=3;%3*theta(11);%3;% autocatalytic polymerization
    p.aut_popz1=p.aut_popz;
    p.alpha_PopZPodJ=50*theta(12);%1;%PodJ on PopZ
    p.depol_popz=0.5*theta(13);%0.1; % deplymerization 
    % diffussion rates
    p.D_popzm =835;%750; % PopZ monomer dffusion  
p.D_popzp =0.0005;%40;%0.0005;  % PopZ polymer diffusion
   %% SpmX
   p.syn_spmx =6*theta(15);%0.04;
   p.deg_spmx = 0.1*theta(16);%0.1;
     p.dnv_spmx =1*theta(18);%0.1;%
   p.depol_spmx = 0.36*p.dnv_spmx;%theta(18);%
   p.aut_spmx=1*theta(18);
   p.Ja_SpmXCtrA=1;
   p.alpha_SpmXPopZ=50;
p.D_spmx =200;%um^2/min
p.D_spmxp=0.0005;
%% CtrA
p.syn_ctrA1=0.008;%0.0083;%Shenghua2008-exp; 0.026;%Murray-fitted
p.Ji_CtrACtrA=3;
p.syn_ctrA2=0.073;%0.073;%Shenghua2008-exp
p.Ja_CtrACtrA=3;
p.deg_ctrA1=0.0038;%0.0038;%exp-calculated
p.deg_ctrA2=0.12;%0.1 exp-estimated
p.Jd_CtrA=0.2;
p.dephoCtrA=100000*log(2)/8;%Murray-estimated~0.14
p.phoCtrA=100000*0.03;%7/93*p.dephoCtrA/0.5%~0.021;%Murray-estimated
p.D_CtrA=427;
p.D_CtrAP=p.D_CtrA;
%% PleC
p.syn_pleC=0.01;%0.053%Bronson
p.deg_pleC=0.02;%0.028;%Bronson
p.fb_PleC=10*0.1;
p.alpha_PleCPodJ=10;%0.5;
p.bf_PleC=10*0.05;
p.D_PleC=71;
%% DivJ
p.syn_divJ=0.008*2;%0.005;%Bronson
p.deg_divJ=0.035*2;%0.035;%Bronson
p.fb_DivJ=1;
p.alpha_DivJSpmX=5;
p.bf_DivJ=0.5;
p.D_DivJ=108;
%% DivK
p.syn_divK1=0.0004;%0.0004;%Bronson
p.syn_divK2=0.05*3;%0.125;%Bronson
p.Ja_DivKCtrA=0.7;%1.7/2;%Bronson
p.deg_divK=0.014*3;%Bronson; 0.002;%Shenghua2008-exp
p.phoDivK=0.01;
p.alpha_DivKDivJ=50;
p.dephoDivK=1;
p.D_DivK=1319;
p.D_DivKP=p.D_DivK;%p.D_DivK;

p.fb_DivKP=10000;
p.bf_DivKP=0.1;
%% PerP
p.syn_perP=0.5;%0.0428;%Bronson
p.Ja_PerPCtrA=1/2;%2.8/2;%Bronon
p.deg_perP=0.04;%0.04;%Bronson
p.D_PerP=853;
elseif ver==4  %%%%%%%%%%%%for mutant analysis
      p.D_podJm =100;%for PopZ only & 100 in manuscript;  
p.D_podJL =0.0005;%40;%0.0005;  

p.syn_podJ = 0.01;%0.2*theta(1);   % monomer synthesis 0.2
p.syn_podJ2 = 0.006;
p.Ji_PodJCtrA = 0.5;
p.n_PodJCtrA = 4;
% p.deg_podJm =theta(2);% 0.05;     0.046 paper116
% p.deg_podJp =theta(2);% 0.05;       
p.deg_podJm =0.3;% 0.05;    0.01?- 0.046 paper116
p.deg_podJp =p.deg_podJm;% 0.05;    
p.deg_podJ1 =0.007;%0.007;%0.007;%paper116

p.dnv_podJ =0.1*theta(3);% 15              % denovo polymerization 12 in dissertation for table C.2; 1 for table C.1
p.aut1_podJ =30*theta(4);%10;               % autocatalytic polymerization (pole)
p.aut1_podJ1=0*0.1*p.aut1_podJ;%central compartment autcatalytic polymerization
% p.depol_basal=0.2;
p.depol_podJ =0.2*theta(5);%0.5;             % deplymerization 

p.podj=2;
p.deg_s=0.5*theta(6);%0.05  %0.01 paper116ref
   p.alpha_PodJSpmX =100*theta(17);%2; %SpmX on PodJ dnv
%% PopZ
p.syn_popz=0.1*theta(8);%1.5;% popZ monomer synthesis
    p.deg_popzm=0.01;%theta(9);%0.05 % popZ monomer degradation
    p.deg_popzp=p.deg_popzm;%theta(9);%0.05;
    p.dnv_popz=0.7*theta(10);%60;% denovo polymerization
    p.aut_popz=12;%3*theta(11);%3;% autocatalytic polymerization
     p.aut_popz1=0.7*p.aut_popz;
    p.alpha_PopZPodJ=50*theta(12);%1;%PodJ on PopZ
    p.depol_popz=0.5*theta(13);%0.1; % deplymerization 
    % diffussion rates
    p.D_popzm =835;%750; % PopZ monomer dffusion  
p.D_popzp =0.0005;%40;%0.0005;  % PopZ polymer diffusion
   %% SpmX
   p.syn_spmx =6*theta(15);%0.04;
   p.deg_spmx = 0.1*theta(16);%0.1;
     p.dnv_spmx =1*theta(18);%0.1;%
   p.depol_spmx = 0.36*p.dnv_spmx;%theta(18);%
   p.aut_spmx=1*theta(18);
   p.Ja_SpmXCtrA=1;
   p.alpha_SpmXPopZ=50;
p.D_spmx =200;%um^2/min
p.D_spmxp=0.0005;
%% CtrA
p.syn_ctrA1=0.008;%0.0083;%Shenghua2008-exp; 0.026;%Murray-fitted
p.Ji_CtrACtrA=3;
p.syn_ctrA2=0.073;%0.073;%Shenghua2008-exp
p.Ja_CtrACtrA=3;
p.deg_ctrA1=0.0038;%0.0038;%exp-calculated
p.deg_ctrA2=0.15;%0.1 exp-estimated
p.Jd_CtrA=0.2;
p.dephoCtrA=100000*log(2)/8;%Murray-estimated~0.14
p.phoCtrA=100000*0.03;%7/93*p.dephoCtrA/0.5%~0.021;%Murray-estimated
p.D_CtrA=427;
p.D_CtrAP=p.D_CtrA;
%% PleC
p.syn_pleC=0.01;%0.053%Bronson
p.deg_pleC=0.02;%0.028;%Bronson
p.fb_PleC=10*0.1;
p.alpha_PleCPodJ=10;%0.5;
p.bf_PleC=10*0.05;
p.D_PleC=71;
%% DivJ
p.syn_divJ=0.008*2;%0.005;%Bronson
p.deg_divJ=0.035*2;%0.035;%Bronson
p.fb_DivJ=1;
p.alpha_DivJSpmX=5;
p.bf_DivJ=0.5;
p.D_DivJ=108;
%% DivK
p.syn_divK1=0.0004;%0.0004;%Bronson
p.syn_divK2=0.05*3;%0.125;%Bronson
p.Ja_DivKCtrA=0.7;%1.7/2;%Bronson
p.deg_divK=0.014*3;%Bronson; 0.002;%Shenghua2008-exp
p.phoDivK=0.01;
p.alpha_DivKDivJ=50;
p.dephoDivK=1;
p.D_DivK=1319;
p.D_DivKP=p.D_DivK;%p.D_DivK;

p.fb_DivKP=10000;
p.bf_DivKP=0.1;
%% PerP
p.syn_perP=0.5;%0.0428;%Bronson
p.Ja_PerPCtrA=1/2;%2.8/2;%Bronon
p.deg_perP=0.04;%0.04;%Bronson
p.D_PerP=853;
elseif ver==5
       p.D_podJm =100;%for PopZ only & 100 in manuscript;  
p.D_podJL =0.0005;%40;%0.0005;  

p.syn_podJ = 0.01;%0.2*theta(1);   % monomer synthesis 0.2
p.syn_podJ2 = 0.006;
p.Ji_PodJCtrA = 0.5;
p.n_PodJCtrA = 4;
% p.deg_podJm =theta(2);% 0.05;     0.046 paper116
% p.deg_podJp =theta(2);% 0.05;       
p.deg_podJm =0.3;% 0.05;    0.01?- 0.046 paper116
p.deg_podJp =p.deg_podJm;% 0.05;    
p.deg_podJ1 =0.007;%0.007;%0.007;%paper116

p.dnv_podJ =0.1*theta(3);% 15              % denovo polymerization 12 in dissertation for table C.2; 1 for table C.1
p.aut1_podJ =30*theta(4);%10;               % autocatalytic polymerization (pole)
p.aut1_podJ1=0*0.1*p.aut1_podJ;%central compartment autcatalytic polymerization
% p.depol_basal=0.2;
p.depol_podJ =0.2*theta(5);%0.5;             % deplymerization 

p.podj=2;
p.deg_s=0.5*theta(6);%0.05  %0.01 paper116ref
   p.alpha_PodJSpmX =100*theta(17);%2; %SpmX on PodJ dnv
%% PopZ
p.syn_popz=0.1*theta(8);%1.5;% popZ monomer synthesis
    p.deg_popzm=0.01;%theta(9);%0.05 % popZ monomer degradation
    p.deg_popzp=p.deg_popzm;%theta(9);%0.05;
    p.dnv_popz=0.3*theta(10);%60;% denovo polymerization
    p.aut_popz=0.6;%3*theta(11);%3;% autocatalytic polymerization
     p.aut_popz1=0.7*p.aut_popz;
    p.alpha_PopZPodJ=60*theta(12);%1;%PodJ on PopZ
    p.depol_popz=0.5*theta(13);%0.1; % deplymerization 
    % diffussion rates
    p.D_popzm =835;%750; % PopZ monomer dffusion  
p.D_popzp =0.0005;%40;%0.0005;  % PopZ polymer diffusion
   %% SpmX
   p.syn_spmx =6*theta(15);%0.04;
   p.deg_spmx = 0.1*theta(16);%0.1;
     p.dnv_spmx =0.8*theta(18);%0.1;%
   p.depol_spmx =p.dnv_spmx;%theta(18);%
   p.aut_spmx=1.2*theta(18);
   p.Ja_SpmXCtrA=1;
   p.alpha_SpmXPopZ=200;
p.D_spmx =200;%um^2/min
p.D_spmxp=0.0005;
%% CtrA
p.syn_ctrA1=0.008;%0.0083;%Shenghua2008-exp; 0.026;%Murray-fitted
p.Ji_CtrACtrA=3;
p.syn_ctrA2=0.073;%0.073;%Shenghua2008-exp
p.Ja_CtrACtrA=3;
p.deg_ctrA1=0.0038;%0.0038;%exp-calculated
p.deg_ctrA2=0.15;%0.1 exp-estimated
p.Jd_CtrA=0.2;
p.dephoCtrA=100000*log(2)/8;%Murray-estimated~0.14
p.phoCtrA=100000*0.03;%7/93*p.dephoCtrA/0.5%~0.021;%Murray-estimated
p.D_CtrA=427;
p.D_CtrAP=p.D_CtrA;
%% PleC
p.syn_pleC=0.01;%0.053%Bronson
p.deg_pleC=0.02;%0.028;%Bronson
p.fb_PleC=10*0.1;
p.alpha_PleCPodJ=10;%0.5;
p.bf_PleC=10*0.05;
p.D_PleC=71;
%% DivJ
p.syn_divJ=0.008*2;%0.005;%Bronson
p.deg_divJ=0.035*2;%0.035;%Bronson
p.fb_DivJ=1;
p.alpha_DivJSpmX=5;
p.bf_DivJ=0.5;
p.D_DivJ=108;
%% DivK
p.syn_divK1=0.0004;%0.0004;%Bronson
p.syn_divK2=0.05*3;%0.125;%Bronson
p.Ja_DivKCtrA=0.7;%1.7/2;%Bronson
p.deg_divK=0.014*3;%Bronson; 0.002;%Shenghua2008-exp
p.phoDivK=0.01;
p.alpha_DivKDivJ=50;
p.dephoDivK=1;
p.D_DivK=1319;
p.D_DivKP=p.D_DivK;%p.D_DivK;

p.fb_DivKP=10000;
p.bf_DivKP=0.1;
%% PerP
p.syn_perP=0.5;%0.0428;%Bronson
p.Ja_PerPCtrA=1/2;%2.8/2;%Bronon
p.deg_perP=0.04;%0.04;%Bronson
p.D_PerP=853;
elseif ver==6
      p.D_podJm =100;%for PopZ only & 100 in manuscript;  
p.D_podJL =0.0005;%40;%0.0005;  

p.syn_podJ = 0.01;%0.2*theta(1);   % monomer synthesis 0.2
p.syn_podJ2 = 0.006;
p.Ji_PodJCtrA = 0.5;
p.n_PodJCtrA = 4;
% p.deg_podJm =theta(2);% 0.05;     0.046 paper116
% p.deg_podJp =theta(2);% 0.05;       
p.deg_podJm =0.3;% 0.05;    0.01?- 0.046 paper116
p.deg_podJp =p.deg_podJm;% 0.05;    
p.deg_podJ1 =0.007;%0.007;%0.007;%paper116

p.dnv_podJ =0.1*theta(3);% 15              % denovo polymerization 12 in dissertation for table C.2; 1 for table C.1
p.aut1_podJ =30*theta(4);%10;               % autocatalytic polymerization (pole)
p.aut1_podJ1=0*0.1*p.aut1_podJ;%central compartment autcatalytic polymerization
% p.depol_basal=0.2;
p.depol_podJ =0.2*theta(5);%0.5;             % deplymerization 

p.podj=2;
p.deg_s=0.5*theta(6);%0.05  %0.01 paper116ref
   p.alpha_PodJSpmX =100*theta(17);%2; %SpmX on PodJ dnv
%% PopZ
p.syn_popz=0.1*theta(8);%1.5;% popZ monomer synthesis
    p.deg_popzm=0.02;%theta(9);%0.05 % popZ monomer degradation
    p.deg_popzp=p.deg_popzm;%theta(9);%0.05;
    p.dnv_popz=0.5*theta(10);%60;% denovo polymerization
    p.aut_popz=3;%3*theta(11);%3;% autocatalytic polymerization
     p.aut_popz1=0.7*p.aut_popz;
    p.alpha_PopZPodJ=60*theta(12);%1;%PodJ on PopZ
    p.depol_popz=0.5*theta(13);%0.1; % deplymerization 
    % diffussion rates
    p.D_popzm =835;%750; % PopZ monomer dffusion  
p.D_popzp =0.0005;%40;%0.0005;  % PopZ polymer diffusion
   %% SpmX
   p.syn_spmx =6*theta(15);%0.04;
   p.deg_spmx = 0.1*theta(16);%0.1;
     p.dnv_spmx =0.8*theta(18);%0.1;%
   p.depol_spmx =1.6*p.dnv_spmx;%theta(18);%
   p.aut_spmx=1.2*theta(18);
   p.Ja_SpmXCtrA=1;
   p.alpha_SpmXPopZ=5;
p.D_spmx =200;%um^2/min
p.D_spmxp=0.0005;
%% CtrA
p.syn_ctrA1=0.008;%0.0083;%Shenghua2008-exp; 0.026;%Murray-fitted
p.Ji_CtrACtrA=3;
p.syn_ctrA2=0.073;%0.073;%Shenghua2008-exp
p.Ja_CtrACtrA=3;
p.deg_ctrA1=0.0038;%0.0038;%exp-calculated
p.deg_ctrA2=0.15;%0.1 exp-estimated
p.Jd_CtrA=0.2;
p.dephoCtrA=100000*log(2)/8;%Murray-estimated~0.14
p.phoCtrA=100000*0.03;%7/93*p.dephoCtrA/0.5%~0.021;%Murray-estimated
p.D_CtrA=427;
p.D_CtrAP=p.D_CtrA;
%% PleC
p.syn_pleC=0.01;%0.053%Bronson
p.deg_pleC=0.02;%0.028;%Bronson
p.fb_PleC=10*0.1;
p.alpha_PleCPodJ=10;%0.5;
p.bf_PleC=10*0.05;
p.D_PleC=71;
%% DivJ
p.syn_divJ=0.008*2;%0.005;%Bronson
p.deg_divJ=0.035*2;%0.035;%Bronson
p.fb_DivJ=1;
p.alpha_DivJSpmX=5;
p.bf_DivJ=0.5;
p.D_DivJ=108;
%% DivK
p.syn_divK1=0.0004;%0.0004;%Bronson
p.syn_divK2=0.05*3;%0.125;%Bronson
p.Ja_DivKCtrA=0.7;%1.7/2;%Bronson
p.deg_divK=0.014*3;%Bronson; 0.002;%Shenghua2008-exp
p.phoDivK=0.01;
p.alpha_DivKDivJ=50;
p.dephoDivK=1;
p.D_DivK=1319;
p.D_DivKP=p.D_DivK;%p.D_DivK;

p.fb_DivKP=10000;
p.bf_DivKP=0.1;
%% PerP
p.syn_perP=0.5;%0.0428;%Bronson
p.Ja_PerPCtrA=1/2;%2.8/2;%Bronon
p.deg_perP=0.04;%0.04;%Bronson
p.D_PerP=853;
elseif ver==7
     p.D_podJm =100;%for PopZ only & 100 in manuscript;  
p.D_podJL =0.0005;%40;%0.0005;  

p.syn_podJ = 0.01;%0.2*theta(1);   % monomer synthesis 0.2
p.syn_podJ2 = 0.006;
p.Ji_PodJCtrA = 0.5;
p.n_PodJCtrA = 4;
% p.deg_podJm =theta(2);% 0.05;     0.046 paper116
% p.deg_podJp =theta(2);% 0.05;       
p.deg_podJm =0.3;% 0.05;    0.01?- 0.046 paper116
p.deg_podJp =p.deg_podJm;% 0.05;    
p.deg_podJ1 =0.007;%0.007;%0.007;%paper116

p.dnv_podJ =0.1*theta(3);% 15              % denovo polymerization 12 in dissertation for table C.2; 1 for table C.1
p.aut1_podJ =30*theta(4);%10;               % autocatalytic polymerization (pole)
p.aut1_podJ1=0*0.1*p.aut1_podJ;%central compartment autcatalytic polymerization
% p.depol_basal=0.2;
p.depol_podJ =0.2*theta(5);%0.5;             % deplymerization 

p.podj=2;
p.deg_s=0.5*theta(6);%0.05  %0.01 paper116ref
   p.alpha_PodJSpmX =100*theta(17);%2; %SpmX on PodJ dnv
%% PopZ
p.syn_popz=0.1*theta(8);%1.5;% popZ monomer synthesis
    p.deg_popzm=0.02;%theta(9);%0.05 % popZ monomer degradation
    p.deg_popzp=p.deg_popzm;%theta(9);%0.05;
    p.dnv_popz=0.5*theta(10);%60;% denovo polymerization
    p.aut_popz=5;%3*theta(11);%3;% autocatalytic polymerization
     p.aut_popz1=0.7*p.aut_popz;
    p.alpha_PopZPodJ=60*theta(12);%1;%PodJ on PopZ
    p.depol_popz=0.5*theta(13);%0.1; % deplymerization 
    % diffussion rates
    p.D_popzm =835;%750; % PopZ monomer dffusion  
p.D_popzp =0.0005;%40;%0.0005;  % PopZ polymer diffusion
   %% SpmX
   p.syn_spmx =6*theta(15);%0.04;
   p.deg_spmx = 0.1*theta(16);%0.1;
     p.dnv_spmx =0.8*theta(18);%0.1;%
   p.depol_spmx =1.6*p.dnv_spmx;%theta(18);%
   p.aut_spmx=1.2*theta(18);
   p.Ja_SpmXCtrA=1;
   p.alpha_SpmXPopZ=5;
p.D_spmx =200;%um^2/min
p.D_spmxp=0.0005;
%% CtrA
p.syn_ctrA1=0.008;%0.0083;%Shenghua2008-exp; 0.026;%Murray-fitted
p.Ji_CtrACtrA=3;
p.syn_ctrA2=0.073;%0.073;%Shenghua2008-exp
p.Ja_CtrACtrA=3;
p.deg_ctrA1=0.0038;%0.0038;%exp-calculated
p.deg_ctrA2=0.15;%0.1 exp-estimated
p.Jd_CtrA=0.2;
p.dephoCtrA=100000*log(2)/8;%Murray-estimated~0.14
p.phoCtrA=100000*0.03;%7/93*p.dephoCtrA/0.5%~0.021;%Murray-estimated
p.D_CtrA=427;
p.D_CtrAP=p.D_CtrA;
%% PleC
p.syn_pleC=0.01;%0.053%Bronson
p.deg_pleC=0.02;%0.028;%Bronson
p.fb_PleC=10*0.1;
p.alpha_PleCPodJ=10;%0.5;
p.bf_PleC=10*0.05;
p.D_PleC=71;
%% DivJ
p.syn_divJ=0.008*2;%0.005;%Bronson
p.deg_divJ=0.035*2;%0.035;%Bronson
p.fb_DivJ=1;
p.alpha_DivJSpmX=5;
p.bf_DivJ=0.5;
p.D_DivJ=108;
%% DivK
p.syn_divK1=0.0004;%0.0004;%Bronson
p.syn_divK2=0.05*3;%0.125;%Bronson
p.Ja_DivKCtrA=0.7;%1.7/2;%Bronson
p.deg_divK=0.014*3;%Bronson; 0.002;%Shenghua2008-exp
p.phoDivK=0.01;
p.alpha_DivKDivJ=50;
p.dephoDivK=1;
p.D_DivK=1319;
p.D_DivKP=p.D_DivK;%p.D_DivK;

p.fb_DivKP=10000;
p.bf_DivKP=0.1;
%% PerP
p.syn_perP=0.5;%0.0428;%Bronson
p.Ja_PerPCtrA=1/2;%2.8/2;%Bronon
p.deg_perP=0.04;%0.04;%Bronson
p.D_PerP=853;
end
if strcmp(mutant,'deltaPodJ')
    p.syn_podJ = 0;
p.syn_podJ2 = 0;
elseif strcmp(mutant,'PodJ+')
    p.syn_podJ = p.syn_podJ*2;
p.syn_podJ2 = p.syn_podJ2*2;
elseif strcmp(mutant,'deltaSpmX')
    p.syn_spmx=0;
elseif strcmp(mutant,'deltaPopZ')
    p.syn_popz=0;
end
end
