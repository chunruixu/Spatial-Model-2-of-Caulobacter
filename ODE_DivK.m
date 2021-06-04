%6/2/2021
function dydt = ODE_DivK(t,y)

global p;
   
%% Equations

dydt=zeros(107,1);
m=0.01;%0.15;
mplec=0.15;%0.15;
%% PodJL_m
m_podj=((1-m)*y(103)+m);
% Bin 1
dydt(1) = 0*p.syn_podJ*m_podj+0*p.syn_podJ2*p.Ji_PodJCtrA^p.n_PodJCtrA/(p.Ji_PodJCtrA^p.n_PodJCtrA+y(1+48)^p.n_PodJCtrA)...
    - (p.deg_podJ1+p.deg_podJm*y(1+96))*y(1)... % synthesis & degradation*PerP (97-102); 
    - p.dnv_podJ*y(1)...% denovo polymerization; SpmX (19-24&25-30)
    - p.aut1_podJ/(1+p.alpha_PodJSpmX*(y(1+18)+y(1+24)))*y(1+6)^p.podj*y(1)...% autocatalytic polymerization at poles
+ p.depol_podJ*y(1+6)...%
    + p.D_podJm*(y(2)-y(1))/(y(107)^2)-p.mu*y(1); %
% Bin 2 - 5
for v1=2:5
    dydt(v1) = p.syn_podJ*m_podj+p.syn_podJ2*p.Ji_PodJCtrA^p.n_PodJCtrA/(p.Ji_PodJCtrA^p.n_PodJCtrA+y(v1+48)^p.n_PodJCtrA)...%CtrAP 49-54
        - (p.deg_podJ1+p.deg_podJm*y(v1+96))*y(v1)... % synthesis & degradation*PerP (65-68); y(61)=SpodJ 
    - p.dnv_podJ*y(v1)...% denovo polymerization; SpmX (13-16&17-20)
    - p.aut1_podJ1/(1+p.alpha_PodJSpmX*(y(v1+18)+y(v1+24)))*y(v1+6)^p.podj*y(v1)...% autocatalytic polymerization at centers
+ p.depol_podJ*y(v1+6)...%
    +p.D_podJm*(y(v1-1)-2*y(v1)+y(v1+1))/(y(107)^2)- p.mu*y(v1);
end
% Bin 6
    dydt(6) = 0*p.syn_podJ*m_podj+0*p.syn_podJ2*p.Ji_PodJCtrA^p.n_PodJCtrA/(p.Ji_PodJCtrA^p.n_PodJCtrA+y(6+48)^p.n_PodJCtrA)...%CtrAP 49-54
        - (p.deg_podJ1+p.deg_podJm*y(6+96))*y(6)... % synthesis & degradation*PerP (65-68); y(61)=SpodJ 
    - p.dnv_podJ*y(6)...% denovo polymerization; SpmX (13-16&17-20)
    - p.aut1_podJ/(1+p.alpha_PodJSpmX*(y(6+18)+y(6+24)))*y(6+6)^p.podj*y(6)...% autocatalytic polymerization at centers
+ p.depol_podJ*y(6+6)...%
    +p.D_podJm*(y(5)-y(6))/(y(107)^2)- p.mu*y(6);

%% PodJL_p
% Bin 7
dydt(7) =-(p.deg_podJ1+p.deg_podJp*y(7+90))*y(7)...%PerP (97-102)
    + p.dnv_podJ*y(7-6)...%SpmX (19-24&25-30)
    + p.aut1_podJ/(1+p.alpha_PodJSpmX*(y(7+12)+y(7+18)))*y(7)^p.podj*y(7-6)...% autocatalytic polymerization at poles
    - p.depol_podJ*y(7)...%
    + p.D_podJL*(y(8)-y(7))/(y(107)^2)- p.mu*y(7); % polymer diffusion % dilution
% Bin 8-11
for v2=8:11
    dydt(v2) =-(p.deg_podJ1+p.deg_podJp*y(v2+90))*y(v2)...%PerP (97-102)
    + p.dnv_podJ*y(v2-6)...%SpmX (19-24&25-30)
    + p.aut1_podJ1/(1+p.alpha_PodJSpmX*(y(v2+12)+y(v2+18)))*y(v2)^p.podj*y(v2-6)...% autocatalytic polymerization at poles
    - p.depol_podJ*y(v2)...%
    +p.D_podJL*(y(v2-1)-2*y(v2)+y(v2+1))/(y(107)^2)- p.mu*y(v2); % polymer diffusion % dilution
end
%Bin 12
    dydt(12) =-(p.deg_podJ1+p.deg_podJp*y(12+90))*y(12)...%PerP (97-102)
    + p.dnv_podJ*y(12-6)...%SpmX (19-24&25-30)
    + p.aut1_podJ/(1+p.alpha_PodJSpmX*(y(12+12)+y(12+18)))*y(12)^p.podj*y(12-6)...% autocatalytic polymerization at poles
    - p.depol_podJ*y(12)...%
    +p.D_podJL*(y(11)-y(12))/(y(107)^2)- p.mu*y(12); % polymer diffusion % dilution

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PodJS
% no diffussion assumed
for i=13:18
dydt(i)=(p.deg_podJ1+p.deg_podJm*y(i+84))*(y(i-12)+y(i-6))...%%PerP (97-102)
    -(p.mu+p.deg_s)*y(i);
end

%% SpmXm
%Bin 19
q=1;
dydt(19) = 0*p.syn_spmx*y(19+30)^q/(y(19+30)^q+p.Ja_SpmXCtrA^q)-(p.deg_spmx+p.mu)*y(19)...%CtrAP 49-54
    -p.dnv_spmx*y(19)+p.depol_spmx*y(19+6)...
    -p.aut_spmx*(1+p.alpha_SpmXPopZ*(y(19+12)+y(19+18)))*y(19)*y(19+6)^2 ...%PopZ: 31-36&37-42
+p.D_spmx*(y(20)-y(19))/(y(107)^2);

%Bin 20-23
for v3=20:23
dydt(v3) = p.syn_spmx*y(v3+30)^q/(y(v3+30)^q+p.Ja_SpmXCtrA^q)-(p.deg_spmx+p.mu)*y(v3)...%CtrAP 33-36 
    -p.dnv_spmx*y(v3)+p.depol_spmx*y(v3+6)...
  -p.aut_spmx*(1+p.alpha_SpmXPopZ*(y(v3+12)+y(v3+18)))*y(v3)*y(v3+6)^2 ...%PopZp 25-28
    +p.D_spmx*(y(v3-1)-2*y(v3)+y(v3+1))/(y(107)^2);
end
%Bin 24
dydt(24) = 0*p.syn_spmx*y(24+30)^q/(y(24+30)^q+p.Ja_SpmXCtrA^q)-(p.deg_spmx+p.mu)*y(24)...%CtrAP 49-54
    -p.dnv_spmx*y(24)+p.depol_spmx*y(24+6)...
    -p.aut_spmx*(1+p.alpha_SpmXPopZ*(y(24+12)+y(24+18)))*y(24)*y(24+6)^2 ...%PopZ: 31-36&37-42
+p.D_spmx*(y(23)-y(24))/(y(107)^2);


%% SpmXp
%Bin 25
dydt(25)=-(p.deg_spmx+p.mu)*y(25)...
    +p.dnv_spmx*y(25-6)-p.depol_spmx*y(25)...
      +p.aut_spmx*(1+p.alpha_SpmXPopZ*(y(25+6)+y(25+12)))*y(25-6)*y(25)^2 ...%PopZ: 31-36&37-42
     +p.D_spmxp*(y(26)-y(25))/(y(107)^2);
 %Bin 26-29
 for v4=26:29
     dydt(v4)=-(p.deg_spmx+p.mu)*y(v4)...
    +p.dnv_spmx*y(v4-6)-p.depol_spmx*y(v4)...
      +p.aut_spmx*(1+p.alpha_SpmXPopZ*(y(v4+6)+y(v4+12)))*y(v4-6)*y(v4)^2 ...%PopZp 25-28
      +p.D_spmxp*(y(v4-1)-2*y(v4)+y(v4+1))/(y(107)^2);
 end
%Bin 30
dydt(30)=-(p.deg_spmx+p.mu)*y(30)...
    +p.dnv_spmx*y(30-6)-p.depol_spmx*y(30)...
      +p.aut_spmx*(1+p.alpha_SpmXPopZ*(y(30+6)+y(30+12)))*y(30-6)*y(30)^2 ...%PopZ: 31-36&37-42
     +p.D_spmxp*(y(29)-y(30))/(y(107)^2);


%% PopZm
%Bin 31
dydt(31) =0*p.syn_popz  - (p.mu+p.deg_popzm)*y(31)...
    -p.dnv_popz*y(31)-p.aut_popz*(1+p.alpha_PopZPodJ*(y(31-30)+y(31-24)))*y(31+6)^2*y(31)...%PodJp 7-12
    +p.depol_popz*y(31+6)...
    +p.D_popzm*(y(32)-y(31))/(y(107)^2);
%Bin 32-35
for v5=32:35
    dydt(v5) =p.syn_popz  - (p.mu+p.deg_popzm)*y(v5)...
    -p.dnv_popz*y(v5)-p.aut_popz1*(1+p.alpha_PopZPodJ*(y(v5-30)+y(v5-24)))*y(v5+6)^2*y(v5)...%PodJp 7-12
    +p.depol_popz*y(v5+6)...
+p.D_popzm*(y(v5-1)-2*y(v5)+y(v5+1))/(y(107)^2);
end
%Bin 36
dydt(36) =0*p.syn_popz  - (p.mu+p.deg_popzm)*y(36)...
    -p.dnv_popz*y(36)-p.aut_popz*(1+p.alpha_PopZPodJ*(y(36-30)+y(36-24)))*y(36+6)^2*y(36)...%PodJp 7-12
    +p.depol_popz*y(36+6)...
    +p.D_popzm*(y(35)-y(36))/(y(107)^2);


%% PopZp
%Bin 37
dydt(37) = -(p.mu+p.deg_popzp)*y(37)...
    +p.dnv_popz*y(37-6)+p.aut_popz*(1+p.alpha_PopZPodJ*(y(37-36)+y(37-30)))*y(37)^2*y(37-6)...%PodJp 7-12
    -p.depol_popz*y(37)...
 +p.D_popzp*(y(38)-y(37))/(y(107)^2);
%Bin 38-41
for v6=38:41
    dydt(v6) = -(p.mu+p.deg_popzp)*y(v6)...
    +p.dnv_popz*y(v6-6)+p.aut_popz1*(1+p.alpha_PopZPodJ*(y(v6-36)+y(v6-30)))*y(v6)^2*y(v6-6)...%PodJp 7-12
    -p.depol_popz*y(v6)...
   +p.D_popzp*(y(v6-1)-2*y(v6)+y(v6+1))/(y(107)^2);
end
%Bin 42
dydt(42) = -(p.mu+p.deg_popzp)*y(42)...
    +p.dnv_popz*y(42-6)+p.aut_popz*(1+p.alpha_PopZPodJ*(y(42-36)+y(42-30)))*y(42)^2*y(42-6)...%PodJp 7-12
    -p.depol_popz*y(42)...
 +p.D_popzp*(y(41)-y(42))/(y(107)^2);


%% CtrA
m_ctrA = (1-m)*y(104)+m; %S_ctrA=y(104)
%Bin 43
dydt(43)=0*p.syn_ctrA1*(1-y(43+6)/(y(43+6)+y(43)+p.Ji_CtrACtrA))*m_ctrA...
    +0*p.syn_ctrA2*y(43+6)/(y(43+6)+y(43)+p.Ja_CtrACtrA) ...
    -(p.mu+p.deg_ctrA1+p.deg_ctrA2*(y(43+42)+y(43+48))^2/(p.Jd_CtrA^2+(y(43+42)+y(43+48))^2))*y(43) ...%DivKPT (85~90)+(91~96)
    -p.phoCtrA*y(43)+p.dephoCtrA*(y(43+42)+y(43+48))*y(43+6) ...
+p.D_CtrA*(y(44)-y(43))/(y(107)^2);
%Bin 44-47
for v7=44:47
    dydt(v7)=p.syn_ctrA1*(1-y(v7+6)/(y(v7+6)+y(v7)+p.Ji_CtrACtrA))*m_ctrA...
    +p.syn_ctrA2*y(v7+6)/(y(v7+6)+y(v7)+p.Ja_CtrACtrA) ...
    -(p.mu+p.deg_ctrA1+p.deg_ctrA2*(y(v7+42)+y(v7+48))^2/(p.Jd_CtrA^2+(y(v7+42)+y(v7+48))^2))*y(v7) ...%DivKPT (85~90)+(91~96)
    -p.phoCtrA*y(v7)+p.dephoCtrA*(y(v7+42)+y(v7+48))*y(v7+6) ...
+p.D_CtrA*(y(v7-1)-2*y(v7)+y(v7+1))/(y(107)^2);
end
%Bin 48
dydt(48)=0*p.syn_ctrA1*(1-y(48+6)/(y(48+6)+y(48)+p.Ji_CtrACtrA))*m_ctrA...
    +0*p.syn_ctrA2*y(48+6)/(y(48+6)+y(48)+p.Ja_CtrACtrA) ...
    -(p.mu+p.deg_ctrA1+p.deg_ctrA2*(y(48+42)+y(48+48))^2/(p.Jd_CtrA^2+(y(48+42)+y(48+48))^2))*y(48) ...%DivKPT (85~90)+(91~96)
    -p.phoCtrA*y(48)+p.dephoCtrA*(y(48+42)+y(48+48))*y(48+6) ...
+p.D_CtrA*(y(47)-y(48))/(y(107)^2);


%% CtrAP
%Bin 49
 dydt(49)=-(p.mu+p.deg_ctrA1+p.deg_ctrA2*(y(49+36)+y(49+42))^2/(p.Jd_CtrA^2+(y(49+36)+y(49+42))^2))*y(49) ...%DivKPT (85~90)+(91~96)
    +p.phoCtrA*y(49-6)-p.dephoCtrA*(y(49+36)+y(49+42))*y(49) ...
    + p.D_CtrAP*(y(50)-y(49))/(y(107)^2);
%Bin 50-53
for v8=50:53
    dydt(v8)=-(p.mu+p.deg_ctrA1+p.deg_ctrA2*(y(v8+36)+y(v8+42))^2/(p.Jd_CtrA^2+(y(v8+36)+y(v8+42))^2))*y(v8) ...%
    +p.phoCtrA*y(v8-6)-p.dephoCtrA*(y(v8+36)+y(v8+42))*y(v8) ...
    +p.D_CtrAP*(y(v8-1)-2*y(v8)+y(v8+1))/(y(107)^2);
end
%Bin 54
 dydt(54)=-(p.mu+p.deg_ctrA1+p.deg_ctrA2*(y(54+36)+y(54+42))^2/(p.Jd_CtrA^2+(y(54+36)+y(54+42))^2))*y(54) ...%DivKPT (85~90)+(91~96)
    +p.phoCtrA*y(54-6)-p.dephoCtrA*(y(54+36)+y(54+42))*y(54) ...
    + p.D_CtrAP*(y(53)-y(54))/(y(107)^2);

%% PleCf
m_pleC = (1-mplec)*y(105)+mplec;%S_pleC=y(105)
%Bin 55
dydt(55)=0*p.syn_pleC*m_pleC -(p.mu+p.deg_pleC)*y(55) ...
    -p.fb_PleC*(1+p.alpha_PleCPodJ*(y(55-54)+y(55-48)))*y(55)+p.bf_PleC*y(55+6) ...%PodJp 7-12
    + p.D_PleC*(y(56)-y(55))/(y(107)^2);
%Bin 56-59
for v9=56:59
    dydt(v9)=p.syn_pleC*m_pleC -(p.mu+p.deg_pleC)*y(v9) ...
    -p.fb_PleC*(1+p.alpha_PleCPodJ*(y(v9-54)+y(v9-48)))*y(v9)+p.bf_PleC*y(v9+6) ...%
    + p.D_PleC*(y(v9-1)-2*y(v9)+y(v9+1))/(y(107)^2);
end
%Bin 60
dydt(60)=0*p.syn_pleC*m_pleC -(p.mu+p.deg_pleC)*y(60) ...
    -p.fb_PleC*(1+p.alpha_PleCPodJ*(y(60-54)+y(60-48)))*y(60)+p.bf_PleC*y(60+6) ...%PodJp 7-12
    + p.D_PleC*(y(59)-y(60))/(y(107)^2);


%% PleCb
for i=61:66
    dydt(i)=-(p.mu+p.deg_pleC)*y(i) ...
        +p.fb_PleC*(1+p.alpha_PleCPodJ*(y(i-60)+y(i-54)))*y(i-6)-p.bf_PleC*y(i) ;%PodJp 7-12
end

%% DivJf
%Bin 67
dydt(67)=0*p.syn_divJ -(p.mu+p.deg_divJ)*y(67) ...
    -p.fb_DivJ*(1+p.alpha_DivJSpmX*(y(67-48)+y(67-42)))*y(67)+p.bf_DivJ*y(67+6) ...%SpmXp 19-24&25-30
    + p.D_DivJ*(y(68)-y(67))/(y(107)^2);
%Bin 68-71
for v10=68:71
    dydt(v10)=p.syn_divJ -(p.mu+p.deg_divJ)*y(v10) ...
    -p.fb_DivJ*(1+p.alpha_DivJSpmX*(y(v10-48)+y(v10-42)))*y(v10)+p.bf_DivJ*y(v10+6) ...%
    +p.D_DivJ*(y(v10-1)-2*y(v10)+y(v10+1))/(y(107)^2);
end
%Bin 72
dydt(72)=0*p.syn_divJ -(p.mu+p.deg_divJ)*y(72) ...
    -p.fb_DivJ*(1+p.alpha_DivJSpmX*(y(72-48)+y(72-42)))*y(72)+p.bf_DivJ*y(72+6) ...%SpmXp 19-24&25-30
    + p.D_DivJ*(y(71)-y(72))/(y(107)^2);

%% DivJb
for i=73:78
    dydt(i)=-(p.mu+p.deg_divJ)*y(i) ...
        +p.fb_DivJ*(1+p.alpha_DivJSpmX*(y(i-54)+y(i-48)))*y(i-6)-p.bf_DivJ*y(i) ;%SpmXp 19-24&25-30
end
%% DivK
%Bin 79
dydt(79)=0*p.syn_divK1+0*p.syn_divK2*y(79-30)^2/(y(79-30)^2+p.Ja_DivKCtrA^2) ...%CtrAP 49-54
    -(p.mu+p.deg_divK)*y(79) ...
    -p.phoDivK*(p.alpha_DivKDivJ*y(79-6)+1)*y(79)+p.dephoDivK*(y(79-24)+y(79-18))*(y(79+6)+y(79+12)) ...%DivJb 73-78; PleC 55-60&61-66
+p.D_DivK*(y(80)-y(79))/(y(107)^2);
%Bin 80-83
for v11=80:83
    dydt(v11)=p.syn_divK1+p.syn_divK2*y(v11-30)^2/(y(v11-30)^2+p.Ja_DivKCtrA^2) ...%CtrAP 
    -(p.mu+p.deg_divK)*y(v11) ...
    -p.phoDivK*(p.alpha_DivKDivJ*y(v11-6)+1)*y(v11)+p.dephoDivK*(y(v11-24)+y(v11-18))*(y(v11+6)+y(v11+12)) ...%
+p.D_DivK*(y(v11-1)-2*y(v11)+y(v11+1))/(y(107)^2);
end
%Bin 84
dydt(84)=0*p.syn_divK1+0*p.syn_divK2*y(84-30)^2/(y(84-30)^2+p.Ja_DivKCtrA^2) ...%CtrAP 49-54
    -(p.mu+p.deg_divK)*y(84) ...
    -p.phoDivK*(p.alpha_DivKDivJ*y(84-6)+1)*y(84)+p.dephoDivK*(y(84-24)+y(84-18))*(y(84+6)+y(84+12)) ...%DivJb 73-78; PleC 55-60&61-66
+p.D_DivK*(y(83)-y(84))/(y(107)^2);

%% DivKPf
%Bin 85
dydt(85)= -(p.mu+p.deg_divK)*y(85) ...
    +p.phoDivK*(p.alpha_DivKDivJ*y(85-12)+1)*y(85-6)-p.dephoDivK*(y(85-30)+y(85-24))*y(85) ...%DivJb 73-78; PleC 55-60&61-66
-p.fb_DivKP*y(85)+p.bf_DivKP*y(85+6)+p.D_DivKP*(y(86)-y(85))/(y(107)^2);
%Bin 86-89
for v12=86:89
    dydt(v12)=(p.mu+p.deg_divK)*y(v12) ...
    +p.phoDivK*(p.alpha_DivKDivJ*y(v12-12)+1)*y(v12-6)-p.dephoDivK*(y(v12-30)+y(v12-34))*y(v12) ...%DivJ 45-48&49-52; PleC 37-40&41-44
-p.fb_DivKP*y(v12)+p.bf_DivKP*y(v12+6)+p.D_DivKP*(y(v12-1)-2*y(v12)+y(v12+1))/(y(107)^2);
end
%Bin 90
dydt(90)= -(p.mu+p.deg_divK)*y(90) ...
    +p.phoDivK*(p.alpha_DivKDivJ*y(90-12)+1)*y(90-6)-p.dephoDivK*(y(90-30)+y(90-24))*y(90) ...%DivJb 73-78; PleC 55-60&61-66
-p.fb_DivKP*y(90)+p.bf_DivKP*y(90+6)+p.D_DivKP*(y(89)-y(90))/(y(107)^2);

%% DivKPb
for i=91:96
dydt(i)=-(p.mu+p.deg_divK)*y(i)+p.fb_DivKP*y(i-6)-p.bf_DivKP*y(i)-p.dephoDivK*(y(i-36)+y(i-30))*y(i);
end
%% PerP
m_perP=((1-m)*y(106)+m);
%Bin 97
dydt(97)=0*p.syn_perP*m_perP*y(97-48)^2/(y(97-48)^2+p.Ja_PerPCtrA^2) ...%CtrAP 49-54
    -(p.mu+p.deg_perP)*y(97) ...
    +p.D_PerP*(y(98)-y(97))/(y(107)^2);
%Bin 98-101
for v13=98:101
    dydt(v13)=p.syn_perP*m_perP*y(v13-48)^2/(y(v13-48)^2+p.Ja_PerPCtrA^2) ...%CtrAP
    -(p.mu+p.deg_perP)*y(v13) ...
  +p.D_PerP*(y(v13-1)-2*y(v13)+y(v13+1))/(y(107)^2);
end
%Bin 102
dydt(102)=0*p.syn_perP*m_perP*y(102-48)^2/(y(102-48)^2+p.Ja_PerPCtrA^2) ...%CtrAP 49-54
    -(p.mu+p.deg_perP)*y(102) ...
    +p.D_PerP*(y(101)-y(102))/(y(107)^2);

%% S
for i=103:106
% dydt(i)=0;
dydt(i)=0;
end%S - PodJ, CtrA, PleC, PerP
%% cell growth equation
dydt(107)=p.mu*y(107);








