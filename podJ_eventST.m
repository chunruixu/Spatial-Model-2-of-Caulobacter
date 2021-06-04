
function [value,isterminal,direction] = podJ_eventST(t,y)

global T_e1
% T_e1=4;
CtrAP=sum(y(33:36));
% T_e1=34-30;%30-25;%T_ini;
% T_term=115-34;%120;
% T_e2=(T_term-T_e1)*0.37+T_e1;%S_CtrA changed from 0 to 1; the fork passes ctrA
% T_e3=(T_term-T_e1)*0.65+T_e1;%S_pleC changed from 0 to 1; the fork passes pleC
% T_e4=(T_term-T_e1)*0.74+T_e1;%S_perP changed from 0 to 1; the fork passes perP
% T_e5=(T_term-T_e1)*0.87+T_e1;%S_podJ changed from 0 to 1; the fork passes podJ
T1=max((34-T_e1),0);
T_Sphase=90;
T_term=min((T1+T_Sphase),125-34);
T_e2=(T_term-T1)*0.37+T1;%S_CtrA changed from 0 to 1; the fork passes ctrA
T_e3=(T_term-T1)*0.65+T1;%S_pleC changed from 0 to 1; the fork passes pleC
T_e4=(T_term-T1)*0.74+T1;%S_perP changed from 0 to 1; the fork passes perP
T_e5=(T_term-T1)*0.87+T1;%S_podJ changed from 0 to 1; the fork passes podJ


value = [ sign(0.25-CtrAP);sign(t - T_e2); sign(t - T_e3); sign(t - T_e4);sign(t - T_e5);sign(t-T_term)];
isterminal = [1; 1; 1;1;1;1];
direction = [+1; +1; +1;+1;+1;+1];

end