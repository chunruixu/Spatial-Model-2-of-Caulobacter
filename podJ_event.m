% function [value,isterminal,direction] = podJ_event(t,y)
% %% for main_PopZ1; Ncf: 75-125
% PerP_end = 30;
% SpmX_start = 34;
% 
% value = [sign(t - PerP_end); sign(t-SpmX_start)];
% isterminal = [1; 1];
% direction = [+1; +1];
% 
% 
% 
% end

function [value,isterminal,direction] = podJ_event(t,y)
global T_e1
%% for sw-pd and pd-end cell simulation (0-120; 120:150min)
CtrAP=sum(y(33:36));
% T_e1=30.6;%T_ini=22+-3;
% T_term=115;%120;
% T_e2=(T_term-T_e1)*0.37+T_e1;%S_CtrA changed from 0 to 1; the fork passes ctrA
% T_e3=(T_term-T_e1)*0.65+T_e1;%S_pleC changed from 0 to 1; the fork passes pleC
% T_e4=(T_term-T_e1)*0.74+T_e1;%S_perP changed from 0 to 1; the fork passes perP
% T_e5=(T_term-T_e1)*0.87+T_e1;%S_podJ changed from 0 to 1; the fork passes podJ
T_Sphase=90;
T_term=min((T_e1+T_Sphase),125);
T_e2=(T_term-T_e1)*0.37+T_e1;%S_CtrA changed from 0 to 1; the fork passes ctrA
T_e3=(T_term-T_e1)*0.65+T_e1;%S_pleC changed from 0 to 1; the fork passes pleC
T_e4=(T_term-T_e1)*0.74+T_e1;%S_perP changed from 0 to 1; the fork passes perP
T_e5=(T_term-T_e1)*0.87+T_e1;%S_podJ changed from 0 to 1; the fork passes podJ

% T_e2=T_Sphase*0.37+T_e1;%S_CtrA changed from 0 to 1; the fork passes ctrA
% T_e3=T_Sphase*0.65+T_e1;%S_pleC changed from 0 to 1; the fork passes pleC
% T_e4=T_Sphase*0.74+T_e1;%S_perP changed from 0 to 1; the fork passes perP
% T_e5=T_Sphase*0.87+T_e1;%S_podJ changed from 0 to 1; the fork passes podJ
% 
%max*0.2;
value = [sign(0.25-CtrAP);sign(t - T_e2); sign(t - T_e3); sign(t - T_e4);sign(t - T_e5); sign(t-T_term)];
isterminal = [1; 1; 1;1;1;1];
direction = [+1; +1; +1;+1;+1;+1];
end