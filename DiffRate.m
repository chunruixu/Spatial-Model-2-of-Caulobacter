function [D, D_] = DiffRate(M)
%KDa
alpha=4.3*10^3;
D0=0.65;%um^2/s
D_=alpha/M^2+D0;
D=D_*60;%um^2/min
end