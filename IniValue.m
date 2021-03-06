function y0=IniValue(yout,celltype)
%6 compartment %even bins
Y0=yout(:,end);
% L=Y0(107)+2*Y0(108);%length of half-cell at the end time point
% a=Y0(107)-0.2*L; b=Y0(108)-0.5*L;
if strcmp(celltype,'SW')
for i=1:17
y0_(i*6-6+1)=Y0(i*6-6+3);
y0_(i*6-6+2)=Y0(i*6-6+3);
y0_(i*6-6+3)=Y0(i*6-6+2);
y0_(i*6-6+4)=Y0(i*6-6+2);
y0_(i*6-6+5)=Y0(i*6-6+1);
y0_(i*6-6+6)=Y0(i*6-6+1);
end
y0_(107)=0.02*100/6; 
elseif strcmp(celltype,'ST')
    for i=1:17
y0_(i*6-6+1)=Y0(i*6-6+4);
y0_(i*6-6+2)=Y0(i*6-6+4);
y0_(i*6-6+3)=Y0(i*6-6+5);
y0_(i*6-6+4)=Y0(i*6-6+5);
y0_(i*6-6+5)=Y0(i*6-6+6);
y0_(i*6-6+6)=Y0(i*6-6+6);
    end
y0_(107)=0.026*100/6; 
end
y0_(103:106)=0;
y0=y0_';


% %6 compartment 20%+15%+15%
% Y0=yout(:,end);
% L=Y0(107)+2*Y0(108);%length of half-cell at the end time point
% % a=Y0(107)-0.2*L; b=Y0(108)-0.5*L;
% if strcmp(celltype,'SW')
% for i=1:17
% y0_(i*6-6+1)=Y0(i*6-6+3);
% y0_(i*6-6+2)=1/3*Y0(i*6-6+2)+2/3*Y0(i*6-6+3);
% y0_(i*6-6+3)=Y0(i*6-6+2);
% y0_(i*6-6+4)=1/3*Y0(i*6-6+1)+2/3*Y0(i*6-6+2);
% y0_(i*6-6+5)=Y0(i*6-6+1);
% y0_(i*6-6+6)=Y0(i*6-6+1);
% end
% y0_(107)=0.02*20; y0_(108)=0.02*15;
% elseif strcmp(celltype,'ST')
%     for i=1:17
% y0_(i*6-6+1)=Y0(i*6-6+4);
% y0_(i*6-6+2)=2/3*Y0(i*6-6+4)+1/3*Y0(i*6-6+5);
% y0_(i*6-6+3)=Y0(i*6-6+5);
% y0_(i*6-6+4)=2/3*Y0(i*6-6+5)+1/3*Y0(i*6-6+6);
% y0_(i*6-6+5)=Y0(i*6-6+6);
% y0_(i*6-6+6)=Y0(i*6-6+6);
%     end
% y0_(107)=0.026*20; y0_(108)=0.026*15;
% end
% y0_(103:106)=0;
% y0=y0_';