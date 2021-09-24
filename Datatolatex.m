clear all; close all; clc
syms alpha1 h1op beta1 He uop 
syms s real 
syms h1 h2 u
% syms He beta1 h1op adis b q_1 q_2 q_3 q_0 h1opi
% alpha1 = 0.0004949; 

beta1 = 0.15;
He = 0.52;
V=[];
h1op=[];
A=[];
B=[];
G=sym(zeros(100, 1))';
%Gl={};%sym(zeros(100, 1))';
Gl={};
% Parameters
% adis=50;
% q_3=0.0000000702626295;
% q_2=-0.000000123003337;
% q_1=0.0000000734410693;
% q_0=0.00000000600146569;
% b=1.82800065;
% for i=1:10   %i=1:100
%     h1op(i)=i/20; %h1op(i)=i/200;
%     uop=adis^b*(q_3*h1op(i)^3+q_2*h1op(i)^2+q_1*h1op(i)+q_0);
%     A(i)=-adis^b*(3*h1op(i)^2*q_3+2*h1op(i)*q_2+q_1)/(beta1*sqrt(1-(He-h1op(i))^2/He^2))-(uop-adis^b*(h1op(i)^3*q_3+h1op(i)^2*q_2+h1op(i)*q_1+q_0))*(He-h1op(i))/(beta1*(1-(He-h1op(i))^2/He^2)^(3/2)*He^2);
%     B(i)=1/(beta1*sqrt(1-((He-h1op(i))/(He))^2));
%     G(i)=vpa(B(i)/(s-A(i)),6); %vpa(B(i)/(s-A(i)),4)
% end


h1opx=0.05:0.05:1;
udisx=100*ones(1,10);
uopx=[];
Fu1=(u-(0.0000000702626295*h1^3 - 0.000000123003337*h1^2 +...
+0.0000000734410693*h1+ 0.00000000600146569)*udis^1.82800065)/(beta1*sqrt(1-((He-h1)/(He))^2));
for k=1:10 
    uopx(k)=(0.0000000702626295*h1opx(k).^3 - 0.000000123003337*h1opx(k).^2 + ...
    0.0000000734410693*h1opx(k) + 0.0000000060014657)*udisx(k).^1.82800065;
    Ax(k)=double(subs(diff(Fu1,h1),{h1 u},{h1opx(k) uopx(k)}));
    Bx(k)=double(subs(diff(Fu1,u),{h1 u},{h1opx(k) uopx(k)}));
end

%% Converting G(s) to equation type eneviroment in latex
for k=1:10
   Gl{k}=latex((vpa(G(k),6))); %latex((vpa(G(k),4)));
   Gl1(k)=strrep(Gl(k),'\','\(\displaystyle\' );
   Gl2(k)=strcat(Gl1(k),'\)');
end
%% Converting numbers to equation type enviroment in latex
for k=1:10
   h11op{k}=latex((vpa(h1op(k),3)));
   hl1op(k)=strcat('$',h11op(k),'$');
   Al{k}=latex((vpa(A(k),5)));  %Al{k}=latex((vpa(A(k),5)));
   Al1(k)=strcat('$',Al(k),'$');
   Bl{k}=latex((vpa(B(k),6)));   %latex((vpa(B(k),4)));
   Bl1(k)=strcat('$',Bl(k),'$');
end
%
T=table(hl1op',Al1',Bl1',Gl2','VariableNames',{'h1op','A','B','Gs'});
Ta=T(1:10,:);
%Tb=T(6:10,:);
Tb.Properties.VariableNames = {'h2op','A1','B1','Gs1'};
%tablejoin=[Ta,Tb];
%% table2latex(T(1:5,:))
%table2latex(tablejoin)
table2latex(Ta)
