%********************************************************************%
%       Computing of thermodynamic equilibrium                       %
%           variables:                                               %
%               * x: etanol concentration                            %
%               * P: Column pressure                                 %
%********************************************************************%
function [y,T]=gammaphi(x,P)

% Initial conditions
% Teperature initial condition
T_ini=[351.3126 351.4011 351.5037 351.6212 351.7574 351.9183 352.1147 352.3651 352.7020 353.1893 353.9689 355.3891 358.8806 368.2881];
% Vapor concentration initial condition
yyyy=[0.8267 0.8060 0.7875 0.7699 0.7527 0.7352 0.7165 0.6959 0.6721 0.6431 0.6054 0.5524 0.4520 0.1733];

%Parámetros ecuación de Antoine para etanol grados K y presion en bar 
A1=5.2567;    
B1=1598.673;
C1=-46.424;       
%Parámetros ecuación de Antoine para water grados K y presion en bar 
A2=5.0768;    
B2=1659.793;
C2=-45.854; 
%************************************************************************
%Coeficientes margules para hallar gamma
A12=1.6022;
A21=0.7947;
%coeficientes van-laar
AA12=1.6789;
AA21=0.9227;
%calculado por margules
gamma1=exp((A12+(2.*x.*(A21-A12))).*(1-x));
gamma2=exp((A21+(2.*(1-x).*(A12-A21))).*(x));
%calculado por van laar
gamma11=exp(AA12.*(AA21.*(1-x)./((AA12.*x)+(AA21.*(1-x)))).^2);
gamma22=exp(AA21.*(AA12.*(x)./((AA12.*x)+(AA21.*(1-x)))).^2);
T0=300;
% Compute bubble temperatures by margules coeficients
for i=1:14
Temp(i)=fzero(@(T) (P-(x(i)*gamma1(i)*10^(A1-(B1/(C1+T))))-((1-x(i))*gamma2(i)*10^(A2-(B2/(C2+T))))),T0);
psat1(i)=10^(A1-(B1/(C1+Temp(i))));
p1(i)=psat1(i)*gamma1(i)*x(i);
y1(i)=p1(i)/P;
%Compute bubble temperature by van laar coeficients
Temp1(i)=fzero(@(T) (P-(x(i)*gamma11(i)*10^(A1-(B1/(C1+T))))-((1-x(i))*gamma22(i)*10^(A2-(B2/(C2+T))))),T0);
psat11(i)=10^(A1-(B1/(C1+Temp1(i))));
p11(i)=psat11(i)*gamma11(i)*x(i);
y11(i)=p11(i)/P;
end
T=Temp1;
y=y11;
end