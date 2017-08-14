%*******************************************************************%
% script to obtain Binary Distillation Column results               %
% Created by Jorge Cote                                             %
%*******************************************************************%

global u

% Steady State Initial Conditions for the States
x_ss(1:14)=[0.8060 0.7758 0.7473 0.7192 0.6905 0.66 0.6262 0.5871 0.5394 0.4771 0.3891 0.2574 0.1025 0.0159];
T_ini=[351.3126 351.4011 351.5037 351.6212 351.7574 351.9183 352.1147 352.3651 352.7020 353.1893 353.9689 355.3891 358.8806 368.2881];
y_ss=[0.8267 0.8060 0.7875 0.7699 0.7527 0.7352 0.7165 0.6959 0.6721 0.6431 0.6054 0.5524 0.4520 0.1733];
P=1.01;              %Presion total columna en bar
x_ss = x_ss';
%Inizialization variables
% Reflux Ratio
u.rr = 3;
% Distillation Flowrate
u.D=260;
% Feed flowrate
u.Feed=600;
% Feed molar composition
u.x_Feed=0.4;

% Final Time (hr)
tf = 1;

[t,x] = ode15s('distill',0:0.002:tf,x_ss);
for i=1:size(t,1)
    [y Temp]=gammaphi(x(i,:),P);
    y_(i,:)=y;
    Temp_(i,:)=Temp;
end
% Parse out the distillate composition
x1 = x(:,1);
x2 = x(:,14);
datax_d=[datax_d x1];
datax_b=[datax_b x2];
%Plot the results
figure(1);
 plot(t,x1);
figure(2);
plot(t,x2);
 figure(3);
 plot(T_ini);
 hold on
 plot(Temp_(1,:));