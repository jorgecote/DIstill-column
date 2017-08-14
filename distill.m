%*******************************************************************%
% Binary Distillation Column model for water-etanol mix             %
%       variables:                                                  %
%  t --> time                                                       %
%  x --> mole fraction of A at each stage                           %
%  xdot --> state derivatives                                       %
%*******************************************************************%
function xdot = distill(t,x)

global u

% Inputs (1):
% rr => Reflux Ratio (L/D)
rr=u.rr;

% Stages (14):
% x(1) - Reflux Drum Liquid Mole Fraction of Component A
% x(2) - Tray 1 - Liquid Mole Fraction of Component A
% .
% .
% .
% x(12) - Tray 11 - Liquid Mole Fraction of Component A (Feed Location)
% .
% .
% .
% x(13) - Tray 12 - Liquid Mole Fraction of Component A
% x(14) - Reboiler Liquid Mole Fraction of Component A

% Parameters
% Feed Flowrate (Kgmol/hr)
Feed = u.Feed;
% Mole Fraction of Feed
x_Feed = u.x_Feed;
% Distillate Flowrate (Kgmol/hr)
D=u.D;
% Flowrate of the Liquid in the Rectification Section (Kgmol/hr)
L=rr*D;
u.L=L;
% Vapor Flowrate in the Column (Kgmol/hr)
V=L+D;
u.V=V;
% Flowrate of the Liquid in the Stripping Section (Kgmol/hr)
FL=Feed+L;
u.FL=FL;
%Bottom Flowrate
B=Feed-D;
u.B=B;
% Total Molar Holdup in the Condenser (initial conditions)
Mtn=[10 1.3219 1.3383 1.3552 1.3732 1.3936 1.4176 1.4480 1.4895 1.5535 1.6730 2.3508 3.6787 10];


% Vapor Mole Fractions of Component A
% From the equilibrium assumption and mole balances
[y,T]=gammaphi(x,1.01);


% Compute xdot
xdot(1) = 1/Mtn(1)*V*(y(2)-x(1));
xdot(2) = 1/Mtn(2)*(L*(x(1)-x(2))-V*(y(2)-y(3)));
xdot(3) = 1/Mtn(3)*(L*(x(2)-x(3))-V*(y(3)-y(4)));
xdot(4) = 1/Mtn(4)*(L*(x(3)-x(4))-V*(y(4)-y(5)));
xdot(5) = 1/Mtn(5)*(L*(x(4)-x(5))-V*(y(5)-y(6)));
xdot(6) = 1/Mtn(6)*(L*(x(5)-x(6))-V*(y(6)-y(7)));
xdot(7) = 1/Mtn(7)*(L*(x(6)-x(7))-V*(y(7)-y(8)));
xdot(8) = 1/Mtn(8)*(L*(x(7)-x(8))-V*(y(8)-y(9)));
xdot(9) = 1/Mtn(9)*(L*(x(8)-x(9))-V*(y(9)-y(10)));
xdot(10) = 1/Mtn(10)*(L*(x(9)-x(10))-V*(y(10)-y(11)));
xdot(11) = 1/Mtn(11)*(L*(x(10)-x(11))-V*(y(11)-y(12)));
xdot(12) = 1/Mtn(12)*(Feed*x_Feed+L*x(11)-FL*x(12)-V*(y(12)-y(13)));
xdot(13) = 1/Mtn(13)*(FL*(x(12)-x(13))-V*(y(13)-y(14)));
xdot(14) = 1/Mtn(14)*(FL*x(13)-B*x(14)-V*y(14));

xdot = xdot';  % xdot must be a column vector
