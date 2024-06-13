%% Three-dimensional Skater Model
% Eline van der Kruk
% versie 2.0
% 28-01-2014

%% File information
% This file contains the Runge Kutta solver for the three-dimensional model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% classical Runge Kutta 4 step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sprintf('Rk4 odetmt start')
function [y,Labdas,Forces, Q, Fi, F_wrijv, Fwrb, ydotu, thetab] = rk4(b0,VAR,h,freqLPM,m_skater,alpha,fric_coef,mu,skate, VAR2)

[k1 Labdas Forces Q Fi F_wrijv Fwrb thetab] = odetmt(b0,VAR(1,:),m_skater,alpha,fric_coef,mu,skate,freqLPM);
b0new = b0+(h/2)*k1'*(1/freqLPM);
disp(k1)
[k2] = odetmt(b0new,VAR(2,:),m_skater,alpha,fric_coef,mu,skate,freqLPM);
b0new = b0+(h/2)*k2'*(1/freqLPM);

[k3] = odetmt(b0new,VAR(2,:),m_skater,alpha,fric_coef,mu,skate,freqLPM);
b0new = b0+h*k3'*(1/freqLPM);

[k4] = odetmt(b0new,VAR(3,:),m_skater,alpha,fric_coef,mu,skate,freqLPM);
% fprintf('k1=  : %d\n', k1);
% fprintf('k2=  : %d\n', k2);
% fprintf('k3=  : %d\n', k3);
% fprintf('k4=  : %d\n', k4);
y = b0 + (h/6)*(k1'+2*k2'+2*k3'+k4')*(1/freqLPM);
% fprintf('y=  : %d\n', y);
K = ((k1'+2*k2'+2*k3'+k4')/6);
% fprintf('K=  : %d\n', K);
% Rename variables (sequence VAR:
% US,VS,THETAS,US,VS,THETAS)
XB1 = y(1);
YB1 = y(2);
ZB1 = y(3);
dXB1 = y(4);
dYB1 = y(5);
dZB1 = y(6);
ddXB1 = K(4);
ddYB1 = K(5);%K(7);
ddZB1 = K(6);%K(8);

mm =3;

if skate == 'LS'
US1 = VAR(mm,1);
VS1 = VAR(mm,2);
WS1 = VAR(mm,3);
THETA_S1 = VAR(mm,4);

dUS1 = VAR(mm,9);
dVS1 = VAR(mm,10);
dWS1 = VAR(mm,11);
dTHETA_S1 = VAR(mm,12);

ddUS1 = VAR(mm,17);
ddVS1 = VAR(mm,18);
ddWS1 = VAR(mm,19);
ddTHETA_S1 = VAR(mm,20);

else
US1 = VAR(mm,5);
VS1 = VAR(mm,6);
WS1 = VAR(mm,7);
THETA_S1 = VAR(mm,8);    
 
dUS1 = VAR(mm,13);
dVS1 = VAR(mm,14);
dWS1 = VAR(mm,15);
dTHETA_S1 = VAR(mm,16);

ddUS1 = VAR(mm,21);
ddVS1 = VAR(mm,22);
ddWS1 = VAR(mm,23);
ddTHETA_S1 = VAR(mm,24);
end

% Coordinate Projection method
% Constraint skate is on the ice:
% sprintf('coordinate projection method')

% if strcmp(skate, 'LS')
% EPS1 = ZB1-WS1;
% iterat = 1;
% tol = 1e-12;
% maxiterat = 100;
%  
% while   max(abs(EPS1))>tol && maxiterat>iterat 
%         DEPS = [ 0, 0, 1];
%         delta = -DEPS'*inv(DEPS*DEPS')*EPS1';
%     
%          y(1:3) = y(1:3)+delta';
%     
%         XB1 = y(1);
%         YB1 = y(2);
%         ZB1 = y(3);
%         
%         iterat = iterat +1;
% end
% end
%     
if strcmp(skate, 'LS')
      EPS = -sin(THETA_S1)*(-dUS1*cos(THETA_S1) + dYB1 - dVS1*sin(THETA_S1) - dTHETA_S1*VS1*cos(THETA_S1) + dTHETA_S1*US1*sin(THETA_S1)) - cos(THETA_S1)*(dXB1 - dVS1*cos(THETA_S1) + dUS1*sin(THETA_S1) + dTHETA_S1*US1*cos(THETA_S1) + dTHETA_S1*VS1*sin(THETA_S1));
        D = ([ -cos(THETA_S1), -sin(THETA_S1), 0]); %(anders dan Danique?!)(minnen aangepast)
elseif strcmp(skate,'RS')
    EPS = -sin(THETA_S1)*(-dUS1*cos(THETA_S1) + dYB1 - dVS1*sin(THETA_S1) - dTHETA_S1*VS1*cos(THETA_S1) + dTHETA_S1*US1*sin(THETA_S1)) +cos(THETA_S1)*(-dUS1*sin(THETA_S1) + dVS1*cos(THETA_S1) + dXB1 - dTHETA_S1*US1*cos(THETA_S1) - dTHETA_S1*VS1*sin(THETA_S1));
    D = ( [ cos(THETA_S1), -sin(THETA_S1), 0]);
end



% Moore-Penrose pseudo invers
deltaX = D'*inv(D*D.')*-EPS;
deltaX;

% New position and velocities of the center of mass
y = y+[0,0,0,deltaX'];



% [ydotu Labdas Forces Q Fi F_wrijv Fwrb thetab] = odetmt(y,VAR2,m_skater,alpha,fric_coef,mu,skate,freqLPM);
[ydotu] = odetmt(y,VAR2,m_skater,alpha,fric_coef,mu,skate,freqLPM);


global g
