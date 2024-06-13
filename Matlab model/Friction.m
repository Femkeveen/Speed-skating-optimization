%% Skater Model 
% Feb 2013, E van der Kruk


function [Fb_friction F_skate THETA_B] = Friction(qi,dqi,freqLPM,m,alpha,...
    k1,mu,skate,dxs,dys)
g = 9.81;

global gain

mb = (1-(alpha))*m;
ms = (alpha)*m;

%Body
xb = qi(1);
yb = qi(2);
dxb = dqi(1);
dyb = dqi(2);

vb = norm([dxb;dyb]);
Fb_friction = k1*(vb)^2; %aangepast op 11-12-13
% Fb_friction = (vb)^2; 
% Fb_friction = (1-alpha)*k1*(vb)^2; %aangepast op 11-12-13
% Fb_friction = (1-gain)*k1*(vb)^2; %aangepast op 11-12-13


THETA_B = -asin(dxb/sqrt(((dyb)^2)+((dxb)^2)));


% vs = norm([dxs;dys]);
F_skate = (mu*m*g); %geen luchtwrijving op de schaats
% F_skate = (mu*m*g)+alpha*k1*(vs)^2 ;
% F_skate = (mu*m*g)+gain*k1*(vs)^2 ;

end


