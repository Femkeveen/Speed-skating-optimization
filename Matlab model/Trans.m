%% Tranformation functions
% Eline van der Kruk
% 25-02-2014

function [T,gcon, Cktot, Cki] = Trans(qi,dqi,kk);

% variables
ub      = qi(1);
vb      = qi(2);
ws      = qi(3);
us      = qi(4);
vs      = qi(5);
thetas  = qi(6);

dub      = dqi(1);
dvb      = dqi(2);
dws      = dqi(3);
dus      = dqi(4);
dvs      = dqi(5);
dthetas  = dqi(6);

% transformation matrix
T = [   1, 0, 0,              0,               0,                                     0;
        0, 1, 0,              0,               0,                                     0;
        0, 0, 1,              0,               0,                                     0;
        1, 0, 0, kk*sin(thetas), -kk*cos(thetas), kk*us*cos(thetas) + kk*vs*sin(thetas);
        0, 1, 0,   -cos(thetas),    -sin(thetas),       us*sin(thetas) - vs*cos(thetas);
        0, 0, 0,              0,               0,                                    kk];
% Convective acceleration terms
gcon = [                                                                                                                                                     0;
                                                                                                                                                             0;
                                                                                                                                                             0;
 dthetas*(dthetas*(kk*vs*cos(thetas) - kk*us*sin(thetas)) + dus*kk*cos(thetas) + dvs*kk*sin(thetas)) + dthetas*dus*kk*cos(thetas) + dthetas*dvs*kk*sin(thetas);
                   dthetas*(dthetas*(us*cos(thetas) + vs*sin(thetas)) - dvs*cos(thetas) + dus*sin(thetas)) - dthetas*dvs*cos(thetas) + dthetas*dus*sin(thetas);
                                                                                                                                                             0];
                                                                                                                                                         
Cktot = [ -kk*cos(thetas), -sin(thetas), 0, - cos(thetas)*sin(thetas)*kk^2 + cos(thetas)*sin(thetas), kk^2*cos(thetas)^2 + sin(thetas)^2, sin(thetas)*(vs*cos(thetas) - us*sin(thetas)) - kk*cos(thetas)*(kk*us*cos(thetas) + kk*vs*sin(thetas))];
 
% Cki = dthetas*(dvs*(- 2*cos(thetas)*sin(thetas)*kk^2 + 2*cos(thetas)*sin(thetas)) - dvb*cos(thetas) + dthetas*(cos(thetas)*(vs*cos(thetas) - us*sin(thetas)) - sin(thetas)*(us*cos(thetas) + vs*sin(thetas)) - kk*cos(thetas)*(kk*vs*cos(thetas) - kk*us*sin(thetas)) + kk*sin(thetas)*(kk*us*cos(thetas) + kk*vs*sin(thetas))) + dus*(cos(thetas)^2 - sin(thetas)^2 - kk^2*cos(thetas)^2 + kk^2*sin(thetas)^2) + dub*kk*sin(thetas)) - dthetas*dus*(kk^2*cos(thetas)^2 + sin(thetas)^2) + dthetas*dvs*(- cos(thetas)*sin(thetas)*kk^2 + cos(thetas)*sin(thetas));

Cki = dthetas* (dvs*(- 2*cos(thetas)*sin(thetas)*kk^2 + 2*cos(thetas)*sin(thetas)) ...
                - dvb*cos(thetas) + ...
                 dthetas*(cos(thetas)*(vs*cos(thetas)- us*sin(thetas)) ...
                    - sin(thetas)*(us*cos(thetas) + vs*sin(thetas)) ...
                    - kk*cos(thetas)*(kk*vs*cos(thetas) - kk*us*sin(thetas)) ...
                    + kk*sin(thetas)*(kk*us*cos(thetas) + kk*vs*sin(thetas))) ...
                 + dus*(cos(thetas)^2 - sin(thetas)^2 - kk^2*cos(thetas)^2 + kk^2*sin(thetas)^2) ...
                 + dub*kk*sin(thetas)) ...
        - dthetas*dus*(kk^2*cos(thetas)^2 + sin(thetas)^2) ...
       + dthetas*dvs*(- cos(thetas)*sin(thetas)*kk^2 + cos(thetas)*sin(thetas));

end

