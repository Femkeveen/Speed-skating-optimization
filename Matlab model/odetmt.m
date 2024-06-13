%% Three-dimensional Skater Model
% Eline van der Kruk
% versie 2.0
% 28-01-2014

%% File information
% This file contains contains the data from the file 'Equations of Motion'
% and finds the solution to the differential equation.

function [ydot Labda FSkate Q Fi Fs Fb thetab] = odetmt(b0,q,m_skater,alpha,fric_coef,mu,skate,freqLPM)
% sprintf('start ode tmt')

mb = (1-(alpha))*m_skater;
ms = (alpha)*m_skater;
Ms = 0;
Ib = 0;%1.26;%0.044;% %0.044;% [Ackland, Blanksby, Bloomfield]% 1.26 %http://www.phys.washington.edu/users/jeff/courses/ken_young_webs/208A/body_segment_mass.txt
Is = 0;%(1/12)*mrs*(269e-3)^2;% 0.0038;

locd = 1:3;
loco = 4:6;


%% extract parameters
xb = b0(1);
yb = b0(2);
zb = b0(3);
ub = b0(1);
vb = b0(2);
% ws = b0(3);

dxb = b0(4);
dyb = b0(5);
dzb = b0(6);
dub = b0(4);
dvb = b0(5);
% dws = b0(6);

if skate == 'LS'
us = q(1);
vs = q(2);
ws = q(3);
thetas = q(4);


dus = q(9);
dvs =q(10);
dws = q(11);
dthetas = q(12);

ddus = q(17) ;
ddvs = q(18);
ddws = q(19);
ddthetas = q(20);
end

if skate == 'RS'
us = q(5);
vs = q(6);
ws = q(7);
thetas = q(8);

dus = q(13);
dvs = q(14);
dws = q(15);
dthetas = q(16);

ddus = q(21);
ddvs = q(22);
ddws = q(23);
ddthetas = q(24);
end


%% Define the vectors

% variables
qd = [ub vb ws]';%b0(1:3)';
dqd = [dub dvb dws]';%b0(4:6)';
qo = [us vs thetas]';
dqo = [dus dvs dthetas]';

qi = [qd;qo];
dqi = [dqd;dqo];
% ddqi = [ddub ddvb ddwb ddus ddvs ddthetas];


%% find the solution
if skate == 'LS';
    kk = 1;
    aa = 0;
else kk = -1;
    aa = 4;
end

dxs = b0(4)+kk*q(12+aa)*sin(q(4+aa))*q(2+aa)-...
    kk*cos(q(4+aa))*q(10+aa)+kk*q(12+aa)*cos(q(4+aa))*q(1+aa)...
    +kk*q(9+aa)*sin(q(4+aa)); %alleen als er luchtwrijving is op de schaats?


dys = b0(5)+sin(q(4+aa))*(-q(10+aa)+q(12+aa)*q(1+aa))...
    -cos(q(4+aa))*(q(12+aa)*q(2+aa)+q(13+aa));%alleen als er luchtwrijving is op de schaats?


% syms Qus Qvs QMs real

Mij = (diag([mb mb mb ms ms Is])); %mass matrix
[Fb Fs thetab] = Friction(qi,dqi,freqLPM,m_skater,alpha,fric_coef,mu,skate,dxs,dys);
fi =([ Fb*sin(thetab), -Fb*cos(thetab), -(981*mb)/100, Fs*kk*sin(thetas), -Fs*cos(thetas), kk*Ms]);
[T,gcon, Cktot, Cki] = Trans(qi,dqi,kk);

T = (T);
gcon = (gcon);
Cktot = (Cktot);
Cki = (Cki);
% 
% if skate == 'RS' %test om deze weg te laten
%     Q = [0 0 0 (q(23)+9.81)*m_skater 0 0];
%     else Q = [0 0 0 (q(19)+9.81)*m_skater 0 0];
%     end
%  

% Reduced mass matrix and force matrix
Mredtot = (T'*Mij*T);  
fredtot = (T'*(fi'-Mij*gcon));%+Q'; 
Fi = (T'*(fi'-Mij*gcon)); %

% solutions
Ltot = ([Mredtot Cktot';
        Cktot zeros(1,1)]);
Rtot = (-[fredtot; -Cki]);
    
    Mdd = (Mredtot(locd,locd));
    Mdo =(Mredtot(locd,loco));
    Mod = (Mredtot(loco,locd));
    Moo = (Mredtot(loco,loco));
    
    Ckd = (Cktot(1,locd));
    Cko = (Cktot(1,loco));

    dqd = (dqi(locd)');
    dqo = (dqi(loco)');
    ddqo = ([ddus ddvs ddthetas]');%ddqi(4:6)';
        
    Fd = (fredtot(locd));
    Fo = (fredtot(loco));
 
% ---------------------------------------------------------------------%
%  sol1 = inv([Mdd Ckd'; Ckd zeros(1,1)])*[Fd-Mdo*ddqo ;-Cki-Cko*ddqo];
sol1 = ([Mdd Ckd'; Ckd zeros(1,1)])\[Fd-Mdo*ddqo ;-Cki-Cko*ddqo];

    ddqd = sol1(locd);
    if skate == 'RS' %test om deze weg te laten
    ddqd(3) = q(23);
    else ddqd(3) = q(19);
    end
        
 sol2 = [Mod Moo Cko']*[ddqd;ddqo;sol1(4)];
%------------------------------------------------------------------------%
% ddqd = sol1(1:3);
ddxb = sol1(1);
ddyb = sol1(2);
ddzb = sol1(3);

ddqo = [ddus ddvs ddthetas]';
Labda = (sol1(4));
                             
FSkate = (sol2);

Fi = (Fi);
Q = FSkate-Fi(loco);
ydot = ([dqd';ddqd]);


end






