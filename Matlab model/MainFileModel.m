%% Simple 3D skater model
% I have first processed the data and saved it. Now just loading that data
% and commented all the processing steps to keep some consistency in the
% results. 

clc
clear all
close all

participant = 'P7'
file = 'Lap Fast 3';

alpha = 0.1;
fric_coef = 0.15;
mu = 0.003;

load('info.mat');
for i = 1:length(info.names)
    if strcmp(info.names{i},participant)>0
        m_skater = info.meetschaats(i);
    end
end

gv = 9.81;
        load('processed_data_P4.mat')
        load('measured_forces_P4.mat')
        load('Qualysis_data_P4.mat')
%%
% Eline van der Kruk
% 31-03-2015

% last modified 14-04-2015
% updated 23-04-2015
% clc
% clear all
% close all
% 
% participant = 'Paul-Yme'
% file = 'Lap Fast 3';
% 
% alpha = 0.1;
% fric_coef = 0.15;
% mu = 0.003;
% 
% %% Load data files
% 
% load('info.mat');
% for i = 1:length(info.names)
%     if strcmp(info.names{i},participant)>0
%         m_skater = info.meetschaats(i);
%     end
% end
% 
% gv = 9.81;
% %%
% [Visual3D, Meetschaats, QTM] = loadData(participant,file);
% 
% %% knip force data zodat alleen het stuk gelijk met V3D overblijft
% LTime = Meetschaats.LTime;
% RTime = Meetschaats.RTime;
% VTime = Visual3D.TIME.pX;
% 
% freqMS = 1/mean(diff(LTime));
% freqQual = 1/mean(diff(VTime));
% 
% [v startLV3D] = min(abs(LTime-Visual3D.newTime(1)));
% [v startRV3D] = min(abs(RTime-Visual3D.newTime(1)));
% [v endLV3D] = min(abs(LTime-Visual3D.newTime(end)));
% [v endRV3D] = min(abs(RTime-Visual3D.newTime(end)));
% 
% % Force data (instrumented skate)
% RSFn = Meetschaats.RFZ(startRV3D:endRV3D);
% LSFn = Meetschaats.LFZ(startLV3D:endLV3D);
% RSFl = -Meetschaats.RFX(startRV3D:endRV3D);
% LSFl = Meetschaats.LFX(startLV3D:endLV3D);
% 
% % dRSFn = 
% 
% ForceLS = [LSFl LSFn];
% ForceRS = [RSFl RSFn];
% 
% if length(ForceLS(:,2))>length(ForceRS(:,2))
%     Lorce = length(ForceRS(:,2));
% else Lorce = length(ForceLS(:,2));
% end

%% Define variables
% and resample naar 100Hz (Meetschaats)
% 
% for i = 1:length(RSFn)
%     ii = (startRV3D-1)+i;
%     [v kk] = min(abs(Visual3D.newTime-LTime(ii)));
%     
%     NVTime(i,:)= Visual3D.newTime(kk);
% 
%     % position data (Qualisys)
%     XB(i,:)  = -Visual3D.RTH.CGPos.pY(kk);
%     YB(i,:)  = Visual3D.RTH.CGPos.pX(kk);
%     ZB(i,:)  = Visual3D.RTH.CGPos.pZ(kk);
% 
%     XLS(i,:) = -Visual3D.LFT.CGPos.pY(kk);%Visual3D.L_FAL.pY(kk);%
%     YLS(i,:) = Visual3D.LFT.CGPos.pX(kk);%Visual3D.L_FAL.pX(kk);%
%     ZLS(i,:) = Visual3D.LFT.CGPos.pZ(kk);%Visual3D.L_FAL.pZ(kk);%
% 
%     XRS(i,:) = -Visual3D.RFT.CGPos.pY(kk);%Visual3D.R_FAL.pY(kk);%
%     YRS(i,:) = Visual3D.RFT.CGPos.pX(kk);%Visual3D.R_FAL.pX(kk);%
%     ZRS(i,:) = Visual3D.RFT.CGPos.pZ(kk);%Visual3D.R_FAL.pZ(kk);%
%     
%     % velocity data
%     VDXB(i,:)  = -Visual3D.RTH.CGVel.pY(kk);
%     VDYB(i,:)  = Visual3D.RTH.CGVel.pX(kk);
%     VDZB(i,:)  = Visual3D.RTH.CGVel.pZ(kk);
% 
%     VDXLS(i,:) = -Visual3D.LFT.CGVel.pY(kk);%Visual3D.L_FAL.pY(kk);%
%     VDYLS(i,:) = Visual3D.LFT.CGVel.pX(kk);%Visual3D.L_FAL.pX(kk);%
%     VDZLS(i,:) = Visual3D.LFT.CGVel.pZ(kk);%Visual3D.L_FAL.pZ(kk);%
% 
%     VDXRS(i,:) = -Visual3D.RFT.CGVel.pY(kk);%Visual3D.R_FAL.pY(kk);%
%     VDYRS(i,:) = Visual3D.RFT.CGVel.pX(kk);%Visual3D.R_FAL.pX(kk);%
%     VDZRS(i,:) = Visual3D.RFT.CGVel.pZ(kk);%Visual3D.R_FAL.pZ(kk);%
%     
%     % acceleration data
%     VDDXB(i,:)  = -Visual3D.RTH.CGAcc.pY(kk);
%     VDDYB(i,:)  = Visual3D.RTH.CGAcc.pX(kk);
%     VDDZB(i,:)  = Visual3D.RTH.CGAcc.pZ(kk);
% 
%     VDDXLS(i,:) = -Visual3D.LFT.CGAcc.pY(kk);%Visual3D.L_FAL.pY(kk);%
%     VDDYLS(i,:) = Visual3D.LFT.CGAcc.pX(kk);%Visual3D.L_FAL.pX(kk);%
%     VDDZLS(i,:) = Visual3D.LFT.CGAcc.pZ(kk);%Visual3D.L_FAL.pZ(kk);%
% 
%     VDDXRS(i,:) = -Visual3D.RFT.CGAcc.pY(kk);%Visual3D.R_FAL.pY(kk);%
%     VDDYRS(i,:) = Visual3D.RFT.CGAcc.pX(kk);%Visual3D.R_FAL.pX(kk);%
%     VDDZRS(i,:) = Visual3D.RFT.CGAcc.pZ(kk);%Visual3D.R_FAL.pZ(kk);%
% 
%     % orientation data (Qualisys)
%     LSroll(i,:) = deg2rad(Visual3D.ANG_GLOB_LF_ZXY.pY(kk));
%     RSroll(i,:) = deg2rad(-Visual3D.ANG_GLOB_RF_ZXY.pY(kk));
% 
%     LSyaw(i,:) = deg2rad(Visual3D.ANG_GLOB_LF_ZXY.pZ(kk));
%     RSyaw(i,:) = deg2rad(Visual3D.ANG_GLOB_RF_ZXY.pZ(kk));
%     
% end
%%
% Differentiate the data
% dt = 1/freqMS;
% 
%     % velocity
%     DXB  = [VDXB(1); diff(XB)/dt];
%     DYB  = [VDYB(1); diff(YB)/dt];
%     DZB  = [VDZB(1); diff(ZB)/dt];
% 
%     DXLS = [VDXLS(1); diff(XLS)/dt];
%     DYLS = [VDYLS(1); diff(YLS)/dt];
%     DZLS = [VDZLS(1); diff(ZLS)/dt];
% 
%     DXRS = [VDXRS(1); diff(XRS)/dt];
%     DYRS = [VDYRS(1); diff(YRS)/dt];
%     DZRS = [VDZRS(1); diff(ZRS)/dt];
% 
%     % acceleration
%     DDXB  = [VDDXB(1); diff(DXB)/dt];
%     DDYB  = [VDDYB(1); diff(DYB)/dt];
%     DDZB =  [VDDZB(1); diff(DZB)/dt];
% 
%     DDXLS = [VDDXLS(1); diff(DXLS)/dt];
%     DDYLS = [VDDYLS(1); diff(DYLS)/dt];
%     DDZLS = [VDDZLS(1); diff(DZLS)/dt];
% 
%     DDXRS = [VDDXRS(1); diff(DXRS)/dt];
%     DDYRS = [VDDYRS(1); diff(DYRS)/dt];
%     DDZRS = [VDDZRS(1); diff(DZRS)/dt];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Prepare data for the model%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Save processed data %%%%

% %% -------------------seperate the strokes-----------------------------%%
% 
% for i = 1:Lorce
%     if ForceLS(i,2)>ForceRS(i,2);
%         skate_Array(i) = 1;
%     else skate_Array(i) = 2;
%     end
% end
% 
% sc = skate_Array(1);
% ii = 1;
% Stroke(1,1) = 1;
% for i = 1:Lorce-1
%    if skate_Array(i) == sc
%        
%    else Stroke(ii,2) = i-1;
%         ii = ii+1;
%         Stroke(ii,1) = i;
%         sc = skate_Array(i);
%    end
%     
% end
% Stroke(end,2) = Lorce;
% 
%  figure
% a(1) = subplot(211)
%  plot(RSFn,'r');hold on
%  plot(LSFn,'b'); hold on
%  grid minor
% a(2) =  subplot(212)
%  plot(skate_Array,'k','Linewidth',2);
% grid minor
% linkaxes(a,'x')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fit variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % variables which have to be fitted
% aaa = 30;
% for ik = 1:length(Stroke(:,1));
%     try
%             T2 = Stroke(ik,1)-aaa:1:Stroke(ik,2)+aaa;
%             input_variables=[XB(T2),YB(T2),ZB(T2),XLS(T2),YLS(T2),ZLS(T2),XRS(T2),YRS(T2),ZRS(T2)];
%             s1 = aaa+1;
%             s2 = length(T2)-aaa;
%     catch
%         try
%                 T2 = Stroke(ik,1)-aaa:1:Stroke(ik,2);
%                 input_variables=[XB(T2),YB(T2),ZB(T2),XLS(T2),YLS(T2),ZLS(T2),XRS(T2),YRS(T2),ZRS(T2)];
%                 s1 = aaa+1;
%                 s2 = length(T2);
%         catch
%                 T2 = Stroke(ik,1):1:Stroke(ik,2)+aaa;
%                 input_variables=[XB(T2),YB(T2),ZB(T2),XLS(T2),YLS(T2),ZLS(T2),XRS(T2),YRS(T2),ZRS(T2)];
%                 s1 = 1;
%                 s2 = length(T2)-aaa;
%         end
%     end
% 
%     T = length(T2);
%  
%     [POS3,VEL3,ACC3,coeffval3]= fitfunction_new(input_variables,T,T2,freqMS);
%     POS(:,Stroke(ik,1):Stroke(ik,2)) = POS3(:,s1:s2);
%     VEL(:,Stroke(ik,1):Stroke(ik,2)) = VEL3(:,s1:s2);
%     ACC(:,Stroke(ik,1):Stroke(ik,2)) = ACC3(:,s1:s2);
%     coeffval.(char(strcat('n',(num2str(ik))))) = coeffval3;
% end
% 
%     T2 = 1:1:Lorce;
%     input_variables=[XB(T2),YB(T2),ZB(T2),XLS(T2),YLS(T2),ZLS(T2),XRS(T2),YRS(T2),ZRS(T2)];
%     rPOS = input_variables';
%     rVEL = [DXB(T2),DYB(T2),DZB(T2),DXLS(T2),DYLS(T2),DZLS(T2),DXRS(T2),DYRS(T2),DZRS(T2)]';
%     rACC = [DDXB(T2),DDYB(T2),DDZB(T2),DDXLS(T2),DDYLS(T2),DDZLS(T2),DDXRS(T2),DDYRS(T2),DDZRS(T2)]';
% 
% % mean fit error
% for ii = 1:length(POS(:,1))
%     error(ii) = sqrt(mean((input_variables(:,ii) - POS(ii,:)').^2));
% end
% 
% ii=2;
% figure
% subplot(211)
% plot(RSFn(T2),'r');hold on
% plot(LSFn(T2),'b');hold on
% grid minor
% subplot(212)
% % plot(rACC(ii,:),'k','Linewidth',2);hold on
% plot(ACC(ii,:),'k','Linewidth',2)
% % plot(input_variables(:,ii+3),'m');hold on
% % plot(POS(ii+3,:),'b')
% grid minor
% 
% % figure
% % plot(input_variables(:,ii)-POS(ii,:)');
% 
T1 = 1:1:length(POS(4,:));
% 
% SA = skate_Array(T2);
% 
% %% check data
% 
% iii = 1;
% figure
% for iii = 1:length(Stroke(:,1))
% plot(abs(coeffval.(char(strcat('n',(num2str(iii)))))(:,2:end)'),'.','Linewidth',2);hold on
% legend('XB','YB','ZB','XLS','YLS','ZLS','XRS','YRS','ZRS');
% end
% title('amplitude plot');xlabel('coefficient');ylabel('absolute amplitude')

%% Express positions in the generalized coordinates.

%%%%LET OP HIER WORDT POS RPOS etc.%%%%%%
% POS = rPOS;
% VEL = rVEL;
% ACC = rACC;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T2 = 1:1:420; 

        THETA_B     = zeros(1,length(T2));
        THETA_LS    = zeros(1,length(T2));
        VLS         = zeros(1,length(T2));
        ULS         = zeros(1,length(T2));
        WLS         = zeros(1,length(T2));
        THETA_RS    = zeros(1,length(T2));
        VRS         = zeros(1,length(T2));
        URS         = zeros(1,length(T2));
        WRS         = zeros(1,length(T2));
        
 for k = 1:length(T2)
            xb = POS(1,k);
            yb = POS(2,k);
            zb = POS(3,k);
            xls = POS(4,k);
            yls = POS(5,k);
            zls = POS(6,k);
            xrs = POS(7,k);
            yrs = POS(8,k);
            zrs = POS(9,k);
            
            dxb = VEL(1,k);
            dyb = VEL(2,k);
            dzb = VEL(3,k);
            dxls = VEL(4,k);
            dyls = VEL(5,k);
            dzls = VEL(6,k);
            dxrs = VEL(7,k);
            dyrs = VEL(8,k);
            dzrs = VEL(9,k);
            
            % NOTE: hoeken nog niet veranderd naar daadwerkelijke hoeken ->
            % uitzoeken of dit nodig is!
            
            % angle theta of body
            %THETA_B(k) = atan(dyb/dxb);
            THETA_B(k) = -asin(dxb/sqrt(((dyb)^2)+((dxb)^2))); %sin veranderd naar cos
            
            % generalized coordinates of the left skate
            THETA_LS(k) = -asin(dxls/sqrt(((dyls)^2)+((dxls)^2))); %- veranderen sin veranderd naar cos
            VLS(k) = (sin(THETA_LS(k))*(yb-yls)+cos(THETA_LS(k))*(xb-xls));
            ULS(k) = ((cos(THETA_LS(k))*(yb-yls)-sin(THETA_LS(k))*(xb-xls)));
            WLS(k) = zb-zls;
            
            % generalized coordinates of the right skate
            THETA_RS(k) = asin(dxrs/sqrt(((dyrs)^2)+((dxrs)^2)));%- veranderd
            VRS(k) = (sin(THETA_RS(k))*(yb-yrs)-cos(THETA_RS(k))*(xb-xrs));
            URS(k) = (cos(THETA_RS(k))*(yb-yrs)+sin(THETA_RS(k))*(xb-xrs));
            WRS(k) = zb-zrs;
   end


   %%
    
    
%     [POS1,VEL1,ACC1]=fitfunction1(input_variables1,T,T2,freqMS);
%     
%     aaa = 30;
% for ik = 1:length(Stroke(:,1));
%     try
%             T2 = Stroke(ik,1)-aaa:1:Stroke(ik,2)+aaa;
%             input_variables1=[THETA_B(T2)',ULS(T2)',VLS(T2)',WLS(T2)',THETA_LS(T2)',URS(T2)',VRS(T2)',WRS(T2)',THETA_RS(T2)'];
%             s1 = aaa+1;
%             s2 = length(T2)-aaa;
%     catch
%         try
%                 T2 = Stroke(ik,1)-aaa:1:Stroke(ik,2);
%                 input_variables1=[THETA_B(T2)',ULS(T2)',VLS(T2)',WLS(T2)',THETA_LS(T2)',URS(T2)',VRS(T2)',WRS(T2)',THETA_RS(T2)'];
%                 s1 = aaa+1;
%                 s2 = length(T2);
%         catch
%                 T2 = Stroke(ik,1):1:Stroke(ik,2)+aaa;
%                 input_variables1=[THETA_B(T2)',ULS(T2)',VLS(T2)',WLS(T2)',THETA_LS(T2)',URS(T2)',VRS(T2)',WRS(T2)',THETA_RS(T2)'];
%                 s1 = 1;
%                 s2 = length(T2)-aaa;
%         end
%     end
% 
%     T = length(T2);
%  
%     [POS3,VEL3,ACC3,coeffval3]= fitfunction_new(input_variables1,T,T2,freqMS);
%     POS1(:,Stroke(ik,1):Stroke(ik,2)) = POS3(:,s1:s2);
%     VEL1(:,Stroke(ik,1):Stroke(ik,2)) = VEL3(:,s1:s2);
%     ACC1(:,Stroke(ik,1):Stroke(ik,2)) = ACC3(:,s1:s2);
%     coeffval1.(char(strcat('n',(num2str(ik))))) = coeffval3;
% end
% 
% 
% figure
% for iii = 1:length(Stroke(:,1))
% plot(abs(coeffval1.(char(strcat('n',(num2str(iii)))))(:,2:end)'),'.','Linewidth',2);hold on
% legend('XB','YB','ZB','XLS','YLS','ZLS','XRS','YRS','ZRS');
% end
% title('amplitude plot');xlabel('coefficient');ylabel('absolute amplitude')
% 
% 
%  input_variables1 = [THETA_B',ULS',VLS',WLS',THETA_LS',URS',VRS',WRS',THETA_RS'];
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     input_variables1 = inpaint_nans(input_variables1);%deze functie nog checken!!!
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   %%  
%     % velocity
%     DTHETA_B     = [VEL1(1), diff(THETA_B)/dt]; %sin veranderd naar cos
%     
%     DTHETA_LS    = [VEL1(2), diff(THETA_LS)/dt]; %- veranderen sin veranderd naar cos
%     DVLS         = [VEL1(3), diff(VLS)/dt];
%     DULS         = [VEL1(4), diff(ULS)/dt];
%     DWLS         = [VEL1(5), diff(WLS)/dt];%-zls;
%     
%     DTHETA_RS    = [VEL1(6), diff(THETA_RS)/dt];%- veranderd
%     DVRS         = [VEL1(7), diff(VRS)/dt];
%     DURS         = [VEL1(8), diff(URS)/dt];
%     DWRS         = [VEL1(9), diff(WRS)/dt];%-zrs;
%     
%     % acceleration
%     DDTHETA_B     = [ACC1(1), diff(DTHETA_B)/dt]; %sin veranderd naar cos
%     
%     DDTHETA_LS    = [ACC1(2), diff(DTHETA_LS)/dt]; %- veranderen sin veranderd naar cos
%     DDVLS         = [ACC1(3), diff(DVLS)/dt];
%     DDULS         = [ACC1(4), diff(DULS)/dt];
%     DDWLS         = [ACC1(5), diff(DWLS)/dt];%-zls;
%     
%     DDTHETA_RS    = [ACC1(6), diff(DTHETA_RS)/dt];%- veranderd
%     DDVRS         = [ACC1(7), diff(DVRS)/dt];
%     DDURS         = [ACC1(8), diff(DURS)/dt];
%     DDWRS         = [ACC1(9), diff(DWRS)/dt];%-zrs;    
%%    
%     rPOS1 = input_variables1';
%     rVEL1 = [DTHETA_B',DULS',DVLS',DWLS',DTHETA_LS',DURS',DVRS',DWRS',DTHETA_RS'];
%     rACC1 = [DDTHETA_B',DDULS',DDVLS',DDWLS',DDTHETA_LS',DDURS',DDVRS',DDWRS',DDTHETA_RS'];
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LET OP!%%%%%%%%%%%%%%%%%%%%%%
%     %%%% POS1 wordt hier rPOS1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     POS1 = rPOS1;
% %     VEL1 = rVEL1';
% %     ACC1 = rACC1';
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     T2 = 1:1:Lorce;
%     
%     % mean fit error
% for ii = 1:length(POS(:,1))
%     error1(ii) = sqrt(mean((input_variables(:,ii) - POS(ii,:)').^2));
% end
%     
% ii=1;
% figure
% subplot(211)
% plot(RSFn(T2),'r');hold on
% plot(LSFn(T2),'b');hold on
% grid minor
% subplot(212)
% plot(input_variables1(:,ii),'k');hold on
% plot(POS1(ii,:),'r')
% % plot(input_variables(:,ii+3),'m');hold on
% % plot(POS(ii+3,:),'b')
% grid minor

    
%     %% Comparison Force measured by the instrumented skate
%  for i = 1:Lorce
% %     LRot_s2g = [-cos(LSroll(i)),-sin(LSroll(i));-sin(LSroll(i)),cos(LSroll(i));] ;% skate to global
%      LRot_s2g = [cos(LSroll(i)),-sin(LSroll(i));sin(LSroll(i)),cos(LSroll(i));] ;% skate to global
%      RRot_s2g = [-cos(RSroll(i)),sin(RSroll(i));-sin(RSroll(i)),cos(RSroll(i));]  ;% skate to global
% %     RRot_s2g = [-cos(RSroll(i)),0;-sin(RSroll(i)),cos(RSroll(i));]  ;% skate to global
%    
%     GForceLS(i,:) = LRot_s2g*ForceLS(i,:)';%[Fx Fz]
%     GForceRS(i,:) = RRot_s2g*ForceRS(i,:)';% [Fx Fz]
%  end


% %% Integrate the Upper Body Position
%       global kkk  KTHETAB GTHETA_B
%       GTHETA_B = THETA_B;
%       kkk = 1;
%     
% 
%         skate_Array = skate_Array(T2);
%         
%         q = [POS1(2,:)',POS1(3,:)',POS1(4,:)',POS1(5,:)',POS1(6,:)',POS1(7,:)',POS1(8,:)',POS1(9,:)',...
%             VEL1(2,:)',VEL1(3,:)',VEL1(4,:)',VEL1(5,:)',VEL1(6,:)',VEL1(7,:)',VEL1(8,:)',VEL1(9,:)',...
%             ACC1(2,:)',ACC1(3,:)',ACC1(4,:)',ACC1(5,:)',ACC1(6,:)',ACC1(7,:)',ACC1(8,:)',ACC1(9,:)'];
                                                      % [ULS VLS WLS THETALS URS VRS WRS THETARS]

        %%

        %%
        % Constants
        h = 1;                               % stepsize
        %T = length(T2);                      % timeperiod in 1/100s
        T = length(skate_Array);
        N = T/h;                             % number of steps
        
        iter = 2;
        
        % Interpolate data for runge kutta
       % q1 = interp1(q,1:(1/iter)*h:T,'spline');   % interpolate data with spline function
                %%
%         save('processed_data_P4.mat', 'POS', 'VEL', 'freqMS', 'q', 'skate_Array', 'q1')
%         save('measured_forces_P4.mat', 'GForceLS', 'GForceRS', 'LSFl', 'LSFn', 'RSFl', 'RSFn')

        %save('processed_data.mat', 'POS', 'VEL', 'skate_Array', 'freqMS', 'q', 'q1'); % List all variables you want to save
        
        b0= zeros(length(N-1),6);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Classical Runge-Kutta 4th Order
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Initial conditions
        b0(1,:) = [POS(1,1) POS(2,1) POS(3,1) VEL(1,1) VEL(2,1) VEL(3,1)];
   %%
        iter2 = 3;
        k = 1;
       for i = 1:N-1
            q2 = q1(k:k+iter2-1,:);
            i;
            if skate_Array(i)== 1
                skate = 'LS';
            elseif skate_Array(i)==2
                skate = 'RS';
            elseif isempty(skate)
                 continue
            end
                     
            [b0(i+1,:),SkateLabdas(i,:),SkateForces(i,:), Q(i,:),Fi(i,:),...
                F_wrijv(i,:), F_wrijvb(i,:), b0dot(i+1,:),thetab(i,:)] = ...
                rk4(b0(i,:),q2,h,freqMS,m_skater,alpha,fric_coef,mu,skate,q(i+1,:));
            k = k+iter;
        end
        
        
        %% Vermogen berekeningen
        F_wrijvy = -F_wrijvb.*cos(thetab)- F_wrijv.*cos(thetab);
        F_wrijvx = F_wrijvb.*sin(thetab)+ F_wrijv.*sin(thetab);;
               
%         Pwy = F_wrijvy.*DYB(2:end);
%         Pwx = F_wrijvx.*DXB(2:end);
        
        % human power globaal assenstelsel
        for i = 1:length(THETA_LS)
            if skate_Array(i) == 1
                Fhy(i) = -GForceLS(i,1).*sin(THETA_LS(i));
                Fhx(i) = -GForceLS(i,1).*cos(THETA_LS(i));
                Fhz(i) = GForceLS(i,2);
            else
                Fhy(i) = GForceRS(i,1).*sin(THETA_RS(i));
                Fhx(i) = -GForceRS(i,1).*cos(THETA_RS(i));
                Fhz(i) = GForceRS(i,2);
            end
        end
%         
%         Phy = Fhy'.*DYB(T2);
%         Phx = Fhx'.*DXB(T2);
%         Phz = Fhz'.*DZB(T2);
% %         
%         % human power schaats stelsel
%         for i = 1:length(THETA_LS)
%             if skate_Array(i) == 1
%                 Fhv(i) = -GForceLS(i,1);
%             else
%                 Fhv(i) = GForceRS(i,1);
% %             end
%         end
%         for i = 1:length(THETA_LS)
%             if skate_Array(i) == 1
%                 Phv(i) = -GForceLS(i,1)*DVLS(i);
%             else
%                 Phv(i) = GForceRS(i,1)*DVRS(i);
%             end
%         end
        
%         figure
%         ab(1) = subplot(211)
%         plot(RSFn(T2),'r');hold on
%         plot(LSFn(T2),'b');hold on
%         plot(RSFl(T2),'k');hold on
%         plot(LSFl(T2),'--k');hold on
%         grid minor
%         ab(2) = subplot(212)
%         plot(Pwx,'--k','Linewidth',2);hold on
%         plot(Pwy,'k','Linewidth',2);hold on
%         plot(Pwx+Pwy,'b','Linewidth',2);hold on
%         legend('Pwx','Pwy','Ptotal');ylabel('Friction Power')
%%  
        figure
        ab(1) = subplot(311)
        plot(RSFn(T2),'r');hold on
        plot(LSFn(T2),'b');hold on
        plot(RSFl(T2),'k');hold on
        plot(LSFl(T2),'--k');hold on
        grid minor
        ab(2) = subplot(312)
        plot(Fhx,'r','Linewidth',2);hold on
        plot(Fhy,'b','Linewidth',2);hold on
        plot(Fhz,'g','Linewidth',2);hold on
        legend('Fhx','Fhy','Fhz');ylabel('Human Force Directions')
        grid minor
        ab(3) = subplot(313)
%         plot(Phx,'r','Linewidth',2);hold on
%         plot(Phy,'b','Linewidth',2);hold on
%         plot(Phz,'g','Linewidth',2);hold on
% %         legend('Phx','Phy','Phz');ylabel('Human Power Directions')
% %         grid minor
%         linkaxes(ab,'x')
%%         
%         figure
%         ab(1) = subplot(311)
%         plot(RSFn(T2),'r');hold on
%         plot(LSFn(T2),'b');hold on
%         plot(RSFl(T2),'k');hold on
%         plot(LSFl(T2),'--k');hold on
%         grid minor
%         ab(2) = subplot(312)
%         plot(Fhz,'r','Linewidth',2);hold on
%         plot(Fhv,'b','Linewidth',2);hold on
%         legend('Fhz','Fhv');ylabel('Human Force Directions')
%         grid minor
%         ab(3) = subplot(313)
%         plot(Phz,'r','Linewidth',2);hold on
%         plot(Phv,'b','Linewidth',2);hold on
%         plot(Phz+Phv','g','Linewidth',2);hold on
%         legend('Phz','Phv','Phz+Phv');ylabel('Human Power (skate frame)')
%         grid minor
%         linkaxes(ab,'x')
% %         
%         figure
%         ab(1) = subplot(411)
%         plot(RSFn(T2),'r');hold on
%         plot(LSFn(T2),'b');hold on
%         plot(RSFl(T2),'k');hold on
%         plot(LSFl(T2),'--k');hold on
%         grid minor
%         ab(2) = subplot(412)
%         plot(Fhx,'r','Linewidth',2);hold on
%         plot(Fhy,'b','Linewidth',2);hold on
%         plot(Fhz,'g','Linewidth',2);hold on
%         legend('Fhx','Fhy','Fhz');ylabel('Human Force Directions')
%         grid minor
%         ab(3) = subplot(413)
%         plot(DXB(T2),'r','Linewidth',2);hold on
%         plot(DYB(T2),'b','Linewidth',2);hold on
%         plot(DZB(T2),'g','Linewidth',2);hold on
%         legend('dx','dy','dz');ylabel('Upper body Velocity Directions')
%         grid minor
%         linkaxes(ab,'x')
%          ab(4) = subplot(414)
%         plot(POS(1,:),'k','Linewidth',2);hold on
% for i = 1:length(POS(4,:));
%     if skate_Array(i) == 1
%         plot(T1(i),POS(4,i),'.b','Linewidth',2);hold on
%     else
%         plot(T1(i),POS(7,i),'.r','Linewidth',2);hold on
%     end
% end
%         grid minor; ylabel('Lateral Position');xlabel('frame')
%         
%         linkaxes(ab,'x')
% %         
%         figure
%         av1(1) = subplot(511)
%         plot(RSFn(T2),'r');hold on
%         plot(LSFn(T2),'b');hold on
%         plot(RSFl(T2),'k');hold on
%         plot(LSFl(T2),'--k');hold on
%         grid minor
%         av1(5) = subplot(512)
%         plot(rad2deg(RSroll(T2)),'r','Linewidth',2);hold on
%         plot(rad2deg(LSroll(T2)),'b','Linewidth',2);hold on
%         legend('roll rs','roll ls');ylabel('Upper body Velocity Directions')
%         grid minor
%         av1(2) = subplot(513)
%         plot(Fhx,'r','Linewidth',2);hold on
%         plot(Fhy,'b','Linewidth',2);hold on
%         plot(Fhz,'g','Linewidth',2);hold on
%         legend('Fhx','Fhy','Fhz');ylabel('Human Force Directions')
%         grid minor
%         av1(3) = subplot(514)
%         plot(THETA_RS,'r','Linewidth',2);hold on
%         plot(THETA_LS,'b','Linewidth',2);hold on
%         plot(THETA_B,'k','Linewidth',2);hold on
%         legend('theta rs','theta ls','theta b');ylabel('Upper body Velocity Directions')
%         grid minor
%         linkaxes(ab,'x')
%          av1(4) = subplot(515)
%         plot(POS(1,:),'k','Linewidth',2);hold on
% for i = 1:length(POS(4,:));
%     if skate_Array(i) == 1
%         plot(T1(i),POS(4,i),'.b','Linewidth',2);hold on
%     else
%         plot(T1(i),POS(7,i),'.r','Linewidth',2);hold on
%     end
% end
%         grid minor; ylabel('Lateral Position');xlabel('frame')
%         
%         linkaxes(av1,'x')
%         
%         figure
%         ab(1) = subplot(411)
%         plot(RSFn(T2),'r');hold on
%         plot(LSFn(T2),'b');hold on
%         plot(RSFl(T2),'k');hold on
%         plot(LSFl(T2),'--k');hold on
%         grid minor; ylabel('Normal and Lateral Force [N]')
%         ab(2) = subplot(412)
%         plot(Phy+Phx+Phz,'k','Linewidth',2);hold on
%         plot(Phz+Phv','.b','Linewidth',2);hold on
%         plot(abs(Pwx+Pwy),'--k','Linewidth',2);hold on
%         legend('Human power','Human power skateframe','Power Friction');ylabel('Human power')
%         grid minor
%         ab(3) = subplot(413)
%         for i = 1:length(skate_Array)
%             if skate_Array(i) == 2
%                 plot(T1(i), THETA_RS(i),'.r','Linewidth',2);hold on
%             else
%                 plot(T1(i), THETA_LS(i),'.b','Linewidth',2);hold on
%             end
%         end
%         legend('Steer Angle');ylabel('Steer Angle')
%         grid minor
%         
%         ab(4) = subplot(414)
%         plot(POS(1,:),'k','Linewidth',2);hold on
% for i = 1:length(POS(4,:));
%     if skate_Array(i) == 1
%         plot(T1(i),POS(4,i),'.b','Linewidth',2);hold on
%     else
%         plot(T1(i),POS(7,i),'.r','Linewidth',2);hold on
%     end
% end
%         grid minor; ylabel('Lateral Position');xlabel('frame')
                
        %% Test plots
        
T2 = length(skate_Array)      
%save('Qualysis_data_P4.mat', 'XB', 'YB', 'ZB', 'XLS', 'YLS', 'ZLS', 'XRS', 'YRS', 'ZRS', ...
%     'VDXB', 'VDYB', 'VDZB', 'VDXLS', 'VDYLS', 'VDZLS', 'VDXRS', 'VDYRS', 'VDZRS', ...
%     'VDDXB','VDDYB', 'VDDZB', 'VDDXLS', 'VDDYLS', 'VDDZLS', 'VDDXRS', 'VDDYRS', 'VDDZRS',  ...
%     'LSroll', 'RSroll', 'LSyaw', 'RSyaw');
% Figure for plotting the position
figure
subplot(311)
plot(XB,'k','Linewidth',2);hold on
plot(POS(1,:),'--r','Linewidth',2);hold on
plot(b0(:,1),'--b','Linewidth',2);hold on
legend('Qualisys','Fitted','Modelled')
ylabel('Lateral Position [m]');xlabel('Frames [-]')
subplot(312)
plot(YB,'k','Linewidth',2);hold on
plot(POS(2,:),'--r','Linewidth',2);hold on
plot(b0(:,2),'--b','Linewidth',2);hold on
legend('Qualisys','Fitted','Modelled')
ylabel('Longitudinal Position [m]');xlabel('Frames [-]')
subplot(313)
plot(ZB,'k','Linewidth',2);hold on
plot(POS(3,:),'--r','Linewidth',2);hold on
%plot(POS1(4,:),'--m','Linewidth',2);hold on
%plot(POS1(8,:),'--g','Linewidth',2);hold on
plot(b0(:,3),'--b','Linewidth',2);hold on
legend('Qualisys','Fitted','Modelled')
ylabel('Upward Position [m]');xlabel('Frames [-]')
%% Velocity and force 
figure
ab(1) = subplot(411)
plot(RSFn,'r');hold on
plot(LSFn,'b');hold on
plot(RSFl,'k');hold on
plot(LSFl,'--k');hold on
legend('RSFn', 'LSFn', 'RSFl', 'LSFl')
ylabel('Force')
grid minor
ab(2) =subplot(412)
plot(VEL(1,:),'--r','Linewidth',2);hold on
plot(b0(:,4),'--b','Linewidth',2);hold on
legend('Fitted','Modelled')
ylabel('Lateral Velocity[m/s]');xlabel('Frames [-]')
ab(3) =subplot(413)
plot(VEL(2,:),'--r','Linewidth',2);hold on
plot(b0(:,5),'--b','Linewidth',2);hold on
legend('Fitted','Modelled')
ylabel('Longitudinal Velocity [m/s]');xlabel('Frames [-]')
ab(4) =subplot(414)
plot(VEL(3,:),'-r','Linewidth',2);hold on
plot(b0(:,6),'--b','Linewidth',2);hold on
legend('Fitted','Modelled')
ylabel('Upward Velocity [m/s]');xlabel('Frames [-]')

linkaxes(ab,'x')
%%
figure
a(1)= subplot(311)
plot(RSFn,'r');hold on
plot(LSFn,'b');hold on
plot(RSFl,'k');hold on
plot(LSFl,'--k');hold on
a(2)= subplot(312)
plot(-GForceLS(:,1),'b');hold on
plot(GForceRS(:,1),'r');hold on
plot(SkateLabdas,'k');hold on
a(3)= subplot(313)
plot(LSroll(T2),'b');hold on
plot(RSroll(T2),'r');hold on
linkaxes(a,'x')
%%
hold off
% 
% 
% 
%%
% figure
% ac(1) = subplot(311)
% plot(RSFn,'r');hold on
% plot(LSFn,'b');hold on
% plot(RSFl,'k');hold on
% plot(LSFl,'--k');hold on
% legend('RSFn', 'LSFn', 'RSFl', 'LSFl')
% grid minor
% ac(2) = subplot(312)
% plot(abs(F_wrijvy),'m','Linewidth',2);hold on
% for i = 1:length(skate_Array)-1;
%     if skate_Array(i) == 1
%     plot(SkateLabdas(i).*sin(THETA_LS(i+1))','.b');hold on
%     plot((-GForceLS(i,1)).*sin(THETA_LS(i))','.k');hold on
%     else
%     plot(SkateLabdas(i).*sin(THETA_RS(i+1))','.r');hold on
%     plot(GForceRS(i,1).*sin(THETA_RS(i))','.k');hold on
%     end
% end
% grid minor; title ('voorwaarte kracht (kracht in y)')
% legend('Luchtwrijving','Modelled', 'Measured'); 
% ac(3) = subplot(313)
% plot(abs(F_wrijvx),'m','Linewidth',2);hold on
% for i = 1:length(skate_Array)-1;
%     if skate_Array(i) == 1
%     plot(SkateLabdas(i).*cos(THETA_LS(i+1))','.b');hold on
%     plot(-GForceLS(i,1).*cos(THETA_LS(i))','.k');hold on
%     else
%     plot(SkateLabdas(i).*cos(THETA_RS(i+1))','.r');hold on
%     plot(GForceRS(i,1).*cos(THETA_RS(i))','.k');hold on
%     end
% end
% grid minor; title ('zijwaartse kracht (kracht in x)')
% legend('Luchtwrijving','Modelled', 'Measured'); 
% linkaxes(ac,'x')
% hold off

% figure
% plot((THETA_LS(2:end))','k');hold on
% plot((THETA_RS(2:end))','--k');hold on



%% test fitfunctions
% 
% tVAR = WLS;
% dt = 1/freqMS;
% Ttime = 0:dt:T1(end)*dt;
% 
% tdVAR = diff(tVAR)/dt;
% 
% [B,A]=butter(2,10/50,'low');
% ftdVAR = filtfilt(B,A,tdVAR);
% 
% tddVAR = diff(tdVAR)/dt;
% 
% figure
% ac(1) = subplot(411)
% plot(RSFn(T2),'r');hold on
% plot(LSFn(T2),'b');hold on
% plot(RSFl(T2),'k');hold on
% plot(LSFl(T2),'--k');hold on
% grid minor
% ac(2) = subplot(412)
% plot(T2,tVAR,'k');
% ylabel('position [m]')
% ac(3) = subplot(413)
% plot(T2(2:end),tdVAR,'b');hold on
% plot(T2(2:end),ftdVAR,'m');
% % plot(T2(2:end),DXB,'r');
% ylabel('velocity [m/s]')
% ac(4) = subplot(414)
% plot(T2(3:end),tddVAR,'r');hold on
% % plot(T2,DDXB(T2),'m')
% ylabel('acceleration [m/s^2]'); 
% 
% linkaxes(ac,'x')
% 
% % figure
% % subplot(311)
% % plot(T2,tVAR(T2),'k');
% % ylabel('position [m]')
% % subplot(312)
% % plot(T2(2:end),tdVAR(T2(2:end)),'b');hold on
% % plot(T2(2:end),ftdVAR(T2(2:end)),'m');
% % plot(T2(2:end),DXB(T2(2:end)),'r');
% % ylabel('velocity [m/s]')
% % subplot(313)
% % plot(T2(2:end),tddVAR(T2(2:end)),'r');hold on
% % plot(T2,DDXB(T2),'m')
% % ylabel('acceleration [m/s^2]'); 
% 
% 
