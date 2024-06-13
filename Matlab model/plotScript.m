% Plot script

%% Leg Extension opgedeeld per slag
T2 = 1:1:length(skate_Array);
T1 = T2
figure
a(1) = subplot(311)
plot(RSFn(T2),'r','Linewidth',2);hold on
plot(LSFn(T2),'b','Linewidth',2);hold on
grid minor;title('Leg extension')
a(2) = subplot(312)
for i = 1:length(ULS)
    if skate_Array(i)==1
        plot(T1(i),ULS(i),'.b','Linewidth',2);hold on
    else
        plot(T1(i),URS(i),'.r','Linewidth',2);hold on
    end
end
    legend('U');ylabel('extension [m]')
    grid minor
a(2) = subplot(313)
for i = 1:length(ULS)
    if skate_Array(i)==1
        plot(T1(i),VLS(i),'.b','Linewidth',2);hold on
    else
        plot(T1(i),VRS(i),'.r','Linewidth',2);hold on
    end
end
    legend('V');ylabel('extension [m]')
    grid minor;linkaxes(a,'x')

    %% Force en skate array
figure
a(1) = subplot(211)
plot(RSFn(T2),'r','Linewidth',2);hold on
plot(LSFn(T2),'b','Linewidth',2);hold on
grid minor;title('Leg extension')
a(2) = subplot(212)
plot(skate_Array)

%% Leg extension geheel
figure
a(1) = subplot(311)
plot(RSFn(T2),'r','Linewidth',2);hold on
plot(LSFn(T2),'b','Linewidth',2);hold on
grid minor;title('Leg extension')
a(2) = subplot(312)
plot(ULS,'b','Linewidth',2);hold on
plot(URS,'r','Linewidth',2);hold on
legend('ULS','URS');ylabel('extension [m]')
grid minor
a(3) = subplot(313)
plot(VLS,'b','Linewidth',2);hold on
plot(VRS,'r','Linewidth',2);hold on
legend('VLS','VRS');ylabel('extension [m]')
grid minor
linkaxes(a,'x')

%% Check magnitude of forces
% MFL = [(ForceLS(:,1).^2+ForceLS(:,2).^2)];
% MFR = [(ForceRS(:,1).^2+ForceRS(:,2).^2)];
GMFL = [(GForceLS(:,1).^2+GForceLS(:,2).^2)];
GMFR = [(GForceRS(:,1).^2+GForceRS(:,2).^2)];

figure
% plot(MFL,'b','Linewidth',2);hold on
% plot(MFR,'r','Linewidth',2);hold on
plot(GMFL,'.k','Linewidth',2);hold on
plot(GMFR,'.k','Linewidth',2);hold on

%% check the upper body angle

figure
a5(1)= subplot(211)
plot(RSFn(T2),'r');hold on
plot(LSFn(T2),'b');hold on
plot(RSFl(T2),'k');hold on
plot(LSFl(T2),'--k');hold on
grid minor
a5(2) = subplot(212)
plot(rad2deg(THETA_B),'b','Linewidth',2);hold on
plot(rad2deg(thetab),'r','Linewidth',2)
grid minor

linkaxes(a5,'x')

%% check steer angle
 figure
a(1) = subplot(311)
plot(RSFn(T2),'r');hold on
plot(LSFn(T2),'b');hold on
grid minor;title('Steer Angle Comparison')
a(2) =subplot(312)
plot(90+LSyaw,'b');hold on
plot(rad2deg(THETA_LS),'k');hold on
grid minor;ylim([-30 30]);legend('steer angle Visual3D','Angle velocity')
a(3) = subplot(313)
plot(90+RSyaw,'r');hold on
plot(-rad2deg(THETA_RS),'k');hold on
grid minor;ylim([-30 30]);legend('steer angle Visual3D','Angle velocity')
linkaxes(a,'x')

%% detrend plot

figure
ac(1) = subplot(211)
plot(RSFn(T2),'r');hold on
plot(LSFn(T2),'b');hold on
plot(RSFl(T2),'k');hold on
plot(LSFl(T2),'--k');hold on
grid minor
ac(1) = subplot(212)
plot(ZLS,'b');hold on
plot(ZRS,'r');hold on

linkaxes(ac,'x')

LTrend = detrend(ZLS);
RTrend = detrend(ZRS);

figure
ac(1) = subplot(211)
plot(RSFn(T2),'r');hold on
plot(LSFn(T2),'b');hold on
plot(RSFl(T2),'k');hold on
plot(LSFl(T2),'--k');hold on
grid minor
ac(2) = subplot(212)
plot(XLS(T2),'b');hold on
% plot(LTrend(T2)+0.2,'k');hold on
plot(XRS(T2),'r');hold on
% plot(RTrend(T2)+0.2,'--k');hold on

 %% Forces etc.
 
 figure
plot(F_wrijvx,'b');hold on
plot(F_wrijvy,'r');hold on
plot(abs(F_wrijvy+F_wrijvx),'b','Linewidth',2);hold on
plot(SkateLabdas.*sin(THETA_LS(2:end))','k');hold on
plot(SkateLabdas.*sin(THETA_RS(2:end))','--k');hold on
plot(-GForceLS(T2,1).*sin(THETA_LS)','m');hold on
plot(GForceRS(T2,1).*sin(THETA_RS)','--c');hold on
legend('FwrijvX','FwrijvY', 'F_wrijvy+F_wrijvx', 'LF','RF','LMF','RMF')

%% check fitdata
T1 = 1:1:length(POS(4,:));

SA = skate_Array(T2);

figure
a(1) = subplot(411)
plot(RSFn(T2),'r','Linewidth',2);hold on
plot(LSFn(T2),'b','Linewidth',2); hold on
grid minor
a(2)=subplot(412)
plot(POS(1,:),'k','Linewidth',2);hold on
for i = 1:length(POS(4,:));
    if SA(i)==1
        plot(T1(i),POS(4,i),'.b','Linewidth',2);hold on
%         plot(T1(i),rPOS(4,i),'.k','Linewidth',2);hold on
    else
        plot(T1(i),POS(7,i),'.r','Linewidth',2);hold on
%          plot(T1(i),rPOS(7,i),'.k','Linewidth',2);hold on
    end
end
grid minor;ylabel('X-position')
a(3)=subplot(413)
plot(POS(2,:),'k','Linewidth',2);hold on
for i = 1:length(POS(4,:));
    if SA(i)==1
        plot(T1(i),POS(5,i),'.b','Linewidth',2);hold on
%         plot(T1(i),rPOS(5,i),'.k','Linewidth',2);hold on
    else
        plot(T1(i),POS(8,i),'.r','Linewidth',2);hold on
%         plot(T1(i),rPOS(8,i),'.k','Linewidth',2);hold on
    end
end
grid minor;ylabel('Y-position')
a(4)=subplot(414)
plot(POS(3,:),'k','Linewidth',2);hold on
for i = 1:length(POS(4,:));
    if SA(i)==1
        plot(T1(i),POS(6,i),'.b','Linewidth',2);hold on
%         plot(T1(i),rPOS(6,i),'.k','Linewidth',2);hold on
    else
        plot(T1(i),POS(9,i),'.r','Linewidth',2);hold on
%         plot(T1(i),rPOS(9,i),'.k','Linewidth',2);hold on
    end
end
grid minor;ylabel('Z-position');linkaxes(a,'x');ylim([0 1.5])

figure
a(1) = subplot(411)
plot(RSFn(T2),'r','Linewidth',2);hold on
plot(LSFn(T2),'b','Linewidth',2); hold on
grid minor
a(2)=subplot(412)
plot(VEL(1,:),'k','Linewidth',2);hold on
for i = 1:length(POS(4,:));
    if SA(i)==1
        plot(T1(i),VEL(4,i),'.b','Linewidth',2);hold on
        plot(T1(i),rVEL(4,i),'.k','Linewidth',2);hold on
    else
        plot(T1(i),VEL(7,i),'.r','Linewidth',2);hold on
        plot(T1(i),rVEL(7,i),'.k','Linewidth',2);hold on
    end
end
grid minor;ylabel('X-velocity')
a(3)=subplot(413)
plot(VEL(2,:),'k','Linewidth',2);hold on
for i = 1:length(POS(4,:));
    if SA(i)==1
        plot(T1(i),VEL(5,i),'.b','Linewidth',2);hold on
        plot(T1(i),rVEL(5,i),'.k','Linewidth',2);hold on
    else
        plot(T1(i),VEL(8,i),'.r','Linewidth',2);hold on
        plot(T1(i),rVEL(8,i),'.k','Linewidth',2);hold on
    end
end
grid minor;ylabel('Y-velocity')
a(4)= subplot(414)
plot(VEL(3,:),'k','Linewidth',2);hold on
for i = 1:length(POS(4,:));
    if SA(i)==1
        plot(T1(i),VEL(6,i),'.b','Linewidth',2);hold on
        plot(T1(i),rVEL(6,i),'.k','Linewidth',2);hold on
    else
        plot(T1(i),VEL(9,i),'.r','Linewidth',2);hold on
        plot(T1(i),rVEL(9,i),'.k','Linewidth',2);hold on
    end
end
grid minor;ylabel('Z-velocity');
linkaxes(a,'x')

figure
a(1) = subplot(411)
plot(RSFn(T2),'r','Linewidth',2);hold on
plot(LSFn(T2),'b','Linewidth',2); hold on
grid minor
a(2)=subplot(412)
plot(ACC(1,:),'k','Linewidth',2);hold on
for i = 1:length(POS(4,:));
    if SA(i)==1
        plot(T1(i),ACC(4,i),'.b','Linewidth',2);hold on
        plot(T1(i),rACC(4,i),'.k','Linewidth',2);hold on
    else
        plot(T1(i),ACC(7,i),'.r','Linewidth',2);hold on
        plot(T1(i),rACC(7,i),'.k','Linewidth',2);hold on
    end
end
grid minor;ylabel('X-acceleratie')
a(3)=subplot(413)
plot(ACC(2,:),'k','Linewidth',2);hold on
for i = 1:length(POS(4,:));
    if SA(i)==1
        plot(T1(i),ACC(5,i),'.b','Linewidth',2);hold on
        plot(T1(i),rACC(5,i),'.k','Linewidth',2);hold on
    else
        plot(T1(i),ACC(8,i),'.r','Linewidth',2);hold on
        plot(T1(i),rACC(8,i),'.k','Linewidth',2);hold on
    end
end
grid minor;ylabel('Y-acceleratie')
a(4)= subplot(414)
plot(ACC(3,:),'k','Linewidth',2);hold on
for i = 1:length(POS(4,:));
    if SA(i)==1
        plot(T1(i),ACC(6,i),'.b','Linewidth',2);hold on
        plot(T1(i),rACC(6,i),'.k','Linewidth',2);hold on
    else
        plot(T1(i),ACC(9,i),'.r','Linewidth',2);hold on
         plot(T1(i),rACC(9,i),'.k','Linewidth',2);hold on
    end
end
grid minor;ylabel('Z-acceleratie');
linkaxes(a,'x')
hold off
