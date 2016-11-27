%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2016 Bart Grasmeijer
%       grasmeijerb
%
%       bartgrasmeijer@gmail.com
%       Meppel
%
%   This library is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This library is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this library.  If not, see <http://www.gnu.org/licenses/>.
%   --------------------------------------------------------------------


clearvars;
close all;
clc;

% Solving the shallow water equations in shallow basins and estuaries.
% Based on Friedrichs (2011), Chapter 3 in Contemporary Issues in Estuarine
% Physics.

%**************************************************************************
%**************************************************************************
%*              Paremeter settings
%**************************************************************************
%**************************************************************************


deltaT=60;              % time step in seconds. Choose appropriate time step yourself based on Courant number.
deltaX=1000;             % spatial step in meters
% Lbasin=2.15e5;           % Length of the basin or estuary in meters
Lbasin=290e3;           % Length of the basin or estuary in meters
% Lb=2e4;                  % e-folding length scale for width.
Lb = 200e4;                  % e-folding length scale for width.
Ld = 100e4;                  % e-folding length scale for depth.
% B0=5e3;                  % Width of the basin in meters at seaward side.
B0=500;                  % Width of the basin in meters at seaward side.
% H0=5.8;                  % Depth of basin.
H0=15;                  % Depth of basin.
M2amp=0.75;                 % Amplitude of M2 tide at seaward side.
discharge=1500;             % Constant river discharge at landward boundary.
% Cd=2.5e-3;               % Drag coefficient
Cd= 10e-3;

Do = 2500;              % Dispersion coefficient at mouth (range of 500 to 1500 )
So = 25;                % Salinity at mouth (kg/m3)
Tlws = 600;               % Time at which low water slack occurs (see tidal velocity plot)

g = 9.81;               % acceleration of gravity (m/s^2)

computetransport = false;

%**************************************************************************
%**************************************************************************
%*                 End of parameter setting
%**************************************************************************
%**************************************************************************


Tm2=12*3600+25*60;      % M2 tidal period in seconds
% time=0:deltaT:15*Tm2;    % time in seconds
time=0:deltaT:16*Tm2;    % time in seconds
Nt=length(time);

% %% load sediment fractions
% xls_sed_fname = 'c:\Users\grasmeijerb\Documents\C03041.002054_Nieuwe_Waterweg\01_data-in\20150225_sedimenteigenschappen\Sediment waterweg.xlsx';
% [xls_sed_num,xls_sed_txt,xls_sed_raw] = xlsread(xls_sed_fname,'Sheet1','c4:i15');
% 
% sed_km = xls_sed_num(:,1);
% sed_plut = xls_sed_num(:,3)./100;
% sed_psilt = xls_sed_num(:,4)./100;
% sed_psand = xls_sed_num(:,5)./100;
% sed_pgravel = xls_sed_num(:,6)./100;
% 
% sed_km_i = [1009:1035];
% sed_plut_i = interp1(sed_km,sed_plut,sed_km_i,'linear','extrap');
% sed_psilt_i = interp1(sed_km,sed_psilt,sed_km_i,'linear','extrap');
% sed_psand_i = interp1(sed_km,sed_psand,sed_km_i,'linear','extrap');
% sed_pgravel_i = interp1(sed_km,sed_pgravel,sed_km_i,'linear','extrap');
% 
% sed_plut_i = smooth(sed_km_i,sed_plut_i,.95,sed_km_i);
% sed_psilt_i = smooth(sed_km_i,sed_psilt_i,.95,sed_km_i);
% sed_psand_i = smooth(sed_km_i,sed_psand_i,.95,sed_km_i);
% sed_pgravel_i = smooth(sed_km_i,sed_pgravel_i,.95,sed_km_i);
% 
% sed_p_all_i = sed_plut_i+sed_psilt_i+sed_psand_i+sed_pgravel_i;
% 
% mykm = [0:26];
% mypmud = fliplr(sed_psilt_i);
% mypsand = fliplr(sed_psand_i);
% mypgravel = fliplr(sed_pgravel_i);
% 
% mykm(28) = 100;
% mypmud(28) = mypmud(27);
% mypsand(28) = mypsand(27);
% mypgravel(28) = mypgravel(27);
% 
% mypall = mypmud+mypsand+mypgravel
% 
% save('fractions','mykm','mypmud','mypsand','mypgravel');
% 
% %%
% figure;
% plot(mykm,mypmud);
% hold on;
% plot(mykm,mypsand);
% plot(mykm,mypgravel);
% legend('mud','sand','gravel')

%Define frequencies to be analysed. To determine the amplitudes and phase
%use the code you designed in the first Matlab exercise.

global wn

wn(1)=2*pi/Tm2;
wn(2)=2*wn(1);
wn(3)=3*wn(1);

x=0:deltaX:Lbasin;
Nx=length(x);

load('fractions.mat');
pmud = interp1(mykm.*1000,mypmud,x,'linear','extrap');
psand = interp1(mykm.*1000,mypsand,x,'linear','extrap');
pgravel = interp1(mykm.*1000,mypgravel,x,'linear','extrap');

pall = pmud+psand+pgravel;

disp('overwriting pmud, psand and pgravel...')
pmud = 0.30.*ones(size(pmud));
psand = 0.60.*ones(size(pmud));
pgravel = 0.10.*ones(size(pmud));

%%
figure;
plot(x./1000,pmud);
hold on;
plot(x./1000,psand);
plot(x./1000,pgravel);
legend('mud','sand','gravel')
xlim([0 70]);
%%

% B(1:Nx)=B0*exp(-x/Lb);
B(1:Nx)=B0;                % when basin width has to be constant.

%% initial depth
% H(1:Nx)=H0*exp(-x/Ld);
H(1:Nx)=H0;
H(x<20000) = 16;
H(x>=20000&x<30000) = 14.5;
H(x>=30000&x<35000) = 12;
H(x>=35000&x<40000) = 8;
% H(x>=40000) = 8;%8-(x(x>=40000)-40000).*0.5e-4;
H(x>=40000) = 8-(x(x>=40000)-40000).*0.3e-4;
% H(1:Nx)=H0;

%% new depth
H2(1:Nx)=H0;
H2(x<20000) = 16;
H2(x>=20000&x<30000) = 16;
H2(x>=30000&x<35000) = 12;
H2(x>=35000&x<40000) = 8;
% H(x>=40000) = 8;%8-(x(x>=40000)-40000).*0.5e-4;
H2(x>=40000) = 8-(x(x>=40000)-40000).*0.3e-4;


Z=zeros(Nx-1,Nt);           % Z points shifted half a grid point to the right with respect to Q points. we start with Q point, therefore 1 Z point less than Q points.
Q=zeros(Nx,Nt);
A=(B.*H)'*ones(1,Nt);       % A at Q points
P=B'*ones(1,Nt);            % Wetted perimeter at Q points.
Inertia=zeros(Nx,Nt);       % Initalize Inertia, Pressure Gradient and Friction for further analysis later on.
PG=zeros(Nx,Nt);
Fric=zeros(Nx,Nt);

myHm0(1:Nt) = 0.001;          % wave height;

salmean = zeros(Nx,1);
cmudmean_top = zeros(Nx,1);
cmudmean_bed = zeros(Nx,1);

%%
figure;
subplot(2,1,1);
plot(x,-H);
grid on;
xlim([0 50000]);
subplot(2,1,2)
plot(x,A(:,1));
xlim([0 50000]);
%%


% Boundary conditions
% Z(1,:)=M2amp*sin(2*pi*time/Tm2);          % prescribed water levels
% disp('only M2')
% Z1 = M2amp*sin(2*pi*time/Tm2);          % prescribed water levels

disp('tide from rws...')
load('c:\Users\grasmeijerb\Documents\anaconda\CoastalEngineeringTools\TidalChannel\tidal_predictions_from_rws_getij_nl.mat');
mytime = (tide.datenum-tide.datenum(1)).*24*3600;
myz = tide.zwl;
Z1 = interp1(mytime,myz,time);

Q(Nx,:)=-discharge;                      % river discharge; most often river discharge is taken zero.

[Z2,U2,A2,P2] = TidesInChannel('H',H,'Z1',Z1,'time',time);
[Z3,U3,A3,P3] = TidesInChannel('H',H2,'Z1',Z1,'time',time);


%%
ix = find(x==10000);
figure;
plot(time,Z2(ix,:));
hold on;
plot(time,Z3(ix,:));
plot(time,Z1);
legend('ori','new','Z1')
%%

courant=sqrt(9.8*max(H))*deltaT/deltaX;

% For numerical part, follow thesis of Speer (1984). Staggered grid. Z points shifted half deltaX to
% the right of U points:
% Q1 Z1 Q2 Z2 ........ ZN QN+1

% solve Bdz/dt + dQ/dx=0
% solve dQ/dt + d/dx Q^2/A = -gA dZ/dx - Cd Q|Q|P/A*A

% Z is water level, B=estuary width, Q is discharge, A is cross-sectional
% area, P is wetted perimeter (~channel width when H/B is small)

% First start simple: rectangular basin, no advection, no river flow.
% Numerical scheme from Speer (1984).

% for pt=1:Nt-1
%     for px=2:Nx-1
%         Z(px,pt+1)=Z(px,pt)-(deltaT/(0.5*(B(px)+B(px+1))))*(Q(px+1,pt)-Q(px,pt))/deltaX;
%     end
%     for px=2:Nx-1
%         A(px,pt+1)=B(px)*(H(px)+0.5*Z(px,pt+1)+0.5*Z(px-1,pt+1));           % A at Q points
%         P(px,pt+1)=B(px)+2*H(px)+Z(px,pt+1) +Z(px-1,pt+1);                  % P at Q points
%     end
%     for px=2:Nx-1
%         Q(px,pt+1)=Q(px,pt) ...                                            % Inertia.
%             -9.81*A(px,pt+1)*(deltaT/deltaX)*(Z(px,pt+1)-Z(px-1,pt+1)) ...  % Pressure gradient
%             -Cd*deltaT*abs(Q(px,pt))*Q(px,pt)*P(px,pt)/(A(px,pt)*A(px,pt)); % Friction
%         Inertia(px,pt+1)=(Q(px,pt+1)-Q(px,pt))/deltaT;
%         PG(px,pt+1)=-9.81*A(px,pt+1)*(1/deltaX)*(Z(px,pt+1)-Z(px-1,pt+1));
%         Fric(px,pt+1)=-Cd*abs(Q(px,pt))*Q(px,pt)*P(px,pt)/(A(px,pt)*A(px,pt));
%     end
%     Q(1,pt+1)=Q(2,pt+1)+B(1)*deltaX*(Z(1,pt+1)-Z(1,pt))/deltaT;
% end
% 
% U=Q./A;         % Flow velocity in m/s

Upeak = maxmax(U2);
Upeakx = max(U2(:,1000:end),[],2);
Upeak3 = maxmax(U3);
Upeakx3 = max(U3(:,1000:end),[],2);


Le = Upeak.*Tm2./3.18;    % Tidal excursion (m)

Uriver = discharge./(B0.*H0);     % (+=landward; - = seaward)

if Uriver>0
    Lsmin = Do./(0.5.*abs(Uriver));
    Lsmax = Le+Lsmin;
else
    Lsmin = Le;
    Lsmax = Le;
end


% Analyse last tidal period only. For example determine M2, M4, M6 and mean
% of water level. Design you own code here. I used my code (harmfit). You
% can determine HW, LW, moments of LW and HW, propagation speed of LW wave
% and HW wave...

multiWaitbar('CloseAll');

Nsteps=floor(Tm2/deltaT);
for px=1:Nx-1
    multiWaitbar('Looping through x...',px./(Nx-1));
%     coefin=[0.1, 1, 0.2, 0.1, 1, 0.2, 0.1];
%     coefout=nlinfit(time(end-Nsteps:end),Z(px,end-Nsteps:end),@harmfit,coefin);
%     Z0(px)=coefout(1);
%     ZM2(px)=sqrt(coefout(2)^2 + coefout(5)^2);
%     ZM4(px)=sqrt(coefout(3)^2 + coefout(6)^2);
%     ZM6(px)=sqrt(coefout(4)^2 + coefout(7)^2);
%     phaseZM2(px)=atan2(coefout(2),coefout(5));
%     phaseZM4(px)=atan2(coefout(3),coefout(6));
%     phaseZM6(px)=atan2(coefout(4),coefout(7));
%     coefin=[0.1, 1, 0.2, 0.1, 1, 0.2, 0.1];
%     coefout=nlinfit(time(end-Nsteps:end),U(px,end-Nsteps:end),@harmfit,coefin);
%     U0(px)=coefout(1);
%     UM2(px)=sqrt(coefout(2)^2 + coefout(5)^2);
%     UM4(px)=sqrt(coefout(3)^2 + coefout(6)^2);
%     UM6(px)=sqrt(coefout(4)^2 + coefout(7)^2);
%     phaseUM2(px)=atan2(coefout(2),coefout(5));
%     phaseUM4(px)=atan2(coefout(3),coefout(6));
%     phaseUM6(px)=atan2(coefout(4),coefout(7));
    
    %     disp('computing horizontal salinity gradient...')
    h = (H(px)+Z(px,:));                                                    % depth (m)
    if x(px)<0
        sal = So;
    else if x(px)-0.5*Le+0.5*Le*cos(2.*pi./Tm2.*(time-Tlws))-Lsmin<0
            sal = So.*((1-(1./(Lsmin+500.*h(px))).*(x(px)+150.*h(px)-0.5.*Le+0.5.*Le.*cos(2.*pi./Tm2.*(time-Tlws)))).^2);
            sal(sal>So) = So;
        else
            sal = 0;
        end
    end
    
    if computetransport
        
        %% initialize tsand input
        mydatenum = time./(24.*3600);
        a = 0.10;                                                               % reference level (m)
        D50 = 0.00025;
        D90 = 0.0020;
        nrofsigmalevels = 20;
        zbedt = 0-h;                                                            % bed level relative to reference
        salt = 25.*ones(size(h));                                               % salinity (ppt);
        Up = U(px,:);                                                           % velocity (m/s)
        kscvel = 1;                                                             % current-related roughness (effective roughness) for velocity profile (m)
        
%         psand = 0.90;
%         pmud = 0.10;
%         pgravel = 0;
        
        
        zt = NaN(length(mydatenum),nrofsigmalevels);
        zi = 1:1:nrofsigmalevels;                                                   % sigma levels
        disp('creating z levels (sigma)...')
        zt(:,1) = a;
        blup1 = a.*(h./a);
        blup2 = (zi(2:end)./(length(zi)));
        blup3 = blup1'.*blup2;
        zt(:,2:end) = blup3;
%         zt(:,2:end) = a.*(h./a).^(zi(2:end)'./(length(zi)));
%         for i = 1:length(h)
%             multiWaitbar( 'creating z (sigma)', i/length(h) );
%             zt(i,1) = a;
%             zt(i,2:end) = a.*(h(i)./a).^(zi(2:end)./(length(zi)));
%             dz(i,:) = gradient(zt(i,:));
%         end
        
        disp('Creating velocity profiles...')
        C = 18.*log10(12.*h./kscvel);
        fc = 8.*g./C.^2;
        for i = 1:length(Up)
            U1t(i,:) = Up(i)./(-1+log(h(i)/(0.033*kscvel))).*log(zt(i,:)./(0.033.*kscvel));
        end
        
        zvelt = -zt;
        V1t = zeros(size(U1t));
        
        disp('running tsand...')
        [qtotsandx,qtotsandy, qssandxVRadjust, qssandyVRadjust, z, Uz, Vz, csand, qssandx, qssandy, cmud, qmudxtot] = tsand('Times',mydatenum','Hs',myHm0','U1',U1t,...
            'V1',V1t,'zvel',zvelt,...
            'zbedt',zbedt','d',h',...
            'salt',sal','D50',D50,'D90',D90,...
            'psand',psand(px),'pmud',pmud(px),'pgravel',pgravel(px),...
            'wsmud',0.0002);
        
        %     figure;
        %     plot(time,U(px,:))
        %
        %     figure;
        %     plot(time,sal);
        %     grid on;
        %
        %     figure;
        %     plot(time,qtotsandx);
    end
    salmean(px) = mean(sal);
    if computetransport
%     cmudmean_top(px) = mean(cmud(end-3:end));
%     cmudmean_bed(px) = mean(cmud(1:3));
    end
end


TP=trapz(time(end-Nsteps:end),abs(Q(2,end-Nsteps:end)))/2;      % this calculates the the tidal prism.


%% WITH NEW DEPTH

Nsteps=floor(Tm2/deltaT);
for px=1:Nx-1
    multiWaitbar('Looping through x...',px./(Nx-1));
%     coefin3=[0.1, 1, 0.2, 0.1, 1, 0.2, 0.1];
%     coefout3=nlinfit(time(end-Nsteps:end),Z3(px,end-Nsteps:end),@harmfit,coefin3);
%     Z03(px)=coefout3(1);
%     ZM23(px)=sqrt(coefout3(2)^2 + coefout3(5)^2);
%     ZM43(px)=sqrt(coefout3(3)^2 + coefout3(6)^2);
%     ZM63(px)=sqrt(coefout3(4)^2 + coefout3(7)^2);
%     phaseZM23(px)=atan2(coefout3(2),coefout3(5));
%     phaseZM43(px)=atan2(coefout3(3),coefout3(6));
%     phaseZM63(px)=atan2(coefout3(4),coefout3(7));
%     coefin3=[0.1, 1, 0.2, 0.1, 1, 0.2, 0.1];
%     coefout3=nlinfit(time(end-Nsteps:end),U3(px,end-Nsteps:end),@harmfit,coefin3);
%     U03(px)=coefout3(1);
%     UM23(px)=sqrt(coefout3(2)^2 + coefout3(5)^2);
%     UM43(px)=sqrt(coefout3(3)^2 + coefout3(6)^2);
%     UM63(px)=sqrt(coefout3(4)^2 + coefout3(7)^2);
%     phaseUM23(px)=atan2(coefout3(2),coefout3(5));
%     phaseUM43(px)=atan2(coefout3(3),coefout3(6));
%     phaseUM63(px)=atan2(coefout3(4),coefout3(7));
    
    %     disp('computing horizontal salinity gradient...')
    h = (H(px)+Z(px,:));                                                    % depth (m)
    if x(px)<0
        sal = So;
    else if x(px)-0.5*Le+0.5*Le*cos(2.*pi./Tm2.*(time-Tlws))-Lsmin<0
            sal = So.*((1-(1./(Lsmin+500.*h(px))).*(x(px)+150.*h(px)-0.5.*Le+0.5.*Le.*cos(2.*pi./Tm2.*(time-Tlws)))).^2);
            sal(sal>So) = So;
        else
            sal = 0;
        end
    end
    
    if computetransport
        
        %% initialize tsand input
        mydatenum = time./(24.*3600);
        a = 0.10;                                                               % reference level (m)
        D50 = 0.00025;
        D90 = 0.0020;
        nrofsigmalevels = 20;
        zbedt = 0-h;                                                            % bed level relative to reference
        salt = 25.*ones(size(h));                                               % salinity (ppt);
        Up = U(px,:);                                                           % velocity (m/s)
        kscvel = 1;                                                             % current-related roughness (effective roughness) for velocity profile (m)
        
%         psand = 0.98;
%         pmud = 0.02;
%         pgravel = 0;
        
        zt(:,1) = a;
        blup1 = a.*(h./a);
        blup2 = (zi(2:end)./(length(zi)));
        blup3 = blup1'.*blup2;
        zt(:,2:end) = blup3;

        
%         zt = NaN(length(mydatenum),nrofsigmalevels);
%         zi = 1:1:nrofsigmalevels;                                                   % sigma levels
%         disp('creating z levels (sigma)...')
%         for i = 1:length(h)
%             multiWaitbar( 'creating z (sigma)', i/length(h) );
%             zt(i,1) = a;
%             zt(i,2:end) = a.*(h(i)./a).^(zi(2:end)./(length(zi)));
%             dz(i,:) = gradient(zt(i,:));
%         end
        
        disp('Creating velocity profiles...')
        C = 18.*log10(12.*h./kscvel);
        fc = 8.*g./C.^2;
        for i = 1:length(Up)
            U1t(i,:) = Up(i)./(-1+log(h(i)/(0.033*kscvel))).*log(zt(i,:)./(0.033.*kscvel));
        end
        
        zvelt = -zt;
        V1t = zeros(size(U1t));
        
        disp('running tsand...')
        [qtotsandx,qtotsandy, qssandxVRadjust, qssandyVRadjust, z, Uz, Vz, csand, qssandx, qssandy, cmud, qmudxtot] = tsand('Times',mydatenum','Hs',myHm0','U1',U1t,...
            'V1',V1t,'zvel',zvelt,...
            'zbedt',zbedt','d',h',...
            'salt',sal','D50',D50,'D90',D90,...
            'psand',psand(px),'pmud',pmud(px),'pgravel',pgravel(px),...
            'wsmud',0.0002);
        
        %     figure;
        %     plot(time,U(px,:))
        %
        %     figure;
        %     plot(time,sal);
        %     grid on;
        %
        %     figure;
        %     plot(time,qtotsandx);
    end
    salmean(px) = mean(sal);
    
end


TP=trapz(time(end-Nsteps:end),abs(Q(2,end-Nsteps:end)))/2;      % this calculates the the tidal prism.

%%
disp('Analaysing tides...')
for px = 1:Nx-1
    multiWaitbar('Analysing tides...',px./(Nx-1));
    myzwl = Z2(px,:);
    [NAME,FREQ,TIDECON,XOUT] = t_tide(myzwl','interval',1/60,'latitude',51.45,'output','none');
    for j = 1:length(NAME)
        if strcmp(NAME(j,:),'M2  ')
            M2_obs.amp(px) = TIDECON(j,1);
            M2_obs.amp95(px) = TIDECON(j,2);
            M2_obs.phase(px) = TIDECON(j,3);
        end
        if strcmp(NAME(j,:),'M4  ')
            M4_obs.amp(px) = TIDECON(j,1);
            M4_obs_amp95(px) = TIDECON(j,2);
            M4_obs.phase(px) = TIDECON(j,3);
        end
        if strcmp(NAME(j,:),'M6  ')
            M6_obs.amp(px) = TIDECON(j,1);
            M6_obs.amp95(px) = TIDECON(j,2);
            M6_obs.phase(px) = TIDECON(j,3);
        end
    end
    myzwl = Z3(px,:);
    [NAME,FREQ,TIDECON,XOUT] = t_tide(myzwl','interval',1/60,'latitude',51.45,'output','none');
    for j = 1:length(NAME)
        if strcmp(NAME(j,:),'M2  ')
            M2_obs3.amp(px) = TIDECON(j,1);
            M2_obs3.amp95(px) = TIDECON(j,2);
            M2_obs3.phase(px) = TIDECON(j,3);
        end
        if strcmp(NAME(j,:),'M4  ')
            M4_obs3.amp(px) = TIDECON(j,1);
            M4_obs3_amp95(px) = TIDECON(j,2);
            M4_obs3.phase(px) = TIDECON(j,3);
        end
        if strcmp(NAME(j,:),'M6  ')
            M6_obs3.amp(px) = TIDECON(j,1);
            M6_obs3.amp95(px) = TIDECON(j,2);
            M6_obs3.phase(px) = TIDECON(j,3);
        end
    end
    
    
end

figure;plot(x(1:end-1),M2_obs.amp);hold on; plot(x(1:end-1),M2_obs3.amp)



%%
close all;
xmaxplot = 70;
figure;
figsize = [0 0 5.9 7];     % set figure size to 5.9 inch x 2.95 inch = 15 cm x 7.5 cm
set(gcf,'PaperOrientation','portrait','PaperUnits','inches' ,'PaperPosition',figsize);
subplot(8,1,1);
plot(x./1000,0.5.*B,'b-')
hold on;
plot(x./1000,-0.5.*B,'b-')
grid on;
xlim([0 xmaxplot]);
ylim([-300 300]);
set(gca,'fontsize',7);
title('Channel layout')
subplot(8,1,2);
plot(x./1000,-H,'b-')
hold on;
plot(x./1000,-H2,'r--')
grid on;
legend('present','deepened','location','best')
xlim([0 xmaxplot]);
set(gca,'fontsize',7);
title('Bed level');
subplot(8,1,3);
% plot(x(1:end-1)./1000,Z0)
hold on;
grid on;
xlim([0 xmaxplot]);
set(gca,'fontsize',7);
title('Mean water level');
subplot(8,1,4);
plot(x(1:end-1)./1000,M2_obs.amp,'b-');
hold on;
plot(x(1:end-1)./1000,M2_obs3.amp,'r--');
set(gca,'fontsize',7);
leg = legend('M2 present','M2 deepened','location','southwest');
set(leg,'fontsize',5);
plot(x(1:end-1)./1000,M4_obs.amp,'b-');
plot(x(1:end-1)./1000,M4_obs3.amp,'r--');
grid on;
xlim([0 xmaxplot]);
title('Amplitude M2 and M4 tide');
subplot(8,1,5);
plot(x(1:end-1)./1000,M2_obs3.amp-M2_obs.amp);
hold on;
plot(x(1:end-1)./1000,M4_obs3.amp-M4_obs.amp);
set(gca,'fontsize',7);
ylim([-0.02 0.02])
leg = legend('effect M2','effect M4','location','southwest');
set(leg,'fontsize',5);
title('Effect of deepening on tidal amplitude');
grid on;
xlim([0 xmaxplot]);
subplot(8,1,6);
plot(x(1:end-1)./1000,2.*M2_obs.phase-M4_obs.phase,'b-');
hold on;
plot(x(1:end-1)./1000,2.*M2_obs3.phase-M4_obs3.phase,'r--');
set(gca,'fontsize',7);
leg = legend('present','deepened','location','southwest');
set(leg,'fontsize',5);
title('Effect of deepening on tidal asymmetry (2\phi_{M2}-\phi_{M4})');
grid on;
xlim([0 xmaxplot]);
subplot(8,1,7);
plot(x./1000,Upeakx,'b-');
hold on;
plot(x./1000,Upeakx3,'r--');
xlim([0 xmaxplot]);
set(gca,'fontsize',7);
title('Upeak');
subplot(8,1,8);
plot(x./1000,salmean);
hold on;
grid on;
xlim([0 xmaxplot]);
set(gca,'fontsize',7);
title('Salinity');
xlabel('distance along Rotterdam Waterway (km)');
text(1,0,['\copyright Utrecht University ',datestr(now,10)],'fontsize',5,'rotation',90,'unit','n','ver','t');  % add ARCADIS copyright
annotation('textbox',[1,0.0,0,0],'string',[addslash([mfilename])],'fontsize',4,'horizontalalignment','right','verticalalignment','baseline','color',[0.5 0.5 0.5]);  % add script name
print('-dpng','-r300',['TidalChannel_Le',num2str(Lb)])  % print figure at 300 dpi


figure;plot(x./1000,1000.*cmudmean_bed);xlim([0 xmaxplot]);

