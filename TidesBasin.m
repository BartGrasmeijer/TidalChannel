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


deltaT=100;              % time step in seconds. Choose appropriate time step yourself based on Courant number. 
deltaX=1000;             % spatial step in meters
% Lbasin=2.15e5;           % Length of the basin or estuary in meters
Lbasin=50e3;           % Length of the basin or estuary in meters
Lb=4e4;                  % e-folding length scale for width.
B0=5e3;                  % Width of the basin in meters at seaward side.
H0=5.8;                  % Depth of basin.
M2amp=1;                 % Amplitude of M2 tide at seaward side.
discharge=0;             % Constant river discharge at landward boundary. 
Cd=2.5e-3;               % Drag coefficient

Do = 1500;              % Dispersion coefficient at mouth (range of 500 to 1500 )
So = 25;                % Salinity at mouth (kg/m3)
Tlws = 0;               % Time at which low water slack occurs (see tidal velocity plot)

g = 9.81;               % acceleration of gravity (m/s^2)

%**************************************************************************
%**************************************************************************
%*                 End of parameter setting
%**************************************************************************
%**************************************************************************


Tm2=12*3600+25*60;      % M2 tidal period in seconds
time=0:deltaT:15*Tm2;    % time in seconds
Nt=length(time);



%Define frequencies to be analysed. To determine the amplitudes and phase
%use the code you designed in the first Matlab exercise. 

global wn

wn(1)=2*pi/Tm2;
wn(2)=2*wn(1);
wn(3)=3*wn(1);

x=0:deltaX:Lbasin;
Nx=length(x);

B(1:Nx)=B0*exp(-x/Lb);
%B(1:Nx)=B0;                % when basin width has to be constant.
H(1:Nx)=H0;


Z=zeros(Nx-1,Nt);           % Z points shifted half a grid point to the right with respect to Q points. we start with Q point, therefore 1 Z point less than Q points.      
Q=zeros(Nx,Nt);
A=(B.*H)'*ones(1,Nt);       % A at Q points
P=B'*ones(1,Nt);            % Wetted perimeter at Q points.
Inertia=zeros(Nx,Nt);       % Initalize Inertia, Pressure Gradient and Friction for further analysis later on.
PG=zeros(Nx,Nt);
Fric=zeros(Nx,Nt);

myHm0(1:Nt) = 0.001;          % wave height;


% Boundary conditions
Z(1,:)=M2amp*sin(2*pi*time/Tm2);          % prescribed water levels
Q(Nx,:)=-discharge;                      % river discharge; most often river discharge is taken zero.

%%
figure;
plot(time,Z(1,:));
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

for pt=1:Nt-1
    for px=2:Nx-1
        Z(px,pt+1)=Z(px,pt)-(deltaT/(0.5*(B(px)+B(px+1))))*(Q(px+1,pt)-Q(px,pt))/deltaX;
    end
    for px=2:Nx-1
        A(px,pt+1)=B(px)*(H(px)+0.5*Z(px,pt+1)+0.5*Z(px-1,pt+1));           % A at Q points
        P(px,pt+1)=B(px)+2*H(px)+Z(px,pt+1) +Z(px-1,pt+1);                  % P at Q points
    end
    for px=2:Nx-1
        Q(px,pt+1)=Q(px,pt) ...                                            % Inertia.
            -9.81*A(px,pt+1)*(deltaT/deltaX)*(Z(px,pt+1)-Z(px-1,pt+1)) ...  % Pressure gradient
            -Cd*deltaT*abs(Q(px,pt))*Q(px,pt)*P(px,pt)/(A(px,pt)*A(px,pt)); % Friction
        Inertia(px,pt+1)=(Q(px,pt+1)-Q(px,pt))/deltaT;
        PG(px,pt+1)=-9.81*A(px,pt+1)*(1/deltaX)*(Z(px,pt+1)-Z(px-1,pt+1));
        Fric(px,pt+1)=-Cd*abs(Q(px,pt))*Q(px,pt)*P(px,pt)/(A(px,pt)*A(px,pt));
    end
    Q(1,pt+1)=Q(2,pt+1)+B(1)*deltaX*(Z(1,pt+1)-Z(1,pt))/deltaT;
end

U=Q./A;         % Flow velocity in m/s

Upeak = maxmax(U);

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
    coefin=[0.1, 1, 0.2, 0.1, 1, 0.2, 0.1];
    coefout=nlinfit(time(end-Nsteps:end),Z(px,end-Nsteps:end),@harmfit,coefin);
    Z0(px)=coefout(1);
    ZM2(px)=sqrt(coefout(2)^2 + coefout(5)^2);
    ZM4(px)=sqrt(coefout(3)^2 + coefout(6)^2);
    ZM6(px)=sqrt(coefout(4)^2 + coefout(7)^2);
    phaseZM2(px)=atan2(coefout(2),coefout(5));
    phaseZM4(px)=atan2(coefout(3),coefout(6));
    phaseZM6(px)=atan2(coefout(4),coefout(7));
    coefin=[0.1, 1, 0.2, 0.1, 1, 0.2, 0.1];
    coefout=nlinfit(time(end-Nsteps:end),U(px,end-Nsteps:end),@harmfit,coefin);
    U0(px)=coefout(1);
    UM2(px)=sqrt(coefout(2)^2 + coefout(5)^2);
    UM4(px)=sqrt(coefout(3)^2 + coefout(6)^2);
    UM6(px)=sqrt(coefout(4)^2 + coefout(7)^2);
    phaseUM2(px)=atan2(coefout(2),coefout(5));
    phaseUM4(px)=atan2(coefout(3),coefout(6));
    phaseUM6(px)=atan2(coefout(4),coefout(7));
    

    %% initialize tsand input
    mydatenum = time./(24.*3600);
    a = 0.10;                                                                  % reference level (m)
    D50 = 0.00025;
    D90 = 0.0020;
    nrofsigmalevels = 20;
    h = (H(px)+Z(px,:));                                 % depth (m)
    zbedt = 0-h;                                                               % bed level relative to reference
    salt = 25.*ones(size(h));                                                  % salinity (ppt);
    Up = U(px,:);                                             % velocity (m/s)
    kscvel = 1;                                                                 % current-related roughness (effective roughness) for velocity profile (m)
    
    psand = 0.98;
    pmud = 0.02;
    pgravel = 0;
    

    C = 18.*log10(12.*h./kscvel);
    fc = 8.*g./C.^2;
    if x(px)<0
        sal = So;
    else if x(px)-0.5*Le+0.5*Le*cos(2.*pi./Tm2.*(time-Tlws))-Lsmin<0
            sal = So.*((1-(1./(Lsmin+500.*h(px))).*(x(px)+150.*h(px)-0.5.*Le+0.5.*Le.*cos(2.*pi./Tm2.*(time-Tlws)))).^2);
            sal(sal>So) = So;
        else
            sal = 0;
        end
    end
    
    zt = NaN(length(mydatenum),nrofsigmalevels);
    zi = 1:1:nrofsigmalevels;                                                   % sigma levels
    disp('creating z...')
    for i = 1:length(h)
        multiWaitbar( 'creating z', i/length(h) );
        zt(i,1) = a;
        zt(i,2:end) = a.*(h(i)./a).^(zi(2:end)./(length(zi)));
        dz(i,:) = gradient(zt(i,:));
    end
    
    disp('Creating velocity profiles...')
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
        'psand',psand,'pmud',pmud,'pgravel',pgravel);
    
%     figure;
%     plot(time,U(px,:))
%     
%     figure;
%     plot(time,sal);
%     grid on;
%     
%     figure;
%     plot(time,qtotsandx);

    salmean(px) = mean(sal);

end

    
TP=trapz(time(end-Nsteps:end),abs(Q(2,end-Nsteps:end)))/2;      % this calculates the the tidal prism.

%%
figure;
figsize = [0 0 5.9 3];     % set figure size to 5.9 inch x 2.95 inch = 15 cm x 7.5 cm
set(gcf,'PaperOrientation','portrait','PaperUnits','inches' ,'PaperPosition',figsize);
subplot(4,1,1);
plot(x,0.5.*B,'b-')
hold on;
plot(x,-0.5.*B,'b-')
grid on;
title('Channel layout')
subplot(4,1,2);
plot(x(1:end-1),Z0)
hold on;
grid on;
title('Mean water level')
subplot(4,1,3);
plot(x(1:end-1),ZM2)
hold on;
grid on;
title('Amplitude M2 tide')
subplot(4,1,4);
plot(x(1:end-1),salmean)
hold on;
grid on;
title('Salinity')