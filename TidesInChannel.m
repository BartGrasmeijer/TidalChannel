function [Z1,U,A,P] = TidesInChannel(varargin)
%TIDESINCHANNEL is a 1D model that computes tidal propagation in a tidal
% channel
%
%   More detailed description goes here.
%
%   Syntax:
%   [Z2,U2,A2,P2] = TidesInChannel('H',H,'Z1',Z1);
%
%   Input: For <keyword,value> pairs call Untitled() without arguments.
%   varargin  =
%
%   Output:
%   varargout =
%
%   Example
%   [Z2,U2,A2,P2] = TidesInChannel('H',H,'Z1',Z1);
%
%   See also t_tide

%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2016 Bart Grasmeijer
%       Bart Grasmeijer
%
%       bartgrasmeijer@gmail.com
%
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

% This tool is part of <a href="http://www.OpenEarth.eu">OpenEarthTools</a>.
% OpenEarthTools is an online collaboration to share and manage data and
% programming tools in an open source, version controlled environment.
% Sign up to recieve regular updates of this function, and to contribute
% your own tools.

%% Version <http://svnbook.red-bean.com/en/1.5/svn.advanced.props.special.keywords.html>
% Created: 22 Nov 2016
% Created with Matlab version: 9.1.0.441655 (R2016b)

% $Id: $
% $Date: $
% $Author: $
% $Revision: $
% $HeadURL: $
% $Keywords: $

%%
OPT.deltaT=60;              % time step in seconds. Choose appropriate time step yourself based on Courant number.
OPT.deltaX=1000;             % spatial step in meters
OPT.Lbasin=290e3;           % Length of the basin or estuary in meters
OPT.Lb = 200e4;                  % e-folding length scale for width.
OPT.Ld = 100e4;                  % e-folding length scale for depth.
OPT.B0=500;                  % Width of the basin in meters at seaward side.
OPT.H0=15;                  % Depth of basin.
OPT.M2amp=0.75;                 % Amplitude of M2 tide at seaward side.
discharge=1500;             % Constant river discharge at landward boundary.
OPT.Cd= 10e-3;

OPT.Do = 2500;              % Dispersion coefficient at mouth (range of 500 to 1500 )
OPT.So = 25;                % Salinity at mouth (kg/m3)
OPT.Tlws = 600;               % Time at which low water slack occurs (see tidal velocity plot)
OPT.Tm2=12*3600+25*60;      % M2 tidal period in seconds
OPT.time=0:OPT.deltaT:16*OPT.Tm2;    % time in seconds

OPT.Z1(1,:)=OPT.M2amp*sin(2*pi*OPT.time/OPT.Tm2);          % prescribed water levels

x=0:OPT.deltaX:OPT.Lbasin;
Nx = length(x);

% B(1:Nx)=B0*exp(-x/Lb);
OPT.B(1:Nx)=OPT.B0;                % when basin width has to be constant.
% H(1:Nx)=H0*exp(-x/Ld);
OPT.H(1:Nx)=OPT.H0;
OPT.H(x<20000) = 16;
OPT.H(x>=20000&x<30000) = 14.5;
OPT.H(x>=30000&x<35000) = 12;
OPT.H(x>=35000&x<40000) = 8;
OPT.H(x>=40000) = 8-(x(x>=40000)-40000).*0.3e-4;

% return defaults (aka introspection)
% if nargin==0;
%     varargout = {OPT};
%     return
% end
% overwrite defaults with user arguments
OPT = setproperty(OPT, varargin);
%% code

g = 9.81;

%% start running

global wn

wn(1)=2*pi/OPT.Tm2;
wn(2)=2*wn(1);
wn(3)=3*wn(1);

time = OPT.time;
deltaT = OPT.deltaT;
deltaX = OPT.deltaX;
Lbasin = OPT.Lbasin;
B = OPT.B;
H = OPT.H;
B0 = OPT.B0;
H0 = OPT.H0;
Tm2 = OPT.Tm2;
Cd = OPT.Cd;
Do = OPT.Do;
So = OPT.So;


x=0:deltaX:Lbasin;
Nx=length(x);
Nt=length(time);

Z1=zeros(Nx-1,Nt);           % Z points shifted half a grid point to the right with respect to Q points. we start with Q point, therefore 1 Z point less than Q points.
Q=zeros(Nx,Nt);
A=(OPT.B.*OPT.H)'*ones(1,Nt);       % A at Q points
P=OPT.B'*ones(1,Nt);            % Wetted perimeter at Q points.
Inertia=zeros(Nx,Nt);       % Initalize Inertia, Pressure Gradient and Friction for further analysis later on.
PG=zeros(Nx,Nt);
Fric=zeros(Nx,Nt);

myHm0(1:Nt) = 0.001;          % wave height;
Z1(1,:)=OPT.Z1(1,:);

courant=sqrt(9.8*max(OPT.H))*OPT.deltaT/OPT.deltaX;

for pt=1:Nt-1
    for px=2:Nx-1
        Z1(px,pt+1)=Z1(px,pt)-(deltaT/(0.5*(B(px)+B(px+1))))*(Q(px+1,pt)-Q(px,pt))/deltaX;
    end
    for px=2:Nx-1
        A(px,pt+1)=B(px)*(H(px)+0.5*Z1(px,pt+1)+0.5*Z1(px-1,pt+1));           % A at Q points
        P(px,pt+1)=B(px)+2*H(px)+Z1(px,pt+1) +Z1(px-1,pt+1);                  % P at Q points
    end
    for px=2:Nx-1
        Q(px,pt+1)=Q(px,pt) ...                                            % Inertia.
            -9.81*A(px,pt+1)*(deltaT/deltaX)*(Z1(px,pt+1)-Z1(px-1,pt+1)) ...  % Pressure gradient
            -Cd*deltaT*abs(Q(px,pt))*Q(px,pt)*P(px,pt)/(A(px,pt)*A(px,pt)); % Friction
        Inertia(px,pt+1)=(Q(px,pt+1)-Q(px,pt))/deltaT;
        PG(px,pt+1)=-9.81*A(px,pt+1)*(1/deltaX)*(Z1(px,pt+1)-Z1(px-1,pt+1));
        Fric(px,pt+1)=-Cd*abs(Q(px,pt))*Q(px,pt)*P(px,pt)/(A(px,pt)*A(px,pt));
    end
    Q(1,pt+1)=Q(2,pt+1)+B(1)*deltaX*(Z1(1,pt+1)-Z1(1,pt))/deltaT;
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

% figure;
% plot(time,Z(1,:));
% hold on;
% plot(time,Z2(1,:));
