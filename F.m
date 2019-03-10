%% Linearized INS 15-State Error Equations

function dxdot = F(t,dx) %#ok<INUSL>
% ========================================================================
% This function represents the 15-State linearized INS state error equation
%
% Linearized States (dx) in NED (North-East-Down) Navigation frame are
% defined as below.
%
% Position Errors dP    = [del_lat; del_long; ; del_alt];
%                       = [dphi; dlambda; dh];
%
% Velocity Errors dVn   = [del_vN; del_vE; del_vD];
%                       = [dvN; dvE; dvD];
%
% Attitude Errors dPSIn = [del_psiN; del_psiE; del_psiD];
%                       = [dpN; dpE; dpD];
%
% NOTE: RHO = [psiN; psiE; psiD] are positively defined small-angle
% rotations about the navigation-frame axis to align the navigation frame
% with the computed navigation-frame. epsN and epsD are referred to as tilt
% errors and epsD is referred to as the heading, yaw or azimuth error.
%
% Additional Forcing (6-States) Can be provided by Accelerometer Errors
% (dae) Gyro Errors (dge)
%
% References:
% [1] W. S. Widnall and P. A. Grundy, Inertial Navigation System Error
% Models, Technical Report TR-03-73, Intermetrics, Cambridge, MA, 1973.
%
% [2] Farrell, J. and Barth, M., 1999. The global positioning system and
% inertial navigation (Vol. 61). New York: Mcgraw-hill.
% ========================================================================

% Author: Jyot R. Buch
% Ph.D. Student in Aerospace Engineering and Mechanics
% University of Minnesota, Twin Cities
% Copyright (c) 2019, All Rights Reserved.
% Course Project : Navigation and Attitude Determination Systems

% Vehicle is nominally located at 45 deg North Latitude and 45 deg East
% Longitude and sea-level i.e. h = 0 Altitude
phi = deg2rad(45);
h = 0;

% Earth Radii of curvature
[Rphi,Rlambda] = earthrad(phi);
Re = sqrt(Rlambda*Rphi); % m

% Gravity Model
g = norm(glocal(phi,h));

% Earth Rate
omega_n_ie = earthrate(phi);
omega_ie = norm(omega_n_ie); % rad/sec

% Stationary Error Analysis
fE = 0;
fN = 0;
fD = -g;
vD = 0;
vN = 0;
vE = 0;

% Table 6.1 [2]
OMEGAN = omega_ie * cos(phi);
OMEGAD = -omega_ie * sin(phi);
rhoN = vE/Re;
rhoE = -vN/Re;
rhoD = -vE*tan(phi)/Re;
omegaN = OMEGAN + rhoN;
omegaE = rhoE;
omegaD = OMEGAD + rhoD;

% Table 6.2 [2]
kD = vD/Re;
F41 = -2*OMEGAN*vE - rhoN*vE/(cos(phi)^2);
F43 = rhoE*kD - rhoN*rhoD;
F51 = 2*(OMEGAN*vN + OMEGAD*vD) + rhoN*vN/(cos(phi)^2);
F53 = -rhoE*rhoD - kD*rhoN;
F55 = kD - rhoE*tan(phi);
F63 = rhoN^2+rhoE^2 - (2*g/Re);

% F terms
F13 = rhoE/Re;
F14 = 1/Re;
F21 = -rhoD/cos(phi);
F23 = -rhoN/(Re*cos(phi));
F25 = 1/(Re*cos(phi));
F36 = -1;
F44 = kD;
F45 = 2*omegaD;
F46 = -rhoE;
F48 = fD;
F49 = -fE;
F54 = -(omegaD+OMEGAD);
F56 = omegaN + OMEGAN;
F57 = -fD;
F59 = fN;
F61 = -2*vE*OMEGAD;
F64 = 2*rhoE;
F65 = -2*omegaN;
F67 = fE;
F68 = -fN;
F71 = -OMEGAD;
F73 = rhoN/Re;
F75 = -1/Re;
F78 = omegaD;
F79 = -omegaE;
F83 = rhoE/Re;
F84 = 1/Re;
F87 = -omegaD;
F89 = omegaN;
F91 = OMEGAN + (rhoN/(cos(phi)^2));
F93 = rhoD/Re;
F95 = tan(phi)/Re;
F97 = omegaE;
F98 = -omegaN;

% dxdot = F(t)*dx(t) + f(t);
F = [...
    0   0   F13 F14 0   0   0   0   0   ;
    F21 0   F23 0   F25 0   0   0   0   ;
    0   0   0   0   0   F36 0   0   0   ;
    F41 0   F43 F44 F45 F46 0   F48 F49 ;
    F51 0   F53 F54 F55 F56 F57 0   F59 ;
    F61 0   F63 F64 F65 0   F67 F68 0   ;
    F71 0   F73 0   F75 0   0   F78 F79 ;
    0   0   F83 F84 0   0   F87 0   F89 ;
    F91 0   F93 0   F95 0   F97 F98 0   ;...
    ];

% Remove 3rd and 6th row and column because vertical dynamics are unstable
% They will cause other states to blow up and oscillations will not be
% observed. First remove 6th and then 3rd so no ambiguity in matrix
% indexes.
F(6,:) = [];
F(:,6) = [];
F(3,:) = [];
F(:,3) = [];

% Extract gyro and Accelerometer error States
dae = dx(10:12);
dge = dx(13:15);

% 7 states
derrors = [0;0;dae(1:2);dge]; 
dx = [dx(1:2); dx(4:5); dx(7:9)];

% Return forced dxdot, i.e. with gyro and Accelerometer bias
% Accelerometer Errors affect accelration equation i.e. velocity components
% Rate Gyro Errors affect small-angle rotations
dxdot = F*dx + derrors;
dxdot = [dxdot(1:2); 0; dxdot(3:4); 0; dxdot(5:7); zeros(6,1)];
end