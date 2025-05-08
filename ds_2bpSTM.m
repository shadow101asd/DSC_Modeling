function [x1_C, STM_M] = ds_2bpSTM(x0_C,dt,GM,x1_C)
%
% [x1_C, STM_M] = ds_2bpSTM(x0_C,dt,GM, x1_C)
%
% Compute the final state and/or state transition matrix using universal
% variables. Sources: Battin and Vallado (see details).
% NOTE This is an old routine, which I currnely use for the Caltech Class
% Possible calls:
%
% x1_C          = ds_2bpSTM(x0_C,dt,GM)
% [x1_C, STM_M] = ds_2bpSTM(x0_C,dt,GM)
% [~   , STM_M] = ds_2bpSTM(x0_C,dt,GM, x1_C)
%
% INPUT:
% name      Dim     Description
% x0_C      6x1     Initial state
% dt        1x1     Time interval
% *** optional ***
% x1_C      6x1     Final state (you can pass it as input to compute STM only)
%
% OUTPUT:
% name      Dim     Description
% x1_C      6x1     Final state
% STM_M     6x6     State Transition Matrix
%
% --------------------------last updated SC, 2012 -------------------
% if dt ==0
%     x1_C    = x0_C;
%     STM_M   = eye(6);
%     return
% end
% NORMALIZATION
LSF = norm(x0_C(1:3));
VSF = norm(x0_C(4:6));
XSF = [LSF LSF LSF VSF VSF VSF].';
TSF = LSF/VSF;
GMSF = LSF^3/TSF^2;
x0_C = x0_C./XSF;
dt = dt/TSF;
GM = GM/GMSF;

%% INPUT
FLAG_compute_X1_C = 0; FLAG_compute_STM = 0;
if nargin   == 3
    FLAG_compute_X1_C   = 1;% x1_C must be computed
else
    x1_C = x1_C./XSF;
end
if nargout  == 2
    FLAG_compute_STM    = 1;% STM must be computed
end

r0_C    = x0_C(1:3);
v0_C    = x0_C(4:6);
r0      = work_norm(r0_C.');
v0      = work_norm(v0_C.');

%% COMPUTE SIGMA0, ALPHA0
% Battin, pp 174/175
si0     = work_dot(r0_C.',v0_C.')/sqrt(GM);
alpha   = 2/r0 - v0^2/GM;

%% COMPUTE XSI, R1
toll    = 1e-12;
if FLAG_compute_X1_C
    [r1,xsi] = aux_xsi(alpha, GM, si0, r0, v0, dt, toll);
else
    % Battin, pp 179
    r1_C    = x1_C(1:3);
    v1_C    = x1_C(4:6);
    r1      = work_norm(r1_C.');
    si1     = work_dot(r1_C.',v1_C.')/sqrt(GM);
    xsi     = alpha*sqrt(GM)*dt + si1 - si0;
    
end
%% Compute U0, U1, U2,U3, U4, U5
[~,U1,U2,~,U4,U5] = aux_univ(alpha, xsi);

%% COMPUTE F, G, Ft, Gt, X1_C(if not provided), and return (if no STM is needed)
% Otherwise use
% r = r0*U0+si0*U1+U2;
% Battin, pp 179
F       = 1 - U2/r0;
G       = r0*U1/sqrt(GM)+U2*si0/sqrt(GM);
Ft      = -sqrt(GM)*U1/r1/r0;
Gt      = 1-U2/r1;

if FLAG_compute_X1_C
    r1_C = F*r0_C    + G*v0_C;
    v1_C = Ft*r0_C   + Gt*v0_C;
    
    if ~FLAG_compute_STM
        x1_C     = [r1_C; v1_C];
        x1_C    = x1_C.*XSF;
        return
    end
end
%% Compute STM
% Battin pp 466, Eq 9.74
C       = (3*U5-xsi*U4)/sqrt(GM)-dt*U2;
% Battin, pp.467
dv_C    = v1_C-v0_C;
dr_C    = r1_C-r0_C;


drdr0_M = r1/GM*(dv_C*dv_C.')...
    + (r0*(1-F)*r1_C*r0_C.' + C*v1_C*r0_C.') /r0^3 ...
    + F*eye(3);% Rtilde

drdv0_M = r0/GM*(1-F)*(dr_C*v0_C.' - dv_C*r0_C.')...
    +C/GM*v1_C*v0_C.'   +G*eye(3); % R

dvdr0_M = - dv_C*r0_C.'/r0^2 - r1_C*dv_C.'/r1^2 ...
    - Ft*(r1_C*r1_C.')/r1^2 + Ft*(r1_C*v1_C.' - v1_C*r1_C.')*(r1_C*dv_C.')/GM/r1...%;%...
    - GM*C/r1^3/r0^3*r1_C*r0_C.';% Vtilde
dvdr0_M = dvdr0_M+Ft*eye(3);

dvdv0_M = (dv_C*dv_C.')*r0/GM ...
    +  (r0*(1-F)*r1_C*r0_C.'  - C* r1_C*v0_C.') /r1^3 ...
    + Gt*eye(3);% Rtilde

%% SCALE BACK THE OUTPUT
x1_C  = [r1_C; v1_C];
STM_M = [ drdr0_M, drdv0_M
    dvdr0_M,  dvdv0_M];

% DIMENSIONALIZE
x1_C    = x1_C.*XSF;
STM_M   = [drdr0_M drdv0_M*LSF/VSF;
    dvdr0_M*VSF/LSF dvdv0_M];
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UNIVERSAL FUNCTIONS
function [U0,U1,U2,U3,U4,U5] = aux_univ(alpha, xsi)
if abs(alpha)<1e-10
    % PARABOLA
    %  Use Battin, pp176, Eq. 4.76 with alpha = 0;
    U0 = 1;
    U1 = xsi;
    U2 = xsi^2/2;
    U3 = xsi^3/6;
    U4 = xsi^4/24;
    U5 = xsi^5/120;
else
    % ELLIPSE/HYPERBOLA
    % Use Battin, pp180, and pp176, Eq. 4.76 with alpha = 0;
    if alpha>=1e-10 % Ellipse
        U0 = cos(sqrt(alpha)*xsi);
        U1 = sin(sqrt(alpha)*xsi)/sqrt(alpha);
    else % Hyperbola
        U0 = cosh(sqrt(-alpha)*xsi);
        U1 = sinh(sqrt(-alpha)*xsi)/sqrt(-alpha);
    end
    
    U2 = (1       - U0)/alpha;
    U3 = (xsi     - U1)/alpha;
    U4 = (xsi^2/2 - U2)/alpha;
    U5 = (xsi^3/6 - U3)/alpha;
end
end
%% KEPLER'S EQUATION
% My notes on Battin p 179; See also Vallado 70 and 101
function [r1,xsi] = aux_xsi(alpha, GM, si0, r0, v0, dt, toll)
%  GENERATE THE FIRST GUESS
if abs(alpha)<1e-7  % PARABOLA
    % Eq. for p from a) + b) + c)
    % a) si0^2 = r0^2*v0^2*cos(theta)/GM
    % b) h^2   = r0^2*v0^2*sin^2(theta) =  r0^2*v0^2*[1 - cos^2(theta)]
    % c) p     = h^2/GM
    p = r0^2*v0^2/GM-si0^2;
    s = atan( sqrt(p^3/GM) /3 / dt);
    w = atan( tan(s)^(1/3));
    xsi = 2*sqrt(p)*cot(2*w);
elseif alpha>=1e-7  % ELLIPSE
    xsi = sqrt(GM)*dt*alpha;
else                % HYPERBOLA
    a = 1/alpha;
    xsi = sign(dt)*sqrt(-a)*log( -2*alpha*dt*sqrt(GM)/...
        (si0+sign(dt)*sqrt(-a)*(1-r0*alpha)));
end

xsiold = 1e6;
% FIND THE ZERO OF KEPLER'S EQUATION USING NEWTON'S METHOD
while abs(xsi-xsiold)>toll
    
    [U0,U1,U2,U3] = aux_univ(alpha, xsi);
    
    r1      = r0*U0+si0*U1+U2;% = -dKepler/dxsi
    Kepler  = sqrt(GM)*dt - (r0*U1+si0*U2+U3);
    % Newton's method iteration
    xsiold  = xsi;
    xsi     = xsiold + Kepler/r1;
    
end


end

function y_C = work_norm(x_M)
% y_C = work_norm(x_M)
% Compute the N norms of the N row vectors 
%         x_M(i,:), i=1,N
if isempty(x_M),y_C = [];return;end
y_C = (work_dot(x_M,x_M)).^(0.5);
% faster than sqrt!
end

function z_C = work_dot(x_M,y_M)
%  z_C =work_dot(x_M,y_M)
%
% Compute the N dot products of the N row vectors
%         x_M(i,:), i=1,N
% with the N row vectors
%         y_M(i,:), i=1,N
%

[~,ncol]=size(x_M);
z_C = 0;
for i=1:ncol
    z_C  = z_C + x_M(:,i).*y_M(:,i);
end
end







