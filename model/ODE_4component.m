function dydt = ODE_4component(t,y,tvec,BMP,Act,opts)
%t = current time
%y = level of each transciption factor at time t
%components of y:
%y(1) = X = AP2C; y(2) = Y = SOX17; y(3) = Z = ME TF; y(4) = W = AM TF
%tvec = vector of time points for which BMP and Act levels are defined
%BMP, Act = BMP and Act levels at each time in tvec
%opts = structure array of parameters for the ODE model

%interpolate BMP and Activin levels at the given time point
if numel(tvec) > 1
    bmp = interp1(tvec,BMP,t,'nearest'); % Interpolate at time t
    act = interp1(tvec,Act,t,'nearest'); % Interpolate at time t
elseif numel(tvec) == 1 %if only one level is provided, no need to interpolate
    bmp = BMP;
    act = Act;
end

%protein concentrations
X = y(1); Y = y(2); Z = y(3); W = y(4);
%%%parameters%%%
%Hill coefficient (assume the same for all activation/inhibition functions)
n = opts.n; ns = opts.ns;
%production rates
betaX = opts.betaX; betaY = opts.betaY;
betaAZ = opts.betaAZ; betaZZ = opts.betaZZ; betaBW = opts.betaBW; betaWW = opts.betaWW;
%degradation rates
alphaX = opts.alphaX; alphaY = opts.alphaY; alphaZ = opts.alphaZ; alphaW = opts.alphaW;
%signaling activation thresholds
Kbx = opts.Kbx; Kaby = opts.Kaby; Kaz = opts.Kaz; Kbw = opts.Kbw;
%autoactivation thresholds
Kxy = opts.Kxy; Kzz = opts.Kzz; Kww = opts.Kww;
%inhibition thresholds
Kzx = opts.Kzx; Kwy = opts.Kwy; Kxz = opts.Kxz; Kyw = opts.Kyw;
%association/dissociation constant for [XY]
Kd = opts.Kd;

%separation of timescales -> set [XY] to equilibrium value (for current 
%total X and Y) each time dydt is evaluated
p = [1, -(X + Y + 1/Kd), X*Y];
XY = min(roots(p));
% Xfree = X - XY;
% Yfree = Y - XY;

dydt = NaN(4,1);
%dX/dt (AP2C)
dydt(1) = betaX*((bmp/Kbx)^ns + (XY/Kxy)^n)/...
    (1 + (bmp/Kbx)^ns + (XY/Kxy)^n)*1/(1 + (Z/Kzx)^n) - alphaX*X;
%dY/dt (SOX17)
dydt(2) = betaY*((bmp*act/Kaby)^ns + (XY/Kxy)^n)/...
    (1 + (bmp*act/Kaby)^ns + (XY/Kxy)^n)*1/(1 + (W/Kwy)^n) - alphaY*Y;
%dZ/dt (EOMES)
% dydt(3) = (betaAZ*(act/Kaz)^n + betaZZ*(Z/Kzz)^n)/...
%     (1 + (act/Kaz)^n + (Z/Kzz)^n)*1/(1 + (X/Kxz)^n) - alphaZ*Z;
dydt(3) = (betaAZ*(act/Kaz)^ns/(1 + (act/Kaz)^ns) +...
    betaZZ*(Z/Kzz)^n/(1 + (Z/Kzz)^n))/(1 + (X/Kxz)^n) - alphaZ*Z;
%dW/dt (ISL1)
% dydt(4) = (betaBW*(bmp/Kbw)^n + betaWW*(W/Kww)^n)/...
%     (1 + (bmp/Kbw)^n + (W/Kww)^n)*1/(1 + (Y/Kyw)^n) - alphaW*W;
dydt(4) = (betaBW*(bmp/Kbw)^ns/(1 + (bmp/Kbw)^ns) +...
    betaWW*(W/Kww)^n/(1 + (W/Kww)^n))/(1 + (Y/Kyw)^n) - alphaW*W;




end