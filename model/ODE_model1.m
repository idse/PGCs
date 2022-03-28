function dydt = ODE_model1(t,y,tvec,BMP,Act,opts)
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
X = y(1); Y = y(2);
%%%parameters%%%
%Hill coefficient (assume the same for all activation/inhibition functions)
n = opts.n; ns = opts.ns;
%production rates
betaX = opts.betaX; betaY = opts.betaY;
%degradation rates
alphaX = opts.alphaX; alphaY = opts.alphaY;
%activation thresholds
Kbx = opts.Kbx; Kay = opts.Kay; Kxy = opts.Kxy;
%association/dissociation constant for [XY]
Kd = opts.Kd;

%separation of timescales -> set [XY] to equilibrium value (for current 
%total X and Y) each time dydt is evaluated
p = [1, -(X + Y + 1/Kd), X*Y];
XY = min(roots(p));

dydt = zeros(4,1);
%dX/dt (AP2C)
dydt(1) = betaX*((bmp/Kbx)^ns + (XY/Kxy)^n)/...
    (1 + (bmp/Kbx)^ns + (XY/Kxy)^n) - alphaX*X;
%dY/dt (SOX17)
dydt(2) = betaY*((act/Kay)^ns + (XY/Kxy)^n)/...
    (1 + (act/Kay)^ns + (XY/Kxy)^n) - alphaY*Y;



end