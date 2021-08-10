clear; close all; clc

%% setup
Time = 42; %hours
tspan = [0,Time];
tvec = linspace(0,Time,1000); %vector of time points for which to define BMP and Activin/NODAL
ic = [0;0;0;0]; %initial conditions

R = 350; %colony radius (um)
nbins = 34;

imsize = [round(2.2*R) round(2.2*R)];
DIST = zeros(imsize); DIST(round(imsize(1)/2),round(imsize(2)/2)) = 1;
DIST = bwdist(DIST);

edges = linspace(0,R,nbins+1); %bin edges
d = 0.5*(edges(1:end-1) + edges(2:end)); 
t = linspace(0,Time,1000); %times
[D,T] = meshgrid(d,t);

%activin rescue doses
ndoses = 6;
doses = linspace(0,1,ndoses);

%BMP conditions:
BMP = cell(2,1);
%control
BMP{1} = dynamicBMP(D,T);
%LDN at 24 hours
BMP{2} = dynamicBMP(D,T);
BMP{2}(T>24) = 0;

%NODAL/Activin conditions:
NODAL = cell(9,1);
%WT BMP4-treated (control)
NODAL{1} = dynamicNodal(D,T);
%Dose response (sustained treatment)
for idx = 1:ndoses
    NODAL{idx+1} = rescueActivin(D,T,doses(idx));
end
%max dose only in or after first 24 hours
NODAL{ndoses+2} = rescueActivin(D,T,1);
NODAL{ndoses+2}(T>24) = 0;
NODAL{ndoses+3} = rescueActivin(D,T,1);
NODAL{ndoses+3}(T<24) = 0;

treats = [[ones(9,1) (1:9)'];[2,1]]; %all relevant treatment conditions
nconds = size(treats,1);
pies = {[1,10],2:7,[2,7,8,9]};

%threshold level to define AP2C/SOX17 as on or off
thresh = 0.5;

%% define parameters
%choose which model to use
%model 1 = simplest, model 4 = final, most complex
model = 4;

n = 2; %Hill function exponent
ns = 1; %Hill exponent for BMP and Act
t_half = 4; %set the time scale
alpha = log(2)/t_half; %use the time scale to determine alpha
beta = alpha; %set beta = alpha so steady-state levels are 1

%stuff for Z and W
ss_lowW = 0.1; %signaling-mediated steady-state
ss_lowZ = 0.2;
lambdaZ = 0.22;
lambdaW = 1.3;

thW = lambdaW*ss_lowW; %threshold
thZ = lambdaZ*ss_lowZ; %threshold

opts = struct(...
    'n',n,'ns',ns,... %determines shape of sigmoid activation/inhibition curves
    'betaX',beta,'betaY',beta,...
    'betaZZ',alpha*(thZ+1),'betaAZ',alpha*ss_lowZ,...
    'betaWW',alpha*(thW+1),'betaBW',alpha*ss_lowW,... %production rates
    'alphaX',alpha,'alphaY',alpha,'alphaZ',alpha,'alphaW',alpha,... %degradation rates
    'Kbx',0.5,'Kaby',0.25,'Kaz',0.5,'Kbw',0.9,... %BMP/Act activation thresholds
    'Kzx',0.3,'Kwy',0.08,'Kxz',0.43,'Kyw',0.45,... %inhibition thresholds
    'Kxy',0.5,'Kzz',sqrt(thZ),'Kww',sqrt(thW),...%autoactivation thresholds
    'Kd',100); %disassociation constant for XZ

fn = fieldnames(opts);

%graphical options
fs = 18; %font size
lw = 2; %line width

%folder to save figures
savedir = 'G:\My Drive\Research\Heemskerk lab\HeemskerkLabFiles\Idse\PGC\modelSimulation\figures';
%protein names
channels = {'AP2C','SOX17','ME','AM'};


%test different models
if model == 1
    %minimal model, just need to specify alphaX, alphaY, betaX, betaY, Kbx,
    %Kay, Kxy, and Kd
    opts.Kay = opts.Kaby; %maybe need to tune this
elseif model == 2
    %get rid of production of Z and W
    opts.betaBW = 0;
    opts.betaAZ = 0;
elseif model == 3
    %get rid of production of W
    opts.betaBW = 0;
else
    %no change if model == 4
end

%% compare conditions for given parameters
close all
order = [1 2 3];

%iterate over treatment conditions
V = cell(nconds,1);
parfor idx = 1:nconds
    ii = treats(idx,1); jj = treats(idx,2);
    vals = zeros(nbins,4);
    %iterate over bins within the treatment condition
    for bi = 1:nbins
        bmp = BMP{ii}(:,bi)';
        act = NODAL{jj}(:,bi)';
        if model == 1
            [~,y] = ode45(@(t,y) ODE_model1(t,y,tvec,bmp,act,opts), tspan, ic);
        else
            [~,y] = ode45(@(t,y) ODE_4component(t,y,tvec,bmp,act,opts), tspan, ic);
        end
        vals(bi,1) = y(end,1); vals(bi,2) = y(end,2);
        vals(bi,3) = y(end,3); vals(bi,4) = y(end,4);
    end
    V{idx} = vals;
end

P = cell(3,1);
imgs = vals2ims(V,DIST,edges);
for ii = 1:3
    subplot_tight(1,3,ii)
    P{ii} = pieSlices(imgs(pies{ii},:),R,order,channels);
end
set(gcf,'Position',[167 305 1724 965])
% pause
% close all

%% save pie plots
pref = ['model',num2str(model),'_',channels{order(1)},'_',...
    channels{order(2)},'_',channels{order(3)},'_'];
savenames = {'WT_CNTRLvLDN','NODALKO_doseResponse','NODALKO_ActDurations'};

for ii = 1:3
    savename = fullfile(savedir,[pref,savenames{ii},'.png']);
    imwrite(P{ii},savename)
end

%% solve the ode numerically and plot wrt time
%helps to get a handle on the actual dynamics of each protein in different
%conditions and at different radial bins
close all
f = figure('WindowState','maximized');
for dist = 0:10:120%350 %distance from colony edge
    dose = 0.2;
    bmp = dynamicBMP(dist,tvec);
    act = dynamicNodal(dist,tvec);
%     act = rescueActivin(dist,tvec,dose);
%     act(tvec > 24) = 0;
    
    if model == 1
        [t,y] = ode45(@(t,y) ODE_model1(t,y,tvec,bmp,act,opts), tspan, ic);
    else
        [t,y] = ode45(@(t,y) ODE_4component(t,y,tvec,bmp,act,opts), tspan, ic);
    end
    
    subplot(2,1,1)
    cla
    hold on
    plot(tvec,bmp,'-','LineWidth',lw)
    plot(tvec,act,'--','LineWidth',lw)
    hold off
    ylabel('Signaling levels')
    legend('BMP','Act','Location','southeast')
    cleanSubplot(fs)
    ylim([0,1])
    xlim([0,Time])
    title(strcat("dist = ",num2str(dist)))
    
    subplot(2,1,2)
    cla
    ls = {'-','--','-.',':'};
    hold on
    for ii = 1:4
        plot(t,y(:,ii),ls{ii},'LineWidth',lw)
    end
    yline(opts.Kwy);
    hold off
    ylabel('Protein concentration (au)')
    xlabel('Time (au)')
    legend('AP2C','SOX17','ME','AM','Location','southeast')
    cleanSubplot(fs)
    ylim([0,1])
    xlim([0,Time])
    pause
end
close(f)

%% Make a phase diagram of fate given specific BMP, Act conditions

%for evenly spaced BMP and Act between 0 and 1, simulate time course of
%differentiation
Time = 200;
tspan = [0,Time];
tvec = linspace(0,Time,1000);
bmpl = linspace(0,1,100);
actl = linspace(0,1,100);

tic
X = zeros(length(bmpl)*length(actl),1); Y = X;
% A = X; B = X; C = X;
M = zeros(length(bmpl)*length(actl),4);
xx = zeros(size(M,1),1); yy = xx;
idx = 1;
for ii = 1:length(bmpl)
    fprintf('.')
    if mod(ii,45) == 0
        fprintf('\n')
    end
    for jj = 1:length(actl)
        bmp = bmpl(ii)*ones(size(tvec));
        act = actl(jj)*ones(size(tvec));
        if model == 1
            [~,y] = ode45(@(t,y) ODE_model1(t,y,tvec,bmp,act,opts), tspan, ic);
        else
            [~,y] = ode45(@(t,y) ODE_4component(t,y,tvec,bmp,act,opts), tspan, ic);
        end
        X(idx) = bmpl(ii); Y(idx) = actl(jj);
        M(idx,:) = y(end,:);
        xx(idx) = ii; yy(idx) = jj;
        idx = idx + 1;
    end
end
fprintf('\n')
toc

%% plot phase diagram image
im = zeros(length(bmpl),length(actl),3);

order = [1 2 3];
for idx = 1:length(X)
    ii = xx(idx); jj = yy(idx);
    im(ii,jj,1) = M(idx,order(1));
    im(ii,jj,2) = M(idx,order(2));
    im(ii,jj,3) = M(idx,order(3));
end

figure
imagesc([min(actl),max(actl)],[min(bmpl),max(bmpl)],im)
axis square
set(gca,'YDir','normal') 
set(gca,'TickLength',[0,0])
set(gca,'Box','off')
cleanSubplot
xlabel('NODAL/Activin')
ylabel('BMP')
title(strcat("\color{red}", channels{order(1)},...
    " ", "\color{green}", channels{order(2)},...
    " ", "\color{blue}", channels{order(3)}))

%% save figure
title('')
savename = ['model',num2str(model),'_','phaseDiagram_',channels{order(1)},'_',channels{order(2)},...
    '_',channels{order(3)}];
saveas(gcf,fullfile(savedir,[savename,'.png']))
saveas(gcf,fullfile(savedir,'matlab_fig_files',[savename,'.fig']))

%% Kymographs of BMP and NODAL dynamics in WT BMP4-treated micropatterns

d = linspace(0,R,1000); %edge distances
t = linspace(0,Time,100); %times
[D,T] = meshgrid(d,t);
bmp = dynamicBMP(D,T); %compute BMP wrt to time and distance
nodal = dynamicNodal(D,T); %compute NODAL wrt to time and distance
act = rescueActivin(D,T,1);

%quantities to plot wrt to time and edge distance
plots = {bmp,nodal,min((nodal + bmp),1),act};
names = {'BMP signaling','WT NODAL signaling','Combined signaling',...
    'NODALKO + Activin'};
snames = {'BMP_WT','NODAL_WT','Combined_WT','Activin_NODALKO'};

close all
figure
for ii = 1:length(plots)
    clf
    surf(T,D,plots{ii},'LineStyle','none')
    xlabel('time (hr)')
    ylabel('edge distance (um)')
    view(2)
    cleanSubplot
    ylim([min(d),max(d)])
    xlim([min(t),max(t)])
    colormap jet
    title(names{ii})
    axis square
    
    savename = ['kymograph',snames{ii}];
    saveas(gcf,fullfile(savedir,[savename,'.png']))
    saveas(gcf,fullfile(savedir,'matlab_fig_files',[savename,'.fig']))
    pause(0.5)
end


%% local functions
function BMP = dynamicBMP(d,t,sigma,tau)

if ~exist('sigma','var')
    sigma = 30;
end

if ~exist('tau','var')
    tau = 12;
end

BMP = heaviside(tau - t) + heaviside(t - tau).*exp(-(d/sigma).^2);

end

function Nodal = dynamicNodal(d,t,v,tau)

if ~exist('v','var')
    v = 10.5;
end

if ~exist('tau','var')
    tau = 24;
end

Nodal = heaviside(v*(t-tau) - d);

end

function Act = rescueActivin(d,t,dose,sigma)

if ~exist('lambda','var')
    sigma = 90;
end

Act = ones(size(t)).*dose.*exp(-(d/sigma).^2);

end

function imgs = vals2ims(V,DIST,edges)

nconds = length(V);
imgs = cell(nconds,4);
imsize = size(DIST);
B = zeros(imsize);
bmask = DIST >= imsize(1)/2;

for ii = 1:nconds
    for ci = 1:4
        imgs{ii,ci} = B;
        for bi = 1:length(edges)-1
            mask = (DIST >= edges(bi)) & (DIST < edges(bi+1));
            imgs{ii,ci}(mask) = V{ii}(end+1-bi,ci);
        end
        imgs{ii,ci}(bmask) = 1;
    end
end

end

function pie = pieSlices(imgs,radius,order,channels)

Rmax = uint16(radius + 25); %25 micron margin
[X,Y] = meshgrid(1:size(imgs{1},2),1:size(imgs{1},1));
center = round([size(imgs{1},1)/2,size(imgs{1},2)/2]);
R = sqrt((X - center(1)).^2 + (Y - center(2)).^2);
F = atan2((X - center(1)),-(Y - center(2)));
disk = R > Rmax;

N = size(imgs,1);
mask = cell(N,1); 
for ii = 1:N
    mask{ii} = F < -pi + (ii-1)*2*pi/N | F > -pi + ii*2*pi/N;
    mask{ii} = imdilate(mask{ii},strel('disk',10));
end
lines = mask{1};
for ii = 2:N
    lines = lines & mask{ii};
end

pie = cell(1,1,3);
for ci = 1:3
    pie{ci} = zeros(size(disk));
end

for ii = 1:N
    for cii = 1:3
        pie{cii}(~mask{ii}) = imgs{ii,order(cii)}(~mask{ii});
        pie{cii}(disk) = 1;
        pie{cii}(lines) = 1;
    end
end

pie = cell2mat(pie);

imshow(pie,'InitialMagnification','fit')
cleanSubplot
% set(gcf,'Position',[950 650 660 684])
title(strcat("\color{red}", channels{order(1)},...
    " ", "\color{green}", channels{order(2)},...
    " ", "\color{blue}", channels{order(3)}))

end

function cleanSubplot(varargin)

if nargin == 0
    fs = 20;
    lw = 2;
elseif nargin == 1
    fs = varargin{1};
    lw = 2;
elseif nargin == 2
    fs = varargin{1};
    lw = varargin{2};
end
fgc = 'k';
bgc = 'w';
graphbgc = 1*[1 1 1];

set(gcf,'color',bgc);
set(gca, 'LineWidth', lw);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
set(gca,'XColor',fgc);
set(gca,'YColor',fgc);
set(gca,'Color',graphbgc);

end
