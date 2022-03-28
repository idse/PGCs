clear; close all; clc

%% setup
%load ground truth data
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
dataDir = scriptPath;
load(fullfile(dataDir,'ground_truth.mat'))
groundTruth = cell(size(targetProfiles));

%define space and time coordinate for simulation
Time = 42; %hours
tspan = [0,Time];
tvec = linspace(0,Time,1000); %vector of time points for which to define BMP and Activin/NODAL
ic = [0;0;0;0]; %initial conditions

R = 350; %colony radius (um)
nbins = 34;

imsize = [round(2.2*R) round(2.2*R)];
DIST = zeros(imsize); DIST(round(imsize(1)/2),round(imsize(2)/2)) = 1;
DIST = bwdist(DIST); %use this to make colony simulation visualizations

edges = linspace(0,R,nbins+1); %bin edges
d = 0.5*(edges(1:end-1) + edges(2:end)); 
t = linspace(0,Time,1000); %times
[D,T] = meshgrid(d,t);

%interpolate target gene expression profiles to the same radial bins used
%in simulation
for ii = 1:size(groundTruth,1)
    for jj = 1:size(targetProfiles{ii},2)
        groundTruth{ii}(:,jj) = interp1(rad,targetProfiles{ii}(:,jj),d)';
    end
end

%activin rescue doses
ndoses = 6;
doses = linspace(0,1,ndoses);
ndoses = length(doses);

%BMP conditions:
BMP = cell(2,1);
%control
BMP{1} = dynamicBMP(D,T); %WT BMP spatiotemporal dynamics
%LDN at 24 hours
BMP{2} = dynamicBMP(D,T);
BMP{2}(T>20) = 0; %WT BMP dynamics with LDN at 20 hours

BMPlabels = {'BMP50','LDN@20h'}; %labels for BMP treatment conditions

%NODAL/Activin conditions:
NODAL = cell(ndoses+3,1);
%WT BMP4-treated (control)
NODAL{1} = dynamicNodal(D,T);
%Dose response (sustained treatment)
for idx = 1:ndoses
    NODAL{idx+1} = rescueActivin(D,T,doses(idx));
end
%max dose only in or after first 24 hours (simulate addition of SB to
%completely cut off signaling)
NODAL{ndoses+2} = rescueActivin(D,T,1);
NODAL{ndoses+2}(T<24) = 0;
NODAL{ndoses+3} = rescueActivin(D,T,1);
NODAL{ndoses+3}(T>24) = 0;
%labels for Activin/Nodal treatment conditions
NODALlabels = cellfun(@(x) strrep(x,' ',''),strcat('A',cellstr(num2str(doses')))','UniformOutput',false);
NODALlabels = [{'WT Nodal'},NODALlabels,{'ACT after 24h','ACT first 24h'}];

%all relevant treatment conditions -> should correspod to labels for 
%target data; each row is [BMP index, NODAL index]
treats = [[ones(ndoses,1),(2:ndoses+1)'];[1,1];[2,1];[ones(2,1),((ndoses+2):(ndoses+3))']];
test = strcat(BMPlabels(treats(:,1))',repmat({'_'},size(treats,1),1),NODALlabels(treats(:,2))');
%make sure target and simulation conditions match
disp([test,labels])
nconds = size(treats,1);

%% define starting parameters
n = 2; %Hill function coefficient
ns = 2; %Hill coefficient for BMP and Act
t_half = 4; %set the time scale
alpha = log(2)/t_half; %use the time scale to determine alpha
beta = alpha; %set beta = alpha so steady-state levels are 1

%stuff for Z and W
ss_lowW = 0.1; ss_lowZ = 0.2; lambdaZ = 0.2; lambdaW = 1.2;
thW = lambdaW*ss_lowW; thZ = lambdaZ*ss_lowZ; %threshold

opts = struct(...
    'n',n,'ns',ns,... %determines shape of sigmoid activation/inhibition curves
    'betaX',beta,'betaY',beta,...
    'betaZZ',alpha*(thZ+1),'betaAZ',alpha*ss_lowZ,...
    'betaWW',alpha*(thW+1),'betaBW',alpha*ss_lowW,... %production rates
    'alphaX',alpha,'alphaY',alpha,'alphaZ',alpha,'alphaW',alpha,... %degradation rates
    'Kbx',0.3,'Kaby',0.18,'Kaz',0.8,'Kbw',0.95,... %BMP/Act activation thresholds
    'Kzx',0.1,'Kwy',0.08,'Kxz',0.4,'Kyw',0.85,... %inhibition thresholds
    'Kxy',0.5,'Kzz',sqrt(thZ),'Kww',sqrt(thW),...%autoactivation thresholds
    'Kd',100); %disassociation constant for XZ

fn = fieldnames(opts);
channels = {'TFAP2C','SOX17','ME','AM'};

% model parameters over which to optimize:
fields = {...
    'thZ', 'thW',...
    'ss_lowZ', 'ss_lowW',...
    'Kbx', 'Kaby', 'Kaz', 'Kbw',...
    'Kzx', 'Kwy', 'Kxz', 'Kyw'};
%use local functions to update theta from opts and vice-versa
theta = opts2theta(opts,fields);
nt = length(theta);

%% run the optimization
%create a folder to save results, labeled with the current date
c = clock;
savedir = sprintf('%d%.2d%.2d',c(1),c(2),c(3));
savedir = fullfile(dataDir,savedir(3:end));
if ~exist(savedir,'dir'), mkdir(savedir); end

%set hyperparameters
bmax = 20; %for efficiency, don't need to simulate every radial bin, behavior of interest is at the edge
mu = zeros(nt,1); %zero vector for mean of gaussian rv
sigma = 5e-4; %variance determines the step size for each parameter 
Sigma = diag(sigma*[0.2,0.2,0.2,0.2,1,1,1,1,1,1,1,1]); %diagonal covariance matrix
tempmax = 0.1; %'temperature' determines the acceptance criterion
niters = 1500; %number of iterations per optimization run
nrepeats = 40; %number of independent repeats of the optimzation procedure

for ri = 1:nrepeats
    %randomly initialize parameters
    theta = rand(nt,1);
    newopts = theta2opts(theta,opts,fields);
    V = cell(nconds,1);
    %simulate once to determine initial error
    parfor idx = 1:nconds
        ii = treats(idx,1); jj = treats(idx,2);
        vals = zeros(nbins,2);
        for bi = 1:bmax
            bmp = BMP{ii}(:,bi)';
            act = NODAL{jj}(:,bi)';
            [~,y] = ode45(@(t,y) ODE_4component(t,y,tvec,bmp,act,newopts), tspan, ic);
            vals(bi,:) = y(end,1:2);
        end
        V{idx} = vals;
    end

    errold = mean(cellfun(@(x,y) sum((x(1:bmax,:) - y(1:bmax,:)).^2,'all'),V,groundTruth));
    fprintf('MSE = %g\n',errold)
    
    %initialize variables to store intermediate cost and parameter values
    thetas = zeros(nt,niters);
    costs = zeros(niters,1);
    oldcosts = zeros(niters,1);
    naccept = 0;
    close all
    tic
    %show the error over time
    figure('Position',[1000 602 1250 750])
    for ii = 1:niters
        %update parameters with random step
        thetap = theta + mvnrnd(mu,Sigma,1)';
        %if any parameters stepped below zero, set them below zero and
        %their initial value instead
        halfp = 0.5*theta;
        thetap(thetap <= 0) = halfp(thetap <= 0);
        newopts = theta2opts(thetap,opts,fields);
        %run the simulation over treatment conditions in parallel
        V = cell(nconds,1);
        parfor idx = 1:nconds
            ll = treats(idx,1); jj = treats(idx,2);
            vals = zeros(nbins,2);
            for bi = 1:bmax
                bmp = BMP{ll}(:,bi)';
                act = NODAL{jj}(:,bi)';
                [~,y] = ode45(@(t,y) ODE_4component(t,y,tvec,bmp,act,newopts), tspan, ic);
                vals(bi,:) = real(y(end,1:2));
            end
            V{idx} = vals;
        end
        %error with new parameter values
        errnew = mean(cellfun(@(x,y) sum((x(1:bmax,:) - y(1:bmax,:)).^2,'all'),V,groundTruth));
        
        costs(ii) = errnew;
        oldcosts(ii) = errold;
        thetas(:,ii) = theta';

        %acceptance criterion
        temp = tempmax*(1 - ii/niters);
        if exp((errold - errnew)/temp) > rand(1)
            theta = thetap;
            errold = errnew;
            naccept = naccept + 1;
        end
        %update the graph of costs every 10 iterations
        if mod(ii,10) == 0
            fprintf('.') %progress indicator
            clf
            plot(costs(1:ii),'LineWidth',2)
            hold on
            plot(oldcosts(1:ii),'LineWidth',2)
            hold off
            cleanSubplot
            xlabel('iterations'); ylabel('mean squared error')
            axis square
            drawnow
        end
        if mod(ii,500) == 0
            fprintf('\n')
        end
    end
    fprintf('\nEfficiency = %g\n',naccept/niters)
    fprintf('MSE = %g\n',errold(end))
    toc
    
    %for each run, save a graph of proposed and accepted costs over
    %iterations and the associated parameter and cost values as mat files
    close all
    figure('Position',[1000 602 1250 750])
    plot(costs,'LineWidth',2); hold on
    plot(oldcosts,'LineWidth',2); hold off
    xlabel('iterations'); ylabel('mean squared error')
    cleanSubplot; axis square
    ylim([0,4])
    legend('proposed cost', 'accepted cost')
    saveas(gcf,fullfile(savedir,sprintf('costs_run%.2d.png',ri)))
    
    
    newopts = theta2opts(theta,opts,fields);
    save(fullfile(savedir,sprintf('paramValues_run%.2d',ri)),'newopts','thetas','fields','costs','oldcosts')

    % compare target and calculated profiles
    V = cell(nconds,1);
    for idx = 1:nconds
        ii = treats(idx,1); jj = treats(idx,2);
        vals = zeros(nbins,4);
        for bi = 1:nbins
            bmp = BMP{ii}(:,bi)';
            act = NODAL{jj}(:,bi)';
            [~,y] = ode45(@(t,y) ODE_4component(t,y,tvec,bmp,act,newopts), tspan, ic);
            vals(bi,:) = y(end,:);
        end
        V{idx} = vals;
    end
    errold = mean(cellfun(@(x,y) sum((x(1:bmax,1:2) - y(1:bmax,1:2)).^2,'all'),V,groundTruth));
    fprintf('MSE = %g\n',errold)
    %plot and save simulated and target radial profiles for TFAP2C and
    %SOX17 in each treatment condition
    margin = 0.02;
    close all
    figure('Position',[10 550 2500 700])
    idx = 1;
    for jj = 1:2
        for ii = 1:nconds
            subplot_tight(2,nconds,idx,margin)
            plot(d,V{ii}(:,jj),'LineWidth',2)
            hold on
            plot(d,groundTruth{ii}(:,jj),'LineWidth',2)
            cleanSubplot(12); axis square
            ylim([0,1]); xlim([min(d),max(d)])
            legend('calculated','target')
            idx = idx + 1;
            if ii == 1
                ylabel(channels{jj})
            end
            if jj == 1
                title(labels{ii},'Interpreter','none')
            end
        end
    end
    saveas(gcf,fullfile(savedir,sprintf('profiles_run%.2d.png',ri)))
    
    % plot and save pie plots of calculated profiles
    channels = {'TFAP2C','SOX17','ME','AM'};
    order = [1 2 3];
    pies = {[7,8],1:ndoses,[1,ndoses,10,9]};
    imC = vals2ims(V,DIST,edges);
    
    close all
    figure('WindowState','maximized')
    for ii = 1:length(pies)
        subplot(1,length(pies),ii)
        pieSlices(imC(pies{ii},:),R,order,channels);
        idx = idx + 1;
    end
    saveas(gcf,fullfile(savedir,sprintf('pies_run%.2d',ri)))
end

%% pie plots of target profiles
channels = {'TFAP2C','SOX17','ME','AM'};
order = [1 2 3];
%determine which conditions to put in the same pie plots
pies = {[7,8],1:ndoses,[1,ndoses,10,9]};

VT = cell(size(groundTruth,1),1);
for ii = 1:length(VT)
    VT{ii} = [groundTruth{ii}, zeros(size(groundTruth{ii},1),2)];
end
imT = vals2ims(VT,DIST,edges);

close all
figure('Position',[50 250 2500 750])
P = cell(1,length(pies));
idx = 1;
for ii = 1:length(pies)
    subplot(1,length(pies),idx)
    P{ii} = pieSlices(imT(pies{ii},:),R,order,channels);
    idx = idx + 1;
end

saveas(gcf,fullfile(savedir,'groundTruthPies.png'))

pref = ['groundTruth_',channels{order(1)},'_',...
    channels{order(2)},'_',channels{order(3)},'_'];
savenames = {'WT_CNTRLvLDN','NODALKO_doseResponse','NODALKO_ActDurations'};

for ii = 1:3
    savename = fullfile(savedir,[pref,savenames{ii},'.png']);
    imwrite(P{ii},savename)
end

%% local functions

function opts = theta2opts(theta,oldopts,fields)

opts = oldopts;
for ii = 1:length(fields)
    opts.(fields{ii}) = theta(ii);
end

opts.betaZZ = opts.alphaZ*(opts.thZ + 1);
opts.betaWW = opts.alphaW*(opts.thW + 1);
opts.betaAZ = opts.alphaZ*opts.ss_lowZ;
opts.betaBW = opts.alphaW*opts.ss_lowW;
opts.Kzz = sqrt(opts.thZ); opts.Kww = sqrt(opts.thW);

end

function theta = opts2theta(opts,fields)

opts.thZ = opts.betaZZ/opts.alphaZ - 1;
opts.thW = opts.betaWW/opts.alphaW - 1;
opts.ss_lowZ = opts.betaAZ/opts.alphaZ;
opts.ss_lowW = opts.betaBW/opts.alphaW;

nfields = length(fields);
theta = NaN(nfields,1);
for ii = 1:nfields
    theta(ii) = opts.(fields{ii});
end

end

function BMP = dynamicBMP(d,t,sigma,tau)

if ~exist('sigma','var')
    %spatial extent of BMP signaling gradient in microns
    sigma = 60;
end

if ~exist('tau','var')
    %time that BMP signaling is restricted to the edge in hours
    tau = 12;
end

BMP = heaviside(tau - t) + heaviside(t - tau).*exp(-(d/sigma).^2);

end

function Nodal = dynamicNodal(d,t,v,tau)

if ~exist('v','var')
    %velocity of the traveling nodal wave
    v = 10.5;
end

if ~exist('tau','var')
    %time at which the nodal wave starts traveling in from the colony edge
    tau = 24;
end

Nodal = heaviside(v*(t-tau) - d);

end

function Act = rescueActivin(d,t,dose,sigma)

if ~exist('sigma','var')
    %spatial extent of Activin signaling gradient in microns
    sigma = 100;
end

Act = ones(size(t)).*dose.*exp(-(d/sigma).^2);

end

function imgs = vals2ims(V,DIST,edges)

nconds = length(V);
nc = size(V{1},2);
imgs = cell(nconds,nc);
imsize = size(DIST);
B = zeros(imsize);
bmask = DIST >= imsize(1)/2;

for ii = 1:nconds
    for ci = 1:nc
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
    for cii = 1:length(order)
        pie{cii}(~mask{ii}) = imgs{ii,order(cii)}(~mask{ii});
    end
end

for ii = 1:N
    for cii = 1:3
        pie{cii}(disk) = 1;
        pie{cii}(lines) = 1;
    end
end

pie = cell2mat(pie);

imshow(pie,'InitialMagnification','fit')
cleanSubplot

if length(order) == 2
    title(strcat("\color{red}", channels{order(1)},...
        " ", "\color{green}", channels{order(2)}))
elseif length(order) == 3
    title(strcat("\color{red}", channels{order(1)},...
        " ", "\color{green}", channels{order(2)},...
        " ", "\color{blue}", channels{order(3)}))
end

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
