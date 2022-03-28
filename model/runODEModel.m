clear; close all; clc

%% setup
%load ground truth data
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
dataDir = scriptPath;

%set a folder to which to save results
c = clock;
savedir = sprintf('%d%.2d%.2d',c(1),c(2),c(3));
savedir = fullfile(dataDir,savedir(3:end));
if ~exist(savedir,'dir'), mkdir(savedir); end

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

%activin rescue doses -> doses are evenly spaced on a log scale
ndoses = 6;
doses = linspace(0,1,ndoses);
ndoses = length(doses);

%BMP conditions:
BMP = cell(2,1);
%control (WT BMP dynamics with 50 ng/mL treatment)
BMP{1} = dynamicBMP(D,T);
%LDN at 24 hours
BMP{2} = dynamicBMP(D,T);
BMP{2}(T>20) = 0; %LDN at 20 hours

BMPlabels = {'BMP50','LDN@20h'};

%NODAL/Activin conditions:
NODAL = cell(ndoses+3,1);
%WT BMP4-treated (control)
NODAL{1} = dynamicNodal(D,T);
%Dose response (sustained treatment)
for idx = 1:ndoses
    NODAL{idx+1} = rescueActivin(D,T,doses(idx));
end
%max dose only in or after first 24 hours
NODAL{ndoses+2} = rescueActivin(D,T,1);
NODAL{ndoses+2}(T<24) = 0;
NODAL{ndoses+3} = rescueActivin(D,T,1);
NODAL{ndoses+3}(T>24) = 0;

NODALlabels = cellfun(@(x) strrep(x,' ',''),strcat('A',cellstr(num2str(doses')))','UniformOutput',false);
NODALlabels = [{'WT Nodal'},NODALlabels,{'ACT after 24h','ACT first 24h'}];

%all relevant treatment conditions; each row is [BMP index, NODAL index]
treats = [[ones(ndoses,1),(2:ndoses+1)'];[2,1];[1,1];[ones(2,1),((ndoses+2):(ndoses+3))']];
nconds = size(treats,1);
labels = strcat(BMPlabels(treats(:,1))',repmat({'_'},size(treats,1),1),NODALlabels(treats(:,2))');
%group conditions for plotting as pie slices
pies = {[8,7],1:ndoses,[1,ndoses,10,9]};
for ii = 1:length(pies)
    fprintf('Pie %d:\n',ii)
    disp(labels(pies{ii}))
end

%% define parameters
%choose which model to use
%model 1 = simplest, model 4 = final, most complex
model = 4;
%specify the file containing parameter values
paramfile = 'paramValues';
load(fullfile(dataDir,paramfile))
opts = theta2opts(thetas(:,end),newopts,fields);

%protein names
channels = {'TFAP2C','SOX17','ME','AM'};
%adjust parameters depending on the chosen model
if model == 1
    %minimal model, just need to specify alphaX, alphaY, betaX, betaY, Kbx,
    %Kay, Kxy, and Kd
    opts.Kay = opts.Kaby;
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
order = [1 2 3]; %channel indices for [red green blue] components in plots

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
        vals(bi,:) = y(end,:);
    end
    V{idx} = vals;
end

%scale colors for display
maxvals = [0.9, 0.88, 0.7, 0.85];
V = cellfun(@(x) x./maxvals,V,'UniformOutput',false);

figure
P = cell(3,1);
imgs = vals2ims(V,DIST,edges);
for ii = 1:3
    subplot_tight(1,3,ii)
    P{ii} = pieSlices(imgs(pies{ii},:),R,order,channels);
end
set(gcf,'Position',[167 305 1724 965])

%% save pie plots
pref = ['model',num2str(model),'_',channels{order(1)},'_',...
    channels{order(2)},'_',channels{order(3)},'_'];
savenames = {'WT_CNTRLvLDN','NODALKO_doseResponse','NODALKO_ActDurations'};

for ii = 1:3
    savename = fullfile(savedir,[pref,savenames{ii},'.png']);
    imwrite(P{ii},savename)
end

%% Make a phase diagram of fate given specific BMP, Act conditions
%for evenly spaced BMP and Act between 0 and 1, simulate time course of
%differentiation 
Time = 200; %let run for an arbitrarily long time to reach steady state
tspan = [0,Time];
tvec = linspace(0,Time,1000);

bmpl = linspace(0,1,100);
actl = linspace(0,1,100);

nb = length(bmpl);
na = length(actl);

tic
X = zeros(nb*na,1); Y = X;
M = zeros(nb*na,4);
xx = zeros(size(M,1),1); yy = xx;
idx = 1;
for ii = 1:nb
    fprintf('.')
    if mod(ii,50) == 0
        fprintf('\n')
    end
    for jj = 1:na
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
im = zeros(nb,na,3);

order = [1 2 3];
for idx = 1:length(X)
    ii = xx(idx); jj = yy(idx);
    for ci = 1:length(order)
        im(ii,jj,ci) = M(idx,order(ci));
    end
end

figure
imagesc([0,1],[0,1],im)
axis square
set(gca,'YDir','normal') 
set(gca,'TickLength',[0,0])
set(gca,'Box','off')
cleanSubplot
xlabel('NODAL/Activin')
ylabel('BMP')
if length(order) == 2
    title(strcat("\color{red}", channels{order(1)},...
        " ", "\color{green}", channels{order(2)}))
elseif length(order) == 3
    title(strcat("\color{red}", channels{order(1)},...
        " ", "\color{green}", channels{order(2)},...
        " ", "\color{blue}", channels{order(3)}))
end


ylabs = [0 10 50]; ny = length(ylabs);
xlabs = [0 3 10 30 100]; nx = length(xlabs);
yticks = log(1 + [0 10 50])/log(51);
xticks = linspace(0,1,nx);
set(gca, 'YTick', yticks, 'YTickLabel', ylabs)
set(gca, 'XTick', xticks, 'XTickLabel', xlabs)

%% save figure
title('')
suffix = [];
for ci = 1:length(order)
    suffix = [suffix,'_',channels{order(ci)}]; %#ok<AGROW>
end
savename = ['model',num2str(model),'_phaseDiagram',suffix];
disp(savename)
saveas(gcf,fullfile(savedir,[savename,'.png']))

%% Kymographs of BMP and NODAL dynamics in WT BMP4-treated micropatterns
Time = 42;
d = linspace(0,R,1000); %edge distances
t = linspace(0,Time,100); %times
[D,T] = meshgrid(d,t);
bmp = dynamicBMP(D,T); %compute BMP wrt to time and distance
nodal = dynamicNodal(D,T); %compute NODAL wrt to time and distance
act = rescueActivin(D,T,1);

%quantities to plot wrt to time and edge distance
plots = {bmp,nodal,min((nodal + bmp),1),act};
nplots = length(plots);
names = {'BMP signaling','WT NODAL signaling','Combined signaling',...
    'NODALKO + Activin'};
snames = {'BMP_WT','NODAL_WT','Combined_WT','Activin_NODALKO'};

close all
figure
for ii = 1:nplots
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
    
    if ii == nplots
        h = colorbar;
        h.Label.String = 'Signal level';
    end
    
    savename = ['kymograph',snames{ii}];
    saveas(gcf,fullfile(savedir,[savename,'.png']))
    pause(0.5)
end

%% Simulation of EXP37 conditions - varying BMP, Act rescue doses
Time = 42;
tspan = [0,Time];
tvec = linspace(0,Time,1000);

edges = linspace(0,R,nbins+1); %bin edges
d = 0.5*(edges(1:end-1) + edges(2:end)); 
t = linspace(0,Time,1000); %times
[D,T] = meshgrid(d,t);

%simulated BMP and Activin doses on a log scale (reflects Smad activation
%levels)
BMPdoses = log(1 + [0, 10, 50, 100, 300])/log(51);
ACTdoses = [0, 0.4, 0.6, 0.8, 1]; %this experiment did not include 1 ng/mL activin
%actual BMP4 and Activin A concentrations in ng/mL
BMPlabels = strcat('B',{'0','10','50','100','300'});
ACTlabels = strcat('A',{'0','3','10','30','100'});

nbmp = length(BMPdoses);
nact = length(ACTdoses);
%generate relevant BMP and Activin spatiotemporal dynamics
BMP = cell(1,nbmp); ACT = cell(1,nact);
baseBMP = dynamicBMP(D,T);
for ii = 1:nbmp
    BMP{ii} = baseBMP*BMPdoses(ii);
end
for ii = 1:nact
    ACT{ii} = rescueActivin(D,T,ACTdoses(ii));
end
treats = [kron((1:5)',ones(5,1)),kron(ones(5,1),(1:5)')];
nconds = size(treats,1);

%iterate over treatment conditions in parallel
V = cell(nconds,1);
parfor idx = 1:nconds
    ii = treats(idx,1); jj = treats(idx,2);
    vals = zeros(nbins,4);
    %iterate over bins within the treatment condition
    for bi = 1:nbins
        bmp = BMP{ii}(:,bi)';
        act = ACT{jj}(:,bi)';
        if model == 1
            [~,y] = ode45(@(t,y) ODE_model1(t,y,tvec,bmp,act,opts), tspan, ic);
        else
            [~,y] = ode45(@(t,y) ODE_4component(t,y,tvec,bmp,act,opts), tspan, ic);
        end
        vals(bi,:) = y(end,:);
    end
    V{idx} = vals;
end

imC = vals2ims(V,DIST,edges);

%% Make EXP37 calculated pie plots and approx phase diagram
fs = 48;

close all
channels = {'TFAP2C','SOX17','ME','AM'};
order = [1 2 3]; %channels{order(1)} = red, channels{order(2)} = green, channels{order(3)} = blue

%make pie plots
BMPpies = cell(1,nbmp);
ACTpies = cell(1,nact);
for ii = 1:nbmp
    BMPpies{ii} = (ii-1)*nact+1:ii*nact;
end
for ii = 1:nact
    ACTpies{ii} = ii:nbmp:nconds;
end

allpies = {BMPpies,ACTpies};
alllabels = {BMPlabels,...
    ACTlabels};
figtitles = {'FixedBMP','FixedACT'};
P = cell(1,length(allpies));

for jj = 1:length(allpies)
    pref = ['EXP37_calculated_',figtitles{jj}];
    figure('Position',[50 250 2500 750])
    npies = length(allpies{jj});
    P = cell(1,npies);
    idx = 1;
    for ii = 1:npies
        subplot(1,npies,idx)
        P{ii} = pieSlices(imC(allpies{jj}{ii},:),R,order,channels);
        idx = idx + 1;
        xlabel(alllabels{jj}{ii})
        
        savename = fullfile(savedir,[pref,'_',alllabels{jj}{ii},sprintf('%dhrs.png',Time)]);
        imwrite(P{ii},savename)
    end
    saveas(gcf,[figtitles{jj},sprintf('_calculated_pies_%dhrs',Time)],gcf)
end

% make approx phase diagram based on average expression within 100 um of
% the edge in each condition
dmax = 100; %cutoff for distance from the colony edge
dbins = find(d < dmax);
dedges = edges([dbins,dbins(end)+1]);
im = zeros(nbmp,nact,3);

order = [1 2];
for idx = 1:nconds
    ii = treats(idx,1); jj = treats(idx,2);
    for ci = 1:length(order)
        %do a weighted average based on the area contained within each
        %radial bin to most closely match how we calculated this for real
        %data
        weight = 0; val = 0;
        for di = 1:length(dbins)
            w = (R - dedges(di+1))^2 - (R - dedges(di))^2;
            weight = weight + w;
            val = val + w*V{idx}(dbins(di),order(ci));
        end
        im(ii,jj,ci) = val/weight;
    end
end

%save plots with and without scaling between minimum and maximum values
for doscaling = [true false]
    test = im;
    if doscaling
        suff = '_scaled';
        for ii = 1:size(im,3)
            test(:,:,ii) = test(:,:,ii)/max(test(:,:,ii),[],'all');
        end
    else
        suff = '';
    end
    figure('Position',[975 400 950 800])
    imagesc([0,1],[0,1],test)
    axis square
    set(gca,'YDir','normal') 
    cleanSubplot(fs)
    xlabel('Activin')
    ylabel('BMP')
    if length(order) == 2
        title(strcat("\color{red}", channels{order(1)},...
            " ", "\color{green}", channels{order(2)}))
    elseif length(order) == 3
        title(strcat("\color{red}", channels{order(1)},...
            " ", "\color{green}", channels{order(2)},...
            " ", "\color{blue}", channels{order(3)}))
    end


    ylabs = [0 10 50 100 300]; ny = length(ylabs);
    xlabs = [0 3 10 30 100]; nx = length(xlabs);
    yticks = linspace(0,1,ny);
    xticks = linspace(0,1,nx);
    set(gca, 'YTick', yticks, 'YTickLabel', ylabs)
    set(gca, 'XTick', xticks, 'XTickLabel', xlabs)
    
    savename = fullfile(savedir,['avgEdgeExpression_TFAP2C_SOX17',suff,'.png']);
    saveas(gcf,fname);
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
%     sigma = 30;
    sigma = 60;
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

if ~exist('sigma','var')
%     sigma = 90;
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
% set(gcf,'Position',[950 650 660 684])

if length(order) == 2
    title(strcat("\color{red}", channels{order(1)},...
        " ", "\color{green}", channels{order(2)}))
elseif length(order) == 3
    title(strcat("\color{red}", channels{order(1)},...
        " ", "\color{green}", channels{order(2)},...
        " ", "\color{blue}", channels{order(3)}))
end

end
