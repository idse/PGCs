function res = plotRadialProfiles(stats, meta, options)
    % average profile with equal number of datapoints per 'bin'
    
    if ~isfield(options,'normalize')
        options.normalize = false;
    end
    if ~isfield(options,'nucChannels') && ~isfield(options,'cytChannels')
        error('please specify channels');
    end
    if ~isfield(options,'nucChannels')
        options.nucChannels = [];
    end
    if ~isfield(options,'cytChannels')
        options.cytChannels = [];
    end
    if ~isfield(options,'legend')
        options.legend = true;
    end
    if ~isfield(options,'colors')
        options.colors = lines(5);
        options.colors = options.colors([2 5 1 3 4],:);
    end
    if ~isfield(options,'std') %|| options.normalize == true
        options.std = false;
    end
    if ~isfield(options,'individualColonies')
        options.individualColonies = false;
    end
    if ~isfield(options,'Ilimmargin')
        options.Ilimmargin = 3;
    end
    if isfield(options,'FontSize') 
        fs = options.FontSize;
    else
        fs = 20;
    end
    if isfield(options,'conditionIdx')
        condi = options.conditionIdx;
    else
        error('please specify options.conditionIdx');
    end
    if isfield(options,'pointsPerBin') 
        ptsperbin = options.pointsPerBin;
    else
        ptsperbin = 400;
    end
    if ~isfield(options,'useTrueRadius') 
        options.useTrueRadius = true;
    end
    interpmethod = 'linear';
    
    % average with equal number of datapoints per 'bin' for all colonies
    % combined
    % this reduces artefacts due to small numbers of junk points on the
    % colony edge
    R = sqrt(sum(stats.XY{condi}.^2,2))*meta.xres;
    nL = stats.nucLevel{condi};
    cL = stats.cytLevel{condi};
    
    % set up a regularly spaced grid for interpolation
    dr = 3; % step in micron
    nominalradius = stats.radiusMicron{condi};
    Ngrid = round(nominalradius/dr);
    radialgrid = linspace(0, nominalradius, Ngrid);
    
    % calculate profiles of individual colonies
    %------------------------------------------------------------------
    nuc_profiles = {};
    nuc_profile_stds = {};
    cyt_profiles = {};
    cyt_profile_stds = {};
    
    for si = unique(stats.sample{condi})'
    
        sbi = stats.sample{condi} == si;
        Rs = R(sbi);
        nLs = nL(sbi,:);
        cLs = cL(sbi,:);

        % estimate cell density and colony radius (which may differ
        % slightly from nominal radius of micropattern)
        PI = 3.1415;

        bw = 5;
        [f,xi] = ksdensity(Rs,'Bandwidth',bw,'BoundaryCorrection','reflection');

        dx = xi(2)-xi(1);
        Ncells = numel(Rs);
        celldensity{si} = Ncells*f./(2*PI*xi);

        nominalradius = stats.radiusMicron{condi};
        xrangeinside = xi > bw & xi < 0.9*nominalradius;
        avgdensityinside = sum(2*PI*celldensity{si}(xrangeinside).*xi(xrangeinside)*dx)/sum(2*PI*xi(xrangeinside)*dx);
        % define the radius as the point where the density drops to 10% of
        % the mean inside
        fullxrange = xi > bw & celldensity{si} > 0.1*avgdensityinside;
        trueradius = max(xi(fullxrange));

        if options.useTrueRadius
            radius = trueradius;
            disp('adjusting for true radius of colony');
        else
            radius = nominalradius;
        end

        % calculate radial profiles of markers for individual colonies
        ptsperbinindividualcolonies = 10;
        [~,I] = sort(Rs);
        N = round(numel(I)/ptsperbinindividualcolonies);
        edges = round(linspace(1, numel(I), N));
        overlap = ptsperbinindividualcolonies; % overlap bins for smoothness

        r_tmp = zeros([N-1 1]);
        nuc_profile_tmp = zeros([N-1 4]);
        nuc_profile_std_tmp = zeros([N-1 4]);
        cyt_profile_tmp = zeros([N-1 4]);
        cyt_profile_std_tmp = zeros([N-1 4]);

        for i = 1:N-1
            ptidx = max(min(edges),edges(i)-overlap):min(edges(i+1)+overlap,max(edges));
            r_tmp(i) = mean(Rs(I(ptidx)));
            nuc_profile_tmp(i,:) = nanmean(nLs(I(ptidx),:));
            nuc_profile_std_tmp(i,:) = nanstd(nLs(I(ptidx),:));%/sqrt(ptsperbin);
            cyt_profile_tmp(i,:) = nanmean(cLs(I(ptidx),:));
            cyt_profile_std_tmp(i,:) = nanstd(cLs(I(ptidx),:));%/sqrt(ptsperbin);
        end

        % convert radius to edge distance
        r_tmp = radius - r_tmp;
        R(sbi) = radius - R(sbi);

        nuc_profiles{si} = zeros([Ngrid 4]);
        nuc_profile_stds{si} = zeros([Ngrid 4]);
        cyt_profiles{si} = zeros([Ngrid 4]);
        cyt_profile_stds{si} = zeros([Ngrid 4]);

        % interpolate everything on a fixed radial grid
        for ci = 1:meta.nChannels
            nuc_profiles{si}(:,ci) = interp1(r_tmp, nuc_profile_tmp(:,ci)', radialgrid,interpmethod,'extrap')';
            cyt_profiles{si}(:,ci) = interp1(r_tmp, cyt_profile_tmp(:,ci)', radialgrid,interpmethod,'extrap')';
            nuc_profile_stds{si}(:,ci) = interp1(r_tmp, nuc_profile_std_tmp(:,ci)', radialgrid,interpmethod,'extrap')';
            cyt_profile_stds{si}(:,ci) = interp1(r_tmp, cyt_profile_std_tmp(:,ci)', radialgrid,interpmethod,'extrap')';
        end
        celldensity{si} = interp1(xi, celldensity{si}, radialgrid,interpmethod,'extrap');
    end
    
    % average with equal number of datapoints per 'bin'
    %------------------------------------------------------------------
    [~,I] = sort(R);

    N = round(numel(I)/ptsperbin);
    edges = round(linspace(1, numel(I), N));
    overlap = round(ptsperbin/2);

    r_tmp = zeros([N-1 1]);
    nuc_profile_tmp = zeros([N-1 4]);
    nuc_profile_std_tmp = zeros([N-1 4]);
    cyt_profile_tmp = zeros([N-1 4]);
    cyt_profile_std_tmp = zeros([N-1 4]);

    for i = 1:N-1
        ptidx = max(min(edges),edges(i)-overlap):min(edges(i+1)+overlap,max(edges));
        r_tmp(i) = mean(R(I(ptidx)));
        nuc_profile_tmp(i,:) = nanmean(nL(I(ptidx),:));
        nuc_profile_std_tmp(i,:) = nanstd(nL(I(ptidx),:));%/sqrt(ptsperbin);
        cyt_profile_tmp(i,:) = nanmean(cL(I(ptidx),:));
        cyt_profile_std_tmp(i,:) = nanstd(cL(I(ptidx),:));%/sqrt(ptsperbin);
    end
    
    % interpolate on grid
    nuc_profile = zeros([Ngrid 4]);
    nuc_profile_std = zeros([Ngrid 4]);
    cyt_profile = zeros([Ngrid 4]);
    cyt_profile_std = zeros([Ngrid 4]);
    for ci = 1:meta.nChannels
        nuc_profile(:,ci) = interp1(r_tmp, nuc_profile_tmp(:,ci)', radialgrid, interpmethod, 'extrap')';
        cyt_profile(:,ci) = interp1(r_tmp, cyt_profile_tmp(:,ci)', radialgrid, interpmethod,'extrap')';
        nuc_profile_std(:,ci) = interp1(r_tmp, nuc_profile_std_tmp(:,ci)', radialgrid,interpmethod,'extrap')';
        cyt_profile_std(:,ci) = interp1(r_tmp, cyt_profile_std_tmp(:,ci)', radialgrid,interpmethod,'extrap')';
    end
    
    % alternative error bar: std of colonies
    % since std from 4 colonies is very noisy but true std should vary
    % smoothly in neighboring points we average the estimate of the std
    % between neighbors
    nuc_profile_colstd = imfilter(std(cat(3, nuc_profiles{:}),[],3),[1 1 1 1 1]'/5, 'replicate');
    cyt_profile_colstd = imfilter(std(cat(3, cyt_profiles{:}),[],3),[1 1 1 1 1]'/5, 'replicate');
    
    % determine limits for normalization
    if isfield(options,'nuclimits')
        nuclimits = options.nuclimits;
    else
        nuclimits = [min(nuc_profile(1:end,:))' max(nuc_profile(1:end,:))'];
    end
    if isfield(options,'cytlimits')
        cytlimits = options.cytlimits;
    else
        cytlimits = [min(cyt_profile)' max(cyt_profile)'];
    end

    if options.normalize
        for ci = 1:size(nuc_profile,2)
            nuc_profile(:,ci) = (nuc_profile(:,ci) - nuclimits(ci,1))./(nuclimits(ci,2) - nuclimits(ci,1));
            cyt_profile(:,ci) = (cyt_profile(:,ci) - cytlimits(ci,1))./(cytlimits(ci,2) - cytlimits(ci,1));
            nuc_profile_std(:,ci) = nuc_profile_std(:,ci)./(nuclimits(ci,2) - nuclimits(ci,1));
            cyt_profile_std(:,ci) = cyt_profile_std(:,ci)./(cytlimits(ci,2) - cytlimits(ci,1));
            nuc_profile_colstd(:,ci) = nuc_profile_colstd(:,ci)./(nuclimits(ci,2) - nuclimits(ci,1));
            cyt_profile_colstd(:,ci) = cyt_profile_colstd(:,ci)./(cytlimits(ci,2) - cytlimits(ci,1));
        end
    end

    r = radialgrid;
    
    legendentries = {};
    hold on
    if ~isempty(options.nucChannels)
        for cii = numel(options.nucChannels):-1:1
            ci = options.nucChannels(cii) + 1;
            
            plot(r, nuc_profile(:,ci),'LineWidth',3, 'Color', options.colors(cii,:))

            if isempty(options.cytChannels)
                prefix = [];
            else
                prefix = 'nuc ';
            end
            legendentries = [legendentries, [prefix meta.channelLabel{ci}]];
        end
        for cii = numel(options.nucChannels):-1:1
            ci = options.nucChannels(cii) + 1;
            if options.std
                errorbar(r, nuc_profile(:,ci), nuc_profile_colstd(:,ci),'LineWidth',1, 'Color', options.colors(cii,:)); 
            end
        end
        
        if options.individualColonies
            for pi = 1:numel(positions)
                for cii = 1:numel(options.nucChannels)
                    ci = options.nucChannels(cii) + 1;
                    plot(r(1:end-binmarg), positions(pi).radialProfile.NucAvgSeg(1:end-binmarg,ci), 'Color', options.colors(cii,:)); 
                end
            end
        end
    end
        
    if ~isempty(options.cytChannels)
        
        for cii = 1:numel(options.cytChannels)
            ci = options.cytChannels(cii) + 1;
            if options.std
                errorbar(r, cyt_profile(:,ci), cyt_profile_colstd(:,ci),'--','LineWidth',2, 'Color', options.colors(cii,:)); 
            else
                plot(r, cyt_profile(:,ci),'--','LineWidth',2, 'Color', options.colors(cii,:))
            end
            legendentries = [legendentries, ['cyt ' meta.channelLabel{ci}]];
        end
    end
    hold off
    
    if options.legend
        legend(legendentries, 'Location','NorthEast');
    end

    % make it pretty
    fgc = 'k';
    bgc = 'w';
    graphbgc = 1*[1 1 1]; 

    xlabel('edge distance ( um )', 'FontSize',fs,'FontWeight','Bold','Color',fgc)
    ylabel('intensity (a.u.)', 'FontSize',fs,'FontWeight','Bold','Color',fgc);

    set(gcf,'color',bgc);
    set(gca, 'LineWidth', 2);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    set(gca,'XColor',fgc);
    set(gca,'YColor',fgc);
    set(gca,'Color',graphbgc);
    
    res = struct('r',r, 'nuc_profile',nuc_profile, 'nuc_profile_std',nuc_profile_std,...
        'cyt_profile', cyt_profile, 'cyt_profile_std', cyt_profile_std,...
        'celldensity',{celldensity},'nuc_profile_colstd',nuc_profile_colstd,...
        'cyt_profile_colstd',cyt_profile_colstd);
end