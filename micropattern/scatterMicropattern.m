function scatterMicropattern(stats, meta, dataDir, options)

if ~isfield(options,'channelCombos')
    channelCombos = {[3 2], [4 3], [4 2]};
else
    channelCombos = options.channelCombos;
end
channelThresholds = options.channelThresholds;

% potential different norm for different conditions (necessary when
% combining different datasets)
if ~isfield(options,'normIdx')
    normIdx = zeros([meta.nWells 1],'uint16') + 1;
else
    normIdx = options.normIdx;
end

if isfield(options, 'log1p')
    log1p = options.log1p;
else
    log1p = false;
end

if isfield(options, 'conditionsCombined')
    conditionsCombined = options.conditionsCombined;
else
    conditionsCombined = true;
end

cases = [1 2]; % case 1: color for edge distance, case 2: color for 3rd gene

for casei = cases
for conditionIdx = options.conditionIdx
    
    % norm by mean
    %norm = mean(stats.nucLevel{normIdx}); 

    % norm by mean of positive cells (more robust to changes in the mean due to
    % colony size)
    ni = normIdx(conditionIdx);
    posidx = stats.nucLevel{ni} < options.channelThresholds(ni,:);
    %norm = sum(stats.nucLevel{ni}.*posidx,1)./sum(posidx,1);
    norm = options.channelThresholds(ni,:)/(exp(1)-1); 
    %norm = mean(stats.nucLevel{ni}); 
    
    for comboi = 1:numel(channelCombos)
        
        figure('Position',[0 0 700 500])
        
        channelIdx = channelCombos{comboi};
        cidx3 = setdiff(2:4,channelIdx);
        if isfield(options, 'channelMax')
            xmax = options.channelMax(channelIdx(1));
            ymax = options.channelMax(channelIdx(2));
            zmax = options.channelMax(cidx3);
        else
            xmax = stats.lim{channelIdx(1)}(2)/norm(channelIdx(1));
            ymax = stats.lim{channelIdx(2)}(2)/norm(channelIdx(2));
            zmax = stats.lim{cidx3}(2)/norm(cidx3);
        end
        if size(channelThresholds,1) > 1
            xthresh = channelThresholds(conditionIdx,channelIdx(1))/norm(channelIdx(1));
            ythresh = channelThresholds(conditionIdx,channelIdx(2))/norm(channelIdx(2));
        else
            xthresh = channelThresholds(channelIdx(1))/norm(channelIdx(1));
            ythresh = channelThresholds(channelIdx(2))/norm(channelIdx(2));
        end

        XY = stats.XY{conditionIdx};
        dist = options.radiusMicron - sqrt(pdist2(XY,[0 0],'squaredeuclidean'))*meta.xres;

        %hist(dist,50)

        x = stats.nucLevel{conditionIdx}(:,channelIdx(1))/norm(channelIdx(1));
        y = stats.nucLevel{conditionIdx}(:,channelIdx(2))/norm(channelIdx(2));
        z = stats.nucLevel{conditionIdx}(:,cidx3)/norm(cidx3);
        
        %z = stats.nucLevel{conditionIdx}(:,4);
        %cmax = 200;
        %c(c > cmax) = cmax;
        %c = stats.sample{conditionIdx};
        
        xlims = stats.lim{channelIdx(1)}/norm(channelIdx(1));
        ylims = stats.lim{channelIdx(2)}/norm(channelIdx(2));
        xlims(2) = xmax;
        ylims(2) = ymax;
        
        if log1p
            x = log(1+x); 
            y = log(1+y);
            z = log(1+z);
            xmax = log(1+xmax);
            ymax = log(1+ymax);
            zmax = log(1+zmax);
            z(z>zmax) = zmax;
            xthresh = log(1+xthresh);
            ythresh = log(1+ythresh);
            xlims = log(1+xlims);
            ylims = log(1+ylims);
        else
            zcut = 6;
            z(z<0) = 0;
            z(z>zcut) = zcut;    
        end
        
        if casei == 1
            c = z;
            %c = stats.nucLevel{conditionIdx}(:,cidx3) > options.channelThresholds(ni,cidx3);
            suffix = 'scatter_z';
            collabel = meta.channelLabel{cidx3};
            
        elseif casei == 2
            c = dist;
            suffix = 'scatter_dist';
            collabel = 'distance from edge';
        end
        
        
        lw = 3;
        fs = 36;
        pfs = 26;
        s = 20;

        h = scatter(x,y, s, c,'filled');
        
        colormap(turbo)
        c = colorbar;
        c.Label.String = collabel;
        
        xlim(xlims)
        ylim(ylims)

        xlabel([meta.channelLabel{channelIdx(1)}])% ' (xm)'])
        ylabel([meta.channelLabel{channelIdx(2)}])% ' (xm)'])

        set(gca, 'LineWidth', 2);
        set(gca,'FontSize', fs)
        set(gca,'FontWeight', 'bold')
        set(gca,'Color','w');
        axis square;

        xticks([0 1 2 3 4]);
        yticks([0 1 2 3 4]);
        %ax = ancestor(h, 'axes');
        %ax.XAxis.Exponent = 3;
        %ax.YAxis.Exponent = 3;

        pm = 100*sum(y > ythresh & x <= xthresh)/numel(y);
        pp = 100*sum(y > ythresh & x > xthresh)/numel(y);
        mp = 100*sum(y <= ythresh & x > xthresh)/numel(y);
        mm = 100*sum(y <= ythresh & x <= xthresh)/numel(y);
        disp(['y+x-: ' num2str(pm,3) '%']);
        disp(['y+x+: ' num2str(pp,3) '%']);
        disp(['y-x+: ' num2str(mp,3) '%']);
        disp(['y-x-: ' num2str(mm,3) '%']);

        line(xthresh*[1 1],[0 ymax],'LineStyle','--','LineWidth',lw,'Color','k');
        
        fc = 'k';
        text(xthresh/2, ythresh/2, [num2str(mm,3) '%'],'FontSize',pfs,'Color',fc,'FontWeight', 'bold','Background','w');
        text(xthresh/2, ythresh + (ymax - ythresh)/2, [num2str(pm,3) '%'],'FontSize',pfs,'Color',fc,'FontWeight', 'bold','Background','w')
        text(xthresh + (xmax - xthresh)/2, ythresh/2, [num2str(mp,3) '%'],'FontSize',pfs,'Color',fc,'FontWeight', 'bold','Background','w')
        text(xthresh + (xmax - xthresh)/2, ythresh + (ymax - ythresh)/2, [num2str(pp,3) '%'],'FontSize',pfs,'Color',fc,'FontWeight', 'bold','Background','w')

        N = size(stats.nucLevel{conditionIdx},1);
        text(0.05*xmax, ymax*0.95, ['N=' num2str(N)],'FontSize',pfs,'Color',fc,'FontWeight', 'bold','Background','w')
        
        line([0 xmax],ythresh*[1 1],'LineStyle','--','LineWidth',lw,'Color','k');
        %title(meta.conditions{conditionIdx});

        saveas(gcf, fullfile(dataDir, [meta.channelLabel{channelIdx(1)} '_' meta.channelLabel{channelIdx(2)} '_' meta.conditions{conditionIdx} '_' suffix '.png'])); 
        close;
    end
end
end

if conditionsCombined

%------------------------------------------------------
% combine conditions in one plot and color for condition
%------------------------------------------------------

colors = lines(7);
colors = colors([1:3 5 4 6 7],:);
%colors = lines(numel(options.conditionIdx));

for comboi = 1:numel(channelCombos)
    
    figure('Position',[0 0 600 500])
    hold on;
    xall = [];
    yall = [];
    call = [];
    
    for i = 1:numel(options.conditionIdx)
        
        conditionIdx = options.conditionIdx(end+1-i);
        
        ni = normIdx(conditionIdx);
        posidx = stats.nucLevel{ni} < options.channelThresholds(ni,:);
        %norm = sum(stats.nucLevel{ni}.*posidx,1)./sum(posidx,1);
        norm = options.channelThresholds(ni,:)/(exp(1)-1);
        
        channelIdx = channelCombos{comboi};
        
        if isfield(options, 'channelMax')
            xmax = options.channelMax(channelIdx(1));
            ymax = options.channelMax(channelIdx(2));
        else
            xmax = stats.lim{channelIdx(1)}(2)/norm(channelIdx(1));
            ymax = stats.lim{channelIdx(2)}(2)/norm(channelIdx(2));
        end
        
        if size(channelThresholds,1) > 1
            xthresh = channelThresholds(conditionIdx,channelIdx(1))/norm(channelIdx(1));
            ythresh = channelThresholds(conditionIdx,channelIdx(2))/norm(channelIdx(2));
        else
            xthresh = channelThresholds(channelIdx(1))/norm(channelIdx(1));
            ythresh = channelThresholds(channelIdx(2))/norm(channelIdx(2));
        end

        XY = stats.XY{conditionIdx};
        dist = options.radiusMicron - sqrt(pdist2(XY,[0 0],'squaredeuclidean'))*meta.xres;

        %hist(dist,50)

        x = stats.nucLevel{conditionIdx}(:,channelIdx(1))/norm(channelIdx(1));
        y = stats.nucLevel{conditionIdx}(:,channelIdx(2))/norm(channelIdx(2));
        cidx3 = setdiff(2:4,channelIdx);
        z = stats.nucLevel{conditionIdx}(:,cidx3)/norm(cidx3);
        %z = stats.nucLevel{conditionIdx}(:,4);
        c = repmat(colors(i,:),[size(x,1) 1]);
        %cmax = 200;
        %c(c > cmax) = cmax;
        %c = stats.sample{conditionIdx};
        
        % this is for the code below to randomize points from different
        % sets (not used right now)
        xall = cat(1,xall,x);
        yall = cat(1,yall,y);
        call = cat(1,call,c);
        
        xlims = stats.lim{channelIdx(1)}/norm(channelIdx(1));
        ylims = stats.lim{channelIdx(2)}/norm(channelIdx(2));
        xlims(2) = xmax;
        ylims(2) = ymax;
        
        if log1p
            x = log(1+x); 
            y = log(1+y);
            xmax = log(1+xmax);
            ymax = log(1+ymax);
            xthresh = log(1+xthresh);
            ythresh = log(1+ythresh);
            xlims = log(1+xlims);
            ylims = log(1+ylims);
        end

        s = 10;
        ma = 1;
        h = scatter(x,y, s, c,'filled','MarkerFaceAlpha',ma);
    end
    
%     shuffleidx = randperm(size(xall,1));
%     xall = xall(shuffleidx);
%     yall = yall(shuffleidx);
%     call = call(shuffleidx,:);
%     
%     if log1p
%         xall = log(1+xall); 
%         yall = log(1+yall);
%     end
%     h = scatter(xall,yall, s, call,'filled','MarkerFaceAlpha',0.5);

    colormap(turbo)
    %c = colorbar;
    %c.Label.String = 'distance from edge';

    xlim(xlims)
    ylim(ylims)

    xlabel([meta.channelLabel{channelIdx(1)}])% ' (xm)'])
    ylabel([meta.channelLabel{channelIdx(2)}])% ' (xm)'])

    fs = 36;
    lw = 3;
    
    set(gca, 'LineWidth', 2);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    set(gca,'Color','w');
    axis square;
        
    line(xthresh*[1 1],[0 ymax],'LineStyle','--','LineWidth',lw,'Color','k');

    fc = 'k';
%         text(xthresh/2, ythresh/2, [num2str(mm,3) '%'],'FontSize',pfs,'Color',fc,'FontWeight', 'bold','Background','w');
%         text(xthresh/2, ythresh + (ymax - ythresh)/2, [num2str(pm,3) '%'],'FontSize',pfs,'Color',fc,'FontWeight', 'bold','Background','w')
%         text(xthresh + (xmax - xthresh)/2, ythresh/2, [num2str(mp,3) '%'],'FontSize',pfs,'Color',fc,'FontWeight', 'bold','Background','w')
%         text(xthresh + (xmax - xthresh)/2, ythresh + (ymax - ythresh)/2, [num2str(pp,3) '%'],'FontSize',pfs,'Color',fc,'FontWeight', 'bold','Background','w')

%         N = size(stats.nucLevel{conditionIdx},1);
%         text(0, ymax*0.9, ['N=' num2str(N)],'FontSize',pfs,'Color',fc,'FontWeight', 'bold','Background','w')

    line([0 xmax],ythresh*[1 1],'LineStyle','--','LineWidth',lw,'Color','k');
    %title(meta.conditions{conditionIdx});

    lfs = 20;
    if isfield(options,'legendstr')
        legend(options.legendstr,'FontSize',lfs,'Location','NorthWest');
    else
        legend(meta.conditions(options.conditionIdx(end:-1:1)),'FontSize',lfs,'Location','NorthWest');
    end
    hold off
    saveas(gcf, fullfile(dataDir, [meta.channelLabel{channelIdx(1)} '_' meta.channelLabel{channelIdx(2)} '_' num2str(options.conditionIdx) '_scatter_conditions.png'])); 
    close;
end

%------------------------------------------------------
% combine conditions in one plot and color for remaining marker or
% condition with shuffled points
%------------------------------------------------------

colors = lines(numel(options.conditionIdx));

for casei = cases

for comboi = 1:numel(channelCombos)
    
    figure('Position',[0 0 600 500])
    hold on;
    xall = [];
    yall = [];
    call = [];
    
    for i = 1:numel(options.conditionIdx)
        
        conditionIdx = options.conditionIdx(i);%end+1-i
        
        ni = normIdx(conditionIdx);
        %posidx = stats.nucLevel{ni} < options.channelThresholds(ni,:);
        %norm = sum(stats.nucLevel{ni}.*posidx,1)./sum(posidx,1);
        norm = options.channelThresholds(ni,:)/(exp(1)-1);
        
        channelIdx = channelCombos{comboi};
        if isfield(options, 'channelMax')
            xmax = options.channelMax(channelIdx(1));
            ymax = options.channelMax(channelIdx(2));
        else
            xmax = stats.lim{channelIdx(1)}(2)/norm(channelIdx(1));
            ymax = stats.lim{channelIdx(2)}(2)/norm(channelIdx(2));
        end
        if size(channelThresholds,1) > 1
            xthresh = channelThresholds(conditionIdx,channelIdx(1))/norm(channelIdx(1));
            ythresh = channelThresholds(conditionIdx,channelIdx(2))/norm(channelIdx(2));
        else
            xthresh = channelThresholds(channelIdx(1))/norm(channelIdx(1));
            ythresh = channelThresholds(channelIdx(2))/norm(channelIdx(2));
        end

        XY = stats.XY{conditionIdx};
        dist = options.radiusMicron - sqrt(pdist2(XY,[0 0],'squaredeuclidean'))*meta.xres;

        %hist(dist,50)

        x = stats.nucLevel{conditionIdx}(:,channelIdx(1))/norm(channelIdx(1));
        y = stats.nucLevel{conditionIdx}(:,channelIdx(2))/norm(channelIdx(2));
        cidx3 = setdiff(2:4,channelIdx);
        
        z = stats.nucLevel{conditionIdx}(:,cidx3)/norm(cidx3);
        zcut = 6;
        z(z<0) = 0;
        z(z>zcut) = zcut;
        
        %zscaled = mat2gray(z);
        %zcol = imadjust(zscaled, stretchlim(zscaled, [0 0.999]));
        
        if casei == 1
            c = z;
            suffix = 'scatter3_conditions';
        elseif casei == 2
            c = repmat(colors(i,:),[size(x,1) 1]);
            disp([meta.conditions{conditionIdx} num2str(colors(i,:))]);
            suffix = 'scatter_conditions_shuffle';
        end
        %cmax = 200;
        %c(c > cmax) = cmax;
        %c = stats.sample{conditionIdx};
        
        % this is for the code below to randomize points from different
        % sets 
        xall = cat(1,xall,x);
        yall = cat(1,yall,y);
        call = cat(1,call,c);
        
        xlims = stats.lim{channelIdx(1)}/norm(channelIdx(1));
        ylims = stats.lim{channelIdx(2)}/norm(channelIdx(2));
        xlims(2) = xmax;
        ylims(2) = ymax;
    end
    
    shuffleidx = randperm(size(xall,1));
    xall = xall(shuffleidx);
    yall = yall(shuffleidx);
    call = call(shuffleidx,:);
    
    if log1p
        xall = log(1+xall); 
        yall = log(1+yall);
        xmax = log(1+xmax);
        ymax = log(1+ymax);
        xthresh = log(1+xthresh);
        ythresh = log(1+ythresh);
        xlims = log(1+xlims);
        ylims = log(1+ylims);
    end
  
    h = scatter(xall,yall, s, call,'filled','MarkerFaceAlpha',ma);

    xlim(xlims)
    ylim(ylims)

    xlabel([meta.channelLabel{channelIdx(1)}])% ' (xm)'])
    ylabel([meta.channelLabel{channelIdx(2)}])% ' (xm)'])

    fs = 36;
    
    set(gca, 'LineWidth', 2);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    set(gca,'Color','w');
    axis square;

    if casei == 1
        colormap(turbo)
        c = colorbar;
        c.Label.String = meta.channelLabel{cidx3};
    end
    %title(meta.conditions{conditionIdx});

    line(xthresh*[1 1],[0 ymax],'LineStyle','--','LineWidth',lw,'Color','k');
    line([0 xmax],ythresh*[1 1],'LineStyle','--','LineWidth',lw,'Color','k');
    
    hold off
    saveas(gcf, fullfile(dataDir, [meta.channelLabel{channelIdx(1)} '_' meta.channelLabel{channelIdx(2)} '_' num2str(options.conditionIdx) '_' suffix '.png'])); 
    close;
end
end
end

end