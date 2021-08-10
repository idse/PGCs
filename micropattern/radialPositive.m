function [P,xi,Pstd] = radialPositive(stats, positions, condi, meta, combos,order)
    % radial probability profile
    %
    % [P,xi] = radialProbability(stats, positions, condi, meta, combos)
    % 
    % condi : index of condition/well
    % combos = {[2 3],..} for expression two or more markers (more not
    % implemented)

    bw = 25; % bandwidth for kde

    if ~exist('order','var')
        order = 1:numel(stats.markerChannels);
    end
    channels = stats.markerChannels(order);
    
    colRadius = positions((condi-1)*meta.posPerCondition + 1).radiusMicron;
    
    Psample = [];
    sampleid = unique(stats.sample{condi});
    Nsamples = numel(sampleid);
        
    for si = 1:Nsamples+1 % last iteration combines all samples
        
        if si <= Nsamples
            sampleidx = stats.sample{condi} == sampleid(si);
        else
           sampleidx = true(size(stats.sample{condi}));
        end
            
        for ci = channels

            posidx = stats.nucLevel{condi}(:,ci) > stats.thresholds(ci);
            posidx = posidx & sampleidx;
            
            XY = stats.XY{condi};
            posXY = XY(posidx,:);

            R = sqrt(XY(sampleidx,1).^2 + XY(sampleidx,2).^2);
            posR = sqrt(XY(posidx,1).^2 + XY(posidx,2).^2);

            % histogram approach
            %[N,edges] = histcounts(R);
            %[Npos,~] = histcounts(posR,edges);
            %xi = edges(1:end-1) + diff(edges)/2;
            %bar(xi, Npos./N)
            %plot(xi, Npos./N)
            %nansum(Npos./N)

            if ci == 2 && si == 1
                [f,xi] = ksdensity(R,'BoundaryCorrection','reflection','Bandwidth',bw,'Support','positive');
                P = zeros(max(channels), numel(xi));
            else
                [f,~] = ksdensity(R,xi,'BoundaryCorrection','reflection','Bandwidth',bw,'Support','positive');
            end
            %dx = diff(xi);
            %dx = dx(1);
            %intf = dx*sum(f);
            if sum(posidx) > 0
                [fpos,~] = ksdensity(posR,xi,'BoundaryCorrection','reflection','Bandwidth',bw,'Support','positive');
            else
                fpos = xi*0;
                warning(['no positive cells for ' meta.channelLabel{ci}]);
            end
            %intfpos = dx*sum(fpos);

            P(ci,:) = (sum(posidx)/sum(sampleidx))*fpos./f;
        end
        
        if si == 1
            Psample = P;
        else
            Psample = cat(3,Psample,P);
        end
    end

    Pstd = zeros(size(Psample(:,:,1)));
    Pavg = zeros(size(Psample(:,:,1)));
    for ci = 1:4
        Pstd(ci,:) = std(Psample(ci,:,1:4),[],3);
        Pavg(ci,:) = mean(Psample(ci,:,1:4),3);
    end

    if exist('combos','var')
        
        Pc = zeros(numel(combos), numel(xi));

        Pcsample = [];
        sii = unique(stats.sample{condi});
        
        for si = 1:Nsamples+1 % last iteration combines all samples
            
            if si == Nsamples+1
               sampleidx = true(size(stats.sample{condi}));
            else
                sampleidx = stats.sample{condi} == sii(si);
            end
            
            for ci = 1:numel(combos)

                posidx = stats.nucLevel{condi}(:,combos{ci}(1)) > stats.thresholds(combos{ci}(1)) & ...
                            stats.nucLevel{condi}(:,combos{ci}(2)) > stats.thresholds(combos{ci}(2));
                posidx = posidx & sampleidx;
                
                if sum(posidx) > 0
                    XY = stats.XY{condi};
                    posXY = XY(posidx,:);

                    R = sqrt(XY(sampleidx,1).^2 + XY(sampleidx,2).^2);
                    posR = sqrt(XY(posidx,1).^2 + XY(posidx,2).^2);

                    [f,~] = ksdensity(R,xi,'BoundaryCorrection','reflection','Bandwidth',bw,'Support','positive');
                    [fpos,~] = ksdensity(posR,xi,'BoundaryCorrection','reflection','Bandwidth',bw,'Support','positive');

                    Pc(ci,:) = (sum(posidx)/sum(sampleidx))*fpos./f;
                else
                    warning(['no positive cells for combo ' num2str(combos{ci})]);
                end
            end
            
            if si == 1
                Pcsample = Pc;
            else
                Pcsample = cat(3,Pcsample,Pc);
            end
        end
        
        Pcstd = zeros(size(Pcsample(:,:,1)));
        Pcavg = zeros(size(Pcsample(:,:,1)));
        for ci = 1:numel(combos)
            Pcstd(ci,:) = std(Pcsample(ci,:,1:4),[],3);
            Pcavg(ci,:) = mean(Pcsample(ci,:,1:4),3);
        end
    end
    
    colors = lines(6);
    colors = colors([2 5 1 3 4 6],:);

    lw = 3;
    figure,
    hold on
    for i = 1:numel(channels)
        %plot(xi*meta.xres,P(channels(i),:),'LineWidth',lw,'Color',colors(i,:));
        plot(colRadius - xi*meta.xres,Pavg(channels(i),:),'LineWidth',lw,'Color',colors(i,:));
    end
    if exist('combos','var')
        for i = 1:numel(combos)
            %plot(xi*meta.xres,Pc(i,:),'LineWidth',lw,'Color',colors(numel(channels) + i,:));
            plot(colRadius - xi*meta.xres,Pcavg(i,:),'LineWidth',lw,'Color',colors(numel(channels) + i,:));
        end
    end
    % error bars
    for i = 1:numel(channels)
        %errorbar(xi*meta.xres,P(channels(i),:),Pstd(channels(i),:),'LineWidth',lw,'Color',colors(i,:));
        good = ~isnan(P(channels(i),:));
        fill([colRadius - meta.xres*xi(good),fliplr(colRadius-meta.xres*xi(good))],[P(channels(i),good) + Pstd(channels(i),good), fliplr(P(channels(i),good) - Pstd(channels(i),good))],colors(i,:),'FaceAlpha',0.2,'EdgeColor','none');
    end
    if exist('combos','var')
        for i = 1:numel(combos)
            good = ~isnan(Pc(i,:));
            fill([colRadius - meta.xres*xi(good),fliplr(colRadius-meta.xres*xi(good))],[Pc(i,good) + Pcstd(i,good), fliplr(Pc(i,good) - Pcstd(i,good))],colors(i,:),'FaceAlpha',0.2,'EdgeColor','none');
        end
    end
    hold off
    
    legendstr = meta.channelLabel(channels);
    for ci = 1:numel(combos)
        x = meta.channelLabel(combos{ci});
        legendstr = [legendstr {[x{1} ' & ' x{2}]}];
    end
    fs = 30;
    legend(legendstr,'FontSize',20,'Location','NorthEast')
    axis square
    fgc = 'k';
    bgc = 'w';
    graphbgc = 1*[1 1 1]; 

    xlim([0 colRadius]);
    xlabel('edge distance ( um )', 'FontSize',fs,'FontWeight','Bold','Color',fgc)
    ylabel('positive fraction ', 'FontSize',fs,'FontWeight','Bold','Color',fgc);

    set(gcf,'color',bgc);
    set(gca, 'LineWidth', 2);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    set(gca,'XColor',fgc);
    set(gca,'YColor',fgc);
    set(gca,'Color',graphbgc);
end