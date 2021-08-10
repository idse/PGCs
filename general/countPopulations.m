function counts = countPopulations(positions, meta, stats, dataDir, combo, conditionsidx)

    positive = {};

    Ncells = zeros([meta.nPositions 1]);
    NcellsTotal = zeros([meta.nWells 1]);
    NcellsStd = zeros([meta.nWells 1]);
    
    positivefraction = zeros([meta.nPositions meta.nChannels]);
    positivefractionavg = zeros([meta.nWells meta.nChannels]);
    positivefractionstd = zeros([meta.nWells meta.nChannels]);

    positivemean = zeros([meta.nPositions meta.nChannels]);
    positivemedian = zeros([meta.nPositions meta.nChannels]);
    positivemax = zeros([meta.nPositions meta.nChannels]);
    
    positivefraction_combo = zeros([meta.nPositions 2 2 2]);
    positivefractionavg_combo = zeros([meta.nWells 2 2 2]);
    positivefractionstd_combo = zeros([meta.nWells 2 2 2]);
    
    for condi = 1:meta.nWells

        condPos = meta.posPerCondition*(condi-1)+1:meta.posPerCondition*condi;

        for pi = condPos
            nucLevel = positions(pi).cellData.nucLevel;
            Ncells(pi) = size(nucLevel,1);
            background = positions(pi).cellData.background;    
            nucLevel = nucLevel - background;
            for ci = 1:meta.nChannels
                positive{pi, ci} = nucLevel(:,ci) > stats.thresholds(ci);
                positivemean(pi,ci) = mean(nucLevel(positive{pi, ci},ci));
                positivemedian(pi,ci) = median(nucLevel(positive{pi, ci},ci));
                m = max(nucLevel(positive{pi, ci},ci));
                if isempty(m)
                    m = 0;
                end
                positivemax(pi,ci) = m;
                positivefraction(pi,ci) = sum(positive{pi, ci})/numel(positive{pi, ci});
            end

            positivefraction_combo(pi, 1, 1, 1) = sum(positive{pi, 2} & positive{pi, 3} & positive{pi, 4})/numel(positive{pi, combo(1)});
            positivefraction_combo(pi, 1, 1, 2) = sum(positive{pi, 2} & positive{pi, 3} & ~positive{pi, 4})/numel(positive{pi, combo(1)});
            positivefraction_combo(pi, 1, 2, 1) = sum(positive{pi, 2} & ~positive{pi, 3} & positive{pi, 4})/numel(positive{pi, combo(1)});
            positivefraction_combo(pi, 1, 2, 2) = sum(positive{pi, 2} & ~positive{pi, 3} & ~positive{pi, 4})/numel(positive{pi, combo(1)});
            positivefraction_combo(pi, 2, 1, 1) = sum(~positive{pi, 2} & positive{pi, 3} & positive{pi, 4})/numel(positive{pi, combo(1)});
            positivefraction_combo(pi, 2, 1, 2) = sum(~positive{pi, 2} & positive{pi, 3} & ~positive{pi, 4})/numel(positive{pi, combo(1)});
            positivefraction_combo(pi, 2, 2, 1) = sum(~positive{pi, 2} & ~positive{pi, 3} & positive{pi, 4})/numel(positive{pi, combo(1)});
            positivefraction_combo(pi, 2, 2, 2) = sum(~positive{pi, 2} & ~positive{pi, 3} & ~positive{pi, 4})/numel(positive{pi, combo(1)});
        end
    end

    for condi = 1:meta.nWells

        condPos = meta.posPerCondition*(condi-1)+1:meta.posPerCondition*condi;

        NcellsTotal(condi) = sum(Ncells(condPos));
        NcellsStd(condi) = std(Ncells(condPos));
        
        positivefractionavg(condi,:) = mean(positivefraction(condPos,:),1);
        positivefractionstd(condi,:) = std(positivefraction(condPos,:),1);

        positivefractionavg_combo(condi,:,:,:) = mean(positivefraction_combo(condPos,:,:,:),1);
        positivefractionstd_combo(condi,:,:,:) = std(positivefraction_combo(condPos,:,:,:),1);
    end

    counts = struct(    'Ncells',Ncells,...
                        'positivemean',positivemean,...
                        'positivemedian',positivemedian,...
                        'positivemax',positivemax,...
                        'positivefraction',positivefraction,...
                        'positivefractionavg',positivefractionavg,...
                        'positivefractionstd',positivefractionstd,...
                        'positivefraction_combo',positivefraction_combo,...
                        'positivefractionavg_combo',positivefractionavg_combo,...
                        'positivefractionstd_combo',positivefractionstd_combo);

    % VISUALIZE
    fgc = 'k';
    bgc = 'w';
    fs = 24;
    cases = 1:5;
    
    for casei = cases
    
        if casei == 1 % FRACTIONS OF EACH MARKER
            
            fnameprefix = 'positiveFractions_bar';
            legendstr = meta.channelLabel(2:4);
            vals = positivefractionavg(conditionsidx,2:4)'*100;
            errs = positivefractionstd(conditionsidx,2:4)'*100;
            ylabelstr = '+% of all cells';
            titlestr = [];
            
        elseif casei == 2 % NUMBERS OF EACH MARKER
            
            fnameprefix = 'positive_bar';
            legendstr = meta.channelLabel(2:4);
            vals = (diag(NcellsTotal(conditionsidx))*positivefractionavg(conditionsidx,2:4)/meta.posPerCondition)';
            errs = (diag(NcellsTotal(conditionsidx))*positivefractionstd(conditionsidx,2:4)/meta.posPerCondition)';
            ylabelstr = '+cells / colony';
            titlestr = [];
        
        elseif casei == 3 || casei == 4 % COMBO=SUBSET FRACTION || NUMBERS
            
            if numel(combo) == 3
                
                combosub = zeros([4 3]);
                combosub(:,combo(1)-1) = 1; % shift because we're not including DAPI here so 2->1
                combosub(:,combo(2)-1) = [1 1 2 2]; % 1 = +, 2 = -
                combosub(:,combo(3)-1) = [1 2 1 2];

                vals = [];
                errs = [];
                for i = 1:4 % +-, --, ++, -+

                    comboi = sub2ind([2 2 2],combosub(i,1),combosub(i,2),combosub(i,3));
                    dvals = positivefractionavg_combo(conditionsidx, comboi)';
                    derrs = positivefractionstd_combo(conditionsidx, comboi)';
                    
                    if casei == 4
                        dvals = NcellsTotal(conditionsidx)'.*dvals/meta.posPerCondition;
                        derrs = NcellsTotal(conditionsidx)'.*derrs/meta.posPerCondition;
                    elseif casei == 3
                        dvals = dvals*100;
                        derrs = derrs*100;
                    end
                    
                    vals = cat(1,vals, dvals);
                    errs = cat(1,errs, derrs);
                end

                legendstr = {[meta.channelLabel{combo(2)} '+' meta.channelLabel{combo(3)} '+'],...
                    [meta.channelLabel{combo(2)} '+' meta.channelLabel{combo(3)} '-'],...
                    [meta.channelLabel{combo(2)} '-' meta.channelLabel{combo(3)} '+'],...
                    [meta.channelLabel{combo(2)} '-' meta.channelLabel{combo(3)} '-']};
                
            elseif numel(combo)==2

                legendstr = {[meta.channelLabel{combo(2)} '+'],...
                    [meta.channelLabel{combo(2)} '-']};
                error('implement this');
            end
            
            if casei == 3
                fnameprefix = ['positiveFractionsCombo_' num2str(combo) '_bar'];
                ylabelstr = '% of all cells';
                titlestr = ['breakdown of ' meta.channelLabel{combo(1)} '+'];
            
            elseif casei == 4
                fnameprefix = ['positiveCombo_' num2str(combo) '_bar'];
                ylabelstr = 'cells / colony';
                titlestr = ['breakdown of ' meta.channelLabel{combo(1)} '+'];
            end
        
        elseif casei == 5 % CELL NUMBERS
            fnameprefix = 'cellnumbers_bar';
            vals = NcellsTotal(conditionsidx)'/meta.posPerCondition;
            errs = NcellsStd(conditionsidx)';
            ylabelstr = 'total cells / colony';
            legendstr = [];
        end

        %errorbar(vals', errs', 'LineWidth',2);
        errorbar_groups(vals, errs,...
                    'bar_names', meta.conditions(conditionsidx),...
                    'bar_width',0.75,'errorbar_width',0.5,...
                    'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',2});
        ylim([0 max(vals(:))+max(errs(:))]);
        legend(legendstr,'Location','NorthWest','FontSize',fs-10);
        if ~isempty(titlestr)
            title(titlestr);
        end
        
        ylabel(ylabelstr);
        set(gcf,'color',bgc);
        set(gca, 'LineWidth', 2);
        set(gca,'FontSize', fs)
        set(gca,'FontWeight', 'bold')
        set(gca,'XColor',fgc);
        set(gca,'YColor',fgc);
        set(gca,'Color',bgc);
        axis square
        saveas(gcf, fullfile(dataDir, [fnameprefix num2str(conditionsidx) '.png'])); 
        close;
    end
    
end