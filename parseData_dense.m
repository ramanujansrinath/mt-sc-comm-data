function parsed = parseData_dense(data,zScore,eqTrials,rateMatch)
    if ~exist('zScore','var'); zScore = 0; end
    if ~exist('eqTrials','var'); eqTrials = 0; end
    if ~exist('rateMatch','var'); rateMatch = 0; end
    
    % for data
    base		= data.base;
    areacode	= data.areacode;
    first		= data.first;
    targ		= data.targ;
    targshort	= data.targshort;
    prevshort	= data.prevshort;
    
    % for trials
    attend		= data.attend;
    change		= data.change;
    outcome		= data.outcome;
    valid		= data.valid;
    instruct    = data.instruct;
    
    % for dense extract
    trialnum    = data.trialnum;
    stimnum     = data.stimnum;
    counts      = data.counts;
    alloutcome  = data.alloutcome;
    
    % for visual/visuomotor classification (vis +, mot -)
    visMotIdx = data.SCclassify.A1vmPIindex;
    
    %% neurons
    respratio=1.1; %inclusion threshold
    % NaN out channels in MT that dont respond more to the first stim than baseline
    MTbaselinecounts = nanmean(base(areacode==1,:),2)./5; %divide because baseline was stored as a rate, not a count (oops)
    MTratios = nanmean(first(areacode==1,:),2)./MTbaselinecounts;
    MTidx = MTratios>=respratio;
    % NaN out channels in SC that dont either respond more to first stim than
    % baseline, or target stim than baseline
    SCbaselinecounts = nanmean(base(areacode==2,:),2)./5;
    SCratios1 = nanmean(first(areacode==2,:),2)./SCbaselinecounts;
    SCidx1 = SCratios1>=respratio;
    SCratios2 = nanmean(targ(areacode==2,attend==1),2)./SCbaselinecounts;
    SCidx2 = SCratios2>=respratio;
    SCidx = SCidx1 | SCidx2;

    scInclude=[SCidx; zeros(length(MTidx),1)]==1; % + (areacode==2)')>=2);
    mtInclude=[zeros(length(SCidx),1); MTidx]==1; % + (areacode==1)')>=2);
    
    scVisIdx = visMotIdx(scInclude); % this is fine because the sc neurons are always before the mt neurons

    %% trials/stim
    outcomeToUse = 100;
    % stims = alloutcome == outcomeToUse;
    
    instruct = instruct(trialnum);
    attend = attend(trialnum);
    
    % good = alloutcome==outcomeToUse & instruct==0 & attend==1;
    good = stimnum>1 & instruct==0 & attend==1;
    % disp(['in ' num2str(sum(good))]);
    MTcountsIn=counts(mtInclude,good);
    SCcountsIn=counts(scInclude,good);
    
    % good = alloutcome==outcomeToUse & instruct==0 & attend==2;
    good = stimnum>1 & instruct==0 & attend==2;
    % disp(['out ' num2str(sum(good))]);
    MTcountsOut=counts(mtInclude,good);
    SCcountsOut=counts(scInclude,good);
    
    bad = isnan(MTcountsIn(1,:)) | isnan(SCcountsIn(1,:));
    MTcountsIn(:,bad) = [];
    SCcountsIn(:,bad) = [];
    
    bad = isnan(MTcountsOut(1,:)) | isnan(SCcountsOut(1,:));
    MTcountsOut(:,bad) = [];
    SCcountsOut(:,bad) = [];
    
    parsed.mt_in = MTcountsIn;
    parsed.sc_in = SCcountsIn;

    parsed.mt_out = MTcountsOut;
    parsed.sc_out = SCcountsOut;
    
    parsed.sc_visUnits = scVisIdx;
    
    nTrials = [size(parsed.mt_in,2),size(parsed.mt_out,2)];
    if eqTrials
        nTrials = min(nTrials);
        
        inTrials = sort(randperm(size(parsed.mt_in,2),nTrials));
        parsed.mt_in = parsed.mt_in(:,inTrials);
        parsed.sc_in = parsed.sc_in(:,inTrials);
        
        outTrials = sort(randperm(size(parsed.mt_out,2),nTrials));
        parsed.mt_out = parsed.mt_out(:,outTrials);
        parsed.sc_out = parsed.sc_out(:,outTrials);
    end
    
    disp(['... nTrials: ' num2str(nTrials)])
    
    if rateMatch
        parsed = matchRates(parsed);
    end
    
    if zScore
        parsed.mt_in = zScoreData(parsed.mt_in);
        parsed.sc_in = zScoreData(parsed.sc_in);
        parsed.mt_out = zScoreData(parsed.mt_out);
        parsed.sc_out = zScoreData(parsed.sc_out);
    end
end

function dat = zScoreData(dat)
    dat = (dat - repmat(mean(dat,2),1,size(dat,2))) ./ repmat(std(dat,[],2),1,size(dat,2));
end

function parsed = matchRates(parsed)
    demean_mean = @(x) mean(x - repmat(mean(x,2),1,size(x,2)));
    mt_in = demean_mean(parsed.mt_in);
    sc_in = demean_mean(parsed.sc_in);
    mt_out = demean_mean(parsed.mt_out);
    sc_out = demean_mean(parsed.sc_out);
    
    edges = linspace(min([mt_in sc_in mt_out sc_out]),max([mt_in sc_in mt_out sc_out]),101);
    binwidth = edges(2)-edges(1);
    
    % figure; subplot(121);
    % histogram(mt_in,edges(1:4:end),'DisplayStyle','stairs'); hold on;
    % histogram(sc_in,edges(1:4:end),'DisplayStyle','stairs');
    % histogram(mt_out,edges(1:4:end),'DisplayStyle','stairs');
    % histogram(sc_out,edges(1:4:end),'DisplayStyle','stairs');
    
    [inda,indb] = matchdistributions4({mt_in sc_in mt_out sc_out},edges,binwidth);
    
    parsed.mt_in = parsed.mt_in(:,inda);
    parsed.sc_in = parsed.sc_in(:,inda);
    parsed.mt_out = parsed.mt_out(:,indb);
    parsed.sc_out = parsed.sc_out(:,indb);
    
    % mt_in = demean_mean(parsed.mt_in);
    % sc_in = demean_mean(parsed.sc_in);
    % mt_out = demean_mean(parsed.mt_out);
    % sc_out = demean_mean(parsed.sc_out);
    % 
    % subplot(122);
    % histogram(mt_in,edges(1:4:end),'DisplayStyle','stairs'); hold on;
    % histogram(sc_in,edges(1:4:end),'DisplayStyle','stairs');
    % histogram(mt_out,edges(1:4:end),'DisplayStyle','stairs');
    % histogram(sc_out,edges(1:4:end),'DisplayStyle','stairs');
end

function [inda,indb] = matchdistributions4(dist,edges,binwidth)
    [~,~,n1,~,i1,~]=bindata(dist{1},dist{1},edges,binwidth);
    [~,~,n2,~,i2,~]=bindata(dist{2},dist{2},edges,binwidth);
    [~,~,n3,~,i3,~]=bindata(dist{3},dist{3},edges,binwidth);
    [~,~,n4,~,i4,~]=bindata(dist{4},dist{4},edges,binwidth);
    
    matchn = nanmin([n1;n2;n3;n4]);
    id1 = []; id2 = []; id3 = []; id4 = [];
    for ii=1:length(edges)-1
        id1 = [id1 datasample(i1{ii},matchn(ii),'Replace',false)];
        id2 = [id2 datasample(i2{ii},matchn(ii),'Replace',false)];
        id3 = [id3 datasample(i3{ii},matchn(ii),'Replace',false)];
        id4 = [id4 datasample(i4{ii},matchn(ii),'Replace',false)];
    end
    
    inda = intersect(id1,id2);
    indb = intersect(id3,id4);
    
    inda = sort(datasample(inda,min(length(inda),length(indb)),'Replace',false));
    indb = sort(datasample(indb,min(length(inda),length(indb)),'Replace',false));
    
    % [~,ind1,~,ind2] = matchdistributions(dist(:,1),dist(:,2),edges,binwidth);
    % inda = intersect(ind1,ind2);
    % [~,ind1,~,ind2] = matchdistributions(dist(:,3),dist(:,4),edges,binwidth);
    % indb = intersect(ind1,ind2);
    % [~,inda,~,indb] = matchdistributions(dist(:,inda),dist(:,indb),edges,binwidth);
end