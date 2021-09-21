%% clear and load data 
clear; clc; close all;
addpath(genpath('dep'));
cols = [1 0.5 0.2; 0.2 0.5 1; 0.5 0.8 0.2; 0.8 0.1 0.8];

%load data
nRands = 20;
if ~exist(['data/randSplitsOnline_' num2str(nRands) '.mat'],'file')
    basePath = 'data/raw/';
    files = dir(basePath);
    files = {files(cellfun(@(x) contains(x,'reduced.mat'),{files.name})).name};
    files = cellfun(@(x) x(1:strfind(x,'_')-1),files,'UniformOutput',false);
    for ii=1:length(files)
        % get data
        disp([num2str(ii) ': ' files{ii}])
        data = load([basePath files{ii} '_reduced.mat']);
        data.sessionName = files{ii};
        parsed = parseData_dense(data,0,1,0);
        nMT = size(parsed.mt_in,1);
        nSC = size(parsed.sc_in,1);
        nTrials = size(parsed.mt_in,2);
        
        disp(['... nTrials: ' num2str(nTrials)])
        disp(['... nMT: ' num2str(nMT)])
        disp(['... nSC: ' num2str(nSC)])

        % iterate and plot iterations
        fprintf('... split: %4d/%4d',1,nRands);
        for kk=1:nRands
            fprintf('\b\b\b\b\b\b\b\b\b%4d/%4d',kk,nRands);
            rng('default'); rng('shuffle');
            
            % split mt population
            nMtSrcNeurons = floor(nMT/2);
            mtSrcNeuronIds = sort(randperm(nMT,nMtSrcNeurons));
            mtTargNeuronIds = find(~ismember(1:nMT,mtSrcNeuronIds));
            
            % split sc population
            nScSrcNeurons = floor(nSC/2);
            scSrcNeuronIds = sort(randperm(nSC,nScSrcNeurons));
            scTargNeuronIds = find(~ismember(1:nSC,scSrcNeuronIds));

            mt_sc(ii).in(kk) = subspace_analyse(parsed.sc_in(scTargNeuronIds,:),parsed.mt_in(mtSrcNeuronIds,:)); %#ok<*SAGROW>
            mt_mt(ii).in(kk) = subspace_analyse(parsed.mt_in(mtTargNeuronIds,:),parsed.mt_in(mtSrcNeuronIds,:));
            sc_mt(ii).in(kk) = subspace_analyse(parsed.mt_in(mtTargNeuronIds,:),parsed.sc_in(scSrcNeuronIds,:));
            sc_sc(ii).in(kk) = subspace_analyse(parsed.sc_in(scTargNeuronIds,:),parsed.sc_in(scSrcNeuronIds,:));
            
            mt_sc(ii).out(kk) = subspace_analyse(parsed.sc_out(scTargNeuronIds,:),parsed.mt_out(mtSrcNeuronIds,:));
            mt_mt(ii).out(kk) = subspace_analyse(parsed.mt_out(mtTargNeuronIds,:),parsed.mt_out(mtSrcNeuronIds,:));
            sc_mt(ii).out(kk) = subspace_analyse(parsed.mt_out(mtTargNeuronIds,:),parsed.sc_out(scSrcNeuronIds,:));
            sc_sc(ii).out(kk) = subspace_analyse(parsed.sc_out(scTargNeuronIds,:),parsed.sc_out(scSrcNeuronIds,:));
        end
        
        % calculate averages
        av(ii).mt_mt_in_perf = nanmean([mt_mt(ii).in.optPerf]);
        av(ii).sc_sc_in_perf = nanmean([sc_sc(ii).in.optPerf]);
        av(ii).mt_sc_in_perf = nanmean([mt_sc(ii).in.optPerf]);
        av(ii).sc_mt_in_perf = nanmean([sc_mt(ii).in.optPerf]);
        
        av(ii).mt_mt_out_perf = nanmean([mt_mt(ii).out.optPerf]);
        av(ii).sc_sc_out_perf = nanmean([sc_sc(ii).out.optPerf]);
        av(ii).mt_sc_out_perf = nanmean([mt_sc(ii).out.optPerf]);
        av(ii).sc_mt_out_perf = nanmean([sc_mt(ii).out.optPerf]);
        
        av(ii).mt_mt_in_perfdim = nanmean([mt_mt(ii).in.optDim]);
        av(ii).sc_sc_in_perfdim = nanmean([sc_sc(ii).in.optDim]);
        av(ii).mt_sc_in_perfdim = nanmean([mt_sc(ii).in.optDim]);
        av(ii).sc_mt_in_perfdim = nanmean([sc_mt(ii).in.optDim]);
        
        av(ii).mt_mt_out_perfdim = nanmean([mt_mt(ii).out.optDim]);
        av(ii).sc_sc_out_perfdim = nanmean([sc_sc(ii).out.optDim]);
        av(ii).mt_sc_out_perfdim = nanmean([mt_sc(ii).out.optDim]);
        av(ii).sc_mt_out_perfdim = nanmean([sc_mt(ii).out.optDim]);
        
        av(ii).mt_in_dim = nanmean([mt_mt(ii).in.srcDim]);
        av(ii).sc_in_dim = nanmean([sc_sc(ii).in.srcDim]);
        av(ii).mt_out_dim = nanmean([mt_mt(ii).out.srcDim]);
        av(ii).sc_out_dim = nanmean([sc_sc(ii).out.srcDim]);
        
        av(ii).mt_in_rsc = mean(nanmean([mt_mt(ii).in.srcRsc]));
        av(ii).sc_in_rsc = mean(nanmean([sc_sc(ii).in.srcRsc]));
        av(ii).mt_out_rsc = mean(nanmean([mt_mt(ii).out.srcRsc]));
        av(ii).sc_out_rsc = mean(nanmean([sc_sc(ii).out.srcRsc]));
        av(ii).in_rscSh = mean(nanmean([mt_sc(ii).in.sharedRsc]));
        av(ii).out_rscSh = mean(nanmean([mt_sc(ii).out.sharedRsc]));
        
        av(ii).mt_in_rate = mean(nanmean([mt_mt(ii).in.srcRate]));
        av(ii).sc_in_rate = mean(nanmean([sc_sc(ii).in.srcRate]));
        av(ii).mt_out_rate = mean(nanmean([mt_mt(ii).out.srcRate]));
        av(ii).sc_out_rate = mean(nanmean([sc_sc(ii).out.srcRate]));
        
        fprintf('\n')
    end
    save(['data/randSplitsOnline_' num2str(nRands) '.mat'],'mt_sc','mt_mt','sc_sc','sc_mt','av')
else
    load(['data/randSplitsOnline_' num2str(nRands) '.mat'],'mt_sc','mt_mt','sc_sc','sc_mt','av')
end

%% a little post processing 
for ii=1:15
    ratio.perf_mt_sc(ii) = mean([mt_sc(ii).in.optPerf])./mean([mt_sc(ii).out.optPerf]);
    ratio.perf_mt_mt(ii) = mean([mt_mt(ii).in.optPerf])./mean([mt_mt(ii).out.optPerf]);
    ratio.perf_sc_sc(ii) = mean([sc_sc(ii).in.optPerf])./mean([sc_sc(ii).out.optPerf]);
    ratio.perf_sc_mt(ii) = mean([sc_mt(ii).in.optPerf])./mean([sc_mt(ii).out.optPerf]);
    
    ratio.dim_mt_sc(ii) = mean([mt_sc(ii).in.optDim])./mean([mt_sc(ii).out.optDim]);
    ratio.dim_mt_mt(ii) = mean([mt_mt(ii).in.optDim])./mean([mt_mt(ii).out.optDim]);
    ratio.dim_sc_sc(ii) = mean([sc_sc(ii).in.optDim])./mean([sc_sc(ii).out.optDim]);
    ratio.dim_sc_mt(ii) = mean([sc_mt(ii).in.optDim])./mean([sc_mt(ii).out.optDim]);
    
    ratio.perf_mt_sc_perSes(ii) = mean([mt_sc(ii).in.optPerf]./[mt_sc(ii).out.optPerf]);
    ratio.perf_mt_mt_perSes(ii) = mean([mt_mt(ii).in.optPerf]./[mt_mt(ii).out.optPerf]);
    ratio.perf_sc_sc_perSes(ii) = mean([sc_sc(ii).in.optPerf]./[sc_sc(ii).out.optPerf]);
    ratio.perf_sc_mt_perSes(ii) = mean([sc_mt(ii).in.optPerf]./[sc_mt(ii).out.optPerf]);
    ratio.dim_mt_sc_perSes(ii)  = mean([mt_sc(ii).in.optDim]./[mt_sc(ii).out.optDim]);
    ratio.dim_mt_mt_perSes(ii)  = mean([mt_mt(ii).in.optDim]./[mt_mt(ii).out.optDim]);
    ratio.dim_sc_sc_perSes(ii)  = mean([sc_sc(ii).in.optDim]./[sc_sc(ii).out.optDim]);
    ratio.dim_sc_mt_perSes(ii)  = mean([sc_mt(ii).in.optDim]./[sc_mt(ii).out.optDim]);
    
    stat.perfIdx_mt_sc(ii) = (mean([mt_sc(ii).in.optPerf])-mean([mt_sc(ii).out.optPerf]))./(mean([mt_sc(ii).in.optPerf])+mean([mt_sc(ii).out.optPerf]));
    stat.perfIdx_mt_mt(ii) = (mean([mt_mt(ii).in.optPerf])-mean([mt_mt(ii).out.optPerf]))./(mean([mt_mt(ii).in.optPerf])+mean([mt_mt(ii).out.optPerf]));
    stat.perfIdx_sc_sc(ii) = (mean([sc_sc(ii).in.optPerf])-mean([sc_sc(ii).out.optPerf]))./(mean([sc_sc(ii).in.optPerf])+mean([sc_sc(ii).out.optPerf]));
    stat.perfIdx_sc_mt(ii) = (mean([sc_mt(ii).in.optPerf])-mean([sc_mt(ii).out.optPerf]))./(mean([sc_mt(ii).in.optPerf])+mean([sc_mt(ii).out.optPerf]));
    
    stat.perfIdx_mt_sc_perSes(ii,:) = ([mt_sc(ii).in.optPerf]-[mt_sc(ii).out.optPerf])./([mt_sc(ii).in.optPerf]+[mt_sc(ii).out.optPerf]);
    stat.perfIdx_mt_mt_perSes(ii,:) = ([mt_mt(ii).in.optPerf]-[mt_mt(ii).out.optPerf])./([mt_mt(ii).in.optPerf]+[mt_mt(ii).out.optPerf]);
    stat.perfIdx_sc_sc_perSes(ii,:) = ([sc_sc(ii).in.optPerf]-[sc_sc(ii).out.optPerf])./([sc_sc(ii).in.optPerf]+[sc_sc(ii).out.optPerf]);
    stat.perfIdx_sc_mt_perSes(ii,:) = ([sc_mt(ii).in.optPerf]-[sc_mt(ii).out.optPerf])./([sc_mt(ii).in.optPerf]+[sc_mt(ii).out.optPerf]);
    
    stat.dimIdx_mt_sc_perSes(ii,:) = ([mt_sc(ii).in.optDim]-[mt_sc(ii).out.optDim])./([mt_sc(ii).in.optDim]+[mt_sc(ii).out.optDim]);
    stat.dimIdx_mt_mt_perSes(ii,:) = ([mt_mt(ii).in.optDim]-[mt_mt(ii).out.optDim])./([mt_mt(ii).in.optDim]+[mt_mt(ii).out.optDim]);
    stat.dimIdx_sc_sc_perSes(ii,:) = ([sc_sc(ii).in.optDim]-[sc_sc(ii).out.optDim])./([sc_sc(ii).in.optDim]+[sc_sc(ii).out.optDim]);
    stat.dimIdx_sc_mt_perSes(ii,:) = ([sc_mt(ii).in.optDim]-[sc_mt(ii).out.optDim])./([sc_mt(ii).in.optDim]+[sc_mt(ii).out.optDim]);

    % stat.perf_mt_sc(ii) = signrank([mt_sc(ii).in.optPerf],[mt_sc(ii).out.optPerf]);
    % stat.perf_mt_mt(ii) = signrank([mt_mt(ii).in.optPerf],[mt_mt(ii).out.optPerf]);
    % stat.perf_sc_sc(ii) = signrank([sc_sc(ii).in.optPerf],[sc_sc(ii).out.optPerf]);
    % stat.perf_sc_mt(ii) = signrank([sc_mt(ii).in.optPerf],[sc_mt(ii).out.optPerf]);
    % stat.dim_mt_sc(ii)  = signrank([mt_sc(ii).in.optDim],[mt_sc(ii).out.optDim]);
    % stat.dim_mt_mt(ii)  = signrank([mt_mt(ii).in.optDim],[mt_mt(ii).out.optDim]);
    % stat.dim_sc_sc(ii)  = signrank([sc_sc(ii).in.optDim],[sc_sc(ii).out.optDim]);
    % stat.dim_sc_mt(ii)  = signrank([sc_mt(ii).in.optDim],[sc_mt(ii).out.optDim]);
end

%% fig 4 - overlapping 
clf; set(gcf,'color','w','pos',[53,1,857,704]);
subplot(231); hold on;
plot([0 1],[0 1],'k--','LineWidth',2);
plot([av.mt_mt_in_perf],[av.mt_mt_out_perf],'.','markersize',20,'color',cols(2,:));
plot([av.sc_sc_in_perf],[av.sc_sc_out_perf],'.','markersize',20,'color',cols(3,:));
plot([av.sc_mt_in_perf],[av.sc_mt_out_perf],'.','markersize',20,'color',cols(4,:));
plot([av.mt_sc_in_perf],[av.mt_sc_out_perf],'.','markersize',20,'color',cols(1,:));
fixPlot(gca,[0 1],[0 1],'attend in','attend out',0:0.25:1,0:0.25:1,'pred acc')

subplot(232); hold on;
plot([0 6],[0 6],'k--','LineWidth',2);
plot([av.mt_mt_in_perfdim],[av.mt_mt_out_perfdim],'.','markersize',20,'color',cols(2,:));
plot([av.sc_sc_in_perfdim],[av.sc_sc_out_perfdim],'.','markersize',20,'color',cols(3,:));
plot([av.sc_mt_in_perfdim],[av.sc_mt_out_perfdim],'.','markersize',20,'color',cols(4,:));
plot([av.mt_sc_in_perfdim],[av.mt_sc_out_perfdim],'.','markersize',20,'color',cols(1,:));
fixPlot(gca,[0 6],[0 6],'attend in','attend out',0:2:6,0:2:6,'pred dim')

subplot(234); hold on;
plot([0 4],[0 4],'k--','LineWidth',2);
plot([av.mt_mt_in_perf]./[av.mt_mt_out_perf],[av.mt_mt_in_perfdim]./[av.mt_mt_out_perfdim],'.','markersize',20,'color',cols(2,:));
plot([av.sc_sc_in_perf]./[av.sc_sc_out_perf],[av.sc_sc_in_perfdim]./[av.sc_sc_out_perfdim],'.','markersize',20,'color',cols(3,:));
plot([av.sc_mt_in_perf]./[av.sc_mt_out_perf],[av.sc_mt_in_perfdim]./[av.sc_mt_out_perfdim],'.','markersize',20,'color',cols(4,:));
plot([av.mt_sc_in_perf]./[av.mt_sc_out_perf],[av.mt_sc_in_perfdim]./[av.mt_sc_out_perfdim],'.','markersize',20,'color',cols(1,:));
fixPlot(gca,[0 4],[0 4],'pred (in/out)','preddim (in/out)',1:3,1:3,'x: pred; y: outdim')

subplot(235); hold on;
histogram([av.mt_mt_in_perf]./[av.mt_mt_out_perf],linspace(0,3,15),'DisplayStyle','stairs','EdgeColor',cols(2,:),'LineWidth',2)
histogram([av.sc_sc_in_perf]./[av.sc_sc_out_perf],linspace(0,3,15),'DisplayStyle','stairs','EdgeColor',cols(3,:),'LineWidth',2)
histogram([av.sc_mt_in_perf]./[av.sc_mt_out_perf],linspace(0,3,15),'DisplayStyle','stairs','EdgeColor',cols(4,:),'LineWidth',2)
histogram([av.mt_sc_in_perf]./[av.mt_sc_out_perf],linspace(0,3,15),'DisplayStyle','stairs','EdgeColor',cols(1,:),'LineWidth',2)
fixPlot(gca,[0 4],[0 12],'pred','count',1:3,0:5:12,'pred')

subplot(236); hold on;
histogram([av.mt_mt_in_perfdim]./[av.mt_mt_out_perfdim],linspace(0,3,15),'DisplayStyle','stairs','EdgeColor',cols(2,:),'LineWidth',2)
histogram([av.sc_sc_in_perfdim]./[av.sc_sc_out_perfdim],linspace(0,3,15),'DisplayStyle','stairs','EdgeColor',cols(3,:),'LineWidth',2)
histogram([av.sc_mt_in_perfdim]./[av.sc_mt_out_perfdim],linspace(0,3,15),'DisplayStyle','stairs','EdgeColor',cols(4,:),'LineWidth',2)
histogram([av.mt_sc_in_perfdim]./[av.mt_sc_out_perfdim],linspace(0,3,15),'DisplayStyle','stairs','EdgeColor',cols(1,:),'LineWidth',2)
fixPlot(gca,[0 4],[0 12],'pred','count',1:3,0:5:12,'preddim')

%% fig 4 - pred dims and pred accuracy as ratios of means 
clf; set(gcf,'color','w','pos',[53,1,857,704]);
subplot(221)
plot(ratio.perf_mt_sc,ratio.dim_mt_sc,'.','markersize',15,'color',cols(1,:)); hold on
plot(mean(ratio.perf_mt_sc),mean(ratio.dim_mt_sc),'k+','markersize',15,'LineWidth',2);
fixPlot(gca,[0 3],[0 3],'pred accuracy ratio','pred dims ratio',0:4,0:3,'mt => sc')

subplot(222)
plot(ratio.perf_sc_mt,ratio.dim_sc_mt,'.','markersize',15,'color',cols(4,:)); hold on
plot(mean(ratio.perf_sc_mt),mean(ratio.dim_sc_mt),'k+','markersize',15,'LineWidth',2);
fixPlot(gca,[0 3],[0 3],'pred accuracy ratio','pred dims ratio',0:4,0:3,'sc => mt')

subplot(223)
plot(ratio.perf_mt_mt,ratio.dim_mt_mt,'.','markersize',15,'color',cols(2,:)); hold on
plot(mean(ratio.perf_mt_mt),mean(ratio.dim_mt_mt),'k+','markersize',15,'LineWidth',2);
fixPlot(gca,[0 3],[0 3],'pred accuracy ratio','pred dims ratio',0:4,0:3,'mt => mt')

subplot(224)
plot(ratio.perf_sc_sc,ratio.dim_sc_sc,'.','markersize',15,'color',cols(3,:)); hold on
plot(mean(ratio.perf_sc_sc),mean(ratio.dim_sc_sc),'k+','markersize',15,'LineWidth',2);
fixPlot(gca,[0 3],[0 3],'pred accuracy ratio','pred dims ratio',0:4,0:3,'sc => sc')

%% supp fig 4a - pred accuracy 
clf; set(gcf,'color','w','pos',[53,1,857,704]);
CI95 = repmat(tinv([0.025 0.975],nRands-1)',1,15); 
subplot(221)
line([0 6],[0 6],[0 6],'color','k','linestyle','--','linewidth',2); hold on
me_in  = arrayfun(@(ii) mean([mt_sc(ii).in.optPerf]),1:15);
me_out = arrayfun(@(ii) mean([mt_sc(ii).out.optPerf]),1:15);
st_in  = repmat(arrayfun(@(ii) std([mt_sc(ii).in.optPerf])./sqrt(nRands),1:15),2,1) .* CI95;
st_out = repmat(arrayfun(@(ii) std([mt_sc(ii).out.optPerf])./sqrt(nRands),1:15),2,1) .* CI95;
arrayfun(@(ii) line([me_in(ii)+st_in(1,ii) me_in(ii)+st_in(2,ii)],[me_out(ii) me_out(ii)],'color',cols(1,:),'linewidth',2),1:15);
arrayfun(@(ii) line([me_in(ii) me_in(ii)],[me_out(ii)+st_out(1,ii) me_out(ii)+st_out(2,ii)],'color',cols(1,:),'linewidth',2),1:15);
plot(me_in,me_out,'o','markersize',8,'color',cols(1,:),'MarkerFaceColor','w','LineWidth',2)
plot(mean(me_in),mean(me_out),'k+','markersize',15,'LineWidth',2)
fixPlot(gca,[0 0.3],[0 0.3],'attend in','attend out',0:0.1:1,0:0.1:1,'mt => sc')

subplot(222)
line([0 6],[0 6],[0 6],'color','k','linestyle','--','linewidth',2); hold on
me_in  = arrayfun(@(ii) mean([sc_mt(ii).in.optPerf]),1:15);
me_out = arrayfun(@(ii) mean([sc_mt(ii).out.optPerf]),1:15);
st_in  = repmat(arrayfun(@(ii) std([sc_mt(ii).in.optPerf])./sqrt(nRands),1:15),2,1) .* CI95;
st_out = repmat(arrayfun(@(ii) std([sc_mt(ii).out.optPerf])./sqrt(nRands),1:15),2,1) .* CI95;
arrayfun(@(ii) line([me_in(ii)+st_in(1,ii) me_in(ii)+st_in(2,ii)],[me_out(ii) me_out(ii)],'color',cols(4,:),'linewidth',2),1:15);
arrayfun(@(ii) line([me_in(ii) me_in(ii)],[me_out(ii)+st_out(1,ii) me_out(ii)+st_out(2,ii)],'color',cols(4,:),'linewidth',2),1:15);
plot(me_in,me_out,'o','markersize',8,'color',cols(4,:),'MarkerFaceColor','w','LineWidth',2)
plot(mean(me_in),mean(me_out),'k+','markersize',15,'LineWidth',2)
fixPlot(gca,[0 0.3],[0 0.3],'attend in','attend out',0:0.1:1,0:0.1:1,'sc => mt')

subplot(223)
line([0 6],[0 6],[0 6],'color','k','linestyle','--','linewidth',2); hold on
me_in  = arrayfun(@(ii) mean([mt_mt(ii).in.optPerf]),1:15);
me_out = arrayfun(@(ii) mean([mt_mt(ii).out.optPerf]),1:15);
st_in  = repmat(arrayfun(@(ii) std([mt_mt(ii).in.optPerf])./sqrt(nRands),1:15),2,1) .* CI95;
st_out = repmat(arrayfun(@(ii) std([mt_mt(ii).out.optPerf])./sqrt(nRands),1:15),2,1) .* CI95;
arrayfun(@(ii) line([me_in(ii)+st_in(1,ii) me_in(ii)+st_in(2,ii)],[me_out(ii) me_out(ii)],'color',cols(2,:),'linewidth',2),1:15);
arrayfun(@(ii) line([me_in(ii) me_in(ii)],[me_out(ii)+st_out(1,ii) me_out(ii)+st_out(2,ii)],'color',cols(2,:),'linewidth',2),1:15);
plot(me_in,me_out,'o','markersize',8,'color',cols(2,:),'MarkerFaceColor','w','LineWidth',2)
plot(mean(me_in),mean(me_out),'k+','markersize',15,'LineWidth',2)
fixPlot(gca,[0 1],[0 1],'attend in','attend out',0:0.25:1,0:0.25:1,'mt => mt')

subplot(224)
line([0 6],[0 6],[0 6],'color','k','linestyle','--','linewidth',2); hold on
me_in  = arrayfun(@(ii) mean([sc_sc(ii).in.optPerf]),1:15);
me_out = arrayfun(@(ii) mean([sc_sc(ii).out.optPerf]),1:15);
st_in  = repmat(arrayfun(@(ii) std([sc_sc(ii).in.optPerf])./sqrt(nRands),1:15),2,1) .* CI95;
st_out = repmat(arrayfun(@(ii) std([sc_sc(ii).out.optPerf])./sqrt(nRands),1:15),2,1) .* CI95;
arrayfun(@(ii) line([me_in(ii)+st_in(1,ii) me_in(ii)+st_in(2,ii)],[me_out(ii) me_out(ii)],'color',cols(3,:),'linewidth',2),1:15);
arrayfun(@(ii) line([me_in(ii) me_in(ii)],[me_out(ii)+st_out(1,ii) me_out(ii)+st_out(2,ii)],'color',cols(3,:),'linewidth',2),1:15);
plot(me_in,me_out,'o','markersize',8,'color',cols(3,:),'MarkerFaceColor','w','LineWidth',2)
plot(mean(me_in),mean(me_out),'k+','markersize',15,'LineWidth',2)
fixPlot(gca,[0 1],[0 1],'attend in','attend out',0:0.25:1,0:0.25:1,'sc => sc')

%% supp fig 4b - rsc vs pred ratio 
clf; set(gcf,'color','w','pos',[53,13,1260,784]);
subplot(351)
plot([av.mt_sc_in_perf]./[av.mt_sc_out_perf],[av.mt_in_rsc]-[av.mt_out_rsc],'.','color',cols(1,:),'markersize',15); hold on
plot(mean([av.mt_sc_in_perf]./[av.mt_sc_out_perf]),mean([av.mt_in_rsc]-[av.mt_out_rsc]),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 3],[-0.1 0.2],'pred accuracy ratio','mt rsc difference',0:4,0:0.1:0.2,'mt => sc vs mt')
subplot(356)
plot([av.mt_sc_in_perf]./[av.mt_sc_out_perf],[av.sc_in_rsc]-[av.sc_out_rsc],'.','color',cols(1,:),'markersize',15); hold on
plot(mean([av.mt_sc_in_perf]./[av.mt_sc_out_perf]),mean([av.sc_in_rsc]-[av.sc_out_rsc]),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 3],[-0.1 0.2],'pred accuracy ratio','sc rsc difference',0:4,0:0.1:0.2,'mt => sc vs sc')

subplot(352)
plot([av.mt_mt_in_perf]./[av.mt_mt_out_perf],[av.mt_in_rsc]-[av.mt_out_rsc],'.','color',cols(2,:),'markersize',15); hold on
plot(mean([av.mt_mt_in_perf]./[av.mt_mt_out_perf]),mean([av.mt_in_rsc]-[av.mt_out_rsc]),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 3],[-0.1 0.2],'pred accuracy ratio','mt rsc difference',0:4,0:0.1:0.2,'mt => mt vs mt')
subplot(357)
plot([av.mt_mt_in_perf]./[av.mt_mt_out_perf],[av.sc_in_rsc]-[av.sc_out_rsc],'.','color',cols(2,:),'markersize',15); hold on
plot(mean([av.mt_mt_in_perf]./[av.mt_mt_out_perf]),mean([av.sc_in_rsc]-[av.sc_out_rsc]),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 3],[-0.1 0.2],'pred accuracy ratio','sc rsc difference',0:4,0:0.1:0.2,'mt => mt vs sc')

subplot(353)
plot([av.sc_sc_in_perf]./[av.sc_sc_out_perf],[av.mt_in_rsc]-[av.mt_out_rsc],'.','color',cols(3,:),'markersize',15); hold on
plot(mean([av.sc_sc_in_perf]./[av.sc_sc_out_perf]),mean([av.mt_in_rsc]-[av.mt_out_rsc]),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 3],[-0.1 0.2],'pred accuracy ratio','mt rsc difference',0:4,0:0.1:0.2,'sc => sc vs mt')
subplot(358)
plot([av.sc_sc_in_perf]./[av.sc_sc_out_perf],[av.sc_in_rsc]-[av.sc_out_rsc],'.','color',cols(3,:),'markersize',15); hold on
plot(mean([av.sc_sc_in_perf]./[av.sc_sc_out_perf]),mean([av.sc_in_rsc]-[av.sc_out_rsc]),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 3],[-0.1 0.2],'pred accuracy ratio','sc rsc difference',0:4,0:0.1:0.2,'sc => sc vs sc')

subplot(354)
plot([av.sc_mt_in_perf]./[av.sc_mt_out_perf],[av.mt_in_rsc]-[av.mt_out_rsc],'.','color',cols(4,:),'markersize',15); hold on
plot(mean([av.sc_mt_in_perf]./[av.sc_mt_out_perf]),mean([av.mt_in_rsc]-[av.mt_out_rsc]),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 3],[-0.1 0.2],'pred accuracy ratio','mt rsc difference',0:4,0:0.1:0.2,'sc => mt vs mt')
subplot(359)
plot([av.sc_mt_in_perf]./[av.sc_mt_out_perf],[av.sc_in_rsc]-[av.sc_out_rsc],'.','color',cols(4,:),'markersize',15); hold on
plot(mean([av.sc_mt_in_perf]./[av.sc_mt_out_perf]),mean([av.sc_in_rsc]-[av.sc_out_rsc]),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 3],[-0.1 0.2],'pred accuracy ratio','sc rsc difference',0:4,0:0.1:0.2,'sc => mt vs sc')

subplot(355); hold on;
histogram([av.mt_in_rsc]-[av.mt_out_rsc],linspace(-0.1,0.2,13),'DisplayStyle','stairs','EdgeColor','r','LineWidth',2)
line([mean([av.mt_in_rsc]-[av.mt_out_rsc]) mean([av.mt_in_rsc]-[av.mt_out_rsc])],[0 6],'linestyle','--','color','r','linewidth',2)
fixPlot(gca,[-0.1 0.2],[0 6],'mt rsc','count',0:0.1:0.2,0:2:8,'mt rsc difference')

subplot(3,5,10); hold on;
histogram([av.sc_in_rsc]-[av.sc_out_rsc],linspace(-0.1,0.2,13),'DisplayStyle','stairs','EdgeColor','b','LineWidth',2)
line([mean([av.sc_in_rsc]-[av.sc_out_rsc]) mean([av.sc_in_rsc]-[av.sc_out_rsc])],[0 6],'linestyle','--','color','b','linewidth',2)
fixPlot(gca,[-0.1 0.2],[0 6],'sc rsc','count',0:0.1:0.2,0:2:8,'sc rsc difference')

subplot(3,5,11);
plot([av.mt_sc_in_perf]./[av.mt_sc_out_perf],[av.in_rscSh]-[av.out_rscSh],'.','color',cols(1,:),'markersize',15); hold on
plot(mean([av.mt_sc_in_perf]./[av.mt_sc_out_perf]),mean([av.in_rscSh]-[av.out_rscSh]),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 3],[-0.1 0.2],'pred accuracy ratio','shared rsc difference',0:4,0:0.1:0.2,'mt => sc vs mt')

subplot(3,5,12);
plot([av.mt_mt_in_perf]./[av.mt_mt_out_perf],[av.in_rscSh]-[av.out_rscSh],'.','color',cols(2,:),'markersize',15); hold on
plot(mean([av.mt_mt_in_perf]./[av.mt_mt_out_perf]),mean([av.in_rscSh]-[av.out_rscSh]),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 3],[-0.1 0.2],'pred accuracy ratio','shared rsc difference',0:4,0:0.1:0.2,'mt => mt vs mt')

subplot(3,5,13);
plot([av.sc_sc_in_perf]./[av.sc_sc_out_perf],[av.in_rscSh]-[av.out_rscSh],'.','color',cols(3,:),'markersize',15); hold on
plot(mean([av.sc_sc_in_perf]./[av.sc_sc_out_perf]),mean([av.in_rscSh]-[av.out_rscSh]),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 3],[-0.1 0.2],'pred accuracy ratio','shared rsc difference',0:4,0:0.1:0.2,'sc => sc vs mt')

subplot(3,5,14);
plot([av.sc_mt_in_perf]./[av.sc_mt_out_perf],[av.in_rscSh]-[av.out_rscSh],'.','color',cols(4,:),'markersize',15); hold on
plot(mean([av.sc_mt_in_perf]./[av.sc_mt_out_perf]),mean([av.in_rscSh]-[av.out_rscSh]),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 3],[-0.1 0.2],'pred accuracy ratio','shared rsc difference',0:4,0:0.1:0.2,'sc => mt vs mt')

subplot(3,5,15); hold on;
histogram([av.in_rscSh]-[av.out_rscSh],linspace(-0.1,0.2,13),'DisplayStyle','stairs','EdgeColor','k','LineWidth',2)
line([mean([av.in_rscSh]-[av.out_rscSh]) mean([av.in_rscSh]-[av.out_rscSh])],[0 6],'linestyle','--','color','k','linewidth',2)
fixPlot(gca,[-0.1 0.2],[0 6],'mt rsc','count',0:0.1:0.2,0:2:8,'shared rsc difference')

%% supp fig 4b - attn index 
clf; set(gcf,'color','w','pos',[53,114,1228,591]);
subplot(251)

mt_attnidx = ([av.mt_in_rate]-[av.mt_out_rate])./([av.mt_in_rate]+[av.mt_out_rate]);
sc_attnidx = ([av.sc_in_rate]-[av.sc_out_rate])./([av.sc_in_rate]+[av.sc_out_rate]);

plot([av.mt_sc_in_perf]./[av.mt_sc_out_perf],mt_attnidx,'.','color',cols(1,:),'markersize',15); hold on
plot(mean([av.mt_sc_in_perf]./[av.mt_sc_out_perf]),mean(mt_attnidx),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 3],[-0.15 0.15],'pred accuracy ratio','mt attn idx',0:4,[-0.1 0 0.1],'mt => sc vs mt')
subplot(256)
plot([av.mt_sc_in_perf]./[av.mt_sc_out_perf],sc_attnidx,'.','color',cols(1,:),'markersize',15); hold on
plot(mean([av.mt_sc_in_perf]./[av.mt_sc_out_perf]),mean(sc_attnidx),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 3],[-0.15 0.15],'pred accuracy ratio','sc attn idx',0:4,[-0.1 0 0.1],'mt => sc vs sc')

subplot(252)
plot([av.mt_mt_in_perf]./[av.mt_mt_out_perf],mt_attnidx,'.','color',cols(2,:),'markersize',15); hold on
plot(mean([av.mt_mt_in_perf]./[av.mt_mt_out_perf]),mean(mt_attnidx),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 3],[-0.15 0.15],'pred accuracy ratio','mt attn idx',0:4,[-0.1 0 0.1],'mt => mt vs mt')
subplot(257)
plot([av.mt_mt_in_perf]./[av.mt_mt_out_perf],sc_attnidx,'.','color',cols(2,:),'markersize',15); hold on
plot(mean([av.mt_mt_in_perf]./[av.mt_mt_out_perf]),mean(sc_attnidx),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 3],[-0.15 0.15],'pred accuracy ratio','sc attn idx',0:4,[-0.1 0 0.1],'mt => mt vs sc')

subplot(253)
plot([av.sc_sc_in_perf]./[av.sc_sc_out_perf],mt_attnidx,'.','color',cols(3,:),'markersize',15); hold on
plot(mean([av.sc_sc_in_perf]./[av.sc_sc_out_perf]),mean(mt_attnidx),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 3],[-0.15 0.15],'pred accuracy ratio','mt attn idx',0:4,[-0.1 0 0.1],'sc => sc vs mt')
subplot(258)
plot([av.sc_sc_in_perf]./[av.sc_sc_out_perf],sc_attnidx,'.','color',cols(3,:),'markersize',15); hold on
plot(mean([av.sc_sc_in_perf]./[av.sc_sc_out_perf]),mean(sc_attnidx),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 3],[-0.15 0.15],'pred accuracy ratio','sc attn idx',0:4,[-0.1 0 0.1],'sc => sc vs sc')

subplot(254)
plot([av.sc_mt_in_perf]./[av.sc_mt_out_perf],mt_attnidx,'.','color',cols(4,:),'markersize',15); hold on
plot(mean([av.sc_mt_in_perf]./[av.sc_mt_out_perf]),mean(mt_attnidx),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 3],[-0.15 0.15],'pred accuracy ratio','mt attn idx',0:4,[-0.1 0 0.1],'sc => mt vs mt')
subplot(259)
plot([av.sc_mt_in_perf]./[av.sc_mt_out_perf],sc_attnidx,'.','color',cols(4,:),'markersize',15); hold on
plot(mean([av.sc_mt_in_perf]./[av.sc_mt_out_perf]),mean(sc_attnidx),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 3],[-0.15 0.15],'pred accuracy ratio','sc attn idx',0:4,[-0.1 0 0.1],'sc => mt vs sc')

subplot(255); hold on;
histogram(mt_attnidx,linspace(-0.15,0.15,9),'DisplayStyle','stairs','EdgeColor','r','LineWidth',2)
line([mean(mt_attnidx) mean(mt_attnidx)],[0 6],'linestyle','--','color','r','linewidth',2)
fixPlot(gca,[-0.15 0.15],[0 6],'mt rate','count',[-0.1 0 0.1],0:2:8,'mt attn idx')

subplot(2,5,10); hold on;
histogram(sc_attnidx,linspace(-0.15,0.15,9),'DisplayStyle','stairs','EdgeColor','b','LineWidth',2)
line([mean(sc_attnidx) mean(sc_attnidx)],[0 6],'linestyle','--','color','b','linewidth',2)
fixPlot(gca,[-0.15 0.15],[0 6],'sc rate','count',[-0.1 0 0.1],0:2:8,'sc attn idx')

%% fig 6 - semedo replication 
% dim vs perfdim
set(gcf,'color','w','pos',[75,68,1102,624]); clf;
subplot(241)
plot([av.mt_in_dim],[av.mt_mt_in_perfdim],'.','markersize',15,'color',[0.2 0.5 1]); hold on
plot([av.sc_in_dim],[av.mt_sc_in_perfdim],'.','markersize',15,'color',[1 0.5 0.2]); hold on
plot(mean([av.mt_in_dim]),mean([av.mt_mt_in_perfdim]),'+','markersize',15,'LineWidth',2,'color',[0.2 0.5 1]); hold on
plot(mean([av.sc_in_dim]),mean([av.mt_sc_in_perfdim]),'+','markersize',15,'LineWidth',2,'color',[1 0.5 0.2]); hold on
line([0 6],[0 6],[0 6],'color','k','linestyle','--','linewidth',2)
fixPlot(gca,[0 6],[0 6],'target dims','pred dims',0:2:6,0:2:6,{'target vs pred dims' 'attend in'},{'mt->mt' 'mt->sc'})

subplot(242)
plot([av.mt_out_dim],[av.mt_mt_out_perfdim],'.','markersize',15,'color',[0.2 0.5 1]); hold on
plot([av.sc_out_dim],[av.mt_sc_out_perfdim],'.','markersize',15,'color',[1 0.5 0.2]); hold on
plot(mean([av.mt_out_dim]),mean([av.mt_mt_out_perfdim]),'+','markersize',15,'LineWidth',2,'color',[0.2 0.5 1]); hold on
plot(mean([av.sc_out_dim]),mean([av.mt_sc_out_perfdim]),'+','markersize',15,'LineWidth',2,'color',[1 0.5 0.2]); hold on
line([0 6],[0 6],[0 6],'color','k','linestyle','--','linewidth',2)
fixPlot(gca,[0 6],[0 6],'target dims','pred dims',0:2:6,0:2:6,{'target vs pred dims' 'attend out'},{'mt->mt' 'mt->sc'})

subplot(245)
plot([av.mt_in_dim]./[av.mt_out_dim],[av.mt_mt_in_perfdim]./[av.mt_mt_out_perfdim],'.','markersize',15,'color',[0.2 0.5 1]); hold on
plot([av.sc_in_dim]./[av.sc_out_dim],[av.mt_sc_in_perfdim]./[av.mt_sc_out_perfdim],'.','markersize',15,'color',[1 0.5 0.2]); hold on
plot(mean([av.mt_in_dim]./[av.mt_out_dim]),mean([av.mt_mt_in_perfdim]./[av.mt_mt_out_perfdim]),'+','markersize',15,'LineWidth',2,'color',[0.2 0.5 1]); hold on
plot(mean([av.sc_in_dim]./[av.sc_out_dim]),mean([av.mt_sc_in_perfdim]./[av.mt_sc_out_perfdim]),'+','markersize',15,'LineWidth',2,'color',[1 0.5 0.2]); hold on
fixPlot(gca,[0 2],[0 2],{'target dims' 'in/out'},{'pred dims' 'in/out'},0:2,0:2,'target vs pred dims',{'mt->mt' 'mt->sc'})

subplot(243)
plot([av.sc_in_dim],[av.sc_sc_in_perfdim],'.','markersize',15,'color',[0.2 0.7 0.3]); hold on
plot([av.mt_in_dim],[av.sc_mt_in_perfdim],'.','markersize',15,'color',[0.8 0.3 0.8]); hold on
plot(mean([av.sc_in_dim]),mean([av.sc_sc_in_perfdim]),'+','markersize',15,'LineWidth',2,'color',[0.2 0.7 0.3]); hold on
plot(mean([av.mt_in_dim]),mean([av.sc_mt_in_perfdim]),'+','markersize',15,'LineWidth',2,'color',[0.8 0.3 0.8]); hold on
line([0 6],[0 6],[0 6],'color','k','linestyle','--','linewidth',2)
fixPlot(gca,[0 6],[0 6],'target dims','pred dims',0:2:6,0:2:6,{'target vs pred dims' 'attend in'},{'sc->sc' 'sc->mt'})

subplot(244)
plot([av.sc_out_dim],[av.sc_sc_out_perfdim],'.','markersize',15,'color',[0.2 0.7 0.3]); hold on
plot([av.mt_out_dim],[av.sc_mt_out_perfdim],'.','markersize',15,'color',[0.8 0.3 0.8]); hold on
plot(mean([av.sc_out_dim]),mean([av.sc_sc_out_perfdim]),'+','markersize',15,'LineWidth',2,'color',[0.2 0.7 0.3]); hold on
plot(mean([av.mt_out_dim]),mean([av.sc_mt_out_perfdim]),'+','markersize',15,'LineWidth',2,'color',[0.8 0.3 0.8]); hold on
line([0 6],[0 6],[0 6],'color','k','linestyle','--','linewidth',2)
fixPlot(gca,[0 6],[0 6],'target dims','pred dims',0:2:6,0:2:6,{'target vs pred dims' 'attend out'},{'sc->sc' 'sc->mt'})

subplot(247)
plot([av.sc_in_dim]./[av.sc_out_dim],[av.sc_sc_in_perfdim]./[av.sc_sc_out_perfdim],'.','markersize',15,'color',[0.2 0.7 0.3]); hold on
plot([av.mt_in_dim]./[av.mt_out_dim],[av.sc_mt_in_perfdim]./[av.sc_mt_out_perfdim],'.','markersize',15,'color',[0.8 0.3 0.8]); hold on
plot(mean([av.mt_in_dim]./[av.mt_out_dim]),mean([av.mt_mt_in_perfdim]./[av.mt_mt_out_perfdim]),'+','markersize',15,'LineWidth',2,'color',[0.2 0.7 0.3]); hold on
plot(mean([av.sc_in_dim]./[av.sc_out_dim]),mean([av.mt_sc_in_perfdim]./[av.mt_sc_out_perfdim]),'+','markersize',15,'LineWidth',2,'color',[0.8 0.3 0.8]); hold on
fixPlot(gca,[0 2],[0 2],{'target dims' 'in/out'},{'pred dims' 'in/out'},0:2,0:2,'target vs pred dims',{'sc->sc' 'sc->mt'})

%% supp fig 6a - out dim and pred dim 
clf; set(gcf,'color','w','pos',[75,68,1102,624]);
subplot(241)
plot([av.mt_sc_in_perfdim],[av.mt_sc_out_perfdim],'.','markersize',15,'color',cols(1,:)); hold on
plot(mean([av.mt_sc_in_perfdim]),mean([av.mt_sc_out_perfdim]),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 3],[0 3],'perf dim in','perf dim out',0:3,0:3,'mt => sc')

subplot(242)
plot([av.mt_mt_in_perfdim],[av.mt_mt_out_perfdim],'.','markersize',15,'color',cols(2,:)); hold on
plot(mean([av.mt_mt_in_perfdim]),mean([av.mt_mt_out_perfdim]),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 6],[0 6],'perf dim in','perf dim out',0:2:6,0:2:6,'mt => mt')

subplot(243)
plot([av.sc_sc_in_perfdim],[av.sc_sc_out_perfdim],'.','markersize',15,'color',cols(3,:)); hold on
plot(mean([av.sc_sc_in_perfdim]),mean([av.sc_sc_out_perfdim]),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 6],[0 6],'perf dim in','perf dim out',0:2:6,0:2:6,'sc => sc')

subplot(244)
plot([av.sc_mt_in_perfdim],[av.sc_mt_out_perfdim],'.','markersize',15,'color',cols(4,:)); hold on
plot(mean([av.sc_mt_in_perfdim]),mean([av.sc_mt_out_perfdim]),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 3],[0 3],'perf dim in','perf dim out',0:3,0:3,'sc => mt')

subplot(245)
plot([av.sc_in_dim],[av.sc_out_dim],'.','markersize',15,'color',[0 0.3 0.8]); hold on
plot(mean([av.sc_in_dim]),mean([av.sc_out_dim]),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 6],[0 6],'dim in','dim out',0:2:6,0:2:6,'sc')

subplot(246)
plot([av.mt_in_dim],[av.mt_out_dim],'.','markersize',15,'color',[0.8 0.3 0]); hold on
plot(mean([av.mt_in_dim]),mean([av.mt_out_dim]),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 6],[0 6],'dim in','dim out',0:2:6,0:2:6,'mt')

%% supp fig 6b - per condition 
% dim vs dim
% perfdim vs perfdim
set(gcf,'color','w','pos',[75,68,1102,624]); clf;

subplot(231)
plot([av.mt_in_dim],[av.sc_in_dim],'.','markersize',15,'color',[0 0 0]); hold on
plot([av.mt_out_dim],[av.sc_out_dim],'.','markersize',15,'color',[0.7 0.7 0.7]); hold on
plot(mean([av.mt_in_dim]),mean([av.sc_in_dim]),'k+','markersize',15,'LineWidth',2); hold on
plot(mean([av.mt_out_dim]),mean([av.sc_out_dim]),'+','markersize',15,'LineWidth',2,'color',[0.7 0.7 0.7]); hold on
line([0 6],[0 6],[0 6],'color','k','linestyle','--','linewidth',2)
fixPlot(gca,[0 6],[0 6],'mt dim','sc dim',0:2:6,0:2:6,'mt vs sc dims')

subplot(234)
plot([av.mt_in_dim]./[av.mt_out_dim],[av.sc_in_dim]./[av.sc_out_dim],'.','markersize',15,'color',[0 0 0]); hold on
plot(mean([av.mt_in_dim]./[av.mt_out_dim]),mean([av.sc_in_dim]./[av.sc_out_dim]),'k+','markersize',15,'LineWidth',2); hold on
fixPlot(gca,[0 2],[0 2],{'mt dim' 'in/out'},{'sc dim' 'in/out'},0:2,0:2,'mt vs sc dims')

subplot(232)
plot([av.mt_mt_in_perfdim],[av.mt_sc_in_perfdim],'.','markersize',15,'color',[0 0 0]); hold on
plot([av.mt_mt_out_perfdim],[av.mt_sc_out_perfdim],'.','markersize',15,'color',[0.7 0.7 0.7]); hold on
plot(mean([av.mt_mt_in_perfdim]),mean([av.mt_sc_in_perfdim]),'k+','markersize',15,'LineWidth',2); hold on
plot(mean([av.mt_mt_out_perfdim]),mean([av.mt_sc_out_perfdim]),'+','markersize',15,'LineWidth',2,'color',[0.7 0.7 0.7]); hold on
line([0 6],[0 6],[0 6],'color','k','linestyle','--','linewidth',2)
fixPlot(gca,[0 6],[0 6],'mt->mt perfdim','mt->sc perfdim',0:2:6,0:2:6,'mt->mt vs mt->sc')

subplot(235)
plot([av.mt_mt_in_perfdim]./[av.mt_mt_out_perfdim],[av.mt_sc_in_perfdim]./[av.mt_sc_out_perfdim],'.','markersize',15,'color',[0 0 0]); hold on
plot(mean([av.mt_mt_in_perfdim]./[av.mt_mt_out_perfdim]),mean([av.mt_sc_in_perfdim]./[av.mt_sc_out_perfdim]),'+','markersize',15,'LineWidth',2,'color',[0 0 0]); hold on
fixPlot(gca,[0 2],[0 2],{'mt->mt perfdim' 'in/out'},{'mt->sc perfdim' 'in/out'},0:2,0:2,'mt->mt vs mt->sc')

subplot(233)
plot([av.sc_sc_in_perfdim],[av.sc_mt_in_perfdim],'.','markersize',15,'color',[0 0 0]); hold on
plot([av.sc_sc_out_perfdim],[av.sc_mt_out_perfdim],'.','markersize',15,'color',[0.7 0.7 0.7]); hold on
plot(mean([av.sc_sc_in_perfdim]),mean([av.sc_mt_in_perfdim]),'k+','markersize',15,'LineWidth',2); hold on
plot(mean([av.sc_sc_out_perfdim]),mean([av.sc_mt_out_perfdim]),'+','markersize',15,'LineWidth',2,'color',[0.7 0.7 0.7]); hold on
line([0 6],[0 6],[0 6],'color','k','linestyle','--','linewidth',2)
fixPlot(gca,[0 6],[0 6],'sc->sc perfdim','sc->mt perfdim',0:2:6,0:2:6,'sc->sc vs sc->mt')

subplot(236)
plot([av.sc_sc_in_perfdim]./[av.sc_sc_out_perfdim],[av.sc_mt_in_perfdim]./[av.sc_mt_out_perfdim],'.','markersize',15,'color',[0 0 0]); hold on
plot(mean([av.sc_sc_in_perfdim]./[av.sc_sc_out_perfdim]),mean([av.sc_mt_in_perfdim]./[av.sc_mt_out_perfdim]),'+','markersize',15,'LineWidth',2,'color',[0 0 0]); hold on
fixPlot(gca,[0 2],[0 2],{'sc->sc perfdim' 'in/out'},{'sc->mt perfdim' 'in/out'},0:2,0:2,'sc->sc vs sc->mt')

%% EXTRA FIGURES:

%% non overlapping histograms from fig 4 
nSds = 2;
me_in  = arrayfun(@(ii) mean([mt_sc(ii).in.optPerf]),1:15);
me_out = arrayfun(@(ii) mean([mt_sc(ii).out.optPerf]),1:15);
st_in  = arrayfun(@(ii) std ([mt_sc(ii).in.optPerf]),1:15) .* nSds;
st_out = arrayfun(@(ii) std ([mt_sc(ii).out.optPerf]),1:15) .* nSds;
gc = sum([me_in-st_in/2;me_in+st_in/2]-repmat(me_out,2,1)>0)~=1 & sum([me_out-st_out/2;me_out+st_out/2]-repmat(me_in,2,1)>0)~=1;
subplot(241); cla;
histogram(me_in./me_out,linspace(0,3,11),'DisplayStyle','stairs','edgecolor',cols(1,:),'linewidth',2); hold on;
histogram(me_in(gc)./me_out(gc),linspace(0,3,11),'edgecolor',cols(1,:),'linewidth',2,'facecolor',cols(1,:),'facealpha',1)
fixPlot(gca,[0 4],[0 12],'pred ratio','count',1:3,0:5:12,'mt => sc')

me_in  = arrayfun(@(ii) mean([mt_mt(ii).in.optPerf]),1:15);
me_out = arrayfun(@(ii) mean([mt_mt(ii).out.optPerf]),1:15);
st_in  = arrayfun(@(ii) std ([mt_mt(ii).in.optPerf]),1:15) .* nSds;
st_out = arrayfun(@(ii) std ([mt_mt(ii).out.optPerf]),1:15) .* nSds;
gc = sum([me_in-st_in/2;me_in+st_in/2]-repmat(me_out,2,1)>0)~=1 & sum([me_out-st_out/2;me_out+st_out/2]-repmat(me_in,2,1)>0)~=1;
subplot(242); cla;
histogram(me_in./me_out,linspace(0,3,11),'DisplayStyle','stairs','edgecolor',cols(2,:),'linewidth',2); hold on;
histogram(me_in(gc)./me_out(gc),linspace(0,3,11),'edgecolor',cols(2,:),'linewidth',2,'facecolor',cols(2,:),'facealpha',1)
fixPlot(gca,[0 4],[0 12],'pred ratio','count',1:3,0:5:12,'mt => mt')

me_in  = arrayfun(@(ii) mean([sc_mt(ii).in.optPerf]),1:15);
me_out = arrayfun(@(ii) mean([sc_mt(ii).out.optPerf]),1:15);
st_in  = arrayfun(@(ii) std ([sc_mt(ii).in.optPerf]),1:15) .* nSds;
st_out = arrayfun(@(ii) std ([sc_mt(ii).out.optPerf]),1:15) .* nSds;
gc = sum([me_in-st_in/2;me_in+st_in/2]-repmat(me_out,2,1)>0)~=1 & sum([me_out-st_out/2;me_out+st_out/2]-repmat(me_in,2,1)>0)~=1;
subplot(243); cla;
histogram(me_in./me_out,linspace(0,3,11),'DisplayStyle','stairs','edgecolor',cols(4,:),'linewidth',2); hold on;
histogram(me_in(gc)./me_out(gc),linspace(0,3,11),'edgecolor',cols(4,:),'linewidth',2,'facecolor',cols(4,:),'facealpha',1)
fixPlot(gca,[0 4],[0 12],'pred ratio','count',1:3,0:5:12,'sc => mt')

me_in  = arrayfun(@(ii) mean([sc_sc(ii).in.optPerf]),1:15);
me_out = arrayfun(@(ii) mean([sc_sc(ii).out.optPerf]),1:15);
st_in  = arrayfun(@(ii) std ([sc_sc(ii).in.optPerf]),1:15) .* nSds;
st_out = arrayfun(@(ii) std ([sc_sc(ii).out.optPerf]),1:15) .* nSds;
gc = sum([me_in-st_in/2;me_in+st_in/2]-repmat(me_out,2,1)>0)~=1 & sum([me_out-st_out/2;me_out+st_out/2]-repmat(me_in,2,1)>0)~=1;
subplot(244); cla;
histogram(me_in./me_out,linspace(0,3,11),'DisplayStyle','stairs','edgecolor',cols(3,:),'linewidth',2); hold on;
histogram(me_in(gc)./me_out(gc),linspace(0,3,11),'edgecolor',cols(3,:),'linewidth',2,'facecolor',cols(3,:),'facealpha',1)
fixPlot(gca,[0 4],[0 12],'pred ratio','count',1:3,0:5:12,'sc => sc')

me_in  = arrayfun(@(ii) mean([mt_sc(ii).in.optDim]),1:15);
me_out = arrayfun(@(ii) mean([mt_sc(ii).out.optDim]),1:15);
st_in  = arrayfun(@(ii) std ([mt_sc(ii).in.optDim]),1:15) .* nSds;
st_out = arrayfun(@(ii) std ([mt_sc(ii).out.optDim]),1:15) .* nSds;
gc = sum([me_in-st_in/2;me_in+st_in/2]-repmat(me_out,2,1)>0)~=1 & sum([me_out-st_out/2;me_out+st_out/2]-repmat(me_in,2,1)>0)~=1;
subplot(245); cla;
histogram(me_in./me_out,linspace(0,3,11),'DisplayStyle','stairs','edgecolor',cols(1,:),'linewidth',2); hold on;
histogram(me_in(gc)./me_out(gc),linspace(0,3,11),'edgecolor',cols(1,:),'linewidth',2,'facecolor',cols(1,:),'facealpha',1)
fixPlot(gca,[0 4],[0 12],'pred ratio','count',1:3,0:5:12,'mt => sc')

me_in  = arrayfun(@(ii) mean([mt_mt(ii).in.optDim]),1:15);
me_out = arrayfun(@(ii) mean([mt_mt(ii).out.optDim]),1:15);
st_in  = arrayfun(@(ii) std ([mt_mt(ii).in.optDim]),1:15) .* nSds;
st_out = arrayfun(@(ii) std ([mt_mt(ii).out.optDim]),1:15) .* nSds;
gc = sum([me_in-st_in/2;me_in+st_in/2]-repmat(me_out,2,1)>0)~=1 & sum([me_out-st_out/2;me_out+st_out/2]-repmat(me_in,2,1)>0)~=1;
subplot(246); cla;
histogram(me_in./me_out,linspace(0,3,11),'DisplayStyle','stairs','edgecolor',cols(2,:),'linewidth',2); hold on;
histogram(me_in(gc)./me_out(gc),linspace(0,3,11),'edgecolor',cols(2,:),'linewidth',2,'facecolor',cols(2,:),'facealpha',1)
fixPlot(gca,[0 4],[0 12],'pred ratio','count',1:3,0:5:12,'mt => mt')

me_in  = arrayfun(@(ii) mean([sc_mt(ii).in.optDim]),1:15);
me_out = arrayfun(@(ii) mean([sc_mt(ii).out.optDim]),1:15);
st_in  = arrayfun(@(ii) std ([sc_mt(ii).in.optDim]),1:15) .* nSds;
st_out = arrayfun(@(ii) std ([sc_mt(ii).out.optDim]),1:15) .* nSds;
gc = sum([me_in-st_in/2;me_in+st_in/2]-repmat(me_out,2,1)>0)~=1 & sum([me_out-st_out/2;me_out+st_out/2]-repmat(me_in,2,1)>0)~=1;
subplot(247); cla;
histogram(me_in./me_out,linspace(0,3,11),'DisplayStyle','stairs','edgecolor',cols(4,:),'linewidth',2); hold on;
histogram(me_in(gc)./me_out(gc),linspace(0,3,11),'edgecolor',cols(4,:),'linewidth',2,'facecolor',cols(4,:),'facealpha',1)
fixPlot(gca,[0 4],[0 12],'pred ratio','count',1:3,0:5:12,'sc => mt')

me_in  = arrayfun(@(ii) mean([sc_sc(ii).in.optDim]),1:15);
me_out = arrayfun(@(ii) mean([sc_sc(ii).out.optDim]),1:15);
st_in  = arrayfun(@(ii) std ([sc_sc(ii).in.optDim]),1:15) .* nSds;
st_out = arrayfun(@(ii) std ([sc_sc(ii).out.optDim]),1:15) .* nSds;
gc = sum([me_in-st_in/2;me_in+st_in/2]-repmat(me_out,2,1)>0)~=1 & sum([me_out-st_out/2;me_out+st_out/2]-repmat(me_in,2,1)>0)~=1;
subplot(248); cla;
histogram(me_in./me_out,linspace(0,3,11),'DisplayStyle','stairs','edgecolor',cols(3,:),'linewidth',2); hold on;
histogram(me_in(gc)./me_out(gc),linspace(0,3,11),'edgecolor',cols(3,:),'linewidth',2,'facecolor',cols(3,:),'facealpha',1)
fixPlot(gca,[0 4],[0 12],'pred ratio','count',1:3,0:5:12,'sc => sc')

%% pred dims and pred accuracy as means of ratios 
clf; set(gcf,'color','w','pos',[53,1,857,704]);
subplot(221)
plot(ratio.perf_mt_sc_perSes,ratio.dim_mt_sc_perSes,'.','color',cols(1,:),'markersize',15); hold on
plot(mean(ratio.perf_mt_sc_perSes),mean(ratio.dim_mt_sc_perSes),'k+','markersize',15,'LineWidth',2);
fixPlot(gca,[0 3],[0 3],'pred accuracy ratio','pred dims ratio',0:4,0:3,'mt => sc')

subplot(222)
plot(ratio.perf_sc_mt_perSes,ratio.dim_sc_mt_perSes,'.','color',cols(4,:),'markersize',15); hold on
plot(mean(ratio.perf_sc_mt_perSes),mean(ratio.dim_sc_mt_perSes),'k+','markersize',15,'LineWidth',2);
fixPlot(gca,[0 3],[0 3],'pred accuracy ratio','pred dims ratio',0:4,0:3,'sc => mt')

subplot(223)
plot(ratio.perf_mt_mt_perSes,ratio.dim_mt_mt_perSes,'.','color',cols(2,:),'markersize',15); hold on
plot(mean(ratio.perf_mt_mt_perSes),mean(ratio.dim_mt_mt_perSes),'k+','markersize',15,'LineWidth',2);
fixPlot(gca,[0 3],[0 3],'pred accuracy ratio','pred dims ratio',0:4,0:3,'mt => mt')

subplot(224)
plot(ratio.perf_sc_sc_perSes,ratio.dim_sc_sc_perSes,'.','color',cols(3,:),'markersize',15); hold on
plot(mean(ratio.perf_sc_sc_perSes),mean(ratio.dim_sc_sc_perSes),'k+','markersize',15,'LineWidth',2);
fixPlot(gca,[0 3],[0 3],'pred accuracy ratio','pred dims ratio',0:4,0:3,'sc => sc')

%% histogram of attention mod index for pred acc 
clf; set(gcf,'color','w','pos',[53,1,857,704]);
arrayfun(@(ii) line(subplot(2,2,ii),[0 0],[0 12],'color','k','linewidth',2,'linestyle','--'),1:4)
arrayfun(@(ii) hold(subplot(2,2,ii),'on'),1:4)
histogram(mean(stat.perfIdx_mt_sc_perSes,2),-0.6:0.1:0.6,'edgecolor',cols(1,:),'DisplayStyle','stairs','linewidth',2,'parent',subplot(221));
histogram(mean(stat.perfIdx_mt_mt_perSes,2),-0.6:0.1:0.6,'edgecolor',cols(2,:),'DisplayStyle','stairs','linewidth',2,'parent',subplot(222));
histogram(mean(stat.perfIdx_sc_sc_perSes,2),-0.6:0.1:0.6,'edgecolor',cols(3,:),'DisplayStyle','stairs','linewidth',2,'parent',subplot(223));
histogram(mean(stat.perfIdx_sc_mt_perSes,2),-0.6:0.1:0.6,'edgecolor',cols(4,:),'DisplayStyle','stairs','linewidth',2,'parent',subplot(224));
line(subplot(221),[median(mean(stat.perfIdx_mt_sc_perSes,2)) median(mean(stat.perfIdx_mt_sc_perSes,2))],[0 12],'color',cols(1,:),'linewidth',2,'linestyle','--');
line(subplot(222),[median(mean(stat.perfIdx_mt_mt_perSes,2)) median(mean(stat.perfIdx_mt_mt_perSes,2))],[0 12],'color',cols(2,:),'linewidth',2,'linestyle','--');
line(subplot(223),[median(mean(stat.perfIdx_sc_sc_perSes,2)) median(mean(stat.perfIdx_sc_sc_perSes,2))],[0 12],'color',cols(3,:),'linewidth',2,'linestyle','--');
line(subplot(224),[median(mean(stat.perfIdx_sc_mt_perSes,2)) median(mean(stat.perfIdx_sc_mt_perSes,2))],[0 12],'color',cols(4,:),'linewidth',2,'linestyle','--');
arrayfun(@(ii) fixPlot(subplot(2,2,ii),[-0.6 0.6],[0 12],'attn mod index pred acc','count',-0.6:0.2:0.8,0:2:12),1:4)

%% pred dims 
clf; set(gcf,'color','w','pos',[53,1,857,704]);
nRands = 1000;
CI95 = repmat(tinv([0.025 0.975], nRands-1)',1,15);
subplot(221)
line([0 6],[0 6],[0 6],'color','k','linestyle','--','linewidth',2); hold on
me_in  = arrayfun(@(ii) mean([mt_sc(ii).in.optDim]),1:15);
me_out = arrayfun(@(ii) mean([mt_sc(ii).out.optDim]),1:15);
st_in  = repmat(arrayfun(@(ii) std([mt_sc(ii).in.optDim])./sqrt(nRands),1:15),2,1) .* CI95;
st_out = repmat(arrayfun(@(ii) std([mt_sc(ii).out.optDim])./sqrt(nRands),1:15),2,1) .* CI95;
arrayfun(@(ii) line([me_in(ii)+st_in(1,ii) me_in(ii)+st_in(2,ii)],[me_out(ii) me_out(ii)],'color',cols(1,:),'linewidth',2),1:15);
arrayfun(@(ii) line([me_in(ii) me_in(ii)],[me_out(ii)+st_out(1,ii) me_out(ii)+st_out(2,ii)],'color',cols(1,:),'linewidth',2),1:15);
plot(me_in,me_out,'o','markersize',8,'color',cols(1,:),'MarkerFaceColor','w','LineWidth',2)
plot(mean(me_in),mean(me_out),'k+','markersize',15,'LineWidth',2)
fixPlot(gca,[0 4],[0 4],'attend in','attend out',0:1:10,0:1:10,'mt => sc')

subplot(222)
line([0 6],[0 6],[0 6],'color','k','linestyle','--','linewidth',2); hold on
me_in  = arrayfun(@(ii) mean([sc_mt(ii).in.optDim]),1:15);
me_out = arrayfun(@(ii) mean([sc_mt(ii).out.optDim]),1:15);
st_in  = repmat(arrayfun(@(ii) std([sc_mt(ii).in.optDim])./sqrt(nRands),1:15),2,1) .* CI95;
st_out = repmat(arrayfun(@(ii) std([sc_mt(ii).out.optDim])./sqrt(nRands),1:15),2,1) .* CI95;
arrayfun(@(ii) line([me_in(ii)+st_in(1,ii) me_in(ii)+st_in(2,ii)],[me_out(ii) me_out(ii)],'color',cols(4,:),'linewidth',2),1:15);
arrayfun(@(ii) line([me_in(ii) me_in(ii)],[me_out(ii)+st_out(1,ii) me_out(ii)+st_out(2,ii)],'color',cols(4,:),'linewidth',2),1:15);
plot(me_in,me_out,'ko','markersize',8,'color',cols(4,:),'MarkerFaceColor','w','LineWidth',2)
plot(mean(me_in),mean(me_out),'k+','markersize',15,'LineWidth',2)
fixPlot(gca,[0 4],[0 4],'attend in','attend out',0:1:10,0:1:10,'sc => mt')

subplot(223)
line([0 6],[0 6],[0 6],'color','k','linestyle','--','linewidth',2); hold on
me_in  = arrayfun(@(ii) mean([mt_mt(ii).in.optDim]),1:15);
me_out = arrayfun(@(ii) mean([mt_mt(ii).out.optDim]),1:15);
st_in  = repmat(arrayfun(@(ii) std([mt_mt(ii).in.optDim])./sqrt(nRands),1:15),2,1) .* CI95;
st_out = repmat(arrayfun(@(ii) std([mt_mt(ii).out.optDim])./sqrt(nRands),1:15),2,1) .* CI95;
arrayfun(@(ii) line([me_in(ii)+st_in(1,ii) me_in(ii)+st_in(2,ii)],[me_out(ii) me_out(ii)],'color',cols(2,:),'linewidth',2),1:15);
arrayfun(@(ii) line([me_in(ii) me_in(ii)],[me_out(ii)+st_out(1,ii) me_out(ii)+st_out(2,ii)],'color',cols(2,:),'linewidth',2),1:15);
plot(me_in,me_out,'ko','markersize',8,'color',cols(2,:),'MarkerFaceColor','w','LineWidth',2)
plot(mean(me_in),mean(me_out),'k+','markersize',15,'LineWidth',2)
fixPlot(gca,[0 6],[0 6],'attend in','attend out',0:2:6,0:2:6,'mt => mt')

subplot(224)
line([0 6],[0 6],[0 6],'color','k','linestyle','--','linewidth',2); hold on
me_in  = arrayfun(@(ii) mean([sc_sc(ii).in.optDim]),1:15);
me_out = arrayfun(@(ii) mean([sc_sc(ii).out.optDim]),1:15);
st_in  = repmat(arrayfun(@(ii) std([sc_sc(ii).in.optDim])./sqrt(nRands),1:15),2,1) .* CI95;
st_out = repmat(arrayfun(@(ii) std([sc_sc(ii).out.optDim])./sqrt(nRands),1:15),2,1) .* CI95;
arrayfun(@(ii) line([me_in(ii)+st_in(1,ii) me_in(ii)+st_in(2,ii)],[me_out(ii) me_out(ii)],'color',cols(3,:),'linewidth',2),1:15);
arrayfun(@(ii) line([me_in(ii) me_in(ii)],[me_out(ii)+st_out(1,ii) me_out(ii)+st_out(2,ii)],'color',cols(3,:),'linewidth',2),1:15);
plot(me_in,me_out,'o','markersize',8,'color',cols(3,:),'MarkerFaceColor','w','LineWidth',2)
plot(mean(me_in),mean(me_out),'k+','markersize',15,'LineWidth',2)
fixPlot(gca,[0 6],[0 6],'attend in','attend out',0:2:6,0:2:6,'sc => sc')

%% histogram of attention mod index for pred dims 
clf; set(gcf,'color','w','pos',[53,1,857,704]);
arrayfun(@(ii) line(subplot(2,2,ii),[0 0],[0 12],'color','k','linewidth',2,'linestyle','--'),1:4)
arrayfun(@(ii) hold(subplot(2,2,ii),'on'),1:4)
histogram(mean(stat.dimIdx_mt_sc_perSes,2),-0.4:0.08:0.4,'edgecolor',cols(1,:),'DisplayStyle','stairs','linewidth',2,'parent',subplot(221));
histogram(mean(stat.dimIdx_mt_mt_perSes,2),-0.4:0.08:0.4,'edgecolor',cols(2,:),'DisplayStyle','stairs','linewidth',2,'parent',subplot(222));
histogram(mean(stat.dimIdx_sc_sc_perSes,2),-0.4:0.08:0.4,'edgecolor',cols(3,:),'DisplayStyle','stairs','linewidth',2,'parent',subplot(223));
histogram(mean(stat.dimIdx_sc_mt_perSes,2),-0.4:0.08:0.4,'edgecolor',cols(4,:),'DisplayStyle','stairs','linewidth',2,'parent',subplot(224));
line(subplot(221),[median(mean(stat.dimIdx_mt_sc_perSes,2)) median(mean(stat.dimIdx_mt_sc_perSes,2))],[0 12],'color',cols(1,:),'linewidth',2,'linestyle','--');
line(subplot(222),[median(mean(stat.dimIdx_mt_mt_perSes,2)) median(mean(stat.dimIdx_mt_mt_perSes,2))],[0 12],'color',cols(2,:),'linewidth',2,'linestyle','--');
line(subplot(223),[median(mean(stat.dimIdx_sc_sc_perSes,2)) median(mean(stat.dimIdx_sc_sc_perSes,2))],[0 12],'color',cols(3,:),'linewidth',2,'linestyle','--');
line(subplot(224),[median(mean(stat.dimIdx_sc_mt_perSes,2)) median(mean(stat.dimIdx_sc_mt_perSes,2))],[0 12],'color',cols(4,:),'linewidth',2,'linestyle','--');
arrayfun(@(ii) fixPlot(subplot(2,2,ii),[-0.4 0.4],[0 12],'attn mod index pred dims','count',-0.4:0.2:0.4,0:2:12),1:4)