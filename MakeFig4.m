% MakeFig4
%
% Same as ModelSimRun but configured to generate Fig 4.
%
% Uses the setup generated in ModelSimSetup and runs and plots desired
% simulations.
%
% See also:
% GenSymVars, MakeMats, ModelSim, ModelSimSetup
%
% ..............................................................................
%
% Created: July 19, 2011 by Vasco Curdia
% Updated: April 18, 2014 by Vasco Curdia
%
% Copyright 2011-2014 by Vasco Curdia

%% -----------------------------------------------------------------------------

%% Preamble
clear all
tic

%% Choose specification
Spec2Use = 'BLMVB';

%% Setup
load([Spec2Use,'Sim'])
if strcmp(Spec2Use,'BLMVB')
    FileName.Estimation = 'RPRatio';
    FileName.MCMCDrawsRedux = 'RPRatioMCMCDrawsUpdate2Redux';
    clear ShockSize
    ShockSize.eBL = -0.25;
%     ShockSize.em = -0.25/400;
else
    error('Wrong spec...')
end

%% Settings
nSteps = 29;
nSteps2Show = 25;
UseDist = 'postdraws'; % Options: 'postdraws','priordraws','postp500',...
UseEstimation = 1;
nDraws = 1000; %1000
isExtendedZLB = 0;
isExtendedLSAP = 0;
isRedrawOmega = 0;
isLowPriceRig = 0;
isIsolateOmega = 0;
isIsolateZeta = 0;
isNoNomRig = 0;
isCheckZLB = 1;
PlotCompareBands = [16,84];
isSaveData = 0;
isPlotReferenceData = 1;
Bands2Show = [50,60,70,80,90];

%% shocks
Shocks2Show = fieldnames(ShockSize);

%% Simpath
clear SimPath
SimPath.eBL.Path = ShockSize.eBL*...
                     [(1:4)/4,ones(1,8+isExtendedLSAP*8),(7:-1:0)/8];
SimPath.eBL.EqIdx = nStateEq-2;
SimPath.eBL.StateVarIdx = ismember(StateVar,'BLMVz');

%% ZLB and constant debt
nZLB = 4+isExtendedZLB*1;
nCtDebt = length(SimPath.eBL.Path);

%% Plot Options
% Plots2Show = {'IRFNoZLB','IRFZLB','IRFZLBComp','SimPathNoZLB','SimPathZLB',...
%               'SimPathZLBComp'};
Plots2Show = {'SimPathNoZLB'};
% FigPanels2Show = {'Paper','Slides','RealVars','BLMV'};
FigPanels2Show = {'Paper'};
PlotSupTitle = 0;
PlotPrint = 1;
PlotEpsDelete = 1;
PlotSuffix = ['_',UseDist];
if isRedrawOmega
    PlotSuffix = [PlotSuffix,'_HighSegmentation'];
end
if isLowPriceRig
    PlotSuffix = [PlotSuffix,'_LowPriceRig'];
end
if isExtendedZLB
    PlotSuffix = [PlotSuffix,'_ExtendedZLB'];
end
if isExtendedLSAP
    PlotSuffix = [PlotSuffix,'_ExtendedLSAP'];
end
if isIsolateOmega
    PlotSuffix = [PlotSuffix,'_IsolateOmega'];
end
if isIsolateZeta
    PlotSuffix = [PlotSuffix,'_IsolateZeta'];
end
if isNoNomRig
    PlotSuffix = [PlotSuffix,'_NoNomRig'];
end
if isCheckZLB
    PlotSuffix = [PlotSuffix,'_CheckZLB'];
end
clear FigPanel
if ismember('Paper',FigPanels2Show)
    FigPanel.Paper.Vars2Show = {...
        'dY','Yz',...
        'pi','r',...
        'rL','RP',...
        };
    FigPanel.Paper.Vars2ShowPretty = {...
        'Output Growth','Output Level',...
        'Inflation','FFR',...
        '10Yr Yield','Risk Premium',...
        };
    FigPanel.Paper.Scale2Show = 100*[...
        4,1,...
        4,4,...
        4,4,...
        ];
    FigPanel.Paper.FigShape = {3,2};
end
if ismember('Slides',FigPanels2Show)
    FigPanel.Slides.Vars2Show = {...
        'dY','pi','r',...
        'Yz','rL','RP',...
        };
    FigPanel.Slides.Vars2ShowPretty = {...
        'Output Growth','Inflation','FFR',...
        'Output Level','10Yr Yield','Risk Premium',...
        };
    FigPanel.Slides.Scale2Show = 100*[...
        4,4,4,...
        1,4,4,...
        ];
    FigPanel.Slides.FigShape = {2,3};
end
if ismember('RealVars',FigPanels2Show)
    FigPanel.RealVars.Vars2Show = {...
        'Yz','Czu','Czr','L',...
        'u','Kzbar','Iz','q',...
        'Tz','wz','dw','BLMVz',...
        };
    FigPanel.RealVars.Vars2ShowPretty = {...
        'Yz','Czu','Czr','L',...
        'u','Kzbar','Iz','q',...
        'Tz','wz','dw','BLMVz',...
        };
    FigPanel.RealVars.Scale2Show = 100*[...
        1,1,1,1,...
        1,1,1,1,...
        1,1,1,1,...
        ];
    FigPanel.RealVars.FigShape = {3,4};
end
if ismember('BLMV',FigPanels2Show)
    FigPanel.BLMV.Vars2Show = {...
        'BLMVz',...
        };
    FigPanel.BLMV.Vars2ShowPretty = {...
        'Long term debt (MV)',...
        };
    FigPanel.BLMV.Scale2Show = 100*[...
        1,...
        ];
    FigPanel.BLMV.FigShape = {1,1};
end

LineWidth = [1.5,1.5];
ShadeColor = 1.3*[0.45,0.45,0.5];
LineColor = vcColorScheme;
LineColor = LineColor([3,1],:);

%% -----------------------------------------------------------------------------

%% Set parameters
SetParam = struct;

if isLowPriceRig
    SetParam.zetap = 0.75;
end

if isNoNomRig
    SetParam.zetap = 0.001;
    SetParam.zetaw = 0.001;
end


%% -----------------------------------------------------------------------------

%% Load Parameters from estimation and make sure they coincide with model setup
if UseEstimation
    ParamsOld = Params;
    ParamNamesOld = {ParamsOld(:).name};
    npOld = np;
    load(FileName.Estimation,'Params','np')
    ParamNames = {Params(:).name};
    idxParamMatch = (np==npOld);
    if idxParamMatch, idxParamMatch = strcmp(ParamNamesOld,ParamNames); end
    if ~all(idxParamMatch==1)
        fprintf(2,['Parameters in the estimation file differ from the ',...
                   'parameters in model setup.\n']);
        fprintf(2,['See below the list of parameters in the estimation ',...
                   '(left) and in model setup (right):\n']);
        fprintf(2,'     %20s %20s\n','estimation','model setup');
        for j=1:max(np,npOld)
            if j<=np,Paramj=ParamNames{j};else Paramj = '';end
            if j<=npOld,ParamOldj=ParamNamesOld{j};else ParamOldj = '';end
            fprintf(2,'%3.0f: %20s %20s\n',j,Paramj,ParamOldj);
        end
    end
end

%% Prepare Draws
if strcmp(UseDist,'priordraws')
    xd = feval(FileName.GenPriorDraw,nDraws);
elseif strcmp(UseDist,'postdraws')
    load(FileName.MCMCDrawsRedux,'xd')
    nxd = size(xd,2);
    if nxd<nDraws
        fprintf(2,...
            'Number of MCMC draws available is less than desired draws.\n');
        fprintf(2,'Proceeding with those available.\n');
    else
        nThinningSim = max(floor(nxd/nDraws),1);
        xd = xd(:,1:nThinningSim:end);
    end
else
    if ~isfield(Params,UseDist)
        fprintf(2,'Did not recognize distribution to use. Cannot proceed.\n');
    end
    xd = [Params(:).(UseDist)]';
end
nDraws = size(xd,2);

%% Redraw omegau
if isRedrawOmega && nDraws>1
    idx_omegau = ismember(ParamNames,'omegau');
    omegaud = xd(idx_omegau,:);
    omegaud = omegaud(omegaud<=prctile(omegaud,50));
    xd(idx_omegau,:) = omegaud(randi(length(omegaud),1,nDraws));
end %isRedrawOmega

%% Isolate omega
if isIsolateOmega && nDraws>1
    idx_omegau = ismember(ParamNames,'omegau');
    omegaud = xd(idx_omegau,:);
    xd = repmat([Params(:).postp500]',1,nDraws);
    xd(idx_omegau,:) = omegaud;
end %isIsolateOmega

%% Isolate zeta
if isIsolateZeta && nDraws>1
    idx_dzetap = ismember(ParamNames,'dzetap');
    dzetapd = xd(idx_dzetap,:);
    xd = repmat([Params(:).postp500]',1,nDraws);
    xd(idx_dzetap,:) = dzetapd;
end %isIsolateZeta

%% Set parameters (if any)
ParNames = {Params(:).name};
SetParamNames = fieldnames(SetParam);
for jp=1:length(SetParamNames)
    xd(ismember(ParNames,SetParamNames{jp}),:) = SetParam.(SetParamNames{jp});
end

%% Prepare Simulations
if ~exist('Shocks2Show','var') || strcmp(Shocks2Show{1},'All')
    Shocks2Show = ShockVar;
end
nShocks2Show = length(Shocks2Show);
ShockSizeTmp = ones(nShocks2Show,1);
ShockSizeList = fieldnames(ShockSize);
for jS=1:length(ShockSizeList),Sj = ShockSizeList{jS};
    idxS = ismember(Shocks2Show,Sj);
    ShockSizeTmp(idxS) = ShockSize.(Sj);
end %jS
ShockSize = ShockSizeTmp;
clear ShockSizeTmp

ZLBList = {'NoZLB','ZLB'};

%% Prepare options for simulations
clear o
o.ShockSize = ShockSize;
o.Shocks2Show = Shocks2Show;
o.nShocks2Show = nShocks2Show;
o.StateVar = StateVar;
o.ShockVar = ShockVar;
o.nShockVar = nShockVar;
o.idxNoZLB = 1;
o.idxZLB = 2;
o.idxNoZLBCtDebt = 3;
o.idxZLBCtDebt = 4;
o.RegIdx.NoZLB = [repmat(o.idxNoZLBCtDebt,1,nCtDebt),o.idxNoZLB];
o.RegIdx.ZLB = [repmat(o.idxZLBCtDebt,1,min(nZLB,nCtDebt)),...
    repmat(o.idxNoZLBCtDebt,1,nCtDebt-nZLB),repmat(o.idxZLB,1,nZLB-nCtDebt),...
    o.idxNoZLB];
TReg.NoZLB = length(o.RegIdx.NoZLB);
TReg.ZLB = length(o.RegIdx.ZLB);
o.TReg = max(TReg.NoZLB,TReg.ZLB);
o.ZLBList = ZLBList;
o.nSteps = nSteps;
o.SimPath = SimPath;
o.SimPathList = fieldnames(SimPath);
o.nSimPath = length(o.SimPathList);
for jSP=1:o.nSimPath,SPj = o.SimPathList{jSP};
    o.SimPath.(SPj).ShockVarIdx = ismember(ShockVar,SPj);
    TReg.(SPj) = size(SimPath.(SPj).Path,2);
    o.TReg = max(o.TReg,TReg.(SPj));
end
o.RegIdx.NoZLB = [o.RegIdx.NoZLB,repmat(o.RegIdx.NoZLB(end),1,...
    o.TReg-TReg.NoZLB)];
o.RegIdx.ZLB = [o.RegIdx.ZLB,repmat(o.RegIdx.ZLB(end),1,o.TReg-TReg.ZLB)];
for jSP=1:o.nSimPath,SPj = o.SimPathList{jSP};
    o.SimPath.(SPj).Path = [o.SimPath.(SPj).Path,...
        repmat(o.SimPath.(SPj).Path(:,end),1,o.TReg-TReg.(SPj))];
end
o.isCheckZLB = isCheckZLB;

%% Run simulations for alternative parameter draws
fprintf('Generating simulations for alternative parameter values...\n')
clear Out
parfor jd=1:nDraws
    Out(jd) = ModelSimRunFcn(xd(:,jd),FileName,o);
end
for jd=nDraws:-1:1
    if ~all(Out(jd).eu(:)==1)
        Out(jd) = [];
    end
end
nDraws = length(Out);

%% -----------------------------------------------------------------------------

%% Plot results

%% preliminary stuff
tid = 0:nSteps2Show-1;
FigPanelList = fieldnames(FigPanel);
nFigPanel = length(FigPanelList);
nPlotCompareBands = length(PlotCompareBands);
if isPlotReferenceData
    load('ReferencePlotData')
end
for jFig=1:nFigPanel, Figj = FigPanelList{jFig};
    FigPanelFields = fieldnames(FigPanel.(Figj));
    for j=1:length(FigPanelFields)
        eval(sprintf('%1$s = FigPanel.(Figj).%1$s;',FigPanelFields{j}))
    end
    nVars2Show = length(Vars2Show);
    PlotDir = ['Plots_',Figj,'/'];
    if~isdir(PlotDir),mkdir(PlotDir),end
    
    %% Plot IRFs with bands
    for jZLB=1:2,ZLBj = ZLBList{jZLB};
        if ~ismember(['IRF',ZLBj],Plots2Show), continue, end
        for jS=1:nShocks2Show
            Shockj = Shocks2Show{jS};
            figure
            for jV=1:nVars2Show
                Varj = Vars2Show{jV};
                idxV = ismember(StateVar,Varj);
                hsubp(jV) = subplot(FigShape{:},jV);
                PlotData = zeros(nDraws,nSteps2Show);
                for jD=1:nDraws
                    PlotData(jD,:) = Scale2Show(jV).*Out(jD).IRF.(ZLBj)(idxV,...
                        1:nSteps2Show,jS);
                end
                if nDraws>1
                    if isPlotReferenceData...
                            && ismember(Varj,fieldnames(PlotData2Save))
                        AltData = PlotData2Save.(Varj);
                    else 
                        AltData = [];
                    end
                    vcPlotDistBands(tid,PlotData,'AltData',AltData,...
                        'LineColor',LineColor,'ShadeColor',ShadeColor,...
                        'LineWidth',LineWidth)
                else
                    plot(tid,PlotData,'Color',MedianColor,'LineWidth',LineWidth)
                    if isPlotReferenceData...
                            && ismember(Varj,fieldnames(PlotData2Save))
                        hold on
                        plot(tid,PlotData2Save.(Varj),'--','Color',...
                            AltColor,'LineWidth',LineWidth)
                    end
                end
                title(Vars2ShowPretty{jV})
                axis tight
                set(gca,'XTick',0:4:nSteps)
                yBounds = ylim;
                yBounds = [min([yBounds(1);-0.0125]),max([yBounds(2);0.0125])];
                ylim(yBounds)
            end %jV
            if PlotSupTitle
                suptitle([Shocks2Show{jS},' - ',ZLBj])
            end
            if PlotPrint
                PlotName = [PlotDir,Spec2Use,'_IRF_',Shockj,'_',ZLBj,...
                    '_',Figj,PlotSuffix,'.eps'];
                print('-depsc2',PlotName)
                eval(['!epstopdf ',PlotName])
                if PlotEpsDelete,delete(PlotName),end
            end
            clear hsubp
        end %jS
    end
    
    %% Compare IRFs with ZLB vs noZLB
    if ismember('IRFZLBComp',Plots2Show)
        for jS=1:nShocks2Show
            Shockj = Shocks2Show{jS};
            figure
            for jV=1:nVars2Show
                Varj = Vars2Show{jV};
                idxV = ismember(StateVar,Varj);
                hsubp(jV) = subplot(FigShape{:},jV);
                PlotDataNoZLB = zeros(nDraws,nSteps2Show);
                PlotDataZLB = zeros(nDraws,nSteps2Show);
                for jD=1:nDraws
                    PlotDataNoZLB(jD,:) = Scale2Show(jV).*...
                        Out(jD).IRF.NoZLB(idxV,1:nSteps2Show,jS);
                    PlotDataZLB(jD,:) = Scale2Show(jV).*...
                        Out(jD).IRF.ZLB(idxV,1:nSteps2Show,jS);
                end
                if nDraws>1
                    PlotDataNoZLB = prctile(PlotDataNoZLB,[50,PlotCompareBands]);
                    PlotDataZLB = prctile(PlotDataZLB,[50,PlotCompareBands]);
                end
                plot(tid,PlotDataNoZLB(1,:),'-b',tid,PlotDataZLB(1,:),'--r',...
                    'LineWidth',LineWidth)
                hold on
                plot(tid,zeros(size(tid)),'k:')
                for jBand = 1:nPlotCompareBands
                    plot(tid,PlotDataNoZLB(1+jBand,:),':b',...
                        tid,PlotDataZLB(1+jBand,:),':r','LineWidth',1)
                end
                title(Vars2ShowPretty{jV})
                axis tight
                yBounds = ylim;
                yBounds = [min([yBounds(1);-0.0125]),max([yBounds(2);0.0125])];
                ylim(yBounds)
            end %jV
            if PlotSupTitle
                suptitle(Shocks2Show{jS})
            end
            if all([FigShape{:}]==1)
                legend('Ignoring ZLB','Accounting for ZLB',...
                    'Orientation','horizontal','Location','SO');
            else
                hleg = legend('Ignoring ZLB','Accounting for ZLB',...
                    'Orientation','horizontal');
                legPos = get(hleg,'Position');
                xL = get(hsubp((FigShape{1}-1)*FigShape{2}+1),'Position');
                xR = get(hsubp((FigShape{1}-1)*FigShape{2}),'Position');
                legPos(1) = xL(1)+(xR(1)-xL(1))/2+(xL(3)-legPos(3))/2;
                legPos(2) = 0;
                set(hleg,'Position',legPos)
            end
            % print fig
            if PlotPrint
                PlotName = [PlotDir,Spec2Use,'_IRF_',Shockj,'_ZLBComp_',Figj,...
                    PlotSuffix,'.eps'];
                print('-depsc2',PlotName)
                eval(['!epstopdf ',PlotName])
                if PlotEpsDelete,delete(PlotName),end
            end
            clear hsubp
        end %jS
    end %
    
    %% Plot SimPath with bands
    for jZLB=1:2,ZLBj = ZLBList{jZLB};
        if ~ismember(['SimPath',ZLBj],Plots2Show), continue, end
        for jSP=1:o.nSimPath
            SPj = o.SimPathList{jSP};
            figure
            for jV=1:nVars2Show
                Varj = Vars2Show{jV};
                idxV = ismember(StateVar,Varj);
                hsubp(jV) = subplot(FigShape{:},jV);
                PlotData = zeros(nDraws,nSteps2Show);
                for jD=1:nDraws
                    PlotData(jD,:) = Scale2Show(jV).*...
                        Out(jD).SimPath.(SPj).(ZLBj)(idxV,1:nSteps2Show);
                end
                if nDraws>1
                    if isPlotReferenceData...
                            && ismember(Varj,fieldnames(PlotData2Save))
                        AltData = PlotData2Save.(Varj);
                    else 
                        AltData = [];
                    end
                    vcPlotDistBands(tid,PlotData,'Bands2Show',Bands2Show,...
                        'AltData',AltData,...
                        'LineColor',LineColor,'ShadeColor',ShadeColor,...
                        'LineWidth',LineWidth)
                    PlotData = prctile(PlotData,50,1);
                else
                    plot(tid,PlotData,'Color',MedianColor,'LineWidth',LineWidth)
                    if isPlotReferenceData...
                            && ismember(Varj,fieldnames(PlotData2Save))
                        hold on
                        plot(tid,PlotData2Save.(Varj),'--','Color',...
                            AltColor,'LineWidth',LineWidth)
                    end
                end
                if isSaveData && strcmp(ZLBj,'ZLB') && strcmp(SPj,'eBL')...
                        && ~isRedrawOmega && ~isLowPriceRig
                    PlotData2Save.(Varj) = PlotData;
                end
                title(Vars2ShowPretty{jV})
                axis tight
                set(gca,'XTick',0:4:nSteps)
                yBounds = ylim;
                yBounds = [min([yBounds(1);-0.0125]),max([yBounds(2);0.0125])];
                ylim(yBounds)
            end %jV
            if isSaveData && strcmp(ZLBj,'ZLB') && strcmp(SPj,'eBL')...
                    && ~isRedrawOmega && ~isLowPriceRig
                save ReferencePlotData PlotData2Save
            end
            if PlotSupTitle
                suptitle(['Simulated Path for ',SPj,' - ',ZLBj])
            end
            if PlotPrint
                PlotName = [PlotDir,Spec2Use,'_SimPath_',SPj,'_',ZLBj,...
                    '_',Figj,PlotSuffix,'.eps'];
                print('-depsc2',PlotName)
                eval(['!epstopdf ',PlotName])
                if PlotEpsDelete,delete(PlotName),end
            end
            clear hsubp
        end %jS
    end
    
    %% Compare SimPath with ZLB vs noZLB
    if ismember('SimPathZLBComp',Plots2Show)
        for jSP=1:o.nSimPath
            SPj = o.SimPathList{jSP};
            figure
            for jV=1:nVars2Show
                Varj = Vars2Show{jV};
                idxV = ismember(StateVar,Varj);
                hsubp(jV) = subplot(FigShape{:},jV);
                PlotDataNoZLB = zeros(nDraws,nSteps2Show);
                PlotDataZLB = zeros(nDraws,nSteps2Show);
                for jD=1:nDraws
                    PlotDataNoZLB(jD,:) = Scale2Show(jV).*...
                        Out(jD).SimPath.(SPj).NoZLB(idxV,1:nSteps2Show);
                    PlotDataZLB(jD,:) = Scale2Show(jV).*...
                        Out(jD).SimPath.(SPj).ZLB(idxV,1:nSteps2Show);
                end
                if nDraws>1
                    PlotDataNoZLB = prctile(PlotDataNoZLB,[50,PlotCompareBands]);
                    PlotDataZLB = prctile(PlotDataZLB,[50,PlotCompareBands]);
                end
                plot(tid,PlotDataNoZLB(1,:),'-b',tid,PlotDataZLB(1,:),'--r',...
                    'LineWidth',LineWidth)
                hold on
                plot(tid,zeros(size(tid)),'k:')
                for jBand = 1:nPlotCompareBands*(nDraws>1)
                    plot(tid,PlotDataNoZLB(1+jBand,:),':b',...
                        tid,PlotDataZLB(1+jBand,:),':r','LineWidth',1)
                end
                title(Vars2ShowPretty{jV})
                axis tight
                yBounds = ylim;
                yBounds = [min([yBounds(1);-0.0125]),max([yBounds(2);0.0125])];
                ylim(yBounds)
            end %jV
            if PlotSupTitle
                suptitle(['Simulated Path for ',SPj])
            end
            if all([FigShape{:}]==1)
                legend('Ignoring ZLB','Accounting for ZLB',...
                    'Orientation','horizontal','Location','SO');
            else
                hleg = legend('Ignoring ZLB','Accounting for ZLB',...
                    'Orientation','horizontal');
                legPos = get(hleg,'Position');
                xL = get(hsubp((FigShape{1}-1)*FigShape{2}+1),'Position');
                xR = get(hsubp((FigShape{1}-1)*FigShape{2}),'Position');
                legPos(1) = xL(1)+(xR(1)-xL(1))/2+(xL(3)-legPos(3))/2;
                legPos(2) = 0;
                set(hleg,'Position',legPos)
            end
            % print fig
            if PlotPrint
                PlotName = [PlotDir,Spec2Use,'_SimPath_',SPj,'_ZLBComp_',...
                    Figj,PlotSuffix,'.eps'];
                print('-depsc2',PlotName)
                eval(['!epstopdf ',PlotName])
                if PlotEpsDelete,delete(PlotName),end
            end
            clear hsubp
        end %jS
    end %
    
end % FigPanel

%% -----------------------------------------------------------------------------

%% elapsed time
fprintf('\n%s\n\n',vctoc)

%% Save environment
% save(FileName.Output)

%% -----------------------------------------------------------------------------

