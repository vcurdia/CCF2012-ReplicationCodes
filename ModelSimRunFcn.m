function Out=ModelSimRunFcn(x,FileName,o)

% ModelSimRunFcn
%
% The commands inside the loop for each parameter draw in ModelSimRun.
%
% ..............................................................................
%
% Created: July 19, 2011  by Vasco Curdia
% Updated: April 14, 2014 by Vasco Curdia
%
% Copyright 2011-2014 by Vasco Curdia

%% -----------------------------------------------------------------------------

%% evaluate options
oList = fieldnames(o);
for jO=1:length(oList)
    eval(sprintf('%1$s = o.%1$s;',oList{jO}))
end

%% -----------------------------------------------------------------------------

%% IRF

%% REE Mats
Mats = feval(FileName.Mats,x,2,2);
idxR = ismember(StateVar,'r');
% normal state
StateEqMat = Mats.StateEq; %assume idxNoZLB=1
% ZLB
StateEqMat(idxZLB) = StateEqMat(idxNoZLB);
StateEqMat(idxZLB).Gamma0(end,:) = idxR;
StateEqMat(idxZLB).GammaBar(end,:) = 0;
StateEqMat(idxZLB).Gamma1(end,:) = 0;
StateEqMat(idxZLB).Gamma2(end,:) = 0;
StateEqMat(idxZLB).Gamma3(end,:) = 0;
% NoZLB, total debt kept constant
StateEqMat(idxNoZLBCtDebt) = StateEqMat(idxNoZLB);
StateEqMat(idxNoZLBCtDebt).Gamma0(end-1,:) = ismember(StateVar,'BTotMVz');
StateEqMat(idxNoZLBCtDebt).GammaBar(end-1,:) = 0;
StateEqMat(idxNoZLBCtDebt).Gamma1(end-1,:) = 0;
StateEqMat(idxNoZLBCtDebt).Gamma2(end-1,:) = 0;
StateEqMat(idxNoZLBCtDebt).Gamma3(end-1,:) = 0;
% ZLB, total debt kept constant
StateEqMat(idxZLBCtDebt) = StateEqMat(idxZLB);
StateEqMat(idxZLBCtDebt).Gamma0(end-1,:) = ismember(StateVar,'BTotMVz');
StateEqMat(idxZLBCtDebt).GammaBar(end-1,:) = 0;
StateEqMat(idxZLBCtDebt).Gamma1(end-1,:) = 0;
StateEqMat(idxZLBCtDebt).Gamma2(end-1,:) = 0;
StateEqMat(idxZLBCtDebt).Gamma3(end-1,:) = 0;
% Perfect Foresight REE
[REE.NoZLB,eu] = PerfectForesightREE(StateEqMat,RegIdx.NoZLB);
[REE.ZLB,eu] = PerfectForesightREE(StateEqMat,RegIdx.ZLB);
Out.eu = eu;

%% Generate IRFs
ShockSizej = ShockSize;
for jS=1:nShocks2Show
    if strcmp(Shocks2Show{jS},'eBL')
        ShockSizej(jS) = ShockSizej(jS)/...
            REE.NoZLB(end).G2(ismember(StateVar,'XiB'),...
            ismember(ShockVar,'eBL'));
    elseif strcmp(Shocks2Show{jS},'ezeta')
        ShockSizej(jS) = ShockSize(jS)/...
            REE.NoZLB(end).G2(ismember(StateVar,'Xizeta'),...
            ismember(ShockVar,'ezeta'));
    elseif strcmp(Shocks2Show{jS},'em')
        ShockSizej(jS) = ShockSize(jS)/...
            REE.NoZLB(end).G2(ismember(StateVar,'Xim'),...
            ismember(ShockVar,'em'));
    end
end
ShockImpact = eye(nShockVar);
ShockImpact = ShockImpact(:,ismember(ShockVar,Shocks2Show))*diag(ShockSizej);
for jZLB=1:2,ZLBj=ZLBList{jZLB};
    REEj = REE.(ZLBj);
    IRF = REEj(1).G2*ShockImpact;
    for t=2:TReg
        IRF(:,:,t) = REEj(t).G1*IRF(:,:,t-1);
    end
    for t=(TReg+1):nSteps
        IRF(:,:,t) = REEj(end).G1*IRF(:,:,t-1);
    end
    Out.IRF.(ZLBj) = permute(IRF,[1,3,2]);
    clear IRF REEj
end

%% -----------------------------------------------------------------------------

%% Simulated path
for jSP=1:nSimPath, SPj = SimPathList{jSP};
    SimPathj = SimPath.(SPj);
    % Change mats
    StateEqMatSP = StateEqMat;
    nRegs = length(StateEqMatSP);
    EqIdx = SimPathj.EqIdx;
    for jR=1:nRegs
        StateEqMatSP(jR).Gamma0(EqIdx,:) = SimPathj.StateVarIdx;
        StateEqMatSP(jR).GammaBar(EqIdx,:) = 0;
        StateEqMatSP(jR).Gamma1(EqIdx,:) = 0;
        StateEqMatSP(jR).Gamma2(EqIdx,:) = SimPathj.ShockVarIdx;
        StateEqMatSP(jR).Gamma3(EqIdx,:) = 0;
    end
    % Perfect Foresight REE
    ShockPathj = SimPathj.ShockVarIdx*SimPathj.Path;
    clear REESP
    % Simulate
    for jZLB=1:2,ZLBj=ZLBList{jZLB};
        RegIdxj = RegIdx.(ZLBj);
        CheckZLB = 0;
        TRegj = TReg;
        ShockPathjj = ShockPathj;
        while ~CheckZLB
            [REEj,eu] = PerfectForesightREE(StateEqMatSP,RegIdxj,ShockPathjj);
            IRF = REEj(1).GBar+REEj(1).G2*ShockPathj(:,1);
            for t=2:TRegj
                IRF(:,t) = REEj(t).GBar+REEj(t).G1*IRF(:,t-1)+...
                    REEj(t).G2*ShockPathjj(:,t);
            end
            for t=(TRegj+1):nSteps
                IRF(:,t) = REEj(end).GBar+REEj(end).G1*IRF(:,t-1);
            end
            if isCheckZLB
                idxViolateZLB = find(IRF(idxR,:)<-1e-5,1,'first');
                if isempty(idxViolateZLB)
                    CheckZLB = 1;
                else
                    if idxViolateZLB>=TRegj
                        RegIdxj((TRegj+1):(idxViolateZLB+1)) = 1;
                        ShockPathjj(:,(TRegj+1):(idxViolateZLB+1)) = 0;
                        TRegj = idxViolateZLB+1;
                    end
                    RegIdxj(idxViolateZLB) = RegIdxj(idxViolateZLB)+1;
                end
            else
                CheckZLB = 1;
            end
        end
        Out.SimPath.(SPj).(ZLBj) = IRF;
        Out.SimPathRegIdx.(SPj).(ZLBj) = RegIdxj;
        clear IRF REEj
        %     Out.SimPath.(SPj).eu.(ZLBj) = eu;
    end
    Out.SimPath.(SPj).eu = eu;
end


%% -----------------------------------------------------------------------------


