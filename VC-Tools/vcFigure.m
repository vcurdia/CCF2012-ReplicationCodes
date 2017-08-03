function [varargout] = vcFigure(y,varargin)

% vcFigure
%
% Creates figure with multiple plots (or just one). Plots can have multiple
% lines in each, but all need to have the same number of lines.
%
% Syntax:
%
%   vcFigure(y)
%   vcFigure(x,y)
%   vcFigure(...,OptionsStructure)
%   vcFigure(...,'PropertyName',PropertyValue,...)
%   h = vcFigure(...)
%
% If no arguments are specified then an example is shown.
%
% See also:
% vcPlot, vcPlotDistBands, vcColorScheme
%
% .........................................................................
%
% Created: October 11, 2013 by Vasco Curdia
% Updated: October 14, 2013 by Vasco Curdia
%   - allows for options to be specified in either a structure or a list of 
%     options (through varargin), or both.
%
% Copyright 2013 by Vasco Curdia

%% ------------------------------------------------------------------------

%% Options

%% Check Inputs
if nargin==0
%   % Example with bands, including alternatives
%   x = 1:25;
%   y = rand(1000,length(x),9);
%   o.Plot.AltData = [...
%     -0.1+x/25;
%     0.4+x/25;
%     1.3-x/25;
%     ];
%   o.PlotBands = 1;
%   o.Plot.Bands2Show = [50,70,90];
  % Example with lines
  x = 1:10;
  y(1,:) = -2+x;
  y(2,:) = 2+log(x);
  y(3,:) = 1+2*log(x);
  y(4,:) = 5-2*log(x);
  y(5,:) = 7-2*log(x);
  y = repmat(y,[1,1,9]);
  o = struct;
elseif nargin==1
  o = struct;
else
  if isnumeric(varargin{1})
    x = y;
    y = varargin{1};
    varargin(1) = [];
  end
  if isstruct(varargin{1})
    o = varargin{1};
    varargin(1) = [];
  end
  for jo=1:(length(varargin)/2) 
    o.(varargin{(jo-1)*2+1}) = varargin{jo*2};
  end
end

%% Some dimensions
[ny,nx,nPlots] = size(y);

%% Default Options
oDefault.NewFig = 1;
oDefault.Visible = 'on';
oDefault.FigShape = [1,1];
oDefault.PlotBands = 0;
oDefault.Plot = struct;
oDefault.AxisTight = 1;
oDefault.YMinScale = 0.1;
oDefault.YSlack = 0.05;
oDefault.ShowLegend = (ny>1);
oDefault.LegendLocation = 'SouthOutside'; % 'SouthOutside','EmptySlot'

%% Apply defaults to unspecified options
oList = fieldnames(oDefault);
for jo=1:length(oList)
  oName = oList{jo};
  if ~isfield(o,oName)
    o.(oName) = oDefault.(oName);
  end
end

%% Update options
if ~isfield(o.Plot,'ShowLegend')
  o.ShowLegend = (ny>1) && (~o.PlotBands);
end
if ~isfield(o.Plot,'LegendOrientation') && nPlots>1 ...
    && strcmp(o.LegendLocation,'SouthOutside')
  o.Plot.LegendOrientation = 'horizontal';
end
if ~(nPlots>1 && (strcmp(o.LegendLocation,'SouthOutside') || ...
    strcmp(o.LegendLocation,'EmptySlot')))
  o.Plot.LegendLocation = o.LegendLocation;
end

%% ------------------------------------------------------------------------

%% Create figure

%% Initiate figure
if o.NewFig
  h.Figure = figure('Visible',o.Visible);
else
  h.Figure = gcf;
  set(gcf,'Visible',o.Visible)
  clf
end

%% Check x
if ~exist('x','var')
  x = 1:nx;
end

%% Check Figure Shape
jDim = 2;
while prod(o.FigShape)<nPlots
  o.FigShape(jDim) = o.FigShape(jDim)+1;
  jDim = ~(jDim-1)+1;
end

%% Plot data
for jPlot=1:nPlots
  h.SubPlot(jPlot) = subplot(o.FigShape(1),o.FigShape(2),jPlot);
  yj = y(:,:,jPlot);
  oPlot = o.Plot;
  if o.ShowLegend && jPlot==nPlots
    oPlot.ShowLegend = 1;
  else
    oPlot.ShowLegend = 0;
  end
  if o.PlotBands
    h.Plot{jPlot} = vcPlotDistBands(x,yj,oPlot);
  else
    h.Plot{jPlot} = vcPlot(x,yj,oPlot);
  end
  if isfield(o,'TitleList')
    hh = title(o.TitleList{jPlot});
    if isfield(o,'TitleFontSize')
      set(hh,'FontSize',o.TitleFontSize)
    end
  end
  if o.AxisTight, axis tight, end
  if isfield(o,'XLim')
    set(gca,'XLim',o.XLim)
  end
  if isfield(o,'XTick')
    set(gca,'XTick',o.XTick)
  end
  if isfield(o,'XTickLabel')
    set(gca,'XTickLabel',o.XTickLabel)
  end
  if isfield(o,'YLim')
    yBounds = o.YLim(1,1,jPlot);
  else
    yMinScale = o.YMinScale(1+(jPlot-1)*(length(o.YMinScale)>1));
    yBounds = ylim;
    yBounds = [min([1+o.YSlack,-o.YSlack]*yBounds',-yMinScale),...
               max([-o.YSlack,1+o.YSlack]*yBounds',yMinScale)];
  end
  ylim(yBounds)
  if isfield(o,'YGrid')
    set(gca,'YGrid',o.YGrid)
  end
  if isfield(o,'AxesFontSize')
    set(gca,'FontSize',o.AxesFontSize)
    if o.ShowLegend && jPlot==nPlots && h.Plot{jPlot}.Options.ShowLight
      set(h.Plot{jPlot}.LegendObj(ny+1:end),'Visible','off')
      for jy=1:ny
        set(h.Plot{jPlot}.LegendLines(:,jy),'XData',...
          get(h.Plot{jPlot}.LegendObj(ny+2*(jy-1)+1),'XData'))
      end
    end
  end
end
if o.ShowLegend && nPlots>1
  legPos = get(h.Plot{jPlot}.Legend,'Position');
  xIdx = (max(1,o.FigShape(1)-1))*o.FigShape(2);
  xR = get(h.SubPlot(xIdx),'Position');
  xL = get(h.SubPlot(min(nPlots,xIdx+1)),'Position');
  if strcmp(o.LegendLocation,'SouthOutside')
    legPos(1) = xL(1)+(xR(1)-xL(1))/2+(xL(3)-legPos(3))/2;
    legPos(2) = 0;
  elseif strcmp(o.LegendLocation,'EmptySlot')
    if nPlots<prod(o.FigShape)
      legPos(1) = xR(1)+(xR(3)-legPos(3))/2;
      legPos(2) = xL(2)+(xL(4)-legPos(4))/2;
    end
  end
  set(h.Plot{jPlot}.Legend,'Position',legPos)
end
if isfield(o,'SupTitle')
  h.SupTitle = suptitle(o.SupTitle);
  if isfield(o,'SupTitleFont')
    set(h.SupTitle,'fontsize',o.SupTitleFont)
  end
end


%% ------------------------------------------------------------------------

%% Exit
if nargout>0
  h.XData = x;
  h.YData = y;
  h.Options = o;
  varargout = {h};
end

%% ------------------------------------------------------------------------
