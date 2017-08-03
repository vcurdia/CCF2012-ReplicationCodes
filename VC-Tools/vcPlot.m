function varargout = vcPlot(y,varargin)

% vcPlot
%
% Plots several lines with optional light effects.
% 
% Syntax:
%   vcPlot(y)
%   vcPlot(x,y)
%   vcPlot(...,OptionsStructure,...)
%   vcPlot(...,'PropertyName',PropertyValue,...)
%   h = vcPlot(...)
%
% If no arguments are specified then an example is shown.
%
% See also:
% vcFigure, vcPlotDistBands, vcColorScheme
%
% .........................................................................
%
% Created: April 24, 2013
% Updated: October 14, 2013
%   - Allows for only one input, assumed to be y.
%   - allows for options to be specified in either a structure or a list of
%     options (through varargin), or both.
% Updated: October 15, 2013
%   - Allows for multiple layers of light.
%   - Simplified code
%   - If there are not enough ticks for the lines with marker only then
%     the code automatically interpolates the data to add more marker ticks
%
% Copyright (2013) by Vasco Rafael da Silva Curdia

%% ------------------------------------------------------------------------

%% Check main input
if nargin==0
  x = 1:10;
  y(1,:) = -2+x;
  y(2,:) = 2+log(x);
  y(3,:) = 1+2*log(x);
  y(4,:) = 5-2*log(x);
  y(5,:) = 7-2*log(x);
  y(6,:) = 1-2*log(x);
  % y(2,1:1) = NaN;
  o = struct;
elseif nargin==1
  x = 1:size(y,2);
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

%% Dimensions
[ny,nx] = size(y);

%% Default options
oDefault.LineColor = vcColorScheme;
oDefault.LineStyle = {'-','--','-.',':',':',':'};
oDefault.LineMarker = {'none','none','none','o','s','^'};
% oDefault.LineStyle = {'-','--','-','-','none'};
% oDefault.LineMarker = {'none','none','o','s','^'};
% oDefault.LineMarkerSize = [2,2,2,2,2,2];
% oDefault.LineWidth = [1.5,1.5,1.5,0.5,0.5,0.5];
oDefault.LineMarkerSize = [1,1,1,3,3,3];
oDefault.LineWidth = [1.5,1.5,1.5,0.5,0.5,0.5];
oDefault.ShowLight = 1;
oDefault.LightLayers = 2;
oDefault.Brightness = 0.3;
oDefault.ShowZeroLine = 1;
oDefault.ZeroLineColor = 'k';
oDefault.ZeroLineStyle = '-';
oDefault.ZeroLineWidth = 0.5;
oDefault.ShowLegend = (ny>1);
oDefault.LegendLocation = 'Best';
oDefault.LegendOrientation = 'vertical'; %'vertical','horizontal'

%% Apply defaults to unspecified options
oList = fieldnames(oDefault);
for jo=1:length(oList)
  oName = oList{jo};
  if ~isfield(o,oName)
    o.(oName) = oDefault.(oName);
  end
end

%% Additional options if not yet set
if ~isfield(o,'MarkerTicks')
  o.MarkerTicks = max(1,round(40/nx));
end

%% Check x
if ~exist('x','var')
  x = 1:nx;
end

%% Some values
CurrentHold = ishold;
h.Lines = zeros(1+o.LightLayers*o.ShowLight,ny);
xx = x(1):1/o.MarkerTicks:x(nx);
nxx = length(xx);

%% Plot zero line
if o.ShowZeroLine
  h.ZeroLine = plot(xx,zeros(1,nxx),o.ZeroLineStyle,...
    'Color',o.ZeroLineColor,'LineWidth',o.ZeroLineWidth);
  set(get(get(h.ZeroLine,'Annotation'),'LegendInformation'),...
      'IconDisplayStyle','off');
  if ~ishold,hold on,end
else
  h.ZeroLine = [];
end

%% Plot Line with light effect
nLayers = 1+o.ShowLight*o.LightLayers;
for j=1:ny
  yj = y(j,:);
  yDiff = diff(yj);
  yy = zeros(1,nxx);
  idx = find(ismember(xx,x));
  yy(idx) = yj;
  for jM=1:(o.MarkerTicks-1)
    idx = idx(1:nx-1)+1;
    yy(idx) = yj(1:nx-1)+yDiff*jM/o.MarkerTicks;
  end
  for jLayer=1:nLayers
    if nLayers>1
      jColor = o.LineColor(j,:)+(1-o.LineColor(j,:))*o.Brightness*...
        (jLayer-1)/(nLayers-1);
    else
      jColor = o.LineColor(j,:);
    end
    jSize = 1-(jLayer-1)/nLayers;
    h.Lines(jLayer,j) = plot(xx,yy,'LineStyle',o.LineStyle{j},...
      'Color',jColor,'MarkerFaceColor',jColor,...
      'Marker',o.LineMarker{j},'MarkerSize',o.LineMarkerSize(j)*jSize,...
      'LineWidth',o.LineWidth(j)*jSize);
    if ~ishold,hold on,end
  end
end

if isfield(o,'FontSize')
  set(gca,'FontSize',o.FontSize)
end

%% Show Legend
if o.ShowLegend
  if ~isfield(o,'LegendString')
    for j=1:ny
      o.LegendString{j} = sprintf('Line %.0f',j);
    end
  end
  h.LegendItems = h.Lines(1,:);
  if isfield(o,'LegendItems')
    h.LegendItems = [h.LegendItems,o.LegendItems];
  end
  nLegendItems = length(h.LegendItems);
  [h.Legend,h.LegendObj,~,~] = legend(h.LegendItems,o.LegendString,...
    'Location',o.LegendLocation,'Orientation',o.LegendOrientation);
  for j=1:ny*o.ShowLight
    idx = nLegendItems+(j-1)*2+1;
    xx = get(h.LegendObj(idx),'XData');
    yy = get(h.LegendObj(idx),'YData');
    set(h.LegendObj(idx),'Visible','off')
    set(h.LegendObj(idx+1),'Visible','off')
    for jLayer=1:nLayers
      if nLayers>1
        jColor = o.LineColor(j,:)+(1-o.LineColor(j,:))*o.Brightness*...
          (jLayer-1)/(nLayers-1);
      else
        jColor = o.LineColor(j,:);
      end
      jSize = 1-(jLayer-1)/nLayers;
      h.LegendLines(jLayer,j) = line('XData',xx,'YData',yy,...
        'Parent',h.Legend,'LineStyle',o.LineStyle{j},...
        'Color',jColor,'MarkerFaceColor',jColor,...
        'Marker',o.LineMarker{j},'MarkerSize',o.LineMarkerSize(j)*jSize,...
        'LineWidth',o.LineWidth(j)*jSize);
    end
  end
else
  h.Legend = [];
  h.LegendObj =[];
end

%% Restore hold status
if ~CurrentHold
  hold off
end

%% end example
if nargin==0
  hold off
  axis tight
end

%% Exit
if nargout>0
  h.Options = o;
  varargout{1} = h;
end

%% ------------------------------------------------------------------------
