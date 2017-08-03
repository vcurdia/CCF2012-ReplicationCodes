function varargout = vcPlotDistBands(y,varargin)

% vcPlotDistBands
%
% Plots median and percentile bands for matrix X.
%
% Usage: 
%   vcPlotDistBands(y)
%   vcPlotDistBands(x,y)
%   vcPlotDistBands(...,OptionsStructure,...)
%   vcPlotDistBands(...,'PropertyName',PropertyValue,...)
%   h = vcPlotDistBands(...)
%
% If no arguments are specified then an example is shown.
%
% Required input argument:
%
%   y
%   Matrix containing data. Percentiles computed along first dimension.
%
% Options:
%
%   Bands2Show
%   Percent intervals to be shown in the plots, centered around median.
%   Default: [50,60,70,80,90]
%
%   MedianColor
%   Color of median line.
%   Default: [0,0,0.7]
%
%   ShadeColor
%   Base color for bands.
%   Default: [0.2,0.6,0.5]
%
%   LineWidth
%   Width of median line.
%   Default: 1.5
%
%   isZeroLine
%   If 1 plots the zero line. If 0 it does not plot a zero line.
%   Default: 1
%
%   ZeroLineColor
%   Color for the zero line.
%   Default: 'k'
%
%   ZeroLineStyle
%   Style for the zero line.
%   Default: ':'
%
%   tid
%   x axis values.
%   Default: 1:T
%
% See also:
% vcPlot, vcFigure, vcColorScheme
%
% ..............................................................................
%
% Created: October 30, 2008 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% Updated: October 14, 2013 by Vasco Curdia
%   - allows for options to be specified in either a structure or a list of
%     options (through varargin), or both.
% Updated: October 15, 2013 by Vasco Curdia
%   - uses vcPlot to present lines and legend if desired.
% Updated: October 17, 2013 by Vasco Curdia
%   - progression of shades allows now more control to caller.
%   - if shade color is not set in the options, then it is set in a mix 
%     between mostly grey and a hint of the color of the first line.
% 
% Copyright 2008-2013 by Vasco Curdia

%% ------------------------------------------------------------------------

%% Check Inputs
if nargin==0
  % Example
  x = 1:25;
  y = rand(1000,length(x));
  o.AltData = [...
    -0.1+x/25;
    0.4+x/25;
%     1.3-x/25;
    ];
  o.LegendLocation = 'SE';
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

%% default options
oDefault.Bands2Show = [50,70,90]; %[50,60,70,80,90]
oDefault.AltData = [];
oDefault.ShowLegend = 0;
oDefault.LegendLocation = 'Best';
oDefault.LegendOrientation = 'vertical'; %'vertical','horizontal'
oDefault.LegendWithBands = 0;
oDefault.LineColor = vcColorScheme;
% oDefault.LineColor = oDefault.LineColor([3,1,2,4:end],:); 
  % Assumes that 3rd line color is red and first is blue
% oDefault.ShadeColor = [0.2,0.6,0.5];
% oDefault.ShadeColor = [0.45,0.45,0.5]; 
% oDefault.ShadeColor = [0.585,0.585,0.65]; 
oDefault.ShadeColor = [0.72,0.77,0.82]*0.95;
% oDefault.ShadeColor = [0.15,0.25,0.75];
oDefault.ShadeColorBrightness = 0.9;
% oDefault.ShadeFactors = [0.1,0.65]% shade factors at 50 and 90%
% oDefault.ShadeFactors = [0.2,0.7]; % shade factors at 50 and 90%
oDefault.ShadeFactors = [0.1,0.65]; % shade factors at 50 and 90%
% oDefault.ShadeFactors = [0.1,0.7]; % shade factors at 50 and 90%
% MedianColor = [0,0,0.7];
% ShadeColor = [0.2,0.6,0.5];

%% Apply defaults to unspecified options
oList = fieldnames(oDefault);
for jo=1:length(oList)
  oName = oList{jo};
  if ~isfield(o,oName)
    o.(oName) = oDefault.(oName);
  end
end

%% Additional options, if not set before
% o.LineColor = [0.3,0.75,0.3;o.LineColor];
if ~isfield(o,'ShadeColor')
  if ~isfield(o,'ShadeColorBase')
  %   o.ShadeColorBase = [0.70,0.70,0.7]; % grey
  %   o.ShadeColorBase = [0.70,0.80,0.85]; % light blue 
  %   o.ShadeColorBase = [0.60,0.70,0.75]; % grey/Blue
    o.ShadeColorBase = [0.15,0.25,0.75]; % Blue
  end
  if ~isfield(o,'ShadeCLWeight')
    o.ShadeCLWeight = 0.0;
  end
%   o.ShadeColor = [0.65,0.75,0.8]*0.95+0.05*o.LineColor(1,:);
  o.ShadeColor = o.ShadeColorBrightness*o.ShadeColorBase*(1-o.ShadeCLWeight)+...
    o.ShadeCLWeight*o.LineColor(1,:);
end

%% Check x
nx = size(y,2);
if ~exist('x','var')
  x = 1:nx;
end

%% Plot bands
o.Bands2Show = sort(o.Bands2Show,'descend');
nBands = length(o.Bands2Show);
InitHold = ishold;
BandsData = zeros(nBands*2,nx);
BandColorSlope = [-1,1]*o.ShadeFactors'/([1,-1]*o.Bands2Show([1,nBands])');
BandColorCt = o.ShadeFactors(1)-BandColorSlope*o.Bands2Show(nBands);
for jB=1:nBands
  Band = o.Bands2Show(jB);
  BandPath = prctile(y,50+Band/2*[-1,+1]);
%   [Band,BandColorCt+BandColorSlope*Band]
%   BandColor = o.ShadeColor+(1-o.ShadeColor)*(-0.5+(0.7-(-0.5))*Band/100);
  BandColor = o.ShadeColor+...
    (1-o.ShadeColor)*(BandColorCt+BandColorSlope*Band);
  h.Bands(jB) = fill([x,x(end:-1:1)],[BandPath(1,:),BandPath(2,end:-1:1)],...
    BandColor,'EdgeColor',BandColor);
  hold on
  BandsData((jB-1)*2+[1,2],:) = BandPath;
end
h.BandsData = BandsData;

%% Plot zero line
% if o.isZeroLine
%   h.ZeroLine = plot(x,zeros(1,nx),o.ZeroLineStyle,'Color',o.ZeroLineColor,...
%     'LineWidth',o.ZeroLineWidth);
% end

%% Plot lines
YData = [prctile(y,50);o.AltData];
ny = size(YData,1);
if o.ShowLegend
  if~isfield(o,'LegendString')
    o.LegendString{1} = 'Median';
    for j=2:ny
      o.LegendString{j} = sprintf('Alt %.0f',j-1);
    end
    for j=1:nBands*o.LegendWithBands
      o.LegendString{ny+j} = sprintf('%.0f%%',o.Bands2Show(j));
    end
  end
  if o.LegendWithBands
    o.LegendItems = h.Bands;
  end
end
h.Lines = vcPlot(x,YData,o);
o.ShowLight = h.Lines.Options.ShowLight;
h.LegendObj = h.Lines.LegendObj;
if o.ShowLegend && o.ShowLight
  h.LegendLines = h.Lines.LegendLines;
end

%% some more stuff
if ~InitHold
    hold off
end
xlim([x(1),x(end)])


%% Exit
h.XData = x;
h.YData = YData;
h.Options = o;
if nargout==1,varargout = {h};end

%% ------------------------------------------------------------------------
