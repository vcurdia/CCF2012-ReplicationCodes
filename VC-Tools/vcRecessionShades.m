function vcRecessionShades(Rec,varargin)

% vcRecessionShades
%
% Adds recession shades to a plot. 
%
% Note: This function should be used when the plot is completely finished and no 
% further axis resizing will take place.
%
% Usage:
%
%   vcRecessionShades(Rec)
%   Where Rec is a series with 1 at the positions of the recessions, and 0
%   otherwise.
%
%   vcRecessionShades(...,'TimeIdx',TimeIdx,...)
%   Sets the values time index in the plot. Default: 1,2,3,...
%
%   vcRecessionShades(...,'EdgeColor',value,...)
%   Sets the color of the edge of the shades to value. Default: 'none'
%
%   vcRecessionShades(...,'FaceColor',value,...)
%   Sets the color of the face of the shades to value. Default: .85*[1,1,1]
%
%   vcRecessionShades(...,'LayerPos',value,...)
%   If value is set to 'top' then the axis lines show up on top of the
%   figure. If instead value is set to 'bottom' then the axis lines show up
%   below everything in the figure, including the shades. Default: 'top'
%
% ..............................................................................
% 
% Created: April 26, 2011 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2011 by Vasco Curdia

%% -----------------------------------------------------------------------------

%% Defaults
EdgeColor = 'none';
FaceColor = .85*[1,1,1];
LayerPos = 'top';
TimeIdx = 1:length(Rec);
ShowBaseline = 'off';

%% Check options
if ~isempty(varargin)
    nOptions = length(varargin);
    if mod(nOptions,2), error('Incorrect number of optional arguments.'), end
    for jO=1:nOptions/2
        eval(sprintf('%s = varargin{%.0f};',varargin{(jO-1)*2+1},jO*2))
    end
end

%% apply shades
yBounds = get(gca,'YLim');
xBounds = get(gca,'XLim');
hold on
h1 = bar(TimeIdx,yBounds(2)*Rec,1,...
    'EdgeColor',EdgeColor,'FaceColor',FaceColor,'ShowBaseline',ShowBaseline);
h2 = bar(TimeIdx,yBounds(1)*Rec,1,...
    'EdgeColor',EdgeColor,'FaceColor',FaceColor,'ShowBaseline',ShowBaseline);
uistack(h1,'bottom')
uistack(h2,'bottom')
hold off
ylim(yBounds)
xlim(xBounds)
set(gca,'Layer',LayerPos)

%% -----------------------------------------------------------------------------

