function TimeIdx = TimeIdxCreate(TimeStart,TimeEnd)

% TimeIdxCreate
%
% Creates a time index
%
% Usage:
%   TimeIdxCreate(TimeStart,TimeEnd)
%
% Inputs:
%
%   TimeStart (string)
%   Start time period. 
%
%   TimeEnd (string)
%   Last time period.
%
% Output:
%
%   TimeIdx (cell array of strings)
%   Time index.
%
% Convention used:
% Dates in format of '####q#' or '####m##' for quarterly and monthly data,
% respectively. For monthly frequency, for single digit periods they can be show
% up in TimeStart and TimeEnd as either # or ## (e.g. 1 or 01) but the index
% will be in the format ##.
%
% In order to search the time index can use the ismember function. Below a
% couple of examples.
%
% 1. Create TimeIdx, from 1986q3 to 2010q2
%   TimeIdx = TimeIdxCreate('1986q3','2010q2');
%
% 2. Search for the position of 2005q1
%   ismember(TimeIdx,'2005q1')
% 
% 3. Search for the position of 1993q1 through 2010q1
%   ismember(TimeIdx,TimeIdxCreate('1993q1','2010q1'))
% 
% 4. Search for the position of 1993q1 through end of TimeIdx
%   ismember(TimeIdx,TimeIdxCreate('1993q1',TimeIdx{end}))
%
% (These searches give a logic array with 1 where the samples match. If need the
% position values only for the common elements (like 3,4,5,...) then can use the
% find function on top of ismember.)
% 
% See also
% ismember, find
%
% ..............................................................................
% 
% Created: October 18, 2010 by Vasco Curdia
% Updated: November 3, 2010 by Vasco Curdia
% Updated: November 16, 2010 by Ging Cee Ng
% Updated: November 16, 2010 by Vasco Curdia
%          Changed  index to cell array of strings.
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2010-2011 by Vasco Curdia

%% -----------------------------------------------------------------------------

%% Find frequency
if ismember('q',TimeStart)
    nPeriods = 4;
    PerStr = 'q';
    nPerStr = 1;
elseif ismember('m',TimeStart)
    nPeriods = 12;
    PerStr = 'm';
    nPerStr = 2;
else
    error('Frequency could not be detected.')
end

%% Identify limiting dates
StartYear = eval(TimeStart(1:4));
StartPer = eval(TimeStart(6:end));
EndYear = eval(TimeEnd(1:4));
EndPer = eval(TimeEnd(6:end));

%% Create index
TimeIdx={};
for yr=StartYear:EndYear
    for per=1:nPeriods
        if yr==StartYear && per<StartPer
            continue
        elseif yr==EndYear && per>EndPer
            break
        end
        TimeIdx{end+1} = sprintf(['%04.0f%s%0',int2str(nPerStr),'.0f'],yr,PerStr,per);
    end
end

%% -----------------------------------------------------------------------------


