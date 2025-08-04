function ylimfac(varargin)
%ylimfac(_) function that nicely handles the y-limits of the axes
% usages:
%       ylimfac         : Applies on the current axis and puts some extra 
%                         whitespace of 10% of ymax-ymin to the axes
%       ylimfac(ax)     : Applies 10% on the axis specified by ax
%       ylimfac(ax,fac) : Applies 100*fac% on the axis specified by ax if
%                         ax is [], ax = gca;
%

narginchk(0,2)
ax = gca;
fac = 0.1;

for ii = 1:nargin
    switch nargin
        case 1
            if isa(varargin{1},'matlab.graphics.axis.Axes')
                ax = varargin{1};
            else 
                ax = gca; 
            end
        case 2
            fac = varargin{2};
    end
end

% find largest y values
maxy = -Inf;
miny = Inf;
objs = findobj(ax);
for ii = 1:numel(objs)
    if strcmp(objs(ii).Type,'line')
        if max(objs(ii).YData) > maxy
            maxy = max(objs(ii).YData);
        end
        if min(objs(ii).YData) < miny
            miny = min(objs(ii).YData);
        end
    end
    if strcmp(objs(ii).Type,'constantline')
        if objs(ii).InterceptAxis == 'y'
            if max(objs(ii).Value) > maxy
                maxy = max(objs(ii).Value);
            end
            if min(objs(ii).Value) < miny
                miny = min(objs(ii).Value);
            end
        end
    end
end
dy = fac*(maxy-miny);
ylim(ax,[miny-dy, maxy+dy])
end
