
% -------------------------------------------------------------------------
% Signal Hunter - electrophysiological signal analysis  
% Copyright (C) 2013, 2013-2016  University of Sao Paulo
% 
% Homepage:  http://df.ffclrp.usp.br/biomaglab
% Contact:   biomaglab@gmail.com
% License:   GNU - GPL 3 (LICENSE.txt)
% 
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or any later
% version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
% 
% The full GNU General Public License can be accessed at file LICENSE.txt
% or at <http://www.gnu.org/licenses/>.
% 
% -------------------------------------------------------------------------
% 
% Author: Victor Hugo Souza
% Date: 13.11.2016


function [hamp, hlat, hend] = plot_emggait(ax, trigger, signal, xs, amp, lat)

axes(ax);
hold on
plot(xs, signal, 'LineWidth', 2);

% if sum(signal - trigger) ~= 0
%     norm_trigger_off = (trigger/max(trigger))*max(signal) - mean((trigger/max(trigger))*max(signal));
%     plot(xs, norm_trigger_off, 'Color', [153 153 153]/255);
% end

% amp(1) = 0;
% lat(1) = 0;

yl = get(ax, 'YLim');
xl = get(ax, 'XLim');
% % if amp(1) ~=0
%     hpmin = plot(pmin(1), pmin(2), 'xr', 'MarkerSize', 15, 'LineWidth', 2);
for id = 1:size(amp, 1)
    hamp = plot(amp(id, 1), amp(id, 2), 'x', 'MarkerSize', 15,...
        'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560]);
    htrig = plot([trigger(id, 1) trigger(id, 1)], [yl(1) yl(2)],...
        'LineWidth', 2, 'Color', [0. 1. 0.]);
    % % else
    % % %     hpmin = nan;
    % %     hamp = nan;
    % % end
    % %
    % % if lat(1) ~= 0
%     hlat = plot([amp(id, 1) amp(id, 1)], yl, '--',...
%         'MarkerSize', 15, 'LineWidth', 2, 'Color', [0.4660 0.6740 0.1880]);
    
%     hstart = plot([lat(id, 1) lat(id, 1)], yl, 'r--',...
%         'MarkerSize', 15, 'LineWidth', 2);
%     hend = plot([lat(id, 2) lat(id, 2)], yl, 'k--',...
%         'MarkerSize', 15, 'LineWidth', 2);
    
    if mod(id, 2)
        hf = fill([lat(id, 1) lat(id, 2) lat(id, 2) lat(id, 1)], [yl(1) yl(1) yl(2) yl(2)], 'w');
        set(hf, 'FaceColor', [1, 0, 0], 'FaceAlpha', 0.1)
    else
        hf = fill([lat(id, 1) lat(id, 2) lat(id, 2) lat(id, 1)], [yl(1) yl(1) yl(2) yl(2)], 'w');
        set(hf, 'FaceColor', 0.44*[1, 1, 1], 'FaceAlpha', 0.1)
    end
    % % else
    % %     hlat = nan;
    % %     hend = nan;
    % % end
    
    % set(hamp, 'Visible', 'off');
    % set(hlat, 'Visible', 'off');
%     set(hend, 'Visible', 'off');
hlat = 0;
hend = 0;
end

hold off
