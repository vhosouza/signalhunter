
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


function output_reader = reader_emganalysis

% loading signal and configuration data
[filename, pathname, fid] = uigetfile({'*.txt','Text files (*.txt)';...
    '*.xlsx;*.xls','MS Excel Files (*.xlsx, *xls)'}, 'Select the signal file',...
    'MultiSelect', 'on');
% data_aux = load([pathname filename]);

if fid == 1
    if iscell(filename)
        for nf = 1:size(filename, 2)
            data_aux = readtable([pathname filename{nf}], 'Delimiter', '\t',...
                'ReadVariableNames',false);
            if nf == 1
                data_aux = table2array(data_aux(:,1:2));
                data = data_aux;
            else
                data_aux = table2array(data_aux(:,2:end));
                data = horzcat(data, data_aux);
            end
        end
    else
        data = readtable([pathname filename], 'Delimiter', '\t',...
            'ReadVariableNames',false);
        data = table2array(data(:,1:end));
    %     data_aux = strrep(data_aux,',','.');
    %     data_aux = str2double(data_aux);
    end
elseif fid == 2
    data = xlsread([pathname filename]);
%     data_aux = strrep(data_aux,',','.');
%     data_aux = str2double(data_aux);
end

titles_part = 'DATA';
fig_titles = cell(1, size(data, 2)-1);
for nf = 1:size(filename, 2)
    fig_titles{nf} = [titles_part num2str(nf)];
end

fs = 1/(data(3,1) - data(2,1));
signal = data(:,2:end);
baseline = signal(1:1*round(fs),:);
offset = repmat(mean(baseline), [size(signal,1),1]);

output_reader.n_channels = size(data,2)-1;
output_reader.signal = signal - offset;
output_reader.xs = data(:,1);
output_reader.force = data(:,2);
output_reader.baseline = data(1:1*round(fs),2:end);
output_reader.fs = fs;
output_reader.offset = mean(output_reader.baseline);

output_reader.fig_titles = fig_titles;
% output_reader.fig_titles = {'FORCE', 'EMG1', 'EMG2', 'EMG3', 'EMG4'};

