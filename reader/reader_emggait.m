
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
% Date: 15.04.2018


function output_reader = reader_emggait

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
%             data_aux = readtable([pathname filename{nf}], 'Delimiter', '\t',...
%                 'ReadVariableNames',false, 'Format', '%s%s');  
%             data_aux = table2cell(data_aux);
%             data_aux = strrep(data_aux,',','.');
%             data_aux = str2double(data_aux);
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

% Unfinished attempt to check for filenames and reorganize them
% titles_part = {'Footswitch', 'GluteoD', 'GluteoE', 'IliocostalD',...
%     'IliocostalE', 'LonguissimoD', 'LonguissimoE'};
% 
% for nf = 1:size(filename, 2)
%     for nt = 1:size(titles, 2)
%         id_t = strfind(lower(filename), lower(titles{nt}));
%     if ~find(~cellfun('isempty', id_t))
%         k = strfind(lower(filename{1}), 'ld')
%         fig_titles{nf} = [titles_part num2str(nf)];
%     end
% end

n_channels = size(data,2)-1;
fs = 1/(data(3,1) - data(2,1));
signal = data(:,2:end);
signal_f = signal;

signal(:, 1:end) = detrend(signal(:, 1:end));
signal_f(:, 1) = signal(:, 1);
signal_f(:, 2:end) = lowpass(abs(signal(:, 2:end)), 10, fs);

if iscell(filename)
    fig_titles = cell(1, size(filename, 2));
    for nf = 1:n_channels
        title_aux = strsplit(filename{nf}, '.');
        fig_titles{nf} = title_aux{1};
    end
else
    fig_titles = cell(1, n_channels);
    titles_part = 'DATA';
    for nf = 1:n_channels
        fig_titles{nf} = [titles_part num2str(nf)];
    end
end

output_reader.n_channels = n_channels;
output_reader.signal = signal;
output_reader.signal_f = signal_f;
output_reader.xs = data(:,1);
output_reader.fs = fs;

output_reader.fig_titles = fig_titles;

