
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
% 
answer = questdlg('Process which task?', ...
    'Task selection', 'Gait', 'Sit-Stand', 'Sit-Stand');
switch answer
    case 'Gait'
        task_id = 0;
    case 'Sit-Stand'
        task_id = 1;
end


% loading signal and configuration data
[filename, pathname, fid] = uigetfile({'*.txt','Text files (*.txt)';...
    '*.xlsx;*.xls','MS Excel Files (*.xlsx, *xls)'}, 'Select the signal file',...
    'MultiSelect', 'on');

% data_aux = load([pathname filename]);
% save('cassia_sl.mat','fid', 'filename', 'pathname');

% task_id = 1;
% sub = 4;
% % 
% if sub == 1
%     dat = load('ana.mat');
% elseif sub == 2
%     dat = load('cassia.mat');
% elseif sub == 3
%     dat = load('ana_sl.mat');
% elseif sub == 4
%     dat = load('cassia_sl.mat');
% end
% 
% fid = dat.fid;
% filename = dat.filename;
% pathname = dat.pathname;

data_aux = cell(1, size(filename, 2));
data_len = zeros(1, size(filename, 2));

if fid == 1
    if iscell(filename)
            for nf = 1:size(filename, 2)
                data_aux{nf} = readtable([pathname filename{nf}], 'Delimiter', '\t',...
                    'ReadVariableNames',false);
                data_len(nf) = size(data_aux{nf}, 1);
                minlen = min(data_len);
%             data_aux = readtable([pathname filename{nf}], 'Delimiter', '\t',...
%                 'ReadVariableNames',false, 'Format', '%s%s');
%             data_aux = table2cell(data_aux);
%             data_aux = strrep(data_aux,',','.');
%             data_aux = str2double(data_aux);
            end
            try
            for nf = 1:size(filename, 2)
                if nf == 1
                    data_arr = table2array(data_aux{nf}(1:minlen,1:2));
                    data_emg = data_arr;
                else
                    data_arr = table2array(data_aux{nf}(1:minlen,2:end));
                    data_emg = horzcat(data_emg, data_arr);
                end
            end
            catch
                hwarn = warndlg('Please check if decimal point is a dot and not a comma.',...
                    'Atention!');
                uiwait(hwarn);
        end
    else
        data_emg = readtable([pathname filename], 'Delimiter', '\t',...
            'ReadVariableNames',false);
        data_emg = table2array(data_emg(:,1:end));
    %     data_aux = strrep(data_aux,',','.');
    %     data_aux = str2double(data_aux);
    end
elseif fid == 2
    data_emg = xlsread([pathname filename]);
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

n_channels = size(data_emg,2) - 1;
fs = 1/(data_emg(3,1) - data_emg(2,1));
signal = data_emg(:,2:end);
signal_f = signal;

signal(:, 1) = detrend(signal(:, 1));
% signal(:, 1) = abs(detrend(signal(:, 1)));
signal(:, 2:end) = detrend(signal(:, 2:end));
signal_f(:, 1) = signal(:, 1);
if task_id
    signal_f(:, 1) = signal_f(:, 1) - mean(signal_f(1:fs, 1));
end
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
output_reader.xs = data_emg(:,1);
output_reader.fs = fs;
output_reader.task_id = task_id;

output_reader.fig_titles = fig_titles;

