
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


function processed = process_emggait(reader)
%PROCESS_EMG Process and calculate EMG start, end, mean and median
% frequency and RMS amplitude.
% 
% INPUT:
% 
% reader: structure with read variables
% 
% OUTPUT:
%
% processed: structure with data processed
% 

signal = reader.signal;
signal_f = reader.signal_f;
n_channels = reader.n_channels;
fs = reader.fs;

% this for sitting
% minpdist = 1.*fs;
% [pks_foot, locs_foot] = findpeaks(signal(:,1), 'MinPeakHeight', 0.1,...
%     'MinPeakDistance', minpdist, 'MinPeakProminence', 0.1);

% this for walking
% minpdist = 0.80*fs;
% [pks_foot, locs_foot] = findpeaks(signal(:,1),'MinPeakHeight',0.001,...
%     'MinPeakDistance', minpdist, 'MinPeakProminence',0.1);

[pks_foot, locs_foot_I] = findpeaks(signal_f(:, 1), 'MinPeakHeight', 0.1,...
    'MinPeakProminence', 0.1);

% Custom made block to remove consecutive peaks distant less than x
% seconds. The automatic findpeaks does not work properly for this data
% mindist = 1.5*fs; % sitting-standing
mindist = 0.8*fs;
rm_id = [];
for ip = 2:numel(pks_foot)
    if locs_foot_I(ip) - locs_foot_I(ip-1) < mindist
        rm_id = [rm_id; ip];
    end
end
pks_foot(rm_id) = [];
locs_foot_I(rm_id) = [];

pks_data = cell(1, n_channels);
locs_data = cell(1, n_channels);
locs_data_I = cell(1, n_channels);

% emg window limits, where with I indicates index and without the time
% instant
emg_start_I = cell(1, n_channels);
emg_start = cell(1, n_channels);
emg_end_I = cell(1, n_channels);
emg_end = cell(1, n_channels);
diff_t = cell(1, n_channels);

% peaks location and amplitude, 1 is for footswitch/trigger and others the
% emg channels
pks_data{1} = pks_foot;
locs_data{1} = locs_foot_I./fs;
locs_data_I{1} = locs_foot_I;

emg_start_I{1} = locs_data_I{1} - 0.4*fs;
emg_start{1} = locs_data{1} - 0.4;
emg_end_I{1} = locs_data_I{1} + 0.6*fs;
emg_end{1} = locs_data{1} + 0.6;

% remove negative indices found in emg_start
neg_ind = emg_start_I{1}(:) < 0;
emg_end_I{1}(neg_ind) = [];
locs_data_I{1}(neg_ind) = [];
locs_data{1}(neg_ind) = [];
pks_data{1}(neg_ind) = [];
emg_start{1}(neg_ind) = [];
emg_end{1}(neg_ind) = [];
emg_start_I{1}(neg_ind) = [];

% remove indices bigger than signal length found in emg_end
big_ind = emg_end_I{1}(:) > length(signal_f(:,1));
emg_start_I{1}(big_ind) = [];
locs_data_I{1}(big_ind) = [];
locs_data{1}(big_ind) = [];
pks_data{1}(big_ind) = [];
emg_start{1}(big_ind) = [];
emg_end{1}(big_ind) = [];
emg_end_I{1}(big_ind) = [];

%     [pks_data{id}, locs_data{id}] = findpeaks(-signal(:, id),...
%         'MinPeakDistance', minpdist, 'MinPeakProminence', 10);

amp_med = zeros(1, n_channels);
diff_t_med = zeros(1, n_channels);
drop_peak = zeros(1, n_channels);
n_peaks = zeros(1, n_channels);

% scale is used to compute scaled median absolute deviation to identify
% outliers, check matlab new documentation
scale = -1/(sqrt(2)*erfcinv(3/2));
for id = 2:n_channels    
    for ip = 1:length(emg_start_I{1})
        [pks_data{id}(ip, 1), locs_aux_I] = findpeaks(signal_f(emg_start_I{1}(ip):emg_end_I{1}(ip), id),...
            'SortStr', 'descend', 'NPeaks', 1);

        locs_data{id}(ip, 1) = (emg_start_I{1}(ip) + locs_aux_I)./fs;
        locs_data_I{id}(ip, 1) = emg_start_I{1}(ip) + locs_aux_I;

        emg_start_I{id}(ip, 1) = emg_start_I{1}(ip);
        emg_start{id}(ip, 1) = emg_start{1}(ip);
        emg_end_I{id}(ip, 1) = emg_end_I{1}(ip);
        emg_end{id}(ip, 1) = emg_end{1}(ip);
        
        diff_t{id}(ip, 1) = locs_data{id}(ip, 1) - locs_data{1}(ip);
    end
    
    % check for outliers in peaks (possible problem in emg contact)    
    out_lim = 3*scale*median(abs(pks_data{id}-median(pks_data{id})));
    out_I = find(pks_data{id} > median(pks_data{id}) + out_lim);
    
    if out_I
        locs_data{id}(out_I) = [];
        locs_data_I{id}(out_I) = [];
        emg_start_I{id}(out_I) = [];
        emg_start{id}(out_I) = [];
        emg_end_I{id}(out_I) = [];
        emg_end{id}(out_I) = [];
        pks_data{id}(out_I) = [];
        diff_t{id}(out_I) = [];
    end
    
    n_peaks(id) = length(emg_start_I{id});
    drop_peak(id) = numel(out_I);
    amp_med(id) = median(pks_data{id});
    diff_t_med(id) = median(diff_t{id});
    
end

fmed = cell(1, n_channels);
fmean = cell(1, n_channels);
rms = cell(1, n_channels);

rms_med = zeros(1, n_channels);
fmed_med = zeros(1, n_channels);
fmean_med = zeros(1, n_channels);

for id = 1:n_channels
    for ip = 1:n_peaks(id)
        [fmed{id}(ip), rms{id}(ip), fmean{id}(ip)]= fmed_rms(signal(emg_start_I{id}(ip):emg_end_I{id}(ip)-1, id),...
            fs, emg_end_I{id}(ip)-emg_start_I{id}(ip));
    end    
    rms_med(id) = median(rms{id});
    fmed_med(id) = median(fmed{id});
    fmean_med(id) = median(fmean{id});
    
end

processed.fmed = fmed;
processed.fmean = fmean;
processed.rms = rms;

processed.rms_med = rms_med;
processed.fmed_med = fmed_med;
processed.rms_med = rms_med;
 
processed.emg_start_I = emg_start_I;
processed.emg_end_I = emg_end_I;
processed.emg_start = emg_start;
processed.emg_end = emg_end;

processed.pmax_I = locs_data;
processed.amp_max = pks_data;

end
