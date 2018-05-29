
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

task_id = reader.task_id;

if task_id
    
    % signal = reader.signal;
    % signal_f = reader.signal_f;
    n_channels = reader.n_channels + 6;
    fs = reader.fs;
    
    % signals from from last 4 channels are divided in two, pre and post the
    % footswitch is activated, visual analysis of two peaks of activation in
    % those time windows
    signal_f = zeros(size(reader.signal_f));
    signal_f(:,1) = reader.signal_f(:,1);
    signal_f(:,2) = reader.signal_f(:,2);
    signal_f(:,3) = reader.signal_f(:,2);
    signal_f(:,4) = reader.signal_f(:,3);
    signal_f(:,5) = reader.signal_f(:,3);
    signal_f(:,6) = reader.signal_f(:,4);
    signal_f(:,7) = reader.signal_f(:,4);
    signal_f(:,8) = reader.signal_f(:,5);
    signal_f(:,9) = reader.signal_f(:,5);
    signal_f(:,10) = reader.signal_f(:,6);
    signal_f(:,11) = reader.signal_f(:,6);
    signal_f(:,12) = reader.signal_f(:,7);
    signal_f(:,13) = reader.signal_f(:,7);
    
    signal = zeros(size(reader.signal));
    signal(:,1) = reader.signal(:,1);
    signal(:,2) = reader.signal(:,2);
    signal(:,3) = reader.signal(:,2);
    signal(:,4) = reader.signal(:,3);
    signal(:,5) = reader.signal(:,3);
    signal(:,6) = reader.signal(:,4);
    signal(:,7) = reader.signal(:,4);
    signal(:,8) = reader.signal(:,5);
    signal(:,9) = reader.signal(:,5);
    signal(:,10) = reader.signal(:,6);
    signal(:,11) = reader.signal(:,6);
    signal(:,12) = reader.signal(:,7);
    signal(:,13) = reader.signal(:,7);
    
    % Custom made block to remove consecutive peaks distant less than x
    % seconds. The automatic findpeaks does not work properly for this data
    % this adjusted for SL
    [pks_foot, locs_foot_I] = findpeaks(signal_f(:, 1), 'MinPeakHeight', 0.1,...
        'MinPeakProminence', 0.1);
    if isempty(pks_foot)
        [~, locs_foot_I] = findpeaks(abs(signal_f(:, 1)), 'MinPeakHeight', 0.1,...
            'MinPeakProminence', 0.1);
        pks_foot = signal_f(locs_foot_I, 1);
    end
    mindist = 1.5*fs; % sit-stand
    
    rm_id = [];
    for ip = 2:numel(pks_foot)
        if locs_foot_I(ip) - locs_foot_I(ip-1) < mindist
            rm_id = [rm_id; ip];
        end
    end
    pks_foot(rm_id) = [];
    locs_foot_I(rm_id) = [];
    
    pks_trigger = cell(1, n_channels);
    locs_trigger = cell(1, n_channels);
    locs_trigger_I = cell(1, n_channels);
    
    % emg window limits, where with I indicates index and without the time
    % instant
    emg_start_I = cell(1, n_channels);
    emg_start = cell(1, n_channels);
    emg_end_I = cell(1, n_channels);
    emg_end = cell(1, n_channels);
    diff_t = cell(1, n_channels-1);
    peak_t = cell(1, n_channels-1);
    peak_I = cell(1, n_channels-1);
    signal_split = cell(1, n_channels-1);
    
    % peaks location and amplitude, 1 is for footswitch/trigger and others the
    % emg channels
    pks_trigger{1} = pks_foot;
    locs_trigger{1} = locs_foot_I./fs;
    locs_trigger_I{1} = locs_foot_I;
    
    emg_start_I{1} = locs_trigger_I{1} - 1.5*fs;
    emg_start{1} = locs_trigger{1} - 1.5;
    emg_end_I{1} = locs_trigger_I{1} + 1.5*fs;
    emg_end{1} = locs_trigger{1} + 1.5;
    
    % remove negative indices found in emg_start
    neg_ind = emg_start_I{1}(:) < 0;
    emg_end_I{1}(neg_ind) = [];
    locs_trigger_I{1}(neg_ind) = [];
    locs_trigger{1}(neg_ind) = [];
    pks_trigger{1}(neg_ind) = [];
    emg_start{1}(neg_ind) = [];
    emg_end{1}(neg_ind) = [];
    emg_start_I{1}(neg_ind) = [];
    
    % remove indices bigger than signal length found in emg_end
    big_ind = emg_end_I{1}(:) > length(signal_f(:,1));
    emg_start_I{1}(big_ind) = [];
    locs_trigger_I{1}(big_ind) = [];
    locs_trigger{1}(big_ind) = [];
    pks_trigger{1}(big_ind) = [];
    emg_start{1}(big_ind) = [];
    emg_end{1}(big_ind) = [];
    emg_end_I{1}(big_ind) = [];
    
    amp_med = zeros(1, n_channels);
    n_peaks = zeros(1, n_channels);
    
    % this applies the emg window start and end to all channels based on the
    % trigger
    emg_start_aux_I = cell(1, n_channels);
    emg_end_aux_I = cell(1, n_channels);
    pulses = cell(1, n_channels);
    for in = 1:n_channels
        emg_start_aux_I{in} = emg_start_I{1};
        emg_end_aux_I{in} = emg_end_I{1};
        pulses{in} = (1:length(emg_start_I{1}))';
    end
    
    % id of the channels with unique emg signal
    ch_id = [2 4 6 8 10 12];
    % scale is used to compute scaled median absolute deviation to identify
    % outliers, check matlab new documentation
    scale = -1/(sqrt(2)*erfcinv(3/2));
    
    for id = ch_id
        for ip = 1:length(emg_start_I{1})
            % other channels have two peaks (pre/post) each step
            % limits the time window between emg_start and the trigger
            % instant for one channel and trigger instant and emg_end for
            % consecutive channel, building a pre-post pair
            [pks_trigger{id}(ip, 1), ~] = findpeaks(signal_f(emg_start_I{1}(ip):locs_trigger_I{1}(ip), id),...
                'SortStr', 'descend', 'NPeaks', 1);
            locs_trigger_I{id}(ip, 1) = locs_trigger_I{1}(ip);
            locs_trigger{id}(ip, 1) = (locs_trigger_I{1}(ip))./fs;
            
            [pks_trigger{id+1}(ip, 1), ~] = findpeaks(signal_f(locs_trigger_I{1}(ip):emg_end_I{1}(ip), id+1),...
                'SortStr', 'descend', 'NPeaks', 1);
            locs_trigger_I{id+1}(ip, 1) = locs_trigger_I{1}(ip);
            locs_trigger{id+1}(ip, 1) = (locs_trigger_I{1}(ip))./fs;
            
            % difference between trigger peak and emg peak in miliseconds
            %         diff_t{id-1}(ip, 1) = 1000.0*(locs_data{id}(ip, 1) - locs_data{1}(ip));
        end
        
        delt = 1.*fs;
        id_2 = id + 1;
        
        out_lim = 10*scale*median(abs(pks_trigger{id}-median(pks_trigger{id})));
        out_I = find(pks_trigger{id} > median(pks_trigger{id}) + out_lim);
        if out_I
            emg_start_aux_I{id}(out_I) = [];
            emg_end_aux_I{id}(out_I) = [];
            pks_trigger{id}(out_I) = [];
            locs_trigger{id}(out_I) = [];
            locs_trigger_I{id}(out_I) = [];
            pulses{id}(out_I) = [];
        end
        
        n_peaks(id) = length(emg_start_aux_I{id});
        %         drop_peak(id) = numel(out_I);
        
        signal_split{id-1} = [];
        signal_split{id-1} = zeros(locs_trigger_I{id}(1)-emg_start_aux_I{id}(1), n_peaks(id));
        
        for ip = 1:n_peaks(id)
            signal_split{id-1}(:, ip) = signal_f(emg_start_aux_I{id}(ip):locs_trigger_I{id}(ip)-1, id);
        end
        signal_avg = mean(signal_split{id-1}, 2);
        [~, locs_avg_I] = findpeaks(signal_avg, 'SortStr', 'descend',...
            'NPeaks', 1);
        
        output = compute_start_end(signal_f(:, id), delt, pulses{id},...
            emg_start_aux_I{id}, locs_avg_I, locs_trigger{id}, locs_trigger_I{id}, fs);
        
        signal_split{id-1} = output.signal_split;
        emg_start_I{id} = output.emg_start_I;
        emg_end_I{id} = output.emg_end_I;
        emg_start{id} = output.emg_start;
        emg_end{id} = output.emg_end;
        pks_trigger{id} = output.peak;
        peak_t{id-1} = output.peak_t;
        peak_I{id-1} = output.peak_I;
        locs_trigger{id} = output.locs_trigger;
        locs_trigger_I{id} = output.locs_trigger_I;
        pulses{id} = output.pulses;
        %         drop_peak(id) = drop_peak(id) + output.ndrops;
        n_peaks(id) = size(signal_split{id-1}, 2);
        
        % start the id_2
            out_lim = 10*scale*median(abs(pks_trigger{id_2}-median(pks_trigger{id_2})));
            out_I = find(pks_trigger{id_2} > median(pks_trigger{id_2}) + out_lim);
            if out_I
                emg_start_aux_I{id_2}(out_I) = [];
                emg_end_aux_I{id_2}(out_I) = [];
                pks_trigger{id_2}(out_I) = [];
                locs_trigger{id_2}(out_I) = [];
                locs_trigger_I{id_2}(out_I) = [];
                pulses{id_2}(out_I) = [];
            end
            
            n_peaks(id_2) = length(emg_start_aux_I{id_2});
            %         drop_peak(id_2) = numel(out_I);
            
            signal_split{id_2-1} = [];
            signal_split{id_2-1} = zeros(emg_end_aux_I{id_2}(1)-locs_trigger_I{id_2}(1), n_peaks(id_2));
            
            for ip = 1:n_peaks(id_2)
                signal_split{id_2-1}(:, ip) = signal_f(locs_trigger_I{id_2}(ip):emg_end_aux_I{id_2}(ip)-1, id_2);
            end
            signal_avg = mean(signal_split{id_2-1}, 2);
            [~, locs_avg_I] = findpeaks(signal_avg, 'SortStr', 'descend',...
                'NPeaks', 1);
            
            output = compute_start_end(signal_f(:, id_2), delt, pulses{id_2},...
                locs_trigger_I{id_2}, locs_avg_I, locs_trigger{id_2}, locs_trigger_I{id_2}, fs);
            
            signal_split{id_2-1} = output.signal_split;
            emg_start_I{id_2} = output.emg_start_I;
            emg_end_I{id_2} = output.emg_end_I;
            emg_start{id_2} = output.emg_start;
            emg_end{id_2} = output.emg_end;
            pks_trigger{id_2} = output.peak;
            peak_t{id_2-1} = output.peak_t;
            peak_I{id_2-1} = output.peak_I;
            locs_trigger{id_2} = output.locs_trigger;
            locs_trigger_I{id_2} = output.locs_trigger_I;
            pulses{id_2} = output.pulses;
            %         drop_peak(id_2) = drop_peak(id_2) + output.ndrops;
            n_peaks(id_2) = size(signal_split{id_2-1}, 2);
                      
    end
    
    fin_pulses = pulses{1};
    for k = 2:length(pulses)
        fin_pulses = intersect(fin_pulses, pulses{k});
    end
    
    n0_peaks = numel(emg_start_I{1});
    peak_drop = zeros(1, n_channels);
    for k = 1:length(pulses)
        [~, ~, fin_inter] = intersect(fin_pulses, pulses{k});
        mask = (1:numel(pulses{k}));
        rm_pulses = setxor(mask, fin_inter);
        emg_start_I{k}(rm_pulses) = [];
        emg_end_I{k}(rm_pulses) = [];
        emg_start{k}(rm_pulses) = [];
        emg_end{k}(rm_pulses) = [];
        pks_trigger{k}(rm_pulses) = [];
        locs_trigger{k}(rm_pulses) = [];
        locs_trigger_I{k}(rm_pulses) = [];
        pulses{k}(rm_pulses) = [];
        n_peaks(k) = numel(pulses{k});
        peak_drop(k) = n0_peaks - n_peaks(k);
        
        if k > 1
            signal_split{k-1}(:, rm_pulses) = [];
            peak_t{k-1}(rm_pulses) = [];
            peak_I{k-1}(rm_pulses) = [];
        end
    end
    
    fmed = cell(1, n_channels-1);
    fmean = cell(1, n_channels-1);
    rms = cell(1, n_channels-1);
    arv = cell(1, n_channels-1);
    
    amp_med = zeros(3, n_channels-1);
    diff_t_med = zeros(3, n_channels-1);
    rms_med = zeros(3, n_channels-1);
    arv_med = zeros(3, n_channels-1);
    fmed_med = zeros(3, n_channels-1);
    fmean_med = zeros(3, n_channels-1);
    n_total = 0;
    
    musc = {};
    musc_med = {};
    musc_code = [1 1 1 1 2 2 2 2 3 3 3 3];
    musc_label = {'GLU', 'ILI', 'LON'};
    
    side = {};
    side_med = {};
    side_code = [1 1 2 2 1 1 2 2 1 1 2 2];
    side_label = ['R', 'L'];
    
    inst = {};
    inst_med = {};
    inst_code = [2 3 2 3 2 3 2 3 2 3 2 3];
    inst_label = {'NA', 'PRE', 'POS'};
    
    for id = 1:length(signal_split)
        % ignore channel 1 that is trigger
        ch = id + 1;
        for ip = 1:n_peaks(ch)
            [fmed{id}(ip, 1), rms{id}(ip, 1), fmean{id}(ip, 1)] = fmed_rms(signal(emg_start_I{ch}(ip):emg_end_I{ch}(ip)-1, ch),...
                fs, emg_end_I{ch}(ip)-emg_start_I{ch}(ip));
            
            diff_t{id}(ip, 1) = locs_trigger_I{ch}(ip) - peak_I{id}(ip);
            arv{id}(ip, 1) = mean(abs(signal(emg_start_I{ch}(ip):emg_end_I{ch}(ip)-1)));
            
        end
        %     m_label = {['M' num2str(id)]};
        m_label = {musc_label{musc_code(id)}};
        musc = [musc; repmat(m_label, [numel(fmed{id}), 1])];
        musc_med = [musc_med; m_label];
        
        s_label = {side_label(side_code(id))};
        side = [side; repmat(s_label, [numel(fmed{id}), 1])];
        side_med = [side_med; s_label];
        
        i_label = {inst_label{inst_code(id)}};
        inst = [inst; repmat(i_label, [numel(fmed{id}), 1])];
        inst_med = [inst_med; i_label];
        
        amp_med(:, id) = prctile(pks_trigger{ch}, [25 50 75]);
        arv_med(:, id) = prctile(arv{id}, [25 50 75]);
        rms_med(:, id) = prctile(rms{id}, [25 50 75]);
        fmed_med(:, id) = prctile(fmed{id}, [25 50 75]);
        fmean_med(:, id) = prctile(fmean{id}, [25 50 75]);
        diff_t_med(:, id) = prctile(diff_t{id}, [25 50 75]);
        n_total = n_total + numel(fmed{id});
        
        %     amp_med(id) = median(pks_trigger{ch});
        %     arv_med(id) = median(arv{id});
        %     rms_med(id) = median(rms{id});
        %     fmed_med(id) = median(fmed{id});
        %     fmean_med(id) = median(fmean{id});
        %     diff_t_med(id) = median(diff_t{id});
    end
    
    % musc = {};
    % musc_med = {};
    % for id = 1:n_channels-1
    %     m_label = {['M' num2str(id)]};
    %     musc = [musc; repmat(m_label, [numel(fmed{id}), 1])];
    %     musc_med = [musc_med; m_label];
    % end
    
    
    emg_start_all = cell(1, 7);
    emg_end_all = cell(1, 7);
    locs_trigger_all = cell(1, 7);
    pks_trigger_all = cell(1, 7);
    
    emg_start_all{1} = emg_start{1};
    emg_end_all{1} = emg_end{1};
    pks_trigger_all{1} = pks_trigger{1};
    locs_trigger_all{1} = locs_trigger{1};
    
    id = 2;
    for ch = 2:2:numel(emg_start)
        emg_start_all{id} = zeros(1*size(emg_start{ch}, 1), 1);
        emg_start_all{id} = emg_start{ch};
        %     emg_start_all{id}(2:2:end) = emg_start{ch+1};
        
        emg_end_all{id} = zeros(1*size(emg_start{ch}, 1), 1);
        %     emg_end_all{id}(1:2:end) = emg_end{ch};
        emg_end_all{id} = emg_end{ch+1};
        
        locs_trigger_all{id} = zeros(1*size(emg_start{ch}, 1), 1);
        locs_trigger_all{id} = locs_trigger{ch};
        %     locs_trigger_all{id}(2:2:end) = locs_trigger{ch+1};
        
        pks_trigger_all{id} = zeros(1*size(emg_start{ch}, 1), 1);
        %     pks_trigger_all{id}(1:2:end) = pks_trigger{ch};
        pks_trigger_all{id} = pks_trigger{ch+1};
        id = id + 1;
    end
    
else
    
    % signal = reader.signal;
    % signal_f = reader.signal_f;
    n_channels = reader.n_channels + 4;
    fs = reader.fs;
    
    % signals from from last 4 channels are divided in two, pre and post the
    % footswitch is activated, visual analysis of two peaks of activation in
    % those time windows
    signal_f = zeros(size(reader.signal_f));
    signal_f(:,1) = reader.signal_f(:,1);
    signal_f(:,2) = reader.signal_f(:,2);
    signal_f(:,3) = reader.signal_f(:,3);
    signal_f(:,4) = reader.signal_f(:,4);
    signal_f(:,5) = reader.signal_f(:,4);
    signal_f(:,6) = reader.signal_f(:,5);
    signal_f(:,7) = reader.signal_f(:,5);
    signal_f(:,8) = reader.signal_f(:,6);
    signal_f(:,9) = reader.signal_f(:,6);
    signal_f(:,10) = reader.signal_f(:,7);
    signal_f(:,11) = reader.signal_f(:,7);
    
    signal = zeros(size(reader.signal));
    signal(:,1) = reader.signal(:,1);
    signal(:,2) = reader.signal(:,2);
    signal(:,3) = reader.signal(:,3);
    signal(:,4) = reader.signal(:,4);
    signal(:,5) = reader.signal(:,4);
    signal(:,6) = reader.signal(:,5);
    signal(:,7) = reader.signal(:,5);
    signal(:,8) = reader.signal(:,6);
    signal(:,9) = reader.signal(:,6);
    signal(:,10) = reader.signal(:,7);
    signal(:,11) = reader.signal(:,7);
    
    % Custom made block to remove consecutive peaks distant less than x
    % seconds. The automatic findpeaks does not work properly for this data
    if task_id
        % this adjusted for SL
        [pks_foot, locs_foot_I] = findpeaks(signal_f(:, 1), 'MinPeakHeight', 0.1,...
            'MinPeakProminence', 0.1);
        if isempty(pks_foot)
            [~, locs_foot_I] = findpeaks(abs(signal_f(:, 1)), 'MinPeakHeight', 0.1,...
                'MinPeakProminence', 0.1);
            pks_foot = signal_f(locs_foot_I, 1);
        end
        mindist = 1.5*fs; % sit-stand
    else
        % this works for gait
        [pks_foot, locs_foot_I] = findpeaks(signal_f(:, 1), 'MinPeakHeight', 0.1,...
            'MinPeakProminence', 0.1);
        mindist = 0.8*fs;
    end
    
    rm_id = [];
    for ip = 2:numel(pks_foot)
        if locs_foot_I(ip) - locs_foot_I(ip-1) < mindist
            rm_id = [rm_id; ip];
        end
    end
    pks_foot(rm_id) = [];
    locs_foot_I(rm_id) = [];
    
    pks_trigger = cell(1, n_channels);
    locs_trigger = cell(1, n_channels);
    locs_trigger_I = cell(1, n_channels);
    
    % emg window limits, where with I indicates index and without the time
    % instant
    emg_start_I = cell(1, n_channels);
    emg_start = cell(1, n_channels);
    emg_end_I = cell(1, n_channels);
    emg_end = cell(1, n_channels);
    diff_t = cell(1, n_channels-1);
    peak_t = cell(1, n_channels-1);
    peak_I = cell(1, n_channels-1);
    signal_split = cell(1, n_channels-1);
    
    % peaks location and amplitude, 1 is for footswitch/trigger and others the
    % emg channels
    pks_trigger{1} = pks_foot;
    locs_trigger{1} = locs_foot_I./fs;
    locs_trigger_I{1} = locs_foot_I;
    
    emg_start_I{1} = locs_trigger_I{1} - 0.5*fs;
    emg_start{1} = locs_trigger{1} - 0.5;
    emg_end_I{1} = locs_trigger_I{1} + 0.5*fs;
    emg_end{1} = locs_trigger{1} + 0.5;
    
    % remove negative indices found in emg_start
    neg_ind = emg_start_I{1}(:) < 0;
    emg_end_I{1}(neg_ind) = [];
    locs_trigger_I{1}(neg_ind) = [];
    locs_trigger{1}(neg_ind) = [];
    pks_trigger{1}(neg_ind) = [];
    emg_start{1}(neg_ind) = [];
    emg_end{1}(neg_ind) = [];
    emg_start_I{1}(neg_ind) = [];
    
    % remove indices bigger than signal length found in emg_end
    big_ind = emg_end_I{1}(:) > length(signal_f(:,1));
    emg_start_I{1}(big_ind) = [];
    locs_trigger_I{1}(big_ind) = [];
    locs_trigger{1}(big_ind) = [];
    pks_trigger{1}(big_ind) = [];
    emg_start{1}(big_ind) = [];
    emg_end{1}(big_ind) = [];
    emg_end_I{1}(big_ind) = [];
    
    amp_med = zeros(1, n_channels);
    n_peaks = zeros(1, n_channels);
    
    % this applies the emg window start and end to all channels based on the
    % trigger
    emg_start_aux_I = cell(1, n_channels);
    emg_end_aux_I = cell(1, n_channels);
    pulses = cell(1, n_channels);
    for in = 1:n_channels
        emg_start_aux_I{in} = emg_start_I{1};
        emg_end_aux_I{in} = emg_end_I{1};
        pulses{in} = (1:length(emg_start_I{1}))';
    end
    
    % id of the channels with unique emg signal
    ch_id = [2 3 4 6 8 10];
    % scale is used to compute scaled median absolute deviation to identify
    % outliers, check matlab new documentation
    scale = -1/(sqrt(2)*erfcinv(3/2));
    
    for id = ch_id
        for ip = 1:length(emg_start_I{1})
            if id < 4
                % first 2 channels have only one peak for each step in gait
                [pks_trigger{id}(ip, 1), ~] = findpeaks(signal_f(emg_start_I{1}(ip):emg_end_I{1}(ip), id),...
                    'SortStr', 'descend', 'NPeaks', 1);
                locs_trigger_I{id}(ip, 1) = locs_trigger_I{1}(ip);
                locs_trigger{id}(ip, 1) = (locs_trigger_I{1}(ip))./fs;
            else
                % other channels have two peaks (pre/post) each step
                % limits the time window between emg_start and the trigger
                % instant for one channel and trigger instant and emg_end for
                % consecutive channel, building a pre-post pair
                [pks_trigger{id}(ip, 1), ~] = findpeaks(signal_f(emg_start_I{1}(ip):locs_trigger_I{1}(ip), id),...
                    'SortStr', 'descend', 'NPeaks', 1);
                locs_trigger_I{id}(ip, 1) = locs_trigger_I{1}(ip);
                locs_trigger{id}(ip, 1) = (locs_trigger_I{1}(ip))./fs;
                
                [pks_trigger{id+1}(ip, 1), ~] = findpeaks(signal_f(locs_trigger_I{1}(ip):emg_end_I{1}(ip), id+1),...
                    'SortStr', 'descend', 'NPeaks', 1);
                locs_trigger_I{id+1}(ip, 1) = locs_trigger_I{1}(ip);
                locs_trigger{id+1}(ip, 1) = (locs_trigger_I{1}(ip))./fs;
            end
            
            % difference between trigger peak and emg peak in miliseconds
            %         diff_t{id-1}(ip, 1) = 1000.0*(locs_data{id}(ip, 1) - locs_data{1}(ip));
        end
        
        if id < 4
            delt = 0.45*fs;
            
            % check for outliers in peaks (possible problem in emg contact)
            out_lim = 10*scale*median(abs(pks_trigger{id}-median(pks_trigger{id})));
            out_I = find(pks_trigger{id} > median(pks_trigger{id}) + out_lim);
            
            if out_I
                emg_start_aux_I{id}(out_I) = [];
                emg_end_aux_I{id}(out_I) = [];
                locs_trigger{id}(out_I) = [];
                locs_trigger_I{id}(out_I) = [];
                pks_trigger{id}(out_I) = [];
                pulses{id}(out_I) = [];
            end
            
            n_peaks(id) = length(emg_start_aux_I{id});
            %         drop_peak(id) = numel(out_I);
            
            % split all steps in the continuous recording into columns of array
            signal_aux = zeros(emg_end_aux_I{id}(1)-emg_start_aux_I{id}(1),...
                n_peaks(id));
            for ip = 1:n_peaks(id)
                signal_aux(:, ip) = signal_f(emg_start_aux_I{id}(ip):emg_end_aux_I{id}(ip)-1, id);
            end
            
            % find the peak of activation in the avereaged signal across steps
            % and in the initial window, next center the emg start and end
            % window relative to the peak
            signal_avg = mean(signal_aux, 2);
            [~, locs_avg_I] = findpeaks(signal_avg, 'SortStr', 'descend',...
                'NPeaks', 1);
            
            % create the array of split signals centered in the peak of
            % averaged signal
            
            output = compute_start_end(signal_f(:, id), delt, pulses{id},...
                emg_start_aux_I{id}, locs_avg_I, locs_trigger{id}, locs_trigger_I{id}, fs);
            
            signal_split{id-1} = output.signal_split;
            emg_start_I{id} = output.emg_start_I;
            emg_end_I{id} = output.emg_end_I;
            emg_start{id} = output.emg_start;
            emg_end{id} = output.emg_end;
            pks_trigger{id} = output.peak;
            peak_t{id-1} = output.peak_t;
            peak_I{id-1} = output.peak_I;
            locs_trigger{id} = output.locs_trigger;
            locs_trigger_I{id} = output.locs_trigger_I;
            pulses{id} = output.pulses;
            
        else
            delt = 0.3*fs;
            id_2 = id + 1;
            
            out_lim = 10*scale*median(abs(pks_trigger{id}-median(pks_trigger{id})));
            out_I = find(pks_trigger{id} > median(pks_trigger{id}) + out_lim);
            if out_I
                emg_start_aux_I{id}(out_I) = [];
                emg_end_aux_I{id}(out_I) = [];
                pks_trigger{id}(out_I) = [];
                locs_trigger{id}(out_I) = [];
                locs_trigger_I{id}(out_I) = [];
                pulses{id}(out_I) = [];
            end
            
            n_peaks(id) = length(emg_start_aux_I{id});
            %         drop_peak(id) = numel(out_I);
            
            signal_split{id-1} = [];
            signal_split{id-1} = zeros(locs_trigger_I{id}(1)-emg_start_aux_I{id}(1), n_peaks(id));
            
            for ip = 1:n_peaks(id)
                signal_split{id-1}(:, ip) = signal_f(emg_start_aux_I{id}(ip):locs_trigger_I{id}(ip)-1, id);
            end
            signal_avg = mean(signal_split{id-1}, 2);
            [~, locs_avg_I] = findpeaks(signal_avg, 'SortStr', 'descend',...
                'NPeaks', 1);
            
            output = compute_start_end(signal_f(:, id), delt, pulses{id},...
                emg_start_aux_I{id}, locs_avg_I, locs_trigger{id}, locs_trigger_I{id}, fs);
            
            signal_split{id-1} = output.signal_split;
            emg_start_I{id} = output.emg_start_I;
            emg_end_I{id} = output.emg_end_I;
            emg_start{id} = output.emg_start;
            emg_end{id} = output.emg_end;
            pks_trigger{id} = output.peak;
            peak_t{id-1} = output.peak_t;
            peak_I{id-1} = output.peak_I;
            locs_trigger{id} = output.locs_trigger;
            locs_trigger_I{id} = output.locs_trigger_I;
            pulses{id} = output.pulses;
            %         drop_peak(id) = drop_peak(id) + output.ndrops;
            n_peaks(id) = size(signal_split{id-1}, 2);
            
            % start the id_2
            out_lim = 10*scale*median(abs(pks_trigger{id_2}-median(pks_trigger{id_2})));
            out_I = find(pks_trigger{id_2} > median(pks_trigger{id_2}) + out_lim);
            if out_I
                emg_start_aux_I{id_2}(out_I) = [];
                emg_end_aux_I{id_2}(out_I) = [];
                pks_trigger{id_2}(out_I) = [];
                locs_trigger{id_2}(out_I) = [];
                locs_trigger_I{id_2}(out_I) = [];
                pulses{id_2}(out_I) = [];
            end
            
            n_peaks(id_2) = length(emg_start_aux_I{id_2});
            %         drop_peak(id_2) = numel(out_I);
            
            signal_split{id_2-1} = [];
            signal_split{id_2-1} = zeros(emg_end_aux_I{id_2}(1)-locs_trigger_I{id_2}(1), n_peaks(id_2));
            
            for ip = 1:n_peaks(id_2)
                signal_split{id_2-1}(:, ip) = signal_f(locs_trigger_I{id_2}(ip):emg_end_aux_I{id_2}(ip)-1, id_2);
            end
            signal_avg = mean(signal_split{id_2-1}, 2);
            [~, locs_avg_I] = findpeaks(signal_avg, 'SortStr', 'descend',...
                'NPeaks', 1);
            
            output = compute_start_end(signal_f(:, id_2), delt, pulses{id_2},...
                locs_trigger_I{id_2}, locs_avg_I, locs_trigger{id_2}, locs_trigger_I{id_2}, fs);
            
            signal_split{id_2-1} = output.signal_split;
            emg_start_I{id_2} = output.emg_start_I;
            emg_end_I{id_2} = output.emg_end_I;
            emg_start{id_2} = output.emg_start;
            emg_end{id_2} = output.emg_end;
            pks_trigger{id_2} = output.peak;
            peak_t{id_2-1} = output.peak_t;
            peak_I{id_2-1} = output.peak_I;
            locs_trigger{id_2} = output.locs_trigger;
            locs_trigger_I{id_2} = output.locs_trigger_I;
            pulses{id_2} = output.pulses;
            %         drop_peak(id_2) = drop_peak(id_2) + output.ndrops;
            n_peaks(id_2) = size(signal_split{id_2-1}, 2);
            
        end
    end
    
    fin_pulses = pulses{1};
    for k = 2:length(pulses)
        fin_pulses = intersect(fin_pulses, pulses{k});
    end
    
    n0_peaks = numel(emg_start_I{1});
    peak_drop = zeros(1, n_channels);
    for k = 1:length(pulses)
        [~, ~, fin_inter] = intersect(fin_pulses, pulses{k});
        mask = (1:numel(pulses{k}));
        rm_pulses = setxor(mask, fin_inter);
        emg_start_I{k}(rm_pulses) = [];
        emg_end_I{k}(rm_pulses) = [];
        emg_start{k}(rm_pulses) = [];
        emg_end{k}(rm_pulses) = [];
        pks_trigger{k}(rm_pulses) = [];
        locs_trigger{k}(rm_pulses) = [];
        locs_trigger_I{k}(rm_pulses) = [];
        pulses{k}(rm_pulses) = [];
        n_peaks(k) = numel(pulses{k});
        peak_drop(k) = n0_peaks - n_peaks(k);
        
        if k > 1
            signal_split{k-1}(:, rm_pulses) = [];
            peak_t{k-1}(rm_pulses) = [];
            peak_I{k-1}(rm_pulses) = [];
        end
    end
    
    fmed = cell(1, n_channels-1);
    fmean = cell(1, n_channels-1);
    rms = cell(1, n_channels-1);
    arv = cell(1, n_channels-1);
    
    amp_med = zeros(3, n_channels-1);
    diff_t_med = zeros(3, n_channels-1);
    rms_med = zeros(3, n_channels-1);
    arv_med = zeros(3, n_channels-1);
    fmed_med = zeros(3, n_channels-1);
    fmean_med = zeros(3, n_channels-1);
    n_total = 0;
    
    musc = {};
    musc_med = {};
    musc_code = [1 1 2 2 2 2 3 3 3 3];
    musc_label = {'GLU', 'ILI', 'LON'};
    
    side = {};
    side_med = {};
    side_code = [1 2 1 1 2 2 1 1 2 2];
    side_label = ['R', 'L'];
    
    inst = {};
    inst_med = {};
    inst_code = [1 1 2 3 2 3 2 3 2 3];
    inst_label = {'NA', 'PRE', 'POS'};
    
    for id = 1:length(signal_split)
        % ignore channel 1 that is trigger
        ch = id + 1;
        for ip = 1:n_peaks(ch)
            [fmed{id}(ip, 1), rms{id}(ip, 1), fmean{id}(ip, 1)] = fmed_rms(signal(emg_start_I{ch}(ip):emg_end_I{ch}(ip)-1, ch),...
                fs, emg_end_I{ch}(ip)-emg_start_I{ch}(ip));
            
            diff_t{id}(ip, 1) = locs_trigger_I{ch}(ip) - peak_I{id}(ip);
            arv{id}(ip, 1) = mean(abs(signal(emg_start_I{ch}(ip):emg_end_I{ch}(ip)-1)));
            
        end
        %     m_label = {['M' num2str(id)]};
        m_label = {musc_label{musc_code(id)}};
        musc = [musc; repmat(m_label, [numel(fmed{id}), 1])];
        musc_med = [musc_med; m_label];
        
        s_label = {side_label(side_code(id))};
        side = [side; repmat(s_label, [numel(fmed{id}), 1])];
        side_med = [side_med; s_label];
        
        i_label = {inst_label{inst_code(id)}};
        inst = [inst; repmat(i_label, [numel(fmed{id}), 1])];
        inst_med = [inst_med; i_label];
        
        amp_med(:, id) = prctile(pks_trigger{ch}, [25 50 75]);
        arv_med(:, id) = prctile(arv{id}, [25 50 75]);
        rms_med(:, id) = prctile(rms{id}, [25 50 75]);
        fmed_med(:, id) = prctile(fmed{id}, [25 50 75]);
        fmean_med(:, id) = prctile(fmean{id}, [25 50 75]);
        diff_t_med(:, id) = prctile(diff_t{id}, [25 50 75]);
        n_total = n_total + numel(fmed{id});
        
        %     amp_med(id) = median(pks_trigger{ch});
        %     arv_med(id) = median(arv{id});
        %     rms_med(id) = median(rms{id});
        %     fmed_med(id) = median(fmed{id});
        %     fmean_med(id) = median(fmean{id});
        %     diff_t_med(id) = median(diff_t{id});
    end
    
    % musc = {};
    % musc_med = {};
    % for id = 1:n_channels-1
    %     m_label = {['M' num2str(id)]};
    %     musc = [musc; repmat(m_label, [numel(fmed{id}), 1])];
    %     musc_med = [musc_med; m_label];
    % end
    
    
    emg_start_all = cell(1, 7);
    emg_end_all = cell(1, 7);
    locs_trigger_all = cell(1, 7);
    pks_trigger_all = cell(1, 7);
    
    for id = 1:3
        emg_start_all{id} = emg_start{id};
        emg_end_all{id} = emg_end{id};
        pks_trigger_all{id} = pks_trigger{id};
        locs_trigger_all{id} = locs_trigger{id};
    end
    
    id = 4;
    for ch = 4:2:numel(emg_start)
        emg_start_all{id} = zeros(1*size(emg_start{ch}, 1), 1);
        emg_start_all{id} = emg_start{ch};
        %     emg_start_all{id}(2:2:end) = emg_start{ch+1};
        
        emg_end_all{id} = zeros(1*size(emg_start{ch}, 1), 1);
        %     emg_end_all{id}(1:2:end) = emg_end{ch};
        emg_end_all{id} = emg_end{ch+1};
        
        locs_trigger_all{id} = zeros(1*size(emg_start{ch}, 1), 1);
        locs_trigger_all{id} = locs_trigger{ch};
        %     locs_trigger_all{id}(2:2:end) = locs_trigger{ch+1};
        
        pks_trigger_all{id} = zeros(1*size(emg_start{ch}, 1), 1);
        %     pks_trigger_all{id}(1:2:end) = pks_trigger{ch};
        pks_trigger_all{id} = pks_trigger{ch+1};
        id = id + 1;
    end
    
end
processed.fmed = fmed;
processed.fmean = fmean;
processed.rms = rms;
processed.arv = arv;
processed.diff_t = diff_t;

processed.rms_med = rms_med;
processed.arv_med = arv_med;
processed.fmed_med = fmed_med;
processed.fmean_med = fmean_med;
processed.diff_t_med = diff_t_med;

processed.emg_start_I = emg_start_I;
processed.emg_end_I = emg_end_I;

processed.emg_start = emg_start_all;
processed.emg_end = emg_end_all;
processed.pmax_I = locs_trigger_all;
processed.amp_max = pks_trigger_all;

% processed.emg_start = emg_start;
% processed.emg_end = emg_end;
% processed.pmax_I = locs_trigger;
% processed.amp_max = pks_trigger;

processed.n_total = n_total;
processed.musc = musc;
processed.musc_med = musc_med;
processed.side = side;
processed.side_med = side_med;
processed.inst = inst;
processed.inst_med = inst_med;

if n_channels == 11 || n_channels == 13
    hf_all = figure(2);
    plot_emggait_all(hf_all, signal_split, n_channels)
else
    
end

end

function output = compute_start_end(signal, delt, pulses, start_I, peak_ref_I, locs_trigger, locs_trigger_I, fs)

n_peaks = numel(pulses);

signal_split = zeros(2*delt, n_peaks);
emg_start_I = zeros(n_peaks, 1);
emg_end_I = zeros(n_peaks, 1);
emg_start = zeros(n_peaks, 1);
emg_end = zeros(n_peaks, 1);
peak = zeros(n_peaks, 1);
peak_t = zeros(n_peaks, 1);
peak_I = zeros(n_peaks, 1);
ndrops = 0;

for ip = 1:n_peaks
    new_start_I = start_I(ip) + (peak_ref_I - delt);
    new_end_I = start_I(ip) + (peak_ref_I + delt);
    
    % check if window does not exceed array limits
    if new_start_I >= 1 && new_end_I <= size(signal, 1)
        
        [peak(ip), peak_aux] = findpeaks(signal(new_start_I:new_end_I-1),...
            'SortStr', 'descend', 'NPeaks', 1);
        
        signal_split(:, ip) = signal(new_start_I:new_end_I-1);
        peak_I(ip) = new_start_I + peak_aux;
        peak_t(ip) = (new_start_I + peak_aux)./fs;
        
        emg_start_I(ip) = new_start_I;
        emg_start(ip) = new_start_I/fs;
        emg_end_I(ip) = new_end_I;
        emg_end(ip) = new_end_I/fs;
        
    end
end

if ~isempty(find(emg_start_I == 0, 1))
    ind = find(emg_start_I == 0);
    signal_split(:, ind) = [];
    emg_start_I(ind) = [];
    emg_end_I(ind) = [];
    emg_start(ind) = [];
    emg_end(ind) = [];
    peak(ind) = [];
    peak_t(ind) = [];
    peak_I(ind) = [];
    locs_trigger(ind) = [];
    locs_trigger_I(ind) = [];
    pulses(ind) = [];
    ndrops = numel(ind);
end

output = [];
output.signal_split = signal_split;
output.emg_start_I = emg_start_I;
output.emg_end_I = emg_end_I;
output.emg_start = emg_start;
output.emg_end = emg_end;
output.peak = peak;
output.peak_t = peak_t;
output.peak_I = peak_I;
output.locs_trigger = locs_trigger;
output.locs_trigger_I = locs_trigger_I;
output.pulses = pulses;
output.ndrops = ndrops;

end