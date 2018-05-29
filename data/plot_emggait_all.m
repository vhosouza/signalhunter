function plot_emggait_all(hf, signal_split, n_channels)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

clf;
if n_channels == 11
    set(hf, 'Position', [100. 100. 896. 480.], 'Color', 'w',...
        'Name', 'Split contractions');
    % ha = tight_subplot(2,5,[.01 .03],[.1 .01],[.01 .01]);
    ha = tight_subplot(2, 5, [0.06 .03], [.06 0.06], [0.06 .01]);
    % ha = tight_subplot(2,5);
    
    fig_titles = {'GluteoD', 'GluteoE', 'IliocostalD Pre',...
        'IliocostalD Pos', 'IliocostalE Pre', 'IliocostalE Pos', ...
        'LonguissimoD Pre', 'LonguissimoD Pos','LonguissimoE Pre',...
        'LonguissimoE Pos'};
    
    max_all = @(x) max(max(x));
    yl = 1.1*max(cellfun(max_all, signal_split));
    
    kk = 1;
    for i1 = 1:2:10
        axes(ha(kk));
        avg = mean(signal_split{i1}, 2);
        st = std(signal_split{i1}, [], 2);
        %     hs = subplot(2, 5, kk);
        hold on
        for is = 1:size(signal_split{i1}, 2);
            plot(signal_split{i1}(:, is), 'Color', [153 153 153]/255);
        end
        plot(avg, 'k', 'LineWidth', 1.5);
        plot(avg - st, 'r--', 'LineWidth', 1.5);
        plot(avg + st, 'r--', 'LineWidth', 1.5);
        hold off
        %     ha(kk).YLim = [0 yl];
        ha(kk).Title = text('String', fig_titles{i1});
        
        kk2 = kk + 5;
        axes(ha(kk2));
        avg = mean(signal_split{i1+1}, 2);
        st = std(signal_split{i1+1}, [], 2);
        %     hs = subplot(2, 5, kk);
        hold on
        for is = 1:size(signal_split{i1+1}, 2);
            plot(signal_split{i1+1}(:, is), 'Color', [153 153 153]/255);
        end
        plot(avg, 'k', 'LineWidth', 1.5);
        plot(avg - st, 'r--', 'LineWidth', 1.5);
        plot(avg + st, 'r--', 'LineWidth', 1.5);
        hold off
        %     ha(kk2).YLim = [0 yl];
        ha(kk2).Title = text('String', fig_titles{i1+1});
        set(ha(kk2),'XLabel', text('String', 'Time (a.u.)'));
        
        kk = kk + 1;
    end
    
    set(ha(1),'YLabel', text('String', 'Amplitude (a.u.)'));
    set(ha(6),'YLabel', text('String', 'Amplitude (a.u.)'));
    set(ha(2:5),'YLabel', text('String', ''));
    set(ha(7:end),'YLabel', text('String', ''));
    
    set(ha(:),'XTickLabel', '');
    
    % set(ha(2:5),'YTickLabel', '');
    % set(ha(7:end),'YTickLabel', '');
    
else
    
    set(hf, 'Position', [100. 100. 896. 480.], 'Color', 'w',...
        'Name', 'Split contractions');
    % ha = tight_subplot(2,5,[.01 .03],[.1 .01],[.01 .01]);
    ha = tight_subplot(2, 6, [0.06 .03], [.06 0.06], [0.06 .01]);
    % ha = tight_subplot(2,5);
    
    fig_titles = {'GluteoD Pre', 'GluteoD Pos', 'GluteoE Pre', 'GluteoE Pos',...
        'IliocostalD Pre', 'IliocostalD Pos', 'IliocostalE Pre', 'IliocostalE Pos', ...
        'LonguissimoD Pre', 'LonguissimoD Pos','LonguissimoE Pre',...
        'LonguissimoE Pos'};
    
    max_all = @(x) max(max(x));
    yl = 1.1*max(cellfun(max_all, signal_split));
    
    kk = 1;
    for i1 = 1:2:12
        axes(ha(kk));
        avg = mean(signal_split{i1}, 2);
        st = std(signal_split{i1}, [], 2);
        %     hs = subplot(2, 5, kk);
        hold on
        for is = 1:size(signal_split{i1}, 2);
            plot(signal_split{i1}(:, is), 'Color', [153 153 153]/255);
        end
        plot(avg, 'k', 'LineWidth', 1.5);
        plot(avg - st, 'r--', 'LineWidth', 1.5);
        plot(avg + st, 'r--', 'LineWidth', 1.5);
        hold off
        %     ha(kk).YLim = [0 yl];
        ha(kk).Title = text('String', fig_titles{i1});
        
        kk2 = kk + 6;
        axes(ha(kk2));
        avg = mean(signal_split{i1+1}, 2);
        st = std(signal_split{i1+1}, [], 2);
        %     hs = subplot(2, 5, kk);
        hold on
        for is = 1:size(signal_split{i1+1}, 2);
            plot(signal_split{i1+1}(:, is), 'Color', [153 153 153]/255);
        end
        plot(avg, 'k', 'LineWidth', 1.5);
        plot(avg - st, 'r--', 'LineWidth', 1.5);
        plot(avg + st, 'r--', 'LineWidth', 1.5);
        hold off
        %     ha(kk2).YLim = [0 yl];
        ha(kk2).Title = text('String', fig_titles{i1+1});
        set(ha(kk2),'XLabel', text('String', 'Time (a.u.)'));
        
        kk = kk + 1;
    end
    
    set(ha(1),'YLabel', text('String', 'Amplitude (a.u.)'));
    set(ha(7),'YLabel', text('String', 'Amplitude (a.u.)'));
    set(ha(2:6),'YLabel', text('String', ''));
    set(ha(8:end),'YLabel', text('String', ''));
    
    set(ha(:),'XTickLabel', '');
    
    % set(ha(2:5),'YTickLabel', '');
    % set(ha(7:end),'YTickLabel', '');
end

