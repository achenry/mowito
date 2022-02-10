% clear all; close all; clc;
clc; close all;

addpath('data_functions')

gs_nom = 480; % Region 3 rated generator speed [rpm]

% Juan experiments are cold-started with Delta Beta set to 0

if ~exist('./data', 'dir')
    mkdir('./data');
end

if ~exist('./results', 'dir')
    mkdir('./results');
end

switch_names = {'dat_N10_FB', 'dat_N15_FB', 'dat_N20_FB', 'dat_N10_FF', 'dat_N15_FF', 'dat_N20_FF'};
wpass = 0.005;
for switch_number = [1, 2, 3, 4, 5, 6]
    filepath = ['./data/', switch_names{switch_number}, '.mat'];
    if exist(filepath, 'file')
        load(filepath);
    else
        disp(['Running Auswertung for ', switch_names{switch_number}])
        Auswertung_v2;
        data_struct = unpack_Auswertung_data_mpc(t1, data);
        
        if switch_number == 1
            % fetch first and last value of wind speed
            
            data_struct.HW_wind_speed_filt ...
                = lowpass(data_struct.HW_wind_speed, wpass);

            idx_1 = 100;
            idx_2 = length(data_struct.HW_wind_speed(1:end - 100));

            hw_ws_1 = data_struct.HW_wind_speed_filt(idx_1);
            hw_ws_2 = data_struct.HW_wind_speed_filt(idx_2);

            data_struct.HW_wind_speed_filt ...
                = data_struct.HW_wind_speed_filt(idx_1:idx_2, :);
        else
            % filter wind ts, get t1 and t2
            nan_ws_idx = isnan(data_struct.HW_wind_speed);
            t = 1:numel(data_struct.HW_wind_speed);
            data_struct.HW_wind_speed(nan_ws_idx) ... 
                 = interp1(t(~nan_ws_idx), ...
                 data_struct.HW_wind_speed(~nan_ws_idx), ...
                 t(nan_ws_idx));
            nonan_ws_idx = find(~isnan(data_struct.HW_wind_speed));

            idx_1_nonan = nonan_ws_idx(1);
            idx_2_nonan = nonan_ws_idx(end);
            
            data_struct.HW_wind_speed_filt ...
                = lowpass(data_struct.HW_wind_speed(idx_1_nonan:idx_2_nonan, :), wpass);

            hw_ws_1_idx = find(data_struct.HW_wind_speed_filt >= hw_ws_1, 1, 'first');
            hw_ws_2_idx = find(data_struct.HW_wind_speed_filt <= hw_ws_2); % get last value that is less than from end
            hw_ws_2_idx = hw_ws_2_idx(hw_ws_2_idx > (length(data_struct.HW_wind_speed_filt) - 3000));
            hw_ws_2_idx = hw_ws_2_idx(1);

            idx_1 = idx_1_nonan - 1 + hw_ws_1_idx;
            idx_2 = idx_1_nonan - 1 + hw_ws_2_idx;

            data_struct.HW_wind_speed_filt ...
                = data_struct.HW_wind_speed_filt(hw_ws_1_idx:hw_ws_2_idx, :);

        end

        data_struct.time = data_struct.time(idx_1:idx_2, :);
        data_struct.time = data_struct.time - data_struct.time(1);
        data_struct.pitch_command = data_struct.pitch_command(idx_1:idx_2, :);
        data_struct.HW_wind_speed = data_struct.HW_wind_speed(idx_1:idx_2, :);
        data_struct.feedback_only = data_struct.feedback_only(idx_1:idx_2, :);
        data_struct.solve_status = data_struct.solve_status(idx_1:idx_2, :);
        data_struct.solve_time_ms = data_struct.solve_time_ms(idx_1:idx_2, :);
        data_struct.iterations = data_struct.iterations(idx_1:idx_2, :);
        
        data_struct.rpm = resample(RPM_res_f, 100, 5000);
        data_struct.rpm = data_struct.rpm(idx_1:idx_2, :);

        %figure(3)
%         subplot(3, 1, 1);
%         plot(1:length(RPM_res_f), RPM_res_f);
%         subplot(3, 1, 2);
%         plot(1:length(data_struct.rpm), data_struct.rpm);
%         subplot(3, 1, 3);
        %plot(1:length(data_struct.rpm), data_struct.rpm);

        figure(2);
        subplot(length(switch_names), 3, (switch_number * 3) - 2);
        plot(data_struct.time, data_struct.HW_wind_speed_filt);
        if switch_number == 1
            title('Filtered Horizontal Wind_Speed (m/s)');
        end

        if mod(switch_number, 3) == 1
            ylb = ylabel(switch_names(switch_number), 'Rotation', 0, 'Interpreter', 'none');
            ylb.Position(1) = ylb.Position(1) - abs(ylb.Position(1) * 0.1);
        end

        hold on;
        subplot(length(switch_names), 3, (switch_number * 3) - 1);
        plot(data_struct.time, data_struct.pitch_command);

        if switch_number == 1
            title('Blpitch Command (deg)');
        end

        subplot(length(switch_names), 3, (switch_number * 3));
        plot(data_struct.time, data_struct.rpm);

        if switch_number == 1
            title(['Generator Speed (rpm)']);
        end

        data_struct.rpm_error = data_struct.rpm - gs_nom;
        data_struct.cum_rpm_error = cumsum(data_struct.rpm_error);
        data_struct.cum_blpitch_travel = cumsum(abs(data_struct.pitch_command));
        assignin('base', switch_names{switch_number}, struct(data_struct));
        save(filepath, switch_names{switch_number});
    end
end
savefig('./filtered_ts.fig');

%% Select Data
data = struct;
data.horizon_lengths = [10, 15, 20];
data.FB.dataset = {dat_N10_FB, dat_N15_FB, dat_N20_FB};
data.FF.dataset = {dat_N10_FF, dat_N15_FF, dat_N20_FF};

data_labels = ['FB', 'FF'];

%% Plot cumulative generator speed tracking error for different horizon lengths

% tf.FB = min([length(data.FB{1}.time), length(data.FB{2}.time), length(data.FB{3}.time)]);
% tf.FF = min([length(data.FF{1}.time), length(data.FF{2}.rpm), length(data.FF{3}.time)]);

data.FB.cumulative_rpm_error = [];
data.FF.cumulative_rpm_error = [];

data.FB.cumulative_pitch_travel = [];
data.FF.cumulative_pitch_travel = [];

for d = 1:length(data.horizon_lengths)
    data.FB.cumulative_rpm_error = [data.FB.cumulative_rpm_error, nonan(data.FB.dataset{d}.cum_rpm_error(end))];
    data.FF.cumulative_rpm_error = [data.FF.cumulative_rpm_error, nonan(data.FF.dataset{d}.cum_rpm_error(end))];

    data.FB.cumulative_pitch_travel = [data.FB.cumulative_pitch_travel, nonan(data.FB.dataset{d}.cum_blpitch_travel(end))];
    data.FF.cumulative_pitch_travel = [data.FF.cumulative_pitch_travel, sum(abs(nonan(data.FF.dataset{d}.cum_blpitch_travel(end))))];
end

figure
ax1 = subplot(2, 1, 1);
plot(data.horizon_lengths, data.FB.cumulative_rpm_error)
subtitle('FB')

ax2 = subplot(2, 1, 2);
plot(data.horizon_lengths, data.FF.cumulative_rpm_error)
subtitle('FF')

xlabel('$N$', "Interpreter","latex")
ylabel('$\sum_k |\omega_k - \omega_{0}|$', "Interpreter","latex")
sgtitle('Cumulative generator speed tracking error vs N')
linkaxes([ax1, ax2], 'y')
savefig(['./results/', 'CumulativeDeltaOmega_vs_N']);

figure
ax1 = subplot(2, 1, 1);
plot(data.horizon_lengths, data.FB.cumulative_pitch_travel)
subtitle('FB')

ax2 = subplot(2, 1, 2);
plot(data.horizon_lengths, data.FF.cumulative_pitch_travel)
subtitle('FF')

xlabel('$N$', "Interpreter","latex")
ylabel('$\sum_k |\Delta \beta_k|$', "Interpreter","latex")
sgtitle('Cumulative pitch travel vs N')
linkaxes([ax1, ax2], 'y')
savefig(['./results/', 'CumulativePitchTravel_vs_N']);

%% Plot blade pitch command vs time for different horizon lengths
colors = ['b', 'r', 'g'];

figure
axes = [];
for d = 1:length(data.horizon_lengths)
    axes = [axes subplot(2, 3, d)];
    plot(data.FB.dataset{d}.time, data.FB.dataset{d}.cum_rpm_error, colors(d))
    subtitle(['N=', num2str(data.horizon_lengths(d)), ', FB'], 'Interpreter', 'latex');
    axes = [axes subplot(2, 3, 3 + d)];
    plot(data.FF.dataset{d}.time, data.FF.dataset{d}.cum_rpm_error, colors(d))
    subtitle(['N=', num2str(data.horizon_lengths(d)), ', FF'], 'Interpreter', 'latex');
end

linkaxes(axes, 'y');
sgtitle('Cumulative generator speed error vs time for different horizon lengths')
xlabel('$t$', "Interpreter","latex")
ylabel('$\beta(t)$', "Interpreter","latex")
savefig(['./results/', 'CumulativeDeltaOmega_vs_t']);


figure
axes = [];
for d = 1:length(data.horizon_lengths)
    axes = [axes subplot(2, 3, d)];
    plot(data.FB.dataset{d}.time, cumsum(data.FB.dataset{d}.pitch_command), colors(d))
    subtitle(['N=', num2str(data.horizon_lengths(d)), ', FB'], 'Interpreter', 'latex')
    axes = [axes subplot(2, 3, 3 + d)];
    plot(data.FF.dataset{d}.time, cumsum(data.FF.dataset{d}.pitch_command), colors(d))
    subtitle(['N=', num2str(data.horizon_lengths(d)), ', FF'], 'Interpreter', 'latex')
end

linkaxes(axes, 'y');
sgtitle('Cumulative pitch travel vs time for different horizon lengths')
xlabel('$t$', "Interpreter","latex")
ylabel('$\beta(t)$', "Interpreter","latex")
savefig(['./results/', 'CumulativePitchTravel_vs_t']);

%% Plot solve status for different horizon lengths
data.FB.solve_status = [];
data.FF.solve_status = [];

for d = 1:length(data.horizon_lengths)
    data.FB.solve_status = [data.FB.solve_status, sum(data.FB.dataset{d}.solve_status == 1) / length(data.FB.dataset{d}.solve_status)];
    data.FF.solve_status = [data.FF.solve_status, sum(data.FF.dataset{d}.solve_status == 1) / length(data.FF.dataset{d}.solve_status)];
end

figure
axes = [];

axes = [axes subplot(2, 1, 1)];
plot(data.horizon_lengths, data.FB.solve_status);
subtitle('FB');

axes = [axes subplot(2, 1, 2)];
plot(data.horizon_lengths, data.FF.solve_status);
subtitle('FF');

xlabel('$N$', "Interpreter","latex")
ylabel('$\%$', "Interpreter","latex")
sgtitle('Percentage of solved time instances vs N');
linkaxes(axes, 'y');
savefig(['./results/', 'PercentSolved_vs_N']);

%% Plot mean, min, max n-iters vs time for different horizon lengths

data.FB.mean_iters = [];
data.FF.mean_iters = [];

data.FB.min_iters = [];
data.FF.min_iters = [];

data.FB.max_iters = [];
data.FF.max_iters = [];

for d = 1:length(data.horizon_lengths)
    data.FB.mean_iters = [data.FB.mean_iters, mean(data.FB.dataset{d}.iterations(~isnan(data.FB.dataset{d}.iterations)))];
    data.FF.mean_iters = [data.FF.mean_iters, mean(data.FF.dataset{d}.iterations(~isnan(data.FF.dataset{d}.iterations)))];

    data.FB.max_iters = [data.FB.max_iters, max(data.FB.dataset{d}.iterations(~isnan(data.FB.dataset{d}.iterations)))];
    data.FF.max_iters = [data.FF.max_iters, max(data.FF.dataset{d}.iterations(~isnan(data.FF.dataset{d}.iterations)))];

    data.FB.min_iters = [data.FB.min_iters, min(data.FB.dataset{d}.iterations(~isnan(data.FB.dataset{d}.iterations)))];
    data.FF.min_iters = [data.FF.min_iters, min(data.FF.dataset{d}.iterations(~isnan(data.FF.dataset{d}.iterations)))];
end

figure
subplot(2, 1, 1)
plot(data.horizon_lengths, data.FB.mean_iters)
subtitle('Mean Num. Iterations, FB');

subplot(2, 1, 2)
plot(data.horizon_lengths, data.FF.mean_iters)
subtitle('Mean Num. Iterations, FF');

xlabel('$N$', "Interpreter","latex")
sgtitle('Mean num iterations vs time for different horizon lengths')
savefig(['./results/', 'MeanIter_vs_t']);

% subplot(3, 1, 2)
% plot(data.horizon_lengths, max_iters)
% subtitle('Max Num. Iterations', "Interpreter","latex")
% 
% subplot(3, 1, 3)
% plot(data.horizon_lengths, min_iters)
% subtitle('Min Num. Iterations', "Interpreter","latex")


%% Plot cost at each time-step x'Qx + u'Ru TODO

data.FB.stage_cost = [];
data.FF.stage_cost = [];

Q = [1 0; 0 1];
R = 1e5;

axes1 = []; axes2 = []; axes3 = [];

for d = 1:length(data.horizon_lengths)

    % time
    t_FB = data.FB.dataset{d}.time;
    t_FF = data.FF.dataset{d}.time;

    % output states
    y_FB = [data.FB.dataset{d}.rpm_error, data.FB.dataset{d}.cum_rpm_error];
    y_FF = [data.FF.dataset{d}.rpm_error, data.FF.dataset{d}.cum_rpm_error];

    u_FB = data.FB.dataset{d}.pitch_command;
    u_FF = data.FF.dataset{d}.pitch_command;

    sc_FB = [];
    sc_FF = [];
    % TODO cost of each term, y1, y2, u
    % TODO run MPC with real data, recording num iterations
    % TODO set ylim to equal
    for k = 1:length(y_FB)
        sc_FB = [sc_FB; [y_FB(k, 1)^2 * Q(1, 1), y_FB(k, 2)^2 * Q(2, 2), u_FB(k)^2 * R]];
    end

    for k = 1:length(y_FF)
%         sc_FF = [sc_FF; y_FF(k, :) * Q * y_FF(k, :)' + u_FF(k) * R * u_FF(k)];
            sc_FF = [sc_FF; [y_FF(k, 1)^2 * Q(1, 1), y_FF(k, 2)^2 * Q(2, 2), u_FF(k)^2 * R]];
    end

    figure(8) % RotSpeed Error Cost
    axes1 = [axes1 subplot(2, 3, d)];
    plot(t_FB, sc_FB(:, 1));
    subtitle(['N=', num2str(data.horizon_lengths(d)), ', FB']);
    axes1 = [axes1 subplot(2, 3, 3 + d)];
    plot(t_FF, sc_FF(:, 1));
    subtitle(['N=', num2str(data.horizon_lengths(d)), ', FF']);
    xlabel('$t$', "Interpreter","latex")
    % ylabel('$$', "Interpreter","latex")
    sgtitle('RotSpeed Error Stage Cost vs time for different horizon lengths');
    
    figure(9) % Cumulative RotSpeed Error Cost
    axes2 = [axes2 subplot(2, 3, d)];
    plot(t_FB, sc_FB(:, 2));
    subtitle(['N=', num2str(data.horizon_lengths(d)), ', FB']);
    axes2 = [axes2 subplot(2, 3, 3 + d)];
    plot(t_FF, sc_FF(:, 2));
    subtitle(['N=', num2str(data.horizon_lengths(d)), ', FF']);
    xlabel('$t$', "Interpreter","latex")
    % ylabel('$$', "Interpreter","latex")
    sgtitle('Cumulative RotSpeed Error Stage Cost vs time for different horizon lengths');

    figure(10) % Blpitch travel cost
    axes3 = [axes3 subplot(2, 3, d)];
    plot(t_FB, sc_FB(:, 3));
    subtitle(['N=', num2str(data.horizon_lengths(d)), ', FB']);
    axes3 = [axes3 subplot(2, 3, 3 + d)];
    plot(t_FF, sc_FF(:, 3));
    subtitle(['N=', num2str(data.horizon_lengths(d)), ', FF']);
    xlabel('$t$', "Interpreter","latex")
    % ylabel('$$', "Interpreter","latex")
    sgtitle('Blpitch Travel Stage Cost vs time for different horizon lengths');

    data.FB.stage_cost = [data.FB.stage_cost, sum(sc_FB)];
    data.FF.stage_cost = [data.FF.stage_cost, sum(sc_FF)];

end

ylim auto
linkaxes(axes, 'y');
savefig(['./results/', 'StageCost_vs_t']);

%% Plot Total Stage Cost vs N
figure
ax1 = subplot(2, 1, 1);
plot(data.horizon_lengths, data.FB.stage_cost)
subtitle('FB')

ax2 = subplot(2, 1, 2);
plot(data.horizon_lengths, data.FF.stage_cost)
subtitle('FF')

xlabel('$N$', "Interpreter","latex")
ylabel('$\sum_k \text{Stage Cost}$', "Interpreter","latex")
sgtitle('Total Stage Cost vs N')
linkaxes([ax1, ax2], 'y')
savefig(['./results/', 'TotalStageCost_vs_N']);

%% Specify if each case was warm-started or cold-started and whether it converged or was truncated

% ? how do I know if a case was warm-started or cold-started ?
% Juan - using their MPC code w/ different horizon lengths, cold-started
% test-delay, using qpOasis, warm-started, N = 10

% ? how do I know if a case converged or not ?
% "didn't converge" - alg timed out before satisfying convergence criteria
% solve_status, NaN when not running,
% if only 1s, probably meaningless
% qpOases (test_delay) returns complete solve_status
% Juan, Vlaho will know
% MPC runs during gusts, channels returning correct data then
% iterations: small number prob converged, large number may not have
% iterations are usually consistent

% ? what do different blocks correspond to ?
% long-run, split into blocks, concatenated

% ? feedback vs feed-forward ?
% could also compare feedback with feed-forward, MPC takes more iters with
% FF

% ? what does 240... numbers mean ?
% speed of fans in RPM, higher numbers = higher wind speeds

% ? daq ?
% data acquisition

