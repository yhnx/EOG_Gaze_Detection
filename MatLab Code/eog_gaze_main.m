%% BM4112 - Horizontal EOG gaze-angle estimation with a vector Kalman filter
%  Problem : EOG baseline drift corrupts gaze-angle estimates.
%  Method  : (1) Least-squares calibration  EOG = g*angle + c
%            (2) 2-state Kalman filter [gaze angle; drift] with saccade-adaptive process noise
%            (3) ML estimate of measurement-noise variance R from a quiet segment
%  Run     : edit the CONFIG block, then run this script. With no data file it runs on synthetic EOG.
clear; clc; close all;

%% ---------------- CONFIG (edit for your dataset) ----------------
cfg.file    = '';          % e.g. 'eog_data.csv'  ('' -> synthetic demo)
cfg.col     = 1;           % column holding the HORIZONTAL EOG channel
cfg.fs      = 250;         % sampling rate (Hz) - check the dataset's description file
cfg.units   = 1;           % multiply raw data by this to get uV (e.g. 1000 if file is in mV)
cfg.truth   = [];          % optional ground-truth gaze angle (deg), same length as signal
cfg.calSeg  = [];          % calibration fixations IN TIME ORDER: rows [startSample endSample angleDeg]
cfg.g       = 16;          % fallback gain (uV/deg) if no calibration segments are given
cfg.quietSeg= [];          % [startSample endSample] fixation with no saccade, to estimate R (ML)
% Kalman tuning
cfg.qDrift  = 0.5;         % drift process noise (uV^2/sample)  - larger = faster drift tracking
cfg.qSacc   = 400;         % gaze process noise during saccade (deg^2/sample)
cfg.gate    = 4;           % innovation gate (sigmas) that declares a saccade
cfg.holdSec = 0.06;        % saccade-mode hold time (s)

%% ---------------- LOAD DATA ----------------
if isempty(cfg.file) || ~isfile(cfg.file)
    fprintf('No data file found -> using synthetic EOG (known ground truth).\n');
    [y, cfg.truth, driftTrue, cfg.fs] = make_synthetic_eog(1);
    % calibration = first 10 fixations (20 s): known angles, in time order
    cfg.calSeg = zeros(10,3);
    for m = 0:9
        if m == 0, s0 = round(0.5*cfg.fs); else, s0 = round((2*m+0.3)*cfg.fs); end
        e0 = round((2*m+1.9)*cfg.fs);
        cfg.calSeg(m+1,:) = [s0 e0 cfg.truth(round((s0+e0)/2))];
    end
    cfg.quietSeg = [round(0.2*cfg.fs) round(1.9*cfg.fs)];
else
    raw = readmatrix(cfg.file);
    y   = raw(:, cfg.col) * cfg.units;
    driftTrue = [];
end
y = y(:);  fs = cfg.fs;  N = numel(y);  t = (0:N-1)'/fs;

%% ---------------- 1. LEAST-SQUARES CALIBRATION ----------------
if ~isempty(cfg.calSeg)
    lv = zeros(size(cfg.calSeg,1),1);
    for i = 1:size(cfg.calSeg,1)
        lv(i) = mean(y(cfg.calSeg(i,1):cfg.calSeg(i,2)));
    end
    g = ls_calibrate(lv, cfg.calSeg(:,3));
    fprintf('LS calibration: g = %.2f uV/deg\n', g);
    if abs(g) < 1, warning('Calibration gain implausibly small - add more/larger-angle targets.'); end
else
    g = cfg.g;  fprintf('No calibration segments: using g = %.1f uV/deg\n', g);
end

%% ---------------- 2. ML ESTIMATE OF MEASUREMENT NOISE ----------------
if ~isempty(cfg.quietSeg)
    seg = y(cfg.quietSeg(1):cfg.quietSeg(2));
    R = var(diff(seg))/2;          % ML variance of white noise, insensitive to slow drift
else
    R = (0.05*std(y))^2;  warning('No quiet segment given: R set heuristically.');
end
fprintf('Measurement noise: R = %.1f uV^2 (sigma = %.1f uV)\n', R, sqrt(R));

%% ---------------- 3. KALMAN FILTER ----------------
[theta, drift, isSacc, ~] = kalman_gaze_drift(y, g, R, cfg.qDrift, cfg.qSacc, fs, cfg.gate, cfg.holdSec);

% Naive estimate (no drift handling) for comparison: angle = (y - first-sample offset)/g
thetaNaive = (y - y(1)) / g;

%% ---------------- 4. METRICS ----------------
if ~isempty(cfg.truth)
    rmse = @(a,b) sqrt(mean((a-b).^2));
    fixMask = ~(conv(double(isSacc), ones(25,1), 'same') > 0);   % exclude saccade transients
    fprintf('\nRMSE gaze angle (all samples)    : naive %.2f deg  |  Kalman %.2f deg\n', ...
            rmse(thetaNaive,cfg.truth), rmse(theta,cfg.truth));
    fprintf('RMSE gaze angle (fixations only) : naive %.2f deg  |  Kalman %.2f deg\n', ...
            rmse(thetaNaive(fixMask),cfg.truth(fixMask)), rmse(theta(fixMask),cfg.truth(fixMask)));
    if ~isempty(driftTrue)
        fprintf('Drift tracking RMSE              : %.1f uV\n', rmse(drift - drift(1) + driftTrue(1), driftTrue));
    end
else
    fprintf('\nNo ground truth supplied: judge by plots (flat fixations, return to 0 deg at centre).\n');
end
fprintf('Saccade samples flagged: %.1f %%\n', 100*mean(isSacc));

%% ---------------- 5. PLOTS (16:9, large fonts for the video) ----------------
set(0,'DefaultAxesFontSize',16,'DefaultLineLineWidth',1.5);
f = figure('Position',[50 50 1920 1080],'Color','w');
subplot(3,1,1);
plot(t,y,'Color',[.6 .6 .6]); hold on;
plot(t, g*theta + drift, 'b'); plot(t, drift, 'r', 'LineWidth', 2.5);
if ~isempty(driftTrue), plot(t, driftTrue - driftTrue(1) + drift(1), 'k--'); end
ylabel('EOG (\muV)'); title('Raw EOG, estimated drift and filtered signal'); grid on;
legend('Raw EOG','Kalman fit','Estimated drift','True drift (synthetic)','Location','best');
subplot(3,1,2);
plot(t,thetaNaive,'Color',[.7 .7 .7]); hold on; plot(t,theta,'b','LineWidth',2);
if ~isempty(cfg.truth), plot(t,cfg.truth,'k--'); end
ylabel('Gaze angle (deg)'); title('Gaze angle: naive vs Kalman'); grid on;
legend('Naive (no drift removal)','Kalman','Ground truth','Location','best');
subplot(3,1,3);
plot(t, isSacc, 'm'); ylim([-0.1 1.1]); ylabel('Saccade flag'); xlabel('Time (s)');
title('Saccade detection from innovation gating'); grid on;
exportgraphics(f,'result_gaze_kalman.png','Resolution',100);

if ~isempty(cfg.truth)
    f2 = figure('Position',[50 50 1920 1080],'Color','w');
    subplot(1,2,1); hold on;
    plot(t,thetaNaive-cfg.truth,'Color',[.7 .7 .7]); plot(t,theta-cfg.truth,'b','LineWidth',2);
    xlabel('Time (s)'); ylabel('Error (deg)'); title('Gaze-angle error'); legend('Naive','Kalman'); grid on;
    subplot(1,2,2); histogram(theta(fixMask)-cfg.truth(fixMask),40,'FaceColor','b'); hold on;
    histogram(thetaNaive(fixMask)-cfg.truth(fixMask),40,'FaceColor',[.7 .7 .7],'FaceAlpha',.6);
    xlabel('Fixation error (deg)'); title('Error distribution (fixations)'); legend('Kalman','Naive'); grid on;
    exportgraphics(f2,'result_error.png','Resolution',100);
end
