function [theta, drift, isSaccade, S_inn] = kalman_gaze_drift(y, g, R, qDrift, qSacc, fs, gate, holdSec)
%KALMAN_GAZE_DRIFT  2-state (vector) Kalman filter: gaze angle + baseline drift from one EOG channel.
%
%   State vector   x = [theta ; d]      theta = gaze angle (deg), d = baseline drift (uV)
%   Measurement    y_k = g*theta_k + d_k + v_k,   v ~ N(0,R)        ->  H = [g 1]
%   Transition     x_k = x_{k-1} + w_k                              ->  F = I (2x2)
%
%   Two process-noise modes (switched by innovation gating):
%     FIXATION : gaze is held, drift wanders   Q = diag(0 , qDrift)
%     SACCADE  : gaze jumps, drift is frozen   Q = diag(qSacc , 0)
%   A saccade is declared when |innovation| > gate*sqrt(S) under the fixation model,
%   and the saccade mode is held for holdSec seconds so the filter can settle.
%
%   Inputs : y (uV, column), g (uV/deg), R (uV^2), qDrift (uV^2/sample), qSacc (deg^2/sample),
%            fs (Hz), gate (e.g. 4), holdSec (e.g. 0.06)
%   Outputs: theta (deg), drift (uV), isSaccade (logical), S_inn (innovation variance)

y = y(:);  N = numel(y);
H = [g 1];
x = [0; y(1)];                 % start at centre gaze, drift = first sample
P = diag([100, 100]);          % loose initial uncertainty
Qf = diag([0, qDrift]);
Qs = diag([qSacc, 0]);
holdN = round(holdSec*fs);  cnt = 0;

theta = zeros(N,1); drift = zeros(N,1); isSaccade = false(N,1); S_inn = zeros(N,1);

for k = 1:N
    % ---- predict + innovation under the fixation model ----
    Pp = P + Qf;
    S  = H*Pp*H' + R;
    nu = y(k) - H*x;                       % innovation
    % ---- saccade gating ----
    if cnt == 0 && abs(nu) > gate*sqrt(S)
        cnt = holdN;
    end
    if cnt > 0
        Pp = P + Qs;  S = H*Pp*H' + R;     % re-predict with saccade model
        cnt = cnt - 1;  isSaccade(k) = true;
    end
    % ---- update ----
    K = (Pp*H')/S;                         % 2x1 Kalman gain
    x = x + K*nu;
    P = (eye(2) - K*H)*Pp;
    theta(k) = x(1);  drift(k) = x(2);  S_inn(k) = S;
end
end
