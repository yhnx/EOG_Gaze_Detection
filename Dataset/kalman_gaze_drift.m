function [theta, drift] = kalman_gaze_drift(y, g, R, qDrift, qSacc, isSacc, tau, sigTheta, fs)
%KALMAN_GAZE_DRIFT  Vector Kalman filter estimating gaze angle AND baseline drift from one EOG channel.
%
%   State        x_k = [theta_k ; d_k]    theta = gaze angle (deg),  d = baseline drift (uV)
%   Measurement  y_k = g*theta_k + d_k + v_k,   v ~ N(0,R)            H = [g  1]
%   Transition   x_k = F*x_{k-1} + w_k,         F = diag(a, 1)
%                a = 1 - 1/(tau*fs)  : gaze is a weak AR(1) process pulled toward the screen centre
%                                      (prior N(0,sigTheta^2)); drift is a random walk
%   Process noise (switched per sample by the saccade flag isSacc):
%       fixation : Q = diag((1-a^2)*sigTheta^2 , qDrift)   gaze ~held, drift wanders
%       saccade  : Q = diag(qSacc , 0)                     gaze jumps, drift frozen
%
%   Inputs : y (uV), g (uV/deg), R (uV^2), qDrift (uV^2/sample), qSacc (deg^2/sample),
%            isSacc (logical, same length as y), tau (s), sigTheta (deg), fs (Hz)
%   Outputs: theta (deg), drift (uV)

y = y(:);  N = numel(y);
a  = 1 - 1/(tau*fs);
F  = diag([a 1]);
H  = [g 1];
QF = diag([(1-a^2)*sigTheta^2, qDrift]);
QS = diag([qSacc, 0]);
x  = [0; y(1)];              % start at centre gaze, drift = first sample
P  = diag([100 100]);
theta = zeros(N,1);  drift = zeros(N,1);
for k = 1:N
    if isSacc(k), Q = QS; else, Q = QF; end
    x = F*x;                  P = F*P*F' + Q;          % predict
    S = H*P*H' + R;           K = (P*H')/S;            % innovation variance, Kalman gain (2x1)
    x = x + K*(y(k) - H*x);   P = (eye(2) - K*H)*P;    % update
    theta(k) = x(1);  drift(k) = x(2);
end
end
