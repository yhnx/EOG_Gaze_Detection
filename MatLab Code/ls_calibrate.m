function [g, sig2] = ls_calibrate(yLevels, angles)
%LS_CALIBRATE  Least-squares gain from calibration fixations (drift-robust, step-based).
%   yLevels : mean EOG (uV) of each calibration fixation, IN TIME ORDER
%   angles  : known gaze angle (deg) of each fixation
%   Absolute levels are corrupted by baseline drift, so we fit the CHANGES between
%   consecutive fixations (the saccade step heights), where slow drift cancels:
%        dy_i = g * dtheta_i + noise        ->   g = (dtheta' * dy) / (dtheta' * dtheta)
%   Under Gaussian noise this LS solution is also the maximum-likelihood estimate.
yLevels = yLevels(:);  angles = angles(:);
dy  = diff(yLevels);   dth = diff(angles);
g   = (dth' * dy) / (dth' * dth);          % LS slope through the origin
res = dy - g*dth;
sig2 = sum(res.^2) / max(numel(dy)-1, 1);  % residual variance (uV^2)
end
