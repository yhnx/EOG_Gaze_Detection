function isSacc = detect_saccades(y, fs, thr, holdSec)
%DETECT_SACCADES  Velocity-threshold saccade detector (robust to slow drift).
%   Smooth over 6 samples, take the velocity, and flag samples whose speed exceeds
%   thr x (robust MAD sigma of the velocity). Flags are widened by holdSec seconds.
w  = 6;
ys = conv(y(:), ones(w,1)/w, 'same');
v  = zeros(size(ys));  v(w+1:end) = (ys(w+1:end) - ys(1:end-w)) * fs / w;
sig = median(abs(v - median(v))) / 0.6745;
isSacc = abs(v) > thr*sig;
isSacc = conv(double(isSacc), ones(max(1,round(holdSec*fs)),1), 'same') > 0;
end
