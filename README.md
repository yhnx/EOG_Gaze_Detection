# Eyris: Gaze Angle from EOG with Baseline-Drift Removal

Estimates the horizontal gaze angle from electrooculography (EOG) while removing the slow baseline drift that corrupts the raw signal. The gain between voltage and gaze angle is calibrated by least squares, and a vector Kalman filter then tracks the gaze angle and the drift together, sample by sample. Implemented in MATLAB.

Course project for BM4112 Medical Electronics and Instrumentation (University of Moratuwa).

## Problem

EOG voltage follows the eye's rotation, so it can be used to find where a person is looking. This is useful for eye-controlled interfaces, communication aids and fatigue monitoring. The signal also contains a baseline that wanders slowly with time. On subject S1 it moves from about 25 µV to about 940 µV over the recording, which looks like roughly 35° of eye movement. Converting voltage to angle without handling this drift gives large errors.

## Method

Each horizontal EOG sample is modelled as

$$y_k = g\,\theta_k + d_k + v_k, \qquad v_k \sim \mathcal{N}(0, R)$$

where $\theta_k$ is the gaze angle (deg), $d_k$ is the baseline drift (µV), $g$ is the gain (µV/deg) and $v_k$ is measurement noise.

**1. Least-squares calibration of the gain.** The first 20 trials have known target angles. For each saccade, the voltage step $\Delta y_i$ between the settled levels before and after is measured. Drift barely changes across one saccade, so it cancels in the difference and $\Delta y_i \approx g\,\Delta\theta_i$. Minimising $\sum_i (\Delta y_i - g\,\Delta\theta_i)^2$ gives

$$\hat g = \frac{\sum_i \Delta\theta_i\,\Delta y_i}{\sum_i \Delta\theta_i^{2}}$$

Differences are used instead of absolute levels, because absolute levels also contain hundreds of µV of drift, which biases $\hat g$. For S1 this gives $\hat g = 27.24$ µV/deg.

**2. Kalman filter.** The state is $\mathbf{x}_k = [\theta_k,\; d_k]^\top$ with measurement matrix $\mathbf{H} = [\hat g \;\; 1]$ and transition matrix $\mathbf{F} = \mathrm{diag}(a, 1)$, where $a = 1 - 1/(\tau f_s)$. Gaze is a weak autoregressive process that relaxes slowly toward the screen centre, and drift is a random walk. The process noise $\mathbf{Q}_k$ switches with eye movement:

- Fixation: $\mathbf{Q} = \mathrm{diag}\big((1-a^2)\sigma_\theta^2,\; q_d\big)$. Gaze is held, drift is free, so slow voltage changes are absorbed as drift.
- Saccade: $\mathbf{Q} = \mathrm{diag}(q_s,\; 0)$. Gaze is free, drift is frozen, so the fast voltage step is credited to gaze.

Saccades are flagged when the smoothed EOG velocity exceeds a threshold of 8 times a robust (MAD) noise level, and each flag is widened by 0.08 s.

Settings tuned on S1 and used unchanged for all subjects:

| Parameter               | Value            |
| ----------------------- | ---------------- |
| Measurement noise R     | 6 µV²            |
| Drift process noise q_d | 60 µV²/sample    |
| Saccade gaze noise q_s  | 1000 deg²/sample |
| Gaze time constant τ    | 5 s              |
| Gaze prior std σ_θ      | 20°              |
| Saccade threshold       | 8 × robust σ     |
| Flag widening           | 0.08 s           |

## Dataset

University of Malta EOG dataset (Barbara, Camilleri and Camilleri), Dataset 1: zero-centred horizontal and vertical bipolar EOG from six healthy subjects. Each subject has 300 trials of 4 s at 256 Hz, in which the subject makes a saccade from the centre to a target, a return saccade, and a blink, with the head fixed on a chin rest. Each subject folder has `EOG.mat`, `ControlSignal.mat` (labels 1 forward saccade, 2 return saccade, 3 blink) and `TargetGA.mat` (target gaze angles in degrees).

The data is not redistributed here. Download it from the dataset page: https://www.um.edu.mt/cbc/ourprojects/EOG/EOGdataset

## Results

Gaze-angle RMSE against the ground truth, using the target angle at the settled part of each saccade window and 0° at the return and blink windows. "Naive" scales the voltage by the gain and subtracts only the first sample, with no drift handling. The gain was calibrated per subject from the first 20 trials. All other settings are the S1 values above.

| Subject  | Naive     | Kalman   | Reduction     |
| -------- | --------- | -------- | ------------- |
| S1       | 20.27°    | 2.20°    | 9×            |
| S2       | 86.30°    | 2.45°    | 35×           |
| S3       | 28.00°    | 3.60°    | 8×            |
| S4       | 156.33°   | 3.96°    | 40×           |
| S5       | 76.75°    | 3.23°    | 24×           |
| S6       | 22.84°    | 2.22°    | 10×           |
| **Mean** | **65.1°** | **2.9°** | **about 22×** |

On S1 the Kalman error is 2.12° at the centre and 2.73° at the targets, and 6.2 % of samples are flagged as saccades.

## Usage

1. Download the dataset and place a subject folder (for example `S1/`) next to the script.
2. Set the data folder in the config section of the script (`cfg.dataDir = 'S1'`).
3. Run the script. It calibrates the gain, runs the Kalman filter, prints the RMSE for the naive and Kalman estimates, and saves the result plots.

If the data folder is not found, the script runs on a synthetic EOG signal with known gaze and drift. Requires MATLAB R2020a or newer (for `exportgraphics`).

## Repository contents

Adjust this list to match the repository.

- Main script with the configuration, data loading, calibration, filtering, metrics and plots
- `kalman_gaze_drift.m`: vector Kalman filter with saccade-switched process noise
- `detect_saccades.m`: velocity-threshold saccade detector
- `ls_calibrate.m`: least-squares gain calibration from saccade steps
- `make_synthetic_eog.m`: synthetic test signal

## Limitations

- The ground truth is the cue angle shown on screen, not a measured eye position, and only the settled parts of each trial are scored.
- The method assumes gaze returns to the screen centre and that drift is slow compared with saccades.
- The settings were tuned on S1 only. A larger drift process noise lowers the error but makes the drift estimate behave more like a fast baseline follower.
- Blinks and the vertical channel are not modelled, and the head is fixed on a chin rest.
- A new calibration is needed for each user and electrode placement.

## Future work

- Add the vertical channel and blink rejection to estimate full two-dimensional gaze direction.
- Evaluate on the non-stationary head dataset from the same group.

## References

1. N. Barbara, T. A. Camilleri and K. P. Camilleri, "A comparison of EOG baseline drift mitigation techniques," _Biomedical Signal Processing and Control_, vol. 57, 101738, 2020. doi:10.1016/j.bspc.2019.101738
2. N. Barbara, T. A. Camilleri and K. P. Camilleri, "Real-time continuous EOG-based gaze angle estimation with baseline drift compensation under stationary head conditions," _Biomedical Signal Processing and Control_, vol. 86, 105282, 2023.
3. Ľ. Hládek, B. Porr and W. O. Brimijoin, "Real-time estimation of horizontal gaze angle by saccade integration using in-ear electrooculography," _PLOS ONE_, 2018. doi:10.1371/journal.pone.0190420

The state-space formulation with gaze and drift estimated together follows the approach of reference 2. The saccade-switched process noise and the specific parameters were developed for this project.
