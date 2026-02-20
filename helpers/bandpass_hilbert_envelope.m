function env = bandpass_hilbert_envelope(eeg, fs, fLow, fHigh, fs_ds)
%BANDPASS_HILBERT_ENVELOPE
% Computes band-limited power (BLP) using bandpass filtering + Hilbert transform.
%
% INPUTS:
%   eeg   : [nChan x nSamples] EEG signal
%   fs    : original sampling frequency (Hz)
%   fLow  : low cutoff frequency (Hz)
%   fHigh : high cutoff frequency (Hz)
%   fs_ds : target sampling frequency for envelope computation (Hz)
%
% OUTPUT:
%   env   : [nChan x nSamples_ds] Hilbert envelope sampled at fs_ds

    if nargin < 5 || isempty(fs_ds)
        fs_ds = fs;  % default: no downsampling
    end

    if fLow <= 0
        error('Invalid fLow: %g Hz', fLow);
    end
    if fs_ds <= 0 || fs_ds > fs
        error('Invalid fs_ds=%g (must be >0 and <= fs=%g)', fs_ds, fs);
    end
    if fHigh >= fs_ds/2
        error('Invalid frequency band: [%g %g] Hz for fs_ds=%g Hz (Nyquist=%g)', ...
            fLow, fHigh, fs_ds, fs_ds/2);
    end

    eeg = double(eeg);

    % -------- downsample (with anti-aliasing) --------
    if fs_ds < fs
        eeg_ds = resample(eeg', fs_ds, fs)';   % [nChan x nSamples_ds]
    else
        eeg_ds = eeg;
    end

    % -------- bandpass filter design (SOS for stability) --------
    filtOrder = 4;
    [sos, g] = butter(filtOrder, [fLow fHigh] / (fs_ds/2), 'bandpass');

    % -------- apply filter (zero-phase) --------
    eeg_filt = filtfilt(sos, g, eeg_ds')';

    % -------- Hilbert transform and envelope --------
    analytic = hilbert(eeg_filt')';
    env = abs(analytic);
end
