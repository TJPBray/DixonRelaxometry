% Simulated slice profile for a 3 mm slice
% Hamming-windowed sinc RF pulse ─ MATLAB R2025a
% -------------------------------------------------------------------------
clear; clc;

%% ---- User-adjustable parameters ----------------------------------------
slice_thickness = 3e-3;        % 3 mm (metres)
T               = 4e-3;        % RF-pulse duration (4 ms)
N               = 1024;        % Number of time samples (power of two → FFT speed)
tbw             = 4;           % Time-bandwidth product of the envelope
G               = 10e-3;       % Slice-select gradient (10 mT/m)

%% ---- Physical constant -------------------------------------------------
gamma = 42.57747892e6;         % ^1H gyromagnetic ratio (Hz/T)

%% ---- Build the Hamming-windowed sinc pulse -----------------------------

% Specify sample spacing (s) for the pulse
dt   = T / N;                              

% Get vector of sampling times for the pulse, centred on 0
t    = linspace(-T/2, T/2, N).';           

% Specify the RF bandwidth for the pulse (Hz)
bw_rf = tbw / T;                           

% Specify the pulse waveform (use the inbuilt MATLAB sinc function, which is equal to sin(pi*x)/(pi*x))
rf    = sinc(bw_rf * t);                   

%Show the pulse
figure
subplot(1,3,1)
plot(rf)
xlabel('Time (s)')
ylabel('Amplitude')
xlim([0 1024])
title('Sinc pulse')

% Apply Hamming window
rf    = rf .* hamming(N);         
subplot(1,3,2)
plot(rf)
xlabel('Time (s)')
ylabel('Amplitude')
xlim([0 1024])
title('Sinc pulse with Hamming window')

subplot(1,3,3)
plot(hamming(N))
xlabel('Time (s)')
ylabel('Amplitude')
xlim([0 1024])
title('Hamming window')


%Show the pulse


% ↳ (Optional) scale rf so the time-integral gives the desired flip angle.

%% ---- Compute spatial slice profile -------------------------------------
RFspec   = fftshift(fft(rf));              % Complex spectrum (NB fftshift moves the zero frequency component to the centre of the spectrum)
profile  = abs(RFspec) / max(abs(RFspec)); % Normalised magnitude

figure
subplot(1,4,1)
plot(rf)
xlabel('Time')

subplot(1,3,2)
plot(abs(fft(rf)))
xlabel('Frequency')

subplot(1,3,3)
plot(profile)
xlabel('Frequency with centre at 0')


% Frequency axis and spatial mapping
df   = 1 / T;                              % Frequency resolution (Hz)
f    = (-N/2:N/2-1).' * df;                % Frequency vector (Hz)
z    = f / (gamma * G);                    % Position (metres)

%% ---- Plot --------------------------------------------------------------
figure('Color','w');
plot(z*1e3, profile, 'b', 'LineWidth', 1.5); hold on;
xline(-slice_thickness*1e3/2, 'r--', 'Slice edges');
xline( slice_thickness*1e3/2, 'r--');
xlabel('Position along slice-select axis (mm)');
ylabel('Normalised magnetisation');
title('Simulated slice profile (Hamming-windowed sinc, 3 mm slice)');
grid on; xlim([-10 10]);