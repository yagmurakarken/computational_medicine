%% run_atrial.m
clc; clear; close all;
addpath('cellular_models/');

%% Simulation options
options = odeset('MaxStep',1e-3,'InitialStep',2e-5);

%% Initial conditions (from Paci2020)
Y0 = [-0.070, 0.32, 0.0002, 0, 0, 1, 1, 1, 0, 1, 0, ...
      0.75, 0.75, 0, 0.1, 1, 0, 9.2, 0, 0.75, 0.3, 0.9, 0.1];

%% Simulation parameters
phenotype = 'atrial';
stimFlag = 0;
tDrugApplication = 10000;
INaFRedMed = 1; 
INaLRedMed = 1; 
ICaLRedMed = 1;
IKrRedMed = 1; 
IKsRedMed = 1; 
INaCaRedMed = 1; 
IfRedMed = 1;

%% Simulation time
t_fixed = 0:0.001:800;

%% Run simulation
[t, Yc] = ode15s(@Paci2020, t_fixed, Y0, options, stimFlag, tDrugApplication, phenotype, ...
                 INaFRedMed, INaLRedMed, ICaLRedMed, IKrRedMed, IKsRedMed, INaCaRedMed, IfRedMed);

%% Ionic and Calcium Handling Currents Computation
INa = zeros(size(t)); 
ICaL = zeros(size(t)); 
INaCa = zeros(size(t));
If = zeros(size(t)); 
IKr = zeros(size(t)); 
IKs = zeros(size(t)); 
IK1 = zeros(size(t));
Iup = zeros(size(t)); 
Irel = zeros(size(t));

for i = 1:length(t)
    [~, dati] = Paci2020(t(i), Yc(i,:), stimFlag, tDrugApplication, phenotype, ...
                         INaFRedMed, INaLRedMed, ICaLRedMed, IKrRedMed, IKsRedMed, INaCaRedMed, IfRedMed);
    INa(i)   = dati(1);
    If(i)    = dati(2);
    ICaL(i)  = dati(3);
    IKr(i)   = dati(6);
    IKs(i)   = dati(5);
    IK1(i)   = dati(7);
    INaCa(i) = dati(8);
    Irel(i)  = dati(13);
    Iup(i)   = dati(14);
end

%% Define Vm
Vm = Yc(:,1);

%% Save results
save('results/atrial_results.mat', 't', 'Yc', 'Vm', 'INa', 'ICaL', 'INaCa', 'If', 'IKr', 'IKs', 'IK1', 'Iup', 'Irel');

%% Verification Plot
figure;
plot(t, Vm*1e3, 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('Vm (mV)');
title('Atrial-like hPSC-CM Action Potential (Paci2020)');
grid on; xlim([790 800]);

%% Biomarker Analysis (last 20 seconds)
analysis_window = t >= 780 & t <= 800;
t_analysis = t(analysis_window);
Vm_analysis = Vm(analysis_window)*1e3;

peak_indices = find(diff(sign(diff(Vm_analysis))) < 0) + 1;
peak_threshold = 0;
valid_peaks = Vm_analysis(peak_indices) > peak_threshold;
pks = Vm_analysis(peak_indices(valid_peaks));
locs = t_analysis(peak_indices(valid_peaks));

spont_intervals = diff(locs);
mean_spont_interval = mean(spont_intervals);

DMP = zeros(length(locs)-1, 1);
for i = 1:length(locs)-1
    segment = (t_analysis >= locs(i)) & (t_analysis <= locs(i+1));
    DMP(i) = min(Vm_analysis(segment));
end
mean_DMP = mean(DMP);
mean_peak_voltage = mean(pks);

APD90 = zeros(length(locs)-1, 1);
for i = 1:length(locs)-1
    AP_segment_time = t_analysis((t_analysis >= locs(i)) & (t_analysis <= locs(i)+mean_spont_interval));
    AP_segment_Vm = Vm_analysis((t_analysis >= locs(i)) & (t_analysis <= locs(i)+mean_spont_interval));

    peak_voltage = max(AP_segment_Vm);
    dmp_voltage = min(AP_segment_Vm);
    repol_90_voltage = peak_voltage - 0.9*(peak_voltage - dmp_voltage);
    post_peak_indices = find(AP_segment_time > locs(i));
    repol_cross_idx = find(AP_segment_Vm(post_peak_indices) <= repol_90_voltage, 1, 'first');
    if ~isempty(repol_cross_idx)
        APD90(i) = (AP_segment_time(post_peak_indices(repol_cross_idx)) - locs(i))*1000;
    else
        APD90(i) = NaN;
    end
end
mean_APD90 = nanmean(APD90);

%% Display Biomarkers
fprintf('Computed Biomarkers:\n');
fprintf('Diastolic Membrane Potential (DMP): %.2f mV\n', mean_DMP);
fprintf('Peak Voltage: %.2f mV\n', mean_peak_voltage);
fprintf('APD90: %.2f ms\n', mean_APD90);
fprintf('Spontaneous Beating Interval: %.2f s\n', mean_spont_interval);