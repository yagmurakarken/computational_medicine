clc; clear; close all;

model_dir = 'results/individual_models/';
model_files = dir(fullfile(model_dir, 'model_*.mat'));

numModels = length(model_files);

% Pre-allocate arrays for biomarker values
DMP = zeros(numModels,1);
PeakVm = zeros(numModels,1);
APD90 = zeros(numModels,1);
SpontInterval = zeros(numModels,1);

for i = 1:numModels
    % Load model data explicitly
    data = load(fullfile(model_dir, model_files(i).name));

    Vm = data.Vm;
    t_fixed = data.t_fixed;

    % Biomarker calculations
    DMP(i) = min(Vm(end-5000:end));
    PeakVm(i) = max(Vm(end-5000:end));

    % APD90 calculation explicitly
    Vm_norm = (Vm - DMP(i)) / (PeakVm(i) - DMP(i));
    threshold = 0.1;
    crossings = diff(Vm_norm > threshold);
    AP_starts = find(crossings == 1);
    AP_ends = find(crossings == -1);

    if length(AP_ends) < length(AP_starts)
        AP_starts(end) = [];
    end

    APD90_temp = t_fixed(AP_ends) - t_fixed(AP_starts);
    APD90(i) = mean(APD90_temp);

    % Spontaneous interval calculation
    intervals = diff(t_fixed(AP_starts));
    SpontInterval(i) = mean(intervals);

    % Explicit prints
    fprintf('\nModel %d Biomarkers:\n', i);
    fprintf('DMP (mV): %g\n', DMP(i)*1000);
    fprintf('Peak Vm (mV): %g\n', PeakVm(i)*1000);
    fprintf('APD90 (ms): %g\n', APD90(i)*1000);
    fprintf('Spontaneous Interval (s): %g\n', SpontInterval(i));
end

% Save biomarker results explicitly
save('results/biomarkers.mat', 'DMP', 'PeakVm', 'APD90', 'SpontInterval');
