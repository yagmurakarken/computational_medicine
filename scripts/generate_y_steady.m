clc; clear; close all;

addpath('cellular_models/');

% Simulation settings
options = odeset('MaxStep',1e-3,'InitialStep',2e-5);
stimFlag = 0;
tDrugApplication = 10000;
phenotype = 'nodal';
t_long = 0:0.001:2000;  % simulate for 2000s

% Start with default Y0 (ventricular-like default state)
Y0_default = [-0.070, 0.32, 0.0002, 0, 0, 1, 1, 1, 0, 1, 0, ...
              0.75, 0.75, 0, 0.1, 1, 0, 9.2, 0, 0.75, 0.3, 0.9, 0.1];

empty_params = struct();  % Placeholder for optional param argument

[~, Yc] = ode15s(@Paci2020, t_long, Y0_default, options, stimFlag, ...
                 tDrugApplication, phenotype, ...
                 1, 1, 1, 1, 1, 1, 1, empty_params);  % Simulate nodal phenotype

% Save the final state as new initial condition
Y0_nodal = Yc(end, :);
save('Y0_nodal.mat', 'Y0_nodal');
fprintf('Saved nodal-specific steady-state Y0 to Y0_nodal.mat\n');

%% Optional: visualize Vm and Cai to confirm stability
Vm = Yc(:,1);
Cai = Yc(:,3);

% Full trace
figure;
subplot(2,1,1);
plot(t_long, Vm*1e3); ylabel('Vm (mV)');
title('Vm Over Time - Nodal Phenotype (2000s)'); grid on;

subplot(2,1,2);
plot(t_long, Cai*1e6); xlabel('Time (s)'); ylabel('Cai (\muM)');
title('Cai Over Time - Nodal Phenotype (2000s)'); grid on;

% Zoom into last 20s
figure;
subplot(2,1,1);
plot(t_long, Vm*1e3); xlim([1980 2000]);
title('Vm (mV) - Last 20s'); ylabel('Vm'); grid on;

subplot(2,1,2);
plot(t_long, Cai*1e6); xlim([1980 2000]);
title('Cai (\muM) - Last 20s'); xlabel('Time (s)'); ylabel('Cai'); grid on;
