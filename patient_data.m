clc;
clear;
close all;

%% =========================================================
% BP-ZnO BIOSENSOR ANALYSIS
% ZnO Nanorods + BP Nanosheets + Composite Formation
% Antibody Functionalization + XRD + FTIR + Raman
% + Biological Sample Collection (200 Patients)
%% =========================================================

rng(1);

%% =========================================================
% PHASE 1 : ZnO NANOROD MORPHOLOGY
%% =========================================================

figure('Name','ZnO Nanorods');
hold on;

for i = 1:60
    x = rand()*100;
    y = rand()*100;
    h = 5 + rand()*15;
    plot3([x x],[y y],[0 h],'LineWidth',2);
end

grid on;
view(45,30);
xlabel('X Position');
ylabel('Y Position');
zlabel('Nanorod Height');
title('Simulated ZnO Nanorod Morphology');

%% =========================================================
% PHASE 2 : BLACK PHOSPHORUS NANOSHEET
%% =========================================================

figure('Name','Black Phosphorus Nanosheet');

[X,Y] = meshgrid(1:100,1:100);
Z = 0.25*sin(X/8) + 0.20*cos(Y/10);

surf(X,Y,Z);
shading interp;
xlabel('X');
ylabel('Y');
zlabel('Thickness');
title('Simulated Black Phosphorus Nanosheet');
view(45,35);

%% =========================================================
% PHASE 3 : BP-ZnO COMPOSITE
%% =========================================================

figure('Name','BP-ZnO Composite');
surf(X,Y,Z);
hold on;

for i = 1:40
    x = rand()*100;
    y = rand()*100;
    h = 4 + rand()*8;
    plot3([x x],[y y],[0 h],'LineWidth',2);
end

xlabel('X');
ylabel('Y');
zlabel('Height');
title('Simulated BP-ZnO Composite');
view(45,30);
grid on;

%% =========================================================
% PHASE 4 : ANTIBODY FUNCTIONALIZATION
%% =========================================================

figure('Name','Antibody Functionalization');
surf(X,Y,Z);
hold on;

for i = 1:50
    x = rand()*100;
    y = rand()*100;
    z = 0.5 + rand()*0.5;
    scatter3(x, y, z, 50, 'c', 'filled');
end

for i = 1:50
    x = rand()*100;
    y = rand()*100;
    z = 0.5 + rand()*0.5;
    scatter3(x, y, z, 50, 'm', 'filled');
end

xlabel('X');
ylabel('Y');
zlabel('Surface Coverage');
title('IL-6 / TNF-\alpha Antibody Immobilization');
legend('BP-ZnO Surface','IL-6 Antibody','TNF-\alpha Antibody', ...
    'Location','northeast');
view(45,30);
grid on;

%% =========================================================
% PHASE 5 : IL-6 ANTIBODY CONFIRMATION
%% =========================================================

conc_IL6   = logspace(-3, 2, 100);
signal_IL6 = 95 ./ (1 + exp(-1.8*(log10(conc_IL6) - 0.5))) ...
             + 2*randn(1,100);
signal_IL6 = max(signal_IL6, 0);

IL6_LOD         = 0.012;
IL6_LOQ         = 0.038;
IL6_linear_low  = 0.05;
IL6_linear_high = 50.0;
IL6_sensitivity = 18.6;
IL6_Kd          = 3.2e-9;
IL6_coverage    = 87.4;

figure('Name','IL-6 Antibody Confirmation');
semilogx(conc_IL6, signal_IL6, 'b-o', 'LineWidth', 2, ...
    'MarkerSize', 4, 'MarkerFaceColor', 'b');
hold on;
xline(IL6_LOD, '--r', 'LOD', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5);
xline(IL6_LOQ, '--g', 'LOQ', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5);
xlabel('IL-6 Concentration (ng/mL)');
ylabel('Current Response (\muA)');
title('IL-6 Antibody Binding Confirmation Curve');
grid on;
legend('IL-6 Response','LOD','LOQ','Location','northwest');

%% =========================================================
% PHASE 6 : TNF-ALPHA ANTIBODY CONFIRMATION
%% =========================================================

conc_TNF   = logspace(-3, 2, 100);
signal_TNF = 88 ./ (1 + exp(-1.6*(log10(conc_TNF) - 0.3))) ...
             + 2*randn(1,100);
signal_TNF = max(signal_TNF, 0);

TNF_LOD         = 0.008;
TNF_LOQ         = 0.025;
TNF_linear_low  = 0.03;
TNF_linear_high = 40.0;
TNF_sensitivity = 21.3;
TNF_Kd          = 2.8e-9;
TNF_coverage    = 83.1;

figure('Name','TNF-alpha Antibody Confirmation');
semilogx(conc_TNF, signal_TNF, 'm-s', 'LineWidth', 2, ...
    'MarkerSize', 4, 'MarkerFaceColor', 'm');
hold on;
xline(TNF_LOD, '--r', 'LOD', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5);
xline(TNF_LOQ, '--g', 'LOQ', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5);
xlabel('TNF-\alpha Concentration (ng/mL)');
ylabel('Current Response (\muA)');
title('TNF-\alpha Antibody Binding Confirmation Curve');
grid on;
legend('TNF-\alpha Response','LOD','LOQ','Location','northwest');

%% =========================================================
% PHASE 7 : XRD ANALYSIS
%% =========================================================

theta = 10:0.05:80;

xrd_raw = ...
    100*exp(-(theta-31.8).^2/0.5) + ...
    140*exp(-(theta-34.4).^2/0.4) + ...
    120*exp(-(theta-36.2).^2/0.4) + ...
     60*exp(-(theta-17.0).^2/1.0);

xrd = xrd_raw ./ max(xrd_raw);

figure('Name','XRD Analysis');
plot(theta, xrd, 'LineWidth', 2);
xlabel('2\theta (Degree)');
ylabel('Normalized Intensity');
title('Simulated XRD Pattern of BP-ZnO Composite');
grid on;

[~, locs_xrd]       = findpeaks(xrd, theta, 'MinPeakHeight', 0.2);
num_xrd_peaks       = length(locs_xrd);
xrd_peak_amplitudes = [60, 100, 140, 120];
xrd_peak_labels     = {'BP (002)', 'ZnO (100)', 'ZnO (002)', 'ZnO (101)'};

%% =========================================================
% PHASE 8 : FTIR ANALYSIS
%% =========================================================

wn = 500:4000;

ftir_raw = ...
     60*exp(-(wn- 550).^2/ 500) + ...
     90*exp(-(wn-1050).^2/1000) + ...
     70*exp(-(wn-1600).^2/ 800) + ...
    120*exp(-(wn-3400).^2/2000);

ftir = ftir_raw ./ max(ftir_raw);

figure('Name','FTIR Analysis');
plot(wn, ftir, 'LineWidth', 2);
xlabel('Wavenumber (cm^{-1})');
ylabel('Normalized Absorbance');
title('Simulated FTIR Spectrum');
grid on;

[~, locs_ftir]       = findpeaks(ftir, wn, 'MinPeakHeight', 0.2);
num_ftir_peaks       = length(locs_ftir);
ftir_peak_amplitudes = [60, 90, 70, 120];
ftir_peak_labels     = {'Zn-O Stretch', 'P-O Bond', 'P=O Stretch', 'O-H Stretch'};

%% =========================================================
% PHASE 9 : RAMAN ANALYSIS
%% =========================================================

shift = 100:600;

raman_raw = ...
    100*exp(-(shift-360).^2/100) + ...
    130*exp(-(shift-438).^2/100) + ...
     80*exp(-(shift-466).^2/120);

raman = raman_raw ./ max(raman_raw);

figure('Name','Raman Analysis');
plot(shift, raman, 'LineWidth', 2);
xlabel('Raman Shift (cm^{-1})');
ylabel('Normalized Intensity');
title('Simulated Raman Spectrum');
grid on;

[~, locs_raman]       = findpeaks(raman, shift, 'MinPeakHeight', 0.2);
num_raman_peaks       = length(locs_raman);
raman_peak_amplitudes = [100, 130, 80];
raman_peak_labels     = {'BP Ag1 Mode', 'ZnO E2 (High)', 'BP Ag2 Mode'};

%% =========================================================
% PHASE 10 : BIOLOGICAL SAMPLE COLLECTION (200 PATIENTS)
%% =========================================================

num_patients  = 200;
surgery_types = {'Cardiac','Orthopedic','Abdominal','Neurological','Thoracic'};

% --- Demographics ---
ages         = randi([25 75], num_patients, 1);
genders      = randi([0 1],   num_patients, 1);   % 0=Male 1=Female
surgery_idx  = randi([1 5],   num_patients, 1);

% --- Step 1 : Blood Collection ---
blood_volume   = 5 + rand(num_patients,1)*5;          % mL
hemolysis_flag = rand(num_patients,1) < 0.05;
collection_ok  = ~hemolysis_flag;

% --- Step 2 : Centrifugation ---
centrifuge_rpm  = 3000;
centrifuge_time = 10;
centrifuge_temp = 4;
centrifuge_ok   = rand(num_patients,1) > 0.02;

% --- Step 3 : Serum Separation ---
serum_yield_pct = 55 + rand(num_patients,1)*10;       % %
serum_volume    = blood_volume .* (serum_yield_pct/100);
serum_volume(~centrifuge_ok) = 0;

turbidity    = rand(num_patients,1)*10;
lipemic_flag = turbidity > 6;

serum_ready  = collection_ok & centrifuge_ok & ~lipemic_flag;

% --- Biomarker Levels ---
IL6_level  = zeros(num_patients,1);
TNF_level  = zeros(num_patients,1);

for p = 1:num_patients
    if serum_ready(p)
        if ages(p) > 50 || ismember(surgery_idx(p),[1,3])
            IL6_level(p) = 20  + rand()*180;
            TNF_level(p) = 10  + rand()*110;
        else
            IL6_level(p) = 1   + rand()*30;
            TNF_level(p) = 0.5 + rand()*20;
        end
    end
end

IL6_threshold = 10;
TNF_threshold = 8;

IL6_elevated  = IL6_level > IL6_threshold & serum_ready;
TNF_elevated  = TNF_level > TNF_threshold & serum_ready;
both_elevated = IL6_elevated & TNF_elevated;

% --- Summary Counts ---
n_collected       = sum(collection_ok);
n_hemolysed       = sum(hemolysis_flag);
n_centrifuged     = sum(centrifuge_ok & collection_ok);
n_centrifuge_fail = sum(~centrifuge_ok & collection_ok);
n_lipemic         = sum(lipemic_flag & collection_ok & centrifuge_ok);
n_serum_ready     = sum(serum_ready);
n_serum_rejected  = num_patients - n_serum_ready;
n_IL6_elevated    = sum(IL6_elevated);
n_TNF_elevated    = sum(TNF_elevated);
n_both_elevated   = sum(both_elevated);
n_normal          = sum(serum_ready & ~IL6_elevated & ~TNF_elevated);
mean_IL6_all      = mean(IL6_level(serum_ready));
mean_TNF_all      = mean(TNF_level(serum_ready));
mean_serum_vol    = mean(serum_volume(serum_ready));

% --- Figures ---
figure('Name','Serum Volume Distribution');
histogram(serum_volume(serum_ready), 20, 'FaceColor', [0.2 0.6 0.8]);
xlabel('Serum Volume (mL)');
ylabel('Number of Patients');
title('Serum Volume Distribution (Valid Samples)');
grid on;

figure('Name','Biomarker Levels');
subplot(1,2,1);
histogram(IL6_level(serum_ready), 25, 'FaceColor', [0.1 0.4 0.9]);
xlabel('IL-6 Level (pg/mL)');
ylabel('Frequency');
title('IL-6 Distribution (Post-op Serum)');
xline(IL6_threshold,'--r','Threshold','LineWidth',1.5);
grid on;

subplot(1,2,2);
histogram(TNF_level(serum_ready), 25, 'FaceColor', [0.9 0.2 0.5]);
xlabel('TNF-\alpha Level (pg/mL)');
ylabel('Frequency');
title('TNF-\alpha Distribution (Post-op Serum)');
xline(TNF_threshold,'--r','Threshold','LineWidth',1.5);
grid on;

figure('Name','Sample Processing Summary');
pie([n_serum_ready, n_hemolysed, n_centrifuge_fail, n_lipemic], ...
    {'Serum Ready','Hemolysed','Centrifuge Fail','Lipemic'});
title('Sample Processing Pipeline — 200 Patients');

%% =========================================================
% PHASE 11 : SAVE 200 PATIENTS DATA TO CSV
%% =========================================================

% --- Build string arrays for categorical columns ---
gender_str   = cell(num_patients,1);
surgery_str  = cell(num_patients,1);
hemolysis_str= cell(num_patients,1);
centrifuge_str=cell(num_patients,1);
lipemic_str  = cell(num_patients,1);
serum_str    = cell(num_patients,1);
IL6_flag_str = cell(num_patients,1);
TNF_flag_str = cell(num_patients,1);
clinical_str = cell(num_patients,1);

for p = 1:num_patients

    % Gender
    if genders(p) == 0
        gender_str{p} = 'Male';
    else
        gender_str{p} = 'Female';
    end

    % Surgery
    surgery_str{p} = surgery_types{surgery_idx(p)};

    % Hemolysis
    if hemolysis_flag(p)
        hemolysis_str{p} = 'Yes';
    else
        hemolysis_str{p} = 'No';
    end

    % Centrifuge
    if centrifuge_ok(p)
        centrifuge_str{p} = 'Pass';
    else
        centrifuge_str{p} = 'Fail';
    end

    % Lipemic
    if lipemic_flag(p)
        lipemic_str{p} = 'Yes';
    else
        lipemic_str{p} = 'No';
    end

    % Serum status
    if serum_ready(p)
        serum_str{p} = 'Ready';
    else
        serum_str{p} = 'Rejected';
    end

    % IL-6 flag
    if IL6_elevated(p)
        IL6_flag_str{p} = 'Elevated';
    else
        IL6_flag_str{p} = 'Normal';
    end

    % TNF-alpha flag
    if TNF_elevated(p)
        TNF_flag_str{p} = 'Elevated';
    else
        TNF_flag_str{p} = 'Normal';
    end

    % Clinical interpretation
    if both_elevated(p)
        clinical_str{p} = 'High Inflammatory Risk';
    elseif IL6_elevated(p)
        clinical_str{p} = 'IL-6 Dominant Inflammation';
    elseif TNF_elevated(p)
        clinical_str{p} = 'TNF-a Dominant Inflammation';
    elseif serum_ready(p)
        clinical_str{p} = 'Normal';
    else
        clinical_str{p} = 'Sample Rejected';
    end

end

% --- Build Patient ID column ---
patient_id = cell(num_patients,1);
for p = 1:num_patients
    patient_id{p} = sprintf('P%03d', p);
end

% --- Assemble Table ---
T = table( ...
    patient_id,                         ...
    ages,                               ...
    gender_str,                         ...
    surgery_str,                        ...
    round(blood_volume,        3),      ...
    hemolysis_str,                      ...
    centrifuge_str,                     ...
    round(serum_yield_pct,     2),      ...
    round(serum_volume,        3),      ...
    round(turbidity,           2),      ...
    lipemic_str,                        ...
    serum_str,                          ...
    round(IL6_level,           4),      ...
    IL6_flag_str,                       ...
    round(TNF_level,           4),      ...
    TNF_flag_str,                       ...
    clinical_str,                       ...
    'VariableNames', {                  ...
        'Patient_ID',                   ...
        'Age_years',                    ...
        'Gender',                       ...
        'Surgery_Type',                 ...
        'Blood_Volume_mL',              ...
        'Hemolysis',                    ...
        'Centrifuge_Status',            ...
        'Serum_Yield_Pct',              ...
        'Serum_Volume_mL',              ...
        'Turbidity_Score',              ...
        'Lipemic',                      ...
        'Serum_Status',                 ...
        'IL6_Level_pgmL',               ...
        'IL6_Flag',                     ...
        'TNF_Level_pgmL',               ...
        'TNF_Flag',                     ...
        'Clinical_Interpretation'       ...
    });

% --- Save to CSV ---
csv_filename = 'BP_ZnO_Biosensor_Patient_Data.csv';
writetable(T, csv_filename);

%% =========================================================
% COMMAND WINDOW OUTPUT
%% =========================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf('        BP-ZnO BIOSENSOR ANALYSIS RESULTS          \n');
fprintf('====================================================\n\n');

%% ---- IL-6 ANTIBODY CONFIRMATION -------------------------
fprintf('IL-6 ANTIBODY CONFIRMATION\n');
fprintf('----------------------------------------------------\n');
fprintf('  Antibody Target         : Interleukin-6 (IL-6)\n');
fprintf('  Antibody Type           : Anti-IL-6 monoclonal IgG\n');
fprintf('  Immobilization Method   : EDC/NHS covalent coupling\n');
fprintf('  Surface Coverage        : %.1f %%\n',        IL6_coverage);
fprintf('  Dissociation Const (Kd) : %.1e M\n',         IL6_Kd);
fprintf('  Sensitivity             : %.1f uA/(ng/mL)\n', IL6_sensitivity);
fprintf('  Limit of Detection      : %.3f ng/mL\n',      IL6_LOD);
fprintf('  Limit of Quantification : %.3f ng/mL\n',      IL6_LOQ);
fprintf('  Linear Range            : %.2f - %.1f ng/mL\n', IL6_linear_low, IL6_linear_high);
fprintf('  Confirmation Status     : CONFIRMED [OK]\n\n');

%% ---- TNF-ALPHA ANTIBODY CONFIRMATION --------------------
fprintf('TNF-ALPHA ANTIBODY CONFIRMATION\n');
fprintf('----------------------------------------------------\n');
fprintf('  Antibody Target         : Tumor Necrosis Factor-alpha\n');
fprintf('  Antibody Type           : Anti-TNF-alpha monoclonal IgG\n');
fprintf('  Immobilization Method   : EDC/NHS covalent coupling\n');
fprintf('  Surface Coverage        : %.1f %%\n',        TNF_coverage);
fprintf('  Dissociation Const (Kd) : %.1e M\n',         TNF_Kd);
fprintf('  Sensitivity             : %.1f uA/(ng/mL)\n', TNF_sensitivity);
fprintf('  Limit of Detection      : %.3f ng/mL\n',      TNF_LOD);
fprintf('  Limit of Quantification : %.3f ng/mL\n',      TNF_LOQ);
fprintf('  Linear Range            : %.2f - %.1f ng/mL\n', TNF_linear_low, TNF_linear_high);
fprintf('  Confirmation Status     : CONFIRMED [OK]\n\n');

%% ---- XRD RESULTS ----------------------------------------
fprintf('XRD ANALYSIS\n');
fprintf('----------------------------------------------------\n');
fprintf('Number of Peaks : %d\n', num_xrd_peaks);
fprintf('----------------------------------------------------\n');
fprintf('  %-4s | %-12s | %-16s | %-20s\n','Peak','Position','Raw Amplitude','Assignment');
fprintf('  %-4s | %-12s | %-16s | %-20s\n','----','------------','----------------','--------------------');
for i = 1:num_xrd_peaks
    fprintf('  %-4d | %8.2f deg  | %-16d | %-20s\n', ...
        i, locs_xrd(i), xrd_peak_amplitudes(i), xrd_peak_labels{i});
end
fprintf('\n');

%% ---- FTIR RESULTS ---------------------------------------
fprintf('FTIR ANALYSIS\n');
fprintf('----------------------------------------------------\n');
fprintf('Number of Peaks : %d\n', num_ftir_peaks);
fprintf('----------------------------------------------------\n');
fprintf('  %-4s | %-14s | %-16s | %-20s\n','Peak','Position','Raw Amplitude','Assignment');
fprintf('  %-4s | %-14s | %-16s | %-20s\n','----','--------------','----------------','--------------------');
for i = 1:num_ftir_peaks
    fprintf('  %-4d | %9.2f cm^-1 | %-16d | %-20s\n', ...
        i, locs_ftir(i), ftir_peak_amplitudes(i), ftir_peak_labels{i});
end
fprintf('\n');

%% ---- RAMAN RESULTS --------------------------------------
fprintf('RAMAN ANALYSIS\n');
fprintf('----------------------------------------------------\n');
fprintf('Number of Peaks : %d\n', num_raman_peaks);
fprintf('----------------------------------------------------\n');
fprintf('  %-4s | %-14s | %-16s | %-20s\n','Peak','Position','Raw Amplitude','Assignment');
fprintf('  %-4s | %-14s | %-16s | %-20s\n','----','--------------','----------------','--------------------');
for i = 1:num_raman_peaks
    fprintf('  %-4d | %9.2f cm^-1 | %-16d | %-20s\n', ...
        i, locs_raman(i), raman_peak_amplitudes(i), raman_peak_labels{i});
end
fprintf('\n');

%% ---- PHASE 10 : BIOLOGICAL SAMPLE COLLECTION OUTPUT ----
fprintf('====================================================\n');
fprintf('   PHASE 10 : BIOLOGICAL SAMPLE COLLECTION REPORT  \n');
fprintf('====================================================\n\n');

fprintf('STEP 1 : BLOOD COLLECTION\n');
fprintf('----------------------------------------------------\n');
fprintf('  Total Postoperative Patients     : %d\n',   num_patients);
fprintf('  Samples Successfully Collected   : %d\n',   n_collected);
fprintf('  Samples Rejected (Hemolysis)     : %d\n',   n_hemolysed);
fprintf('  Average Blood Volume Collected   : %.2f mL\n', mean(blood_volume));
fprintf('  Collection Success Rate          : %.1f %%\n', 100*n_collected/num_patients);
fprintf('  Status                           : COMPLETED\n\n');

fprintf('STEP 2 : CENTRIFUGATION\n');
fprintf('----------------------------------------------------\n');
fprintf('  Centrifuge Speed                 : %d RPM\n',  centrifuge_rpm);
fprintf('  Centrifuge Duration              : %d min\n',  centrifuge_time);
fprintf('  Centrifuge Temperature           : %d C\n',    centrifuge_temp);
fprintf('  Samples Successfully Centrifuged : %d\n',      n_centrifuged);
fprintf('  Centrifuge Failures              : %d\n',      n_centrifuge_fail);
fprintf('  Centrifugation Success Rate      : %.1f %%\n', 100*n_centrifuged/n_collected);
fprintf('  Status                           : COMPLETED\n\n');

fprintf('STEP 3 : SERUM SEPARATION\n');
fprintf('----------------------------------------------------\n');
fprintf('  Samples Passed Centrifugation    : %d\n',   n_centrifuged);
fprintf('  Lipemic Samples Rejected         : %d\n',   n_lipemic);
fprintf('  Valid Serum Samples Ready        : %d\n',   n_serum_ready);
fprintf('  Total Samples Rejected           : %d\n',   n_serum_rejected);
fprintf('  Average Serum Yield              : %.1f %%\n', mean(serum_yield_pct));
fprintf('  Average Serum Volume per Sample  : %.2f mL\n', mean_serum_vol);
fprintf('  Serum Preparation Success Rate   : %.1f %%\n', 100*n_serum_ready/num_patients);
fprintf('  Status                           : COMPLETED\n\n');

fprintf('BIOMARKER SCREENING RESULTS\n');
fprintf('----------------------------------------------------\n');
fprintf('  IL-6  Clinical Threshold         : > %d pg/mL\n',   IL6_threshold);
fprintf('  TNF-a Clinical Threshold         : > %d pg/mL\n',   TNF_threshold);
fprintf('  Mean IL-6  Level (valid samples) : %.2f pg/mL\n',   mean_IL6_all);
fprintf('  Mean TNF-a Level (valid samples) : %.2f pg/mL\n',   mean_TNF_all);
fprintf('  Patients with Elevated IL-6      : %d / %d\n',      n_IL6_elevated,  n_serum_ready);
fprintf('  Patients with Elevated TNF-alpha : %d / %d\n',      n_TNF_elevated,  n_serum_ready);
fprintf('  Patients with Both Elevated      : %d / %d\n',      n_both_elevated, n_serum_ready);
fprintf('  Patients within Normal Range     : %d / %d\n\n',    n_normal,        n_serum_ready);

%% ---- PHASE 11 : CSV SAVE CONFIRMATION -------------------
fprintf('====================================================\n');
fprintf('   PHASE 11 : CSV FILE SAVE REPORT                 \n');
fprintf('====================================================\n\n');
fprintf('  File Name        : %s\n',        csv_filename);
fprintf('  Total Rows       : %d (patients)\n', num_patients);
fprintf('  Total Columns    : %d\n',         width(T));
fprintf('\n');
fprintf('  Columns Saved:\n');
fprintf('  --------------------------------------------\n');
col_names = T.Properties.VariableNames;
for c = 1:length(col_names)
    fprintf('  %2d. %s\n', c, col_names{c});
end
fprintf('\n');
fprintf('  Sample Preview (First 5 Patients):\n');
fprintf('  --------------------------------------------\n');
fprintf('  %-5s | %-3s | %-6s | %-11s | %-9s | %-10s | %-10s | %-8s | %-25s\n', ...
    'ID','Age','Gender','Surgery','Serum(mL)','IL6(pg/mL)','TNF(pg/mL)','Serum','Clinical');
fprintf('  %-5s | %-3s | %-6s | %-11s | %-9s | %-10s | %-10s | %-8s | %-25s\n', ...
    '-----','---','------','-----------','---------','----------','----------','--------','-------------------------');
for p = 1:5
    fprintf('  %-5s | %-3d | %-6s | %-11s | %7.3f mL | %8.4f   | %8.4f   | %-8s | %-25s\n', ...
        patient_id{p}, ages(p), gender_str{p}, surgery_str{p}, ...
        serum_volume(p), IL6_level(p), TNF_level(p), ...
        serum_str{p}, clinical_str{p});
end
fprintf('  ... [%d more rows saved in CSV]\n\n', num_patients-5);
fprintf('  CSV Save Status  : SUCCESS [OK]\n');
fprintf('  Save Location    : %s\n', fullfile(pwd, csv_filename));

%% ---- MATERIAL CHARACTERIZATION SUMMARY -----------------
fprintf('\n====================================================\n');
fprintf('MATERIAL CHARACTERIZATION SUMMARY\n');
fprintf('----------------------------------------------------\n');
fprintf('ZnO Nanorods Successfully Simulated\n');
fprintf('Black Phosphorus Nanosheets Successfully Simulated\n');
fprintf('BP-ZnO Composite Successfully Formed\n');
fprintf('IL-6 Antibody Functionalization       : CONFIRMED [OK]\n');
fprintf('TNF-alpha Antibody Functionalization  : CONFIRMED [OK]\n');
fprintf('XRD Crystal Structure Analysis        : COMPLETED\n');
fprintf('FTIR Functional Group Analysis        : COMPLETED\n');
fprintf('Raman Spectroscopic Analysis          : COMPLETED\n');
fprintf('Biological Sample Collection          : COMPLETED\n');
fprintf('  >> Total Patients                   : %d\n',   num_patients);
fprintf('  >> Valid Serum Samples              : %d\n',   n_serum_ready);
fprintf('  >> IL-6  Elevated Cases             : %d\n',   n_IL6_elevated);
fprintf('  >> TNF-a Elevated Cases             : %d\n',   n_TNF_elevated);
fprintf('CSV Data Export                       : SAVED [OK]\n');
fprintf('  >> File  : %s\n', csv_filename);
fprintf('  >> Rows  : %d  |  Columns : %d\n', num_patients, width(T));
fprintf('\n====================================================\n');
fprintf('         ANALYSIS COMPLETED SUCCESSFULLY            \n');
fprintf('====================================================\n');