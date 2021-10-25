%% comtrack.m
%% Khilesh Jairamdas August 9, 2021
%% This script loads data from a Kinetisense Excel output into arrays
%% It extracts body landmark coords into vectors
%% It performs calculations to get and plot the COM position over time

clc
clearvars

%Source and output directory
source_dir = 'C:\Users\VR-7\Documents\MATLAB\Example Patient Data'; 
output_dir = 'C:\Users\VR-7\Documents\MATLAB\Example Patient Output';
source_files = dir(fullfile(source_dir, '*.xlsx'));
fileIdx = 1; % Select which file to open

%Read data into table from Excel
dataTable = readtable(fullfile(source_dir, source_files(fileIdx).name));

%Table size has 2 numbers: # of rows and # of columns
%We only want the # rows, so 
tableSize = size(dataTable);
%size(dataArray) ========= [M, N] <---- tableSize
dataPoints = tableSize(1);

prompt = "Please enter centre of mass offset from trunk centre: ";
offsetInput = char(inputdlg(prompt));
offsetDouble = str2double(offsetInput);

if ~isnan(offsetDouble)
    valid = 1;
else
    valid = 0;
end

while (valid == 0)
    prompt = "That is not a valid value for mass offset. Please try again: ";
    offsetInput = char(inputdlg(prompt));
    offsetDouble = str2double(offsetInput);
    if ~isnan(offsetDouble)
        valid = 1;
    else
        valid = 0;
    end
end

%Dropdown menu for user to choose unit convention
conventions = {'cm','mm','in'};
[idx,tf] = listdlg('PromptString',{'Please select unit convention.'},...
    'ListString',conventions,'SelectionMode','single');

%Depending on chosen convention and provided shoe size, calculate
%offset in m
switch idx
    case 1
        offsetVal = offsetDouble / 100;
    case 2
        offsetVal = offsetDouble / 1000;
    case 3
        offsetVal = offsetDouble * 25.4 / 1000;
end


%Extract data columns into vector arrays
leftShoulderX = dataTable.ShoulderLeftX;
leftShoulderY = dataTable.ShoulderLeftY;
leftShoulderZ = dataTable.ShoulderLeftZ;
rightShoulderX = dataTable.ShoulderRightX;
rightShoulderY = dataTable.ShoulderRightY;
rightShoulderZ = dataTable.ShoulderRightZ;
leftHipX = dataTable.HipLeftX;
leftHipY = dataTable.HipLeftY;
leftHipZ = dataTable.HipLeftZ;
rightHipX = dataTable.HipRightX;
rightHipY = dataTable.HipRightY;
rightHipZ = dataTable.HipRightZ;
leftAnkleX = dataTable.AnkleLeftX;
leftAnkleY = dataTable.AnkleLeftY;
leftAnkleZ = dataTable.AnkleLeftZ;
rightAnkleX = dataTable.AnkleRightX;
rightAnkleY = dataTable.AnkleRightY;
rightAnkleZ = dataTable.AnkleRightZ;
timeStamp = dataTable.Timestamp;

%
headX = dataTable.HeadX;
headY = dataTable.HeadY;
headZ = dataTable.HeadZ;
%

[yr, mth, day, hr, mn, s] = datevec(timeStamp);
timeAbs = 3600*hr + 60*mn + s;
time = timeAbs - timeAbs(1);

trunkCentreX = mean([leftShoulderX, rightShoulderX, leftHipX, rightHipX], 2);
trunkCentreY = mean([leftShoulderY, rightShoulderY, leftHipY, rightHipY], 2);
trunkCentreZ = mean([leftShoulderZ, rightShoulderZ, leftHipZ, rightHipZ], 2);

midShoulderX = (leftShoulderX + rightShoulderX)/2;
midShoulderY = (leftShoulderY + rightShoulderY)/2;
midShoulderZ = (leftShoulderZ + rightShoulderZ)/2;
midHipX = (leftHipX + rightHipX)/2;
midHipY = (leftHipY + rightHipY)/2;
midHipZ = (leftHipZ + rightHipZ)/2;

midAnkleX = (leftAnkleX + rightAnkleX)/2;
midAnkleY = (leftAnkleY + rightAnkleY)/2;
midAnkleZ = (leftAnkleZ + rightAnkleZ)/2;

crossVector = [midShoulderX - midHipX, midShoulderY - midHipY, midShoulderZ - midHipZ];
crossVectorMag = vecnorm(crossVector, 2, 2);
crossUnitVector = crossVector ./ crossVectorMag;

comOffsetVector = crossUnitVector * offsetVal;

comX = trunkCentreX + comOffsetVector(:, 1);
comY = trunkCentreY + comOffsetVector(:, 2);
comZ = trunkCentreZ + comOffsetVector(:, 3);

paramDX = comZ - midAnkleZ;
paramHX = comY - midAnkleY;

fullArrayX = [leftShoulderX; rightShoulderX; leftHipX; rightHipX; leftAnkleX; rightAnkleX; headX; comX];
fullArrayY = [leftShoulderY; rightShoulderY; leftHipY; rightHipY; leftAnkleY; rightAnkleY; headY; comY];
fullArrayZ = [leftShoulderZ; rightShoulderZ; leftHipZ; rightHipZ; leftAnkleZ; rightAnkleZ; headY; comZ];

for i = 1:length(time)
    animplot = figure(4);
    plot3(leftShoulderX(i), leftShoulderZ(i), leftShoulderY(i), 'bx', ...
        rightShoulderX(i), rightShoulderZ(i), rightShoulderY(i), 'bx', ...
        leftHipX(i), leftHipZ(i), leftHipY(i), 'bx', ...
        rightHipX(i), rightHipZ(i), rightHipY(i), 'bx', ...
        trunkCentreX(i), trunkCentreZ(i), trunkCentreY(i), 'gx', ...
        midShoulderX(i), midShoulderZ(i), midShoulderY(i), 'rx', ...
        midHipX(i), midHipZ(i), midHipY(i), 'rx', ...
        comX(i), comZ(i), comY(i), 'ro', ...
        leftAnkleX(i), leftAnkleZ(i), leftAnkleY(i), 'bx', ...
        rightAnkleX(i), rightAnkleZ(i), rightAnkleY(i), 'bx', ...
        headX(i), headZ(i), headY(i), 'bo', 'MarkerSize', 10);
    xlabel('x');
    ylabel('z');
    zlabel('y');
    xlim([min(fullArrayX), max(fullArrayX)]);
    ylim([min(fullArrayZ), max(fullArrayZ)]);
    zlim([min(fullArrayY), max(fullArrayY)]);
    title(strcat({'Frame: '}, string(i)));
    %Change plot view
    view(-90, 0); % Front view (0, 0); Side View (90, 0) or (-90, 0); Comment out line for isometric
    %GIF writer
    gifFile = string(strcat(output_dir, {'\'}, source_files(fileIdx).name, {'_animplot3.gif'}));
    im = frame2im(getframe(animplot));
    [A,map] = rgb2ind(im, 256);
    if i == 1
        imwrite(A,map,gifFile,'gif','LoopCount',Inf,'Delay',1/60);
    else
        imwrite(A,map,gifFile,'gif','WriteMode','append','Delay',1/60);
    end
end

%% Signal Processing - Output Parameters

% An estimate of the lean speeds(Gradient of +D)
Lean_Speed = P_PeakD;
t_Lean_Speed = t(P_PeakD_index);
% An estimate of the lean back speeds (Gradient of -D)
Lean_Back_Speed = N_PeakD;
t_Lean_Back_Speed = t(N_PeakD_index);
% An estimate of the rise speed (Gradient of +H)
Rise_Speed = P_PeakH;
t_Rise_Speed = t(P_PeakH_index);
% An estimate of the sit speed (Gradient of -H)
Sit_Speed = N_PeakH;
t_Sit_Speed = t(N_PeakH_index);
% Lean-in to Lean-back (+D) Step Height
P_delta_D = abs(mean(P_FlatD)-mean(N_FlatD));
% Lean-back to Lean-in (-D) Step Height
N_delta_D = abs(mean(N_FlatD)-mean(P_FlatD));
% Sit to Stand (+H) Step Height
P_delta_H = abs(mean(P_FlatH)-mean(N_FlatH));
% Stand to Sit (-H) Step Height
N_delta_H = abs(mean(N_FlatH)-mean(P_FlatH));
% Leaning Variation (+D rms)
Leaning_rms = rms(P_FlatD);
% Leaning-back Variation (-D rms)
Leaning_back_rms = rms(N_FlatD);
% Standing Variation (+H rms)
Standing_rms = rms(P_FlatH);
% Siting Variation (-H)
Sitting_rms = rms(N_FlatH);

% Cycling Period
Ts_D = diff(t(P_PeakD_index));
Ts_H = diff(t(P_PeakH_index));

%% Plot Results

% Display original and filtered signal in time domain
figure(1)
subplot(2,1,1)
hold on
plot(t,paramDX,t,new_D,t_Lean_Speed,new_D(P_PeakD_index),'*r',t_Lean_Back_Speed,new_D(N_PeakD_index),'ob')
plot(t(FlatD_index),FlatD,'.g')
yline(base_D,'--k','Baseline');
hold off
ylim([-1 1]);
title('D Response - Time Domain')
xlabel('Time(s)') 
ylabel('D Position(m)') 
legend({'Original','Filtered'},'Location','southeast')
subplot(2,1,2)
hold on
plot(t,paramHX,t,new_H,t_Rise_Speed,new_H(P_PeakH_index),'*r',t_Sit_Speed,new_H(N_PeakH_index),'ob')
plot(t(FlatH_index),FlatH,'.g')
yline(base_H,'--k','Baseline');
hold off
ylim([-0.5 1.5]);
title('H Response - Time Domain')
xlabel('Time(s)') 
ylabel('H Position(m)') 
legend({'Original','Filtered'},'Location','southeast')

% Display original and filtered signal in frequency domain
figure(2)
subplot(2,1,1)
hold on
plot(f,D_mag,f,Df_mag)
xline(fc_D,'--k','Cut-off Frequency')
hold off
title('D Response - Single Sided Frequency')
xlabel('Frequency(Hz)') 
ylabel('D Magnitude(dB)') 
legend({'Original','Filtered'},'Location','northeast')
subplot(2,1,2)
hold on
plot(f,H_mag,f,Hf_mag)
xline(fc_H,'--k','Cut-off Frequency')
hold off
title('H Response - Single Sided Frequency')
xlabel('Frequency(Hz)') 
ylabel('H Magnitude(dB)') 
legend({'Original','Filtered'},'Location','northeast')

% Plot Gradient Response in time domain
figure(3)
subplot(2,1,1)
plot(t,gradient_D,t_Lean_Speed,Lean_Speed,'*r',t_Lean_Back_Speed,Lean_Back_Speed,'ob')
ylim([-0.1 0.1]);
title('D gradient Response - Time Domain')
xlabel('Time(s)') 
ylabel('Gradient D(m/s)') 
subplot(2,1,2)
plot(t,gradient_H,t_Rise_Speed,Rise_Speed,'*r',t_Sit_Speed,Sit_Speed,'ob')
ylim([-0.1 0.1]);
title('H gradient Response - Time Domain')
xlabel('Time(s)') 
ylabel('Gradient H(m/s)') 

%% Export Outputs

[~, dataSetName, ~] = fileparts(fullfile(source_dir, source_files(fileIdx).name));
outputWrite = string(strcat(output_dir, {'\'}, dataSetName, '_Results.xlsx'));

xlsTitles = ["leanSpeed", "leanSpeedTimestamp", "leanBackSpeed", "leanBackSpeedTimestamp", ...
    "riseSpeed", "riseSpeedTimestamp", "dropSpeed", "dropSpeedTimestamp", "leanDistance", ...
    "leanBackDistance", "riseHeight", "dropHeight", "leanRMS", "leanBackRMS", "riseRMS", ...
    "dropRMS", "cyclePeriod"];

xlswrite(outputWrite, xlsTitles, 1, 'A1');

xlswrite(outputWrite, Lean_Speed, 1, 'A2');
xlswrite(outputWrite, t_Lean_Speed, 1, 'B2');
xlswrite(outputWrite, Lean_Back_Speed, 1, 'C2');
xlswrite(outputWrite, t_Lean_Back_Speed, 1, 'D2');
xlswrite(outputWrite, Rise_Speed, 1, 'E2');
xlswrite(outputWrite, t_Rise_Speed, 1, 'F2');
xlswrite(outputWrite, Sit_Speed, 1, 'G2');
xlswrite(outputWrite, t_Sit_Speed, 1, 'H2');
xlswrite(outputWrite, P_delta_D, 1, 'I2');
xlswrite(outputWrite, N_delta_D, 1, 'J2');
xlswrite(outputWrite, P_delta_H, 1, 'K2');
xlswrite(outputWrite, N_delta_H, 1, 'L2');
xlswrite(outputWrite, Leaning_rms, 1, 'M2');
xlswrite(outputWrite, Leaning_back_rms, 1, 'N2');
xlswrite(outputWrite, Standing_rms, 1, 'O2');
xlswrite(outputWrite, Sitting_rms, 1, 'P2');
xlswrite(outputWrite, Ts_D, 1, 'Q2');
xlswrite(outputWrite, Ts_H, 1, 'R2');

saveas(h(1), string(strcat(output_dir, {'\'}, dataSetName, '_figure1.jpg')));
saveas(h(2), string(strcat(output_dir, {'\'}, dataSetName, '_figure2.jpg')));
saveas(h(3), string(strcat(output_dir, {'\'}, dataSetName, '_figure3.jpg')));

