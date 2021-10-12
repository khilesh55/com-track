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

%Read data into table from Excel
dataTable = readtable(fullfile(source_dir, source_files(1).name));

%Table size has 2 numbers: # of rows and # of columns
%We only want the # rows, so 
tableSize = size(dataTable);
%size(dataArray) ========= [M, N] <---- tableSize
dataPoints = tableSize(1);

%The following matrix is has headings [US, Euro, UK, in, cm)
%From https://www.candefashions.com/about/shoe-size-conversion-chart/
shoeSizes = [6      39      5.5     9.25        23.5;
             6.5	39      6       9.5         24.1;
             7      40      6.5     9.625       24.4;
             7.5	40.5	7       9.75        24.8;
             8      41      7.5     9.9375      25.4;
             8.5	41.5	8       10.125      25.7;
             9      42      8.5     10.25       26;
             9.5	42.5	9       10.4375     26.7;
             10     43      9.5     10.5625     27;
             10.5	43.5	10      10.75       27.3;
             11     44      10.5	10.9375     27.9;
             11.5	44.5	11      11.125      28.3;
             12     45      11.5	11.25       28.6;
             13     46      12.5	11.5625     29.4;
             14     47      13.5	11.875      30.2;
             15     48      14.5	12.1875     31 ];

%Ask user for the patient's shoe size
shoeSizeCM = 0;
valid = 0;
prompt = "Please enter participant's shoe size: ";
shoeSizeInput = char(inputdlg(prompt));

%Check if the given shoe size exists in the matrix
if any(shoeSizes(:) == str2double(shoeSizeInput))
    valid = 1;
else valid = 0;
end

%If valid skip, if invalid ask again until valid
while (valid == 0)
    prompt = "That is not a valid shoe size. Please try again: ";
    shoeSizeInput = char(inputdlg(prompt));
    if any(shoeSizes(:) == str2double(shoeSizeInput))
        valid = 1;
    else valid = 0;
    end
end

%Dropdown menu for user to choose sizing convention
conventions = {'US','European','UK'};
[idx,tf] = listdlg('PromptString',{'Please select size convention.'},...
    'ListString',conventions,'SelectionMode','single');
sizeIdx = 0;

%Depending on chosen convention and provided shoe size, locate
%shoe size in cm from the matrix
switch idx
    case 1
        sizeIdx = find(shoeSizes(:,1) == str2double(shoeSizeInput));
        shoeSizeCM = shoeSizes(sizeIdx, 5);
    case 2
        sizeIdx = find(shoeSizes(:,2) == str2double(shoeSizeInput));
        shoeSizeCM = shoeSizes(sizeIdx, 5);
    case 3
        sizeIdx = find(shoeSizes(:,3) == str2double(shoeSizeInput));
        shoeSizeCM = shoeSizes(sizeIdx, 5);
end

shoeSizeCM

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

%Plot figs
figure(1);
title('Centre of Mass X')
plot(time, comX);
xlabel('Time (s)')
ylabel('Centre of Mass X')
      
figure(2);
title('Centre of Mass Y')
plot(time, comY);
xlabel('Time (s)')
ylabel('Centre of Mass Y')

figure(3);
title('Centre of Mass Z')
plot(time, comZ);
xlabel('Time (s)')
ylabel('Centre of Mass Z')

figure(7);
plot(time, paramDX);
xlabel('Time (s)')
ylabel('Parameter D')

figure(8);
plot(time, paramHX);
xlabel('Time (s)')
ylabel('Parameter H')

fullArrayX = [leftShoulderX; rightShoulderX; leftHipX; rightHipX; leftAnkleX; rightAnkleX; headX; comX];
fullArrayY = [leftShoulderY; rightShoulderY; leftHipY; rightHipY; leftAnkleY; rightAnkleY; headY; comY];
fullArrayZ = [leftShoulderZ; rightShoulderZ; leftHipZ; rightHipZ; leftAnkleZ; rightAnkleZ; headY; comZ];

for i = 1:length(time)/5
    figure(4);
    plot3(leftShoulderX(i), leftShoulderZ(i), leftShoulderY(i), 'bx', ...
        rightShoulderX(i), rightShoulderZ(i), rightShoulderY(i), 'bx', ...
        leftHipX(i), leftHipZ(i), leftHipY(i), 'bx', ...
        rightHipX(i), rightHipZ(i), rightHipY(i), 'bx', ...
        trunkCentreX(i), trunkCentreZ(i), trunkCentreY(i), 'go', ...
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
end

%% Signal Processing - Noise Removal

% Remove spikes/outliers
D = filloutliers(paramDX,'clip','movmedian',8,'SamplePoints',time);
H = filloutliers(paramHX,'clip','movmedian',6,'SamplePoints',time);

% Variable declaration
t = time;
n = length(t);
Fs = 1/(mean(diff(time)));
Fn = Fs/2;

% Fourier Transform
f = linspace(0,1,fix(n/2)+1)*Fn;
i = 1:length(f);
fftD = fft(paramDX)/n;
fftH = fft(paramHX)/n;
D_mag = abs(fftD(i))*2;
H_mag = abs(fftH(i))*2;

% Low pass filter design
[D1 D2] = butter(3, 0.9, 'low');
[H1 H2] = butter(1, 0.5, 'low');

% Filter implementation
D_filtered = filter(D1,D2,D);
H_filtered = filter(H1,H2,H);

% Convert Filtered signal to Frequency domain
fftDf = fft(D_filtered)/n;
Df_mag = abs(fftDf(i))*2;
fftHf = fft(H_filtered)/n;
Hf_mag = abs(fftHf(i))*2;

% Calaculate baseline
LM_D = islocalmin(D_filtered,'MinSeparation',2,'SamplePoints',t);
LM_H = islocalmin(H_filtered,'MinSeparation',10,'SamplePoints',t);
base_D = median(D_filtered(LM_D));
base_H = median(H_filtered(LM_H));

% Offset signal to new baseline
new_D = D_filtered-base_D;
new_H = H_filtered-base_H;

% Calculate Gradient Response
gradient_D = gradient(new_D);
gradient_H = gradient(new_H);


%% Signal Processing - Output Parameters

% Standing rms

% An estimate of the gradient for rise speed
RS_index = islocalmax(gradient_H,'MinSeparation',10,'SamplePoints',t);
RS = gradient_H(RS_index);
t_RS = t(RS_index);
% An estimate of the gradient for to sit speed
SS_index = islocalmin(gradient_H,'MinSeparation',10,'SamplePoints',t);
SS = gradient_H(SS_index);
t_SS = t(SS_index);

% An estimate of the lean speed(Gradient of +D)
LS = gradient_D(RS_index);
t_LS = t(RS_index);
% An estimate of the gradient for lean back speed
LBS = gradient_D(SS_index);
t_LBS = t(SS_index);

% The distance D of the COM during the transition from sitting to standing
sit2standD = max(new_D,[RS_index(1) SS_index(1)])-min(new_D,[3 RS_index(1)]);
% The distance D of the COM during the transition from standing to sitting
stand2sitD = max(new_D,[RS_index(1) SS_index(1)])-min(new_D,[SS_index(1) RS_index(2)]);
% The height H of the COM during the transition from sitting to standing
sit2standH = max(new_H,[RS_index(1) SS_index(1)])-min(new_H,[3 RS_index(1)]);
% The height H of the COM during the transition from standing to sitting
stand2sitD = max(new_H,[RS_index(1) SS_index(1)])-min(new_H,[SS_index(1) RS_index(2)]);
% The average sit to stand sequence duration (time per sit to stand cycle/ number of cycles)
index = islocalmax(gradient_H,'MinSeparation',6,'SamplePoints',t);
Ts = mean(diff(t(index)));
%% Plot Results

% Display original and filtered signal in time domain
figure(1)
subplot(2,1,1)
plot(t,paramDX,t,new_D,t_LS,new_D(RS_index),'*r',t_LBS,new_D(SS_index),'ob')
title('D Response - Time Domain')
xlabel('Time(s)') 
ylabel('D Position(m)') 
legend({'Original','Filtered'},'Location','southeast')
subplot(2,1,2)
plot(t,paramHX,t,new_H,t_RS,new_H(RS_index),'*r',t_SS,new_H(SS_index),'ob')
title('H Response - Time Domain')
xlabel('Time(s)') 
ylabel('H Position(m)') 
legend({'Original','Filtered'},'Location','southeast')

% Display original and filtered signal in frequency domain
figure(2)
subplot(2,1,1)
plot(f,D_mag,f,Df_mag)
title('D Response - Single Sided Frequency')
xlabel('Frequency(Hz)') 
ylabel('D Magnitude(dB)') 
legend({'Original','Filtered'},'Location','northeast')
subplot(2,1,2)
plot(f,H_mag,f,Hf_mag)
title('H Response - Single Sided Frequency')
xlabel('Frequency(Hz)') 
ylabel('H Magnitude(dB)') 
legend({'Original','Filtered'},'Location','northeast')

% Plot Gradient Response in time domain
figure(3)
subplot(2,1,1)
plot(t,gradient_D,t_LS,LS,'*r',t_LBS,LBS,'ob')
ylim([-1 1]);
xlabel('Time(s)') 
ylabel('Gradient D(m/s)') 
subplot(2,1,2)
plot(t,gradient_H,t_RS,RS,'*r',t_SS,SS,'ob')
ylim([-0.1 0.1]);
xlabel('Time(s)') 
ylabel('Gradient H(m/s)') 

%% Export Outputs
 

