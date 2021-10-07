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

%Start GitHub - Khilesh

%Alex - slot in the centroid code instead of my approach

%Ella - look up distance from ankle center to back of heel

%Alex - implement above to create bounding box & alert for if COM
% goes outside of bounding box

%James - extend lookup table to include foot size (extrapolation/ calcs
% may be required

%Check axes - Khilesh to verify (Kinetisense literature)

%Uncertainty radius for COM (from Martin + Emily)
 
%% Signal Processing - Noise Removal

% Remove spikes/outliers
CMx = filloutliers(comX,'clip','movmedian',5,'SamplePoints',time);
CMy = filloutliers(comY,'clip','movmedian',5,'SamplePoints',time);
CMz = filloutliers(comZ,'clip','movmedian',5,'SamplePoints',time);
% Variable declaration
t = time;
n = length(t);
Fs = 1/(mean(diff(time)));
Fn = Fs/2;
f = Fs/n*(1:n);
x_mag = abs(fft(CMx));
y_mag = abs(fft(CMy));
z_mag = abs(fft(CMz));
% Low pass filter design
[b a] = butter(1, 0.3, 'low');
x_filtered = filter(b,a,CMx);
[d c] = butter(1, 0.3, 'low');
y_filtered = filter(d,c,CMy);
[v e] = butter(1, 0.3, 'low');
z_filtered = filter(v,e,CMz);
% Offset signal to new baseline
new_x = x_filtered-mean(x_filtered);
new_y = y_filtered+mean(y_filtered)-0.1-offsetVal;
new_z = z_filtered-mean(z_filtered)+0.1;
% Display original and filtered signal in frequency domain
figure(5)
subplot(3,3,1)
xf_mag = abs(fft(x_filtered));
plot(f,abs(x_mag),f,abs(xf_mag))
title('x Response - Frequency Domain')
xlabel('Frequency(Hz)') 
ylabel('x Magnitude(dB)') 
legend({'Original','Filtered'},'Location','northeast')
subplot(3,3,4)
yf_mag = abs(fft(y_filtered));
plot(f,abs(y_mag),f,abs(yf_mag))
title('y Response - Frequency Domain')
xlabel('Frequency(Hz)') 
ylabel('y Magnitude(dB)') 
legend({'Original','Filtered'},'Location','northeast')
subplot(3,3,7)
zf_mag = abs(fft(z_filtered));
plot(f,abs(z_mag),f,abs(zf_mag))
title('z Response - Frequency Domain')
xlabel('Frequency(Hz)') 
ylabel('z Magnitude(dB)') 
legend({'Original','Filtered'},'Location','northeast')
% Display original and filtered signal in time domain
subplot(3,3,2)
plot(t,comX,t,new_x)
title('x Response - Time Domain')
xlabel('Time(s)') 
ylabel('x Position(m)') 
legend({'Original','Filtered'},'Location','southeast')
subplot(3,3,5)
plot(t,comY,t,new_y)
title('y Response - Time Domain')
xlabel('Time(s)') 
ylabel('y Position(m)') 
legend({'Original','Filtered'},'Location','southeast')
subplot(3,3,8)
plot(t,comZ,t,new_z)
title('z Response - Time Domain')
xlabel('Time(s)') 
ylabel('z Position(m)') 
legend({'Original','Filtered'},'Location','southeast')

%% Signal Processing - Output Parameters

% Area under curve
area_x = trapz(t,new_x)
area_y = trapz(t,new_y)
area_z = trapz(t,new_z)
% Plot area under curve
subplot(3,3,3)
area(t,new_x)
ylim([-0.05 0.05]);
title('Area Under Curve - xt Response')
xlabel('Time(s)') 
ylabel('x Position(m)') 
subplot(3,3,6)
area(t,new_y)
ylim([-0.2 0.6]);
title('Area Under Curve - yt Response')
xlabel('Time(s)') 
ylabel('y Position(m)') 
subplot(3,3,9)
area(t,new_z)
ylim([-1 1]);
title('Area Under Curve - zt Response')
xlabel('Time(s)') 
ylabel('z Position(m)') 
% Calculate rms
rms_x = rms(new_x)
rms_y = rms(new_y)
rms_z = rms(new_z)
% Calculate Gradient
gradient_x = gradient(new_x);
gradient_y = gradient(new_y);
gradient_z = gradient(new_z);
figure(6)
subplot(3,1,1)
plot(t,gradient_x)
subplot(3,1,2)
plot(t,gradient_y)
subplot(3,1,3)
plot(t,gradient_z)

%% Signal Processing - Output Parameter: Area Under the Curve
% Established Baseline
% Assign index begin and end
% Area under the curve (trapz)
% Assign Off-balance limits
% Boolean trigger for whether or not user is off balance
