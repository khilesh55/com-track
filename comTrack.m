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

%Extract data columns into vector arrays
leftShoulderX = table2array(dataTable(2:end, 16));
leftShoulderY = table2array(dataTable(2:end, 17));
leftShoulderZ = table2array(dataTable(2:end, 18));
rightShoulderX = table2array(dataTable(2:end, 19));
rightShoulderY = table2array(dataTable(2:end, 20));
rightShoulderZ = table2array(dataTable(2:end, 21));
leftHipX = table2array(dataTable(2:end, 34));
leftHipY = table2array(dataTable(2:end, 35));
leftHipZ = table2array(dataTable(2:end, 36));
rightHipX = table2array(dataTable(2:end, 37));
rightHipY = table2array(dataTable(2:end, 38));
rightHipZ = table2array(dataTable(2:end, 39));

%Turn this into a matrix
dataArray = [leftShoulderX leftShoulderY leftShoulderZ...
    rightShoulderX rightShoulderY rightShoulderZ...
    leftHipX leftHipY leftHipZ...
    rightHipX rightHipY rightHipZ];

%Initialize loop vars
i=1; 
val1=0;
val2=0;
val3=0;

%Iterate through the matrix to smooth outliers 
for j=1:12
    for i=2:(dataPoints-2)
        val1 = dataArray((i-1),j);
        val2 = dataArray(i,j);
        val3 = dataArray((i+1),j);
        
        if ((abs((val2-val1)/val2)>1) && (abs((val2-val3)/val2)>1))
            val2 = (val1+val3)/2;
            dataArray(i,j) = val2;
        end
    end
end

%Extract data columns into vector arrays
leftShoulderX = dataArray(:,1);
leftShoulderY = dataArray(:,2);
leftShoulderZ = dataArray(:,3);
rightShoulderX = dataArray(:,4);
rightShoulderY = dataArray(:,5);
rightShoulderZ = dataArray(:,6);
leftHipX = dataArray(:,7);
leftHipY = dataArray(:,8);
leftHipZ = dataArray(:,9);
rightHipX = dataArray(:,10);
rightHipY = dataArray(:,11);
rightHipZ = dataArray(:,12);

%Get trunk coords X and Y
trunkCentreY = (leftShoulderY-rightShoulderY)/2;
trunkCentreX = ((leftShoulderX-leftHipX)+(rightShoulderX-rightHipX))/4;

%For a reference to initial stand trunk height
refNum = 30;
hRef = (leftShoulderX(1:refNum)-leftHipX(1:refNum))+...
    (rightShoulderX(1:refNum)-rightHipX(1:refNum))/2;
tRef = mean(trunkCentreX(1:refNum));

%Estimate lean angle
alpha = acos((((rightShoulderX-rightHipX)+...
    (leftShoulderX-leftHipX))/2)/hRef);

%Get trunk coord Z
trunkCentreZ = ((rightHipY+leftHipY)/2)+tRef*sin(alpha);

%Generate row vector for plotting
pointNumVec = linspace(1,dataPoints-1,dataPoints-1);
      
%Plot figs
figure(1);
title('Trunk X')
xlabel('Data Point')
plot(pointNumVec, trunkCentreX);
ylabel('Trunk X')
      
figure(2);
title('Trunk Y')
xlabel('Data Point')
plot(pointNumVec, trunkCentreY);
ylabel('Trunk Y')

figure(3);
title('Trunk Z')
xlabel('Data Point')
plot(pointNumVec, trunkCentreZ);
ylabel('Trunk Z')

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

% Parameters extract from data table
t = dataTable(:,55);
n = size(t);
f = 1/t;
% Compute FFT of the time domain Trunk data
FFTx = fft(trunkCentreX);
FFTy = fft(trunkCentreY);
FFTz = fft(trunkCentreZ);
% Compute two-sided amplitude spectrum of signal
P2x = abs(FFTx/n);
P1x = P2x(1:n/2+1)
P1x(2:end-1) = 2*P1x(2:end-1);
P2y = abs(FFTy/n);
P1y = P2y(1:n/2+1)
P1y(2:end-1) = 2*P1y(2:end-1);
P2z = abs(FFTz/n);
P1z = P2z(1:n/2+1)
P1z(2:end-1) = 2*P1z(2:end-1);
% Plot single-sided amplitude spectrum of signal
subplot(1,3,1)
plot(f,P1x) 
title('Single-Sided Amplitude Spectrum of X Trunk Centre')
xlabel('f (Hz)')
ylabel('|P1x(f)|')
subplot(1,3,2)
plot(f,P1y) 
title('Single-Sided Amplitude Spectrum of Y Trunk Centre')
xlabel('f (Hz)')
ylabel('|P1y(f)|')
subplot(1,3,3)
plot(f,P1z) 
title('Single-Sided Amplitude Spectrum of Z Trunk Centre')
xlabel('f (Hz)')
ylabel('|P1z(f)|')
% Identify noise frequency range - apply band pass filter
a = fir1(40,2*[1 200]/f); 
b = fir1(48,2*[1 200]/f); 
c = fir1(60,2*[1 200]/f); 
newTCx = filter(a,1,trunkCentreX);
newTCy = filter(b,1,trunkCentreY);
newTCz = filter(c,1,trunkCentreZ);
% Plot filtered signal against original
subplot(1,3,1)
hold on
plot(t,trunkCentreX)
plot(t,newTCx)
hold off
title('X Trunk Centre')
xlabel('t (s)')
ylabel('Displacement')
subplot(1,3,2)
hold on
plot(t,trunkCentreY)
plot(t,newTCy)
hold off
title('Y Trunk Centre')
xlabel('t (s)')
ylabel('Displacement')
subplot(1,3,3)
hold on
plot(t,trunkCentreZ)
plot(t,newTCz)
hold off
title('Z Trunk Centre')
xlabel('t (s)')
ylabel('Displacement')

%% Signal Processing - Output Parameter: Area Under the Curve
% Established Baseline
% Assign index begin and end
% Area under the curve (trapz)
% Assign Off-balance limits
% Boolean trigger for whether or not user is off balance

%% Signal Processing - Output Parameter: Peaks Identification
% Findpeaks point
findpeaks(trunkCentreX);
findpeaks(trunkCentreY);
findpeaks(trunkCentreZ);
% Plot Peaks as function of time
% Area under the curve (trapz)
