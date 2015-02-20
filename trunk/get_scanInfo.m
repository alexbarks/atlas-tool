clc, clear, close all

[num,txt,raw] = xlsread('C:\1_Chicago\Data\MIMICS\17_Followup\controls_v2','C79:C136')

size(txt)

directory = 'L:\cv_mri\Aorta-4D_Flow\controls'
%directory = 'L:\cv_mri\Aorta-4D_Flow\BAV'

for n = 1:size(txt,1)
    n
%path = [directory '\' txt(n) '\3dpc']
path = [directory '\' char(txt(n)) '\3dpc']
cd(path)
filename = dir('scanInfo_*.xls');
[NUM,TXT,RAW]=xlsread(filename(1).name,'scanInfo');
sex1 = RAW(2,2);
sex(n,1) = sex1{1};
age1 = RAW(4,2);
age(n,1) = str2num(age1{1}(2:3));
resolution = RAW(13,2);
resolution1(n,1) = resolution{1};
resolution = RAW(14,2);
resolution2(n,1) = resolution{1};
resolution = RAW(15,2);
resolution3(n,1) = resolution{1};
FOV = RAW(11,2);
FOV1(n,1) = FOV{1};
FOV = RAW(12,2);
FOV2(n,1) = FOV{1};
FOV3(n,1) = resolution{1} * 30;
TR(n,1) = cell2mat(RAW(16,2))/8;
TE(n,1) = cell2mat(RAW(17,2));
temp_res(n,1) = TR(n,1)*8;
time_frames(n,1) = cell2mat(RAW(18,2));
flip_angle(n,1) = cell2mat(RAW(19,2));
Venc(n,1) = cell2mat(RAW(22,2));
MRsystem{n} = cell2mat(RAW(21,2));
cd ..
end

% [num,txt,raw] = xlsread('C:\1_Chicago\Data\MIMICS\17_Followup\Followups_v4','B71:B118')
% 
% size(txt)
% 
% directory = 'L:\cv_mri\Aorta-4D_Flow\BAV'
% 
% size1 = n
% pause
% for n = 1:size(txt,1)
%     m = n+size1
% %path = [directory '\' txt(n) '\3dpc']
% path = [directory '\' char(txt(n)) '\3dpc']
% cd(path)
% filename = dir('scanInfo_*.xls');
% [NUM,TXT,RAW]=xlsread(filename(1).name,'scanInfo');
% sex1 = RAW(2,2);
% sex(m,1) = sex1{1};
% age1 = RAW(4,2);
% age(m,1) = str2num(age1{1}(2:3));
% resolution = RAW(13,2);
% resolution1(m,1) = resolution{1};
% resolution = RAW(14,2);
% resolution2(m,1) = resolution{1};
% resolution = RAW(15,2);
% resolution3(m,1) = resolution{1};
% FOV = RAW(11,2);
% FOV1(m,1) = FOV{1};
% FOV = RAW(12,2);
% FOV2(m,1) = FOV{1};
% FOV3(m,1) = resolution{1} * 30;
% TR(m,1) = cell2mat(RAW(16,2))/8;
% TE(m,1) = cell2mat(RAW(17,2));
% temp_res(m,1) = TR(m,1)*8;
% time_frames(m,1) = cell2mat(RAW(18,2));
% flip_angle(m,1) = cell2mat(RAW(19,2));
% Venc(m,1) = cell2mat(RAW(22,2));
% MRsystem{m} = cell2mat(RAW(21,2));
% cd ..
% end
% 
% [num,txt,raw] = xlsread('C:\1_Chicago\Data\MIMICS\17_Followup\Followups_v4','B119:B148')
% 
% size(txt)
% 
% directory = 'L:\cv_mri\Aorta-4D_Flow\Aneurysm'
% 
% size1 = m
% pause
% for n = 1:size(txt,1)
%     m = n+size1
% %path = [directory '\' txt(n) '\3dpc']
% path = [directory '\' char(txt(n)) '\3dpc']
% cd(path)
% filename = dir('scanInfo_*.xls');
% [NUM,TXT,RAW]=xlsread(filename(1).name,'scanInfo');
% sex1 = RAW(2,2);
% sex(m,1) = sex1{1};
% age1 = RAW(4,2);
% age(m,1) = str2num(age1{1}(2:3));
% resolution = RAW(13,2);
% resolution1(m,1) = resolution{1};
% resolution = RAW(14,2);
% resolution2(m,1) = resolution{1};
% resolution = RAW(15,2);
% resolution3(m,1) = resolution{1};
% FOV = RAW(11,2);
% FOV1(m,1) = FOV{1};
% FOV = RAW(12,2);
% FOV2(m,1) = FOV{1};
% FOV3(m,1) = resolution{1} * 30;
% TR(m,1) = cell2mat(RAW(16,2))/8;
% TE(m,1) = cell2mat(RAW(17,2));
% temp_res(m,1) = TR(m,1)*8;
% time_frames(m,1) = cell2mat(RAW(18,2));
% flip_angle(m,1) = cell2mat(RAW(19,2));
% Venc(m,1) = cell2mat(RAW(22,2));
% MRsystem{m} = cell2mat(RAW(21,2));
% cd ..
% end

disp(['Age = ' num2str(mean(age)) ' +/- ' num2str(std(age)) '(range:' num2str(min(age)) '-' num2str(max(age)) 'years old)'])
[I,J] = find(sex == 'M');
disp(['Number of males = ' num2str(size(I,1)) ', number of females = ' num2str(size(txt,1)-size(I,1))])
disp(['Resolution 1 = ' num2str(min(resolution1)) '-' num2str(max(resolution1))])
disp(['Resolution 2 = ' num2str(min(resolution2)) '-' num2str(max(resolution2))])
disp(['Resolution 3 = ' num2str(min(resolution3)) '-' num2str(max(resolution3))])
disp(['FOV 1 = ' num2str(min(FOV1)) '-' num2str(max(FOV1))])
disp(['FOV 2 = ' num2str(min(FOV2)) '-' num2str(max(FOV2))])
disp(['FOV 3 = ' num2str(min(FOV3)) '-' num2str(max(FOV3))])
disp(['TE = ' num2str(min(TE)) '-' num2str(max(TE))])
disp(['TR = ' num2str(min(TR)) '-' num2str(max(TR))])
disp(['TempRes = ' num2str(min(temp_res)) '-' num2str(max(temp_res))])
disp(['TimeFrames = ' num2str(min(time_frames)) '-' num2str(max(time_frames))])
disp(['FlipAngle = ' num2str(min(flip_angle)) '-' num2str(max(flip_angle))])
disp(['Venc = ' num2str(min(Venc)) '-' num2str(max(Venc))])
disp(['MRsystems used: = ' unique(MRsystem)])