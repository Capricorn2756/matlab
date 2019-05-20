clear; clc;
fileID = fopen('Fullerene.xyz', 'r');
temp = textscan(fileID, '%f %f %f %f', 1000);
fclose(fileID);
Element = temp{1,1}; % Chemical element number
numberOfPoints = length(Element);
%% Part of gui
prompt = {'Шаг по времени в фс', 'Количество шагов', 'w_min', 'w_max', 'Tao_max', 'Усреднение по времени'};
title = 'Начальные данные';
dims = [1 40];
defaultInput = {'2', '300', '1', '1000', '300', '50'};
answer = inputdlg(prompt, title, dims, defaultInput);
%% Start conditions
Disp = 10^-5;   % сдвиг по координате для нахождения сил в каждой точке
timeStep = str2double(answer{1});   %шаг по времени в фс
FullTime = str2double(answer{2});   %полное количество шагов
wMin = str2double(answer{3});   
wMax = str2double(answer{4});
Tao = str2double(answer{5});
Ttimes = str2double(answer{6}); %сколько значений будем считать для усреднения по времени
mass = 125000; %масса одного атома углерода
fileID = fopen('MDfullereneTermo.xyz', 'w'); %создаем файл для записи результата термолизации
%% координаты фулерена
XYZ = zeros(numberOfPoints,3);
XYZ(:,1) = temp{1,2};
XYZ(:,2) = temp{1,3};
XYZ(:,3) = temp{1,4};
%% Setting waitbar
f = waitbar(0, sprintf('Step - 1 of %.0f', FullTime), 'Name', 'Calculation...', 'CreateCancelBtn',...
    'setappdata(gcbf,''canceling'',1)');
setappdata(f, 'canceling', 0);

%%
V0 = zeros(numberOfPoints, 3);
V0 = randn(size(V0))/10^4;  %матрица случайных скоростей от -10^-4 до 10^-4
fprintf(fileID,'%.0f\nstep 1\n', numberOfPoints);
for i = 1:numberOfPoints
    fprintf(fileID,'%.0f\t%f\t%f\t%f\n',Element(1), XYZ(i,:));  %начальные координаты
end
Forces = forcesFinder(XYZ(:,1), XYZ(:,2), XYZ(:,3), numberOfPoints, Disp);
for i = 2:FullTime
    if getappdata(f,'canceling')
        break
    end
    XYZnew = XYZ + V0*timeStep + Forces./mass*timeStep^2/2;
    Forcesnew = forcesFinder(XYZnew(:,1), XYZnew(:,2), XYZnew(:,3), numberOfPoints, Disp);
    V0 = V0 + timeStep*(Forces./mass + Forcesnew./mass)/2;
    fprintf(fileID,'%.0f\nstep %.0f\n', numberOfPoints, i);
    for j = 1:numberOfPoints
        fprintf(fileID,'%.0f\t%f\t%f\t%f\n', Element(1), XYZnew(j,:));
    end
    XYZ = XYZnew;
    Forces = Forcesnew;
    waitbar(i/FullTime, f, sprintf('Step - %.0f of %.0f', i , FullTime))
end
delete(f); %close waitbar
fclose(fileID); %закрываем файл MDfullereneTermo.xyz


fileID = fopen('MDfullereneVelocity.xyz', 'w'); %создаем файл для записи скоростей
%% Setting waitbar
FullTime = Tao + Ttimes; 
f = waitbar(0, sprintf('Step - 1 of %.0f', FullTime), 'Name', 'Calculation...', ...
    'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
setappdata(f, 'canceling', 0);

%%
%fprintf(fileID,'step 1\n');
for i = 1:numberOfPoints
    fprintf(fileID,'%.8f\t%.8f\t%.8f\n', V0(i,:));  %начальные координаты
end
Forces = forcesFinder(XYZ(:,1), XYZ(:,2), XYZ(:,3), numberOfPoints, Disp);
step = 2; %шаг по врмени для записи скоростей 2 фс
for i = 2:FullTime
    if getappdata(f,'canceling')
        break
    end
    XYZnew = XYZ + V0*step + Forces./mass*step^2/2;
    Forcesnew = forcesFinder(XYZnew(:,1), XYZnew(:,2), XYZnew(:,3), numberOfPoints, Disp);
    V0 = V0 + step*(Forces./mass + Forcesnew./mass)/2;
    %fprintf(fileID,'step %.0f\n', i);
    for j = 1:numberOfPoints
        fprintf(fileID,'%f\t%f\t%f\n', V0(j,:));
    end
    XYZ = XYZnew;
    Forces = Forcesnew;
    waitbar(i/FullTime, f, sprintf('Step - %.0f of %.0f', i , FullTime))
end
delete(f); %close waitbar
fclose(fileID); %закрываем файл MDfullereneVelocity.xyz
fileID = fopen('MDfullereneVelocity.xyz', 'r'); %откываем файл со скоростями
temp = textscan(fileID, '%f %f %f', numberOfPoints*FullTime);
fclose(fileID);
V = zeros(numberOfPoints*FullTime, 3);
V(:,1) = temp{1,1}; %X
V(:,2) = temp{1,2}; %Y
V(:,3) = temp{1,3}; %Z
TaoSteps = (Tao - 10)/2;
A_F = zeros(TaoSteps,1);  %массив фунуций автокорреляции для тао = 10:2:Tao фс
for i = 1:TaoSteps
    for j = 1:Ttimes %усреднение по времени
        V1 = V(1 + (j - 1)*numberOfPoints:j*numberOfPoints,:);
        V2 = V(1 + (j - 1 + 4 + i)*numberOfPoints:(j + 4 + i)*numberOfPoints,:);    %один шаг по времени в файле 2фс, поэтому для Тао=10фс j+4        
        for k = 1:numberOfPoints
            A_F(i) = A_F(i) + sqrt(sum(V1(k,:).^2)*sum(V2(k,:).^2));
        end
    end
end
A_F = A_F./Ttimes;
syms w;
P = 0;
for i = 1:TaoSteps
    P = P + cos(w*(9 + i))*A_F(i);
end
w1= wMin:0.5:wMax;
y1 = subs(P,w,w1);
plot(w1,y1); %построение спектра