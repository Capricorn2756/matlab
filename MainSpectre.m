clear; clc;
fileID = fopen('Fullerene.xyz', 'r');
temp = textscan(fileID, '%f %f %f %f', 1000);
fclose(fileID);
Element = temp{1,1}; % Chemical element number
numberOfPoints = length(Element);
%% Part of gui
prompt = {'Øàã ïî âðåìåíè â ôñ', 'Êîëè÷åñòâî øàãîâ', 'w_min', 'w_max', 'Tao_max', 'Óñðåäíåíèå ïî âðåìåíè'};
title = 'Íà÷àëüíûå äàííûå';
dims = [1 40];
defaultInput = {'2', '300', '1', '1000', '300', '50'};
answer = inputdlg(prompt, title, dims, defaultInput);
%% Start conditions
Disp = 10^-5;   % ñäâèã ïî êîîðäèíàòå äëÿ íàõîæäåíèÿ ñèë â êàæäîé òî÷êå
timeStep = str2double(answer{1});   %øàã ïî âðåìåíè â ôñ
FullTime = str2double(answer{2});   %ïîëíîå êîëè÷åñòâî øàãîâ
wMin = str2double(answer{3});   
wMax = str2double(answer{4});
Tao = str2double(answer{5});
Ttimes = str2double(answer{6}); %ñêîëüêî çíà÷åíèé áóäåì ñ÷èòàòü äëÿ óñðåäíåíèÿ ïî âðåìåíè
mass = 125000; %ìàññà îäíîãî àòîìà óãëåðîäà
fileID = fopen('MDfullereneTermo.xyz', 'w'); %ñîçäàåì ôàéë äëÿ çàïèñè ðåçóëüòàòà òåðìîëèçàöèè
%% êîîðäèíàòû ôóëåðåíà
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
V0 = randn(size(V0))/10^4;  %ìàòðèöà ñëó÷àéíûõ ñêîðîñòåé îò -10^-4 äî 10^-4
fprintf(fileID,'%.0f\nstep 1\n', numberOfPoints);
for i = 1:numberOfPoints
    fprintf(fileID,'%.0f\t%f\t%f\t%f\n',Element(1), XYZ(i,:));  %íà÷àëüíûå êîîðäèíàòû
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
fclose(fileID); %çàêðûâàåì ôàéë MDfullereneTermo.xyz


fileID = fopen('MDfullereneVelocity.xyz', 'w'); %ñîçäàåì ôàéë äëÿ çàïèñè ñêîðîñòåé
%% Setting waitbar
FullTime = Tao + Ttimes; 
f = waitbar(0, sprintf('Step - 1 of %.0f', FullTime), 'Name', 'Calculation...', ...
    'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
setappdata(f, 'canceling', 0);

%%
%fprintf(fileID,'step 1\n');
for i = 1:numberOfPoints
    fprintf(fileID,'%.8f\t%.8f\t%.8f\n', V0(i,:));  %íà÷àëüíûå êîîðäèíàòû
end
Forces = forcesFinder(XYZ(:,1), XYZ(:,2), XYZ(:,3), numberOfPoints, Disp);
step = 2; %øàã ïî âðìåíè äëÿ çàïèñè ñêîðîñòåé 2 ôñ
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
fclose(fileID); %çàêðûâàåì ôàéë MDfullereneVelocity.xyz
fileID = fopen('MDfullereneVelocity.xyz', 'r'); %îòêûâàåì ôàéë ñî ñêîðîñòÿìè
temp = textscan(fileID, '%f %f %f', numberOfPoints*FullTime);
fclose(fileID);
V = zeros(numberOfPoints*FullTime, 3);
V(:,1) = temp{1,1}; %X
V(:,2) = temp{1,2}; %Y
V(:,3) = temp{1,3}; %Z
TaoSteps = (Tao - 10)/2;
A_F = zeros(TaoSteps,1);  %ìàññèâ ôóíóöèé àâòîêîððåëÿöèè äëÿ òàî = 10:2:Tao ôñ
for i = 1:TaoSteps
    for j = 1:Ttimes %óñðåäíåíèå ïî âðåìåíè
        V1 = V(1 + (j - 1)*numberOfPoints:j*numberOfPoints,:);
        V2 = V(1 + (j - 1 + 4 + i)*numberOfPoints:(j + 4 + i)*numberOfPoints,:);    %îäèí øàã ïî âðåìåíè â ôàéëå 2ôñ, ïîýòîìó äëÿ Òàî=10ôñ j+4        
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
plot(w1,y1); %ïîñòðîåíèå ñïåêòðà
