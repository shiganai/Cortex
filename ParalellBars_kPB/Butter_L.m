clear all

Bar_Data_All = (dlmread('Trimmed_BarL1Hirabayashi1.trc','',6,1))/1000;%単位変換[mm]→[m]

Data_Time = Bar_Data_All(:,1) * 1000;

BarR0 = Bar_Data_All(:, 2:4);
BarR50 = Bar_Data_All(:, 5:7);
BarR100 = Bar_Data_All(:, 8:10);
BarR130 = Bar_Data_All(:, 11:13);
BarR180 = Bar_Data_All(:, 14:16);
BarR230 = Bar_Data_All(:, 17:19);

BarL0 = Bar_Data_All(:, 20:22);
BarL50 = Bar_Data_All(:, 23:25);
BarL100 = Bar_Data_All(:, 26:28);
BarL130 = Bar_Data_All(:, 29:31);
BarL180 = Bar_Data_All(:, 32:34);
BarL230 = Bar_Data_All(:, 35:37);

Axe_X_UnitVec = mean(BarR0 - BarR230, 1);
Axe_X_UnitVec(3) = 0;
Axe_X_UnitVec = Axe_X_UnitVec ./ norm(Axe_X_UnitVec);
Axe_Y_UnitVec = [0, 0, 1];

BarR0_2D = Get_2D_Coordinate(BarR0, Axe_X_UnitVec, Axe_Y_UnitVec);
BarR50_2D = Get_2D_Coordinate(BarR50, Axe_X_UnitVec, Axe_Y_UnitVec);
BarR100_2D = Get_2D_Coordinate(BarR100, Axe_X_UnitVec, Axe_Y_UnitVec);
BarR130_2D = Get_2D_Coordinate(BarR130, Axe_X_UnitVec, Axe_Y_UnitVec);
BarR180_2D = Get_2D_Coordinate(BarR180, Axe_X_UnitVec, Axe_Y_UnitVec);
BarR230_2D = Get_2D_Coordinate(BarR230, Axe_X_UnitVec, Axe_Y_UnitVec);

BarL0_2D = Get_2D_Coordinate(BarL0, Axe_X_UnitVec, Axe_Y_UnitVec);
BarL50_2D = Get_2D_Coordinate(BarL50, Axe_X_UnitVec, Axe_Y_UnitVec);
BarL100_2D = Get_2D_Coordinate(BarL100, Axe_X_UnitVec, Axe_Y_UnitVec);
BarL130_2D = Get_2D_Coordinate(BarL130, Axe_X_UnitVec, Axe_Y_UnitVec);
BarL180_2D = Get_2D_Coordinate(BarL180, Axe_X_UnitVec, Axe_Y_UnitVec);
BarL230_2D = Get_2D_Coordinate(BarL230, Axe_X_UnitVec, Axe_Y_UnitVec);

rPB = zeros(size(Data_Time));
linePBs = zeros(size(Data_Time, 1), 100);

figure(1)
clf('reset')
hold on
for ii = 1:size(Data_Time, 1)
    linePB = spline([BarL0_2D(ii,1), BarL50_2D(ii,1), BarL100_2D(ii,1), BarL130_2D(ii,1), BarL180_2D(ii,1), BarL230_2D(ii,1)],...
        [BarL0_2D(ii,2), BarL50_2D(ii,2), BarL100_2D(ii,2), BarL130_2D(ii,2), BarL180_2D(ii,2), BarL230_2D(ii,2)], linspace(BarL0_2D(ii,1), BarL230_2D(ii, 1), 100));
    
    
    
    linePB = linePB - (BarL0_2D(ii,2) + BarL230_2D(ii,2))/2;
    [~, abs_Max_Index] = max(abs(linePB));
    rPB(ii,1) = linePB(abs_Max_Index);
%     rPB(ii, 1) = min(linePB);

    if false
%     if linePB(end) < -0.1 * 1e-3
        rPB(ii,1) = NaN;
    else
        plot(linspace(BarL0(ii,1), BarL230(ii, 1), 100), linePB)
    end
end
hold off

figure(2)
scatter(Data_Time, rPB)

rPB_Hirabayashi = mean(rPB, 'omitnan');


Bar_Data_All = (dlmread('Trimmed_BarL40kg1.trc','',6,1))/1000;%単位変換[mm]→[m]

Data_Time = Bar_Data_All(:,1) * 1000;

BarR0 = Bar_Data_All(:, 2:4);
BarR50 = Bar_Data_All(:, 5:7);
BarR100 = Bar_Data_All(:, 8:10);
BarR130 = Bar_Data_All(:, 11:13);
BarR180 = Bar_Data_All(:, 14:16);
BarR230 = Bar_Data_All(:, 17:19);

BarL0 = Bar_Data_All(:, 20:22);
BarL50 = Bar_Data_All(:, 23:25);
BarL100 = Bar_Data_All(:, 26:28);
BarL130 = Bar_Data_All(:, 29:31);
BarL180 = Bar_Data_All(:, 32:34);
BarL230 = Bar_Data_All(:, 35:37);

Axe_X_UnitVec = mean(BarR0 - BarR230, 1);
Axe_X_UnitVec(3) = 0;
Axe_X_UnitVec = Axe_X_UnitVec ./ norm(Axe_X_UnitVec);
Axe_Y_UnitVec = [0, 0, 1];

BarR0_2D = Get_2D_Coordinate(BarR0, Axe_X_UnitVec, Axe_Y_UnitVec);
BarR50_2D = Get_2D_Coordinate(BarR50, Axe_X_UnitVec, Axe_Y_UnitVec);
BarR100_2D = Get_2D_Coordinate(BarR100, Axe_X_UnitVec, Axe_Y_UnitVec);
BarR130_2D = Get_2D_Coordinate(BarR130, Axe_X_UnitVec, Axe_Y_UnitVec);
BarR180_2D = Get_2D_Coordinate(BarR180, Axe_X_UnitVec, Axe_Y_UnitVec);
BarR230_2D = Get_2D_Coordinate(BarR230, Axe_X_UnitVec, Axe_Y_UnitVec);

BarL0_2D = Get_2D_Coordinate(BarL0, Axe_X_UnitVec, Axe_Y_UnitVec);
BarL50_2D = Get_2D_Coordinate(BarL50, Axe_X_UnitVec, Axe_Y_UnitVec);
BarL100_2D = Get_2D_Coordinate(BarL100, Axe_X_UnitVec, Axe_Y_UnitVec);
BarL130_2D = Get_2D_Coordinate(BarL130, Axe_X_UnitVec, Axe_Y_UnitVec);
BarL180_2D = Get_2D_Coordinate(BarL180, Axe_X_UnitVec, Axe_Y_UnitVec);
BarL230_2D = Get_2D_Coordinate(BarL230, Axe_X_UnitVec, Axe_Y_UnitVec);

rPB = zeros(size(Data_Time));
linePBs = zeros(size(Data_Time, 1), 100);

figure(3)
clf('reset')
hold on
for ii = 1:size(Data_Time, 1)
    linePB = spline([BarL0_2D(ii,1), BarL50_2D(ii,1), BarL100_2D(ii,1), BarL130_2D(ii,1), BarL180_2D(ii,1), BarL230_2D(ii,1)],...
        [BarL0_2D(ii,2), BarL50_2D(ii,2), BarL100_2D(ii,2), BarL130_2D(ii,2), BarL180_2D(ii,2), BarL230_2D(ii,2)], linspace(BarL0_2D(ii,1), BarL230_2D(ii, 1), 100));
    
    
    
    linePB = linePB - (BarL0_2D(ii,2) + BarL230_2D(ii,2))/2;
    [~, abs_Max_Index] = max(abs(linePB));
    rPB(ii,1) = linePB(abs_Max_Index);
%     rPB(ii, 1) = min(linePB);

    if false
%     if linePB(end) < -0.1 * 1e-3
        rPB(ii,1) = NaN;
    else
        plot(linspace(BarL0(ii,1), BarL230(ii, 1), 100), linePB)
    end
end
hold off

figure(4)
scatter(Data_Time, rPB)

rPB_40kg = mean(rPB, 'omitnan');


Bar_Data_All = (dlmread('Trimmed_BarL20kg1.trc','',6,1))/1000;%単位変換[mm]→[m]

Data_Time = Bar_Data_All(:,1) * 1000;

BarR0 = Bar_Data_All(:, 2:4);
BarR50 = Bar_Data_All(:, 5:7);
BarR100 = Bar_Data_All(:, 8:10);
BarR130 = Bar_Data_All(:, 11:13);
BarR180 = Bar_Data_All(:, 14:16);
BarR230 = Bar_Data_All(:, 17:19);

BarL0 = Bar_Data_All(:, 20:22);
BarL50 = Bar_Data_All(:, 23:25);
BarL100 = Bar_Data_All(:, 26:28);
BarL130 = Bar_Data_All(:, 29:31);
BarL180 = Bar_Data_All(:, 32:34);
BarL230 = Bar_Data_All(:, 35:37);

Axe_X_UnitVec = mean(BarR0 - BarR230, 1);
Axe_X_UnitVec(3) = 0;
Axe_X_UnitVec = Axe_X_UnitVec ./ norm(Axe_X_UnitVec);
Axe_Y_UnitVec = [0, 0, 1];

BarR0_2D = Get_2D_Coordinate(BarR0, Axe_X_UnitVec, Axe_Y_UnitVec);
BarR50_2D = Get_2D_Coordinate(BarR50, Axe_X_UnitVec, Axe_Y_UnitVec);
BarR100_2D = Get_2D_Coordinate(BarR100, Axe_X_UnitVec, Axe_Y_UnitVec);
BarR130_2D = Get_2D_Coordinate(BarR130, Axe_X_UnitVec, Axe_Y_UnitVec);
BarR180_2D = Get_2D_Coordinate(BarR180, Axe_X_UnitVec, Axe_Y_UnitVec);
BarR230_2D = Get_2D_Coordinate(BarR230, Axe_X_UnitVec, Axe_Y_UnitVec);

BarL0_2D = Get_2D_Coordinate(BarL0, Axe_X_UnitVec, Axe_Y_UnitVec);
BarL50_2D = Get_2D_Coordinate(BarL50, Axe_X_UnitVec, Axe_Y_UnitVec);
BarL100_2D = Get_2D_Coordinate(BarL100, Axe_X_UnitVec, Axe_Y_UnitVec);
BarL130_2D = Get_2D_Coordinate(BarL130, Axe_X_UnitVec, Axe_Y_UnitVec);
BarL180_2D = Get_2D_Coordinate(BarL180, Axe_X_UnitVec, Axe_Y_UnitVec);
BarL230_2D = Get_2D_Coordinate(BarL230, Axe_X_UnitVec, Axe_Y_UnitVec);

rPB = zeros(size(Data_Time));
linePBs = zeros(size(Data_Time, 1), 100);

figure(5)
clf('reset')
hold on
for ii = 1:size(Data_Time, 1)
    linePB = spline([BarL0_2D(ii,1), BarL50_2D(ii,1), BarL100_2D(ii,1), BarL130_2D(ii,1), BarL180_2D(ii,1), BarL230_2D(ii,1)],...
        [BarL0_2D(ii,2), BarL50_2D(ii,2), BarL100_2D(ii,2), BarL130_2D(ii,2), BarL180_2D(ii,2), BarL230_2D(ii,2)], linspace(BarL0_2D(ii,1), BarL230_2D(ii, 1), 100));
    
    
    
    linePB = linePB - (BarL0_2D(ii,2) + BarL230_2D(ii,2))/2;
    [~, abs_Max_Index] = max(abs(linePB));
    rPB(ii,1) = linePB(abs_Max_Index);
%     rPB(ii, 1) = min(linePB);

    if false
%     if linePB(end) < -0.1 * 1e-3
        rPB(ii,1) = NaN;
    else
        plot(linspace(BarL0(ii,1), BarL230(ii, 1), 100), linePB)
    end
end
hold off

figure(6)
scatter(Data_Time, rPB)

rPB_20kg = mean(rPB, 'omitnan');

Bar_Data_All = (dlmread('Trimmed_BarR1Hirabayashi1.trc','',6,1))/1000;%単位変換[mm]→[m]

Data_Time = Bar_Data_All(:,1) * 1000;

BarR0 = Bar_Data_All(:, 2:4);
BarR50 = Bar_Data_All(:, 5:7);
BarR100 = Bar_Data_All(:, 8:10);
BarR130 = Bar_Data_All(:, 11:13);
BarR180 = Bar_Data_All(:, 14:16);
BarR230 = Bar_Data_All(:, 17:19);

BarL0 = Bar_Data_All(:, 20:22);
BarL50 = Bar_Data_All(:, 23:25);
BarL100 = Bar_Data_All(:, 26:28);
BarL130 = Bar_Data_All(:, 29:31);
BarL180 = Bar_Data_All(:, 32:34);
BarL230 = Bar_Data_All(:, 35:37);

Axe_X_UnitVec = mean(BarR0 - BarR230, 1);
Axe_X_UnitVec(3) = 0;
Axe_X_UnitVec = Axe_X_UnitVec ./ norm(Axe_X_UnitVec);
Axe_Y_UnitVec = [0, 0, 1];

BarR0_2D = Get_2D_Coordinate(BarR0, Axe_X_UnitVec, Axe_Y_UnitVec);
BarR50_2D = Get_2D_Coordinate(BarR50, Axe_X_UnitVec, Axe_Y_UnitVec);
BarR100_2D = Get_2D_Coordinate(BarR100, Axe_X_UnitVec, Axe_Y_UnitVec);
BarR130_2D = Get_2D_Coordinate(BarR130, Axe_X_UnitVec, Axe_Y_UnitVec);
BarR180_2D = Get_2D_Coordinate(BarR180, Axe_X_UnitVec, Axe_Y_UnitVec);
BarR230_2D = Get_2D_Coordinate(BarR230, Axe_X_UnitVec, Axe_Y_UnitVec);

BarL0_2D = Get_2D_Coordinate(BarL0, Axe_X_UnitVec, Axe_Y_UnitVec);
BarL50_2D = Get_2D_Coordinate(BarL50, Axe_X_UnitVec, Axe_Y_UnitVec);
BarL100_2D = Get_2D_Coordinate(BarL100, Axe_X_UnitVec, Axe_Y_UnitVec);
BarL130_2D = Get_2D_Coordinate(BarL130, Axe_X_UnitVec, Axe_Y_UnitVec);
BarL180_2D = Get_2D_Coordinate(BarL180, Axe_X_UnitVec, Axe_Y_UnitVec);
BarL230_2D = Get_2D_Coordinate(BarL230, Axe_X_UnitVec, Axe_Y_UnitVec);

rPB = zeros(size(Data_Time));
linePBs = zeros(size(Data_Time, 1), 100);

figure(7)
clf('reset')
hold on
for ii = 1:size(Data_Time, 1)
    linePB = spline([BarL0_2D(ii,1), BarL50_2D(ii,1), BarL100_2D(ii,1), BarL130_2D(ii,1), BarL180_2D(ii,1), BarL230_2D(ii,1)],...
        [BarL0_2D(ii,2), BarL50_2D(ii,2), BarL100_2D(ii,2), BarL130_2D(ii,2), BarL180_2D(ii,2), BarL230_2D(ii,2)], linspace(BarL0_2D(ii,1), BarL230_2D(ii, 1), 100));
    
    
    
    linePB = linePB - (BarL0_2D(ii,2) + BarL230_2D(ii,2))/2;
    [~, abs_Max_Index] = max(abs(linePB));
    rPB(ii,1) = linePB(abs_Max_Index);
%     rPB(ii, 1) = min(linePB);

%     if false
    if linePB(end) < -0.1 * 1e-3
        rPB(ii,1) = NaN;
    else
        plot(linspace(BarL0(ii,1), BarL230(ii, 1), 100), linePB)
    end
end
hold off

figure(8)
scatter(Data_Time, rPB)

rPB_0kg = mean(rPB, 'omitnan');

g = 9.8;
% X = -[rPB_Hirabayashi, rPB_40kg, rPB_20kg];
% Y = [81.3, 40, 20] * g;
X = -[rPB_Hirabayashi, rPB_40kg, rPB_20kg, rPB_0kg];
Y = [81.3, 40, 20, 0] * g;

figure(9)
scatter(X, Y)
rPB_polyfit = polyfit(X, Y, 1)
rPB_polyfit_Y = polyval(rPB_polyfit, X);
hold on
plot(X, rPB_polyfit_Y)
hold off

mPB = rPB_polyfit(2) / g















































