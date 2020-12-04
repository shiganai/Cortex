clear all

Bar_Data_All = (dlmread('Trimmed_Trimmed_ShSs3-Bar.trc','',6,1))/1000;%単位変換[mm]→[m]

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

Bar0 = (BarR0_2D + BarL0_2D)/2;
Bar50 = (BarR50_2D + BarL50_2D)/2;
Bar100 = (BarR100_2D + BarL100_2D)/2;
Bar130 = (BarR130_2D + BarL130_2D)/2;
Bar180 = (BarR180_2D + BarL180_2D)/2;
Bar230 = (BarR230_2D + BarL230_2D)/2;

BarRs(:,:,1) = [BarR0_2D(:,1), BarR50_2D(:,1), BarR100_2D(:,1), BarR130_2D(:,1), BarR180_2D(:,1), BarR230_2D(:,1)];
BarRs(:,:,2) = [BarR0_2D(:,2), BarR50_2D(:,2), BarR100_2D(:,2), BarR130_2D(:,2), BarR180_2D(:,2), BarR230_2D(:,2)];
BarLs(:,:,1) = [BarL0_2D(:,1), BarL50_2D(:,1), BarL100_2D(:,1), BarL130_2D(:,1), BarL180_2D(:,1), BarL230_2D(:,1)];
BarLs(:,:,2) = [BarL0_2D(:,2), BarL50_2D(:,2), BarL100_2D(:,2), BarL130_2D(:,2), BarL180_2D(:,2), BarL230_2D(:,2)];


Body_Data_All = (dlmread('Trimmed_Trimmed_ShSs3-51P.trc','',6,2))/1000;%単位変換[mm]→[m]

FileName = 'Trimmed_Trimmed_ShSs3-51P.trc';
[handRI, handRO, wrRI, wrRO, elbRI, elbRO, shRF, shRU, shRB, handLI, handLO, wrLI, wrLO, elbLI, elbLO, shLF, shLU, shLB, toeR, baRI, baRO, heelR, anRI, anRO, knRI, knRO, troR, toeL, baLI, baLO, heelL, anLI, anLO, knLI, knLO, troL, head, earR, earL, clav, c7, ribR, ribL, xiph, t12, ASISR, ASISL, PSISR, PSISL, ThR, ThL]...
    = Get_2D_Coordinate_51P(FileName, Axe_X_UnitVec, Axe_Y_UnitVec);

handR = (handRI + handRO)/2;
wrR = (wrRI + wrRO)/2;
elbR = (elbRI + elbRO)/2;
shR = (shRF + shRU + shRB)/3;
% ribR = ribR;
% troR = troR;
knR = (knRI + knRO)/2;
anR = (anRI + anRO)/2;
% toeR = toeR

BodyR(:,:,1) = [handR(:,1), wrR(:,1), elbR(:,1), shR(:,1), ribR(:,1), troR(:,1), knR(:,1), anR(:,1), toeR(:,1)];
BodyR(:,:,2) = [handR(:,2), wrR(:,2), elbR(:,2), shR(:,2), ribR(:,2), troR(:,2), knR(:,2), anR(:,2), toeR(:,2)];

handL = (handLI + handLO)/2;
wrL = (wrLI + wrLO)/2;
elbL = (elbLI + elbLO)/2;
shL = (shLF + shLU + shLB)/3;
knL = (knLI + knLO)/2;
anL = (anLI + anLO)/2;

BodyL(:,:,1) = [handL(:,1), wrL(:,1), elbL(:,1), shL(:,1), ribL(:,1), troL(:,1), knL(:,1), anL(:,1), toeL(:,1)];
BodyL(:,:,2) = [handL(:,2), wrL(:,2), elbL(:,2), shL(:,2), ribL(:,2), troL(:,2), knL(:,2), anL(:,2), toeL(:,2)];

top = head;
ear = (earR + earL)/2;
hand = (handR + handL)/2;
wr = (wrR + wrL)/2;
elb = (elbR + elbL)/2;
sh = (shR + shL)/2;
rib = (t12 + ribR + ribL)/3;
% rib = (ribR + ribL)/2;
tro = (troR + troL)/2;
kn = (knR + knL)/2;
an = (anR + anL)/2;
toe = (toeR + toeL)/2;

mAll = 81;
r_Toe_Ankle = mean(vecnorm(an - toe, 2, 2));
r_Ankle_Knee = mean(vecnorm(kn - an, 2, 2));
r_Knee_Waist = mean(vecnorm(tro - kn, 2, 2));
r_Waist_Rib = mean(vecnorm(rib - tro, 2, 2));
r_Rib_Shoulder = mean(vecnorm(sh - rib, 2, 2));
r_Shoulder_Elbow = mean(vecnorm(elb - sh, 2, 2));
r_Elbow_Wrist = mean(vecnorm(wr - elb, 2, 2));
r_Wrist_Finger = mean(vecnorm(hand - wr, 2, 2));
r_Shoulder_Ear = mean(vecnorm(ear - sh, 2, 2));
r_Ear_Top = mean(vecnorm(top - ear, 2, 2));


% [mArm, mBody, mLeg, rArm, rBody, rLeg, rArmMCD, rBodyMCD, rLegMCD, InertiaArm, InertiaBody, InertiaLeg] = ...
%     Calc_Parameter_AE(mAll, r_Toe_Ankle, r_Ankle_Knee, r_Knee_Waist, r_Rib_Shoulder, r_Shoulder_Elbow, r_Elbow_Wrist, r_Wrist_Finger, r_Shoulder_Ear, r_Ear_Top);
[mArm, mUBody, mLBody, mLeg, rArm, rUBody, rLBody, rLeg, rArmMCD, rUBodyMCD, rLBodyMCD, rLegMCD, InertiaArm, InertiaUBody, InertiaLBody, InertiaLeg]...
    = Calc_Parameter_AE_4Links(mAll, r_Toe_Ankle, r_Ankle_Knee, r_Knee_Waist, r_Waist_Rib, r_Rib_Shoulder, r_Shoulder_Elbow, r_Elbow_Wrist, r_Wrist_Finger, r_Shoulder_Ear, r_Ear_Top);

% wr_sh = sh - wr;
% sh_tro = tro - sh;
% tro_an = an - tro;

% rPB = hand(:,2) - hand(1,2);
hand_x = hand(:,1);
hand_y = hand(:,2);
thHand = GetThetaFromXY(wr, sh) - 1/2*pi;
thShoulder = GetThetaFromXY(sh, rib) - thHand - 1/2*pi;
thRib = GetThetaFromXY(rib, tro) - thHand - 1/2*pi - thShoulder;
thWaist = GetThetaFromXY(tro, an) - thHand - 1/2*pi - thShoulder - thRib;

thShoulder = thShoulder + 2*pi;
thRib = thRib - 0*pi;
thWaist = thWaist - 2*pi;

rPB = zeros(size(Data_Time));
linePBs = zeros(size(Data_Time, 1), 100);
for ii = 1:size(Data_Time, 1)
    linePB = spline([Bar0(ii,1), Bar50(ii,1), Bar100(ii,1), Bar130(ii,1), Bar180(ii,1), Bar230(ii,1)],...
        [Bar0(ii,2), Bar50(ii,2), Bar100(ii,2), Bar130(ii,2), Bar180(ii,2), Bar230(ii,2)], linspace(Bar0(ii,1), Bar230(ii, 1), 100));
    
    linePB = linePB - (Bar0(ii,2) + Bar230(ii,2))/2;
    [~, abs_Max_Index] = max(abs(linePB));
    rPB(ii,1) = linePB(abs_Max_Index);
%     rPB(ii, 1) = min(linePB);
    
    %{
    if ii ~= 1
        hold on
    else
        figure(1)
    end
    plot(linspace(Bar0(ii,1), Bar230(ii, 1), 100), linePB)
    hold off
    %}
end

dthHand = diff(thHand) / diff(Data_Time(1:2));
ddthHand = diff(dthHand) / diff(Data_Time(1:2));

dthShoulder = diff(thShoulder) / diff(Data_Time(1:2));
ddthShoulder = diff(dthShoulder) / diff(Data_Time(1:2));

dthRib = diff(thRib) / diff(Data_Time(1:2));
ddthRib = diff(dthRib) / diff(Data_Time(1:2));

dthWaist = diff(thWaist) / diff(Data_Time(1:2));
ddthWaist = diff(dthWaist) / diff(Data_Time(1:2));

dhand_x = diff(hand_x) / diff(Data_Time(1:2));
ddhand_x = diff(dhand_x) / diff(Data_Time(1:2));

dhand_y = diff(hand_y) / diff(Data_Time(1:2));
ddhand_y = diff(dhand_y) / diff(Data_Time(1:2));

drPB = diff(rPB) / diff(Data_Time(1:2));
ddrPB = diff(drPB) / diff(Data_Time(1:2));

Cut_Off_Freq_tmp = (0.5:0.5:30)';
Nichest_Freq = 200 / 2;
Cut_Range = find(Data_Time == 3.92):find(Data_Time == 7.95);

Data_Time = Data_Time - Data_Time(Cut_Range(1));

%{/
RMS_ddthHand = zeros(size(Cut_Off_Freq_tmp));

for ii = 1:size(Cut_Off_Freq_tmp, 1)
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
    ddthHand_filt_tmp = filtfilt(Butter_b, Butter_a, ddthHand);
    
    RMS_ddthHand(ii) = (mean((ddthHand_filt_tmp - ddthHand).^2))^(1/2);
end

figure(1)
plot(Cut_Off_Freq_tmp, RMS_ddthHand, 'o-')
title('RMS\_thHand')

Cut_Off_Freq_Line_Range = find(abs(Cut_Off_Freq_tmp - 5) < 1e-2):find(abs(Cut_Off_Freq_tmp - 6.5) < 1e-2);
Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddthHand(Cut_Off_Freq_Line_Range), 1);
Cut_Off_Freq_Line_y = polyval(Cut_Off_Freq_Line, Cut_Off_Freq_tmp);
hold on
plot(Cut_Off_Freq_tmp, Cut_Off_Freq_Line_y)
hold off

Cut_Off_Freq = interp1(RMS_ddthHand(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));

[Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
ddthHand_filt = filtfilt(Butter_b, Butter_a, ddthHand);


figure(2)
plot(Data_Time(Cut_Range), ddthHand(Cut_Range), 'Displayname', 'raw', 'LineStyle', ':', 'LineWidth', 1)
hold on
plot(Data_Time(Cut_Range), ddthHand_filt(Cut_Range), 'Displayname', 'filtfilt')
hold off
legend
title('ddthHand')

ddthHand_spline = spline(Data_Time(Cut_Range), ddthHand_filt(Cut_Range));
% dthHand_spline = fnint(ddthHand_spline, 0);
dthHand_spline = fnint(ddthHand_spline, mean(dthHand(Cut_Range(1:2))) * 0.3);
thHand_spline = fnint(dthHand_spline, mean(thHand(Cut_Range(1:2))));

figure(3)
plot(Data_Time(Cut_Range), dthHand(Cut_Range), 'Displayname', 'raw', 'LineStyle', ':', 'LineWidth', 3)
hold on
plot(Data_Time(Cut_Range), fnval(dthHand_spline, Data_Time(Cut_Range)), 'Displayname', 'filtfilt')
hold off
legend
title('dthHand')

figure(4)
plot(Data_Time(Cut_Range), thHand(Cut_Range), 'Displayname', 'raw', 'LineStyle', ':', 'LineWidth', 3)
hold on
plot(Data_Time(Cut_Range), fnval(thHand_spline, Data_Time(Cut_Range)), 'Displayname', 'filtfilt')
hold off
legend
title('thHand')


RMS_ddthShoulder = zeros(size(Cut_Off_Freq_tmp));

for ii = 1:size(Cut_Off_Freq_tmp, 1)
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
    ddthShoulder_filt_tmp = filtfilt(Butter_b, Butter_a, ddthShoulder);
    
    RMS_ddthShoulder(ii) = (mean((ddthShoulder_filt_tmp - ddthShoulder).^2))^(1/2);
end

figure(5)
plot(Cut_Off_Freq_tmp, RMS_ddthShoulder, 'o-')
title('RMS\_thShoulder')

Cut_Off_Freq_Line_Range = find(abs(Cut_Off_Freq_tmp - 4) < 1e-2):find(abs(Cut_Off_Freq_tmp - 6) < 1e-2);
% Cut_Off_Freq_Line_Range = find(Cut_Off_Freq_tmp == 3.5):find(Cut_Off_Freq_tmp == 9);
Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddthShoulder(Cut_Off_Freq_Line_Range), 1);
Cut_Off_Freq_Line_y = polyval(Cut_Off_Freq_Line, Cut_Off_Freq_tmp);
hold on
plot(Cut_Off_Freq_tmp, Cut_Off_Freq_Line_y)
hold off

Cut_Off_Freq = interp1(RMS_ddthShoulder(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));

[Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
ddthShoulder_filt = filtfilt(Butter_b, Butter_a, ddthShoulder);

% Cut_Range = find(Data_Time == 3.5):find(Data_Time == 8.5);

figure(6)
plot(Data_Time(Cut_Range), ddthShoulder(Cut_Range), 'Displayname', 'raw', 'LineStyle', ':', 'LineWidth', 1)
hold on
plot(Data_Time(Cut_Range), ddthShoulder_filt(Cut_Range), 'Displayname', 'filtfilt')
hold off
legend
title('ddthShoulder')

ddthShoulder_spline = spline(Data_Time(Cut_Range), ddthShoulder_filt(Cut_Range));
% dthShoulder_spline = fnint(ddthShoulder_spline, mean(dthShoulder(Cut_Range(1:2))));
dthShoulder_spline = fnint(ddthShoulder_spline, mean(dthShoulder(Cut_Range(1:2))) * 0.4);
thShoulder_spline = fnint(dthShoulder_spline, mean(thShoulder(Cut_Range(1:2))));

figure(7)
plot(Data_Time(Cut_Range), dthShoulder(Cut_Range), 'Displayname', 'raw', 'LineStyle', ':', 'LineWidth', 3)
hold on
plot(Data_Time(Cut_Range), fnval(dthShoulder_spline, Data_Time(Cut_Range)), 'Displayname', 'filtfilt')
hold off
legend
title('dthShoulder')

figure(8)
plot(Data_Time(Cut_Range), thShoulder(Cut_Range), 'Displayname', 'raw', 'LineStyle', ':', 'LineWidth', 3)
hold on
plot(Data_Time(Cut_Range), fnval(thShoulder_spline, Data_Time(Cut_Range)), 'Displayname', 'filtfilt')
hold off
legend
title('thShoulder')

RMS_ddthWaist = zeros(size(Cut_Off_Freq_tmp));

for ii = 1:size(Cut_Off_Freq_tmp, 1)
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
    ddthWaist_filt_tmp = filtfilt(Butter_b, Butter_a, ddthWaist);
    
    RMS_ddthWaist(ii) = (mean((ddthWaist_filt_tmp - ddthWaist).^2))^(1/2);
end

figure(9)
plot(Cut_Off_Freq_tmp, RMS_ddthWaist, 'o-')
title('RMS\_thWaist')

Cut_Off_Freq_Line_Range = find(abs(Cut_Off_Freq_tmp - 5) < 1e-2):find(abs(Cut_Off_Freq_tmp - 6) < 1e-2);
% Cut_Off_Freq_Line_Range = find(Cut_Off_Freq_tmp == 3.5):find(Cut_Off_Freq_tmp == 9);
Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddthWaist(Cut_Off_Freq_Line_Range), 1);
Cut_Off_Freq_Line_y = polyval(Cut_Off_Freq_Line, Cut_Off_Freq_tmp);
hold on
plot(Cut_Off_Freq_tmp, Cut_Off_Freq_Line_y)
hold off

Cut_Off_Freq = interp1(RMS_ddthWaist(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));

[Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
ddthWaist_filt = filtfilt(Butter_b, Butter_a, ddthWaist);

figure(10)
plot(Data_Time(Cut_Range), ddthWaist(Cut_Range), 'Displayname', 'raw', 'LineStyle', ':', 'LineWidth', 1)
hold on
plot(Data_Time(Cut_Range), ddthWaist_filt(Cut_Range), 'Displayname', 'filtfilt')
hold off
legend
title('ddthWaist')

ddthWaist_spline = spline(Data_Time(Cut_Range), ddthWaist_filt(Cut_Range));
% dthWaist_spline = fnint(ddthWaist_spline, mean(dthWaist(Cut_Range(1:2))));
dthWaist_spline = fnint(ddthWaist_spline, mean(dthWaist(Cut_Range(1:2))) * 0.9);
thWaist_spline = fnint(dthWaist_spline, mean(thWaist(Cut_Range(1:2))));

figure(11)
plot(Data_Time(Cut_Range), dthWaist(Cut_Range), 'Displayname', 'raw', 'LineStyle', ':', 'LineWidth', 3)
hold on
plot(Data_Time(Cut_Range), fnval(dthWaist_spline, Data_Time(Cut_Range)), 'Displayname', 'filtfilt')
hold off
legend
title('dthWaist')

figure(12)
plot(Data_Time(Cut_Range), thWaist(Cut_Range), 'Displayname', 'raw', 'LineStyle', ':', 'LineWidth', 3)
hold on
plot(Data_Time(Cut_Range), fnval(thWaist_spline, Data_Time(Cut_Range)), 'Displayname', 'filtfilt')
hold off
legend
title('thWaist')

RMS_ddhand_x = zeros(size(Cut_Off_Freq_tmp));

for ii = 1:size(Cut_Off_Freq_tmp, 1)
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
    ddhand_x_filt_tmp = filtfilt(Butter_b, Butter_a, ddhand_x);
    
    RMS_ddhand_x(ii) = (mean((ddhand_x_filt_tmp - ddhand_x).^2))^(1/2);
end

figure(13)
plot(Cut_Off_Freq_tmp, RMS_ddhand_x, 'o-')
title('RMS\_hand\_x')

Cut_Off_Freq_Line_Range = find(abs(Cut_Off_Freq_tmp - 5) < 1e-2):find(abs(Cut_Off_Freq_tmp - 10) < 1e-2);
% Cut_Off_Freq_Line_Range = find(Cut_Off_Freq_tmp == 10):find(Cut_Off_Freq_tmp == 30);
Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddhand_x(Cut_Off_Freq_Line_Range), 1);
Cut_Off_Freq_Line_y = polyval(Cut_Off_Freq_Line, Cut_Off_Freq_tmp);
hold on
plot(Cut_Off_Freq_tmp, Cut_Off_Freq_Line_y)
hold off

Cut_Off_Freq = interp1(RMS_ddhand_x(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));

[Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
ddhand_x_filt = filtfilt(Butter_b, Butter_a, ddhand_x);

figure(14)
plot(Data_Time(Cut_Range), ddhand_x(Cut_Range), 'Displayname', 'raw', 'LineStyle', ':', 'LineWidth', 1)
hold on
plot(Data_Time(Cut_Range), ddhand_x_filt(Cut_Range), 'Displayname', 'filtfilt')
hold off
legend
title('ddhand\_x')

ddhand_x_spline = spline(Data_Time(Cut_Range), ddhand_x_filt(Cut_Range));
% dhand_x_spline = fnint(ddhand_x_spline, mean(dhand_x(Cut_Range(1:2))));
dhand_x_spline = fnint(ddhand_x_spline, mean(dhand_x(Cut_Range(1:2))) * 0.5);
% dhand_x_spline = fnint(ddhand_x_spline, dhand_x(Cut_Range(1)) * 1.15);
% dhand_x_spline = fnint(ddhand_x_spline, 0);
hand_x_spline = fnint(dhand_x_spline, mean(hand_x(Cut_Range(1:2))));

figure(15)
plot(Data_Time(Cut_Range), dhand_x(Cut_Range), 'Displayname', 'raw', 'LineStyle', ':', 'LineWidth', 3)
hold on
plot(Data_Time(Cut_Range), fnval(dhand_x_spline, Data_Time(Cut_Range)), 'Displayname', 'filtfilt')
hold off
legend
title('dhand\_x')

figure(16)
plot(Data_Time(Cut_Range), hand_x(Cut_Range), 'Displayname', 'raw', 'LineStyle', ':', 'LineWidth', 3)
hold on
plot(Data_Time(Cut_Range), fnval(hand_x_spline, Data_Time(Cut_Range)), 'Displayname', 'filtfilt')
hold off
legend
title('hand\_x')

RMS_ddhand_y = zeros(size(Cut_Off_Freq_tmp));

for ii = 1:size(Cut_Off_Freq_tmp, 1)
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
    ddhand_y_filt_tmp = filtfilt(Butter_b, Butter_a, ddhand_y);
    
    RMS_ddhand_y(ii) = (mean((ddhand_y_filt_tmp - ddhand_y).^2))^(1/2);
end

figure(17)
plot(Cut_Off_Freq_tmp, RMS_ddhand_y, 'o-')
title('RMS\_hand\_y')

Cut_Off_Freq_Line_Range = find(abs(Cut_Off_Freq_tmp - 10) < 1e-2):find(abs(Cut_Off_Freq_tmp - 15) < 1e-2);
% Cut_Off_Freq_Line_Range = find(Cut_Off_Freq_tmp == 3.5):find(Cut_Off_Freq_tmp == 8.5);
% Cut_Off_Freq_Line_Range = find(Cut_Off_Freq_tmp == 10):find(Cut_Off_Freq_tmp == 30);
Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddhand_y(Cut_Off_Freq_Line_Range), 1);
Cut_Off_Freq_Line_y = polyval(Cut_Off_Freq_Line, Cut_Off_Freq_tmp);
hold on
plot(Cut_Off_Freq_tmp, Cut_Off_Freq_Line_y)
hold off

Cut_Off_Freq = interp1(RMS_ddhand_y(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));

[Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
ddhand_y_filt = filtfilt(Butter_b, Butter_a, ddhand_y);

% Cut_Range = find(Data_Time == 3.5):find(Data_Time == 8.5);

figure(18)
plot(Data_Time(Cut_Range), ddhand_y(Cut_Range), 'Displayname', 'raw', 'LineStyle', ':', 'LineWidth', 1)
hold on
plot(Data_Time(Cut_Range), ddhand_y_filt(Cut_Range), 'Displayname', 'filtfilt')
hold off
legend
title('ddhand\_y')

ddhand_y_spline = spline(Data_Time(Cut_Range), ddhand_y_filt(Cut_Range));
% dhand_y_spline = fnint(ddhand_y_spline, dhand_y(Cut_Range(1)) * 0.90);
% dhand_y_spline = fnint(ddhand_y_spline, mean(dhand_y(Cut_Range(1:2))));
dhand_y_spline = fnint(ddhand_y_spline, mean(dhand_y(Cut_Range(1:2))) * 1.2);
hand_y_spline = fnint(dhand_y_spline, mean(hand_y(Cut_Range(1:2))));

figure(19)
plot(Data_Time(Cut_Range), dhand_y(Cut_Range), 'Displayname', 'raw', 'LineStyle', ':', 'LineWidth', 3)
hold on
plot(Data_Time(Cut_Range), fnval(dhand_y_spline, Data_Time(Cut_Range)), 'Displayname', 'filtfilt')
hold off
legend
title('dhand\_y')

figure(20)
plot(Data_Time(Cut_Range), hand_y(Cut_Range), 'Displayname', 'raw', 'LineStyle', ':', 'LineWidth', 3)
hold on
plot(Data_Time(Cut_Range), fnval(hand_y_spline, Data_Time(Cut_Range)), 'Displayname', 'filtfilt')
hold off
legend
title('hand\_y')


RMS_ddrPB = zeros(size(Cut_Off_Freq_tmp));

for ii = 1:size(Cut_Off_Freq_tmp, 1)
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
    ddrPB_filt_tmp = filtfilt(Butter_b, Butter_a, ddrPB);
    
    RMS_ddrPB(ii) = (mean((ddrPB_filt_tmp - ddrPB).^2))^(1/2);
end

figure(21)
plot(Cut_Off_Freq_tmp, RMS_ddrPB, 'o-')
title('RMS\_rPB')

Cut_Off_Freq_Line_Range = find(abs(Cut_Off_Freq_tmp - 25) < 1e-2):find(abs(Cut_Off_Freq_tmp - 30) < 1e-2);
% Cut_Off_Freq_Line_Range = find(Cut_Off_Freq_tmp == 3.5):find(Cut_Off_Freq_tmp == 8.5);
% Cut_Off_Freq_Line_Range = find(Cut_Off_Freq_tmp == 25):find(Cut_Off_Freq_tmp == 30);
Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddrPB(Cut_Off_Freq_Line_Range), 1);
Cut_Off_Freq_Line_y = polyval(Cut_Off_Freq_Line, Cut_Off_Freq_tmp);
hold on
plot(Cut_Off_Freq_tmp, Cut_Off_Freq_Line_y)
hold off

% Cut_Off_Freq = interp1(flip(RMS_ddrPB(1:Cut_Off_Freq_Line_Range(1))), flip(Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1))), Cut_Off_Freq_Line(2))
Cut_Off_Freq = interp1(RMS_ddrPB(find(Cut_Off_Freq_tmp == 7.5):Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(find(Cut_Off_Freq_tmp == 7.5):Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));

[Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
ddrPB_filt = filtfilt(Butter_b, Butter_a, ddrPB);

figure(22)
plot(Data_Time(Cut_Range), ddrPB(Cut_Range), 'Displayname', 'raw', 'LineStyle', ':', 'LineWidth', 1)
hold on
plot(Data_Time(Cut_Range), ddrPB_filt(Cut_Range), 'Displayname', 'filtfilt')
hold off
legend
title('ddrPB')

ddrPB_spline = spline(Data_Time(Cut_Range), ddrPB_filt(Cut_Range));
% drPB_spline = fnint(ddrPB_spline, drPB(Cut_Range(1)) * 0.95);
% drPB_spline = fnint(ddrPB_spline, mean(drPB(Cut_Range(1:2))));
drPB_spline = fnint(ddrPB_spline, mean(drPB(Cut_Range(1:2))) * 0.4);
rPB_spline = fnint(drPB_spline, mean(rPB(Cut_Range(1:2))));

figure(23)
plot(Data_Time(Cut_Range), drPB(Cut_Range), 'Displayname', 'raw', 'LineStyle', ':', 'LineWidth', 3)
hold on
plot(Data_Time(Cut_Range), fnval(drPB_spline, Data_Time(Cut_Range)), 'Displayname', 'filtfilt')
hold off
legend
title('drPB')

figure(24)
plot(Data_Time(Cut_Range), rPB(Cut_Range), 'Displayname', 'raw', 'LineStyle', ':', 'LineWidth', 3)
hold on
plot(Data_Time(Cut_Range), fnval(rPB_spline, Data_Time(Cut_Range)), 'Displayname', 'filtfilt')
hold off
legend
title('rPB')


RMS_ddthRib = zeros(size(Cut_Off_Freq_tmp));

for ii = 1:size(Cut_Off_Freq_tmp, 1)
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
    ddthRib_filt_tmp = filtfilt(Butter_b, Butter_a, ddthRib);
    
    RMS_ddthRib(ii) = (mean((ddthRib_filt_tmp - ddthRib).^2))^(1/2);
end

figure(25)
plot(Cut_Off_Freq_tmp, RMS_ddthRib, 'o-')
title('RMS\_thRib')

Cut_Off_Freq_Line_Range = find(abs(Cut_Off_Freq_tmp - 4) < 1e-2):find(abs(Cut_Off_Freq_tmp - 5) < 1e-2);
% Cut_Off_Freq_Line_Range = find(Cut_Off_Freq_tmp == 5):find(Cut_Off_Freq_tmp == 6);
Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddthRib(Cut_Off_Freq_Line_Range), 1);
Cut_Off_Freq_Line_y = polyval(Cut_Off_Freq_Line, Cut_Off_Freq_tmp);
hold on
plot(Cut_Off_Freq_tmp, Cut_Off_Freq_Line_y)
hold off

% Cut_Off_Freq = interp1(flip(RMS_ddthRib(1:Cut_Off_Freq_Line_Range(1))), flip(Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1))), Cut_Off_Freq_Line(2))
Cut_Off_Freq = interp1(RMS_ddthRib(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));

[Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
ddthRib_filt = filtfilt(Butter_b, Butter_a, ddthRib);

figure(26)
plot(Data_Time(Cut_Range), ddthRib(Cut_Range), 'Displayname', 'raw', 'LineStyle', ':', 'LineWidth', 1)
hold on
plot(Data_Time(Cut_Range), ddthRib_filt(Cut_Range), 'Displayname', 'filtfilt')
hold off
legend
title('ddthRib')

ddthRib_spline = spline(Data_Time(Cut_Range), ddthRib_filt(Cut_Range));
% dthRib_spline = fnint(ddthRib_spline, dthRib(Cut_Range(1)) * 0.95);
% dthRib_spline = fnint(ddthRib_spline, mean(dthRib(Cut_Range(1:3))) * 0.1);
dthRib_spline = fnint(ddthRib_spline, mean(dthRib(Cut_Range(1:2))) * -2);
thRib_spline = fnint(dthRib_spline, mean(thRib(Cut_Range(1:2))));

figure(27)
plot(Data_Time(Cut_Range), dthRib(Cut_Range), 'Displayname', 'raw', 'LineStyle', ':', 'LineWidth', 3)
hold on
plot(Data_Time(Cut_Range), fnval(dthRib_spline, Data_Time(Cut_Range)), 'Displayname', 'filtfilt')
hold off
legend
title('dthRib')

figure(28)
plot(Data_Time(Cut_Range), thRib(Cut_Range), 'Displayname', 'raw', 'LineStyle', ':', 'LineWidth', 3)
hold on
plot(Data_Time(Cut_Range), fnval(thRib_spline, Data_Time(Cut_Range)), 'Displayname', 'filtfilt')
hold off
legend
title('thRib')

%}

%{/
%{
ddthHand_filt = fnval(ddthHand_spline, Data_Time(Cut_Range));
dthHand_filt = fnval(dthHand_spline, Data_Time(Cut_Range));
thHand_filt = fnval(thHand_spline, Data_Time(Cut_Range));
ddthShoulder_filt = fnval(ddthShoulder_spline, Data_Time(Cut_Range));
dthShoulder_filt = fnval(dthShoulder_spline, Data_Time(Cut_Range));
thShoulder_filt = fnval(thShoulder_spline, Data_Time(Cut_Range));
ddthRib_filt = fnval(ddthRib_spline, Data_Time(Cut_Range));
dthRib_filt = fnval(dthRib_spline, Data_Time(Cut_Range));
thRib_filt = fnval(thRib_spline, Data_Time(Cut_Range));
ddthWaist_filt = fnval(ddthWaist_spline, Data_Time(Cut_Range));
dthWaist_filt = fnval(dthWaist_spline, Data_Time(Cut_Range));
thWaist_filt = fnval(thWaist_spline, Data_Time(Cut_Range));
ddhand_x_filt = fnval(ddhand_x_spline, Data_Time(Cut_Range));
dhand_x_filt = fnval(dhand_x_spline, Data_Time(Cut_Range));
hand_x_filt = fnval(hand_x_spline, Data_Time(Cut_Range));
ddhand_y_filt = fnval(ddhand_y_spline, Data_Time(Cut_Range));
dhand_y_filt = fnval(dhand_y_spline, Data_Time(Cut_Range));
hand_y_filt = fnval(hand_y_spline, Data_Time(Cut_Range));

rPB_filt = fnval(rPB_spline, Data_Time(Cut_Range));
drPB_filt = fnval(drPB_spline, Data_Time(Cut_Range));
ddrPB_filt = fnval(ddrPB_spline, Data_Time(Cut_Range));
%}

Data_Time = Data_Time(Cut_Range);

ddthHand = fnval(ddthHand_spline, Data_Time);
dthHand = fnval(dthHand_spline, Data_Time);
thHand = fnval(thHand_spline, Data_Time);
ddthShoulder = fnval(ddthShoulder_spline, Data_Time);
dthShoulder = fnval(dthShoulder_spline, Data_Time);
thShoulder = fnval(thShoulder_spline, Data_Time);
ddthRib = fnval(ddthRib_spline, Data_Time);
dthRib = fnval(dthRib_spline, Data_Time);
thRib = fnval(thRib_spline, Data_Time);
ddthWaist = fnval(ddthWaist_spline, Data_Time);
dthWaist = fnval(dthWaist_spline, Data_Time);
thWaist = fnval(thWaist_spline, Data_Time);
ddhand_x = fnval(ddhand_x_spline, Data_Time);
dhand_x = fnval(dhand_x_spline, Data_Time);
hand_x = fnval(hand_x_spline, Data_Time);
ddhand_y = fnval(ddhand_y_spline, Data_Time);
dhand_y = fnval(dhand_y_spline, Data_Time);
hand_y = fnval(hand_y_spline, Data_Time);

rPB = fnval(rPB_spline, Data_Time);
drPB = fnval(drPB_spline, Data_Time);
ddrPB = fnval(ddrPB_spline, Data_Time);

% rPB_filt = hand_y_filt;
% drPB_filt = dhand_y_filt;
% ddrPB_filt = ddhand_y_filt;

xHand = hand_x;
yHand = hand_y;
dxHand = dhand_x;
dyHand = dhand_y;
ddxHand = ddhand_x;
ddyHand = ddhand_y;

pHand = [hand_x, hand_y];

pShoulder = pHand + rArm * [cos(thHand+1/2*pi), sin(thHand+1/2*pi)];
pRib = pShoulder + rUBody * [cos(thHand+1/2*pi + thShoulder), sin(thHand+1/2*pi + thShoulder)];
pWaist = pRib + rLBody * [cos(thHand+1/2*pi + thShoulder + thRib), sin(thHand+1/2*pi + thShoulder + thRib)];
pToe = pWaist + rLeg * [cos(thHand+1/2*pi + thShoulder + thRib + thWaist), sin(thHand+1/2*pi + thShoulder + thRib + thWaist)];


Body(:,:,1) = [pHand(:,1), pShoulder(:,1), pRib(:,1), pWaist(:,1), pToe(:,1)];
Body(:,:,2) = [pHand(:,2), pShoulder(:,2), pRib(:,2), pWaist(:,2), pToe(:,2)];

kPB = 1.9831 * 1e4;
cPB = 4.1;
mPB = 2;
g = 9.8;

FrPB = find_FrPB(ddrPB,ddthRib,ddthHand,ddthWaist,ddthShoulder,dthRib,dthHand,dthWaist,dthShoulder,g,mArm,mLBody,mLeg,mPB,mUBody,rArm,rArmMCD,rLBody,rLBodyMCD,rLegMCD,rUBody,rUBodyMCD,thHand,thShoulder,thWaist,thRib);
TauHand = find_TauHand(InertiaLeg,InertiaArm,InertiaLBody,InertiaUBody,ddrPB,ddthRib,ddthHand,ddthWaist,ddthShoulder,dthRib,dthHand,dthWaist,dthShoulder,g,mArm,mLBody,mLeg,mUBody,rArm,rArmMCD,rLBody,rLBodyMCD,rLegMCD,rUBody,rUBodyMCD,thHand,thShoulder,thWaist,thRib);
TauShoulder = find_TauShoulder(InertiaLeg,InertiaLBody,InertiaUBody,ddrPB,ddthRib,ddthHand,ddthWaist,ddthShoulder,dthRib,dthHand,dthWaist,dthShoulder,g,mLBody,mLeg,mUBody,rArm,rLBody,rLBodyMCD,rLegMCD,rUBody,rUBodyMCD,thHand,thShoulder,thWaist,thRib);
TauRib = find_TauRib(InertiaLeg,InertiaLBody,ddrPB,ddthRib,ddthHand,ddthWaist,ddthShoulder,dthHand,dthWaist,dthShoulder,g,mLBody,mLeg,rArm,rLBody,rLBodyMCD,rLegMCD,rUBody,thHand,thShoulder,thWaist,thRib);
TauWaist = find_TauWaist(InertiaLeg,ddrPB,ddthRib,ddthHand,ddthWaist,ddthShoulder,dthRib,dthHand,dthShoulder,g,mLeg,rArm,rLBody,rLegMCD,rUBody,thHand,thShoulder,thWaist,thRib);

pG = find_pG(mArm,mLBody,mLeg,mUBody,rArm,rArmMCD,rLBody,rLBodyMCD,rLegMCD,rUBody,rUBodyMCD,thHand,thShoulder,thWaist,thRib,xHand,yHand);

momentumG = find_momentumG(InertiaLeg,InertiaArm,InertiaLBody,InertiaUBody,dthRib,dthHand,dthWaist,dthShoulder,dxHand,dyHand,mArm,mLBody,mLeg,mUBody,rArm,rArmMCD,rLBody,rLBodyMCD,rLegMCD,rUBody,rUBodyMCD,thHand,thShoulder,thWaist,thRib,xHand,yHand);
dmomentumG = find_dmomentumG(InertiaLeg,InertiaArm,InertiaLBody,InertiaUBody,ddthRib,ddthHand,ddthWaist,ddthShoulder,ddxHand,ddyHand,dthRib,dthHand,dthWaist,dthShoulder,mArm,mLBody,mLeg,mUBody,rArm,rArmMCD,rLBody,rLBodyMCD,rLegMCD,rUBody,rUBodyMCD,thHand,thShoulder,thWaist,thRib,xHand,yHand);

pG_Anime(:,:,1) = pG(:,1);
pG_Anime(:,:,2) = pG(:,2);

%{/
Hand_Para = Hand_Para_Matthew;
Shoulder_Para = Shoulder_Para_Matthew;
Rib_Para = Rib_Para_Matthew;
Waist_Para = Waist_Para_Matthew;

thHand_degree = rad2deg(thHand);
thShoulder_degree = rad2deg(thShoulder);
thRib_degree = rad2deg(thRib);
thWaist_degree = rad2deg(thWaist);
dthHand_degree = rad2deg(dthHand);
dthShoulder_degree = rad2deg(dthShoulder);
dthRib_degree = rad2deg(dthRib);
dthWaist_degree = rad2deg(dthWaist);

thHand_degree_VT = thHand_degree;
thHand_degree_VT(thHand_degree < Hand_Para.theta_PE_1_Ext) = Hand_Para.theta_PE_1_Ext;
thHand_degree_VT(thHand_degree > Hand_Para.theta_PE_1_Flex) = Hand_Para.theta_PE_1_Flex;


TauHand_MaxTorque = VoluntaryTorque(Hand_Para, thHand_degree_VT, dthHand_degree, 1) * 2;
TauHand_MinTorque = VoluntaryTorque(Hand_Para, thHand_degree_VT, dthHand_degree, -1) * 2;
TauHand_PassiveTorque = VoluntaryTorque(Hand_Para, thHand_degree_VT, dthHand_degree, 0) * 2;
TauShoulder_MaxTorque = VoluntaryTorque(Shoulder_Para, thShoulder_degree, dthShoulder_degree, 1) * 2;
TauShoulder_MinTorque = VoluntaryTorque(Shoulder_Para, thShoulder_degree, dthShoulder_degree, -1) * 2;
TauShoulder_PassiveTorque = VoluntaryTorque(Shoulder_Para, thShoulder_degree, dthShoulder_degree, 0) * 2;
TauRib_MaxTorque = VoluntaryTorque(Rib_Para, thRib_degree, dthRib_degree, 1) * 2;
TauRib_MinTorque = VoluntaryTorque(Rib_Para, thRib_degree, dthRib_degree, -1) * 2;
TauRib_PassiveTorque = VoluntaryTorque(Rib_Para, thRib_degree, dthRib_degree, 0) * 2;
TauWaist_MaxTorque = VoluntaryTorque(Waist_Para, thWaist_degree, dthWaist_degree, 1) * 2;
TauWaist_MinTorque = VoluntaryTorque(Waist_Para, thWaist_degree, dthWaist_degree, -1) * 2;
TauWaist_PassiveTorque = VoluntaryTorque(Waist_Para, thWaist_degree, dthWaist_degree, 0) * 2;

TauHand_ActivatedRate_Ext = (TauHand - TauHand_PassiveTorque) ./ (TauHand_MaxTorque - TauHand_PassiveTorque);
TauHand_ActivatedRate_Flex = (TauHand_PassiveTorque - TauHand) ./ (TauHand_PassiveTorque - TauHand_MinTorque);
TauShoulder_ActivatedRate_Ext = (TauShoulder - TauShoulder_PassiveTorque) ./ (TauShoulder_MaxTorque - TauShoulder_PassiveTorque);
TauShoulder_ActivatedRate_Flex = (TauShoulder_PassiveTorque - TauShoulder) ./ (TauShoulder_PassiveTorque - TauShoulder_MinTorque);
TauRib_ActivatedRate_Ext = (TauRib - TauRib_PassiveTorque) ./ (TauRib_MaxTorque - TauRib_PassiveTorque);
TauRib_ActivatedRate_Flex = (TauRib_PassiveTorque - TauRib) ./ (TauRib_PassiveTorque - TauRib_MinTorque);
TauWaist_ActivatedRate_Ext = (TauWaist - TauWaist_PassiveTorque) ./ (TauWaist_MaxTorque - TauWaist_PassiveTorque);
TauWaist_ActivatedRate_Flex = (TauWaist_PassiveTorque - TauWaist) ./ (TauWaist_PassiveTorque - TauWaist_MinTorque);

TauHand_ActivatedRate = TauHand_ActivatedRate_Ext;
TauHand_ActivatedRate(TauHand_ActivatedRate < 0) = -TauHand_ActivatedRate_Flex(TauHand_ActivatedRate < 0);

TauShoulder_ActivatedRate = TauShoulder_ActivatedRate_Ext;
TauShoulder_ActivatedRate(TauShoulder_ActivatedRate < 0) = -TauShoulder_ActivatedRate_Flex(TauShoulder_ActivatedRate < 0);

TauRib_ActivatedRate = TauRib_ActivatedRate_Ext;
TauRib_ActivatedRate(TauRib_ActivatedRate < 0) = -TauRib_ActivatedRate_Flex(TauRib_ActivatedRate < 0);

TauWaist_ActivatedRate = TauWaist_ActivatedRate_Ext;
TauWaist_ActivatedRate(TauWaist_ActivatedRate < 0) = -TauWaist_ActivatedRate_Flex(TauWaist_ActivatedRate < 0);

ActivatiedRate_Limit = 10;

TauHand_ActivatedRate(TauHand_ActivatedRate > ActivatiedRate_Limit) = ActivatiedRate_Limit;
TauShoulder_ActivatedRate(TauShoulder_ActivatedRate > ActivatiedRate_Limit) = ActivatiedRate_Limit;
TauRib_ActivatedRate(TauRib_ActivatedRate > ActivatiedRate_Limit) = ActivatiedRate_Limit;
TauWaist_ActivatedRate(TauWaist_ActivatedRate > ActivatiedRate_Limit) = ActivatiedRate_Limit;

TauHand_ActivatedRate(TauHand_ActivatedRate < -ActivatiedRate_Limit) = -ActivatiedRate_Limit;
TauShoulder_ActivatedRate(TauShoulder_ActivatedRate < -ActivatiedRate_Limit) = -ActivatiedRate_Limit;
TauRib_ActivatedRate(TauRib_ActivatedRate < -ActivatiedRate_Limit) = -ActivatiedRate_Limit;
TauWaist_ActivatedRate(TauWaist_ActivatedRate < -ActivatiedRate_Limit) = -ActivatiedRate_Limit;

%}




Fx_Hand = find_Fx_Hand(ddthRib,ddthHand,ddthWaist,ddthShoulder,dthRib,dthHand,dthWaist,dthShoulder,mArm,mLBody,mLeg,mUBody,rArm,rArmMCD,rLBody,rLBodyMCD,rLegMCD,rUBody,rUBodyMCD,thHand,thShoulder,thWaist,thRib);
Fy_Hand = find_Fy_Hand(ddrPB,ddthRib,ddthHand,ddthWaist,ddthShoulder,dthRib,dthHand,dthWaist,dthShoulder,g,mArm,mLBody,mLeg,mPB,mUBody,rArm,rArmMCD,rLBody,rLBodyMCD,rLegMCD,rUBody,rUBodyMCD,thHand,thShoulder,thWaist,thRib);

Torque_F_Hand = cross([pHand - pG, zeros([size(pG, 1), 1])], [Fx_Hand, Fy_Hand, zeros([size(pG, 1), 1])]);
Torque_F_Hand = Torque_F_Hand(:,3);

figure(29)
plot(Data_Time, momentumG)

figure(30)
plot(Data_Time, dmomentumG, 'DisplayName', '角運動量の微分')
hold on
plot(Data_Time, Torque_F_Hand, 'DisplayName', '外力のトルク')
hold off
legend

figure(31)
plot(Data_Time, dmomentumG - Torque_F_Hand, 'DisplayName', '角運動量の微分 - 外力のトルク')
hold on
plot(Data_Time, TauHand, 'DisplayName', '手首のトルク')
hold off
legend


figure(32)
plot(Data_Time, TauHand_ActivatedRate)
hold on
plot(Data_Time, TauShoulder_ActivatedRate)
plot(Data_Time, TauWaist_ActivatedRate)
plot(Data_Time, TauRib_ActivatedRate)
hold off


%{/

% thDataArray = [TauHand_ActivatedRate, TauShoulder_ActivatedRate, TauWaist_ActivatedRate];
% MDataArray = [FrPB];
% otherDataArray = [TauHand, TauShoulder, TauWaist];

% thDataArray = zeros(size(pHand(:,1)));
% MDataArray = zeros(size(pHand(:,1)));
% otherDataArray = zeros(size(pHand(:,1)));

thDataArray = FrPB;
MDataArray = [TauHand_ActivatedRate, TauShoulder_ActivatedRate, TauWaist_ActivatedRate, TauRib_ActivatedRate];
otherDataArray = zeros(size(pHand(:,1)));

sliderTimeh = 0.1;

Anime_Cortex_filt = Anime_Cortex(Data_Time, ...
    zeros([size(Data_Time, 1), 1, 2]), ...
    zeros([size(Data_Time, 1), 1, 2]), ...
    Body, ...
    pG_Anime, ...
    thDataArray, MDataArray, otherDataArray, sliderTimeh);
Anime_Cortex_filt.UIFigure.Name = 'Anime_Cortex_filt.UIFigure.Name';
% ylim(Anime_Cortex_filt.axThData,[-2, 2])
%}




















