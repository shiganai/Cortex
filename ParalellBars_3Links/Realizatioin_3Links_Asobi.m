


clear all

Bar_Data_All = (dlmread('Trimmed_Trimmed_ShSs3-Bar.trc','',6,1))/1000;%�P�ʕϊ�[mm]��[m]

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


Body_Data_All = (dlmread('Trimmed_Trimmed_ShSs3-51P.trc','',6,2))/1000;%�P�ʕϊ�[mm]��[m]

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
tro = (troR + troL)/2;
kn = (knR + knL)/2;
an = (anR + anL)/2;
toe = (toeR + toeL)/2;

mAll = 81;
r_Ankle_Toe = mean(vecnorm(an - toe, 2, 2));
r_Knee_Ankle = mean(vecnorm(kn - an, 2, 2));
r_Waist_Knee = mean(vecnorm(tro - kn, 2, 2));
r_Shoulder_Waist = mean(vecnorm(sh - tro, 2, 2));
r_Elbow_Shoulder = mean(vecnorm(elb - sh, 2, 2));
r_Wrist_Elbow = mean(vecnorm(wr - elb, 2, 2));
r_Finger_Wrist = mean(vecnorm(hand - wr, 2, 2));
r_Ear_Shoulder = mean(vecnorm(ear - sh, 2, 2));
r_Top_Ear = mean(vecnorm(top - ear, 2, 2));

wr_sh = sh - wr;
sh_tro = tro - sh;
tro_an = an - tro;

% rPB = hand(:,2) - hand(1,2);
hand_x = hand(:,1);
hand_y = hand(:,2);
thHand = GetThetaFromXY(wr, sh) - 1/2*pi;
thShoulder = GetThetaFromXY(sh, tro) - thHand - 1/2*pi;
thWaist = GetThetaFromXY(tro, an) - thHand - 1/2*pi - thShoulder;

thShoulder = thShoulder + 2*pi;
thWaist = thWaist - 2*pi;

thBody = GetThetaFromXY(sh, tro) - 1/2*pi;
thLeg = GetThetaFromXY(tro, an) - 1/2*pi;

thBody = thBody + 2*pi;
thLeg = thLeg + 0*pi;

rPB_raw = zeros(size(Data_Time));
linePBs = zeros(size(Data_Time, 1), 100);
for ii = 1:size(Data_Time, 1)
    linePB = spline([Bar0(ii,1), Bar50(ii,1), Bar100(ii,1), Bar130(ii,1), Bar180(ii,1), Bar230(ii,1)],...
        [Bar0(ii,2), Bar50(ii,2), Bar100(ii,2), Bar130(ii,2), Bar180(ii,2), Bar230(ii,2)], linspace(Bar0(ii,1), Bar230(ii, 1), 100));
    
    linePB = linePB - (Bar0(ii,2) + Bar230(ii,2))/2;
    [~, abs_Max_Index] = max(abs(linePB));
    rPB_raw(ii,1) = linePB(abs_Max_Index);
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

dthWaist = diff(thWaist) / diff(Data_Time(1:2));
ddthWaist = diff(dthWaist) / diff(Data_Time(1:2));

dhand_x = diff(hand_x) / diff(Data_Time(1:2));
ddhand_x = diff(dhand_x) / diff(Data_Time(1:2));

dhand_y = diff(hand_y) / diff(Data_Time(1:2));
ddhand_y = diff(dhand_y) / diff(Data_Time(1:2));

drPB_raw = diff(rPB_raw) / diff(Data_Time(1:2));
ddrPB_raw = diff(drPB_raw) / diff(Data_Time(1:2));

dthBody = diff(thBody) / diff(Data_Time(1:2));
ddthBody = diff(dthBody) / diff(Data_Time(1:2));

dthLeg = diff(thLeg) / diff(Data_Time(1:2));
ddthLeg = diff(dthLeg) / diff(Data_Time(1:2));

Cut_Off_Freq_tmp = (0.5:0.5:30)';
Nichest_Freq = 200 / 2;
% Cut_Range = find(Data_Time == 4.19):find(Data_Time == 7.95);
Cut_Range = find(Data_Time == 4.19):find(Data_Time == 8.5);


% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

RMS_ddthHand = zeros(size(Cut_Off_Freq_tmp));

for ii = 1:size(Cut_Off_Freq_tmp, 1)
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
    ddthHand_filt_tmp = filtfilt(Butter_b, Butter_a, ddthHand);
    
    RMS_ddthHand(ii) = (mean((ddthHand_filt_tmp - ddthHand).^2))^(1/2);
end

Cut_Off_Freq_Line_Range = find(Cut_Off_Freq_tmp == 6):find(Cut_Off_Freq_tmp == 24);
Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddthHand(Cut_Off_Freq_Line_Range), 1);

Cut_Off_Freq = interp1(RMS_ddthHand(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));

[Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
ddthHand_filt = filtfilt(Butter_b, Butter_a, ddthHand);



% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

RMS_ddthShoulder = zeros(size(Cut_Off_Freq_tmp));

for ii = 1:size(Cut_Off_Freq_tmp, 1)
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
    ddthShoulder_filt_tmp = filtfilt(Butter_b, Butter_a, ddthShoulder);
    
    RMS_ddthShoulder(ii) = (mean((ddthShoulder_filt_tmp - ddthShoulder).^2))^(1/2);
end

Cut_Off_Freq_Line_Range = find(Cut_Off_Freq_tmp == 3.5):find(Cut_Off_Freq_tmp == 9);
Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddthShoulder(Cut_Off_Freq_Line_Range), 1);

Cut_Off_Freq = interp1(RMS_ddthShoulder(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));

[Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
ddthShoulder_filt = filtfilt(Butter_b, Butter_a, ddthShoulder);

% Cut_Range = find(Data_Time == 3.5):find(Data_Time == 8.5);


% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

RMS_ddthWaist = zeros(size(Cut_Off_Freq_tmp));

for ii = 1:size(Cut_Off_Freq_tmp, 1)
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
    ddthWaist_filt_tmp = filtfilt(Butter_b, Butter_a, ddthWaist);
    
    RMS_ddthWaist(ii) = (mean((ddthWaist_filt_tmp - ddthWaist).^2))^(1/2);
end

Cut_Off_Freq_Line_Range = find(Cut_Off_Freq_tmp == 3.5):find(Cut_Off_Freq_tmp == 9);
Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddthWaist(Cut_Off_Freq_Line_Range), 1);

Cut_Off_Freq = interp1(RMS_ddthWaist(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));

[Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
ddthWaist_filt = filtfilt(Butter_b, Butter_a, ddthWaist);

% Cut_Range = find(Data_Time == 3.5):find(Data_Time == 8.5);

% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

RMS_ddhand_x = zeros(size(Cut_Off_Freq_tmp));

for ii = 1:size(Cut_Off_Freq_tmp, 1)
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
    ddhand_x_filt_tmp = filtfilt(Butter_b, Butter_a, ddhand_x);
    
    RMS_ddhand_x(ii) = (mean((ddhand_x_filt_tmp - ddhand_x).^2))^(1/2);
end

Cut_Off_Freq_Line_Range = find(Cut_Off_Freq_tmp == 10):find(Cut_Off_Freq_tmp == 30);
Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddhand_x(Cut_Off_Freq_Line_Range), 1);

Cut_Off_Freq = interp1(RMS_ddhand_x(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));

[Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
ddhand_x_filt = filtfilt(Butter_b, Butter_a, ddhand_x);

% Cut_Range = find(Data_Time == 3.5):find(Data_Time == 8.5);

% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

RMS_ddhand_y = zeros(size(Cut_Off_Freq_tmp));

for ii = 1:size(Cut_Off_Freq_tmp, 1)
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
    ddhand_y_filt_tmp = filtfilt(Butter_b, Butter_a, ddhand_y);
    
    RMS_ddhand_y(ii) = (mean((ddhand_y_filt_tmp - ddhand_y).^2))^(1/2);
end

% Cut_Off_Freq_Line_Range = find(Cut_Off_Freq_tmp == 3.5):find(Cut_Off_Freq_tmp == 8.5);
Cut_Off_Freq_Line_Range = find(Cut_Off_Freq_tmp == 10):find(Cut_Off_Freq_tmp == 30);
Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddhand_y(Cut_Off_Freq_Line_Range), 1);
Cut_Off_Freq = interp1(RMS_ddhand_y(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));

[Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
ddhand_y_filt = filtfilt(Butter_b, Butter_a, ddhand_y);

% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

RMS_ddrPB_raw = zeros(size(Cut_Off_Freq_tmp));

for ii = 1:size(Cut_Off_Freq_tmp, 1)
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
    ddrPB_raw_filt_tmp = filtfilt(Butter_b, Butter_a, ddrPB_raw);
    
    RMS_ddrPB_raw(ii) = (mean((ddrPB_raw_filt_tmp - ddrPB_raw).^2))^(1/2);
end

% Cut_Off_Freq_Line_Range = find(Cut_Off_Freq_tmp == 3.5):find(Cut_Off_Freq_tmp == 8.5);
Cut_Off_Freq_Line_Range = find(Cut_Off_Freq_tmp == 25):find(Cut_Off_Freq_tmp == 30);
Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddrPB_raw(Cut_Off_Freq_Line_Range), 1);

% Cut_Off_Freq = interp1(flip(RMS_ddrPB_raw(1:Cut_Off_Freq_Line_Range(1))), flip(Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1))), Cut_Off_Freq_Line(2))
Cut_Off_Freq = interp1(RMS_ddrPB_raw(find(Cut_Off_Freq_tmp == 7.5):Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(find(Cut_Off_Freq_tmp == 7.5):Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));

[Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
ddrPB_raw_filt = filtfilt(Butter_b, Butter_a, ddrPB_raw);

% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

RMS_ddthBody = zeros(size(Cut_Off_Freq_tmp));

for ii = 1:size(Cut_Off_Freq_tmp, 1)
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
    ddthBody_filt_tmp = filtfilt(Butter_b, Butter_a, ddthBody);
    
    RMS_ddthBody(ii) = (mean((ddthBody_filt_tmp - ddthBody).^2))^(1/2);
end

Cut_Off_Freq_Line_Range = find(abs(Cut_Off_Freq_tmp - 5) < 1e-3):find(abs(Cut_Off_Freq_tmp - 6) < 1e-3);
Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddthBody(Cut_Off_Freq_Line_Range), 1);

Cut_Off_Freq = interp1(flip(RMS_ddthBody(1:Cut_Off_Freq_Line_Range(1))), flip(Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1))), Cut_Off_Freq_Line(2));
% Cut_Off_Freq = interp1(RMS_ddthBody(find(Cut_Off_Freq_tmp == 7.5):Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(find(Cut_Off_Freq_tmp == 7.5):Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));

[Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
ddthBody_filt = filtfilt(Butter_b, Butter_a, ddthBody);

ddthBody_spline = spline(Data_Time(Cut_Range), ddthBody_filt(Cut_Range));
% dthBody_spline = fnint(ddthBody_spline, dthBody(Cut_Range(1)) * 0.95);
dthBody_spline = fnint(ddthBody_spline, mean(dthBody(Cut_Range(1:2))));
thBody_spline = fnint(dthBody_spline, mean(thBody(Cut_Range(1:2))));


% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

RMS_ddthLeg = zeros(size(Cut_Off_Freq_tmp));

for ii = 1:size(Cut_Off_Freq_tmp, 1)
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
    ddthLeg_filt_tmp = filtfilt(Butter_b, Butter_a, ddthLeg);
    
    RMS_ddthLeg(ii) = (mean((ddthLeg_filt_tmp - ddthLeg).^2))^(1/2);
end

Cut_Off_Freq_Line_Range = find(Cut_Off_Freq_tmp == 5):find(Cut_Off_Freq_tmp == 6);
Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddthLeg(Cut_Off_Freq_Line_Range), 1);

Cut_Off_Freq = interp1(flip(RMS_ddthLeg(1:Cut_Off_Freq_Line_Range(1))), flip(Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1))), Cut_Off_Freq_Line(2));
% Cut_Off_Freq = interp1(RMS_ddthLeg(find(Cut_Off_Freq_tmp == 7.5):Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(find(Cut_Off_Freq_tmp == 7.5):Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));

[Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
ddthLeg_filt = filtfilt(Butter_b, Butter_a, ddthLeg);

ddthLeg_spline = spline(Data_Time(Cut_Range), ddthLeg_filt(Cut_Range));
% dthLeg_spline = fnint(ddthLeg_spline, dthLeg(Cut_Range(1)) * 0.95);
dthLeg_spline = fnint(ddthLeg_spline, mean(dthLeg(Cut_Range(1:2))));
thLeg_spline = fnint(dthLeg_spline, mean(thLeg(Cut_Range(1:2))));


% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------


% Data_Time = Data_Time(Cut_Range) - Data_Time(Cut_Range(1));
Data_Time = Data_Time(Cut_Range);

ddthHand_spline = spline(Data_Time, ddthHand_filt(Cut_Range));
dthHand_spline = fnint(ddthHand_spline, 0);
% dthHand_spline = fnint(ddthHand_spline, mean(dthHand(Cut_Range(1:2))));
thHand_spline = fnint(dthHand_spline, mean(thHand(Cut_Range(1:2))));

ddthShoulder_spline = spline(Data_Time, ddthShoulder_filt(Cut_Range));
dthShoulder_spline = fnint(ddthShoulder_spline, mean(dthShoulder(Cut_Range(1:2))));
thShoulder_spline = fnint(dthShoulder_spline, mean(thShoulder(Cut_Range(1:2))));

ddthWaist_spline = spline(Data_Time, ddthWaist_filt(Cut_Range));
dthWaist_spline = fnint(ddthWaist_spline, mean(dthWaist(Cut_Range(1:2))));
thWaist_spline = fnint(dthWaist_spline, mean(thWaist(Cut_Range(1:2))));

ddhand_x_spline = spline(Data_Time, ddhand_x_filt(Cut_Range));
% dhand_x_spline = fnint(ddhand_x_spline, mean(dhand_x(Cut_Range(1:2))));
dhand_x_spline = fnint(ddhand_x_spline, dhand_x(Cut_Range(1)) * 1.15);
% dhand_x_spline = fnint(ddhand_x_spline, 0);
hand_x_spline = fnint(dhand_x_spline, mean(hand_x(Cut_Range(1:2))));

ddhand_y_spline = spline(Data_Time, ddhand_y_filt(Cut_Range));
dhand_y_spline = fnint(ddhand_y_spline, dhand_y(Cut_Range(1)) * 0.95);
% dhand_y_spline = fnint(ddhand_y_spline, mean(dhand_y(Cut_Range(1:2))));
hand_y_spline = fnint(dhand_y_spline, mean(hand_y(Cut_Range(1:2))));

ddrPB_raw_spline = spline(Data_Time, ddrPB_raw_filt(Cut_Range));
% drPB_raw_spline = fnint(ddrPB_raw_spline, drPB_raw(Cut_Range(1)) * 0.95);
drPB_raw_spline = fnint(ddrPB_raw_spline, mean(drPB_raw(Cut_Range(1:2))));
rPB_raw_spline = fnint(drPB_raw_spline, mean(rPB_raw(Cut_Range(1:2))));


% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

thShoulder_later = fnval(thBody_spline, Data_Time) - fnval(thHand_spline, Data_Time);
thWaist_later = fnval(thLeg_spline, Data_Time) - fnval(thHand_spline, Data_Time) - thShoulder_later;

% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

thShoulder_spline = spline(Data_Time, thShoulder_later);
dthShoulder_spline = fnder(thShoulder_spline, 1);
ddthShoulder_spline = fnder(thShoulder_spline, 2);

thWaist_spline = spline(Data_Time, thWaist_later);
dthWaist_spline = fnder(thWaist_spline, 1);
ddthWaist_spline = fnder(thWaist_spline, 2);


% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ddthHand_filt = fnval(ddthHand_spline, Data_Time);
dthHand_filt = fnval(dthHand_spline, Data_Time);
thHand_filt = fnval(thHand_spline, Data_Time);
ddthShoulder_filt = fnval(ddthShoulder_spline, Data_Time);
dthShoulder_filt = fnval(dthShoulder_spline, Data_Time);
thShoulder_filt = fnval(thShoulder_spline, Data_Time);
ddthWaist_filt = fnval(ddthWaist_spline, Data_Time);
dthWaist_filt = fnval(dthWaist_spline, Data_Time);
thWaist_filt = fnval(thWaist_spline, Data_Time);
ddhand_x_filt = fnval(ddhand_x_spline, Data_Time);
dhand_x_filt = fnval(dhand_x_spline, Data_Time);
hand_x_filt = fnval(hand_x_spline, Data_Time);
ddhand_y_filt = fnval(ddhand_y_spline, Data_Time);
dhand_y_filt = fnval(dhand_y_spline, Data_Time);
hand_y_filt = fnval(hand_y_spline, Data_Time);

rPB_filt = fnval(rPB_raw_spline, Data_Time);
drPB_filt = fnval(drPB_raw_spline, Data_Time);
ddrPB_filt = fnval(ddrPB_raw_spline, Data_Time);

rPB_filt_Hand = hand_y_filt;
drPB_filt_Hand = dhand_y_filt;
ddrPB_filt_Hand = ddhand_y_filt;


mean(rPB_filt(find(abs(Data_Time - 5) < 1/400):find(abs(Data_Time - 6) < 1/400)))






constants = constants_tmp(mAll, r_Ankle_Toe, r_Knee_Ankle, r_Waist_Knee, r_Shoulder_Waist, r_Elbow_Shoulder, r_Wrist_Elbow, r_Finger_Wrist, r_Ear_Shoulder, r_Top_Ear);

% constants.kPB = 1.9831 * 1e4 * 1.6;
% constants.kPB = 1.9831 * 1e4 * 3.2;

constants.kPB = 2.4160e+04;

% constants.cPB = constants.cPB * 2.4e2;
% constants.cPB = 10;
% constants.mPB = 5.1896;

constants.cPB = 8.29;
constants.mPB = 4.84 * 2;

g = constants.g;
kPB = constants.kPB;
cPB = constants.cPB;
mAll = constants.mAll;
mPB = constants.mPB;
mArm = constants.mArm;
mBody = constants.mBody;
mLeg = constants.mLeg;
rArm = constants.rArm;
rBody = constants.rBody;
rLeg = constants.rLeg;
InertiaModel = constants.InertiaModel;
rArmMCD = constants.rArmMCD;
rBodyMCD = constants.rBodyMCD;
rLegMCD = constants.rLegMCD;
InertiaArm = constants.InertiaArm;
InertiaBody = constants.InertiaBody;
InertiaLeg = constants.InertiaLeg;
InertiaG = constants.InertiaG;
Hand_Para = constants.Hand_Para;
Shoulder_Para = constants.Shoulder_Para;
Waist_Para = constants.Waist_Para;


% t = 6 : 1e-3 : Data_Time(end);
t = Data_Time(1) : 1e-3 : Data_Time(end);

q0_Index = find(abs(Data_Time - t(1)) < 1/400);
q0 = [-(mAll + mPB) * g/(2*kPB), thHand_filt(q0_Index), thShoulder_filt(q0_Index), thWaist_filt(q0_Index), ...
    drPB_filt(q0_Index), dthHand_filt(q0_Index), dthShoulder_filt(q0_Index), dthWaist_filt(q0_Index)]';

% q0 = [-(mAll + mPB) * g/(2*kPB), thHand_filt(1,1), thShoulder_filt(1,1), thWaist_filt(1,1), drPB_filt_PB(1,1), dthHand_filt(1,1), dthShoulder_filt(1,1), dthWaist_filt(1,1)]';

ddpp_time_histories = [ddthHand_spline, ddthShoulder_spline, ddthWaist_spline]';

ode = @(t, q) ddt_Realization(t, q, constants, ddpp_time_histories);
% eventFcn = @(t,q)Events(t,q,constants, ddpp_time_histories);
eventFcn = @(t,q)Events_NonStop(t,q,constants);

% odeOption = odeset('AbsTol', 1e-5);
odeOption = odeset('Events', eventFcn, 'AbsTol', 1e-5);

[time_onbar, q, te, ye, ie] = ode45(ode, t, q0, odeOption);

if ~isempty(ie)
    breakNum = ie(end);
else
    breakNum = 0;
end

Spin_num_Goal = 1.5;
ObjectiveValueData = ObjectiveFcn_OnlySpin(time_onbar, q, constants, Spin_num_Goal, breakNum, 0, Data_Time);
DispObjectiveValueData(ObjectiveValueData, 0)

% ObjectiveValueData = ObjectiveFcn(time_onbar, q, constants, 1.2, breakNum);
% DispObjectiveValueData(ObjectiveValueData, 0)

rPB = q(:,1);
thHand = q(:,2);
thShoulder = q(:,3);
thWaist = q(:,4);
drPB = q(:,5);
dthHand = q(:,6);
dthShoulder = q(:,7);
dthWaist = q(:,8);

yHand = rPB(:,1);
xHand = zeros(size(yHand));
dyHand = drPB(:,1);
dxHand = zeros(size(dyHand));

ddrPB = ppval(fnder(spline(time_onbar, drPB)), time_onbar);
ddthHand = ppval(fnder(spline(time_onbar, dthHand)), time_onbar);
ddthShoulder = ppval(fnder(spline(time_onbar, dthShoulder)), time_onbar);
ddthWaist = ppval(fnder(spline(time_onbar, dthWaist)), time_onbar);

MrPB = find_MrPB(ddrPB,ddthHand,ddthWaist,ddthShoulder,dthHand,dthWaist,dthShoulder,g,mArm,mBody,mLeg,mPB,rArm,rArmMCD,rBody,rBodyMCD,rLegMCD,thHand,thShoulder,thWaist);
MthHand = find_MthHand(InertiaLeg,InertiaArm,InertiaBody,ddrPB,ddthHand,ddthWaist,ddthShoulder,dthHand,dthWaist,dthShoulder,g,mArm,mBody,mLeg,rArm,rArmMCD,rBody,rBodyMCD,rLegMCD,thHand,thShoulder,thWaist);
MthShoulder = find_MthShoulder(InertiaLeg,InertiaBody,ddrPB,ddthHand,ddthWaist,ddthShoulder,dthHand,dthWaist,dthShoulder,g,mBody,mLeg,rArm,rBody,rBodyMCD,rLegMCD,thHand,thShoulder,thWaist);
MthWaist = find_MthWaist(InertiaLeg,ddrPB,ddthHand,ddthWaist,ddthShoulder,dthHand,dthShoulder,g,mLeg,rArm,rBody,rLegMCD,thHand,thShoulder,thWaist);

KE = find_KE(InertiaLeg,InertiaArm,InertiaBody,dthHand,dthWaist,dthShoulder,dxHand,dyHand,g,kPB,mArm,mBody,mLeg,mPB,rArm,rArmMCD,rBody,rBodyMCD,rLegMCD,thHand,thShoulder,thWaist,yHand);
% dKE = find_dKE(InertiaLeg,InertiaArm,InertiaBody,ddthHand,ddthWaist,ddthShoulder,ddxHand,ddyHand,dthHand,dthWaist,dthShoulder,dxHand,dyHand,g,kPB,mArm,mBody,mLeg,mPB,rArm,rArmMCD,rBody,rBodyMCD,rLegMCD,thHand,thShoulder,thWaist,yHand);
dKE = ppval(fnder(spline(time_onbar, KE)), time_onbar);

pG = find_pG(mArm,mBody,mLeg,rArm,rArmMCD,rBody,rBodyMCD,rLegMCD,thHand,thShoulder,thWaist,xHand,yHand);
vG = find_vG(dthHand,dthWaist,dthShoulder,dxHand,dyHand,mArm,mBody,mLeg,rArm,rArmMCD,rBody,rBodyMCD,rLegMCD,thHand,thShoulder,thWaist);
momentumG = find_momentumG(InertiaLeg,InertiaArm,InertiaBody,dthHand,dthWaist,dthShoulder,dxHand,dyHand,mArm,mBody,mLeg,rArm,rArmMCD,rBody,rBodyMCD,rLegMCD,thHand,thShoulder,thWaist,xHand,yHand);


pG_Anime(:,:,1) = pG(:,1);
pG_Anime(:,:,2) = pG(:,2);

height = pG(end,2) + vG(end,2)^2/(2*g) + 1.8 - rLeg;
time_until_landing = vG(end,2)/g + (2*height/g)^0.5; %vG(end,2)���}�C�i�X�������ꍇ�ł����v
omegaG = momentumG(end)/InertiaG;
spinNum = time_until_landing * omegaG /(2*pi)

thHand_degree = rad2deg(thHand);
thShoulder_degree = rad2deg(thShoulder);
thWaist_degree = rad2deg(thWaist);
dthHand_degree = rad2deg(dthHand);
dthShoulder_degree = rad2deg(dthShoulder);
dthWaist_degree = rad2deg(dthWaist);

thHand_degree_VT = thHand_degree;
thHand_degree_VT(thHand_degree < Hand_Para.theta_PE_1_Ext) = Hand_Para.theta_PE_1_Ext;
thHand_degree_VT(thHand_degree > Hand_Para.theta_PE_1_Flex) = Hand_Para.theta_PE_1_Flex;


MthHand_MaxTorque = VoluntaryTorque(Hand_Para, thHand_degree_VT, dthHand_degree, 1) * 2;
MthHand_MinTorque = VoluntaryTorque(Hand_Para, thHand_degree_VT, dthHand_degree, -1) * 2;
MthHand_PassiveTorque = VoluntaryTorque(Hand_Para, thHand_degree_VT, dthHand_degree, 0) * 2;
MthShoulder_MaxTorque = VoluntaryTorque(Shoulder_Para, thShoulder_degree, dthShoulder_degree, 1) * 2;
MthShoulder_MinTorque = VoluntaryTorque(Shoulder_Para, thShoulder_degree, dthShoulder_degree, -1) * 2;
MthShoulder_PassiveTorque = VoluntaryTorque(Shoulder_Para, thShoulder_degree, dthShoulder_degree, 0) * 2;
MthWaist_MaxTorque = VoluntaryTorque(Waist_Para, thWaist_degree, dthWaist_degree, 1) * 2;
MthWaist_MinTorque = VoluntaryTorque(Waist_Para, thWaist_degree, dthWaist_degree, -1) * 2;
MthWaist_PassiveTorque = VoluntaryTorque(Waist_Para, thWaist_degree, dthWaist_degree, 0) * 2;

MthHand_ActivatedRate_Ext = (MthHand - MthHand_PassiveTorque) ./ (MthHand_MaxTorque - MthHand_PassiveTorque);
MthHand_ActivatedRate_Flex = (MthHand_PassiveTorque - MthHand) ./ (MthHand_PassiveTorque - MthHand_MinTorque);
MthShoulder_ActivatedRate_Ext = (MthShoulder - MthShoulder_PassiveTorque) ./ (MthShoulder_MaxTorque - MthShoulder_PassiveTorque);
MthShoulder_ActivatedRate_Flex = (MthShoulder_PassiveTorque - MthShoulder) ./ (MthShoulder_PassiveTorque - MthShoulder_MinTorque);
MthWaist_ActivatedRate_Ext = (MthWaist - MthWaist_PassiveTorque) ./ (MthWaist_MaxTorque - MthWaist_PassiveTorque);
MthWaist_ActivatedRate_Flex = (MthWaist_PassiveTorque - MthWaist) ./ (MthWaist_PassiveTorque - MthWaist_MinTorque);

MthHand_ActivatedRate = MthHand_ActivatedRate_Ext;
MthHand_ActivatedRate(MthHand_ActivatedRate < 0) = -MthHand_ActivatedRate_Flex(MthHand_ActivatedRate < 0);

MthShoulder_ActivatedRate = MthShoulder_ActivatedRate_Ext;
MthShoulder_ActivatedRate(MthShoulder_ActivatedRate < 0) = -MthShoulder_ActivatedRate_Flex(MthShoulder_ActivatedRate < 0);

MthWaist_ActivatedRate = MthWaist_ActivatedRate_Ext;
MthWaist_ActivatedRate(MthWaist_ActivatedRate < 0) = -MthWaist_ActivatedRate_Flex(MthWaist_ActivatedRate < 0);

pHand = [zeros(size(rPB(:,1))), rPB];
vHand = [zeros(size(rPB(:,1))), drPB];

% pHand = [hand_x_filt, hand_y_filt];

pShoulder = pHand + rArm * [cos(thHand+1/2*pi), sin(thHand+1/2*pi)];
pWaist = pShoulder + rBody * [cos(thHand+1/2*pi + thShoulder), sin(thHand+1/2*pi + thShoulder)];
pToe = pWaist + rLeg * [cos(thHand+1/2*pi + thShoulder + thWaist), sin(thHand+1/2*pi + thShoulder + thWaist)];

time = time_onbar;

Body(:,:,1) = [pHand(:,1), pShoulder(:,1), pWaist(:,1), pToe(:,1)];
Body(:,:,2) = [pHand(:,2), pShoulder(:,2), pWaist(:,2), pToe(:,2)];

% Animation
%{
pG_Anime(:,:,1) = pG(:,1);
pG_Anime(:,:,2) = pG(:,2);

% thDataArray = [TauHand_ActivatedRate, TauShoulder_ActivatedRate, TauWaist_ActivatedRate];
% MDataArray = [FrPB];
% otherDataArray = [TauHand, TauShoulder, TauWaist];

% thDataArray = zeros(size(pHand(:,1)));
% MDataArray = zeros(size(pHand(:,1)));
% otherDataArray = zeros(size(pHand(:,1)));

thDataArray = MrPB;
MDataArray = [MthHand, MthShoulder, MthWaist];
otherDataArray = zeros(size(pHand(:,1)));

sliderTimeh = 0.1;

Anime_Cortex_filt = Anime_Cortex(time, ...
    zeros([size(time, 1), 1, 2]), ...
    zeros([size(time, 1), 1, 2]), ...
    Body, ...
    pG_Anime, ...
    thDataArray, MDataArray, otherDataArray, sliderTimeh);
Anime_Cortex_filt.UIFigure.Name = 'Anime_Cortex_filt.UIFigure.Name';
% ylim(Anime_Cortex_filt.axThData,[-2, 2])
xlim(Anime_Cortex_filt.axAnime,[-2.5, 2.5])
ylim(Anime_Cortex_filt.axAnime,[-2.5, 2.5])
%}


figure(31)
plot(Data_Time, rPB_filt - rPB_filt(1) + q0(1), 'DisplayName', 'raw', 'LineStyle', ':', 'LineWidth', 2)
hold on
% plot(Data_Time, rPB_filt_Hand - rPB_filt_Hand(1) + q0(1), 'DisplayName', 'raw\_FromHand', 'LineStyle', '--', 'LineWidth', 2)
plot(time, rPB, 'DisplayName', 'simulate', 'LineWidth', 2)
hold off
legend

figure(32)
plot(time, drPB)

RMS_rPB = mean((spline(Data_Time, rPB_filt - rPB_filt(1) + q0(1), time) - rPB).^2)
% RMS_rPB = mean((rPB_filt_PB - rPB_filt_PB(1) + q0(1) - spline(time, rPB, Data_Time)).^2)


%{
if breakNum == 1
    
    pG = find_pG(mArm,mBody,mLeg,rArm,rArmMCD,rBody,rBodyMCD,rLegMCD,thHand,thShoulder,thWaist,xHand,yHand);
    vG = find_vG(dthHand,dthWaist,dthShoulder,dxHand,dyHand,mArm,mBody,mLeg,rArm,rArmMCD,rBody,rBodyMCD,rLegMCD,thHand,thShoulder,thWaist);
    momentumG = find_momentumG(InertiaLeg,InertiaArm,InertiaBody,dthHand,dthWaist,dthShoulder,dxHand,dyHand,mArm,mBody,mLeg,rArm,rArmMCD,rBody,rBodyMCD,rLegMCD,thHand,thShoulder,thWaist,xHand,yHand);
    
    height = pG(end,2) + vG(end,2)^2/(2*g) + 1.8 - rLeg;
    time_until_landing = vG(end,2)/g + (2*height/g)^0.5; %vG(end,2)���}�C�i�X�������ꍇ�ł����v
    omegaG = momentumG(end)/InertiaG;
    spinNum = time_until_landing * omegaG /(2*pi)
    
    t_freefall = 0:1e-3:10;
    q0_freefall = [0, q(end,1), q(end,2), q(end,3), q(end,4), 0, q(end,5), q(end,6), q(end,7), q(end,8)]';
    
    eventFcn_freefall = @(t,q)Events_freefall(t, q, constants);
    odeOption_freefall = odeset('Events', eventFcn_freefall ,'AbsTol', 1e-5);
    ode_freefall = @(t,q) ddt_freefall(t, q, constants);
    
    [time_freefall, q_freefall] = ode45(ode_freefall, t_freefall, q0_freefall, odeOption_freefall);
    
    time = [time_onbar(:,1); time_freefall(2:end,1) + time_onbar(end,1)];
    
    xHand = [xHand; q_freefall(2:end, 1)];
    yHand = [yHand; q_freefall(2:end, 2)];
    thHand = [thHand; q_freefall(2:end,3)];
    thShoulder = [thShoulder; q_freefall(2:end,4)];
    thWaist = [thWaist; q_freefall(2:end,5)];
    
    dxHand = [dxHand; q_freefall(2:end, 6)];
    dyHand = [dyHand; q_freefall(2:end, 7)];
    dthHand = [dthHand; q_freefall(2:end,8)];
    dthShoulder = [dthShoulder; q_freefall(2:end,9)];
    dthWaist = [dthWaist; q_freefall(2:end,10)];
    
    MrPB = [MrPB; zeros(size(time_freefall(2:end,1)))];
    MthHand = [MthHand; zeros(size(time_freefall(2:end,1)))];
    MthShoulder = [MthShoulder; zeros(size(time_freefall(2:end,1)))];
    MthWaist = [MthWaist; zeros(size(time_freefall(2:end,1)))];
    
    pHand = [xHand, yHand];    
    vHand = [dxHand, dyHand];

    rPB = [rPB(:,1); zeros(size(time_freefall(2:end,1))) + rPB(end,1)];
    drPB = [drPB(:,1); zeros(size(time_freefall(2:end,1))) + drPB(end,1)];
    
    ddrPB = [ddrPB(:,1); zeros(size(time_freefall(2:end,1)))];
    ddthHand = [ddthHand(:,1); zeros(size(time_freefall(2:end,1)))];
    ddthShoulder = [ddthShoulder(:,1); zeros(size(time_freefall(2:end,1)))];
    ddthWaist = [ddthWaist(:,1); zeros(size(time_freefall(2:end,1)))];
    
else
    warning("���莸�s")
    time = time_onbar;
end

pShoulder = pHand + rArm * [cos(thHand+1/2*pi), sin(thHand+1/2*pi)];
pWaist = pShoulder + rBody * [cos(thHand+1/2*pi + thShoulder), sin(thHand+1/2*pi + thShoulder)];
pToe = pWaist + rLeg * [cos(thHand+1/2*pi + thShoulder + thWaist), sin(thHand+1/2*pi + thShoulder + thWaist)];

vShoulder = vHand + rArm * [-sin(thHand+1/2*pi), cos(thHand+1/2*pi)] .* dthHand;
vWaist = vShoulder + rBody * [-sin(thHand+1/2*pi + thShoulder), cos(thHand+1/2*pi + thShoulder)] .* (dthHand + dthShoulder);
vToe = vWaist + rLeg * [-sin(thHand+1/2*pi + thShoulder + thWaist), cos(thHand+1/2*pi + thShoulder + thWaist)] .* (dthHand + dthShoulder + dthWaist);

PrPB = MrPB .* drPB;
PthHand = MthHand .* dthHand;
PthShoulder = MthHand .* dthShoulder;
PthWaist = MthWaist .* dthWaist;

pG = find_pG(mArm,mBody,mLeg,rArm,rArmMCD,rBody,rBodyMCD,rLegMCD,thHand,thShoulder,thWaist,xHand,yHand);
vG = find_vG(dthHand,dthWaist,dthShoulder,dxHand,dyHand,mArm,mBody,mLeg,rArm,rArmMCD,rBody,rBodyMCD,rLegMCD,thHand,thShoulder,thWaist);
momentumG = find_momentumG(InertiaLeg,InertiaArm,InertiaBody,dthHand,dthWaist,dthShoulder,dxHand,dyHand,mArm,mBody,mLeg,rArm,rArmMCD,rBody,rBodyMCD,rLegMCD,thHand,thShoulder,thWaist,xHand,yHand);
KE = [KE; zeros([size(momentumG,1)-size(KE,1),1])];
dKE = [dKE; zeros([size(momentumG,1)-size(dKE,1),1])];

hold on
plot(time, rPB)
hold off
title("���莞�܂ł̕��s�_�̎��R������̕ψ�")
xlabel("time [s]")
ylabel("���R������̕ψ� [m]")
legend("���f�[�^","�V�~�����[�V����")
xlim([0 time_onbar(end)])


MthHand_ActivatedRate = zeros(size(time));
MthShoulder_ActivatedRate = zeros(size(time));
MthWaist_ActivatedRate = zeros(size(time));

Anime_main
breakNum


%}














































