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
rib = (ribR + ribL)/2;
tro = (troR + troL)/2;
kn = (knR + knL)/2;
an = (anR + anL)/2;
toe = (toeR + toeL)/2;

%{
mAll = 81;
r_Toe_Ankle = mean(vecnorm(an - toe, 2, 2));
r_Ankle_Knee = mean(vecnorm(kn - an, 2, 2));
r_Knee_Waist = mean(vecnorm(tro - kn, 2, 2));
r_Waist_Rib = mean(vecnorm(sh - tro, 2, 2));
r_Rib_Shoulder = mean(vecnorm(sh - tro, 2, 2));
r_Shoulder_Elbow = mean(vecnorm(elb - sh, 2, 2));
r_Elbow_Wrist = mean(vecnorm(wr - elb, 2, 2));
r_Wrist_Finger = mean(vecnorm(hand - wr, 2, 2));
r_Shoulder_Ear = mean(vecnorm(ear - sh, 2, 2));
r_Ear_Top = mean(vecnorm(top - ear, 2, 2));


% [mArm, mBody, mLeg, rArm, rBody, rLeg, rArmMCD, rBodyMCD, rLegMCD, InertiaArm, InertiaBody, InertiaLeg] = ...
%     Calc_Parameter_AE(mAll, r_Toe_Ankle, r_Ankle_Knee, r_Knee_Waist, r_Rib_Shoulder, r_Shoulder_Elbow, r_Elbow_Wrist, r_Wrist_Finger, r_Shoulder_Ear, r_Ear_Top);
[mArm, mUBody, mLBody, mLeg, rArm, rUBody, rLBody, rLeg, rArmMCD, rUBodyMCD, rLBodyMCD, rLegMCD, InertiaArm, InertiaUBody, InertiaLBody, InertiaLeg]...
    = Calc_Parameter_AE_4Links(mAll, r_Toe_Ankle, r_Ankle_Knee, r_Knee_Waist, r_Waist_Rib, r_Rib_Shoulder, r_Shoulder_Elbow, r_Elbow_Wrist, r_Wrist_Finger, r_Shoulder_Ear, r_Ear_Top);
%}

% rPB = hand(:,2) - hand(1,2);
hand_x = hand(:,1);
hand_y = hand(:,2);
thHand = GetThetaFromXY(wr, sh) - 1/2*pi;
thShoulder = GetThetaFromXY(sh, tro) - thHand - 1/2*pi;
thWaist = GetThetaFromXY(tro, an) - thHand - 1/2*pi - thShoulder;

thShoulder = thShoulder + 2*pi;
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

% Cut_Range_Start_Index_tmp = (find(Data_Time == 3.5):find(Data_Time == 6))';
Cut_Range_Start_Index_tmp = (find(Data_Time == 5):find(Data_Time == 6))';
Cut_Range_End_Time = 7.995;

V0_Ratio_List = (-2:0.1:2)';
Best_V0_Ratio_List_thHand = zeros(size(Cut_Range_Start_Index_tmp));
Best_V0_Ratio_List_thShoulder = zeros(size(Cut_Range_Start_Index_tmp));
Best_V0_Ratio_List_thWaist = zeros(size(Cut_Range_Start_Index_tmp));
Best_V0_Ratio_List_hand_x = zeros(size(Cut_Range_Start_Index_tmp));
Best_V0_Ratio_List_hand_y = zeros(size(Cut_Range_Start_Index_tmp));
Best_V0_Ratio_List_rPB = zeros(size(Cut_Range_Start_Index_tmp));

tic

% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
RMS_thHand_Cut_Range_Start_Index = zeros(size(Cut_Range_Start_Index_tmp));
Cut_Off_Freq_Line_Range = find(abs(Cut_Off_Freq_tmp - 5) < 1e-2):find(abs(Cut_Off_Freq_tmp - 6.5) < 1e-2);
for Cut_Range_Start_Index_tmp_Index = 1:size(Cut_Range_Start_Index_tmp, 1)
    
%     Cut_Range = Cut_Range_Start_Index_tmp(Cut_Range_Start_Index_tmp_Index):find(Data_Time == 7.95);
    Cut_Range = Cut_Range_Start_Index_tmp(Cut_Range_Start_Index_tmp_Index):find(Data_Time == Cut_Range_End_Time);
    
    RMS_ddthHand_Cut_Off_Freq = zeros(size(Cut_Off_Freq_tmp));
    
    for ii = 1:size(Cut_Off_Freq_tmp, 1)
        [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
        ddthHand_filt_tmp = filtfilt(Butter_b, Butter_a, ddthHand);
        
        RMS_ddthHand_Cut_Off_Freq(ii) = (mean((ddthHand_filt_tmp - ddthHand).^2))^(1/2);
    end
    
    Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddthHand_Cut_Off_Freq(Cut_Off_Freq_Line_Range), 1);
    
    Cut_Off_Freq = interp1(RMS_ddthHand_Cut_Off_Freq(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));
    
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
    ddthHand_filt = filtfilt(Butter_b, Butter_a, ddthHand);
    
    RMS_thHand_V0 = zeros(size(V0_Ratio_List));
    
    for V0_Ratio_Index = 1:size(V0_Ratio_List,1)
        ddthHand_spline = spline(Data_Time(Cut_Range), ddthHand_filt(Cut_Range));
        dthHand_spline = fnint(ddthHand_spline, mean(dthHand(Cut_Range(1:2))) * V0_Ratio_List(V0_Ratio_Index));
        thHand_spline = fnint(dthHand_spline, mean(thHand(Cut_Range(1:2))));
        
        RMS_thHand_V0(V0_Ratio_Index) = mean((fnval(thHand_spline, Data_Time(Cut_Range)) - thHand(Cut_Range)).^2)^(1/2);
        
    end
    
    [RMS_thHand_Cut_Range_Start_Index(Cut_Range_Start_Index_tmp_Index), Best_V0_Ratio_Index] = min(RMS_thHand_V0);
    Best_V0_Ratio_List_thHand(Cut_Range_Start_Index_tmp_Index) = V0_Ratio_List(Best_V0_Ratio_Index);
    
end

[min_value_thHand, min_index_thHand] = min(RMS_thHand_Cut_Range_Start_Index);
Best_Cut_Range_Start_Time_thHand = Data_Time(Cut_Range_Start_Index_tmp(min_index_thHand))
Best_V0_Ratio_thHand = Best_V0_Ratio_List_thHand(min_index_thHand)

Cut_Range = Cut_Range_Start_Index_tmp(min_index_thHand):find(Data_Time == Cut_Range_End_Time);
    
RMS_ddthHand_Cut_Off_Freq = zeros(size(Cut_Off_Freq_tmp));
    
for ii = 1:size(Cut_Off_Freq_tmp, 1)
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
    ddthHand_filt_tmp = filtfilt(Butter_b, Butter_a, ddthHand);
    
    RMS_ddthHand_Cut_Off_Freq(ii) = (mean((ddthHand_filt_tmp - ddthHand).^2))^(1/2);
end

Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddthHand_Cut_Off_Freq(Cut_Off_Freq_Line_Range), 1);

Cut_Off_Freq = interp1(RMS_ddthHand_Cut_Off_Freq(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));

[Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
ddthHand_filt = filtfilt(Butter_b, Butter_a, ddthHand);

ddthHand_spline = spline(Data_Time(Cut_Range), ddthHand_filt(Cut_Range));
dthHand_spline = fnint(ddthHand_spline, mean(dthHand(Cut_Range(1:2))) * Best_V0_Ratio_thHand);
thHand_spline = fnint(dthHand_spline, mean(thHand(Cut_Range(1:2))));

toc

% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

tic
RMS_thShoulder_Cut_Range_Start_Index = zeros(size(Cut_Range_Start_Index_tmp));
Cut_Off_Freq_Line_Range = find(abs(Cut_Off_Freq_tmp - 4) < 1e-2):find(abs(Cut_Off_Freq_tmp - 6) < 1e-2);
for Cut_Range_Start_Index_tmp_Index = 1:size(Cut_Range_Start_Index_tmp, 1)
    
%     Cut_Range = Cut_Range_Start_Index_tmp(Cut_Range_Start_Index_tmp_Index):find(Data_Time == 7.95);
    Cut_Range = Cut_Range_Start_Index_tmp(Cut_Range_Start_Index_tmp_Index):find(Data_Time == Cut_Range_End_Time);
    
    RMS_ddthShoulder_Cut_Off_Freq = zeros(size(Cut_Off_Freq_tmp));
    
    for ii = 1:size(Cut_Off_Freq_tmp, 1)
        [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
        ddthShoulder_filt_tmp = filtfilt(Butter_b, Butter_a, ddthShoulder);
        
        RMS_ddthShoulder_Cut_Off_Freq(ii) = (mean((ddthShoulder_filt_tmp - ddthShoulder).^2))^(1/2);
    end
    
    Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddthShoulder_Cut_Off_Freq(Cut_Off_Freq_Line_Range), 1);
    
    Cut_Off_Freq = interp1(RMS_ddthShoulder_Cut_Off_Freq(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));
    
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
    ddthShoulder_filt = filtfilt(Butter_b, Butter_a, ddthShoulder);
    
    RMS_thShoulder_V0 = zeros(size(V0_Ratio_List));
    
    for V0_Ratio_Index = 1:size(V0_Ratio_List,1)
        ddthShoulder_spline = spline(Data_Time(Cut_Range), ddthShoulder_filt(Cut_Range));
        dthShoulder_spline = fnint(ddthShoulder_spline, mean(dthShoulder(Cut_Range(1:2))) * V0_Ratio_List(V0_Ratio_Index));
        thShoulder_spline = fnint(dthShoulder_spline, mean(thShoulder(Cut_Range(1:2))));
        
        RMS_thShoulder_V0(V0_Ratio_Index) = mean((fnval(thShoulder_spline, Data_Time(Cut_Range)) - thShoulder(Cut_Range)).^2)^(1/2);
        
    end
    
    [RMS_thShoulder_Cut_Range_Start_Index(Cut_Range_Start_Index_tmp_Index), Best_V0_Ratio_Index] = min(RMS_thShoulder_V0);
    Best_V0_Ratio_List_thShoulder(Cut_Range_Start_Index_tmp_Index) = V0_Ratio_List(Best_V0_Ratio_Index);
    
end

[min_value_thShoulder, min_index_thShoulder] = min(RMS_thShoulder_Cut_Range_Start_Index);
Best_Cut_Range_Start_Time_thShoulder = Data_Time(Cut_Range_Start_Index_tmp(min_index_thShoulder))
Best_V0_Ratio_thShoulder = Best_V0_Ratio_List_thShoulder(min_index_thShoulder)

Cut_Range = Cut_Range_Start_Index_tmp(min_index_thShoulder):find(Data_Time == Cut_Range_End_Time);
    
RMS_ddthShoulder_Cut_Off_Freq = zeros(size(Cut_Off_Freq_tmp));
    
for ii = 1:size(Cut_Off_Freq_tmp, 1)
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
    ddthShoulder_filt_tmp = filtfilt(Butter_b, Butter_a, ddthShoulder);
    
    RMS_ddthShoulder_Cut_Off_Freq(ii) = (mean((ddthShoulder_filt_tmp - ddthShoulder).^2))^(1/2);
end

Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddthShoulder_Cut_Off_Freq(Cut_Off_Freq_Line_Range), 1);

Cut_Off_Freq = interp1(RMS_ddthShoulder_Cut_Off_Freq(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));

[Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
ddthShoulder_filt = filtfilt(Butter_b, Butter_a, ddthShoulder);

ddthShoulder_spline = spline(Data_Time(Cut_Range), ddthShoulder_filt(Cut_Range));
dthShoulder_spline = fnint(ddthShoulder_spline, mean(dthShoulder(Cut_Range(1:2))) * Best_V0_Ratio_thShoulder);
thShoulder_spline = fnint(dthShoulder_spline, mean(thShoulder(Cut_Range(1:2))));

toc

% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

tic
RMS_thWaist_Cut_Range_Start_Index = zeros(size(Cut_Range_Start_Index_tmp));
Cut_Off_Freq_Line_Range = find(abs(Cut_Off_Freq_tmp - 5) < 1e-2):find(abs(Cut_Off_Freq_tmp - 6) < 1e-2);
for Cut_Range_Start_Index_tmp_Index = 1:size(Cut_Range_Start_Index_tmp, 1)
    
%     Cut_Range = Cut_Range_Start_Index_tmp(Cut_Range_Start_Index_tmp_Index):find(Data_Time == 7.95);
    Cut_Range = Cut_Range_Start_Index_tmp(Cut_Range_Start_Index_tmp_Index):find(Data_Time == Cut_Range_End_Time);
    
    RMS_ddthWaist_Cut_Off_Freq = zeros(size(Cut_Off_Freq_tmp));
    
    for ii = 1:size(Cut_Off_Freq_tmp, 1)
        [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
        ddthWaist_filt_tmp = filtfilt(Butter_b, Butter_a, ddthWaist);
        
        RMS_ddthWaist_Cut_Off_Freq(ii) = (mean((ddthWaist_filt_tmp - ddthWaist).^2))^(1/2);
    end
    
    Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddthWaist_Cut_Off_Freq(Cut_Off_Freq_Line_Range), 1);
    
    Cut_Off_Freq = interp1(RMS_ddthWaist_Cut_Off_Freq(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));
    
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
    ddthWaist_filt = filtfilt(Butter_b, Butter_a, ddthWaist);
    
    RMS_thWaist_V0 = zeros(size(V0_Ratio_List));
    
    for V0_Ratio_Index = 1:size(V0_Ratio_List,1)
        ddthWaist_spline = spline(Data_Time(Cut_Range), ddthWaist_filt(Cut_Range));
        dthWaist_spline = fnint(ddthWaist_spline, mean(dthWaist(Cut_Range(1:2))) * V0_Ratio_List(V0_Ratio_Index));
        thWaist_spline = fnint(dthWaist_spline, mean(thWaist(Cut_Range(1:2))));
        
        RMS_thWaist_V0(V0_Ratio_Index) = mean((fnval(thWaist_spline, Data_Time(Cut_Range)) - thWaist(Cut_Range)).^2)^(1/2);
        
    end
    
    [RMS_thWaist_Cut_Range_Start_Index(Cut_Range_Start_Index_tmp_Index), Best_V0_Ratio_Index] = min(RMS_thWaist_V0);
    Best_V0_Ratio_List_thWaist(Cut_Range_Start_Index_tmp_Index) = V0_Ratio_List(Best_V0_Ratio_Index);
    
end

[min_value_thWaist, min_index_thWaist] = min(RMS_thWaist_Cut_Range_Start_Index);
Best_Cut_Range_Start_Time_thWaist = Data_Time(Cut_Range_Start_Index_tmp(min_index_thWaist))
Best_V0_Ratio_thWaist = Best_V0_Ratio_List_thWaist(min_index_thWaist)

Cut_Range = Cut_Range_Start_Index_tmp(min_index_thWaist):find(Data_Time == Cut_Range_End_Time);
    
RMS_ddthWaist_Cut_Off_Freq = zeros(size(Cut_Off_Freq_tmp));
    
for ii = 1:size(Cut_Off_Freq_tmp, 1)
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
    ddthWaist_filt_tmp = filtfilt(Butter_b, Butter_a, ddthWaist);
    
    RMS_ddthWaist_Cut_Off_Freq(ii) = (mean((ddthWaist_filt_tmp - ddthWaist).^2))^(1/2);
end

Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddthWaist_Cut_Off_Freq(Cut_Off_Freq_Line_Range), 1);

Cut_Off_Freq = interp1(RMS_ddthWaist_Cut_Off_Freq(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));

[Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
ddthWaist_filt = filtfilt(Butter_b, Butter_a, ddthWaist);

ddthWaist_spline = spline(Data_Time(Cut_Range), ddthWaist_filt(Cut_Range));
dthWaist_spline = fnint(ddthWaist_spline, mean(dthWaist(Cut_Range(1:2))) * Best_V0_Ratio_thWaist);
thWaist_spline = fnint(dthWaist_spline, mean(thWaist(Cut_Range(1:2))));

toc

% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

tic
RMS_hand_x_Cut_Range_Start_Index = zeros(size(Cut_Range_Start_Index_tmp));
Cut_Off_Freq_Line_Range = find(abs(Cut_Off_Freq_tmp - 5) < 1e-2):find(abs(Cut_Off_Freq_tmp - 10) < 1e-2);
for Cut_Range_Start_Index_tmp_Index = 1:size(Cut_Range_Start_Index_tmp, 1)
    
%     Cut_Range = Cut_Range_Start_Index_tmp(Cut_Range_Start_Index_tmp_Index):find(Data_Time == 7.95);
    Cut_Range = Cut_Range_Start_Index_tmp(Cut_Range_Start_Index_tmp_Index):find(Data_Time == Cut_Range_End_Time);
    
    RMS_ddhand_x_Cut_Off_Freq = zeros(size(Cut_Off_Freq_tmp));
    
    for ii = 1:size(Cut_Off_Freq_tmp, 1)
        [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
        ddhand_x_filt_tmp = filtfilt(Butter_b, Butter_a, ddhand_x);
        
        RMS_ddhand_x_Cut_Off_Freq(ii) = (mean((ddhand_x_filt_tmp - ddhand_x).^2))^(1/2);
    end
    
    Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddhand_x_Cut_Off_Freq(Cut_Off_Freq_Line_Range), 1);
    
    Cut_Off_Freq = interp1(RMS_ddhand_x_Cut_Off_Freq(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));
    
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
    ddhand_x_filt = filtfilt(Butter_b, Butter_a, ddhand_x);
    
    RMS_hand_x_V0 = zeros(size(V0_Ratio_List));
    
    for V0_Ratio_Index = 1:size(V0_Ratio_List,1)
        ddhand_x_spline = spline(Data_Time(Cut_Range), ddhand_x_filt(Cut_Range));
        dhand_x_spline = fnint(ddhand_x_spline, mean(dhand_x(Cut_Range(1:2))) * V0_Ratio_List(V0_Ratio_Index));
        hand_x_spline = fnint(dhand_x_spline, mean(hand_x(Cut_Range(1:2))));
        
        RMS_hand_x_V0(V0_Ratio_Index) = mean((fnval(hand_x_spline, Data_Time(Cut_Range)) - hand_x(Cut_Range)).^2)^(1/2);
        
    end
    
    [RMS_hand_x_Cut_Range_Start_Index(Cut_Range_Start_Index_tmp_Index), Best_V0_Ratio_Index] = min(RMS_hand_x_V0);
    Best_V0_Ratio_List_hand_x(Cut_Range_Start_Index_tmp_Index) = V0_Ratio_List(Best_V0_Ratio_Index);
    
end

[min_value_hand_x, min_index_hand_x] = min(RMS_hand_x_Cut_Range_Start_Index);
Best_Cut_Range_Start_Time_hand_x = Data_Time(Cut_Range_Start_Index_tmp(min_index_hand_x))
Best_V0_Ratio_hand_x = Best_V0_Ratio_List_hand_x(min_index_hand_x)

Cut_Range = Cut_Range_Start_Index_tmp(min_index_hand_x):find(Data_Time == Cut_Range_End_Time);
    
RMS_ddhand_x_Cut_Off_Freq = zeros(size(Cut_Off_Freq_tmp));
    
for ii = 1:size(Cut_Off_Freq_tmp, 1)
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
    ddhand_x_filt_tmp = filtfilt(Butter_b, Butter_a, ddhand_x);
    
    RMS_ddhand_x_Cut_Off_Freq(ii) = (mean((ddhand_x_filt_tmp - ddhand_x).^2))^(1/2);
end

Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddhand_x_Cut_Off_Freq(Cut_Off_Freq_Line_Range), 1);

Cut_Off_Freq = interp1(RMS_ddhand_x_Cut_Off_Freq(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));

[Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
ddhand_x_filt = filtfilt(Butter_b, Butter_a, ddhand_x);

ddhand_x_spline = spline(Data_Time(Cut_Range), ddhand_x_filt(Cut_Range));
dhand_x_spline = fnint(ddhand_x_spline, mean(dhand_x(Cut_Range(1:2))) * Best_V0_Ratio_hand_x);
hand_x_spline = fnint(dhand_x_spline, mean(hand_x(Cut_Range(1:2))));

toc

% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

tic
RMS_hand_y_Cut_Range_Start_Index = zeros(size(Cut_Range_Start_Index_tmp));
Cut_Off_Freq_Line_Range = find(abs(Cut_Off_Freq_tmp - 10) < 1e-2):find(abs(Cut_Off_Freq_tmp - 15) < 1e-2);
for Cut_Range_Start_Index_tmp_Index = 1:size(Cut_Range_Start_Index_tmp, 1)
    
%     Cut_Range = Cut_Range_Start_Index_tmp(Cut_Range_Start_Index_tmp_Index):find(Data_Time == 7.95);
    Cut_Range = Cut_Range_Start_Index_tmp(Cut_Range_Start_Index_tmp_Index):find(Data_Time == Cut_Range_End_Time);
    
    RMS_ddhand_y_Cut_Off_Freq = zeros(size(Cut_Off_Freq_tmp));
    
    for ii = 1:size(Cut_Off_Freq_tmp, 1)
        [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
        ddhand_y_filt_tmp = filtfilt(Butter_b, Butter_a, ddhand_y);
        
        RMS_ddhand_y_Cut_Off_Freq(ii) = (mean((ddhand_y_filt_tmp - ddhand_y).^2))^(1/2);
    end
    
    Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddhand_y_Cut_Off_Freq(Cut_Off_Freq_Line_Range), 1);
    
    Cut_Off_Freq = interp1(RMS_ddhand_y_Cut_Off_Freq(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));
    
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
    ddhand_y_filt = filtfilt(Butter_b, Butter_a, ddhand_y);
    
    RMS_hand_y_V0 = zeros(size(V0_Ratio_List));
    
    for V0_Ratio_Index = 1:size(V0_Ratio_List,1)
        ddhand_y_spline = spline(Data_Time(Cut_Range), ddhand_y_filt(Cut_Range));
        dhand_y_spline = fnint(ddhand_y_spline, mean(dhand_y(Cut_Range(1:2))) * V0_Ratio_List(V0_Ratio_Index));
        hand_y_spline = fnint(dhand_y_spline, mean(hand_y(Cut_Range(1:2))));
        
        RMS_hand_y_V0(V0_Ratio_Index) = mean((fnval(hand_y_spline, Data_Time(Cut_Range)) - hand_y(Cut_Range)).^2)^(1/2);
        
    end
    
    [RMS_hand_y_Cut_Range_Start_Index(Cut_Range_Start_Index_tmp_Index), Best_V0_Ratio_Index] = min(RMS_hand_y_V0);
    Best_V0_Ratio_List_hand_y(Cut_Range_Start_Index_tmp_Index) = V0_Ratio_List(Best_V0_Ratio_Index);
    
end

[min_value_hand_y, min_index_hand_y] = min(RMS_hand_y_Cut_Range_Start_Index);
Best_Cut_Range_Start_Time_hand_y = Data_Time(Cut_Range_Start_Index_tmp(min_index_hand_y))
Best_V0_Ratio_hand_y = Best_V0_Ratio_List_hand_y(min_index_hand_y)

Cut_Range = Cut_Range_Start_Index_tmp(min_index_hand_y):find(Data_Time == Cut_Range_End_Time);
    
RMS_ddhand_y_Cut_Off_Freq = zeros(size(Cut_Off_Freq_tmp));
    
for ii = 1:size(Cut_Off_Freq_tmp, 1)
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
    ddhand_y_filt_tmp = filtfilt(Butter_b, Butter_a, ddhand_y);
    
    RMS_ddhand_y_Cut_Off_Freq(ii) = (mean((ddhand_y_filt_tmp - ddhand_y).^2))^(1/2);
end

Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddhand_y_Cut_Off_Freq(Cut_Off_Freq_Line_Range), 1);

Cut_Off_Freq = interp1(RMS_ddhand_y_Cut_Off_Freq(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));

[Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
ddhand_y_filt = filtfilt(Butter_b, Butter_a, ddhand_y);

ddhand_y_spline = spline(Data_Time(Cut_Range), ddhand_y_filt(Cut_Range));
dhand_y_spline = fnint(ddhand_y_spline, mean(dhand_y(Cut_Range(1:2))) * Best_V0_Ratio_hand_y);
hand_y_spline = fnint(dhand_y_spline, mean(hand_y(Cut_Range(1:2))));

toc

% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

tic
RMS_rPB_Cut_Range_Start_Index = zeros(size(Cut_Range_Start_Index_tmp));
Cut_Off_Freq_Line_Range = find(abs(Cut_Off_Freq_tmp - 25) < 1e-2):find(abs(Cut_Off_Freq_tmp - 30) < 1e-2);
for Cut_Range_Start_Index_tmp_Index = 1:size(Cut_Range_Start_Index_tmp, 1)
    
%     Cut_Range = Cut_Range_Start_Index_tmp(Cut_Range_Start_Index_tmp_Index):find(Data_Time == 7.95);
    Cut_Range = Cut_Range_Start_Index_tmp(Cut_Range_Start_Index_tmp_Index):find(Data_Time == Cut_Range_End_Time);
    
    RMS_ddrPB_Cut_Off_Freq = zeros(size(Cut_Off_Freq_tmp));
    
    for ii = 1:size(Cut_Off_Freq_tmp, 1)
        [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
        ddrPB_filt_tmp = filtfilt(Butter_b, Butter_a, ddrPB);
        
        RMS_ddrPB_Cut_Off_Freq(ii) = (mean((ddrPB_filt_tmp - ddrPB).^2))^(1/2);
    end
    
    Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddrPB_Cut_Off_Freq(Cut_Off_Freq_Line_Range), 1);
    
    Cut_Off_Freq = interp1(RMS_ddrPB_Cut_Off_Freq(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));
    
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
    ddrPB_filt = filtfilt(Butter_b, Butter_a, ddrPB);
    
    RMS_rPB_V0 = zeros(size(V0_Ratio_List));
    
    for V0_Ratio_Index = 1:size(V0_Ratio_List,1)
        ddrPB_spline = spline(Data_Time(Cut_Range), ddrPB_filt(Cut_Range));
        drPB_spline = fnint(ddrPB_spline, mean(drPB(Cut_Range(1:2))) * V0_Ratio_List(V0_Ratio_Index));
        rPB_spline = fnint(drPB_spline, mean(rPB(Cut_Range(1:2))));
        
        RMS_rPB_V0(V0_Ratio_Index) = mean((fnval(rPB_spline, Data_Time(Cut_Range)) - rPB(Cut_Range)).^2)^(1/2);
        
    end
    
    [RMS_rPB_Cut_Range_Start_Index(Cut_Range_Start_Index_tmp_Index), Best_V0_Ratio_Index] = min(RMS_rPB_V0);
    Best_V0_Ratio_List_rPB(Cut_Range_Start_Index_tmp_Index) = V0_Ratio_List(Best_V0_Ratio_Index);
    
end

[min_value_rPB, min_index_rPB] = min(RMS_rPB_Cut_Range_Start_Index);
Best_Cut_Range_Start_Time_rPB = Data_Time(Cut_Range_Start_Index_tmp(min_index_rPB))
Best_V0_Ratio_rPB = Best_V0_Ratio_List_rPB(min_index_rPB)

Cut_Range = Cut_Range_Start_Index_tmp(min_index_rPB):find(Data_Time == Cut_Range_End_Time);
    
RMS_ddrPB_Cut_Off_Freq = zeros(size(Cut_Off_Freq_tmp));
    
for ii = 1:size(Cut_Off_Freq_tmp, 1)
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
    ddrPB_filt_tmp = filtfilt(Butter_b, Butter_a, ddrPB);
    
    RMS_ddrPB_Cut_Off_Freq(ii) = (mean((ddrPB_filt_tmp - ddrPB).^2))^(1/2);
end

Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddrPB_Cut_Off_Freq(Cut_Off_Freq_Line_Range), 1);

Cut_Off_Freq = interp1(RMS_ddrPB_Cut_Off_Freq(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));

[Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
ddrPB_filt = filtfilt(Butter_b, Butter_a, ddrPB);

ddrPB_spline = spline(Data_Time(Cut_Range), ddrPB_filt(Cut_Range));
drPB_spline = fnint(ddrPB_spline, mean(drPB(Cut_Range(1:2))) * Best_V0_Ratio_rPB);
rPB_spline = fnint(drPB_spline, mean(rPB(Cut_Range(1:2))));

toc

% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

tic

figure_num = 0;

% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
figure_num = figure_num + 1; 
figure(figure_num)
plot(Data_Time(Cut_Range_Start_Index_tmp), RMS_thHand_Cut_Range_Start_Index)
title('RMS_thHand_Cut_Range_Start_Index')

figure_num = figure_num + 1; 
figure(figure_num)
plot(Data_Time, thHand, 'DisplayName', 'Raw', 'LineStyle', ':', 'LineWidth', 3)
hold on
plot(Data_Time(Cut_Range), fnval(thHand_spline, Data_Time(Cut_Range)), 'DisplayName', 'Raw', 'LineStyle', '-', 'LineWidth', 2)
hold off
title('thHand')

% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
figure_num = figure_num + 1; 
figure(figure_num)
plot(Data_Time(Cut_Range_Start_Index_tmp), RMS_thShoulder_Cut_Range_Start_Index)
title('RMS_thShoulder_Cut_Range_Start_Index')

figure_num = figure_num + 1; 
figure(figure_num)
plot(Data_Time, thShoulder, 'DisplayName', 'Raw', 'LineStyle', ':', 'LineWidth', 3)
hold on
plot(Data_Time(Cut_Range), fnval(thShoulder_spline, Data_Time(Cut_Range)), 'DisplayName', 'Raw', 'LineStyle', '-', 'LineWidth', 2)
hold off
title('thShoulder')

% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
figure_num = figure_num + 1; 
figure(figure_num)
plot(Data_Time(Cut_Range_Start_Index_tmp), RMS_thWaist_Cut_Range_Start_Index)
title('RMS_thWaist_Cut_Range_Start_Index')

figure_num = figure_num + 1; 
figure(figure_num)
plot(Data_Time, thWaist, 'DisplayName', 'Raw', 'LineStyle', ':', 'LineWidth', 3)
hold on
plot(Data_Time(Cut_Range), fnval(thWaist_spline, Data_Time(Cut_Range)), 'DisplayName', 'Raw', 'LineStyle', '-', 'LineWidth', 2)
hold off
title('thWaist')

% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
figure_num = figure_num + 1; 
figure(figure_num)
plot(Data_Time(Cut_Range_Start_Index_tmp), RMS_hand_x_Cut_Range_Start_Index)
title('RMS_hand_x_Cut_Range_Start_Index')

figure_num = figure_num + 1; 
figure(figure_num)
plot(Data_Time, hand_x, 'DisplayName', 'Raw', 'LineStyle', ':', 'LineWidth', 3)
hold on
plot(Data_Time(Cut_Range), fnval(hand_x_spline, Data_Time(Cut_Range)), 'DisplayName', 'Raw', 'LineStyle', '-', 'LineWidth', 2)
hold off
title('hand_x')

% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
figure_num = figure_num + 1; 
figure(figure_num)
plot(Data_Time(Cut_Range_Start_Index_tmp), RMS_hand_y_Cut_Range_Start_Index)
title('RMS_hand_y_Cut_Range_Start_Index')

figure_num = figure_num + 1; 
figure(figure_num)
plot(Data_Time, hand_y, 'DisplayName', 'Raw', 'LineStyle', ':', 'LineWidth', 3)
hold on
plot(Data_Time(Cut_Range), fnval(hand_y_spline, Data_Time(Cut_Range)), 'DisplayName', 'Raw', 'LineStyle', '-', 'LineWidth', 2)
hold off
title('hand_y')

% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
figure_num = figure_num + 1; 
figure(figure_num)
plot(Data_Time(Cut_Range_Start_Index_tmp), RMS_rPB_Cut_Range_Start_Index)
title('RMS_rPB_Cut_Range_Start_Index')

figure_num = figure_num + 1; 
figure(figure_num)
plot(Data_Time, rPB, 'DisplayName', 'Raw', 'LineStyle', ':', 'LineWidth', 3)
hold on
plot(Data_Time(Cut_Range), fnval(rPB_spline, Data_Time(Cut_Range)), 'DisplayName', 'Raw', 'LineStyle', '-', 'LineWidth', 2)
hold off
title('rPB')


toc
    

%{
for Cut_Range_Start_Index_tmp_Index = 1:size(Cut_Range_Start_Index_tmp, 1)
    
    
    RMS_ddthShoulder_Cut_Off_Freq = zeros(size(Cut_Off_Freq_tmp));
    
    for ii = 1:size(Cut_Off_Freq_tmp, 1)
        [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
        ddthShoulder_filt_tmp = filtfilt(Butter_b, Butter_a, ddthShoulder);
        
        RMS_ddthShoulder_Cut_Off_Freq(ii) = (mean((ddthShoulder_filt_tmp - ddthShoulder).^2))^(1/2);
    end
    
    Cut_Off_Freq_Line_Range = find(abs(Cut_Off_Freq_tmp - 4) < 1e-2):find(abs(Cut_Off_Freq_tmp - 6) < 1e-2);
    % Cut_Off_Freq_Line_Range = find(Cut_Off_Freq_tmp == 3.5):find(Cut_Off_Freq_tmp == 9);
    Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddthShoulder_Cut_Off_Freq(Cut_Off_Freq_Line_Range), 1);
    
    Cut_Off_Freq = interp1(RMS_ddthShoulder_Cut_Off_Freq(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));
    
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
    ddthShoulder_filt = filtfilt(Butter_b, Butter_a, ddthShoulder);
    
    RMS_thShoulder_V0 = zeros(size(V0_Ratio_List));
    
    for V0_Ratio_Index = 1:size(V0_Ratio_List,1)
        ddthShoulder_spline = spline(Data_Time(Cut_Range), ddthShoulder_filt(Cut_Range));
%         dthShoulder_spline = fnint(ddthShoulder_spline, 0);
        dthShoulder_spline = fnint(ddthShoulder_spline, mean(dthShoulder(Cut_Range(1:2))) * V0_Ratio_List(V0_Ratio_Index));
        thShoulder_spline = fnint(dthShoulder_spline, mean(thShoulder(Cut_Range(1:2))));
        
        RMS_thShoulder_V0(V0_Ratio_Index) = mean((fnval(thShoulder_spline, Data_Time(Cut_Range)) - thShoulder(Cut_Range)).^2)^(1/2);
        
    end
    
    [RMS_thShoulder_tmp, Best_V0_Ratio_Index] = min(RMS_thShoulder_V0);
    Best_V0_Ratio_List_thShoulder(Cut_Range_Start_Index_tmp_Index) = V0_Ratio_List(Best_V0_Ratio_Index);
    
    
    RMS_ddthWaist_Cut_Off_Freq = zeros(size(Cut_Off_Freq_tmp));
    
    for ii = 1:size(Cut_Off_Freq_tmp, 1)
        [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
        ddthWaist_filt_tmp = filtfilt(Butter_b, Butter_a, ddthWaist);
        
        RMS_ddthWaist_Cut_Off_Freq(ii) = (mean((ddthWaist_filt_tmp - ddthWaist).^2))^(1/2);
    end
    
    Cut_Off_Freq_Line_Range = find(abs(Cut_Off_Freq_tmp - 5) < 1e-2):find(abs(Cut_Off_Freq_tmp - 6) < 1e-2);
    % Cut_Off_Freq_Line_Range = find(Cut_Off_Freq_tmp == 3.5):find(Cut_Off_Freq_tmp == 9);
    Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddthWaist_Cut_Off_Freq(Cut_Off_Freq_Line_Range), 1);
    
    Cut_Off_Freq = interp1(RMS_ddthWaist_Cut_Off_Freq(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));
    
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
    ddthWaist_filt = filtfilt(Butter_b, Butter_a, ddthWaist);
    
    RMS_thWaist_V0 = zeros(size(V0_Ratio_List));
    
    for V0_Ratio_Index = 1:size(V0_Ratio_List,1)
        ddthWaist_spline = spline(Data_Time(Cut_Range), ddthWaist_filt(Cut_Range));
%         dthWaist_spline = fnint(ddthWaist_spline, 0);
        dthWaist_spline = fnint(ddthWaist_spline, mean(dthWaist(Cut_Range(1:2))) * V0_Ratio_List(V0_Ratio_Index));
        thWaist_spline = fnint(dthWaist_spline, mean(thWaist(Cut_Range(1:2))));
        
        RMS_thWaist_V0(V0_Ratio_Index) = mean((fnval(thWaist_spline, Data_Time(Cut_Range)) - thWaist(Cut_Range)).^2)^(1/2);
        
    end
    
    [RMS_thWaist_tmp, Best_V0_Ratio_Index] = min(RMS_thWaist_V0);
    Best_V0_Ratio_List_thWaist(Cut_Range_Start_Index_tmp_Index) = V0_Ratio_List(Best_V0_Ratio_Index);
    
    
    RMS_ddhand_x_Cut_Off_Freq = zeros(size(Cut_Off_Freq_tmp));
    
    for ii = 1:size(Cut_Off_Freq_tmp, 1)
        [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
        ddhand_x_filt_tmp = filtfilt(Butter_b, Butter_a, ddhand_x);
        
        RMS_ddhand_x_Cut_Off_Freq(ii) = (mean((ddhand_x_filt_tmp - ddhand_x).^2))^(1/2);
    end
    
    Cut_Off_Freq_Line_Range = find(abs(Cut_Off_Freq_tmp - 5) < 1e-2):find(abs(Cut_Off_Freq_tmp - 10) < 1e-2);
    % Cut_Off_Freq_Line_Range = find(Cut_Off_Freq_tmp == 10):find(Cut_Off_Freq_tmp == 30);
    Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddhand_x_Cut_Off_Freq(Cut_Off_Freq_Line_Range), 1);
    
    Cut_Off_Freq = interp1(RMS_ddhand_x_Cut_Off_Freq(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));
    
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
    ddhand_x_filt = filtfilt(Butter_b, Butter_a, ddhand_x);
    
    RMS_hand_x_V0 = zeros(size(V0_Ratio_List));
    
    for V0_Ratio_Index = 1:size(V0_Ratio_List,1)
        ddhand_x_spline = spline(Data_Time(Cut_Range), ddhand_x_filt(Cut_Range));
%         dhand_x_spline = fnint(ddhand_x_spline, 0);
        dhand_x_spline = fnint(ddhand_x_spline, mean(dhand_x(Cut_Range(1:2))) * V0_Ratio_List(V0_Ratio_Index));
        hand_x_spline = fnint(dhand_x_spline, mean(hand_x(Cut_Range(1:2))));
        
        RMS_hand_x_V0(V0_Ratio_Index) = mean((fnval(hand_x_spline, Data_Time(Cut_Range)) - hand_x(Cut_Range)).^2)^(1/2);
        
    end
    
    [RMS_hand_x_tmp, Best_V0_Ratio_Index] = min(RMS_hand_x_V0);
    Best_V0_Ratio_List_hand_x(Cut_Range_Start_Index_tmp_Index) = V0_Ratio_List(Best_V0_Ratio_Index);
    
    
    RMS_ddhand_y_Cut_Off_Freq = zeros(size(Cut_Off_Freq_tmp));
    
    for ii = 1:size(Cut_Off_Freq_tmp, 1)
        [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
        ddhand_y_filt_tmp = filtfilt(Butter_b, Butter_a, ddhand_y);
        
        RMS_ddhand_y_Cut_Off_Freq(ii) = (mean((ddhand_y_filt_tmp - ddhand_y).^2))^(1/2);
    end
    
    Cut_Off_Freq_Line_Range = find(abs(Cut_Off_Freq_tmp - 10) < 1e-2):find(abs(Cut_Off_Freq_tmp - 15) < 1e-2);
    % Cut_Off_Freq_Line_Range = find(Cut_Off_Freq_tmp == 3.5):find(Cut_Off_Freq_tmp == 8.5);
    % Cut_Off_Freq_Line_Range = find(Cut_Off_Freq_tmp == 10):find(Cut_Off_Freq_tmp == 30);
    Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddhand_y_Cut_Off_Freq(Cut_Off_Freq_Line_Range), 1);
    
    Cut_Off_Freq = interp1(RMS_ddhand_y_Cut_Off_Freq(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));
    
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
    ddhand_y_filt = filtfilt(Butter_b, Butter_a, ddhand_y);
    
    RMS_hand_y_V0 = zeros(size(V0_Ratio_List));
    
    for V0_Ratio_Index = 1:size(V0_Ratio_List,1)
        ddhand_y_spline = spline(Data_Time(Cut_Range), ddhand_y_filt(Cut_Range));
%         dhand_y_spline = fnint(ddhand_y_spline, 0);
        dhand_y_spline = fnint(ddhand_y_spline, mean(dhand_y(Cut_Range(1:2))) * V0_Ratio_List(V0_Ratio_Index));
        hand_y_spline = fnint(dhand_y_spline, mean(hand_y(Cut_Range(1:2))));
        
        RMS_hand_y_V0(V0_Ratio_Index) = mean((fnval(hand_y_spline, Data_Time(Cut_Range)) - hand_y(Cut_Range)).^2)^(1/2);
        
    end
    
    [RMS_hand_y_tmp, Best_V0_Ratio_Index] = min(RMS_hand_y_V0);
    Best_V0_Ratio_List_hand_y(Cut_Range_Start_Index_tmp_Index) = V0_Ratio_List(Best_V0_Ratio_Index);
    
    
    
    RMS_ddrPB_Cut_Off_Freq = zeros(size(Cut_Off_Freq_tmp));
    
    for ii = 1:size(Cut_Off_Freq_tmp, 1)
        [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
        ddrPB_filt_tmp = filtfilt(Butter_b, Butter_a, ddrPB);
        
        RMS_ddrPB_Cut_Off_Freq(ii) = (mean((ddrPB_filt_tmp - ddrPB).^2))^(1/2);
    end
    
    Cut_Off_Freq_Line_Range = find(abs(Cut_Off_Freq_tmp - 25) < 1e-2):find(abs(Cut_Off_Freq_tmp - 30) < 1e-2);
    % Cut_Off_Freq_Line_Range = find(Cut_Off_Freq_tmp == 3.5):find(Cut_Off_Freq_tmp == 8.5);
    % Cut_Off_Freq_Line_Range = find(Cut_Off_Freq_tmp == 25):find(Cut_Off_Freq_tmp == 30);
    Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddrPB_Cut_Off_Freq(Cut_Off_Freq_Line_Range), 1);
    
    % Cut_Off_Freq = interp1(flip(RMS_ddrPB(1:Cut_Off_Freq_Line_Range(1))), flip(Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1))), Cut_Off_Freq_Line(2))
    Cut_Off_Freq = interp1(RMS_ddrPB_Cut_Off_Freq(find(Cut_Off_Freq_tmp == 7.5):Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(find(Cut_Off_Freq_tmp == 7.5):Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));
    
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
    ddrPB_filt = filtfilt(Butter_b, Butter_a, ddrPB);
    
    RMS_rPB_V0 = zeros(size(V0_Ratio_List));
    
    for V0_Ratio_Index = 1:size(V0_Ratio_List,1)
        ddrPB_spline = spline(Data_Time(Cut_Range), ddrPB_filt(Cut_Range));
%         drPB_spline = fnint(ddrPB_spline, 0);
        drPB_spline = fnint(ddrPB_spline, mean(drPB(Cut_Range(1:2))) * V0_Ratio_List(V0_Ratio_Index));
        rPB_spline = fnint(drPB_spline, mean(rPB(Cut_Range(1:2))));
        
        RMS_rPB_V0(V0_Ratio_Index) = mean((fnval(rPB_spline, Data_Time(Cut_Range)) - rPB(Cut_Range)).^2)^(1/2);
        
    end
    
    [RMS_rPB_tmp, Best_V0_Ratio_Index] = min(RMS_rPB_V0);
    Best_V0_Ratio_List_rPB(Cut_Range_Start_Index_tmp_Index) = V0_Ratio_List(Best_V0_Ratio_Index);
    
    %{
    RMS_ddthRib_Cut_Off_Freq = zeros(size(Cut_Off_Freq_tmp));
    
    for ii = 1:size(Cut_Off_Freq_tmp, 1)
        [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
        ddthRib_filt_tmp = filtfilt(Butter_b, Butter_a, ddthRib);
        
        RMS_ddthRib_Cut_Off_Freq(ii) = (mean((ddthRib_filt_tmp - ddthRib).^2))^(1/2);
    end
    
    Cut_Off_Freq_Line_Range = find(abs(Cut_Off_Freq_tmp - 3) < 1e-2):find(abs(Cut_Off_Freq_tmp - 4) < 1e-2);
    % Cut_Off_Freq_Line_Range = find(Cut_Off_Freq_tmp == 5):find(Cut_Off_Freq_tmp == 6);
    Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddthRib_Cut_Off_Freq(Cut_Off_Freq_Line_Range), 1);
    
    % Cut_Off_Freq = interp1(flip(RMS_ddthRib(1:Cut_Off_Freq_Line_Range(1))), flip(Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1))), Cut_Off_Freq_Line(2))
    Cut_Off_Freq = interp1(RMS_ddthRib_Cut_Off_Freq(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));
    
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
    ddthRib_filt = filtfilt(Butter_b, Butter_a, ddthRib);
    
    RMS_thRib_V0 = zeros(size(V0_Ratio_List));
    
    for V0_Ratio_Index = 1:size(V0_Ratio_List,1)
        ddthRib_spline = spline(Data_Time(Cut_Range), ddthRib_filt(Cut_Range));
%         dthRib_spline = fnint(ddthRib_spline, 0);
        dthRib_spline = fnint(ddthRib_spline, mean(dthRib(Cut_Range(1:2))) * V0_Ratio_List(V0_Ratio_Index));
        thRib_spline = fnint(dthRib_spline, mean(thRib(Cut_Range(1:2))));
        
        RMS_thRib_V0(V0_Ratio_Index) = mean((fnval(thRib_spline, Data_Time(Cut_Range)) - thRib(Cut_Range)).^2)^(1/2);
        
    end
    
    [RMS_thRib_tmp, Best_V0_Ratio_Index] = min(RMS_thRib_V0);
    Best_V0_Ratio_List_thRib(Cut_Range_Start_Index_tmp_Index) = V0_Ratio_List(Best_V0_Ratio_Index);
    %}
    
    RMS_th_level(Cut_Range_Start_Index_tmp_Index)...
        = RMS_thHand_tmp ...
        + RMS_thShoulder_tmp ...
        + RMS_thWaist_tmp ...
        + RMS_hand_x_tmp ...
        + RMS_hand_y_tmp ...
        + RMS_rPB_tmp;
    
%     ddthHand_filt = fnval(ddthHand_spline, Data_Time(Cut_Range));
%     dthHand_filt = fnval(dthHand_spline, Data_Time(Cut_Range));
%     thHand_filt = fnval(thHand_spline, Data_Time(Cut_Range));
%     ddthShoulder_filt = fnval(ddthShoulder_spline, Data_Time(Cut_Range));
%     dthShoulder_filt = fnval(dthShoulder_spline, Data_Time(Cut_Range));
%     thShoulder_filt = fnval(thShoulder_spline, Data_Time(Cut_Range));
%     ddthWaist_filt = fnval(ddthWaist_spline, Data_Time(Cut_Range));
%     dthWaist_filt = fnval(dthWaist_spline, Data_Time(Cut_Range));
%     thWaist_filt = fnval(thWaist_spline, Data_Time(Cut_Range));
%     ddhand_x_filt = fnval(ddhand_x_spline, Data_Time(Cut_Range));
%     dhand_x_filt = fnval(dhand_x_spline, Data_Time(Cut_Range));
%     hand_x_filt = fnval(hand_x_spline, Data_Time(Cut_Range));
%     ddhand_y_filt = fnval(ddhand_y_spline, Data_Time(Cut_Range));
%     dhand_y_filt = fnval(dhand_y_spline, Data_Time(Cut_Range));
%     hand_y_filt = fnval(hand_y_spline, Data_Time(Cut_Range));
%     
%     rPB_filt = fnval(rPB_spline, Data_Time(Cut_Range));
%     drPB_filt = fnval(drPB_spline, Data_Time(Cut_Range));
%     ddrPB_filt = fnval(ddrPB_spline, Data_Time(Cut_Range));
    
%     RMS_th_level(Cut_Range_Start_Index_tmp_Index)...
%         = mean((fnval(thHand_spline, Data_Time(Cut_Range)) - thHand(Cut_Range)).^2)^(1/2) ...
%         + mean((fnval(thShoulder_spline, Data_Time(Cut_Range)) - thShoulder(Cut_Range)).^2)^(1/2) ...
%         + mean((fnval(thWaist_spline, Data_Time(Cut_Range)) - thWaist(Cut_Range)).^2)^(1/2) ...
%         + mean((fnval(hand_x_spline, Data_Time(Cut_Range)) - hand_x(Cut_Range)).^2)^(1/2) ...
%         + mean((fnval(hand_y_spline, Data_Time(Cut_Range)) - hand_y(Cut_Range)).^2)^(1/2);
    
    if mod(Cut_Range_Start_Index_tmp_Index, 50) == 0
        fprintf(strcat("Cut_Range_Start_Index_tmp_Index = ", num2str(Cut_Range_Start_Index_tmp_Index), '\n'))
        toc
        tic
    end
end
toc

[min_value, min_index] = min(RMS_th_level)
Data_Time(Cut_Range_Start_Index_tmp(min_index))

Best_V0_Ratio_thShoulder = Best_V0_Ratio_List_thShoulder(min_index)
Best_V0_Ratio_thWaist = Best_V0_Ratio_List_thWaist(min_index)
Best_V0_Ratio_hand_x = Best_V0_Ratio_List_hand_x(min_index)
Best_V0_Ratio_hand_y = Best_V0_Ratio_List_hand_y(min_index)
Best_V0_Ratio_rPB = Best_V0_Ratio_List_rPB(min_index)

figure(1)
plot(Data_Time(Cut_Range_Start_Index_tmp), RMS_th_level)

%}
%}

%{
RMS_dthHand = zeros(size(Cut_Off_Freq_tmp));

for ii = 1:size(Cut_Off_Freq_tmp, 1)
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
    dthHand_filt_tmp = filtfilt(Butter_b, Butter_a, dthHand);
    
    RMS_dthHand(ii) = (mean((dthHand_filt_tmp - dthHand).^2))^(1/2);
end

figure(1)
plot(Cut_Off_Freq_tmp, RMS_dthHand, 'o-')
ylim([0, 1e-3])

Cut_Off_Freq_Line_Range = 30:70;
Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_dthHand(Cut_Off_Freq_Line_Range), 1)
Cut_Off_Freq_Line_y = polyval(Cut_Off_Freq_Line, Cut_Off_Freq_tmp);
hold on
plot(Cut_Off_Freq_tmp, Cut_Off_Freq_Line_y)
hold off

Cut_Off_Freq = interp1(RMS_dthHand(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2))

[Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
dthHand_filt = filtfilt(Butter_b, Butter_a, dthHand);

figure(2)
plot(Data_Time(1:end-1), dthHand, 'Displayname', 'raw', 'LineStyle', ':', 'LineWidth', 3)
hold on
plot(Data_Time(1:end-1), dthHand_filt, 'Displayname', 'filtfilt')
hold off
legend

dthHand_pp = spline(Data_Time(1:end-1), dthHand);
ddthHand = fnval(fnder(dthHand_pp, 1), Data_Time(1:end-1));

figure(3)
plot(Data_Time(1:end-1), ddthHand)
%}

%{
RMS_thHand = zeros(size(Cut_Off_Freq_tmp));

for ii = 1:size(Cut_Off_Freq_tmp, 1)
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
    thHand_filt_tmp = filtfilt(Butter_b, Butter_a, thHand);
    
    RMS_thHand(ii) = (mean((thHand_filt_tmp - thHand).^2))^(1/2);
end

figure(1)
plot(Cut_Off_Freq_tmp, RMS_thHand, 'o-')
ylim([0, 1e-3])

Cut_Off_Freq_Line_Range = 30:70;
Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_thHand(Cut_Off_Freq_Line_Range), 1)
Cut_Off_Freq_Line_y = polyval(Cut_Off_Freq_Line, Cut_Off_Freq_tmp);
hold on
plot(Cut_Off_Freq_tmp, Cut_Off_Freq_Line_y)
hold off

Cut_Off_Freq = interp1(RMS_thHand(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2))

[Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
thHand_filt = filtfilt(Butter_b, Butter_a, thHand);

figure(2)
plot(Data_Time, thHand, 'Displayname', 'raw', 'LineStyle', ':', 'LineWidth', 3)
hold on
plot(Data_Time, thHand_filt, 'Displayname', 'filtfilt')
hold off
legend

Down_Sampling_Index = 1:5:size(Data_Time,1);

thHand_pp = spline(Data_Time(Down_Sampling_Index), thHand(Down_Sampling_Index));
dthHand = fnval(fnder(thHand_pp, 1), Data_Time);
ddthHand = fnval(fnder(thHand_pp, 2), Data_Time);


figure(3)
plot(Data_Time, dthHand)

figure(4)
plot(Data_Time, ddthHand)

%{
[dthHand_pp, ddthHand_pp] = GetSmoothingDerivates_pp(Data_Time, dthHand);

figure(3)
hold on
plot(Data_Time, fnval(dthHand_pp, Data_Time))
hold off

figure(4)
hold on
plot(Data_Time, fnval(ddthHand_pp, Data_Time))
hold off
%}

%}

%{
ddthHand_filt = fnval(ddthHand_spline, Data_Time(Cut_Range));
dthHand_filt = fnval(dthHand_spline, Data_Time(Cut_Range));
thHand_filt = fnval(thHand_spline, Data_Time(Cut_Range));
ddthShoulder_filt = fnval(ddthShoulder_spline, Data_Time(Cut_Range));
dthShoulder_filt = fnval(dthShoulder_spline, Data_Time(Cut_Range));
thShoulder_filt = fnval(thShoulder_spline, Data_Time(Cut_Range));
ddthWaist_filt = fnval(ddthWaist_spline, Data_Time(Cut_Range));
dthWaist_filt = fnval(dthWaist_spline, Data_Time(Cut_Range));
thWaist_filt = fnval(thWaist_spline, Data_Time(Cut_Range));

rArm = mean(vecnorm(wr_sh, 2, 2));
rBody = mean(vecnorm(sh_tro, 2, 2));
rLeg = mean(vecnorm(tro_an, 2, 2));

pHand = hand(Cut_Range, :);

pShoulder = pHand + rArm * [cos(thHand_filt+1/2*pi), sin(thHand_filt+1/2*pi)];
pWaist = pShoulder + rBody * [cos(thHand_filt+1/2*pi + thShoulder_filt), sin(thHand_filt+1/2*pi + thShoulder_filt)];
pToe = pWaist + rLeg * [cos(thHand_filt+1/2*pi + thShoulder_filt + thWaist_filt), sin(thHand_filt+1/2*pi + thShoulder_filt + thWaist_filt)];


Body(:,:,1) = [pHand(:,1), pShoulder(:,1), pWaist(:,1), pToe(:,1)];
Body(:,:,2) = [pHand(:,2), pShoulder(:,2), pWaist(:,2), pToe(:,2)];

thDataArray = [thHand_filt, thShoulder_filt, thWaist_filt];
MDataArray = [dthHand_filt, dthShoulder_filt, dthWaist_filt];
otherDataArray = [ddthHand_filt, ddthShoulder_filt, ddthWaist_filt];

sliderTimeh = 0.1;

Anime_Cortex(Data_Time(Cut_Range) - Data_Time(Cut_Range(1)), zeros([size(Data_Time(Cut_Range), 1), 1, 2]), zeros([size(Data_Time(Cut_Range), 1), 1, 2]), zeros([size(Data_Time(Cut_Range), 1), 1, 2]), Body, thDataArray, MDataArray, otherDataArray, sliderTimeh)
%}



































