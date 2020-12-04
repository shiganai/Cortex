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

hand = (handR + handL)/2;
wr = (wrR + wrL)/2;
sh = (shR + shL)/2;
tro = (troR + troL)/2;
an = (anR + anL)/2;

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

Cut_Off_Freq_tmp = (0.5:0.5:30)';
Nichest_Freq = 200 / 2;

Cut_Range_Start_Index_tmp = (find(Data_Time == 3.5):find(Data_Time == 6.5))';
RMS_th_level = zeros(size(Cut_Range_Start_Index_tmp));

tic
for Cut_Range_Start_Index_tmp_Index = 1:size(Cut_Range_Start_Index_tmp, 1)
    
    
%     Cut_Range = Cut_Range_Start_Index_tmp(Cut_Range_Start_Index_tmp_Index):find(Data_Time == 7.95);
    Cut_Range = Cut_Range_Start_Index_tmp(Cut_Range_Start_Index_tmp_Index):find(Data_Time == 8.5);
    
    %{/
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
    
    ddthHand_spline = spline(Data_Time(Cut_Range), ddthHand_filt(Cut_Range));
    dthHand_spline = fnint(ddthHand_spline, 0);
    % dthHand_spline = fnint(ddthHand_spline, mean(dthHand(Cut_Range(1:2))));
    thHand_spline = fnint(dthHand_spline, mean(thHand(Cut_Range(1:2))));
    
    
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
    
    ddthShoulder_spline = spline(Data_Time(Cut_Range), ddthShoulder_filt(Cut_Range));
    dthShoulder_spline = fnint(ddthShoulder_spline, mean(dthShoulder(Cut_Range(1:2))));
    thShoulder_spline = fnint(dthShoulder_spline, mean(thShoulder(Cut_Range(1:2))));
    
    
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
    
    ddthWaist_spline = spline(Data_Time(Cut_Range), ddthWaist_filt(Cut_Range));
    dthWaist_spline = fnint(ddthWaist_spline, mean(dthWaist(Cut_Range(1:2))));
    thWaist_spline = fnint(dthWaist_spline, mean(thWaist(Cut_Range(1:2))));
    
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
    
    ddhand_x_spline = spline(Data_Time(Cut_Range), ddhand_x_filt(Cut_Range));
%     dhand_x_spline = fnint(ddhand_x_spline, 0);
    dhand_x_spline = fnint(ddhand_x_spline, mean(dhand_x(Cut_Range(1:2))));
    hand_x_spline = fnint(dhand_x_spline, mean(hand_x(Cut_Range(1:2))));
    
    RMS_ddhand_y = zeros(size(Cut_Off_Freq_tmp));
    
    for ii = 1:size(Cut_Off_Freq_tmp, 1)
        [Butter_b, Butter_a] = butter(2, Cut_Off_Freq_tmp(ii) / Nichest_Freq);
        ddhand_y_filt_tmp = filtfilt(Butter_b, Butter_a, ddhand_y);
        
        RMS_ddhand_y(ii) = (mean((ddhand_y_filt_tmp - ddhand_y).^2))^(1/2);
    end
    
    Cut_Off_Freq_Line_Range = find(Cut_Off_Freq_tmp == 10):find(Cut_Off_Freq_tmp == 30);
    Cut_Off_Freq_Line = polyfit(Cut_Off_Freq_tmp(Cut_Off_Freq_Line_Range), RMS_ddhand_y(Cut_Off_Freq_Line_Range), 1);
    
    Cut_Off_Freq = interp1(RMS_ddhand_y(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_tmp(1:Cut_Off_Freq_Line_Range(1)), Cut_Off_Freq_Line(2));
    
    [Butter_b, Butter_a] = butter(2, Cut_Off_Freq/ Nichest_Freq);
    ddhand_y_filt = filtfilt(Butter_b, Butter_a, ddhand_y);
    
    ddhand_y_spline = spline(Data_Time(Cut_Range), ddhand_y_filt(Cut_Range));
    dhand_y_spline = fnint(ddhand_y_spline, mean(dhand_y(Cut_Range(1:2))));
    hand_y_spline = fnint(dhand_y_spline, mean(hand_y(Cut_Range(1:2))));
    
    RMS_th_level(Cut_Range_Start_Index_tmp_Index)...
        = mean((fnval(thHand_spline, Data_Time(Cut_Range)) - thHand(Cut_Range)).^2)^(1/2) ...
        + mean((fnval(thShoulder_spline, Data_Time(Cut_Range)) - thShoulder(Cut_Range)).^2)^(1/2) ...
        + mean((fnval(thWaist_spline, Data_Time(Cut_Range)) - thWaist(Cut_Range)).^2)^(1/2) ...
        + mean((fnval(hand_x_spline, Data_Time(Cut_Range)) - hand_x(Cut_Range)).^2)^(1/2) ...
        + mean((fnval(hand_y_spline, Data_Time(Cut_Range)) - hand_y(Cut_Range)).^2)^(1/2);
    
    
    if mod(Cut_Range_Start_Index_tmp_Index, 50) == 0
        fprintf(strcat("Cut_Range_Start_Index_tmp_Index = ", num2str(Cut_Range_Start_Index_tmp_Index), '\n'))
        toc
        tic
    end
end
toc

[min_value, min_index] = min(RMS_th_level)
Data_Time(Cut_Range_Start_Index_tmp(min_index))


figure(1)
plot(Data_Time(Cut_Range_Start_Index_tmp), RMS_th_level)

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



































