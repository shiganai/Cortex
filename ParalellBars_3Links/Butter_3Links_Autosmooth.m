clear all

Bar_Data_All = (dlmread('Trimmed_Trimmed_ShSs3-Bar.trc','',6,1))/1000;%単位変換[mm]→[m]
Body_Data_All = (dlmread('Trimmed_Trimmed_ShSs3-51P.trc','',6,2))/1000;%単位変換[mm]→[m]

Data_Time = Bar_Data_All(:,1) * 1000;

% rdata = [Bar_Data_All(:,2:end), Body_Data_All];
nFr = size(Data_Time, 1);
time_int = diff(Data_Time(1:2)); % 1/200
g = 'n';

np = size(Bar_Data_All(:,2:end),2)/3; % 3 * 2 * 2 + 51;
[Bar_Data_All(:,2:end), ocf] = autosmooth(Bar_Data_All(:,2:end),np,nFr,time_int,g);

np = size(Body_Data_All,2)/3; % 3 * 2 * 2 + 51;
[Body_Data_All, ~] = autosmooth(Body_Data_All,np,nFr,time_int,g);

Cut_Range = find(abs(Data_Time - 6.5650) < 1/400):find(abs(Data_Time - 7.965) < 1/400);
Bar_Data_All = Bar_Data_All(Cut_Range, :); % 1598 番目が離手
Body_Data_All = Body_Data_All(Cut_Range, :);

Data_Time = Data_Time(Cut_Range, :);


%----------------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------------

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


%----------------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------------

[handRI, handRO, wrRI, wrRO, elbRI, elbRO, shRF, shRU, shRB, handLI, handLO, wrLI, wrLO, elbLI, elbLO, shLF, shLU, shLB, toeR, baRI, baRO, heelR, anRI, anRO, knRI, knRO, troR, toeL, baLI, baLO, heelL, anLI, anLO, knLI, knLO, troL, head, earR, earL, clav, c7, ribR, ribL, xiph, t12, ASISR, ASISL, PSISR, PSISL, ThR, ThL]...
    = Get_2D_Coordinate_51P_Autosmooth(Body_Data_All, Axe_X_UnitVec, Axe_Y_UnitVec);

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

% とりあえずの左右プロット
%{/
thDataArray = zeros([size(Data_Time, 1), 1]);
MDataArray = zeros([size(Data_Time, 1), 1]);
otherDataArray = zeros([size(Data_Time, 1), 1]);
% otherDataArray = [MthHand, MthShoulder, MthWaist];

sliderTimeh = 0.1;


Anime_Cortex = Anime_Cortex(Data_Time, ...
    BarRs, ...
    BarLs, ...
    BodyR, ...
    BodyL, ...
    thDataArray, MDataArray, otherDataArray, sliderTimeh);
Anime_Cortex.UIFigure.Name = 'Anime_Cortex.UIFigure.Name';
ylim(Anime_Cortex.axThData,[-2, 2])
%}

wr_sh = sh - wr;
sh_tro = tro - sh;
tro_an = an - tro;

% rPB = hand(:,2) - hand(1,2);
hand_x = hand(:,1);
hand_y = hand(:,2);
thHand = GetThetaFromXY(wr, sh) - 1/2*pi;
thShoulder = GetThetaFromXY(sh, tro) - thHand - 1/2*pi;
thWaist = GetThetaFromXY(tro, an) - thHand - 1/2*pi - thShoulder;

thShoulder = thShoulder + 0*pi;
thWaist = thWaist - 0*pi;

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

figure_Number = 0;

figure_Number = figure_Number + 1;
figure(figure_Number)
plot(Data_Time, rPB);
xlabel('時間')
ylabel('rPB\_raw')

figure_Number = figure_Number + 1;
figure(figure_Number)
plot(Data_Time(1:end-2), ddrPB);
xlabel('時間')
ylabel('ddrPB\_raw')

figure_Number = figure_Number + 1;
figure(figure_Number)
plot(Data_Time, thHand);
xlabel('時間')
ylabel('thHand')

figure_Number = figure_Number + 1;
figure(figure_Number)
plot(Data_Time(1:end-2), ddthHand);
xlabel('時間')
ylabel('ddthHand')

figure_Number = figure_Number + 1;
figure(figure_Number)
plot(Data_Time, thShoulder);
xlabel('時間')
ylabel('thShoulder')

figure_Number = figure_Number + 1;
figure(figure_Number)
plot(Data_Time(1:end-2), ddthShoulder);
xlabel('時間')
ylabel('ddthShoulder')

figure_Number = figure_Number + 1;
figure(figure_Number)
plot(Data_Time, thWaist);
xlabel('時間')
ylabel('thWaist')

figure_Number = figure_Number + 1;
figure(figure_Number)
plot(Data_Time(1:end-2), ddthWaist);
xlabel('時間')
ylabel('ddthWaist')

dhand_x = fnval(fnder(spline(Data_Time, hand_x), 1), Data_Time);
ddhand_x = fnval(fnder(spline(Data_Time, hand_x), 2), Data_Time);

dhand_y = fnval(fnder(spline(Data_Time, hand_y), 1), Data_Time);
ddhand_y = fnval(fnder(spline(Data_Time, hand_y), 2), Data_Time);

drPB = fnval(fnder(spline(Data_Time, rPB), 1), Data_Time);
ddrPB = fnval(fnder(spline(Data_Time, rPB), 2), Data_Time);

dthHand = fnval(fnder(spline(Data_Time, thHand), 1), Data_Time);
ddthHand = fnval(fnder(spline(Data_Time, thHand), 2), Data_Time);

dthShoulder = fnval(fnder(spline(Data_Time, thShoulder), 1), Data_Time);
ddthShoulder = fnval(fnder(spline(Data_Time, thShoulder), 2), Data_Time);

dthWaist = fnval(fnder(spline(Data_Time, thWaist), 1), Data_Time);
ddthWaist = fnval(fnder(spline(Data_Time, thWaist), 2), Data_Time);

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

[mArm, mBody, mLeg, rArm, rBody, rLeg, rArmMCD, rBodyMCD, rLegMCD, InertiaArm, InertiaBody, InertiaLeg] = ...
    Calc_Parameter_AE(mAll, r_Ankle_Toe, r_Knee_Ankle, r_Waist_Knee, r_Shoulder_Waist, r_Elbow_Shoulder, r_Wrist_Elbow, r_Finger_Wrist, r_Ear_Shoulder, r_Top_Ear);

xHand = hand_x;
yHand = hand_y;

pHand = [hand_x, hand_y];

pShoulder = pHand + rArm * [cos(thHand+1/2*pi), sin(thHand+1/2*pi)];
pWaist = pShoulder + rBody * [cos(thHand+1/2*pi + thShoulder), sin(thHand+1/2*pi + thShoulder)];
pToe = pWaist + rLeg * [cos(thHand+1/2*pi + thShoulder + thWaist), sin(thHand+1/2*pi + thShoulder + thWaist)];

vHand = [dhand_x, dhand_y];
vShoulder = vHand + rArm * [-sin(thHand+1/2*pi), cos(thHand+1/2*pi)] .* dthHand;
vWaist = vShoulder + rBody * [-sin(thHand+1/2*pi + thShoulder), cos(thHand+1/2*pi + thShoulder)] .* (dthHand + dthShoulder);
vToe = vWaist + rLeg * [-sin(thHand+1/2*pi + thShoulder + thWaist), cos(thHand+1/2*pi + thShoulder + thWaist)] .* (dthHand + dthShoulder + dthWaist);

pArmMCD = pHand + rArmMCD * [cos(thHand+1/2*pi), sin(thHand+1/2*pi)];
pBodyMCD = pShoulder + rBodyMCD * [cos(thHand+1/2*pi + thShoulder), sin(thHand+1/2*pi + thShoulder)];
pLegMCD = pWaist + rLegMCD * [cos(thHand+1/2*pi + thShoulder + thWaist), sin(thHand+1/2*pi + thShoulder + thWaist)];

vArmMCD = vHand + rArmMCD * [-sin(thHand+1/2*pi), cos(thHand+1/2*pi)] .* dthHand;
vBodyMCD = vShoulder + rBodyMCD * [-sin(thHand+1/2*pi + thShoulder), cos(thHand+1/2*pi + thShoulder)] .* (dthHand + dthShoulder);
vLegMCD = vWaist + rLegMCD * [-sin(thHand+1/2*pi + thShoulder + thWaist), cos(thHand+1/2*pi + thShoulder + thWaist)] .* (dthHand + dthShoulder + dthWaist);

pG = find_pG(mArm,mBody,mLeg,rArm,rArmMCD,rBody,rBodyMCD,rLegMCD,thHand,thShoulder,thWaist,xHand,yHand);
pG_Anime(:,:,1) = pG(:,1);
pG_Anime(:,:,2) = pG(:,2);


Body(:,:,1) = [pHand(:,1), pShoulder(:,1), pWaist(:,1), pToe(:,1)];
Body(:,:,2) = [pHand(:,2), pShoulder(:,2), pWaist(:,2), pToe(:,2)];

kPB = 2.4160e+04;
cPB = 8.29;
mPB = 4.84 * 2;

% kPB = 1.9831 * 1e4;
% cPB = 4.1;
% mPB = 2;

g = 9.8;

MrPB = find_MrPB(ddrPB,ddthHand,ddthWaist,ddthShoulder,dthHand,dthWaist,dthShoulder,g,mArm,mBody,mLeg,mPB,rArm,rArmMCD,rBody,rBodyMCD,rLegMCD,thHand,thShoulder,thWaist);
MthHand = find_MthHand(InertiaLeg,InertiaArm,InertiaBody,ddrPB,ddthHand,ddthWaist,ddthShoulder,dthHand,dthWaist,dthShoulder,g,mArm,mBody,mLeg,rArm,rArmMCD,rBody,rBodyMCD,rLegMCD,thHand,thShoulder,thWaist);
MthShoulder = find_MthShoulder(InertiaLeg,InertiaBody,ddrPB,ddthHand,ddthWaist,ddthShoulder,dthHand,dthWaist,dthShoulder,g,mBody,mLeg,rArm,rBody,rBodyMCD,rLegMCD,thHand,thShoulder,thWaist);
MthWaist = find_MthWaist(InertiaLeg,ddrPB,ddthHand,ddthWaist,ddthShoulder,dthHand,dthShoulder,g,mLeg,rArm,rBody,rLegMCD,thHand,thShoulder,thWaist);


Hand_Para = Hand_Para_Matthew;
Shoulder_Para = Shoulder_Para_Matthew;
Waist_Para = Waist_Para_Matthew;

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

ActivatiedRate_Limit = 3;

MthHand_ActivatedRate(MthHand_ActivatedRate > ActivatiedRate_Limit) = ActivatiedRate_Limit;
MthShoulder_ActivatedRate(MthShoulder_ActivatedRate > ActivatiedRate_Limit) = ActivatiedRate_Limit;
MthWaist_ActivatedRate(MthWaist_ActivatedRate > ActivatiedRate_Limit) = ActivatiedRate_Limit;

MthHand_ActivatedRate(MthHand_ActivatedRate < -ActivatiedRate_Limit) = -ActivatiedRate_Limit;
MthShoulder_ActivatedRate(MthShoulder_ActivatedRate < -ActivatiedRate_Limit) = -ActivatiedRate_Limit;
MthWaist_ActivatedRate(MthWaist_ActivatedRate < -ActivatiedRate_Limit) = -ActivatiedRate_Limit;

base_y_rPB = (0.002857 + -0.007649)/2;
base_y_Position = mean(Bar0(:,2));

T = 1/2 * mPB * vecnorm(vHand, 2, 2).^2 ...
    + 1/2 * mArm * vecnorm(vArmMCD, 2, 2).^2 ...
    + 1/2 * mBody * vecnorm(vBodyMCD, 2, 2).^2 ...
    + 1/2 * mLeg * vecnorm(vLegMCD, 2, 2).^2 ...
    + 1/2 * InertiaArm * dthHand.^2 ...
    + 1/2 * InertiaBody * (dthHand + dthShoulder).^2 ...
    + 1/2 * InertiaLeg * (dthHand + dthShoulder + dthWaist).^2;

U = 1/2 * (2*kPB) * (pHand(:,2) - base_y_Position - base_y_rPB).^2 ...
    + mPB * g * (pHand(:,2) - base_y_Position) ...
    + mArm * g * (pArmMCD(:,2) - base_y_Position) ...
    + mBody * g * (pBodyMCD(:,2) - base_y_Position) ...
    + mLeg * g * (pLegMCD(:,2) - base_y_Position);

figure_Number = figure_Number + 1;
figure(figure_Number)
plot(Data_Time, MthHand, 'DisplayName', '手首のトルク')
hold on
plot(Data_Time, MthShoulder, 'DisplayName', '肩のトルク')
plot(Data_Time, MthWaist, 'DisplayName', '腰のトルク')
hold off
legend
xlabel('時間')


figure_Number = figure_Number + 1;
figure(figure_Number)
plot(Data_Time, MthHand_ActivatedRate, 'DisplayName', '手首')
hold on
plot(Data_Time, MthShoulder_ActivatedRate, 'DisplayName', '肩')
plot(Data_Time, MthWaist_ActivatedRate, 'DisplayName', '腰')
hold off
legend
xlabel('時間')

figure_Number = figure_Number + 1;
figure(figure_Number)
plot(Data_Time, [T, U, T + U])
legend('T', 'U', 'T + U')
grid on
xlabel('時間')

%{
thDataArray = [MthHand_ActivatedRate, MthShoulder_ActivatedRate, MthWaist_ActivatedRate];
MDataArray = [MrPB];
otherDataArray = thHand;
% otherDataArray = [MthHand, MthShoulder, MthWaist];

sliderTimeh = 0.1;


Anime_Cortex = Anime_Cortex(Data_Time - Data_Time(Cut_Range(1)), ...
    zeros([size(Data_Time, 1), 1, 2]), ...
    zeros([size(Data_Time, 1), 1, 2]), ...
    Body, ...
    pG_Anime, ...
    thDataArray, MDataArray, otherDataArray, sliderTimeh);
Anime_Cortex.UIFigure.Name = 'Anime_Cortex.UIFigure.Name';
ylim(Anime_Cortex.axThData,[-2, 2])
%}





