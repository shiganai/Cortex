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
thHand = GetThetaFromXY(wr, sh) - 1/2*pi;
thShoulder = GetThetaFromXY(sh, tro) - thHand - 1/2*pi;
thWaist = GetThetaFromXY(tro, an) - thHand - 1/2*pi - thShoulder;

thShoulder = thShoulder + 2*pi;
thWaist = thWaist - 2*pi;

[thHand_pp, dthHand_pp, ddthHand_pp] = GetSmoothingDerivates_pp(Data_Time, thHand);
[thShoulder_pp, dthShoulder_pp, ddthShoulder_pp] = GetSmoothingDerivates_pp(Data_Time, thShoulder);
[thWaist_pp, dthWaist_pp, ddthWaist_pp] = GetSmoothingDerivates_pp(Data_Time, thWaist);

thHand = fnval(thHand_pp, Data_Time);
dthHand = fnval(dthHand_pp, Data_Time);
ddthHand = fnval(ddthHand_pp, Data_Time);

thShoulder = fnval(thShoulder_pp, Data_Time);
dthShoulder = fnval(dthShoulder_pp, Data_Time);
ddthShoulder = fnval(ddthShoulder_pp, Data_Time);

thWaist = fnval(thWaist_pp, Data_Time);
dthWaist = fnval(dthWaist_pp, Data_Time);
ddthWaist = fnval(ddthWaist_pp, Data_Time);

thDataArray = [thHand, thShoulder, thWaist];
MDataArray = [dthHand, dthShoulder, dthWaist];
otherDataArray = [ddthHand, ddthShoulder, ddthWaist];
sliderTimeh = 0.1;

Anime_Cortex(Data_Time, BarRs, BarLs, BodyR, BodyL, thDataArray, MDataArray, otherDataArray, sliderTimeh)

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

figure(2)
plot(Data_Time, rPB)

drPB = diff(rPB) / diff(Data_Time(1:2));
ddrPB = diff(drPB) / diff(Data_Time(1:2));

figure(3)
plot(Data_Time(1:end-1), drPB);

figure(4)
plot(Data_Time(1:end-2), ddrPB);

%{
figure(3)
plot(Data_Time, [vecnorm(wr_sh, 2, 2), vecnorm(sh_tro, 2, 2), vecnorm(tro_an, 2, 2)])

rArm = mean(vecnorm(wr_sh, 2, 2));
rBody = mean(vecnorm(sh_tro, 2, 2));
rLeg = mean(vecnorm(tro_an, 2, 2));

pHand = hand;

pShoulder = pHand + rArm * [cos(thHand+1/2*pi), sin(thHand+1/2*pi)];
pWaist = pShoulder + rBody * [cos(thHand+1/2*pi + thShoulder), sin(thHand+1/2*pi + thShoulder)];
pToe = pWaist + rLeg * [cos(thHand+1/2*pi + thShoulder + thWaist), sin(thHand+1/2*pi + thShoulder + thWaist)];


Body(:,:,1) = [pHand(:,1), pShoulder(:,1), pWaist(:,1), pToe(:,1)];
Body(:,:,2) = [pHand(:,2), pShoulder(:,2), pWaist(:,2), pToe(:,2)];

Anime_Cortex(Data_Time, zeros([size(Data_Time, 1), 1, 2]), zeros([size(Data_Time, 1), 1, 2]), zeros([size(Data_Time, 1), 1, 2]), Body, thDataArray, MDataArray, otherDataArray, sliderTimeh)
%}



































