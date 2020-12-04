function [handRI, handRO, wrRI, wrRO, elbRI, elbRO, shRF, shRU, shRB, handLI, handLO, wrLI, wrLO, elbLI, elbLO, shLF, shLU, shLB, toeR, baRI, baRO, heelR, anRI, anRO, knRI, knRO, troR, toeL, baLI, baLO, heelL, anLI, anLO, knLI, knLO, troL, head, earR, earL, clav, c7, ribR, ribL, xiph, t12, ASISR, ASISL, PSISR, PSISL, ThR, ThL]...
    = Get_2D_Coordinate_51P(FileName, Axe1, Axe2)
%GET_2D_COORDINATE この関数の概要をここに記述
%   詳細説明をここに記述

Body_Data_All = (dlmread(FileName,'',6,2))/1000;%単位変換[mm]→[m]

% handRI = Get_2D_Coordinate(Body_Data_All(:,1:3), Axe1, Axe2);

OutputName = ["handRI", "handRO", "wrRI", "wrRO", "elbRI", "elbRO", "shRF", "shRU", "shRB", "handLI", "handLO", "wrLI", "wrLO", "elbLI", "elbLO", "shLF", "shLU", "shLB", "toeR", "baRI", "baRO", "heelR", "anRI", "anRO", "knRI", "knRO", "troR", "toeL", "baLI", "baLO", "heelL", "anLI", "anLO", "knLI", "knLO", "troL", "head", "earR", "earL", "clav", "c7", "ribR", "ribL", "xiph", "t12", "ASISR", "ASISL", "PSISR", "PSISL", "ThR", "ThL"]';

for ii = 1:51
    eval(strcat(...
        OutputName(ii,1), " = Get_2D_Coordinate(Body_Data_All(:,", num2str(3*ii-2), ":", num2str(3*ii), "), Axe1, Axe2);" ...
        ));
end
end

