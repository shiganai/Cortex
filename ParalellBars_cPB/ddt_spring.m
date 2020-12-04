function dotq = ddt_spring(t, q, kPB, cPB, mPB)
%DDT_SPRIG この関数の概要をここに記述
%   詳細説明をここに記述
rPB = q(1);
drPB = q(2);

ddrPB = 2 * (-kPB * rPB - cPB * drPB) / mPB;

dotq = [drPB, ddrPB]';
end

