function dotq = ddt_spring(t, q, kPB, cPB, mPB)
%DDT_SPRIG ���̊֐��̊T�v�������ɋL�q
%   �ڍא����������ɋL�q
rPB = q(1);
drPB = q(2);

ddrPB = 2 * (-kPB * rPB - cPB * drPB) / mPB;

dotq = [drPB, ddrPB]';
end

