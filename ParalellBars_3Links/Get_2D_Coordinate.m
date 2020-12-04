function Coordinate_2D = Get_2D_Coordinate(Coordinate_3D, Axe1, Axe2)
%GET_2D_COORDINATE ���̊֐��̊T�v�������ɋL�q
%   �ڍא����������ɋL�q

Coordinate_2D = [Coordinate_3D(:,1) * Axe1(1) + Coordinate_3D(:,2) * Axe1(2) + Coordinate_3D(:,3) * Axe1(3),...
    Coordinate_3D(:,1) * Axe2(1) + Coordinate_3D(:,2) * Axe2(2) + Coordinate_3D(:,3) * Axe2(3)];
end

