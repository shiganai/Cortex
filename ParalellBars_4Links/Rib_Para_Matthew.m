classdef Rib_Para_Matthew
    %HIP_PARA ���̃N���X�̊T�v�������ɋL�q
    %   �ڍא����������ɋL�q
    
    properties(Constant)
        % �_���f�[�^�� (Lumbar)
        % �p�x�̐��̕����͔����. 
        % Flex �͔���. 
        % Ext �͋��g�����. 
        
        % �V�~�����[�V�������邩��
        % Ext �͋��g
        % Flex �͔�����Ɏ��.
        
        tau_M_0_Ext = 688 % ��
        theta_A_0_Ext = -89
        theta_A_W_Ext = 66 % ��
        omega_Max_Ext = 458 % ��
        t_V_E_Max_Ext = 1.10
        t_V_C_Half_Ext = 0.15
        theta_PE_0_Ext = -18 % ��
        theta_PE_1_Ext = -78 % ��
        
        tau_M_0_Flex = -212 % �� �����l�Ɛ����t
        theta_A_0_Flex = 0 % �����l�Ɛ����t
        theta_A_W_Flex = 360 % ��
        omega_Max_Flex = -1102 % �� �����l�Ɛ����t
        t_V_E_Max_Flex = 1.10
        t_V_C_Half_Flex = 0.15
        theta_PE_0_Flex = 0 % ��
        theta_PE_1_Flex = 45 % ��
    end
end

