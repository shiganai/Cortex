classdef Rib_Para_Matthew
    %HIP_PARA このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties(Constant)
        % 論文データは (Lumbar)
        % 角度の正の方向は反る方. 
        % Flex は反る. 
        % Ext は屈身する方. 
        
        % シミュレーションするから
        % Ext は屈身
        % Flex は反る方に取る.
        
        tau_M_0_Ext = 688 % 正
        theta_A_0_Ext = -89
        theta_A_W_Ext = 66 % 正
        omega_Max_Ext = 458 % 正
        t_V_E_Max_Ext = 1.10
        t_V_C_Half_Ext = 0.15
        theta_PE_0_Ext = -18 % 負
        theta_PE_1_Ext = -78 % 負
        
        tau_M_0_Flex = -212 % 負 文献値と正負逆
        theta_A_0_Flex = 0 % 文献値と正負逆
        theta_A_W_Flex = 360 % 正
        omega_Max_Flex = -1102 % 負 文献値と正負逆
        t_V_E_Max_Flex = 1.10
        t_V_C_Half_Flex = 0.15
        theta_PE_0_Flex = 0 % 正
        theta_PE_1_Flex = 45 % 正
    end
end

