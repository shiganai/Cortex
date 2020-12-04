classdef constants_tmp
    %CONSTANTS このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    %{
    properties(Constant)
        
        constants_Name = 'constants_Yeadon'
        g = 9.8
        kPB = 1.9831 * 1e4
        cPB = 4.1
        mAll = 80
        mPB = 2
        
        mArm = 10.83
        mBody = 35.47
        mLeg = 13.91 + 7.59
        
        rArm = 0.537;
        rBody = 0.569;
        rLeg = 0.8727;% = 0.374*7/3;
        
%         InertiaModel = InertiaModel_CompleteSticks;
        InertiaModel = struct(InertiaModel_Yeadon);
        
        rArmMCD = constants_HH.rArm / constants_HH.InertiaModel.rArm * constants_HH.InertiaModel.rArmMCD;
        rBodyMCD = constants_HH.rBody / constants_HH.InertiaModel.rBody * constants_HH.InertiaModel.rBodyMCD;
        rLegMCD = constants_HH.rLeg / constants_HH.InertiaModel.rLeg * constants_HH.InertiaModel.rLegMCD;
        
        InertiaArm = constants_HH.mArm / constants_HH.InertiaModel.mArm * (constants_HH.rArm / constants_HH.InertiaModel.rArm)^2 * constants_HH.InertiaModel.InertiaArm;
        InertiaBody = constants_HH.mBody / constants_HH.InertiaModel.mBody * (constants_HH.rBody / constants_HH.InertiaModel.rBody)^2 * constants_HH.InertiaModel.InertiaBody;
        InertiaLeg = constants_HH.mLeg / constants_HH.InertiaModel.mLeg * (constants_HH.rLeg / constants_HH.InertiaModel.rLeg)^2 * constants_HH.InertiaModel.InertiaLeg;
        
        InertiaG =  find_InertiaG(constants_HH.InertiaLeg,constants_HH.InertiaArm,constants_HH.InertiaBody,...
            constants_HH.mArm,constants_HH.mBody,constants_HH.mLeg,...
            constants_HH.rArm,constants_HH.rArmMCD,constants_HH.rBody,constants_HH.rBodyMCD,constants_HH.rLegMCD,...
            pi, 0)
        
        Hand_Para = struct(Hand_Para_Matthew);
        Shoulder_Para = struct(Shoulder_Para_Matthew);
        Waist_Para = struct(Waist_Para_Matthew);
    end
    %}
    
    properties
        constants_Name
        g
        kPB
        cPB
        mAll
        mPB
        
        mArm
        mUBody
        mLBody
        mLeg
        
        rArm
        rUBody
        rLBody
        rLeg
        
        %         InertiaModel = InertiaModel_CompleteSticks;
        InertiaModel
        
        rArmMCD
        rUBodyMCD
        rLBodyMCD
        rLegMCD
        
        InertiaArm
        InertiaUBody
        InertiaLBody
        InertiaLeg
        
        InertiaG
        
        Hand_Para
        Shoulder_Para
        Rib_Para
        Waist_Para
    end
    
    methods
        function obj = constants_tmp(mAll, r_Toe_Ankle, r_Ankle_Knee, r_Knee_Waist, r_Waist_Rib, r_Rib_Shoulder, r_Shoulder_Elbow, r_Elbow_Wrist, r_Wrist_Finger, r_Shoulder_Ear, r_Ear_Top)
            
            obj.constants_Name = 'constants_tmp';
            obj.g = 9.8;
            obj.kPB = 1.9831 * 1e4;
            obj.cPB = 4.1;
            obj.mAll = mAll;
            obj.mPB = 2;
%             obj.rAll = 1.815;
            
            obj.InertiaModel = [];
            
%             r_Ankle_Toe = 22/100;
%             r_Knee_Ankle = 42/100;
%             r_Waist_Knee = 44/100;
%             r_Shoulder_Waist = (30 + 27)/100;
%             r_Elbow_Shoulder = 32/100;
%             r_Wrist_Elbow = 28/100;
%             r_Finger_Wrist = 9.5/100;
%             r_Ear_Shoulder = (31 - 15)/100;
%             r_Top_Ear = 15/100;

            [obj.mArm, obj.mUBody, obj.mLBody, obj.mLeg, obj.rArm, obj.rUBody, obj.rLBody, obj.rLeg, obj.rArmMCD, obj.rUBodyMCD, obj.rLBodyMCD, obj.rLegMCD, obj.InertiaArm, obj.InertiaUBody, obj.InertiaLBody, obj.InertiaLeg]...
                    = Calc_Parameter_AE_4Links(mAll, r_Toe_Ankle, r_Ankle_Knee, r_Knee_Waist, r_Waist_Rib, r_Rib_Shoulder, r_Shoulder_Elbow, r_Elbow_Wrist, r_Wrist_Finger, r_Shoulder_Ear, r_Ear_Top);
            
%             obj.InertiaG =  find_InertiaG(obj.InertiaLeg,obj.InertiaArm,obj.InertiaBody,...
%                 obj.mArm,obj.mBody,obj.mLeg,...
%                 obj.rArm,obj.rArmMCD,obj.rBody,obj.rBodyMCD,obj.rLegMCD,...
%                 pi, 0);
            
            obj.InertiaG = find_InertiaG(obj.InertiaLeg,obj.InertiaArm,obj.InertiaLBody,obj.InertiaUBody,...
                obj.mArm,obj.mLBody,obj.mLeg,obj.mUBody,obj.rArm,obj.rArmMCD,obj.rLBody,obj.rLBodyMCD,obj.rLegMCD,obj.rUBody,obj.rUBodyMCD,...
                0 , pi, 0, 0, 0, 0);
            
            obj.Hand_Para = struct(Hand_Para_Matthew);
            obj.Shoulder_Para = struct(Shoulder_Para_Matthew);
            obj.Rib_Para = struct(Rib_Para_Matthew);
            obj.Waist_Para = struct(Waist_Para_Matthew);
        end
    end
end

