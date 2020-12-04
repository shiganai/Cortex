function [mArm, mUBody, mLBody, mLeg, rArm, rUBody, rLBody, rLeg, rArmMCD, rUBodyMCD, rLBodyMCD, rLegMCD, InertiaArm, InertiaUBody, InertiaLBody, InertiaLeg]...
    = Calc_Parameter_AE_4Links(mAll, r_Toe_Ankle, r_Ankle_Knee, r_Knee_Waist, r_Waist_Rib, r_Rib_Shoulder, r_Shoulder_Elbow, r_Elbow_Wrist, r_Wrist_Finger, r_Shoulder_Ear, r_Ear_Top)
%CALC_PARAMETER_AE この関数の概要をここに記述
%   詳細説明をここに記述

h_Toe = 0;
h_Ankle = h_Toe + r_Toe_Ankle;
h_Knee = h_Ankle + r_Ankle_Knee;
h_Waist = h_Knee + r_Knee_Waist;
h_Rib = h_Waist + r_Waist_Rib;
h_Shoulder = h_Rib + r_Rib_Shoulder;
h_Elbow = h_Shoulder + r_Shoulder_Elbow;
h_Wrist = h_Elbow + r_Elbow_Wrist;
h_Finger = h_Wrist;
h_Ear = h_Shoulder + r_Shoulder_Ear;
h_Top = h_Ear + r_Ear_Top;

% ScatVec = [h_Toe, h_Ankle, h_Knee, h_Waist, h_Shoulder, h_Elbow, h_Wrist, h_Finger];
% scatter(zeros(size(ScatVec)), ScatVec)


m_Toe_Ankle = mAll *  2 * 1.1/100;
m_Ankle_Knee = mAll *  2 * 5.1/100;
m_Knee_Waist = mAll *  2 * 11/100;
m_Waist_Rib= mAll *  18.7/100;
m_Rib_Shoulder = mAll *  30.2/100;
m_Shoulder_Elbow = mAll *  2 * 2.7/100;
m_Elbow_Wrist = mAll *  2 * 1.6/100;
m_Wrist_Finger = mAll * 2 * 0.6/100;
m_Shoulder_Top = mAll *  6.9/100;

hG_Toe_Ankle = h_Toe + (h_Ankle - h_Toe) * 59.5/100;
hG_Ankle_Knee = h_Knee + (h_Ankle - h_Knee) * 40.6/100;
hG_Knee_Waist = h_Waist + (h_Knee - h_Waist) * 47.5/100;
hG_Waist_Rib = h_Rib + (h_Waist - h_Rib) * 60.9/100;
hG_Rib_Shoulder = h_Shoulder + (h_Rib - h_Shoulder) * 42.8/100;
hG_Shoulder_Elbow = h_Shoulder + (h_Elbow - h_Shoulder) * 52.9/100;
hG_Elbow_Wrist = h_Elbow + (h_Wrist - h_Elbow) * 41.5/100;
hG_Wrist_Finger = h_Wrist + (h_Finger - h_Wrist) * 41.5/100;
hG_Shoulder_Top = h_Top + (h_Ear - h_Top) * 82.1/100;

% hold on
% ScatVec = [hG_Toe_Ankle, hG_Ankle_Knee, hG_Knee_Waist, hG_Waist_Shoulder, hG_Shoulder_Elbow, hG_Elbow_Wrist, hG_Wrist_Finger];
% scatter(zeros(size(ScatVec)), ScatVec)
% hold off

Inertia_Toe_Ankle = m_Toe_Ankle * (r_Toe_Ankle * 20.4/100)^2;
Inertia_Ankle_Knee = m_Ankle_Knee * (r_Ankle_Knee * 27.4/100)^2;
Inertia_Knee_Waist = m_Knee_Waist * (r_Knee_Waist * 27.8/100)^2;
Inertia_Waist_Rib = m_Waist_Rib * (r_Waist_Rib * 42.5/100)^2;
Inertia_Rib_Shoulder = m_Rib_Shoulder * (r_Rib_Shoulder * 35.0/100)^2;
Inertia_Shoulder_Elbow = m_Shoulder_Elbow * (r_Shoulder_Elbow * 26.2/100)^2;
Inertia_Elbow_Wrist = m_Elbow_Wrist * (r_Elbow_Wrist * 27.9/100)^2;
Inertia_Wrist_Finger = m_Wrist_Finger * (r_Wrist_Finger * 51.9/100)^2;
Inertia_Shoulder_Top = m_Shoulder_Top * (r_Ear_Top * 47.9/100)^2;


mArm = m_Shoulder_Elbow + m_Elbow_Wrist + m_Wrist_Finger;
mUBody = m_Rib_Shoulder + m_Shoulder_Top;
mLBody = m_Waist_Rib;
mLeg = m_Toe_Ankle + m_Ankle_Knee + m_Knee_Waist;

rArm = h_Finger - h_Shoulder;
rUBody = h_Shoulder - h_Rib;
rLBody = h_Rib - h_Waist;
rLeg = h_Waist - h_Toe;

hG_Arm = (m_Shoulder_Elbow * hG_Shoulder_Elbow + m_Elbow_Wrist * hG_Elbow_Wrist + m_Wrist_Finger * hG_Wrist_Finger) / mArm;
hG_UBody = (m_Rib_Shoulder * hG_Rib_Shoulder + m_Shoulder_Top * hG_Shoulder_Top) / mUBody;
hG_LBody = (m_Waist_Rib * hG_Waist_Rib) / mLBody;
hG_Leg = (m_Toe_Ankle * hG_Toe_Ankle + m_Ankle_Knee * hG_Ankle_Knee + m_Knee_Waist * hG_Knee_Waist) / mLeg;

rArmMCD = hG_Arm - h_Shoulder;
rUBodyMCD = h_Shoulder - hG_UBody;
rLBodyMCD = h_Rib - hG_LBody;
rLegMCD = h_Waist - hG_Leg;

if rArmMCD <= 0
    error('rArmMCD は正であるはずです')
end
if rUBodyMCD <= 0
    error('rUBodyMCD は正であるはずです')
end
if rLBodyMCD <= 0
    error('rLBodyMCD は正であるはずです')
end
if rLegMCD <= 0
    error('rLegMCD は正であるはずです')
end

InertiaArm = Inertia_Shoulder_Elbow + Inertia_Elbow_Wrist + Inertia_Wrist_Finger...
    + m_Shoulder_Elbow * (hG_Shoulder_Elbow - hG_Arm)^2 ...
    + m_Elbow_Wrist * (hG_Elbow_Wrist - hG_Arm)^2 ...
    + m_Wrist_Finger * (hG_Wrist_Finger - hG_Arm)^2 ...
    ;
InertiaUBody = Inertia_Rib_Shoulder + Inertia_Shoulder_Top...
    + m_Rib_Shoulder * (hG_Rib_Shoulder - hG_UBody)^2 ...
    + m_Shoulder_Top * (hG_Shoulder_Top - hG_UBody)^2 ...
    ;
InertiaLBody = Inertia_Waist_Rib...
    + m_Waist_Rib * (hG_Waist_Rib - hG_LBody)^2 ...
    ;
InertiaLeg = Inertia_Toe_Ankle + Inertia_Ankle_Knee + Inertia_Knee_Waist...
    + m_Toe_Ankle * (hG_Toe_Ankle - hG_Leg)^2 ...
    + m_Ankle_Knee * (hG_Ankle_Knee - hG_Leg)^2 ...
    + m_Knee_Waist * (hG_Knee_Waist - hG_Leg)^2 ...
    ;



end

