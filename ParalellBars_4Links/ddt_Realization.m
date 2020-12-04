function dotq = ddt_Realization(t,q,constants, ddpp_TH)
%DDT この関数の概要をここに記述
%   詳細説明をここに記述

g = constants.g;
kPB = constants.kPB;
cPB = constants.cPB;
mAll = constants.mAll;
mPB = constants.mPB;
mArm = constants.mArm;
mUBody = constants.mUBody;
mLBody = constants.mLBody;
mLeg = constants.mLeg;
rArm = constants.rArm;
rUBody = constants.rUBody;
rLBody = constants.rLBody;
rLeg = constants.rLeg;
InertiaModel = constants.InertiaModel;
rArmMCD = constants.rArmMCD;
rUBodyMCD = constants.rUBodyMCD;
rLBodyMCD = constants.rLBodyMCD;
rLegMCD = constants.rLegMCD;
InertiaArm = constants.InertiaArm;
InertiaUBody = constants.InertiaUBody;
InertiaLBody = constants.InertiaLBody;
InertiaLeg = constants.InertiaLeg;
InertiaG = constants.InertiaG;
Hand_Para = constants.Hand_Para;
Shoulder_Para = constants.Shoulder_Para;
Rib_Para = constants.Rib_Para;
Waist_Para = constants.Waist_Para;

rPB = q(1);
thHand = q(2);
thShoulder = q(3);
thRib = q(4);
thWaist = q(5);

drPB = q(6);
dthHand = q(7);
dthShoulder = q(8);
dthRib = q(9);
dthWaist = q(10);

ddthHand_pp = ddpp_TH(1,1);
ddthShoulder_pp = ddpp_TH(2,1);
ddthRib_pp = ddpp_TH(3,1);
ddthWaist_pp = ddpp_TH(4,1);

ddthHand = fnval(ddthHand_pp, t);
ddthShoulder = fnval(ddthShoulder_pp, t);
ddthRib = fnval(ddthRib_pp, t);
ddthWaist = fnval(ddthWaist_pp, t);

FrPB = 2 * (- kPB * rPB - cPB * drPB);

ddrPB = find_ddrPB_only(FrPB,ddthRib,ddthHand,ddthWaist,ddthShoulder,dthRib,dthHand,dthWaist,dthShoulder,g,mArm,mLBody,mLeg,mPB,mUBody,rArm,rArmMCD,rLBody,rLBodyMCD,rLegMCD,rUBody,rUBodyMCD,thHand,thShoulder,thWaist,thRib);
% ddrPB = find_ddrPB_only(MrPB,ddthHand,ddthWaist,ddthShoulder,dthHand,dthWaist,dthShoulder,g,kPB,mArm,mBody,mLeg,mPB,rArm,rArmMCD,rBody,rBodyMCD,rLegMCD,rPB,thHand,thShoulder,thWaist);
% ddrPB = find_ddrPB_only_mex(MrPB,ddthHand,ddthWaist,ddthShoulder,dthHand,dthWaist,dthShoulder,g,kPB,mArm,mBody,mLeg,mPB,rArm,rArmMCD,rBody,rBodyMCD,rLegMCD,rPB,thHand,thShoulder,thWaist);

dotq = [drPB, dthHand, dthShoulder, dthRib, dthWaist, ddrPB, ddthHand, ddthShoulder, ddthRib, ddthWaist]';
end















































