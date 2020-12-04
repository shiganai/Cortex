function vG = find_vG(dthRib,dthHand,dthWaist,dthShoulder,dxHand,dyHand,mArm,mLBody,mLeg,mUBody,rArm,rArmMCD,rLBody,rLBodyMCD,rLegMCD,rUBody,rUBodyMCD,thHand,thShoulder,thWaist,thRib)
%FIND_VG
%    VG = FIND_VG(T,DTHXIP,DTHHAND,DTHWAIST,DTHSHOULDER,DXHAND,DYHAND,MARM,MLBODY,MLEG,MUBODY,RARM,RARMMCD,RLBODY,RLBODYMCD,RLEGMCD,RUBODY,RUBODYMCD,THHAND,THSHOULDER,THWAIST,THXIP)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    13-Nov-2020 17:09:31

t2 = dthHand+dthShoulder;
t5 = -dxHand;
t6 = mArm+mLBody+mLeg+mUBody;
t7 = pi./2.0;
t3 = dthRib+t2;
t4 = dthWaist+t2;
t8 = t7+thHand;
t9 = 1.0./t6;
t10 = cos(t8);
t11 = sin(t8);
t12 = t8+thShoulder;
t13 = cos(t12);
t14 = sin(t12);
t15 = t12+thWaist;
t16 = t12+thRib;
t19 = dthHand.*rArm.*t10;
t20 = dthHand.*rArm.*t11;
t17 = cos(t16);
t18 = sin(t16);
t21 = rUBody.*t2.*t13;
t22 = rUBody.*t2.*t14;
vG = [-t9.*(-mArm.*(dxHand-dthHand.*rArmMCD.*t11)+mLBody.*(t5+t20+t22+rLBodyMCD.*t3.*t18)+mLeg.*(t5+t20+t22+rLegMCD.*t4.*sin(t15)+rLBody.*t3.*t18)+mUBody.*(t5+t20+rUBodyMCD.*t2.*t14)),t9.*(mLBody.*(dyHand+t19+t21+rLBodyMCD.*t3.*t17)+mArm.*(dyHand+dthHand.*rArmMCD.*t10)+mLeg.*(dyHand+t19+t21+rLegMCD.*t4.*cos(t15)+rLBody.*t3.*t17)+mUBody.*(dyHand+t19+rUBodyMCD.*t2.*t13))];
