function TauRib = find_TauRib(InertiaLeg,InertiaLBody,ddrPB,ddthRib,ddthHand,ddthWaist,ddthShoulder,dthHand,dthWaist,dthShoulder,g,mLBody,mLeg,rArm,rLBody,rLBodyMCD,rLegMCD,rUBody,thHand,thShoulder,thWaist,thRib)
%FIND_TAUXIP
%    TAUXIP = FIND_TAUXIP(INERTIALEG,INERTIALBODY,DDRPB,DDTHXIP,DDTHHAND,DDTHWAIST,DDTHSHOULDER,DTHHAND,DTHWAIST,DTHSHOULDER,G,MLBODY,MLEG,RARM,RLBODY,RLBODYMCD,RLEGMCD,RUBODY,THHAND,THSHOULDER,THWAIST,THXIP)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    13-Nov-2020 17:13:19

t2 = dthHand.^2;
t3 = dthWaist.^2;
t4 = dthShoulder.^2;
t5 = rLBody.^2;
t6 = rLBodyMCD.^2;
t7 = pi./2.0;
t8 = t7+thHand;
t9 = cos(t8);
t10 = sin(t8);
t11 = t8+thShoulder;
t12 = cos(t11);
t13 = sin(t11);
t14 = t11+thWaist;
t15 = t11+thRib;
t16 = cos(t14);
t17 = cos(t15);
t18 = sin(t14);
t19 = sin(t15);
t20 = t17.^2;
t21 = t19.^2;
TauRib = InertiaLeg.*ddthRib+InertiaLeg.*ddthHand+InertiaLBody.*ddthRib+InertiaLeg.*ddthWaist+InertiaLBody.*ddthHand+InertiaLeg.*ddthShoulder+InertiaLBody.*ddthShoulder+ddrPB.*mLBody.*rLBodyMCD.*t17+ddrPB.*mLeg.*rLBody.*t17+ddthRib.*mLBody.*t6.*t20+ddthRib.*mLBody.*t6.*t21+ddthHand.*mLBody.*t6.*t20+ddthHand.*mLBody.*t6.*t21+ddthShoulder.*mLBody.*t6.*t20+ddthShoulder.*mLBody.*t6.*t21+ddthRib.*mLeg.*t5.*t20+ddthRib.*mLeg.*t5.*t21+ddthHand.*mLeg.*t5.*t20+ddthHand.*mLeg.*t5.*t21+ddthShoulder.*mLeg.*t5.*t20+ddthShoulder.*mLeg.*t5.*t21+g.*mLBody.*rLBodyMCD.*t17+g.*mLeg.*rLBody.*t17-mLBody.*rArm.*rLBodyMCD.*t2.*t10.*t17+mLBody.*rArm.*rLBodyMCD.*t2.*t9.*t19-mLeg.*rArm.*rLBody.*t2.*t10.*t17+mLeg.*rArm.*rLBody.*t2.*t9.*t19+mLeg.*rLBody.*rLegMCD.*t2.*t16.*t19-mLeg.*rLBody.*rLegMCD.*t2.*t17.*t18+mLeg.*rLBody.*rLegMCD.*t3.*t16.*t19-mLeg.*rLBody.*rLegMCD.*t3.*t17.*t18+mLeg.*rLBody.*rLegMCD.*t4.*t16.*t19-mLeg.*rLBody.*rLegMCD.*t4.*t17.*t18-mLBody.*rLBodyMCD.*rUBody.*t2.*t13.*t17+mLBody.*rLBodyMCD.*rUBody.*t2.*t12.*t19-mLBody.*rLBodyMCD.*rUBody.*t4.*t13.*t17+mLBody.*rLBodyMCD.*rUBody.*t4.*t12.*t19-mLeg.*rLBody.*rUBody.*t2.*t13.*t17+mLeg.*rLBody.*rUBody.*t2.*t12.*t19-mLeg.*rLBody.*rUBody.*t4.*t13.*t17+mLeg.*rLBody.*rUBody.*t4.*t12.*t19+ddthHand.*mLBody.*rArm.*rLBodyMCD.*t9.*t17+ddthHand.*mLBody.*rArm.*rLBodyMCD.*t10.*t19+ddthHand.*mLeg.*rArm.*rLBody.*t9.*t17+ddthHand.*mLeg.*rArm.*rLBody.*t10.*t19+ddthHand.*mLeg.*rLBody.*rLegMCD.*t16.*t17+ddthHand.*mLeg.*rLBody.*rLegMCD.*t18.*t19+ddthWaist.*mLeg.*rLBody.*rLegMCD.*t16.*t17+ddthWaist.*mLeg.*rLBody.*rLegMCD.*t18.*t19+ddthShoulder.*mLeg.*rLBody.*rLegMCD.*t16.*t17+ddthShoulder.*mLeg.*rLBody.*rLegMCD.*t18.*t19+ddthHand.*mLBody.*rLBodyMCD.*rUBody.*t12.*t17+ddthHand.*mLBody.*rLBodyMCD.*rUBody.*t13.*t19+ddthShoulder.*mLBody.*rLBodyMCD.*rUBody.*t12.*t17+ddthShoulder.*mLBody.*rLBodyMCD.*rUBody.*t13.*t19+ddthHand.*mLeg.*rLBody.*rUBody.*t12.*t17+ddthHand.*mLeg.*rLBody.*rUBody.*t13.*t19+ddthShoulder.*mLeg.*rLBody.*rUBody.*t12.*t17+ddthShoulder.*mLeg.*rLBody.*rUBody.*t13.*t19+dthHand.*dthWaist.*mLeg.*rLBody.*rLegMCD.*t16.*t19.*2.0-dthHand.*dthWaist.*mLeg.*rLBody.*rLegMCD.*t17.*t18.*2.0+dthHand.*dthShoulder.*mLeg.*rLBody.*rLegMCD.*t16.*t19.*2.0-dthHand.*dthShoulder.*mLeg.*rLBody.*rLegMCD.*t17.*t18.*2.0+dthWaist.*dthShoulder.*mLeg.*rLBody.*rLegMCD.*t16.*t19.*2.0-dthWaist.*dthShoulder.*mLeg.*rLBody.*rLegMCD.*t17.*t18.*2.0-dthHand.*dthShoulder.*mLBody.*rLBodyMCD.*rUBody.*t13.*t17.*2.0+dthHand.*dthShoulder.*mLBody.*rLBodyMCD.*rUBody.*t12.*t19.*2.0-dthHand.*dthShoulder.*mLeg.*rLBody.*rUBody.*t13.*t17.*2.0+dthHand.*dthShoulder.*mLeg.*rLBody.*rUBody.*t12.*t19.*2.0;
