function InertiaG = find_InertiaG(InertiaLeg,InertiaArm,InertiaLBody,InertiaUBody,mArm,mLBody,mLeg,mUBody,rArm,rArmMCD,rLBody,rLBodyMCD,rLegMCD,rUBody,rUBodyMCD,thHand,thShoulder,thWaist,thRib,xHand,yHand)
%FIND_INERTIAG
%    INERTIAG = FIND_INERTIAG(T,INERTIALEG,INERTIAARM,INERTIALBODY,INERTIAUBODY,MARM,MLBODY,MLEG,MUBODY,RARM,RARMMCD,RLBODY,RLBODYMCD,RLEGMCD,RUBODY,RUBODYMCD,THHAND,THSHOULDER,THWAIST,THXIP,XHAND,YHAND)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    13-Nov-2020 17:09:37

t2 = mArm+mLBody+mLeg+mUBody;
t3 = pi./2.0;
t4 = t3+thHand;
t5 = 1.0./t2;
t6 = cos(t4);
t7 = sin(t4);
t8 = t4+thShoulder;
t9 = cos(t8);
t10 = sin(t8);
t11 = t8+thWaist;
t12 = t8+thRib;
t13 = rArm.*t6;
t14 = rArmMCD.*t6;
t15 = rArm.*t7;
t16 = rArmMCD.*t7;
t17 = cos(t11);
t18 = cos(t12);
t19 = sin(t11);
t20 = sin(t12);
t21 = rUBody.*t9;
t22 = rUBodyMCD.*t9;
t23 = rUBody.*t10;
t24 = rUBodyMCD.*t10;
t25 = t14+xHand;
t26 = t16+yHand;
t27 = rLBody.*t18;
t28 = rLBodyMCD.*t18;
t29 = rLegMCD.*t17;
t30 = rLBody.*t20;
t31 = rLBodyMCD.*t20;
t32 = rLegMCD.*t19;
t33 = mArm.*t25;
t34 = mArm.*t26;
t35 = t13+t22+xHand;
t36 = t15+t24+yHand;
t37 = mUBody.*t35;
t38 = mUBody.*t36;
t39 = t13+t21+t28+xHand;
t40 = t15+t23+t31+yHand;
t43 = t13+t21+t27+t29+xHand;
t44 = t15+t23+t30+t32+yHand;
t41 = mLBody.*t39;
t42 = mLBody.*t40;
t45 = mLeg.*t44;
t46 = mLeg.*t43;
t47 = t33+t37+t41+t46;
t48 = t34+t38+t42+t45;
t49 = t5.*t47;
t50 = t5.*t48;
t51 = -t49;
t52 = -t50;
InertiaG = InertiaLeg+InertiaArm+InertiaLBody+InertiaUBody+mArm.*((t25+t51).^2+(t26+t52).^2)+mLBody.*((t39+t51).^2+(t40+t52).^2)+mLeg.*((t43+t51).^2+(t44+t52).^2)+mUBody.*((t35+t51).^2+(t36+t52).^2);