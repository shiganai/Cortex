function ddthShoulder = find_ddthShoulder(InertiaLeg,InertiaArm,InertiaBody,MrPB,MthHand,MthWaist,MthShoulder,dthHand,dthWaist,dthShoulder,g,mArm,mBody,mLeg,mPB,rArm,rArmMCD,rBody,rBodyMCD,rLegMCD,thHand,thShoulder,thWaist)
%FIND_DDTHSHOULDER
%    DDTHSHOULDER = FIND_DDTHSHOULDER(INERTIALEG,INERTIAARM,INERTIABODY,MRPB,MTHHAND,MTHWAIST,MTHSHOULDER,DTHHAND,DTHWAIST,DTHSHOULDER,G,MARM,MBODY,MLEG,MPB,RARM,RARMMCD,RBODY,RBODYMCD,RLEGMCD,THHAND,THSHOULDER,THWAIST)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    19-Jul-2020 21:38:05

t2 = dthHand.^2;
t3 = dthWaist.^2;
t4 = dthShoulder.^2;
t5 = mArm.^2;
t6 = mBody.^2;
t7 = mLeg.^2;
t8 = rArm.^2;
t9 = rArmMCD.^2;
t10 = rBody.^2;
t11 = rBodyMCD.^2;
t12 = rLegMCD.^2;
t13 = InertiaLeg.*InertiaArm.*mArm;
t14 = InertiaLeg.*InertiaBody.*mArm;
t15 = InertiaLeg.*InertiaArm.*mBody;
t16 = InertiaLeg.*InertiaBody.*mBody;
t17 = InertiaLeg.*InertiaArm.*mLeg;
t18 = InertiaLeg.*InertiaBody.*mLeg;
t19 = InertiaLeg.*InertiaArm.*mPB;
t20 = InertiaLeg.*InertiaBody.*mPB;
t25 = pi./2.0;
t21 = InertiaBody.*t13;
t22 = InertiaBody.*t15;
t23 = InertiaBody.*t17;
t24 = InertiaBody.*t19;
t26 = t25+thHand;
t27 = cos(t26);
t28 = sin(t26);
t29 = t26+thShoulder;
t30 = cos(t29);
t31 = sin(t29);
t32 = t29+thWaist;
t33 = t27.^2;
t34 = t28.^2;
t35 = cos(t32);
t36 = sin(t32);
t37 = t30.^2;
t38 = t31.^2;
t41 = g.*mBody.*rBodyMCD.*t30;
t42 = g.*mLeg.*rBody.*t30;
t44 = InertiaLeg.*mArm.*mBody.*t8.*t33;
t45 = InertiaLeg.*mArm.*mBody.*t9.*t33;
t46 = InertiaLeg.*mArm.*mLeg.*t8.*t33;
t47 = InertiaLeg.*mArm.*mLeg.*t9.*t33;
t48 = InertiaLeg.*mArm.*mPB.*t9.*t33;
t49 = InertiaLeg.*mBody.*mPB.*t8.*t33;
t50 = InertiaLeg.*mLeg.*mPB.*t8.*t33;
t51 = InertiaLeg.*mArm.*mBody.*t8.*t34;
t52 = InertiaLeg.*mArm.*mBody.*t9.*t34;
t53 = InertiaLeg.*mArm.*mLeg.*t8.*t34;
t54 = InertiaLeg.*mArm.*mLeg.*t9.*t34;
t55 = InertiaLeg.*mArm.*mPB.*t9.*t34;
t56 = InertiaLeg.*mBody.*mPB.*t8.*t34;
t57 = InertiaLeg.*mLeg.*mPB.*t8.*t34;
t58 = InertiaLeg.*mArm.*mBody.*rArm.*rArmMCD.*t33.*2.0;
t59 = InertiaLeg.*mArm.*mLeg.*rArm.*rArmMCD.*t33.*2.0;
t73 = InertiaLeg.*t5.*t9.*t34;
t74 = InertiaLeg.*t6.*t8.*t34;
t75 = InertiaLeg.*t7.*t8.*t34;
t76 = InertiaLeg.*mBody.*mLeg.*t8.*t34.*2.0;
t79 = mBody.*t8.*t14.*t33;
t80 = mBody.*t9.*t14.*t33;
t81 = mLeg.*t8.*t14.*t33;
t82 = mLeg.*t9.*t14.*t33;
t83 = mPB.*t9.*t14.*t33;
t84 = mPB.*t8.*t16.*t33;
t85 = mPB.*t8.*t18.*t33;
t86 = mBody.*t8.*t14.*t34;
t87 = mBody.*t9.*t14.*t34;
t88 = mLeg.*t8.*t14.*t34;
t89 = mLeg.*t9.*t14.*t34;
t90 = mPB.*t9.*t14.*t34;
t91 = mPB.*t8.*t16.*t34;
t92 = mPB.*t8.*t18.*t34;
t93 = mBody.*rArm.*rArmMCD.*t14.*t33.*2.0;
t94 = mLeg.*rArm.*rArmMCD.*t14.*t33.*2.0;
t126 = mLeg.*t8.*t16.*t34.*2.0;
t141 = InertiaLeg.*mArm.*mBody.*rArm.*rBodyMCD.*t27.*t30;
t142 = InertiaLeg.*mArm.*mBody.*rArmMCD.*rBodyMCD.*t27.*t30;
t143 = InertiaLeg.*mArm.*mLeg.*rArm.*rBody.*t27.*t30;
t144 = InertiaLeg.*mArm.*mLeg.*rArmMCD.*rBody.*t27.*t30;
t145 = InertiaLeg.*mBody.*mPB.*rArm.*rBodyMCD.*t27.*t30;
t146 = InertiaLeg.*mLeg.*mPB.*rArm.*rBody.*t27.*t30;
t147 = InertiaLeg.*mArm.*mBody.*rArm.*rBodyMCD.*t28.*t31;
t148 = InertiaLeg.*mArm.*mLeg.*rArm.*rBody.*t28.*t31;
t149 = InertiaLeg.*mBody.*mLeg.*rArm.*rBody.*t28.*t31;
t150 = InertiaLeg.*mBody.*mLeg.*rArm.*rBodyMCD.*t28.*t31;
t151 = InertiaLeg.*mBody.*mPB.*rArm.*rBodyMCD.*t28.*t31;
t152 = InertiaLeg.*mLeg.*mPB.*rArm.*rBody.*t28.*t31;
t153 = InertiaLeg.*rArm.*rBodyMCD.*t6.*t28.*t31;
t154 = InertiaLeg.*rArm.*rBody.*t7.*t28.*t31;
t388 = InertiaLeg.*mArm.*mBody.*mLeg.*rBody.*rBodyMCD.*t8.*t27.*t28.*t30.*t31.*4.0;
t389 = InertiaLeg.*mBody.*mLeg.*mPB.*rBody.*rBodyMCD.*t8.*t27.*t28.*t30.*t31.*4.0;
t390 = InertiaLeg.*mArm.*rArm.*rArmMCD.*t6.*t11.*t27.*t28.*t30.*t31.*2.0;
t391 = InertiaLeg.*mArm.*rArm.*rArmMCD.*t7.*t10.*t27.*t28.*t30.*t31.*2.0;
t394 = InertiaLeg.*mArm.*t6.*t8.*t11.*t27.*t28.*t30.*t31.*2.0;
t395 = InertiaLeg.*mArm.*t7.*t8.*t10.*t27.*t28.*t30.*t31.*2.0;
t396 = InertiaLeg.*mPB.*t6.*t8.*t11.*t27.*t28.*t30.*t31.*2.0;
t397 = InertiaLeg.*mPB.*t7.*t8.*t10.*t27.*t28.*t30.*t31.*2.0;
t39 = t35.^2;
t40 = t36.^2;
t43 = g.*mLeg.*rLegMCD.*t35;
t60 = InertiaLeg.*mArm.*mBody.*t11.*t37;
t61 = InertiaLeg.*mArm.*mLeg.*t10.*t37;
t62 = InertiaLeg.*mBody.*mLeg.*t10.*t37;
t63 = InertiaLeg.*mBody.*mLeg.*t11.*t37;
t64 = InertiaLeg.*mBody.*mPB.*t11.*t37;
t65 = InertiaLeg.*mLeg.*mPB.*t10.*t37;
t66 = InertiaLeg.*mArm.*mBody.*t11.*t38;
t67 = InertiaLeg.*mArm.*mLeg.*t10.*t38;
t68 = InertiaLeg.*mBody.*mLeg.*t10.*t38;
t69 = InertiaLeg.*mBody.*mLeg.*t11.*t38;
t70 = InertiaLeg.*mBody.*mPB.*t11.*t38;
t71 = InertiaLeg.*mLeg.*mPB.*t10.*t38;
t72 = InertiaLeg.*mBody.*mLeg.*rBody.*rBodyMCD.*t37.*2.0;
t77 = -t58;
t78 = -t59;
t95 = InertiaLeg.*t6.*t11.*t38;
t96 = InertiaLeg.*t7.*t10.*t38;
t98 = mBody.*t11.*t13.*t37;
t99 = mLeg.*t10.*t13.*t37;
t100 = mLeg.*t10.*t15.*t37;
t101 = mLeg.*t11.*t15.*t37;
t102 = mPB.*t11.*t15.*t37;
t103 = mPB.*t10.*t17.*t37;
t104 = mBody.*t11.*t13.*t38;
t105 = mLeg.*t10.*t13.*t38;
t106 = mLeg.*t10.*t15.*t38;
t107 = mLeg.*t11.*t15.*t38;
t108 = mPB.*t11.*t15.*t38;
t109 = mPB.*t10.*t17.*t38;
t110 = mLeg.*rBody.*rBodyMCD.*t15.*t37.*2.0;
t123 = InertiaBody.*t73;
t124 = InertiaBody.*t74;
t125 = InertiaBody.*t75;
t127 = -t93;
t128 = -t94;
t155 = mLeg.*rArm.*rLegMCD.*t2.*t28.*t35;
t156 = mLeg.*rArm.*rLegMCD.*t2.*t27.*t36;
t157 = -t142;
t158 = -t144;
t159 = mLeg.*rBody.*rLegMCD.*t3.*t30.*t36;
t160 = mLeg.*rBody.*rLegMCD.*t3.*t31.*t35;
t161 = dthHand.*dthWaist.*mLeg.*rBody.*rLegMCD.*t30.*t36.*2.0;
t162 = dthHand.*dthWaist.*mLeg.*rBody.*rLegMCD.*t31.*t35.*2.0;
t163 = dthWaist.*dthShoulder.*mLeg.*rBody.*rLegMCD.*t30.*t36.*2.0;
t164 = dthWaist.*dthShoulder.*mLeg.*rBody.*rLegMCD.*t31.*t35.*2.0;
t168 = InertiaLeg.*mArm.*mBody.*mLeg.*rArm.*rArmMCD.*rBody.*rBodyMCD.*t33.*t37.*4.0;
t169 = mLeg.*t10.*t37.*t44;
t170 = mLeg.*t11.*t37.*t44;
t171 = mLeg.*t10.*t37.*t45;
t172 = mLeg.*t11.*t37.*t45;
t173 = mPB.*t11.*t37.*t45;
t174 = mPB.*t10.*t37.*t47;
t175 = mLeg.*t10.*t37.*t49;
t176 = mLeg.*t11.*t37.*t49;
t177 = mLeg.*t10.*t37.*t51;
t178 = mLeg.*t10.*t38.*t44;
t179 = mLeg.*t11.*t37.*t51;
t180 = mLeg.*t11.*t38.*t44;
t181 = mLeg.*t10.*t37.*t52;
t182 = mLeg.*t10.*t38.*t45;
t183 = mLeg.*t11.*t37.*t52;
t184 = mLeg.*t11.*t38.*t45;
t185 = mPB.*t11.*t37.*t52;
t186 = mPB.*t11.*t38.*t45;
t187 = mPB.*t10.*t37.*t54;
t188 = mPB.*t10.*t38.*t47;
t189 = mLeg.*t10.*t37.*t56;
t190 = mLeg.*t10.*t38.*t49;
t191 = mLeg.*t11.*t37.*t56;
t192 = mLeg.*t11.*t38.*t49;
t193 = mLeg.*t10.*t37.*t58;
t194 = mLeg.*t11.*t37.*t58;
t195 = mLeg.*rBody.*rBodyMCD.*t37.*t44.*2.0;
t196 = mLeg.*rBody.*rBodyMCD.*t37.*t45.*2.0;
t197 = mLeg.*rBody.*rBodyMCD.*t37.*t49.*2.0;
t198 = mLeg.*t10.*t38.*t51;
t199 = mLeg.*t11.*t38.*t51;
t200 = mLeg.*t10.*t38.*t52;
t201 = mLeg.*t11.*t38.*t52;
t202 = mPB.*t11.*t38.*t52;
t203 = mPB.*t10.*t38.*t54;
t204 = mLeg.*t10.*t38.*t56;
t205 = mLeg.*t11.*t38.*t56;
t206 = mLeg.*t10.*t38.*t58;
t207 = mLeg.*t11.*t38.*t58;
t208 = mLeg.*rBody.*rBodyMCD.*t37.*t52.*2.0;
t213 = mLeg.*rBody.*rBodyMCD.*t38.*t51.*2.0;
t214 = mLeg.*rBody.*rBodyMCD.*t38.*t56.*2.0;
t255 = mArm.*t11.*t37.*t74;
t258 = mBody.*t11.*t37.*t73;
t259 = mArm.*t10.*t37.*t75;
t262 = mLeg.*t10.*t37.*t73;
t263 = mBody.*t10.*t37.*t75;
t264 = mLeg.*t10.*t37.*t74;
t265 = mBody.*t11.*t37.*t75;
t266 = mLeg.*t11.*t37.*t74;
t267 = mPB.*t11.*t37.*t74;
t269 = mPB.*t10.*t37.*t75;
t272 = mBody.*t11.*t38.*t73;
t274 = mLeg.*t10.*t38.*t73;
t275 = mBody.*t10.*t38.*t75;
t276 = mLeg.*t10.*t38.*t74;
t277 = mBody.*t11.*t38.*t75;
t278 = mLeg.*t11.*t38.*t74;
t286 = mBody.*rBody.*rBodyMCD.*t37.*t75.*2.0;
t287 = mLeg.*rBody.*rBodyMCD.*t37.*t74.*2.0;
t291 = mBody.*rBody.*rBodyMCD.*t38.*t75.*2.0;
t292 = mLeg.*rBody.*rBodyMCD.*t38.*t74.*2.0;
t387 = mLeg.*rArmMCD.*rBody.*t28.*t31.*t141.*4.0;
t392 = -t388;
t393 = -t389;
t398 = mArm.*rArm.*rBody.*t7.*t12.*t28.*t30.*t35.*t36;
t399 = mArm.*rArm.*rBody.*t7.*t12.*t27.*t31.*t35.*t36;
t400 = mArm.*rArmMCD.*rBody.*t7.*t12.*t27.*t31.*t35.*t36;
t401 = mBody.*rArm.*rBody.*t7.*t12.*t28.*t30.*t35.*t36;
t402 = mBody.*rArm.*rBodyMCD.*t7.*t12.*t28.*t30.*t35.*t36;
t403 = mPB.*rArm.*rBody.*t7.*t12.*t28.*t30.*t35.*t36;
t404 = mPB.*rArm.*rBody.*t7.*t12.*t27.*t31.*t35.*t36;
t412 = mArm.*rArm.*rArmMCD.*t7.*t12.*t27.*t28.*t35.*t36.*2.0;
t416 = -t394;
t417 = -t395;
t418 = -t396;
t419 = -t397;
t439 = mArm.*t7.*t8.*t12.*t27.*t28.*t35.*t36.*2.0;
t440 = mPB.*t7.*t8.*t12.*t27.*t28.*t35.*t36.*2.0;
t491 = mBody.*rBody.*rBodyMCD.*t7.*t12.*t30.*t31.*t35.*t36.*2.0;
t509 = InertiaBody.*mArm.*t7.*t8.*t12.*t27.*t28.*t35.*t36.*-2.0;
t510 = InertiaBody.*mPB.*t7.*t8.*t12.*t27.*t28.*t35.*t36.*-2.0;
t511 = mArm.*t7.*t10.*t12.*t30.*t31.*t35.*t36.*2.0;
t512 = mBody.*t7.*t10.*t12.*t30.*t31.*t35.*t36.*2.0;
t513 = mPB.*t7.*t10.*t12.*t30.*t31.*t35.*t36.*2.0;
t522 = InertiaArm.*mArm.*t7.*t10.*t12.*t30.*t31.*t35.*t36.*-2.0;
t523 = InertiaArm.*mBody.*t7.*t10.*t12.*t30.*t31.*t35.*t36.*-2.0;
t524 = InertiaArm.*mPB.*t7.*t10.*t12.*t30.*t31.*t35.*t36.*-2.0;
t537 = mArm.*mBody.*rArm.*rArmMCD.*rBody.*rBodyMCD.*t7.*t12.*t30.*t31.*t33.*t35.*t36.*4.0;
t550 = mArm.*mBody.*rArm.*rArmMCD.*rBody.*rBodyMCD.*t7.*t12.*t27.*t28.*t35.*t36.*t37.*-2.0;
t551 = mArm.*mBody.*rArm.*rArmMCD.*rBody.*rBodyMCD.*t7.*t12.*t27.*t28.*t35.*t36.*t38.*-2.0;
t552 = mArm.*mBody.*rArm.*rArmMCD.*t7.*t10.*t12.*t30.*t31.*t33.*t35.*t36.*4.0;
t576 = rBody.*rBodyMCD.*t6.*t7.*t8.*t12.*t30.*t31.*t34.*t35.*t36.*4.0;
t577 = mArm.*mBody.*t7.*t8.*t11.*t12.*t27.*t28.*t35.*t36.*t37.*-2.0;
t578 = mBody.*mPB.*t7.*t8.*t11.*t12.*t27.*t28.*t35.*t36.*t37.*-2.0;
t579 = mArm.*mBody.*t7.*t8.*t11.*t12.*t27.*t28.*t35.*t36.*t38.*-2.0;
t580 = mBody.*mPB.*t7.*t8.*t11.*t12.*t27.*t28.*t35.*t36.*t38.*-2.0;
t581 = mArm.*mBody.*t7.*t8.*t10.*t12.*t30.*t31.*t33.*t35.*t36.*-2.0;
t582 = mArm.*mBody.*t7.*t9.*t10.*t12.*t30.*t31.*t33.*t35.*t36.*-2.0;
t583 = mArm.*mPB.*t7.*t9.*t10.*t12.*t30.*t31.*t33.*t35.*t36.*-2.0;
t584 = mBody.*mPB.*t7.*t8.*t10.*t12.*t30.*t31.*t33.*t35.*t36.*-2.0;
t585 = t5.*t7.*t9.*t10.*t12.*t30.*t31.*t34.*t35.*t36.*2.0;
t586 = t6.*t7.*t8.*t10.*t12.*t30.*t31.*t34.*t35.*t36.*2.0;
t587 = t6.*t7.*t8.*t11.*t12.*t30.*t31.*t34.*t35.*t36.*2.0;
t588 = mArm.*mBody.*t7.*t8.*t10.*t12.*t30.*t31.*t34.*t35.*t36.*-2.0;
t589 = mArm.*mBody.*t7.*t9.*t10.*t12.*t30.*t31.*t34.*t35.*t36.*-2.0;
t590 = mArm.*mPB.*t7.*t9.*t10.*t12.*t30.*t31.*t34.*t35.*t36.*-2.0;
t591 = mBody.*mPB.*t7.*t8.*t10.*t12.*t30.*t31.*t34.*t35.*t36.*-2.0;
t97 = -t72;
t111 = InertiaArm.*mArm.*mLeg.*t12.*t39;
t112 = InertiaBody.*mArm.*mLeg.*t12.*t39;
t113 = InertiaArm.*mBody.*mLeg.*t12.*t39;
t114 = InertiaBody.*mBody.*mLeg.*t12.*t39;
t115 = InertiaArm.*mLeg.*mPB.*t12.*t39;
t116 = InertiaBody.*mLeg.*mPB.*t12.*t39;
t117 = InertiaArm.*mArm.*mLeg.*t12.*t40;
t118 = InertiaBody.*mArm.*mLeg.*t12.*t40;
t119 = InertiaArm.*mBody.*mLeg.*t12.*t40;
t120 = InertiaBody.*mBody.*mLeg.*t12.*t40;
t121 = InertiaArm.*mLeg.*mPB.*t12.*t40;
t122 = InertiaBody.*mLeg.*mPB.*t12.*t40;
t129 = InertiaArm.*t95;
t130 = InertiaArm.*t96;
t131 = -t110;
t132 = InertiaArm.*t7.*t12.*t40;
t133 = InertiaBody.*t7.*t12.*t40;
t165 = -t161;
t166 = -t163;
t167 = -t159;
t209 = mArm.*mBody.*mLeg.*t8.*t12.*t33.*t39;
t210 = mArm.*mBody.*mLeg.*t9.*t12.*t33.*t39;
t211 = mArm.*mLeg.*mPB.*t9.*t12.*t33.*t39;
t212 = mBody.*mLeg.*mPB.*t8.*t12.*t33.*t39;
t215 = mArm.*mBody.*mLeg.*t8.*t12.*t34.*t39;
t216 = mArm.*mBody.*mLeg.*t8.*t12.*t33.*t40;
t217 = mArm.*mBody.*mLeg.*t9.*t12.*t34.*t39;
t218 = mArm.*mBody.*mLeg.*t9.*t12.*t33.*t40;
t219 = mArm.*mLeg.*mPB.*t9.*t12.*t34.*t39;
t220 = mArm.*mLeg.*mPB.*t9.*t12.*t33.*t40;
t221 = mBody.*mLeg.*mPB.*t8.*t12.*t34.*t39;
t222 = mBody.*mLeg.*mPB.*t8.*t12.*t33.*t40;
t223 = mArm.*mBody.*mLeg.*rArm.*rArmMCD.*t12.*t33.*t39.*2.0;
t224 = mArm.*mBody.*mLeg.*t8.*t12.*t34.*t40;
t225 = mArm.*mBody.*mLeg.*t9.*t12.*t34.*t40;
t226 = mArm.*mLeg.*mPB.*t9.*t12.*t34.*t40;
t227 = mBody.*mLeg.*mPB.*t8.*t12.*t34.*t40;
t228 = mArm.*mBody.*mLeg.*rArm.*rArmMCD.*t12.*t33.*t40.*2.0;
t247 = mArm.*mBody.*mLeg.*t11.*t12.*t37.*t39;
t248 = mBody.*mLeg.*mPB.*t11.*t12.*t37.*t39;
t249 = mArm.*mBody.*mLeg.*t11.*t12.*t37.*t40;
t250 = mArm.*mBody.*mLeg.*t11.*t12.*t38.*t39;
t251 = mBody.*mLeg.*mPB.*t11.*t12.*t37.*t40;
t252 = mBody.*mLeg.*mPB.*t11.*t12.*t38.*t39;
t253 = mArm.*mBody.*mLeg.*t11.*t12.*t38.*t40;
t254 = mBody.*mLeg.*mPB.*t11.*t12.*t38.*t40;
t256 = mArm.*t8.*t33.*t95;
t257 = mArm.*t9.*t33.*t95;
t260 = mArm.*t8.*t33.*t96;
t261 = mArm.*t9.*t33.*t96;
t268 = mPB.*t8.*t33.*t95;
t270 = mPB.*t8.*t33.*t96;
t271 = mArm.*t9.*t34.*t95;
t273 = mArm.*t9.*t34.*t96;
t279 = mBody.*rArm.*rArmMCD.*t33.*t61.*-2.0;
t280 = mLeg.*rArm.*rArmMCD.*t33.*t60.*-2.0;
t281 = -t195;
t282 = -t196;
t283 = -t197;
t284 = mArm.*rArm.*rArmMCD.*t33.*t95.*2.0;
t285 = mArm.*rArm.*rArmMCD.*t33.*t96.*2.0;
t288 = mBody.*rArm.*rArmMCD.*t33.*t67.*-2.0;
t289 = mLeg.*rArm.*rArmMCD.*t33.*t66.*-2.0;
t290 = -t208;
t293 = mArm.*t7.*t8.*t12.*t34.*t39;
t294 = mArm.*t7.*t8.*t12.*t33.*t40;
t295 = mArm.*t7.*t9.*t12.*t33.*t40;
t296 = mLeg.*t5.*t9.*t12.*t34.*t39;
t297 = mBody.*t7.*t8.*t12.*t34.*t39;
t298 = mLeg.*t6.*t8.*t12.*t34.*t39;
t299 = mPB.*t7.*t8.*t12.*t34.*t39;
t300 = mPB.*t7.*t8.*t12.*t33.*t40;
t301 = -t213;
t302 = -t214;
t303 = mArm.*t7.*t9.*t12.*t34.*t40;
t304 = mLeg.*t5.*t9.*t12.*t34.*t40;
t305 = mBody.*t7.*t8.*t12.*t34.*t40;
t306 = mLeg.*t6.*t8.*t12.*t34.*t40;
t308 = mArm.*rArm.*rArmMCD.*t7.*t12.*t33.*t40.*2.0;
t325 = mArm.*t7.*t10.*t12.*t37.*t40;
t326 = mArm.*t7.*t10.*t12.*t38.*t39;
t327 = mBody.*t7.*t10.*t12.*t37.*t40;
t328 = mBody.*t7.*t10.*t12.*t38.*t39;
t329 = mBody.*t7.*t11.*t12.*t37.*t40;
t330 = mLeg.*t6.*t11.*t12.*t38.*t39;
t331 = mPB.*t7.*t10.*t12.*t37.*t40;
t332 = mPB.*t7.*t10.*t12.*t38.*t39;
t333 = mBody.*t7.*t11.*t12.*t38.*t40;
t334 = mLeg.*t6.*t11.*t12.*t38.*t40;
t335 = mBody.*rBody.*rBodyMCD.*t7.*t12.*t37.*t40.*2.0;
t346 = -t286;
t347 = -t287;
t348 = -t291;
t349 = -t292;
t365 = mArm.*mBody.*mLeg.*rArm.*rBodyMCD.*t12.*t27.*t30.*t39;
t366 = mArm.*mBody.*mLeg.*rArmMCD.*rBodyMCD.*t12.*t27.*t30.*t39;
t367 = mBody.*mLeg.*mPB.*rArm.*rBodyMCD.*t12.*t27.*t30.*t39;
t368 = mArm.*mBody.*mLeg.*rArm.*rBodyMCD.*t12.*t27.*t30.*t40;
t369 = mArm.*mBody.*mLeg.*rArmMCD.*rBodyMCD.*t12.*t27.*t30.*t40;
t370 = mBody.*mLeg.*mPB.*rArm.*rBodyMCD.*t12.*t27.*t30.*t40;
t371 = mArm.*mBody.*mLeg.*rArm.*rBodyMCD.*t12.*t28.*t31.*t39;
t372 = mBody.*mLeg.*mPB.*rArm.*rBodyMCD.*t12.*t28.*t31.*t39;
t373 = mArm.*mBody.*mLeg.*rArm.*rBodyMCD.*t12.*t28.*t31.*t40;
t374 = mBody.*mLeg.*mPB.*rArm.*rBodyMCD.*t12.*t28.*t31.*t40;
t375 = mArm.*rArm.*rBody.*t7.*t12.*t27.*t30.*t40;
t376 = mArm.*rArmMCD.*rBody.*t7.*t12.*t27.*t30.*t40;
t377 = mPB.*rArm.*rBody.*t7.*t12.*t27.*t30.*t40;
t378 = mArm.*rArm.*rBody.*t7.*t12.*t28.*t31.*t39;
t379 = mBody.*rArm.*rBody.*t7.*t12.*t28.*t31.*t39;
t380 = mLeg.*rArm.*rBodyMCD.*t6.*t12.*t28.*t31.*t39;
t381 = mPB.*rArm.*rBody.*t7.*t12.*t28.*t31.*t39;
t382 = mBody.*rArm.*rBodyMCD.*t7.*t12.*t28.*t31.*t40;
t383 = mLeg.*rArm.*rBodyMCD.*t6.*t12.*t28.*t31.*t40;
t413 = mArm.*mBody.*rArm.*rArmMCD.*rBody.*rBodyMCD.*t7.*t12.*t33.*t37.*t40.*4.0;
t415 = InertiaBody.*t412;
t478 = InertiaBody.*t439;
t479 = InertiaBody.*t440;
t480 = -t398;
t481 = -t399;
t482 = -t401;
t483 = -t403;
t484 = -t404;
t485 = t5.*t7.*t9.*t10.*t12.*t34.*t37.*t40;
t486 = t5.*t7.*t9.*t10.*t12.*t34.*t38.*t39;
t487 = t6.*t7.*t8.*t10.*t12.*t34.*t37.*t40;
t488 = t6.*t7.*t8.*t10.*t12.*t34.*t38.*t39;
t489 = t6.*t7.*t8.*t11.*t12.*t34.*t37.*t40;
t490 = t6.*t7.*t8.*t11.*t12.*t34.*t38.*t39;
t492 = -t439;
t493 = -t440;
t501 = rBody.*rBodyMCD.*t6.*t7.*t8.*t12.*t34.*t37.*t40.*2.0;
t502 = rBody.*rBodyMCD.*t6.*t7.*t8.*t12.*t34.*t38.*t39.*2.0;
t508 = InertiaArm.*t491;
t516 = -t511;
t517 = -t512;
t518 = -t513;
t519 = InertiaArm.*t511;
t520 = InertiaArm.*t512;
t521 = InertiaArm.*t513;
t527 = mArm.*mLeg.*rArm.*rArmMCD.*t6.*t11.*t12.*t27.*t28.*t30.*t31.*t39.*2.0;
t528 = mArm.*mBody.*rBody.*rBodyMCD.*t7.*t8.*t12.*t27.*t28.*t30.*t31.*t39.*2.0;
t529 = mBody.*mPB.*rBody.*rBodyMCD.*t7.*t8.*t12.*t27.*t28.*t30.*t31.*t39.*2.0;
t530 = mArm.*mLeg.*rArm.*rArmMCD.*t6.*t11.*t12.*t27.*t28.*t30.*t31.*t40.*2.0;
t531 = mArm.*mBody.*rBody.*rBodyMCD.*t7.*t8.*t12.*t27.*t28.*t30.*t31.*t40.*2.0;
t532 = mBody.*mPB.*rBody.*rBodyMCD.*t7.*t8.*t12.*t27.*t28.*t30.*t31.*t40.*2.0;
t533 = mBody.*rBody.*rBodyMCD.*t37.*t412;
t534 = mBody.*rBody.*rBodyMCD.*t38.*t412;
t535 = mArm.*mLeg.*t6.*t8.*t11.*t12.*t27.*t28.*t30.*t31.*t39.*2.0;
t536 = mLeg.*mPB.*t6.*t8.*t11.*t12.*t27.*t28.*t30.*t31.*t39.*2.0;
t538 = mArm.*mLeg.*t6.*t8.*t11.*t12.*t27.*t28.*t30.*t31.*t40.*2.0;
t539 = mLeg.*mPB.*t6.*t8.*t11.*t12.*t27.*t28.*t30.*t31.*t40.*2.0;
t544 = mBody.*t11.*t37.*t412;
t545 = mBody.*rBody.*rBodyMCD.*t37.*t439;
t546 = mBody.*rBody.*rBodyMCD.*t37.*t440;
t547 = mBody.*t11.*t38.*t412;
t548 = mBody.*rBody.*rBodyMCD.*t38.*t439;
t549 = mBody.*rBody.*rBodyMCD.*t38.*t440;
t553 = mArm.*t8.*t33.*t491;
t554 = mArm.*t9.*t33.*t491;
t555 = mPB.*t8.*t33.*t491;
t558 = mArm.*t8.*t34.*t491;
t559 = mArm.*t9.*t34.*t491;
t560 = mPB.*t8.*t34.*t491;
t561 = -t537;
t564 = mBody.*t11.*t37.*t439;
t565 = mBody.*t11.*t37.*t440;
t566 = mBody.*t11.*t38.*t439;
t567 = mBody.*t11.*t38.*t440;
t568 = mBody.*t8.*t33.*t511;
t569 = mBody.*t9.*t33.*t511;
t570 = mPB.*t9.*t33.*t511;
t571 = mPB.*t8.*t33.*t512;
t572 = mBody.*t8.*t34.*t511;
t573 = mBody.*t9.*t34.*t511;
t574 = mPB.*t9.*t34.*t511;
t575 = mPB.*t8.*t34.*t512;
t592 = -t585;
t593 = -t586;
t594 = -t587;
t134 = InertiaBody.*t111;
t135 = InertiaBody.*t113;
t136 = InertiaBody.*t115;
t137 = InertiaBody.*t117;
t138 = InertiaBody.*t119;
t139 = InertiaBody.*t121;
t140 = InertiaBody.*t132;
t229 = mBody.*t8.*t33.*t112;
t230 = mBody.*t9.*t33.*t112;
t231 = mPB.*t9.*t33.*t112;
t232 = mPB.*t8.*t33.*t114;
t233 = mBody.*t8.*t34.*t112;
t234 = mBody.*t8.*t33.*t118;
t235 = mBody.*t9.*t34.*t112;
t236 = mBody.*t9.*t33.*t118;
t237 = mPB.*t9.*t34.*t112;
t238 = mPB.*t9.*t33.*t118;
t239 = mPB.*t8.*t34.*t114;
t240 = mPB.*t8.*t33.*t120;
t241 = mBody.*rArm.*rArmMCD.*t33.*t112.*2.0;
t242 = mBody.*t8.*t34.*t118;
t243 = mBody.*t9.*t34.*t118;
t244 = mPB.*t9.*t34.*t118;
t245 = mPB.*t8.*t34.*t120;
t246 = mBody.*rArm.*rArmMCD.*t33.*t118.*2.0;
t307 = -t223;
t309 = -t228;
t310 = InertiaBody.*t293;
t311 = mArm.*t8.*t33.*t133;
t312 = mArm.*t9.*t33.*t133;
t313 = InertiaBody.*t296;
t314 = InertiaBody.*t297;
t315 = InertiaBody.*t298;
t316 = InertiaBody.*t299;
t317 = mPB.*t8.*t33.*t133;
t318 = mArm.*t9.*t34.*t133;
t319 = InertiaBody.*t304;
t320 = mBody.*t8.*t34.*t133;
t321 = InertiaBody.*t306;
t323 = mArm.*rArm.*rArmMCD.*t33.*t133.*2.0;
t336 = mBody.*t11.*t37.*t111;
t337 = mPB.*t11.*t37.*t113;
t338 = mBody.*t11.*t37.*t117;
t339 = mBody.*t11.*t38.*t111;
t340 = mPB.*t11.*t37.*t119;
t341 = mPB.*t11.*t38.*t113;
t342 = mBody.*t11.*t38.*t117;
t343 = mPB.*t11.*t38.*t119;
t344 = -t284;
t345 = -t285;
t350 = -t308;
t352 = -t335;
t353 = mArm.*t10.*t37.*t132;
t354 = InertiaArm.*t326;
t355 = mBody.*t10.*t37.*t132;
t356 = InertiaArm.*t328;
t357 = mBody.*t11.*t37.*t132;
t358 = InertiaArm.*t330;
t359 = mPB.*t10.*t37.*t132;
t360 = InertiaArm.*t332;
t361 = mBody.*t11.*t38.*t132;
t362 = InertiaArm.*t334;
t363 = mBody.*rBody.*rBodyMCD.*t37.*t132.*2.0;
t384 = -t366;
t385 = -t369;
t386 = -t376;
t405 = mPB.*t11.*t37.*t210;
t406 = mPB.*t11.*t37.*t217;
t407 = mPB.*t11.*t37.*t218;
t408 = mPB.*t11.*t38.*t210;
t409 = mPB.*t11.*t37.*t225;
t410 = mPB.*t11.*t38.*t217;
t411 = mPB.*t11.*t38.*t218;
t414 = mPB.*t11.*t38.*t225;
t420 = mBody.*t10.*t37.*t294;
t421 = mBody.*t8.*t33.*t326;
t422 = mBody.*t11.*t37.*t293;
t423 = mBody.*t11.*t37.*t294;
t424 = mBody.*t10.*t37.*t295;
t425 = mBody.*t9.*t33.*t326;
t426 = mArm.*t11.*t37.*t298;
t427 = mArm.*t8.*t33.*t330;
t428 = mBody.*t11.*t37.*t295;
t429 = mArm.*t9.*t33.*t330;
t430 = mBody.*t11.*t37.*t296;
t431 = mPB.*t10.*t37.*t295;
t432 = mPB.*t9.*t33.*t326;
t433 = mBody.*t10.*t37.*t300;
t434 = mPB.*t8.*t33.*t328;
t435 = mPB.*t11.*t37.*t297;
t436 = mBody.*t11.*t37.*t300;
t437 = mPB.*t11.*t37.*t298;
t438 = mPB.*t8.*t33.*t330;
t441 = mArm.*t10.*t37.*t305;
t442 = mBody.*t10.*t38.*t293;
t443 = mBody.*t11.*t38.*t293;
t444 = mBody.*t11.*t38.*t294;
t445 = mBody.*t10.*t37.*t303;
t446 = mBody.*t9.*t34.*t326;
t447 = mArm.*t11.*t37.*t306;
t448 = mArm.*t8.*t33.*t334;
t449 = mBody.*t11.*t37.*t303;
t450 = mBody.*t11.*t38.*t295;
t451 = mArm.*t9.*t34.*t330;
t452 = mArm.*t9.*t33.*t334;
t453 = mBody.*t11.*t37.*t304;
t454 = mBody.*t11.*t38.*t296;
t455 = mPB.*t10.*t37.*t303;
t456 = mPB.*t9.*t34.*t326;
t457 = mPB.*t10.*t37.*t305;
t458 = mPB.*t10.*t38.*t297;
t459 = mPB.*t11.*t38.*t297;
t460 = mBody.*t11.*t38.*t300;
t461 = mPB.*t11.*t37.*t306;
t462 = mPB.*t8.*t33.*t334;
t463 = mBody.*t10.*t37.*t308;
t464 = mBody.*rArm.*rArmMCD.*t33.*t326.*2.0;
t465 = mBody.*t11.*t37.*t308;
t466 = mArm.*rArm.*rArmMCD.*t33.*t330.*2.0;
t467 = mBody.*rBody.*rBodyMCD.*t37.*t294.*2.0;
t468 = mBody.*rBody.*rBodyMCD.*t37.*t295.*2.0;
t469 = mBody.*rBody.*rBodyMCD.*t37.*t300.*2.0;
t470 = mBody.*t11.*t38.*t303;
t471 = mArm.*t9.*t34.*t334;
t472 = mBody.*t11.*t38.*t304;
t473 = mBody.*t11.*t38.*t308;
t474 = mArm.*rArm.*rArmMCD.*t33.*t334.*2.0;
t475 = mBody.*rBody.*rBodyMCD.*t38.*t293.*2.0;
t476 = mBody.*rBody.*rBodyMCD.*t37.*t303.*2.0;
t477 = mPB.*rBody.*rBodyMCD.*t38.*t297.*2.0;
t494 = mBody.*rArm.*rArmMCD.*t33.*t325.*-2.0;
t496 = mArm.*rArm.*rArmMCD.*t33.*t329.*-2.0;
t503 = mArm.*rArm.*rArmMCD.*t33.*t333.*-2.0;
t514 = -t501;
t515 = -t502;
t525 = mBody.*rArmMCD.*rBodyMCD.*t27.*t30.*t378.*2.0;
t526 = mBody.*rArmMCD.*rBodyMCD.*t28.*t31.*t375.*2.0;
t540 = -t528;
t541 = -t529;
t542 = -t531;
t543 = -t532;
t556 = -t535;
t557 = -t536;
t562 = -t538;
t563 = -t539;
t322 = -t241;
t324 = -t246;
t351 = -t323;
t364 = -t363;
t495 = -t464;
t497 = -t466;
t498 = -t467;
t499 = -t468;
t500 = -t469;
t504 = -t474;
t505 = -t475;
t506 = -t476;
t507 = -t477;
t595 = t21+t22+t23+t24+t79+t80+t81+t82+t83+t84+t85+t86+t87+t88+t89+t90+t91+t92+t98+t99+t100+t101+t102+t103+t104+t105+t106+t107+t108+t109+t123+t124+t125+t126+t127+t128+t129+t130+t131+t134+t135+t136+t137+t138+t139+t140+t168+t169+t170+t171+t172+t173+t174+t175+t176+t177+t178+t179+t180+t181+t182+t183+t184+t185+t186+t187+t188+t189+t190+t191+t192+t198+t199+t200+t201+t202+t203+t204+t205+t229+t230+t231+t232+t233+t234+t235+t236+t237+t238+t239+t240+t242+t243+t244+t245+t255+t256+t257+t258+t259+t260+t261+t262+t263+t264+t265+t266+t267+t268+t269+t270+t271+t272+t273+t274+t275+t276+t277+t278+t279+t280+t281+t282+t283+t288+t289+t290+t301+t302+t310+t311+t312+t313+t314+t315+t316+t317+t318+t319+t320+t321+t322+t324+t336+t337+t338+t339+t340+t341+t342+t343+t344+t345+t346+t347+t348+t349+t351+t353+t354+t355+t356+t357+t358+t359+t360+t361+t362+t364+t387+t390+t391+t392+t393+t405+t406+t407+t408+t409+t410+t411+t413+t414+t415+t416+t417+t418+t419+t420+t421+t422+t423+t424+t425+t426+t427+t428+t429+t430+t431+t432+t433+t434+t435+t436+t437+t438+t441+t442+t443+t444+t445+t446+t447+t448+t449+t450+t451+t452+t453+t454+t455+t456+t457+t458+t459+t460+t461+t462+t470+t471+t472+t485+t486+t487+t488+t489+t490+t494+t495+t496+t497+t498+t499+t500+t503+t504+t505+t506+t507+t508+t509+t510+t514+t515+t522+t523+t524+t525+t526+t527+t530+t540+t541+t542+t543+t544+t545+t546+t547+t548+t549+t550+t551+t552+t553+t554+t555+t556+t557+t558+t559+t560+t561+t562+t563+t576+t577+t578+t579+t580+t581+t582+t583+t584+t588+t589+t590+t591+t592+t593+t594;
t596 = 1.0./t595;
ddthShoulder = t596.*(MrPB-g.*mArm-g.*mBody-g.*mLeg-g.*mPB+mArm.*rArmMCD.*t2.*t28+mBody.*rArm.*t2.*t28+mBody.*rBodyMCD.*t2.*t31+mBody.*rBodyMCD.*t4.*t31+mLeg.*rArm.*t2.*t28+mLeg.*rBody.*t2.*t31+mLeg.*rBody.*t4.*t31+mLeg.*rLegMCD.*t2.*t36+mLeg.*rLegMCD.*t3.*t36+mLeg.*rLegMCD.*t4.*t36+dthHand.*dthShoulder.*mBody.*rBodyMCD.*t31.*2.0+dthHand.*dthShoulder.*mLeg.*rBody.*t31.*2.0+dthHand.*dthWaist.*mLeg.*rLegMCD.*t36.*2.0+dthHand.*dthShoulder.*mLeg.*rLegMCD.*t36.*2.0+dthWaist.*dthShoulder.*mLeg.*rLegMCD.*t36.*2.0).*(rArm.*t16.*t27+rArm.*t18.*t27+rArm.*t27.*t62+rArm.*t27.*t63+rArm.*t27.*t68+rArm.*t27.*t69+rArm.*t27.*t95+rArm.*t27.*t96+rArm.*t27.*t114+rArm.*t27.*t120+rArm.*t27.*t133+rArmMCD.*t14.*t27+rArmMCD.*t27.*t60+rArmMCD.*t27.*t61+rArmMCD.*t27.*t66+rArmMCD.*t27.*t67+rArmMCD.*t27.*t112+rArm.*t27.*t327+rArm.*t27.*t328+rArm.*t27.*t329+rArm.*t27.*t330+rArmMCD.*t27.*t118+rArm.*t27.*t333+rArm.*t27.*t334+rArmMCD.*t27.*t147+rArmMCD.*t27.*t148+rArmMCD.*t27.*t247+rArmMCD.*t27.*t249+rArmMCD.*t27.*t250+rArmMCD.*t27.*t253+rArm.*t27.*t491+rArmMCD.*t27.*t325+rArmMCD.*t27.*t326+rArmMCD.*t27.*t371+rArmMCD.*t27.*t373+rArmMCD.*t27.*t378+rArmMCD.*t27.*t480-rBody.*t17.*t30-rBody.*t30.*t47-rBody.*t30.*t54-rBody.*t30.*t75-rBody.*t30.*t132-rBodyMCD.*t15.*t30-rBodyMCD.*t30.*t45-rBodyMCD.*t30.*t52-rBodyMCD.*t30.*t74-rBody.*t30.*t295-rBody.*t30.*t303-rBody.*t30.*t305-rBodyMCD.*t30.*t113-rBodyMCD.*t30.*t119-rBodyMCD.*t30.*t149.*2.0-rBody.*t30.*t382-rBodyMCD.*t30.*t210-rBodyMCD.*t30.*t217-rBodyMCD.*t30.*t218-rBodyMCD.*t30.*t225-rBodyMCD.*t30.*t297-rBodyMCD.*t30.*t298-rBodyMCD.*t30.*t306-rBodyMCD.*t30.*t379-InertiaLeg.*mBody.*mLeg.*rBody.*t8.*t30.*t34-InertiaLeg.*mBody.*mLeg.*rBodyMCD.*t8.*t30.*t34-InertiaLeg.*rArm.*t6.*t11.*t28.*t30.*t31-InertiaLeg.*rArm.*t7.*t10.*t28.*t30.*t31-InertiaBody.*rArm.*t7.*t12.*t28.*t35.*t36+InertiaLeg.*rBody.*t7.*t8.*t27.*t28.*t31+InertiaArm.*rBody.*t7.*t12.*t31.*t35.*t36+InertiaLeg.*rBodyMCD.*t6.*t8.*t27.*t28.*t31+InertiaLeg.*mArm.*mBody.*rArm.*rArmMCD.*rBodyMCD.*t30.*t33+InertiaLeg.*mArm.*mLeg.*rArm.*rArmMCD.*rBody.*t30.*t33-InertiaLeg.*mBody.*mLeg.*rArm.*rBody.*rBodyMCD.*t27.*t37.*2.0+InertiaLeg.*mBody.*mLeg.*rBody.*t8.*t27.*t28.*t31+InertiaLeg.*mBody.*mLeg.*rBodyMCD.*t8.*t27.*t28.*t31+mArm.*rArm.*rArmMCD.*rBody.*t7.*t12.*t30.*t33.*t40-mBody.*rArm.*rBody.*rBodyMCD.*t7.*t12.*t27.*t37.*t40.*2.0+mArm.*rBody.*t7.*t9.*t12.*t31.*t33.*t35.*t36+mArm.*rBody.*t7.*t9.*t12.*t31.*t34.*t35.*t36-mBody.*rArm.*t7.*t11.*t12.*t28.*t35.*t36.*t37-mBody.*rArm.*t7.*t11.*t12.*t28.*t35.*t36.*t38+mBody.*rBody.*t7.*t8.*t12.*t27.*t28.*t31.*t39+mBody.*rBody.*t7.*t8.*t12.*t31.*t34.*t35.*t36+mBody.*rBodyMCD.*t7.*t8.*t12.*t27.*t28.*t31.*t40-mBody.*rBodyMCD.*t7.*t8.*t12.*t31.*t34.*t35.*t36-mLeg.*rArm.*t6.*t11.*t12.*t28.*t30.*t31.*t39-mLeg.*rArm.*t6.*t11.*t12.*t28.*t30.*t31.*t40+mLeg.*rBodyMCD.*t6.*t8.*t12.*t27.*t28.*t31.*t39+mLeg.*rBodyMCD.*t6.*t8.*t12.*t27.*t28.*t31.*t40+mArm.*mBody.*mLeg.*rArm.*rArmMCD.*rBodyMCD.*t12.*t30.*t33.*t39+mArm.*mBody.*mLeg.*rArm.*rArmMCD.*rBodyMCD.*t12.*t30.*t33.*t40-mArm.*rArm.*rArmMCD.*rBody.*t7.*t12.*t31.*t33.*t35.*t36+mBody.*rArm.*rBody.*rBodyMCD.*t7.*t12.*t28.*t35.*t36.*t37+mBody.*rArm.*rBody.*rBodyMCD.*t7.*t12.*t28.*t35.*t36.*t38-mArm.*rArmMCD.*t7.*t10.*t12.*t27.*t30.*t31.*t35.*t36.*2.0-mBody.*rArm.*t7.*t10.*t12.*t27.*t30.*t31.*t35.*t36.*2.0-mBody.*rBody.*t7.*t8.*t12.*t27.*t28.*t30.*t35.*t36+mBody.*rBodyMCD.*t7.*t8.*t12.*t27.*t28.*t30.*t35.*t36)+t596.*(-MthHand+t41+t42+t43+t160+t162+t164+t165+t166+t167+g.*mArm.*rArmMCD.*t27+g.*mBody.*rArm.*t27+g.*mLeg.*rArm.*t27-mBody.*rArm.*rBodyMCD.*t4.*t27.*t31+mBody.*rArm.*rBodyMCD.*t4.*t28.*t30-mLeg.*rArm.*rBody.*t4.*t27.*t31+mLeg.*rArm.*rBody.*t4.*t28.*t30-mLeg.*rArm.*rLegMCD.*t3.*t27.*t36+mLeg.*rArm.*rLegMCD.*t3.*t28.*t35-mLeg.*rArm.*rLegMCD.*t4.*t27.*t36+mLeg.*rArm.*rLegMCD.*t4.*t28.*t35-dthHand.*dthShoulder.*mBody.*rArm.*rBodyMCD.*t27.*t31.*2.0+dthHand.*dthShoulder.*mBody.*rArm.*rBodyMCD.*t28.*t30.*2.0-dthHand.*dthShoulder.*mLeg.*rArm.*rBody.*t27.*t31.*2.0+dthHand.*dthShoulder.*mLeg.*rArm.*rBody.*t28.*t30.*2.0-dthHand.*dthWaist.*mLeg.*rArm.*rLegMCD.*t27.*t36.*2.0+dthHand.*dthWaist.*mLeg.*rArm.*rLegMCD.*t28.*t35.*2.0-dthHand.*dthShoulder.*mLeg.*rArm.*rLegMCD.*t27.*t36.*2.0+dthHand.*dthShoulder.*mLeg.*rArm.*rLegMCD.*t28.*t35.*2.0-dthWaist.*dthShoulder.*mLeg.*rArm.*rLegMCD.*t27.*t36.*2.0+dthWaist.*dthShoulder.*mLeg.*rArm.*rLegMCD.*t28.*t35.*2.0).*(t14+t16+t18+t20+t60+t61+t62+t63+t64+t65+t66+t67+t68+t69+t70+t71+t95+t96+t97+t112+t114+t116+t118+t120+t122+t133+t141+t143+t145+t146+t147+t148+t149+t150+t151+t152+t153+t154+t157+t158+t247+t248+t249+t250+t251+t252+t253+t254+t325+t326+t327+t328+t329+t330+t331+t332+t333+t334+t352+t365+t367+t368+t370+t371+t372+t373+t374+t375+t377+t378+t379+t380+t381+t382+t383+t384+t385+t386+t400+t402+t480+t481+t482+t483+t484+t491+t516+t517+t518)-t596.*(-MthShoulder+t41+t42+t43-t155+t156+t160+t162+t164+t165+t166+t167+mBody.*rArm.*rBodyMCD.*t2.*t27.*t31-mBody.*rArm.*rBodyMCD.*t2.*t28.*t30+mLeg.*rArm.*rBody.*t2.*t27.*t31-mLeg.*rArm.*rBody.*t2.*t28.*t30).*(t13+t14+t15+t16+t17+t18+t19+t20+t44+t45+t46+t47+t48+t49+t50+t51+t52+t53+t54+t55+t56+t57+t60+t61+t62+t63+t64+t65+t66+t67+t68+t69+t70+t71+t73+t74+t75+t76+t77+t78+t95+t96+t97+t111+t112+t113+t114+t115+t116+t117+t118+t119+t120+t121+t122+t132+t133+t141.*2.0-t142.*2.0+t143.*2.0-t144.*2.0+t145.*2.0+t146.*2.0+t147.*2.0+t148.*2.0+t149.*2.0+t150.*2.0+t151.*2.0+t152.*2.0+t153.*2.0+t154.*2.0+t209+t210+t211+t212+t215+t216+t217+t218+t219+t220+t221+t222+t224+t225+t226+t227+t247+t248+t249+t250+t251+t252+t253+t254+t293+t294+t295+t296+t297+t298+t299+t300+t303+t304+t305+t306+t307+t309+t325+t326+t327+t328+t329+t330+t331+t332+t333+t334+t350+t352+t365.*2.0-t366.*2.0+t367.*2.0+t368.*2.0-t369.*2.0+t370.*2.0+t371.*2.0+t372.*2.0+t373.*2.0+t374.*2.0+t375.*2.0-t376.*2.0+t377.*2.0+t378.*2.0+t379.*2.0+t380.*2.0+t381.*2.0+t382.*2.0+t383.*2.0-t398.*2.0-t399.*2.0+t400.*2.0-t401.*2.0+t402.*2.0-t403.*2.0-t404.*2.0+t412+t491+t492+t493+t516+t517+t518)-t596.*(MthWaist-t43+t155-t156-mLeg.*rBody.*rLegMCD.*t2.*t30.*t36+mLeg.*rBody.*rLegMCD.*t2.*t31.*t35-mLeg.*rBody.*rLegMCD.*t4.*t30.*t36+mLeg.*rBody.*rLegMCD.*t4.*t31.*t35-dthHand.*dthShoulder.*mLeg.*rBody.*rLegMCD.*t30.*t36.*2.0+dthHand.*dthShoulder.*mLeg.*rBody.*rLegMCD.*t31.*t35.*2.0).*(t13+t15+t17+t19+t44+t45+t46+t47+t48+t49+t50+t51+t52+t53+t54+t55+t56+t57+t73+t74+t75+t76+t77+t78+t111+t113+t115+t117+t119+t121+t132+t141+t143+t145+t146+t147+t148+t149+t150+t151+t152+t153+t154+t157+t158+t209+t210+t211+t212+t215+t216+t217+t218+t219+t220+t221+t222+t224+t225+t226+t227+t293+t294+t295+t296+t297+t298+t299+t300+t303+t304+t305+t306+t307+t309+t350+t365+t367+t368+t370+t371+t372+t373+t374+t375+t377+t378+t379+t380+t381+t382+t383+t384+t385+t386+t400+t402+t412+t480+t481+t482+t483+t484+t492+t493-InertiaBody.*rArm.*rLegMCD.*t7.*t28.*t36+InertiaArm.*rBody.*rLegMCD.*t7.*t31.*t36-InertiaBody.*mArm.*mLeg.*rArm.*rLegMCD.*t27.*t35-InertiaBody.*mArm.*mLeg.*rArm.*rLegMCD.*t28.*t36+InertiaBody.*mArm.*mLeg.*rArmMCD.*rLegMCD.*t27.*t35+InertiaArm.*mArm.*mLeg.*rBody.*rLegMCD.*t30.*t35+InertiaArm.*mArm.*mLeg.*rBody.*rLegMCD.*t31.*t36-InertiaBody.*mBody.*mLeg.*rArm.*rLegMCD.*t28.*t36+InertiaArm.*mBody.*mLeg.*rBody.*rLegMCD.*t30.*t35+InertiaArm.*mBody.*mLeg.*rBody.*rLegMCD.*t31.*t36-InertiaArm.*mBody.*mLeg.*rBodyMCD.*rLegMCD.*t30.*t35-InertiaBody.*mLeg.*mPB.*rArm.*rLegMCD.*t27.*t35-InertiaBody.*mLeg.*mPB.*rArm.*rLegMCD.*t28.*t36+InertiaArm.*mLeg.*mPB.*rBody.*rLegMCD.*t30.*t35+InertiaArm.*mLeg.*mPB.*rBody.*rLegMCD.*t31.*t36-mArm.*rArm.*rLegMCD.*t7.*t10.*t27.*t35.*t38-mArm.*rArm.*rLegMCD.*t7.*t10.*t28.*t36.*t37+mArm.*rArmMCD.*rLegMCD.*t7.*t10.*t27.*t35.*t38+mArm.*rBody.*rLegMCD.*t7.*t8.*t30.*t34.*t35+mArm.*rBody.*rLegMCD.*t7.*t8.*t31.*t33.*t36+mArm.*rBody.*rLegMCD.*t7.*t9.*t31.*t33.*t36+mArm.*rBody.*rLegMCD.*t7.*t9.*t31.*t34.*t36-mBody.*rArm.*rLegMCD.*t7.*t10.*t28.*t36.*t37-mBody.*rArm.*rLegMCD.*t7.*t11.*t28.*t36.*t37-mBody.*rArm.*rLegMCD.*t7.*t11.*t28.*t36.*t38+mBody.*rBody.*rLegMCD.*t7.*t8.*t30.*t34.*t35+mBody.*rBody.*rLegMCD.*t7.*t8.*t31.*t34.*t36-mBody.*rBodyMCD.*rLegMCD.*t7.*t8.*t30.*t34.*t35-mBody.*rBodyMCD.*rLegMCD.*t7.*t8.*t31.*t34.*t36-mLeg.*rArm.*rLegMCD.*t6.*t11.*t28.*t36.*t38+mLeg.*rBody.*rLegMCD.*t5.*t9.*t30.*t34.*t35+mLeg.*rBody.*rLegMCD.*t6.*t8.*t30.*t34.*t35+mLeg.*rBody.*rLegMCD.*t5.*t9.*t31.*t34.*t36+mLeg.*rBody.*rLegMCD.*t6.*t8.*t31.*t34.*t36-mLeg.*rBodyMCD.*rLegMCD.*t6.*t8.*t30.*t34.*t35-mLeg.*rBodyMCD.*rLegMCD.*t6.*t8.*t31.*t34.*t36-mPB.*rArm.*rLegMCD.*t7.*t10.*t27.*t35.*t38-mPB.*rArm.*rLegMCD.*t7.*t10.*t28.*t36.*t37+mPB.*rBody.*rLegMCD.*t7.*t8.*t30.*t34.*t35+mPB.*rBody.*rLegMCD.*t7.*t8.*t31.*t33.*t36-mArm.*mBody.*mLeg.*rArm.*rLegMCD.*t11.*t27.*t35.*t37-mArm.*mBody.*mLeg.*rArm.*rLegMCD.*t11.*t27.*t35.*t38-mArm.*mBody.*mLeg.*rArm.*rLegMCD.*t11.*t28.*t36.*t37-mArm.*mBody.*mLeg.*rArm.*rLegMCD.*t11.*t28.*t36.*t38+mArm.*mBody.*mLeg.*rArmMCD.*rLegMCD.*t11.*t27.*t35.*t37+mArm.*mBody.*mLeg.*rArmMCD.*rLegMCD.*t11.*t27.*t35.*t38+mArm.*mBody.*mLeg.*rBody.*rLegMCD.*t8.*t30.*t33.*t35+mArm.*mBody.*mLeg.*rBody.*rLegMCD.*t8.*t30.*t34.*t35+mArm.*mBody.*mLeg.*rBody.*rLegMCD.*t9.*t30.*t33.*t35+mArm.*mBody.*mLeg.*rBody.*rLegMCD.*t8.*t31.*t33.*t36+mArm.*mBody.*mLeg.*rBody.*rLegMCD.*t9.*t30.*t34.*t35+mArm.*mBody.*mLeg.*rBody.*rLegMCD.*t8.*t31.*t34.*t36+mArm.*mBody.*mLeg.*rBody.*rLegMCD.*t9.*t31.*t33.*t36+mArm.*mBody.*mLeg.*rBody.*rLegMCD.*t9.*t31.*t34.*t36-mArm.*mBody.*mLeg.*rBodyMCD.*rLegMCD.*t8.*t30.*t33.*t35-mArm.*mBody.*mLeg.*rBodyMCD.*rLegMCD.*t9.*t30.*t33.*t35-mArm.*mBody.*mLeg.*rBodyMCD.*rLegMCD.*t9.*t30.*t34.*t35-mArm.*mBody.*mLeg.*rBodyMCD.*rLegMCD.*t8.*t31.*t34.*t36+mArm.*mLeg.*mPB.*rBody.*rLegMCD.*t9.*t30.*t33.*t35+mArm.*mLeg.*mPB.*rBody.*rLegMCD.*t9.*t30.*t34.*t35+mArm.*mLeg.*mPB.*rBody.*rLegMCD.*t9.*t31.*t33.*t36+mArm.*mLeg.*mPB.*rBody.*rLegMCD.*t9.*t31.*t34.*t36-mBody.*mLeg.*mPB.*rArm.*rLegMCD.*t11.*t27.*t35.*t37-mBody.*mLeg.*mPB.*rArm.*rLegMCD.*t11.*t27.*t35.*t38-mBody.*mLeg.*mPB.*rArm.*rLegMCD.*t11.*t28.*t36.*t37-mBody.*mLeg.*mPB.*rArm.*rLegMCD.*t11.*t28.*t36.*t38+mBody.*mLeg.*mPB.*rBody.*rLegMCD.*t8.*t30.*t33.*t35+mBody.*mLeg.*mPB.*rBody.*rLegMCD.*t8.*t30.*t34.*t35+mBody.*mLeg.*mPB.*rBody.*rLegMCD.*t8.*t31.*t33.*t36+mBody.*mLeg.*mPB.*rBody.*rLegMCD.*t8.*t31.*t34.*t36-mBody.*mLeg.*mPB.*rBodyMCD.*rLegMCD.*t8.*t30.*t33.*t35-mBody.*mLeg.*mPB.*rBodyMCD.*rLegMCD.*t8.*t31.*t34.*t36-mArm.*rArm.*rArmMCD.*rBody.*rLegMCD.*t7.*t31.*t33.*t36.*2.0+mBody.*rArm.*rBody.*rBodyMCD.*rLegMCD.*t7.*t28.*t36.*t37.*2.0+mBody.*rArm.*rBody.*rBodyMCD.*rLegMCD.*t7.*t28.*t36.*t38+mLeg.*rArm.*rBody.*rBodyMCD.*rLegMCD.*t6.*t28.*t36.*t38+mArm.*rArm.*rLegMCD.*t7.*t10.*t27.*t30.*t31.*t36+mArm.*rArm.*rLegMCD.*t7.*t10.*t28.*t30.*t31.*t35-mArm.*rArmMCD.*rLegMCD.*t7.*t10.*t27.*t30.*t31.*t36-mArm.*rBody.*rLegMCD.*t7.*t8.*t27.*t28.*t30.*t36-mArm.*rBody.*rLegMCD.*t7.*t8.*t27.*t28.*t31.*t35+mBody.*rArm.*rLegMCD.*t7.*t10.*t28.*t30.*t31.*t35-mLeg.*rArm.*rLegMCD.*t6.*t11.*t28.*t30.*t31.*t35+mPB.*rArm.*rLegMCD.*t7.*t10.*t27.*t30.*t31.*t36+mPB.*rArm.*rLegMCD.*t7.*t10.*t28.*t30.*t31.*t35-mPB.*rBody.*rLegMCD.*t7.*t8.*t27.*t28.*t30.*t36-mPB.*rBody.*rLegMCD.*t7.*t8.*t27.*t28.*t31.*t35-mArm.*mBody.*mLeg.*rArm.*rArmMCD.*rBody.*rLegMCD.*t30.*t33.*t35.*2.0-mArm.*mBody.*mLeg.*rArm.*rArmMCD.*rBody.*rLegMCD.*t31.*t33.*t36.*2.0+mArm.*mBody.*mLeg.*rArm.*rArmMCD.*rBodyMCD.*rLegMCD.*t30.*t33.*t35.*2.0+mArm.*mBody.*mLeg.*rArm.*rBody.*rBodyMCD.*rLegMCD.*t27.*t35.*t37+mArm.*mBody.*mLeg.*rArm.*rBody.*rBodyMCD.*rLegMCD.*t28.*t36.*t38-mArm.*mBody.*mLeg.*rArmMCD.*rBody.*rBodyMCD.*rLegMCD.*t27.*t35.*t37+mBody.*mLeg.*mPB.*rArm.*rBody.*rBodyMCD.*rLegMCD.*t27.*t35.*t37+mBody.*mLeg.*mPB.*rArm.*rBody.*rBodyMCD.*rLegMCD.*t28.*t36.*t38-mArm.*mBody.*mLeg.*rBodyMCD.*rLegMCD.*t8.*t27.*t28.*t30.*t36-mArm.*mBody.*mLeg.*rBodyMCD.*rLegMCD.*t8.*t27.*t28.*t31.*t35-mBody.*mLeg.*mPB.*rBodyMCD.*rLegMCD.*t8.*t27.*t28.*t30.*t36-mBody.*mLeg.*mPB.*rBodyMCD.*rLegMCD.*t8.*t27.*t28.*t31.*t35+mArm.*rArm.*rArmMCD.*rBody.*rLegMCD.*t7.*t27.*t28.*t30.*t36+mArm.*rArm.*rArmMCD.*rBody.*rLegMCD.*t7.*t27.*t28.*t31.*t35-mBody.*rArm.*rBody.*rBodyMCD.*rLegMCD.*t7.*t28.*t30.*t31.*t35+mLeg.*rArm.*rBody.*rBodyMCD.*rLegMCD.*t6.*t28.*t30.*t31.*t35+mArm.*mBody.*mLeg.*rArm.*rArmMCD.*rBodyMCD.*rLegMCD.*t27.*t28.*t30.*t36+mArm.*mBody.*mLeg.*rArm.*rArmMCD.*rBodyMCD.*rLegMCD.*t27.*t28.*t31.*t35+mArm.*mBody.*mLeg.*rArm.*rBody.*rBodyMCD.*rLegMCD.*t27.*t30.*t31.*t36+mArm.*mBody.*mLeg.*rArm.*rBody.*rBodyMCD.*rLegMCD.*t28.*t30.*t31.*t35-mArm.*mBody.*mLeg.*rArmMCD.*rBody.*rBodyMCD.*rLegMCD.*t27.*t30.*t31.*t36+mBody.*mLeg.*mPB.*rArm.*rBody.*rBodyMCD.*rLegMCD.*t27.*t30.*t31.*t36+mBody.*mLeg.*mPB.*rArm.*rBody.*rBodyMCD.*rLegMCD.*t28.*t30.*t31.*t35);