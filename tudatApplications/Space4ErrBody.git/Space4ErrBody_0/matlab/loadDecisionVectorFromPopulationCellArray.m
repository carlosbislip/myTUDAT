function [ population ] = loadDecisionVectorFromPopulationCellArray( populationData )

decisionVector = cell2mat(populationData(:,3:end));

population.nodeIntervalAscent1                    = decisionVector(:,1);
population.nodeIntervalAscent2                    = decisionVector(:,2);
population.nodeIntervalAscent3                    = decisionVector(:,3);
population.nodeIntervalAscent4                    = decisionVector(:,4);
population.nodeIntervalAscent5                    = decisionVector(:,5);
population.nodeIntervalAscent6                    = decisionVector(:,6);
population.nodeIntervalAscent7                    = decisionVector(:,7);
population.nodeIntervalAscent8                    = decisionVector(:,8);
population.nodeIntervalAscent9                    = decisionVector(:,9);

population.angleOfAttackAscent1                   = decisionVector(:,10);
population.angleOfAttackAscent2                   = decisionVector(:,11);
population.angleOfAttackAscent3                   = decisionVector(:,12);
population.angleOfAttackAscent4                   = decisionVector(:,13);
population.angleOfAttackAscent5                   = decisionVector(:,14);
population.angleOfAttackAscent6                   = decisionVector(:,15);
population.angleOfAttackAscent7                   = decisionVector(:,16);
population.angleOfAttackAscent8                   = decisionVector(:,17);
population.angleOfAttackAscent9                   = decisionVector(:,18);
population.angleOfAttackAscent10                  = decisionVector(:,19);

population.bankAngleAscent1                       = decisionVector(:,20);
population.bankAngleAscent2                       = decisionVector(:,21);
population.bankAngleAscent3                       = decisionVector(:,22);
population.bankAngleAscent4                       = decisionVector(:,23);
population.bankAngleAscent5                       = decisionVector(:,24);
population.bankAngleAscent6                       = decisionVector(:,25);
population.bankAngleAscent7                       = decisionVector(:,26);
population.bankAngleAscent8                       = decisionVector(:,27);
population.bankAngleAscent9                       = decisionVector(:,28);
population.bankAngleAscent10                      = decisionVector(:,29);

population.thrustElevationAngleAscent1            = decisionVector(:,30);
population.thrustElevationAngleAscent2            = decisionVector(:,31);
population.thrustElevationAngleAscent3            = decisionVector(:,32);
population.thrustElevationAngleAscent4            = decisionVector(:,33);
population.thrustElevationAngleAscent5            = decisionVector(:,34);
population.thrustElevationAngleAscent6            = decisionVector(:,35);
population.thrustElevationAngleAscent7            = decisionVector(:,36);
population.thrustElevationAngleAscent8            = decisionVector(:,37);
population.thrustElevationAngleAscent9            = decisionVector(:,38);
population.thrustElevationAngleAscent10           = decisionVector(:,39);

population.thrustAzimuthAngleAscent1              = decisionVector(:,40);
population.thrustAzimuthAngleAscent2              = decisionVector(:,41);
population.thrustAzimuthAngleAscent3              = decisionVector(:,42);
population.thrustAzimuthAngleAscent4              = decisionVector(:,43);
population.thrustAzimuthAngleAscent5              = decisionVector(:,44);
population.thrustAzimuthAngleAscent6              = decisionVector(:,45);
population.thrustAzimuthAngleAscent7              = decisionVector(:,46);
population.thrustAzimuthAngleAscent8              = decisionVector(:,47);
population.thrustAzimuthAngleAscent9              = decisionVector(:,48);
population.thrustAzimuthAngleAscent10             = decisionVector(:,49);

population.throttleSettingAscent1                 = decisionVector(:,50);
population.throttleSettingAscent2                 = decisionVector(:,51);
population.throttleSettingAscent3                 = decisionVector(:,52);
population.throttleSettingAscent4                 = decisionVector(:,53);
population.throttleSettingAscent5                 = decisionVector(:,54);
population.throttleSettingAscent6                 = decisionVector(:,55);
population.throttleSettingAscent7                 = decisionVector(:,56);
population.throttleSettingAscent8                 = decisionVector(:,57);
population.throttleSettingAscent9                 = decisionVector(:,58);
population.throttleSettingAscent10                = decisionVector(:,59);

population.initialFlightPathAngle                 = decisionVector(:,60);
population.initialLaunchHeading                   = decisionVector(:,61);
population.initialVelocity                        = decisionVector(:,62);
population.maximumVelocity                        = decisionVector(:,63);
population.maximumHeight                          = decisionVector(:,64);
population.additionalMass                         = decisionVector(:,65);
population.terminationDistanceRatio               = decisionVector(:,66);

population.nodeIntervalDescent1                   = decisionVector(:,67);
population.nodeIntervalDescent2                   = decisionVector(:,68);
population.nodeIntervalDescent3                   = decisionVector(:,69);
population.nodeIntervalDescent4                   = decisionVector(:,70);
population.nodeIntervalDescent5                   = decisionVector(:,71);
population.nodeIntervalDescent6                   = decisionVector(:,72);
population.nodeIntervalDescent7                   = decisionVector(:,73);
population.nodeIntervalDescent8                   = decisionVector(:,74);
population.nodeIntervalDescent9                   = decisionVector(:,75);

population.angleOfAttackDescent1                  = decisionVector(:,76);
population.angleOfAttackDescent2                  = decisionVector(:,77);
population.angleOfAttackDescent3                  = decisionVector(:,78);
population.angleOfAttackDescent4                  = decisionVector(:,79);
population.angleOfAttackDescent5                  = decisionVector(:,80);
population.angleOfAttackDescent6                  = decisionVector(:,81);
population.angleOfAttackDescent7                  = decisionVector(:,82);
population.angleOfAttackDescent8                  = decisionVector(:,83);
population.angleOfAttackDescent9                  = decisionVector(:,84);
population.angleOfAttackDescent10                 = decisionVector(:,85);

population.bankAngleDescent1                      = decisionVector(:,86);
population.bankAngleDescent2                      = decisionVector(:,87);
population.bankAngleDescent3                      = decisionVector(:,88);
population.bankAngleDescent4                      = decisionVector(:,89);
population.bankAngleDescent5                      = decisionVector(:,90);
population.bankAngleDescent6                      = decisionVector(:,91);
population.bankAngleDescent7                      = decisionVector(:,92);
population.bankAngleDescent8                      = decisionVector(:,93);
population.bankAngleDescent9                      = decisionVector(:,94);
population.bankAngleDescent10                     = decisionVector(:,95);

population.thrustElevationAngleDescent1           = decisionVector(:,96);
population.thrustElevationAngleDescent2           = decisionVector(:,97);
population.thrustElevationAngleDescent3           = decisionVector(:,98);
population.thrustElevationAngleDescent4           = decisionVector(:,99);
population.thrustElevationAngleDescent5           = decisionVector(:,100);
population.thrustElevationAngleDescent6           = decisionVector(:,101);
population.thrustElevationAngleDescent7           = decisionVector(:,102);
population.thrustElevationAngleDescent8           = decisionVector(:,103);
population.thrustElevationAngleDescent9           = decisionVector(:,104);
population.thrustElevationAngleDescent10          = decisionVector(:,105);

population.thrustAzimuthAngleDescent1             = decisionVector(:,106);
population.thrustAzimuthAngleDescent2             = decisionVector(:,107);
population.thrustAzimuthAngleDescent3             = decisionVector(:,108);
population.thrustAzimuthAngleDescent4             = decisionVector(:,109);
population.thrustAzimuthAngleDescent5             = decisionVector(:,110);
population.thrustAzimuthAngleDescent6             = decisionVector(:,111);
population.thrustAzimuthAngleDescent7             = decisionVector(:,112);
population.thrustAzimuthAngleDescent8             = decisionVector(:,113);
population.thrustAzimuthAngleDescent9             = decisionVector(:,114);
population.thrustAzimuthAngleDescent10            = decisionVector(:,115);

population.throttleSettingDescent1                = decisionVector(:,116);
population.throttleSettingDescent2                = decisionVector(:,117);
population.throttleSettingDescent3                = decisionVector(:,118);
population.throttleSettingDescent4                = decisionVector(:,119);
population.throttleSettingDescent5                = decisionVector(:,120);
population.throttleSettingDescent6                = decisionVector(:,121);
population.throttleSettingDescent7                = decisionVector(:,122);
population.throttleSettingDescent8                = decisionVector(:,123);
population.throttleSettingDescent9                = decisionVector(:,124);
population.throttleSettingDescent10               = decisionVector(:,125);

population.finalVelocity                          = decisionVector(:,126);
population.skipSuppressionTriggerTime             = decisionVector(:,127);
population.maximumMechanicalLoad                  = decisionVector(:,128);


end