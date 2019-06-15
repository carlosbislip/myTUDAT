function [ fitness ] = loadFitnessVectorFromFitnessCellArray( fitnessData )

fitnessVector = cell2mat(fitnessData(:,3:end));

fitness.nodeIntervalAscent1                    = fitnessVector(:,1);
fitness.nodeIntervalAscent2                    = fitnessVector(:,2);
fitness.nodeIntervalAscent3                    = fitnessVector(:,3);



end