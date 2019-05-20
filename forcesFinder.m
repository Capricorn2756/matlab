function Forces = forcesFinder(X,Y,Z,numberOfPoints,Disp)
tempX = X;
tempY = Y;
tempZ = Z;
Forces = zeros(numberOfPoints,3);
Uc_c = findEnergyOfc_cInteraction(X,Y,Z,numberOfPoints);
for i = 1:numberOfPoints
    tempX(i) = tempX(i) + Disp;
    Forces(i,1) = -(findEnergyOfc_cInteraction(tempX,Y,Z,numberOfPoints) - Uc_c)/Disp;
    tempX(i) = X(i); 
    
    tempY(i) = tempY(i) + Disp;
    Forces(i,2) = -(findEnergyOfc_cInteraction(X,tempY,Z,numberOfPoints) - Uc_c)/Disp;
    tempY(i) = Y(i);
    
    tempZ(i) = tempZ(i) + Disp;
    Forces(i,3) = -(findEnergyOfc_cInteraction(X,Y,tempZ,numberOfPoints) - Uc_c)/Disp;
    tempZ(i) = Z(i);
end