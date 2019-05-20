function Uc_c = findEnergyOfc_cInteraction(X,Y,Z,numberOfPoints)
%% constants
De = 6.325; % eV
S = 1.29;
Beta = 1.5; % 1/?
Re = 1.315; % ?
R1 = 1.7; % ?
R2 = 2.0; % ?
Delta = 0.80469;
a = 0.011304;
c = 19;
d = 2.5;

%%
for i = 1:numberOfPoints
    for j = 1:numberOfPoints
        r(i,j) = sqrt((X(i)-X(j))^2+(Y(i)-Y(j))^2+(Z(i)-Z(j))^2);   %find the distance between each atom
        %find f cut function
        if r(i,j)<=R1
            F(i,j) = 1;
        elseif r(i,j)>R2
            F(i,j) = 0;
        else
            F(i,j) = 1/2 + cos((r(i,j)-R1)/(R2-R1)*pi)/2;
        end
        Ur(i,j) = De/(S-1)*exp(-Beta*sqrt(2*S)*(r(i,j)-Re));    %find U(R) and U(A)
        Ua(i,j) = De*S/(S-1)*exp(-Beta*sqrt(2/S)*(r(i,j)-Re));
    end
end

Uc_c = 0;
summa=0;
for j = 1:numberOfPoints
    for i = 1:numberOfPoints
        for k = 1:numberOfPoints
            if j ~= i && j ~= k && i ~= k
                %find Gijk
                %costeta = ((X(i)-X(j))*(X(k)-X(j))+(Y(i)-Y(j))*(Y(k)-Y(j))+(Z(i)-Z(j))*(Z(k)-Z(j)))/(r(i,j)*r(j,k));
                summa = summa + F(j,k)*a*(1+c^2/d^2-c^2/(d^2+(1+ ...
                    ((X(i)-X(j))*(X(k)-X(j))+(Y(i)-Y(j))*(Y(k)-Y(j))+(Z(i)-Z(j))*(Z(k)-Z(j))) ...
                    /(r(i,j)*r(j,k)))^2));
            end
        end
        if j~=i
            Uc_c = Uc_c + F(i,j)*(Ur(i,j)-(1+summa)^-Delta*Ua(i,j));%final
        end
        summa=0;
    end
end
Uc_c = Uc_c/2;
end
