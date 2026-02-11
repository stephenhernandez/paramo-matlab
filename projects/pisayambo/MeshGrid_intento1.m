nFilas = 61; filasX = NaN(length(X),nFilas); filasY = filasX;dz = 100; newDX = dz*cosd(35)/tand(60);
newDY = dz*sind(35)/tand(60);
for i = 1:length(x)
    X_ = X(i);
    Y_ = Y(i);
    for j = 1:nFilas
        filasX(i,j+1) = X_-newDX; filasY(i,j+1) = Y_ + newDY;
    end
end