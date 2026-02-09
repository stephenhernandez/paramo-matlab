function RM = rotation_matrix(theta)
theta = theta*(pi/180.);
costheta = cos(theta);
sintheta = sin(theta);
RM = [costheta -sintheta; ...
    sintheta costheta];