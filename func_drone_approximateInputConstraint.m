function [points,A,b] = func_drone_approximateInputConstraint(s1,s2,Tmax,epmax)
g=9.81;
R=Tmax*sin(epmax);
v3 =Tmax*cos(epmax)-g ;
alpha=linspace(0,2*pi,s1);
r = linspace(0,R,s2);
points=[0,0,-g];
for i = 1:s1
    for j = 1:s2
        alpha_i = alpha(i);
        r_j = r(j);
        points = [points;
            [r_j*cos(alpha_i),r_j*sin(alpha_i),sqrt(Tmax^2-r_j^2)-g]];
        points = [points; [R*cos(alpha_i),R*sin(alpha_i),v3]];
        
    end
end

H=polytope(points');
A = H.A;
b = H.b;
end