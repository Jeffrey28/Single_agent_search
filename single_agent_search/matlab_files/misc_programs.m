%% used for saving miscellaneous code snippets 

% check if the point-to-line distance constraint is enforced or not
v1 = [76.11;8.782];
v2 = [68.15;14.58];
x0 = 73.333;
y0 = 13.333;
r = 3.333;
a = v2(2)-v1(2);
b = v1(1)-v2(1);
c = v1(2)*v2(1)-v2(2)*v1(1);

val = (x0^2-r^2)*a^2+(y0^2-r^2)*b^2+c^2+2*a*b*x0*y0+2*a*c*x0+2*b*c*y0;
