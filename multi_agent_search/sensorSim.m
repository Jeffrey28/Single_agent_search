function z = sensorSim(x,y,tx,ty,sigmaVal)
x_r = [x;y];
t = [tx;ty];
sigma = sigmaVal*[1 0; 0 1];
prob = exp(-1/2*(t-x_r)'/sigma*(t-x_r));
z = (rand(1,1) < prob);
end