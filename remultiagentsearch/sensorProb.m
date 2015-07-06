function prob = sensorProb(x,y,xlen,ylen,sigmaVal)
sigma = sigmaVal*[1 0; 0 1];
x_r = [x;y];
[ptx,pty] = meshgrid(1:xlen,1:ylen);
pt = [ptx(:)';pty(:)'];
pt = pt - x_r*ones(1,xlen*ylen);
tmp = pt'/sigma*pt;
tmp_diag = diag(tmp);
% prob = exp(-1/2*(t-x_r)'/sigma*(t-x_r));
prob = exp(-1/2*(reshape(tmp_diag,xlen,ylen))');
end