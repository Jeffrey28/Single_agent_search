function particles = genParticles(agent,w,mu,sigma,n_data)
gmm_obj = gmdistribution(mu',sigma,w');
% generate random data from gmm
particles = random(gmm_obj,n_data);
% make the particles a 2-by-n_data array
if size(particles,2) == 2
    particles = particles';
end