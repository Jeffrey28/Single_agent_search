%% Optimization problems
% modify this test file based on CS287

% problem0 is just a QP. Your optimizer should jump to the solution
% right away if the initial trust region size is big enough

clear % clear global variables
close all

simSetup;

zerofunc = @(x) 0;
neginffunc = @(x) -1e5;
hinge = @(x) sum(max(x,0));
abssum = @(x) sum(abs(x));

n = 2;

penalty_coeff = cfg.initial_penalty_coeff; % Coefficient of l1 penalties
trust_box_size = cfg.initial_trust_box_size; % The trust region will be a box around the current iterate x.

problem_default = struct();
problem_default.f = zerofunc;
problem_default.g = neginffunc;
problem_default.h = zerofunc;
problem_default.Q = zeros(n,n);
problem_default.q = zeros(1,n);
problem_default.A_ineq = zeros(0,n);
problem_default.b_ineq = zeros(0,1);
problem_default.A_eq = zeros(0,n);
problem_default.b_eq = zeros(0,1);
problem_default.glin = neginffunc;
problem_default.hlin = zerofunc;
problem_default.hinge = hinge;
problem_default.abssum = abssum;
problem_default.penalty_coeff = penalty_coeff;
problem_default.trust_box_size = trust_box_size;
problem_default.snum = 0;


% f,g (QP)
problem0 = problem_default;
problem0.x0 = [1,1]';
problem0.f = @(x) x(1).^2 + x(2).^2; % Objective
problem0.g = @(x) 3 - x(1) - x(2);   % Inequality constraint (should be <= 0)
problem0.xtrue = [1.5, 1.5]';                 % Optimal solution

% f,g (QP)
problem1 = problem_default;
problem1.x0 = [-2,1]';
problem1.f = @(x) (x(2) - x(1).^2).^2 + (1 - x(1)).^2;
problem1.g = @(x) -1.5 - x(2);
problem1.xtrue = [1,1]';

% f,g (QP)
problem2 = problem_default;
problem2.x0 = [10,1]';
problem2.f = @(x) x(2) + 1e-5 + (x(2) - x(1)).^2;
problem2.g = @(x) -x(2);
problem2.xtrue = [0,0]';

% f (QP), h
problem3 = problem_default;
problem3.x0 = [10,1]';
problem3.f = @(x) (1 - x(1)).^2;
problem3.h = @(x) 10*(x(2) - x(1).^2);
problem3.xtrue = [1,1]';

% f,h
problem4 = problem_default;
problem4.x0 = [2,2]';
problem4.f = @(x) log(1+x(1).^2) - x(2);
problem4.h = @(x) (1+x(1).^2).^2 + x(2).^2 - 4;
problem4.xtrue = [0,sqrt(3)]';

% This is an LP
% q'x, glin, q
angles = (1:6)'*2*pi/6;
problem5 = problem_default;
problem5.x0 = [0,0]';
% problem5.A_ineq = [cos(angles), sin(angles)];
% problem5.b_ineq = ones(length(angles),1);
problem5.glin = @(x) [cos(angles), sin(angles)]*x-ones(length(angles),1);
problem5.q = -[cos(pi/6), sin(pi/6)];
problem5.xtrue = [1, tan(pi/6)]';

% Same LP constraints, but poorly scaled and implemented as nonlinear constraints
problem6 = problem_default;
problem6.x0 = [0,0]';
problem6.Q = .1*eye(2);
problem6.q = -[cos(pi/6), sin(pi/6)];
problem6.g = @(x) .01*problem5.glin(x) %(problem5.A_ineq*x - problem5.b_ineq);
problem6.xtrue = [1, tan(pi/6)]';
% Note that if you set problem6.Q = zeros(2) instead, our SQP algorithm will
% diverge. That's because the penalized problem can be unbounded, though
% the constrained problem is bounded on the feasible set. More sophisticated
% optimization algorithms adjust the penalty coefficient to ensure that
% each step decreases the constraint violations. 

% f, g(linear),h(linear)
problem7 = problem_default;
problem7.x0 = [0,0]';
problem7.f = @(x) x(1).^4 + x(2).^4;
problem7.g = @(x) 3 - x(1) - x(2);
problem7.h = @(x) x(1) - 2*x(2);
problem7.xtrue = [2 1]';

% Kowabunga!
% x'Qx, g
problem8 = problem_default;
problem8.x0 = [5,5]';
problem8.g = @(x) [((x(1)-0).^2 + (x(2)-0).^2 - 4); -((x(1)-1).^2 + (x(2)-1).^2 - .25); -((x(1)+1).^2 + (x(2)-1).^2 - .25); -((x(1)-0).^2 + 7*(x(2)+1-x(1).^2/2).^2 - .8)];
problem8.Q = eye(2);
problem8.xtrue = [0,0]';



problems = {problem7};%{problem0, problem1, problem2, problem3, problem4, problem5, problem6, problem7, problem8};

for i_problem = 1:length(problems)
    figure(1);
    clf
    hold on;    
    legenditems = {};

    
    problem = problems{i_problem};
    
    % plots level sets of f (red), g (green) and h (blue)
    [xplot, yplot]= meshgrid(-10:.1:10, -10:.1:10);
    fplot=zeros(size(xplot));
    gplot=zeros(size(yplot));
    hplot=zeros(size(xplot));
    for i = 1:size(xplot,1)
        for j=1:size(xplot,2)
            x = [xplot(i,j);yplot(i,j)];
            fplot(i,j) = max([problem.f(x) + x'*problem.Q*x + problem.q*x; problem.A_eq*x - problem.b_eq]);
            gplot(i,j) = max([problem.g(x); problem.A_ineq*x - problem.b_ineq]);
            hplot(i,j) = problem.h(x);
        end
    end

    if ~isequal(problem.g, neginffunc) || size(problem.A_ineq,1) > 0
        contourf(xplot,yplot,-gplot,[0 0],'Color','g');
        colormap(repmat([0,1,0],64,1));
        legenditems{end+1} = 'feasible set for ineq constraints';
    end

    if ~isequal(problem.f, zerofunc) || numel(problem.Q)>0 || numel(problem.q)>0
        contour(xplot,yplot,fplot,'Color','r');        
        legenditems{end+1} = 'objective';
    end

    
    if ~isequal(problem.h, zerofunc) || size(problem.A_eq,1) > 0
        contour(xplot,yplot,hplot,[0 0],'Color','b');
        legenditems{end+1} = 'feasible set for eq constraints';
    end
    plot_path_2d(); % Clear persistent variables in that function
    axis('equal');

    
    legend(legenditems{:});



    % Set parameters of the optimizer
%     cfg = struct();
    cfg.min_approx_improve = 1e-8;
    cfg.min_trust_box_size = 1e-5;
    cfg.callback = @plot_path_2d;
    cfg.full_hessian = true;
    
    % x, Q, q, f, A_ineq, b_ineq, A_eq, b_eq, glin, hlin, g, h, hinge, abssum, cfg, penalty_coeff, trust_box_size, snum, fld
    
    % Run the optimizer
    x = rbt.minimize_merit_function(problem.x0, problem.Q, problem.q, problem.f, problem.A_ineq, problem.b_ineq, ...
        problem.A_eq, problem.b_eq, problem.glin, problem.hlin, problem.g, problem.h, hinge, abssum, cfg, penalty_coeff, trust_box_size, [], fld);    
    fprintf('computed: %s. true %s. max error %.3e\n',mat2str(x,3), mat2str(problem.xtrue,3), max(abs(x-problem.xtrue)));
    
    saveas(gcf, sprintf('problem%i_path.png',i_problem));
    

    disp('press enter to continue to the next problem')
    pause();
end