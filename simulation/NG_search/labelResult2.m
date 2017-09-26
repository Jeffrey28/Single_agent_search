function labeledRes = labelResult2(res,fctname,gmm_num,N) %this,
% label each row and column of input res with the name of
% corresponding variables

% this one is for ipopt variables

% s = [z(:),u(:),x(:),xpred(:),P(:),P_pred(:),K(:),Pinverse(:),auxt(:),auxm(:)]

%%% define header of variables
% header consisting of each variable's name
varHeader = {};
% z
for ii = 1:N+1
    varHeader = [varHeader,{sprintf('z(1,%d)',ii-1)},{sprintf('z(2,%d)',ii-1)}...
        ,{sprintf('z(3,%d)',ii-1)},{sprintf('z(4,%d)',ii-1)}];
end
% u
for ii = 1:N
    varHeader = [varHeader,{sprintf('u(1,%d)',ii)},{sprintf('u(2,%d)',ii)}];
end
% x
for ii = 1:N+1
    for jj = 1:gmm_num
        varHeader = [varHeader,{sprintf('x(1,%d,%d)',jj,ii-1)},{sprintf('x(2,%d,%d)',jj,ii-1)}];
    end
end
%x_pred
for ii = 1:N
    for jj = 1:gmm_num
        varHeader = [varHeader,{sprintf('x_pred(1,%d,%d)',jj,ii)},{sprintf('x_pred(2,%d,%d)',jj,ii)}];
    end
end
% P
for ii = 1:N+1
    for jj = 1:gmm_num
        varHeader = [varHeader,{sprintf('P(1,1,%d,%d)',jj,ii-1)},{sprintf('P(2,1,%d,%d)',jj,ii-1)},...
            {sprintf('P(1,2,%d,%d)',jj,ii-1)},{sprintf('P(2,2,%d,%d)',jj,ii-1)}];
    end
end
% P_pred
for ii = 1:N
    for jj = 1:gmm_num
        varHeader = [varHeader,{sprintf('P_pred(1,1,%d,%d)',jj,ii)},{sprintf('P_pred(2,1,%d,%d)',jj,ii)},...
            {sprintf('P_pred(1,2,%d,%d)',jj,ii)},{sprintf('P_pred(2,2,%d,%d)',jj,ii)}];
    end
end
% K
for ii = 1:N
    for jj = 1:gmm_num
        varHeader = [varHeader,{sprintf('K(1,1,%d,%d)',jj,ii)},{sprintf('K(2,1,%d,%d)',jj,ii)},...
            {sprintf('K(1,2,%d,%d)',jj,ii)},{sprintf('K(2,2,%d,%d)',jj,ii)}];
    end
end
% Pinverse
for ii = 1:N+1
    for jj = 1:gmm_num
        varHeader = [varHeader,{sprintf('Pinv(1,1,%d,%d)',jj,ii-1)},{sprintf('Pinv(2,1,%d,%d)',jj,ii-1)},...
            {sprintf('Pinv(1,2,%d,%d)',jj,ii-1)},{sprintf('Pinv(2,2,%d,%d)',jj,ii-1)}];
    end
end

% auxt
for ii = 1:N
    for jj = 1:gmm_num
        varHeader = [varHeader,{sprintf('auxt(%d,%d)',jj,ii)}];
    end
end

% auxm
for ii = 1:N
    for jj = 1:gmm_num
        for ll = 1:gmm_num
            varHeader = [varHeader,{sprintf('auxm(%d,%d,%d)',jj,ll,ii)}];
        end
    end
end

% gamma variable
for ii = 1:N
    for jj = 1:gmm_num
        varHeader = [varHeader,{sprintf('gamvar(%d,%d)',jj,ii)}];
    end
end

%%% turn result into cell
cellRes = num2cell(res);
constrHeader = {};

%%% label result
switch fctname
    case 's'
        % decision variable
        % row: variables in s (note: s is a column vector)
        labeledRes = [varHeader',cellRes];
        
    case 'fgrad'
        % fgrad
        % column: variables in s (note: fgrad is a row vector)
        labeledRes = [varHeader;cellRes];
        
    case 'hjac'
        % hjac
        % jacobian of h
        % column: variables in s
        % row: nonlinear equality constraints
        
        %
        for ii = 1:N
            % z(:,ii+1) - (z(:,ii)+ [z(4,ii)*cos(z(3,ii)); z(4,ii)*sin(z(3,ii));...
            % u(:,ii)]*dt);
            constrHeader = [constrHeader;{sprintf('z(1,%d)-z(1,%d)',ii,ii-1)};...
                {sprintf('z(2,%d)-z(2,%d)',ii,ii-1)};{sprintf('z(3,%d)-z(3,%d)',ii,ii-1)};...
                {sprintf('z(4,%d)-z(4,%d)',ii,ii-1)}];
        end
        
        for ii = 1:N
            % K = P_k+1|k*C_k+1'(C_k+1*P_k+1|k*C_k+1'+R)^-1
            constrHeader = [constrHeader;{sprintf('K(1,1,%d)',ii)};...
                {sprintf('K(2,1,%d)=0',ii)};{sprintf('K(1,2,%d)',ii)};...
                {sprintf('K(2,2,%d)',ii)}]; % K(2,1,%d) is actually not a constraint since I use triu. this one is always 0.
        end
        
        for ii = 1:N
            % P_k+1|k+1 = P_k+1|k-gamma*K*C*P_k+1|k
            constrHeader = [constrHeader;{sprintf('P(1,1,%d)',ii)};...
                {sprintf('P(2,1,%d)',ii)};{sprintf('P(1,2,%d)',ii)};...
                {sprintf('P(2,2,%d)',ii)}];
        end
        
        labeledRes = [{'nothing'},varHeader;constrHeader,cellRes];
        
    case 'hlin'
        % hlin
        % row: value for linear equality constraints
        
        % z(:,1) == this.state;
        % x(:,1) == this.est_pos(:);
        % for jj = 1:this.gmm_num
        %  triu(P(:,:,jj,1)) == triu(this.P{jj});
        % end
        constrHeader = [constrHeader;{'z(1,1)-init'};{'z(2,1)-init'};...
            {'z(3,1)-init'};{'z(4,1)-init'};{'x(1,1)-init'};{'x(2,1)-init'};...
            {'P(1,1,1)-init'};{'P(2,1,1)-init'};{'P(1,2,1)-init'};{'P(2,2,1)-init'}];
        % x_k+1|k = f(x_k)
        for ii = 1:N
            constrHeader = [constrHeader;{sprintf('x_pred(1,%d)-f(x(1,%d))',ii,ii)};...
                {sprintf('x_pred(2,%d)-f(x(2,%d))',ii,ii)}];
        end
        % P_k+1|k=A*P_k*A+Q
        for ii = 1:N
            constrHeader = [constrHeader;{sprintf('P_pred(1,1,%d)-AP(1,1,%d))A',ii,ii)};...
                {sprintf('P_pred(2,1,%d)-AP(2,1,%d))A=0',ii,ii)};{sprintf('P_pred(1,2,%d)-AP(1,2,%d))A',ii,ii)};...
                {sprintf('P_pred(2,2,%d)-AP(2,2,%d))A',ii,ii)}];
        end
        % x_k+1|k+1 = x_k+1|k
        for ii = 1:N
            constrHeader = [constrHeader;{sprintf('x(1,%d)-x_pred(1,%d)',ii+1,ii)};...
                {sprintf('x(2,%d)-x_pred(2,%d)',ii+1,ii)}];
        end
        
        % symmetric constraint of P and Pinv   
        for ii = 1:N+1
            for jj = 1:gmm_num
                constrHeader = [constrHeader;{sprintf('P(1,2,%d,%d)-P(2,1,%d,%d)',jj,ii)};...
                    {sprintf('Pinv(1,2,%d,%d)-Pinv(2,1,%d,%d)',jj,ii)}];
            end
        end
        
        labeledRes = [constrHeader,cellRes];
        
    case 'glin'
        % glin
        % row: value for linear inequality constraints
        % this.w_lb <= u(1,:) <= this.w_ub;
        % this.a_lb <= u(2,:) <= this.a_ub;
        % this.v_lb <= z(4,:) <= this.v_ub;
        for ii = 1:N
            constrHeader = [constrHeader;{sprintf('u(1,%d)-w_ub',ii)}];
        end
        for ii = 1:N
            constrHeader = [constrHeader;{sprintf('u(2,%d)-a_ub',ii)}];
        end
        for ii = 1:N+1
            constrHeader = [constrHeader;{sprintf('z(4,%d)-v_ub',ii)}];
        end
        for ii = 1:N
            constrHeader = [constrHeader;{sprintf('w_lb-u(1,%d)',ii)}];
        end
        for ii = 1:N
            constrHeader = [constrHeader;{sprintf('a_lb-u(2,%d)',ii)}];
        end
        for ii = 1:N+1
            constrHeader = [constrHeader;{sprintf('v_lb-z(4,%d)',ii)}];
        end
        
        % [fld.fld_cor(1);fld.fld_cor(3)]<=z(1:2,ii+1)<=[fld.fld_cor(2);fld.fld_cor(4)];
        for ii = 1:N
            constrHeader = [constrHeader;{sprintf('z(1,%d)-fld_cor(2)',ii+1)};...
                {sprintf('z(2,%d)-fld_cor(4)',ii+1)};{sprintf('fld_cor(1)-z(1,%d)',ii+1)};...
                {sprintf('fld_cor(3)-z(2,%d)',ii+1)}];
        end
        
        % P, Ppred has positive diagonal element
        for ii = 1:N
            for jj = 1:gmm_num                
                constrHeader = [constrHeader;{sprintf('-P(1,1,%d,%d)',jj,ii+1)};...
                    {sprintf('-P(2,2,%d,%d)',jj,ii+1)};...
                    {sprintf('-Ppred(1,1,%d,%d)',jj,ii+1)};...
                    {sprintf('-Ppred(2,2,%d,%d)',jj,ii+1)}];
            end
        end
        
        % auxiliary variable t,m are nonnegative
        % t
        for ii = 1:N
            for jj = 1:gmm_num                
                constrHeader = [constrHeader;{sprintf('-auxt(%d,%d)',jj,ii)}];
            end
        end
        
        % m
        for ii = 1:N
            for jj = 1:gmm_num
                for ll = 1:gmm_num
                    constrHeader = [constrHeader;{sprintf('-auxm(%d,%d,%d)',jj,ll,ii)}];
                end
            end
        end
        
        labeledRes = [constrHeader,cellRes];
        
    case 'h'
        % nonlinear equality constraint
        % h is a column vector
        % row: value for nonlinear equality constraints
        
        for ii = 1:N
            % z(:,ii+1) - (z(:,ii)+ [z(4,ii)*cos(z(3,ii)); z(4,ii)*sin(z(3,ii));...
            % u(:,ii)]*dt);
            constrHeader = [constrHeader;{sprintf('z(1,%d)-z(1,%d)',ii,ii-1)};...
                {sprintf('z(2,%d)-z(2,%d)',ii,ii-1)};{sprintf('z(3,%d)-z(3,%d)',ii,ii-1)};...
                {sprintf('z(4,%d)-z(4,%d)',ii,ii-1)}];
        end
        
        for ii = 1:N
            % K = P_k+1|k*C_k+1'(C_k+1*P_k+1|k*C_k+1'+R)^-1
            constrHeader = [constrHeader;{sprintf('K(1,1,%d)',ii)};...
                {sprintf('K(2,1,%d)=0',ii)};{sprintf('K(1,2,%d)',ii)};...
                {sprintf('K(2,2,%d)',ii)}];
        end
        
        for ii = 1:N
            % P_k+1|k+1 = P_k+1|k-gamma*K*C*P_k+1|k
            constrHeader = [constrHeader;{sprintf('P(1,1,%d)',ii)};...
                {sprintf('P(2,1,%d)=0',ii)};{sprintf('P(1,2,%d)',ii)};...
                {sprintf('P(2,2,%d)',ii)}];
        end
        
        % triu(P(:,:,jj,ii)*Pinverse(:,:,jj,ii))-eye(2)
        for ii = 1:N+1
            for jj = 1:gmm_num
                constrHeader = [constrHeader;{sprintf('[P(:,:,%d,%d)*Pinv(:,:,%d,%d)-1](1,1)',jj,ii)};...
                    {sprintf('[P(:,:,%d,%d)-Pinv(:,:,%d,%d)](2,1)',jj,ii)};...
                    {sprintf('[P(:,:,%d,%d)-Pinv(:,:,%d,%d)](1,2)',jj,ii)};...
                    {sprintf('[P(:,:,%d,%d)-Pinv(:,:,%d,%d)-1](2,2)',jj,ii)}];
            end
        end
        
        % gam_var(jj,ii)*gamDen(jj,ii) = 1
        for ii = 1:N
            for jj = 1:gmm_num
                constrHeader = [constrHeader;{sprintf('gam_var(%d,%d)',jj,ii)}];
            end
        end
        
        labeledRes = [constrHeader,cellRes];
        
    case 'g'
        % (x(jj)-x(ll))'*Pinv(ll)*(x(jj)-x(ll))+2log(2pi/wt(ll))+log(P(ll))+2log(auxm(jj,ll))
        % <= 0
        for ii = 1:N
            for jj = 1:gmm_num
                for ll = 1:gmm_num
                    constrHeader = [constrHeader;{sprintf('mrelation(%d,%d,%d)',jj,ll,ii)}];
                end
            end
        end
        
        % exp(-auxt(jj)/wt(jj))-sum(auxm(jj,:)) <= 0
        for ii = 1:N
            for jj = 1:gmm_num
                constrHeader = [constrHeader;{sprintf('tmrelation(%d,%d)',jj,ii)}];                
            end
        end
        
        labeledRes = [constrHeader,cellRes];
end
end