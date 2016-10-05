% draw entropy
clc
clear
folder_path = ('.\sim_res');
figure('WindowStyle','normal','Position', [500 50 800 700])
for ii = 1:3
    data_name = sprintf('test%d_offset.mat',ii);
    file_name = fullfile(folder_path,data_name);
    load(file_name)
    len = size(rbt.est_pos_hist,2); 
    P_trace = zeros(len,1);
    P_determinant = zeros(len,1);
    for jj = 1:len
        P_trace(jj) = rbt.P_hist(1,2*jj-1)+rbt.P_hist(2,2*jj);
        P_determinant(jj) = det(rbt.P_hist(:,2*jj-1:2*jj));
    end
%     figure(1)
    hold on
    plot(P_trace,'-.','LineWidth',2)
    
%     figure(2)     
%     semilogy(1:len,P_determinant,'-.','LineWidth',2)
%     hold on
end
% axis equal

[~,hdl4] = legend('test1','test2','test3','test4');
textobj = findobj(hdl4,'type','text');
set(textobj,'fontsize',20);

box on
set(gca,'fontsize',30);
title('Trace of Covariance Matrix P','FontSize',25)
xlabel('Time Step','FontSize',30)
ylabel('Trace','FontSize',30)
ylim([0,21])
xlim([1,50])

% figure(2)
% legend('test1','test2','test3','test4')
% box on
% title('Determinant of Covariance Matrix P')
% xlabel('Time Step')
% ylabel('Determinant')
% ylim([0,110])