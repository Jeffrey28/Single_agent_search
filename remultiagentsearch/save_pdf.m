function save_pdf(h,filename)
% plot(1:10);
set(h,'Units','Inches');
pos = get(h,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,filename,'-dpdf','-r0')