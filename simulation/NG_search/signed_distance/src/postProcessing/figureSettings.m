function f=figureSettings(varargin)

% h=axes;
f = gcf;

set(gcf,'DefaultTextInterpreter', 'latex')
set(gcf,'DefaultTextFontName', 'Times');
set(gcf,'DefaultTextUnits','normalized');
set(gcf,'DefaultAxesFontName', 'Times');
set(gcf,'DefaultAxesLineWidth',2);
set(gcf,'DefaultLineLineWidth',2);
set(gcf,'paperPositionMode','auto');
set(gcf,'DefaultLegendInterpreter','latex');
% set(h,'DefaultTextInterpreter', 'latex')

if nargin > 0
    fontSize = varargin{1};
    set(gcf,'DefaultTextFontSize',fontSize);
    set(gcf,'DefaultAxesFontSize',fontSize);
%     set(h,'linewidth',2,'fontname','Times','fontSize',fontSize,'box','on');
%     hold(h,'on');
end

if nargin > 1
    set(gcf,'units',varargin{2},'position',varargin{3});
%     set(gcf,'PaperUnits',varargin{2},'PaperPosition',varargin{3});
end

if nargin > 3
    set(gcf,'DefaultLineLineWidth',varargin{4});
    set(gcf,'DefaultAxesLineWidth',varargin{4});
end
