%% Test script to compute signed distance for given ego and target positions
clc
clear
close all

%% Load data and set paths
pathData = 'results/';
pathFigs = 'figs/';
addpath(genpath('src'));

if ispc
    pathFigs = '';
end

load([pathData, 'sample_data']);

saveFigs = 0;   %sum(strcmp(input('Save figures [y/n]? ','s'),{'y','Y'}));

% set path on Mac for using pdflatex
if ismac
    setenv('PATH', [getenv('PATH') ':/Library/TeX/texbin:/usr/local/bin']);
end

%% Linearized signed distance
egoPosAll   = [[20,0.5*roadParams.width,0]',[28,0.5*roadParams.width,0]'];
tarPosAll   = repmat([30,-0.5*roadParams.width,0]',1,2);
c           = vehParams.c;
dSafe       = 0.5;

for n=1:size(egoPosAll,2)
    figNo = figure(40+n);
    figName = [pathFigs, 'env_LinearSignedDist', num2str(n)];
    h = figureSettings(30,'normalized',[0.2 0.2 0.5 0.25]);
    hold on
    box on
    egoPos = egoPosAll(:,n);
    tarPos = tarPosAll(:,n);
    Obstacles{1} = [tarPos(1)-b, tarPos(2)-c; tarPos(1)-b, tarPos(2)+c; tarPos(1)+a, tarPos(2)+c; tarPos(1)+a, tarPos(2)-c];
    sMin = egoPos(1)-10;
    sMax = egoPos(1)+20;
    eYOff = 0.2;
    eYMin = -roadParams.width-eYOff;
    eYMax = roadParams.width+eYOff;
    makeRobotPoly   = @(x) orientedBoxToPolygon([x(1), x(2), vehParams.a, vehParams.b, 2*vehParams.c, x(3)]);
    [g, gradG]  = g_collisions(egoPos, dSafe, [3,1], makeRobotPoly, Obstacles);
    bU = [sMax; roadParams.width; 2*pi];
    bL = [sMin; -roadParams.width; -2*pi];
    ACons = [-gradG;...
        eye(3);...
        -eye(3)];
    bCons = [g - gradG*egoPos;...
        bU;...
        -bL];
    PCons = Polyhedron(ACons,bCons);
    PCons_xy = PCons.slice(3,egoPos(3));
    PCons_xy.minHRep();
    hSet = patch(PCons_xy.V(:,1),PCons_xy.V(:,2),0.9*[211 211 211]/255,'facealpha',0.2);
    % plot(PCons_xy.getFacet);
    % hSet = plot(PCons_xy,'alpha',0.2);
    span = sMin:sMax;
    roadFig(1)  = plot(span,zeros(length(span),1),'--k','LineWidth',2);
    roadFig(2)  = plot(span,-roadParams.width*ones(length(span),1),'k','LineWidth',2);
    roadFig(3)  = plot(span,roadParams.width*ones(length(span),1),'k','LineWidth',2);
    carFig      = PlotCar2(egoPos, vehParams, 0, [0.8 0.78 0.78]);
    targetFig   = PlotCar2(tarPos, vehParams, 0, [0.8 0.1 0.1]);
    textFontSize = 30;
    text((egoPos(1)-sMin-1.2*vehParams.length)/(sMax - sMin),(egoPos(2)-eYMin+0)/(eYMax-eYMin),'$E$','units','normalized','fontsize',textFontSize)
    text((tarPos(1)-sMin-1.2*vehParams.length)/(sMax - sMin),(tarPos(2)-eYMin+0)/(eYMax-eYMin),'$V$','units','normalized','fontsize',textFontSize)
    grid off
    xlabel('X [m]');
    ylabel('Y [m]');
    % hAnnotation = get(hSet,'Annotation');
    % hLegendEntry = get(hAnnotation,'LegendInformation');
    % set(hLegendEntry,'IconDisplayStyle','on')
    % legend(hSet,{'Unsafe region'},'location','southwest');
    axis([sMin sMax eYMin eYMax])
    if saveFigs
        mlf2pdf(figNo,figName)
    end
end
