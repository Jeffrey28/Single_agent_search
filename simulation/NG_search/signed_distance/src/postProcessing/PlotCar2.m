function handle = PlotCar2(zk, vehicleParam, axisoff, color, plotSensor, sensorColor, sensorParam)

% Plot a car as defined by input params
%   zk              :   Vehicle state: zk = [xpos, ypos, yaw, velocity, yaw_rate, acceleration]
%   vehicleParam    :   Parameters describing the vehicle
%   axisoff         :   Set to one for not adjusting axis in this function
%   color           :   A string with car color values
%
%   R   = [cos(psi), -sin(psi);sin(psi), cos(psi)];]
%   psi = Heading Angle

if nargin < 1
    offset = [0 0];
    R = eye(2);
    axisoff = 0;
    color = 'b';
    plotSensor = 0;
elseif nargin < 2
    R = eye(2);
    axisoff = 0;
    color = 'b';
    plotSensor = 0;
elseif nargin < 3
    axisoff = 0;
    color = 'b';
    plotSensor = 0;
elseif nargin < 4
    color = 'b';
    plotSensor = 0;
elseif nargin < 5
    plotSensor = 0;
end

   
%% Create car edges to draw
y   = linspace(vehicleParam.width/2,(-vehicleParam.width)/2,50);
x   = linspace(-vehicleParam.length/2,vehicleParam.length/2,50);

xf  = sqrt(vehicleParam.frontRadius^2 - (y).^2);
xf  = xf - xf(1) + vehicleParam.length/2;
xr  = -sqrt(vehicleParam.rearRadius^2 - (y).^2);
xr 	= xr - xr(1) - vehicleParam.length/2;

yl  = sqrt(vehicleParam.latRadius^2 - (x).^2);
yl  = yl - yl(1) + vehicleParam.width/2;

yr  = -yl;

%% Create tires
tireLength = 0.5;
tireWidth = 0.25;
frWheel = [vehicleParam.length/2-vehicleParam.distFrontWheel-tireLength/2, -(vehicleParam.width/2-tireWidth);...
           vehicleParam.length/2-vehicleParam.distFrontWheel-tireLength/2, -(vehicleParam.width/2);...
           vehicleParam.length/2-vehicleParam.distFrontWheel+tireLength/2, -(vehicleParam.width/2);...
           vehicleParam.length/2-vehicleParam.distFrontWheel+tireLength/2, -(vehicleParam.width/2-tireWidth)];
flWheel = [vehicleParam.length/2-vehicleParam.distFrontWheel-tireLength/2, (vehicleParam.width/2-tireWidth);...
           vehicleParam.length/2-vehicleParam.distFrontWheel-tireLength/2, (vehicleParam.width/2);...
           vehicleParam.length/2-vehicleParam.distFrontWheel+tireLength/2, (vehicleParam.width/2);...
           vehicleParam.length/2-vehicleParam.distFrontWheel+tireLength/2, (vehicleParam.width/2-tireWidth)];
rrWheel = [vehicleParam.length/2-vehicleParam.distFrontWheel-vehicleParam.wheelBase-tireLength/2, -(vehicleParam.width/2-tireWidth);...
           vehicleParam.length/2-vehicleParam.distFrontWheel-vehicleParam.wheelBase-tireLength/2, -(vehicleParam.width/2);...
           vehicleParam.length/2-vehicleParam.distFrontWheel-vehicleParam.wheelBase+tireLength/2, -(vehicleParam.width/2);...
           vehicleParam.length/2-vehicleParam.distFrontWheel-vehicleParam.wheelBase+tireLength/2, -(vehicleParam.width/2-tireWidth)];
rlWheel = [vehicleParam.length/2-vehicleParam.distFrontWheel-vehicleParam.wheelBase-tireLength/2, (vehicleParam.width/2-tireWidth);...
           vehicleParam.length/2-vehicleParam.distFrontWheel-vehicleParam.wheelBase-tireLength/2, (vehicleParam.width/2);...
           vehicleParam.length/2-vehicleParam.distFrontWheel-vehicleParam.wheelBase+tireLength/2, (vehicleParam.width/2);...
           vehicleParam.length/2-vehicleParam.distFrontWheel-vehicleParam.wheelBase+tireLength/2, (vehicleParam.width/2-tireWidth)];


%% Create Sensors
if plotSensor
    N = 100; %Number of points in sensor sectors
    sensor = zeros(2,N+1,length(sensorParam));

    for k = 1:length(sensorParam);
        theta = linspace(-sensorParam(k).openAngle/2 + sensorParam(k).alignment, sensorParam(k).openAngle/2 + sensorParam(k).alignment, N)';
        sensor(:,1,k) = [sensorParam(k).x, sensorParam(k).y]';
        sensor(:,2:end,k) = [sensorParam(k).maxRange*cos(theta), sensorParam(k).maxRange*sin(theta)]';
    end
end

%% Create outline
outLineOffset = 0.15;
outLineLeft = [x' (yl-outLineOffset)']';
outLineRight = [x' (yr+outLineOffset)']';

%% Create windows
frontWindowStart   = vehicleParam.length/2 - vehicleParam.length/4;
frontWindowEnd     = vehicleParam.length/2 - vehicleParam.length/2;

rearWindowStart   = -vehicleParam.length/2 + vehicleParam.length/8.5;
rearWindowEnd     = -vehicleParam.length/2 + vehicleParam.length/4;

frontWindowFront = [xf'-frontWindowStart, y(end:-1:1)'];
frontWindowFront = frontWindowFront(find(frontWindowFront(:,2) <= outLineLeft(2,min(find(outLineLeft(1,:)>=frontWindowStart))) & ...
                    frontWindowFront(:,2) >= outLineRight(2,min(find(outLineRight(1,:)>=frontWindowStart)))),:);
frontWindow = [frontWindowFront;
            frontWindowEnd, vehicleParam.width/2-0.3*vehicleParam.width/2;
            frontWindowEnd, -(vehicleParam.width/2-0.3*vehicleParam.width/2)];


rearWindowFront = [xr'+ vehicleParam.length/8.5, y(end:-1:1)'];
rearWindowFront = rearWindowFront(find(rearWindowFront(:,2) <= outLineLeft(2,min(find(outLineLeft(1,:)>=rearWindowStart))) & ...
                    rearWindowFront(:,2) >= outLineRight(2,min(find(outLineRight(1,:)>=rearWindowStart)))),:);

rearWindow = [rearWindowFront;
            rearWindowEnd, (vehicleParam.width/2-0.3*vehicleParam.width/2);
            rearWindowEnd, -(vehicleParam.width/2-0.3*vehicleParam.width/2)];

longDistance = 0.2;
latDistance = 0;
sideWindowRight = [outLineRight(:,find(outLineRight(1,:) <= frontWindowStart-longDistance & outLineRight(1,:) >= rearWindowStart+longDistance+0.1))';
                   frontWindowEnd-longDistance-0.15, -(vehicleParam.width/2-0.3*vehicleParam.width/2+latDistance);
                   rearWindowEnd+longDistance, -(vehicleParam.width/2-0.3*vehicleParam.width/2+latDistance)];
sideWindowLeft = [outLineLeft(:,find(outLineLeft(1,:) <= frontWindowStart-longDistance & outLineRight(1,:) >= rearWindowStart+longDistance+0.1))';
                   frontWindowEnd-longDistance-0.15, (vehicleParam.width/2-0.3*vehicleParam.width/2+latDistance);
                   rearWindowEnd+longDistance, (vehicleParam.width/2-0.3*vehicleParam.width/2+latDistance)];
               
        
%% Rotate
psi = zk(3);
R   = [cos(psi), -sin(psi);sin(psi), cos(psi)];

xfy = R * [xf' y']';
xry = R * [xr' y']';
xyl = R * [x' yl']';
xyr = R * [x' yr']';

outLineOffset = 0.15;
outLineLeft = R * [x' (yl-outLineOffset)']';
outLineRight = R * [x' (yr+outLineOffset)']';

frWheel = (R * frWheel')';
flWheel = (R * flWheel')';
rrWheel = (R * rrWheel')';
rlWheel = (R * rlWheel')';

frontWindow = (R * frontWindow')';
rearWindow = (R * rearWindow')';
sideWindowRight = (R * sideWindowRight')';
sideWindowLeft = (R * sideWindowLeft')';

if plotSensor
     for k = 1:length(sensorParam);
        sensor(:,:,k) = R*sensor(:,:,k);
     end
end

%% Translate
offset = [zk(1) zk(2)];

x1=xfy(1,:)+offset(1);
x2=xry(1,:)+offset(1);
x3=xyl(1,:)+offset(1);
x4=xyr(1,:)+offset(1);
y1=xfy(2,:)+offset(2);
y2=xry(2,:)+offset(2);
y3=xyl(2,:)+offset(2);
y4=xyr(2,:)+offset(2);

%% Heading
arrow   = R*[xf(floor(end/2)) xf(floor(end/2))+2; 0 0];

%% Plot vehicle base
handle(1) = fill(x1,y1,color, 'EdgeColor', color);
hold on
handle(2) = fill(x2,y2,color, 'EdgeColor', color);
handle(3) = fill(x3,y3,color, 'EdgeColor', color);
handle(4) = fill(x4,y4,color, 'EdgeColor', color);
handle(5) = fill([x1(1),x1(end),x2(end),x2(1)],[y1(1),y1(end),y2(end),y2(1)],color, 'EdgeColor', 'none');
% alpha(0.75);


%% Plot Tires
% handle(6) = fill(flWheel(:,1)+offset(1), flWheel(:,2)+offset(2),'k');
% handle(7) = fill(frWheel(:,1)+offset(1), frWheel(:,2)+offset(2),'k');
% handle(8) = fill(rrWheel(:,1)+offset(1), rrWheel(:,2)+offset(2),'k');
% handle(9) = fill(rlWheel(:,1)+offset(1), rlWheel(:,2)+offset(2),'k');

%% Plot Windows
windowColor = [0.2 0.2 0.5];
handle(6) = fill(frontWindow(:,1)+offset(1), frontWindow(:,2)+offset(2),windowColor);
handle(7) = fill(rearWindow(:,1)+offset(1), rearWindow(:,2)+offset(2),windowColor);
handle(8) = fill(sideWindowRight(:,1)+offset(1), sideWindowRight(:,2)+offset(2),windowColor);
handle(9) = fill(sideWindowLeft(:,1)+offset(1), sideWindowLeft(:,2)+offset(2),windowColor);

%% Plot heading
%handle(10) = plot(arrow(1,:)+offset(1),arrow(2,:)+offset(2),'k');

%% Plot outline
outLineColor = color-0.2;
outLineColor(find(outLineColor < 0)) = zeros(size(find(outLineColor < 0)));
handle(10) = plot(outLineLeft(1,:)+offset(1),outLineLeft(2,:)+offset(2),'LineWidth',1,'Color',outLineColor);
handle(11) = plot(outLineRight(1,:)+offset(1),outLineRight(2,:)+offset(2),'LineWidth',1,'Color',outLineColor);

%% Plot Sensors
if plotSensor
    for k = 1:length(sensorParam)
        handle(11+k) = fill(sensor(1,:,k)+offset(1), sensor(2,:,k)+offset(2), sensorColor, 'EdgeColor', sensorColor, 'FaceAlpha', 'flat', 'FaceVertexAlphaData', .2);
    end
end

%set(handle,'LineWidth',2);
if ~axisoff
    %1
    %axis(2*[-vehicleParam.width vehicleParam.width -vehicleParam.length vehicleParam.length]);
    %axis equal;
end