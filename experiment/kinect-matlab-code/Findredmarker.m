function [point, vector] = Findredmarker()

% %getting rgb data from kinect
% vid = videoinput('kinect', 1, 'RGB_1280x960');
% src = getselectedsource(vid);
% vid.FramesPerTrigger = 1;
% start(vid);
% stoppreview(vid);
% I = getdata(vid);

% %depth data
% vid = videoinput('kinect', 2, 'Depth_640x480');
% src = getselectedsource(vid);
% vid.FramesPerTrigger = 1;
% src.DepthMode = 'Near';
% start(vid);
% stoppreview(vid);
% depth = getdata(vid);
%----------------------------------------------------------------/\------------------
addpath('Mex');
clear all
close all

% Create Kinect 2 object and initialize it
% Available sources: 'color', 'depth', 'infrared', 'body_index', 'body',
% 'face' and 'HDface'
k2 = Kin2('color','depth','infrared');

%--------------------------------------------------------------------------------

I = k2.getColor;
depth= k2.getDepth;
%seperate into the hue, saturation, and value
hsv=rgb2hsv(I);
h=hsv(:,:,1);
s=hsv(:,:,2);
v=hsv(:,:,3);

%take the dimensions of the matrix
a=numel(h(:,1));
b=numel(h(1,:));
pixel=zeros(a,b);

%iterate through each pixel and compare to the values we are looking for, trying to find red
for i=1:a
    for j=1:b
        if (abs(h(i,j))<.1 || abs(h(i,j)-1)<.1) && s(i,j)>.6 && v(i,j)>.4
            pixel(i,j)=1;
        else
            pixel(i,j)=0;
        end
    end
end

%BW morph dilate and erode narrow in on the pixels by eliminating small patches
BW=bwmorph(pixel,'erode');
BW=bwmorph(BW,'dilate');
BW=bwmorph(BW,'dilate');
BW=bwmorph(pixel,'erode');

%L locates the clumps and gives them a number. We find how many clumps there are
L=bwlabel(BW,4);
n=max(max(L));
u=[1:n];

%iterate through each clump to find the largest
for i=1:n
    [r,c]=find(L==i);
    u(i)=length(r);
end
[x,y]=max(u);    
[r,c]=find(L==y);

%taking the median point for the largest clump (the center of the marker)
x_mid=median(c);
y_mid=median(r);
%finding the depth at that point (divided by 2 to account for the smaller
%pixel range)
z_mid=depth(round(y_mid/2),round(x_mid/2));

%tracing the edges
dim = size(BW);
[mx i]=max(c);
my=r(i);
boundary = bwtraceboundary(BW,[my, mx],'N');

%find the longest distance from the midpoint to the boundary 
[l num]= distance(boundary,[y_mid,x_mid]);
[mlength mnum]=max(l);

%the point the furthest from the midpoint is the tip of the marker
x_tip=boundary(mnum,2);
y_tip=boundary(mnum,1);
%find the depth at that point
z_tip=depth(round(y_tip/2),round(x_tip/2));

%creating a vector from the midpoint to the edge the greatest distance away
vec1=[double(x_tip)-double(x_mid), double(y_tip)-double(y_mid), double(z_tip)-double(z_mid)];


%correcting for if the vector is pointing down to make sure it points upwards
if vec1(2)>0
    vec1=-1.*vec1;
    x_tip=double(x_mid)+vec1(1);
    y_tip=double(y_mid)+vec1(2);
end

%converting to a unit vector
vec1u=vec1/norm(vec1);

point = [x_tip, y_tip, z_tip];
vector=vec1u;
end