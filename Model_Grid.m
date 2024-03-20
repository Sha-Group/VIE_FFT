clc;clear

%  Model data£¨all grids coordinate£©
X=textread('./result/model.txt');

%  draw geometry with a specific material
index=2;  % metal index ***
%rst=find(X(:,4)==index);
X(:,1)=X(:,1).*(X(:,4)==index);
X(:,2)=X(:,2).*(X(:,4)==index);
X(:,3)=X(:,3).*(X(:,4)==index);

% find bounding box
num=length(X);
min_x=min(X(:,1));
max_x=max(X(:,1));
min_y=min(X(:,2));
max_y=max(X(:,2));
min_z=min(X(:,3));
max_z=max(X(:,3));

%  find available coordinates
P=zeros(max_x-min_x+1,max_y-min_y+1,max_z-min_z+1);  %  initial
for m=1:num;
    P(X(m,1)-min_x+1,X(m,2)-min_y+1,X(m,3)-min_z+1)=1;
end

%  find the number of boundary points
index=0;
for x=1:max_x-min_x+1;
   disp('boundary points number£º')
   max_x-min_x+1-x
    for y=1:max_y-min_y+1;
        for z=1:max_z-min_z+1;
            if (x~=1 & y~=1 & z~=1 & x~=max_x-min_x+1 & y~=max_y-min_y+1 & z~=max_z-min_z+1)  %  non-bounding box boundary
                t=abs(P(x-1,y,z)-P(x+1,y,z))+abs(P(x,y+1,z)-P(x,y-1,z))+abs(P(x,y,z+1)-P(x,y,z-1));
                if (t~=0)   %  boundary points of objects
                    index=index+1;
                end
            else   %  bounding box boundary
                if (P(x,y,z)~=0 & ~(x==1 & y==1 & z==1))  %  boundary points of objects
                    index=index+1;
                end
            end
        end
    end
end

%  get the coordinates of boundary points
Q=zeros(index,3);
index=1;
for x=1:max_x-min_x+1;
    disp('boundary points coordinates£º')
    max_x-min_x+1-x
    for y=1:max_y-min_y+1;
        for z=1:max_z-min_z+1;
            if (x~=1 & y~=1 & z~=1 & x~=max_x-min_x+1 & y~=max_y-min_y+1 & z~=max_z-min_z+1)  %  non-bounding box boundary
                t=abs(P(x-1,y,z)-P(x+1,y,z))+abs(P(x,y+1,z)-P(x,y-1,z))+abs(P(x,y,z+1)-P(x,y,z-1));
                if (t~=0)   %  boundary points of objects
                    Q(index,:)=[x,y,z];
                    index=index+1;
                end
            else   %  bounding box boundary
                if (P(x,y,z)~=0 & ~(x==1 & y==1 & z==1))  %  boundary points of objects
                    Q(index,:)=[x,y,z];
                    index=index+1;
                end
            end
        end
    end
end

%  get the cell data
num=length(Q);
index=1;
Quad_x=zeros(4,6*num);
Quad_y=zeros(4,6*num);
Quad_z=zeros(4,6*num);

for m=1:num;
    disp('cell number£º')
    num-m
    %  cell center
    x=Q(m,1);
    y=Q(m,2);
    z=Q(m,3);
    %  6 faces
    Quad_x(:,[index:index+5])=[0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]+x-0.5;
    Quad_y(:,[index:index+5])=[0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]+y-0.5;
    Quad_z(:,[index:index+5])=[0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]+z-0.5;
    index=index+6;
end

%  show staircase model
patch(Quad_x,Quad_y,Quad_z,'r')
axis equal;
view([45,45,45])
title('Grid Generation')