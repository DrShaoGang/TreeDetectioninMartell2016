x=importdata('oak_w_ground.txt ');
data=x';
%figure, plot3(data(1,:),data(2,:),data(3,:),'*')

x_min=min(data(1,:));
y_min=min(data(2,:));
x_max=max(data(1,:));
y_max=max(data(2,:));
z_max=max(data(3,:));
z_min=min(data(3,:));
loc_data=[data(1,:)-x_min;data(2,:)-y_min;data(3,:)];


crown_pts_idx=find(loc_data(3,:)>= 2);
crown_pts=loc_data(:,crown_pts_idx);
%figure, plot3(crown_pts(1,:),crown_pts(2, :), crown_pts(3, :),'*')


%set up resolution
RS=2;

loc_data_ceil=[ceil(crown_pts(1,:)/RS); ceil(crown_pts(2,:)/RS); crown_pts(3,:)];


%loc_data_ceil=[ceil(crown_pts(1,:)); ceil(crown_pts(2,:)); crown_pts(3,:)];
%loc_data_ceil=[ceil(loc_data(1,:)); ceil(loc_data(2,:)); loc_data(3,:)];
plot3(loc_data_ceil(1,:),loc_data_ceil(2, :), loc_data_ceil(3, :),'*')


%create canopy cover model/ pts density map
[m,n]=size(loc_data_ceil);

M=max(loc_data_ceil(2,:));
N=max(loc_data_ceil(1,:));
canopy_surface=zeros(M,N);


for i=1:n
    if loc_data_ceil(2,i)==0
        loc_data_ceil(2,i)=1;
        else if loc_data_ceil(1,i)==0
                loc_data_ceil(1,i)=1;
            end
    end
end


for i=1:n
    if loc_data_ceil(2,i)~=1&&loc_data_ceil(2,i)~=M&&loc_data_ceil(1,i)~=1&&loc_data_ceil(1,i)~=N
        
    canopy_surface(loc_data_ceil(2,i),loc_data_ceil(1,i))=canopy_surface(loc_data_ceil(2,i),loc_data_ceil(1,i))+1;
    canopy_surface(loc_data_ceil(2,i)-1,loc_data_ceil(1,i)-1)=canopy_surface(loc_data_ceil(2,i)-1,loc_data_ceil(1,i)-1)+0.4;
    canopy_surface(loc_data_ceil(2,i)-1,loc_data_ceil(1,i)+1)=canopy_surface(loc_data_ceil(2,i)-1,loc_data_ceil(1,i)+1)+0.4;
    canopy_surface(loc_data_ceil(2,i)+1,loc_data_ceil(1,i)-1)=canopy_surface(loc_data_ceil(2,i)+1,loc_data_ceil(1,i)-1)+0.4;
    canopy_surface(loc_data_ceil(2,i)+1,loc_data_ceil(1,i)+1)=canopy_surface(loc_data_ceil(2,i)+1,loc_data_ceil(1,i)+1)+0.4;
    canopy_surface(loc_data_ceil(2,i)-1,loc_data_ceil(1,i))=canopy_surface(loc_data_ceil(2,i)-1,loc_data_ceil(1,i))+0.6;
    canopy_surface(loc_data_ceil(2,i)+1,loc_data_ceil(1,i))=canopy_surface(loc_data_ceil(2,i)+1,loc_data_ceil(1,i))+0.6;
    canopy_surface(loc_data_ceil(2,i),loc_data_ceil(1,i)-1)=canopy_surface(loc_data_ceil(2,i),loc_data_ceil(1,i)-1)+0.6;
    canopy_surface(loc_data_ceil(2,i),loc_data_ceil(1,i)+1)=canopy_surface(loc_data_ceil(2,i),loc_data_ceil(1,i)+1)+0.6;
    else 
        canopy_surface(loc_data_ceil(2,i),loc_data_ceil(1,i))=canopy_surface(loc_data_ceil(2,i),loc_data_ceil(1,i))+1;
    end
end

figure, image(canopy_surface*10)
figure,imshow(canopy_surface)
figure,surf(canopy_surface)

canopy_surface_flip=zeros(M,N);
for i=1:M
    for j=1:N
    canopy_surface_flip(M-i+1,j)=canopy_surface(i,j);
    end
end
figure, image(canopy_surface_flip*10)



local_m_candi = imregionalmax(canopy_surface_flip);
local_m=(local_m_candi.*canopy_surface_flip)>=4;

tree_candi=zeros(M,N);
tree_candi(local_m)=1;

for i=2:(M-1)
    for j=2:(N-1)
        sum_candi=sum(sum(tree_candi((i-1):(i+1),(j-1):(j+1))));
        if sum_candi>1
            select_candi=tree_candi((i-1):(i+1),(j-1):(j+1)).*canopy_surface_flip((i-1):(i+1),(j-1):(j+1));
            tree_candi((i-1):(i+1),(j-1):(j+1))=select_candi==max(max(select_candi));
        end
    end
end
tree_candi_selected=tree_candi==1;

canopy_surface_m_candi_flip=canopy_surface_flip*7;
%canopy_surface_m_candi_flip(local_m_candi)=255;
canopy_surface_m_candi_flip(tree_candi_selected)=255;
figure, image(canopy_surface_m_candi_flip)
colormap(jet(100))
save('oak_seed.txt','canopy_surface_m_candi_flip','-ascii')

canopy_surface_m_flip=canopy_surface_flip*5;
canopy_surface_m_flip(local_m)=255;
figure, image(canopy_surface_m_flip)



canopy_gap=zeros(M,N);
for i=1:M
    for j=1:N
        if canopy_surface_flip(i,j)<=0
            canopy_gap(i,j)=1;
        end
    end
end

canopy_gap_L= bwareaopen(canopy_gap, 5);
%figure, image(canopy_gap*40)
%figure, image(canopy_gap_L*40)
canopy_cover=~canopy_gap_L;
%figure, image(canopy_cover*40)



%watershed segmentation
%Gradient Magnitude as the Segmentation Function
%hy = fspecial('sobel');
%hx = hy';
%Iy = imfilter(double(canopy_surface*10), hy, 'replicate');
%Ix = imfilter(double(canopy_surface*10), hx, 'replicate');
%gradmag = sqrt(Ix.^2 + Iy.^2);
%figure
%imshow(gradmag,[]), title('Gradient magnitude (gradmag)')

%get inverse image of canopy surface for watershed algorithm
DD=max(max(canopy_surface_flip))-canopy_surface_flip;
figure,image(DD)

%L = watershed(gradmag);
%Lrgb = label2rgb(L);
%figure, imshow(Lrgb), title('Watershed transform of gradient magnitude (Lrgb)')

%separate background
%bw = im2bw(canopy_surface, graythresh(canopy_surface));
%bw = im2bw(canopy_surface, 0.5);
%figure
%imshow(bw), title('Thresholded opening-closing by reconstruction (bw)')

%D = bwdist(~bw);
%D = bwdist(~canopy_cover_flip);

%C=~canopy_cover_flip;
%D(C)=-inf;
%figure, imshow(D)
%DL = watershed(D);
%bgm = DL == 0;
%figure
%imshow(bgm), title('Watershed ridge lines (bgm)')


gradmag2 = imimposemin(DD, local_m_candi);
gradmag2 = imimposemin(DD, tree_candi_selected);
%figure
%imshow(gradmag2,[]), title('Gradient magnitude (gradmag2)')
%figure, image(gradmag2)
L = watershed(gradmag2);
bgm2 = L == 0;
%bgm2(~canopy_cover)=0;
%figure
%imshow(bgm2), title('Watershed ridge lines (bgm)')

L0=L-L;
idx=L(tree_candi_selected);
for i= 1: length(idx)
    tree_idx=L==idx(i);
    L0(tree_idx)=idx(i);
end
bgm2 = L0 == 0;


L0(~canopy_cover_flip)=0;
L0(backgroundT)=0;
Lrgb = label2rgb(L0, 'jet', 'w', 'shuffle');
figure
image(Lrgb)

%figure
%imshow(Lrgb)
%title('Colored watershed label matrix (Lrgb)')

canopy_surface_m1_flip=canopy_surface_flip*10;
%canopy_surface_m1_flip=canopy_height_s_flip;
%canopy_surface_m1_flip=(canopy_surface_flip+5)*7;
%canopy_surface_m1_flip(tree_candi_selected)=255;
canopy_surface_m1_flip(bgm2)=0;
%canopy_surface_m1_flip(~canopy_cover)=0;
backgroundT=canopy_surface_m1_flip<=5;
canopy_surface_m1_flip(backgroundT)=0;

figure, image(canopy_surface_m1_flip)
aix off
colormap(jet(100))

LL=L;
LL(~canopy_cover_flip)=0;
LL(backgroundT)=0;
Lrgb = label2rgb(LL, 'jet', 'w', 'shuffle');
figure
imshow(Lrgb)

figure
image(Lrgb)

figure
image(LL)
title('Colored watershed label matrix (Lrgb)')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%testing canopy height model with similar idea as canopy cover model

[m,n]=size(loc_data_ceil);

M=max(loc_data_ceil(2,:));
N=max(loc_data_ceil(1,:));
canopy_height=zeros(M,N);
canopy_temp=zeros(M,N);

for i=1:n
    if loc_data_ceil(2,i)==0
        loc_data_ceil(2,i)=1;
        else if loc_data_ceil(1,i)==0
                loc_data_ceil(1,i)=1;
            end
    end
end


for i=1:n
    if loc_data_ceil(2,i)~=1&&loc_data_ceil(2,i)~=M&&loc_data_ceil(1,i)~=1&&loc_data_ceil(1,i)~=N
        
    canopy_temp(loc_data_ceil(2,i),loc_data_ceil(1,i))=loc_data_ceil(3,i);
    canopy_temp(loc_data_ceil(2,i)-1,loc_data_ceil(1,i)-1)=loc_data_ceil(3,i)*0.6;
    canopy_temp(loc_data_ceil(2,i)-1,loc_data_ceil(1,i)+1)=loc_data_ceil(3,i)*0.6;
    canopy_temp(loc_data_ceil(2,i)+1,loc_data_ceil(1,i)-1)=loc_data_ceil(3,i)*0.6;
    canopy_temp(loc_data_ceil(2,i)+1,loc_data_ceil(1,i)+1)=loc_data_ceil(3,i)*0.6;
    canopy_temp(loc_data_ceil(2,i)-1,loc_data_ceil(1,i))=loc_data_ceil(3,i)*0.8;
    canopy_temp(loc_data_ceil(2,i)+1,loc_data_ceil(1,i))=loc_data_ceil(3,i)*0.8;
    canopy_temp(loc_data_ceil(2,i),loc_data_ceil(1,i)-1)=loc_data_ceil(3,i)*0.8;
    canopy_temp(loc_data_ceil(2,i),loc_data_ceil(1,i)+1)=loc_data_ceil(3,i)*0.8;
    
    canopy_height(loc_data_ceil(2,i),loc_data_ceil(1,i))=max([canopy_temp(loc_data_ceil(2,i),loc_data_ceil(1,i)) canopy_height(loc_data_ceil(2,i),loc_data_ceil(1,i))]);
    canopy_height(loc_data_ceil(2,i)-1,loc_data_ceil(1,i)-1)=max([canopy_temp(loc_data_ceil(2,i)-1,loc_data_ceil(1,i)-1) canopy_height(loc_data_ceil(2,i)-1,loc_data_ceil(1,i)-1)]);
    canopy_height(loc_data_ceil(2,i)-1,loc_data_ceil(1,i)+1)=max([canopy_temp(loc_data_ceil(2,i)-1,loc_data_ceil(1,i)+1) canopy_height(loc_data_ceil(2,i)-1,loc_data_ceil(1,i)+1)]);
    canopy_height(loc_data_ceil(2,i)+1,loc_data_ceil(1,i)-1)=max([canopy_temp(loc_data_ceil(2,i)+1,loc_data_ceil(1,i)-1) canopy_height(loc_data_ceil(2,i)+1,loc_data_ceil(1,i)-1)]);
    canopy_height(loc_data_ceil(2,i)+1,loc_data_ceil(1,i)+1)=max([canopy_temp(loc_data_ceil(2,i)+1,loc_data_ceil(1,i)+1) canopy_height(loc_data_ceil(2,i)+1,loc_data_ceil(1,i)+1)]);
    canopy_height(loc_data_ceil(2,i)-1,loc_data_ceil(1,i))=max([canopy_temp(loc_data_ceil(2,i)-1,loc_data_ceil(1,i)) canopy_height(loc_data_ceil(2,i)-1,loc_data_ceil(1,i))]);
    canopy_height(loc_data_ceil(2,i)+1,loc_data_ceil(1,i))=max([canopy_temp(loc_data_ceil(2,i)+1,loc_data_ceil(1,i)) canopy_height(loc_data_ceil(2,i)+1,loc_data_ceil(1,i))]);
    canopy_height(loc_data_ceil(2,i),loc_data_ceil(1,i)-1)=max([canopy_temp(loc_data_ceil(2,i),loc_data_ceil(1,i)-1) canopy_height(loc_data_ceil(2,i),loc_data_ceil(1,i)-1)]);
    canopy_height(loc_data_ceil(2,i),loc_data_ceil(1,i)+1)=max([canopy_temp(loc_data_ceil(2,i),loc_data_ceil(1,i)+1) canopy_height(loc_data_ceil(2,i),loc_data_ceil(1,i)+1)]);
    else 
        canopy_height(loc_data_ceil(2,i),loc_data_ceil(1,i))=max([canopy_height(loc_data_ceil(2,i),loc_data_ceil(1,i)) loc_data_ceil(3,i)]);
    end
end

figure, image(canopy_height)
figure,surf(canopy_height)

canopy_height_flip=zeros(M,N);
for i=1:M
    for j=1:N
    canopy_height_flip(M-i+1,j)=canopy_height(i,j);
    end
end
figure, image(canopy_height_flip)
colormap(jet(75))

%finding tree tops
local_max1 = imregionalmax(canopy_height_flip);
figure
imshow(local_max1), title('Regional maxima of opening-closing by reconstruction (fgm)')

canopy_height_m_flip=canopy_height_flip;
canopy_height_m_flip(local_max1)=255;
figure, image(canopy_height_m_flip)
colormap(jet(200))

%smoothing the canopy height model before finding tree tops
h=fspecial('disk',2);
canopy_height_s_flip=imfilter(canopy_height_flip, h, 'replicate');
figure, image(canopy_height_s_flip)
colormap(jet(80))
figure,surf(canopy_height_s_flip)


local_max1_s = imregionalmax(canopy_height_s_flip);

tree_candi1=zeros(M,N);
tree_candi1(local_max1_s)=1;

for i=2:(M-1)
    for j=2:(N-1)
        sum_candi1=sum(sum(tree_candi1((i-1):(i+1),(j-1):(j+1))));
        if sum_candi1>1
            select_candi1=tree_candi1((i-1):(i+1),(j-1):(j+1)).*canopy_height_s_flip((i-1):(i+1),(j-1):(j+1));
            tree_candi1((i-1):(i+1),(j-1):(j+1))=select_candi1==max(max(select_candi1));
        %else
            %tree_candi_selected1((i-1):(i+1),(j-1):(j+1))= tree_candi1((i-1):(i+1),(j-1):(j+1));
        end
    end
end
tree_candi_selected2=tree_candi1==1;

canopy_height_m_s_flip=canopy_height_s_flip;
canopy_height_m_s_flip(tree_candi_selected2)=255;
figure, image(canopy_height_m_s_flip)
colormap(jet(100))

canopy_gap1=zeros(M,N);
for i=1:M
    for j=1:N
        if canopy_height(i,j)==0
            canopy_gap1(i,j)=1;
        end
    end
end


canopy_gap_L1= bwareaopen(canopy_gap1, 5);
figure, image(canopy_gap1*40)
figure, image(canopy_gap_L1*40)
canopy_cover1=~canopy_gap_L1;
figure, image(canopy_cover1*40)

canopy_cover1_flip=zeros(M,N);
for i=1:M
    for j=1:N
    canopy_cover1_flip(M-i+1,j)=canopy_cover1(i,j);
    end
end
figure, image(canopy_cover1_flip*40)
%h=fspecial('disk',1);
%sm_gap=imfilter(gap5, h, 'replicate');


%watershed segmentation
%Gradient Magnitude as the Segmentation Function
%hy = fspecial('sobel');
%hx = hy';
%Iy = imfilter(double(canopy_height_s*10), hy, 'replicate');
%Ix = imfilter(double(canopy_height_s*10), hx, 'replicate');
%gradmag = sqrt(Ix.^2 + Iy.^2);
%figure
%imshow(gradmag,[]), title('Gradient magnitude (gradmag)')

%get inverse image of canopy surface for watershed algorithm
DD1=max(max(canopy_height_s_flip))-canopy_height_s_flip;
figure,image(DD1)


D1 = bwdist(~canopy_cover1_flip);

C1=~canopy_cover1_flip;
D1(C1)=-inf;
figure, imshow(D1)
DL1 = watershed(D1,4);
bgm1 = DL1 == 0;
figure
imshow(bgm1), title('Watershed ridge lines (bgm)')


gradmag3 = imimposemin(DD1, tree_candi_selected2);
figure
imshow(gradmag3,[]), title('Gradient magnitude (gradmag2)')
figure, image(gradmag3)
L1 = watershed(gradmag3);
bgm3 = L1 == 0;
figure
imshow(bgm3), title('Watershed ridge lines (bgm)')

L10=L-L;
idx1=L1(tree_candi_selected2);
for i= 1: length(idx1)
    tree_idx=L1==idx1(i);
    L10(tree_idx)=idx1(i);
end
bgm3 = L10 == 0;


LL1=L10;
LL1(~canopy_cover1_flip)=0;
Lrgb = label2rgb(LL, 'jet', 'w', 'shuffle');
figure
imshow(Lrgb)
title('Colored watershed label matrix (Lrgb)')

canopy_height_m1_flip=canopy_height_s_flip;
%canopy_height_m1_flip(~canopy_cover1_flip)=0
%canopy_surface_m1(local_m)=255;
canopy_height_m1_flip(bgm3)=0;
%canopy_height_m1_flip(tree_candi_selected2)=0;
backgroundT1=canopy_height_m1_flip<=30;
canopy_height_m1_flip(backgroundT1)=0;


figure, image(canopy_height_m1_flip)

colormap(jet(70))


Lrgb1 = label2rgb(L1, 'jet', 'w', 'shuffle');
figure
imshow(Lrgb1)
title('Colored watershed label matrix (Lrgb)')

%watershed use ptsdensity local maxima and canopyheight_cover
gradmag4 = imimposemin(DD1, tree_candi_selected);
figure
imshow(gradmag4,[]), title('Gradient magnitude (gradmag2)')
figure, image(gradmag4)
L2 = watershed(gradmag4);
bgm4 = L2 == 0;
figure
imshow(bgm4), title('Watershed ridge lines (bgm)')

canopy_height_m2_flip=canopy_height_s_flip;
%canopy_surface_m1(local_m)=255;
canopy_height_m2_flip(bgm4)=0;
canopy_height_m2_flip(~canopy_cover1_flip)=0;

figure, image(canopy_height_m2_flip)



%%
%load plot index and tree number
plot_index=[
    108,147,410,460,9; 150,188,425,460,5; 188,220,425,455,3; 218,256,420,460,7;
    256,288,425,460,3; 288,320,425,460,3; 320,360,425,460,8; 356,402,420,464,10;
    105,145,382,414,5; 143,184,388,420,6; 180,216,388,420,4; 216,248,388,420,2;
    248,288,388,425,6; 288,325,388,430,4; 323,360,388,426,8; 358,398,386,424,4;
    109,149,344,388,10; 148,180,348,382,4; 176,218,348,388,8; 218,255,348,388,6;
    252,292,348,388,6; 292,325,348,388,3; 322,364,348,388,11; 108,148,308,348,9;%
    148,180,308,348,3; 176,216,308,348,5; 214,255,310,351,10; 255,292,318,350,4;
    292,325,310,348,3; 322,364,310,348,5; 360,398,315,348,5;  148,180,278,315,5;%
    216,255,278,315,4; 252,292,278,318,9; 292,325,278,315,5;  322,364,278,313,6;
    360,398,278,315,9; 148,180,246,278,4; 176,218,246,278,6;  214,255,246,278,6;
    248,290,246,278,2; 320,360,246,278,8; 360,398,239,278,2;  108,148,210,246,8;%   
    148,180,208,246,6; 176,218,210,246,4; 214,255,210,246,5;  248,293,207,246,6;
    293,325,210,246,1; 320,360,203,246,8; 360,398,200,239,4;  108,150,172,212,9;%
    148,180,172,208,1; 176,218,172,210,8; 214,255,172,210,7;  250,290,170,207,3;
    288,325,172,210,8; 325,363,167,203,3; 360,395,172,205,1;  108,150,140,172,5;%
    146,185,134,174,9; 214,255,138,172,8; 255,290,138,170,2;  288,325,138,172,5;
    325,360,132,167,6; 360,395,128,172,9; 108,150,103,140,11; 148,180,103,138,3;
    176,215,103,138,4; 214,255,103,138,4; 255,290,103,138,4;  288,325,103,135,5;
    325,360,95,130,3;  108,155,70,103,8;  148,180,67,103,7;   180,220,70,103,3;
    255,290,60,100,8;  290,323,55,103,8;  322,360,55,95,6;    360,395,60,95,3;
    108,150,30,70,7;   148,180,30,70,8;   180,220,30,70,10;   255,290,25,60,2;
    322,360,15,58,7;   360,395,20,60,7;
    ];
p_index(:,1:2)=ceil((plot_index(:,1:2)-min(plot_index(:,1))+1)/RS);

p_index(:,3:4)=ceil((max(plot_index(:,4))-plot_index(:,3:4)+1)/RS);

%p_index(:,3:4)=ceil((plot_index(:,3:4)-min(plot_index(:,3))+1)/RS);

p_index(:,5)=plot_index(:,5);

p=62

figure(1), image(canopy_surface_m_flip(p_index(p,4):p_index(p,3),p_index(p,1):p_index(p,2)))
figure(2), image(canopy_height_m_s_flip(p_index(p,4):p_index(p,3),p_index(p,1):p_index(p,2)))
colormap(jet(200))
%figure, image(canopy_surface(p_index(p,3):p_index(p,4),p_index(p,1):p_index(p,2))*10)

%figure, surf(canopy_surface(p_index(p,3):p_index(p,4),p_index(p,1):p_index(p,2))'*10)

figure, image(canopy_surface_m_flip(p_index(p,4):p_index(p,3),p_index(p,1):p_index(p,2)))
figure, image(canopy_surface_m_candi_flip(p_index(p,4):p_index(p,3),p_index(p,1):p_index(p,2)))
figure, image(LL(p_index(p,4):p_index(p,3),p_index(p,1):p_index(p,2)))
figure, image(canopy_cover_flip(p_index(p,4):p_index(p,3),p_index(p,1):p_index(p,2))*40)

height_finder=LL==923;
max(max(height_finder.*canopy_height_flip))*0.3084
crown_dia=(sum(sum(height_finder))*4/pi)^0.5


figure, image(canopy_height_m_s_flip(p_index(p,4):p_index(p,3),p_index(p,1):p_index(p,2)))
colormap(jet(200))
figure, image(LL1(p_index(p,4):p_index(p,3),p_index(p,1):p_index(p,2)))

height_finder1=LL1==674;
max(max(height_finder1.*canopy_height_flip))*0.3084
crown_dia=(sum(sum(height_finder1))*4/pi)^0.5

%%%%%

figure, surf(canopy_surface_m_flip(p_index(p,4):p_index(p,3),p_index(p,1):p_index(p,2)))



figure,surf(canopy_height(p_index(p,4):p_index(p,3),p_index(p,1):p_index(p,2))')

figure, image(canopy_height_s(p_index(p,4):p_index(p,3),p_index(p,1):p_index(p,2)))

figure,surf(canopy_height_s(p_index(p,4):p_index(p,3),p_index(p,1):p_index(p,2))')

figure, image(canopy_height_m(p_index(p,4):p_index(p,3),p_index(p,1):p_index(p,2)))
colormap(jet(200))

figure, image(canopy_height_m_s(p_index(p,4):p_index(p,3),p_index(p,1):p_index(p,2)))
colormap(jet(200))

%get data in sub_plot   %better use code in
      %'Oak_meanshift_auto_moving_plot816_noBoundary' to query points
      [m,n]=size(loc_data_ceil);
      a=1;
      sample=[];
      for j=1:n
          if loc_data_ceil(1,j)<=p_index(p,2)&& loc_data_ceil(1,j)>= p_index(p,1)
              if loc_data_ceil(2,j)<=p_index(p,3) && loc_data_ceil(2,j)>= p_index(p,4)
                  sample(:,a)=loc_data_ceil(:,j);
                  a=a+1;
              end
          end
      end
      
      %
      figure,plot3(sample(1,:),sample(2,:),sample(3,:),'*')
%%      
      
      figure(1), image(canopy_surface_flip*10)
      axis off
      canopy_surface_m_candi_flip(tree_candi_selected)=255;
figure(1), image(canopy_surface_m_candi_flip)
      axis off
      
      figure(1), image(canopy_surface_m1_flip)
      axis off
      
      canopy_surface_m1_flip=canopy_height_s_flip;

canopy_surface_m1_flip(bgm2)=0;
%canopy_surface_m1_flip(~canopy_cover)=0;
backgroundT=canopy_surface_m1_flip<=30;
canopy_surface_m1_flip(backgroundT)=0;

figure(1), image(canopy_surface_m1_flip)
axis off

      figure(1), image(canopy_height_s_flip)
colormap(jet(80))
      axis off
      
      figure(1), image(canopy_height_m_s_flip)
colormap(jet(100))
      axis off
      
      figure(1), image(canopy_height_m1_flip)

      axis off
     
      
      axis off