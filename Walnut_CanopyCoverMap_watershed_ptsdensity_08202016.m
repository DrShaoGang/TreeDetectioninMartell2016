x=importdata('walnut_clip_110415.txt' );
data=x';
%figure, plot3(data(1,:),data(2,:),data(3,:),'*')

x_min=min(data(1,:));
y_min=min(data(2,:));
x_max=max(data(1,:));
y_max=max(data(2,:));
z_max=max(data(3,:));
z_min=min(data(3,:));
loc_data=[data(1,:)-x_min;data(2,:)-y_min;data(3,:)];


crown_pts_idx=find(loc_data(3,:)>= 0);
crown_pts=loc_data(:,crown_pts_idx);
%figure, plot3(crown_pts(1,:),crown_pts(2, :), crown_pts(3, :),'*')


%set up resolution
RS=2;

loc_data_ceil=[ceil(crown_pts(1,:)/RS); ceil(crown_pts(2,:)/RS); crown_pts(3,:)];


%loc_data_ceil=[ceil(crown_pts(1,:)); ceil(crown_pts(2,:)); crown_pts(3,:)];
%loc_data_ceil=[ceil(loc_data(1,:)); ceil(loc_data(2,:)); loc_data(3,:)];
%plot3(loc_data_ceil(1,:),loc_data_ceil(2, :), loc_data_ceil(3, :),'*')


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
%figure,imshow(canopy_surface)
%figure,surf(canopy_surface)



canopy_surface_flip=zeros(M,N);
for i=1:M
    for j=1:N
    canopy_surface_flip(M-i+1,j)=canopy_surface(i,j);
    end
end
figure, image(canopy_surface_flip*10)

h=fspecial('disk',3);
canopy_surface_s_flip=imfilter(canopy_surface_flip, h, 'replicate');
figure, image(canopy_surface_s_flip*10)

local_m_candi = imregionalmax(canopy_surface_s_flip);
local_m=(local_m_candi.*canopy_surface_s_flip)>=1.5;

%fir top fix
local_m(89,11)=0;
local_m(105,11)=0;
%local_m(134,152)=0;

tree_candi1=zeros(M,N);
tree_candi1(local_m)=1;
T=2;
for i=T:(M-T)
    for j=T:(N-T)
        sum_candi=sum(sum(tree_candi1((i-T+1):(i+T),(j-T+1):(j+T))));
        if sum_candi>1
            select_candi1=tree_candi1((i-T+1):(i+T),(j-T+1):(j+T)).*canopy_surface_s_flip((i-T+1):(i+T),(j-T+1):(j+T));
            tree_candi1((i-T+1):(i+T),(j-T+1):(j+T))=select_candi1==max(max(select_candi1));
        end
    end
end
tree_candi_selected1=tree_candi1==1;

%canopy_surface_s_m_candi_flip=canopy_surface_s_flip*5;
%canopy_surface_s_m_candi_flip(local_m_candi)=255;
%figure, image(canopy_surface_s_m_candi_flip)


%canopy_surface_s_m_flip=canopy_surface_s_flip*5;
%canopy_surface_s_m_flip(local_m)=255;
%figure, image(canopy_surface_s_m_flip)
%colormap(jet(70))

canopy_surface_s_m_flip=canopy_surface_s_flip*5;
canopy_surface_s_m_flip(tree_candi_selected1)=255;
figure, image(canopy_surface_s_m_flip)
colormap(jet(70))

save('walnut_seed_pdm.txt','canopy_surface_s_m_flip','-ascii')

canopy_gap=zeros(M,N);
for i=1:M
    for j=1:N
        if canopy_surface_s_m_flip(i,j)==0
            canopy_gap(i,j)=1;
        end
    end
end

canopy_gap_L= bwareaopen(canopy_gap, 5);

figure, image(canopy_gap_L*40)
canopy_cover=~canopy_gap_L;
figure, image(canopy_cover*40)

%get inverse image of canopy surface for watershed algorithm
DD=max(max(canopy_surface_s_flip))-canopy_surface_s_flip;
figure,image(DD)


gradmag2 = imimposemin(DD, tree_candi_selected1);

L = watershed(gradmag2);
bgm2 = L == 0;
bgm2(~canopy_cover)=0;
canopy_surface_m1_flip=canopy_surface_s_flip*10;
%canopy_surface_m1_flip=(canopy_surface_s_flip+1)*10;
%canopy_surface_m1(local_m)=255;
canopy_surface_m1_flip(bgm2)=0;
%canopy_surface_m1_flip(~canopy_cover_flip)=0;

backgroundT=canopy_surface_m1_flip<=5;
canopy_surface_m1_flip(backgroundT)=0;

figure, image(canopy_surface_m1_flip)
figure, surf(canopy_surface_m1_flip')

save('walnut_result_pdm.txt','canopy_surface_m1_flip','-ascii')

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
%figure,surf(canopy_height)

canopy_height_flip=zeros(M,N);
for i=1:M
    for j=1:N
    canopy_height_flip(M-i+1,j)=canopy_height(i,j);
    end
end
figure, image(canopy_height_flip)

%finding tree tops
local_max1 = imregionalmax(canopy_height_flip);

canopy_height_m_flip=canopy_height_flip;
canopy_height_m_flip(local_max1)=255;
figure, image(canopy_height_m_flip)
colormap(jet(200))

%smoothing the canopy height model before finding tree tops
h=fspecial('disk',4);
canopy_height_s_flip=imfilter(canopy_height_flip, h, 'replicate');
figure, image(canopy_height_s_flip)
figure,surf(canopy_height_s_flip)


local_max1_s = imregionalmax(canopy_height_s_flip);


tree_candi=zeros(M,N);
tree_candi(local_max1_s)=1;


tree_candi=zeros(M,N);
tree_candi(local_max1_s)=1;
T=2;%4 by 4 boxes
for i=T:(M-T)
    for j=T:(N-T)
        sum_candi=sum(sum(tree_candi((i-T+1):(i+T),(j-T+1):(j+T))));
        if sum_candi>1
            select_candi=tree_candi((i-T+1):(i+T),(j-T+1):(j+T)).*canopy_height_s_flip((i-T+1):(i+T),(j-T+1):(j+T));
            tree_candi((i-T+1):(i+T),(j-T+1):(j+T))=select_candi==max(max(select_candi));
        end
    end
end
tree_candi_selected=tree_candi==1;

%figure, image(tree_candi.*canopy_surface_flip*100)
%canopy_height_m_s_flip=canopy_height_s_flip;
%canopy_height_m_s_flip(local_max1_s)=255;
%figure, image(canopy_height_m_s_flip)
%colormap(jet(200))

canopy_height_m_s_flip=canopy_height_s_flip;
canopy_height_m_s_flip(tree_candi_selected)=255;
figure, image(canopy_height_m_s_flip)
colormap(jet(200))
save('walnut_seed_chm.txt','canopy_height_m_s_flip','-ascii')


canopy_gap1=zeros(M,N);
for i=1:M
    for j=1:N
        if canopy_height_m_s_flip(i,j)==0
            canopy_gap1(i,j)=1;
        end
    end
end


canopy_gap_L1= bwareaopen(canopy_gap1, 5);
figure, image(canopy_gap1*40)
figure, image(canopy_gap_L1*40)
canopy_cover1=~canopy_gap_L1;
figure, image(canopy_cover1*40)


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





%gradmag3 = imimposemin(DD1, local_max1_s);
gradmag3 = imimposemin(DD1, tree_candi_selected);

L1 = watershed(gradmag3);
bgm3 = L1 == 0;
bgm3(~canopy_cover1)=0;
%figure
%imshow(bgm3), title('Watershed ridge lines (bgm)')



%canopy_height_m1_flip=canopy_height_flip;
canopy_height_m1_flip=canopy_height_s_flip*1.5;
%canopy_surface_m1(local_m)=255;
canopy_height_m1_flip(bgm3)=0;
%canopy_height_m1_flip(~canopy_cover1_flip)=0;

backgroundT1=canopy_height_m1_flip<=20;
canopy_height_m1_flip(backgroundT1)=0;


figure, image(canopy_height_m1_flip)

save('walnut_result_chm.txt','canopy_height_m1_flip','-ascii')




figure(11), image(canopy_surface_s_flip*10)
colormap(jet(70))
axis off
figure(11), image(canopy_surface_s_m_flip)
colormap(jet(60))
axis off
figure(11), image(canopy_surface_m1_flip)
colormap(jet(60))
axis off
figure(11), image(canopy_height_s_flip)
colormap(jet(50))
axis off
figure(11), image(canopy_height_m_s_flip)
colormap(jet(70))
axis off
figure(11), image(canopy_height_m1_flip)
colormap(jet(80))
axis off


%watershed use ptsdensity local maxima and canopyheight_cover
gradmag4 = imimposemin(DD1, local_m);
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
height_finder=LL1==103;
figure,image(LL)
figure,image(LL1)
height_finder=LL==59;
height_finder1=LL1==60;

max(max(height_finder.*canopy_height_flip))*0.3084
crown_dia=(sum(sum(height_finder))*0.371612*4/pi)^0.5
max(max(height_finder1.*canopy_height_flip))*0.3084
crown_dia1=(sum(sum(height_finder1))*0.371612*4/pi)^0.5

