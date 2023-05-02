%read raw lidar data
x=importdata('walnut_clip_110415.txt');
data=x(:,1:3)';

x_min=min(data(1,:));
y_min=min(data(2,:));
x_max=max(data(1,:));
y_max=max(data(2,:));
z_max=max(data(3,:));
z_min=min(data(3,:));
loc_data=[data(1,:)-x_min;data(2,:)-y_min;data(3,:)];
[m,n]=size(data);
figure, plot3(loc_data(1,:),loc_data(2, :), loc_data(3, :),'*')

%define Dy bandwidth

%average point density of the study site
loc_x_max=max(loc_data(1,:))+1;
loc_y_max=max(loc_data(2,:))+1;

Tot_pts_density=n/(loc_x_max*loc_y_max*0.092903);

%creat chm
 chm=zeros(ceil(loc_x_max),ceil(loc_y_max));
  [Q,W]=size(chm);
  grid=5;
  for i=0:grid:Q
      for h=0:grid:W
          a=1;
          sample=[];
          for j=1:n
              if loc_data(1,j)<=i+grid-1 && loc_data(1,j)>= i
                  if loc_data(2,j)<=h+grid-1 && loc_data(2,j)>= h
                      sample(a)=loc_data(3,j);
                      a=a+1;
                  end
              end
          end
          
          if length(sample)~=0
             chm(i+1:i+grid,h+1:h+grid)=max(sample);
          end
      end
  end
  figure, mesh(chm')
  %calculate area meanheight gapfraction at plot y:169-228 x:338-393
  x1=338;
  x2=393;
  y1=169;
  y2=228;
  sub_area=(228-169)*(393-215)*0.092903;
      
      a=1;
      sample=[];
      for j=1:n
          if loc_data(1,j)<=393&& loc_data(1,j)>= 338
              if loc_data(2,j)<=228 && loc_data(2,j)>= 169
                  sample(:,a)=loc_data(:,j);
                  a=a+1;
              end
          end
      end
  figure,scatter3(sample(1,:),sample(2,:),sample(3,:),20,sample(3,:),'*')
  [K,L]=size(sample);
      sub_pts_density=L/sub_area;
      sub_mean_height=mean(sample(3,:));
[E,R]=find(chm(x1:x2,y1:y2)==0);
      sub_gap_fraction=length(E)*0.092903/sub_area;

  D_bandwidth=1.4*((sub_area*(1-sub_gap_fraction)/(3.14*6)).^0.5);
  

  sample=loc_data;
sample1=sample(1:2,:);
%2d MS
bandwidth = D_bandwidth;
tic
[clustCent,point2cluster,clustMembsCell] = MeanShiftCluster(sample1,bandwidth);
toc

numClust = length(clustMembsCell);

    candi=[];
      uncandi=[];
      pts_nu=[];
      b=1;
      c=1;
      for k=1:numClust
          pts_nu(k)=length(clustMembsCell{k});
          if length(clustMembsCell{k})>5
              
                 candi(b)=k;
                
                 b=b+1;
             
          else uncandi(c)=k;
               c=c+1;
          end
      end
    
        length(candi)
      %vertical gap in each candi clusters
uncandi_logic=pts_nu<=5;


clear ratioH
aa=1;
for k = 1:length(candi)
    myMembers = clustMembsCell{candi(k)};
    myMinH(candi(k)) = min(sample(3,myMembers));
    myMaxH(candi(k)) = max(sample(3,myMembers));
    ratioH(aa)=(myMaxH(candi(k))-myMinH(candi(k)))/myMaxH(candi(k));
   
    aa=aa+1;
    
end

tree_index=ratioH>=0.5;
candi_tree=candi(tree_index);
    
uncandi_logic(candi(~tree_index))=1;
candi_crown=find(uncandi_logic==1)
length(candi_tree)





output_center=[clustCent(1,candi_tree); clustCent(2,candi_tree)]';

figure,plot(output_center(:,1),output_center(:,2),'o','MarkerSize',10)



%3d plot candi
figure(15),clf,hold on
cVec = 'bgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmyk';%, cVec = [cVec cVec];
for k = 1:min(length(candi),length(cVec))
    myMembers = clustMembsCell{candi(k)};
    myClustCen = clustCent(:,candi(k));
    plot3(sample1(1,myMembers),sample1(2,myMembers),sample(3,myMembers),[cVec(k) '.'])
    plot(myClustCen(1),myClustCen(2),'o','MarkerEdgeColor','k','MarkerFaceColor',cVec(k), 'MarkerSize',10)
end
title(['no shifting, numClust:' int2str(length(candi))])
      
%2d plot clusters

figure(9),clf,hold on

cVec = 'bgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmyk';%, cVec = [cVec cVec];
for k = 1:min(numClust,length(cVec))
    myMembers = clustMembsCell{k};
    myClustCen = clustCent(:,k);
    plot(sample1(1,myMembers),sample1(2,myMembers),[cVec(k) '.'])
    plot(myClustCen(1),myClustCen(2),'o','MarkerEdgeColor','k','MarkerFaceColor',cVec(k), 'MarkerSize',5)
end
title(['no shifting, numClust:' int2str(numClust)])

%3d plot
figure(12),clf,hold on
cVec = 'bgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmyk';%, cVec = [cVec cVec];
for k = 1:min(numClust,length(cVec))
    myMembers = clustMembsCell{k};
    myClustCen = clustCent(:,k);
    plot3(sample1(1,myMembers),sample1(2,myMembers),sample(3,myMembers),[cVec(k) '.'])
    plot(myClustCen(1),myClustCen(2),'o','MarkerEdgeColor','k','MarkerFaceColor',cVec(k), 'MarkerSize',10)
end
title(['no shifting, numClust:' int2str(numClust)])


      for k = 1:length(candi)
            myMembers = clustMembsCell{candi(k)};
            tree_pts=sample(:,myMembers);
            height(candi(k))=max(tree_pts(3,:));
            crown(candi(k))=(max(tree_pts(1,:))-min(tree_pts(1,:))+max(tree_pts(2,:))-min(tree_pts(2,:)))/2;
            myClustCen(:,candi(k)) = clustCent(:,candi(k));
            output_loc(4,myMembers)=candi(k);
    
      end
clear h_dbh  
h_dbh(1,:)=height;
h_dbh(2,:)=crown;
h_dbh(3:4,:)=myClustCen;

figure, plot3(h_dbh(3,:),h_dbh(4,:),h_dbh(1,:),'o','MarkerSize',10)
figure, plot3(h_dbh(3,:),h_dbh(4,:),h_dbh(2,:),'*','MarkerSize',10)

%sort x y coordinate
S_h1=sort(h_dbh,'descend' );
S_h=sort(S_h1,2,'descend' );


clear height crown myClustCen h_dbh


 
%check ratio of mean and maximum height
aa=1;
for k = 1:length(candi)
    myMembers = clustMembsCell{candi(k)};
    myMeanH = mean(sample(3,myMembers));
    myMaxH = max(sample(3,myMembers));
    ratioH(aa)=myMeanH/myMaxH;
    aa=aa+1;
    
end
 
 %comaprison
  compare=[plot_index(:,5),round(mean(estimate_treenumber,2)),mode(estimate_treenumber,2)];
  compare(:,3)=compare(:,1)-compare(:,2);
  figure, plot(compare(:,1),compare(:,2),'*')
  
  %accuracy test
  aa=find(compare(:,3)>0);
under=compare(aa,:);
  bb=find(compare(:,3)<0);
over=compare(bb,:);
  cc=find(compare(:,3)==0);
correct=compare(cc,:);
  sum(over(:,3))
  sum(under(:,3))
  
  
    kk=find(compare(:,1)<5);
  kk=find(compare(:,1)>4 & compare(:,1)<9);
  kk=find(compare(:,1)>8);
low=compare(kk,:);
    %accuracy test
  aaa=find(low(:,3)>0);
under1=low(aaa,:);
  bbb=find(low(:,3)<0);
over1=low(bbb,:);
  ccc=find(low(:,3)==0);
correct1=low(ccc,:);
  sum(over1(:,3))
  1+sum(over1(:,3))/sum(low(:,1))
  sum(under1(:,3))
   1-sum(under1(:,3))/sum(low(:,1))+sum(over1(:,3))/sum(low(:,1))
   (sum(low(:,1))-sum(under1(:,3)))/sum(low(:,1))
   
  
  xx=compare(:,2);
yy=compare(:,1);
p=polyfit(xx,yy,1);
yyfit = polyval(p,xx);

yyresid = yy - yyfit;
SSresid = sum(yyresid.^2);
SStotal = (length(yy)-1) * var(yy);
rsq = 1 - SSresid/SStotal
figure, plot(xx,yy,'o',xx,xx)

%3d plot
figure(15),clf,hold on
cVec = 'bgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmyk';%, cVec = [cVec cVec];
for k = 1:min(length(candi),length(cVec))
    myMembers = clustMembsCell{candi(k)};
    myClustCen = clustCent(:,candi(k));
    plot3(sample1(1,myMembers),sample1(2,myMembers),sample(3,myMembers),[cVec(k) '.'])
    plot(myClustCen(1),myClustCen(2),'o','MarkerEdgeColor','k','MarkerFaceColor',cVec(k), 'MarkerSize',10)
end
title(['no shifting, numClust:' int2str(length(candi))])
