function [Dexp1,Arc,Xtest,Output] = loadProbs_ClaraData(D3,Dexp)

addpath(genpath('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey'))

%CS5
fid = fopen('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/SmoothSurfaces/CS5/Proj_shots_BLEND_high.csv', 'rt'); 
D0_ = textscan(fid,'%f %f %f %s %s','headerLines', 1,'Delimiter',',');
fclose(fid);
D1.data = cell2mat(D0_(:,1:3));
D1.textdata =D0_{:,4};
Arc = importdata('EmDiscCS5_Arc.mat');
OBJ1=importdata('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/3D_Modellinng/Blender meshes/CS5.mat'); %read_wobj('Medium_resolution_CS5_version10-fordylan-310.obj')
Xex2 = OBJ1.vertices(OBJ1.objects(8).data.vertices,:);
Xem2 = OBJ1.vertices(OBJ1.objects(4).data.vertices,:);
Xtroph2 = OBJ1.vertices(OBJ1.objects(24).data.vertices,:);
Xam2 = OBJ1.vertices(OBJ1.objects(16).data.vertices,:);
Xsys2 = OBJ1.vertices(OBJ1.objects(20).data.vertices,:);
Xve2 = OBJ1.vertices(OBJ1.objects(12).data.vertices,:);
D2 = [Xam2,zeros(size(Xam2,1),1) ; Xem2,ones(size(Xem2,1),1) ; Xsys2,2*ones(size(Xsys2,1),1) ; Xex2,3*ones(size(Xex2,1),1) ; Xve2,4*ones(size(Xve2,1),1) ; Xtroph2,5*ones(size(Xtroph2,1),1) ];

clear idx
locations = D1.textdata(1:end,1);
for i = 1:length(locations)
    try
    idx(i,1) = find(strcmp(D3.textdata(2:end,5) , ['P1_E15B_',strrep(locations{i},'-','_')]));   
    catch
    idx(i,1) = find(strcmp(D3.textdata(2:end,5) , ['P2_E15B_',strrep(locations{i},'-','_')]));    
    end    
end
CellIDs = D3.textdata(idx+1,1);
CellIDsu = D3.textdata(idx+1,6);
D3subs = D3.data(idx,:);
CellIDs2 = D3.textdata(idx+1,11);
CellIDs4 = D3.textdata(idx+1,5);
genes = Dexp.textdata(2:end,1);
headers = Dexp.textdata(1,2:end);
headers = strrep(headers,'"','');
clear indx
for i = 1:length(CellIDs)
    try
    indx(i,1) = find(strcmp(headers , CellIDs{i} ));
    catch
    end
end

CellIDsCS5 = CellIDs(find(indx~=0));
CellIDs2 = CellIDs2(find(indx~=0));
CellIDs4 = CellIDs4(find(indx~=0));
CellIDsu = CellIDsu(find(indx~=0));

XYZ = D1.data(find(indx~=0),:);

Dexp1 = Dexp.data(:,indx(find(indx~=0)));
genes = Dexp.textdata(2:end,1);
XtrainCS5 = XYZ;
XtrainCS5(:,1:3) = XtrainCS5(:,1:3) / 400;
XtestCS5 = D2(:,[1:3]);
XtestCS5(:,1:3) =XtestCS5(:,1:3)/400;

LabBinCS5 = NaN*zeros(length(CellIDs2),1);
LabBinCS5(find(strcmp(CellIDs2,'Am_CS5')==1)) = 0;
LabBinCS5(find(strcmp(CellIDs2,'EmDisc_CS5')==1)) = 1;
LabBinCS5(find(strcmp(CellIDs2,'SYS_CS5')==1)) = 2;
LabBinCS5(find(strcmp(CellIDs2,'ExMes_CS5')==1)) = 3;
LabBinCS5(find(strcmp(CellIDs2,'VE_CS5')==1)) = 4;
LabBinCS5(find(strcmp(CellIDs2,'Tb_CS5')==1)) = 5;

indCS50 = find(D2(:,4)==0);
indCS51 = find(D2(:,4)==1);
indCS52 = find(D2(:,4)==2);
indCS53 = find(D2(:,4)==3);
indCS54 = find(D2(:,4)==4);
indCS55 = find(D2(:,4)==5);


%Load the CS6 datasets
fid = fopen('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/SmoothSurfaces/CS6/Proj_shots_BLEND_high_update.csv', 'rt');
D0_ = textscan(fid,'%f %f %f %s %s','headerLines', 1,'Delimiter',',');
fclose(fid);
D1.data = cell2mat(D0_(:,1:3));
D1.textdata = D0_{:,4};
Arc = importdata('EmDiscCS6_Arc.mat');
OBJ2 = importdata('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/3D_Modellinng/Blender meshes/CS6.mat');
Xex2    = [OBJ2.vertices(OBJ2.objects(28).data.vertices,:)]; 
Xem2    = [OBJ2.vertices(OBJ2.objects(8).data.vertices,:)];
Xtroph2 = [OBJ2.vertices(OBJ2.objects(20).data.vertices,:)];
Xam2    = [OBJ2.vertices(OBJ2.objects(4).data.vertices,:)];
Xsys2   = [OBJ2.vertices(OBJ2.objects(16).data.vertices,:)];
Xve2    = [OBJ2.vertices(OBJ2.objects(24).data.vertices,:)];
Xpgc2   = [OBJ2.vertices(OBJ2.objects(12).data.vertices,:)];
D2 = [Xam2,zeros(size(Xam2,1),1) ; Xem2,ones(size(Xem2,1),1) ; Xsys2,2*ones(size(Xsys2,1),1) ; Xex2,3*ones(size(Xex2,1),1) ; Xve2,4*ones(size(Xve2,1),1) ; Xtroph2,5*ones(size(Xtroph2,1),1); Xpgc2,6*ones(size(Xpgc2,1),1)  ];

locations = D1.textdata(1:end,1);
%idx = NaN*zeros(length(locations),1);
for i = 1:length(locations)
    try
        idx(i,1) = find(strcmp(D3.textdata(2:end,5) , ['E15C_',strrep(locations{i},'-','_')]));   
    end    
end
CellIDs  = D3.textdata(idx+1,1);
D3subs   = D3.data(idx,:);
CellIDsu = D3.textdata(idx+1,6);
CellIDs2 = D3.textdata(idx+1,11);
CellIDs4 = D3.textdata(idx+1,5);

genes   = Dexp.textdata(2:end,1);
headers = Dexp.textdata(1,2:end);
headers = strrep(headers,'"','');
headers = strrep(headers,'X3536STDY','3536STDY');
clear indx
for i = 1:length(CellIDs)
    try
    indx(i,1) = find(strcmp(headers , CellIDs{i} ));
    catch
    end
end

%Are there any that don't exist?
CellIDsCS6 = CellIDs(find(indx~=0));
CellIDs2 = CellIDs2(find(indx~=0));
CellIDs4 = CellIDs4(find(indx~=0));
CellIDsu = CellIDsu(find(indx~=0));
XYZ = D1.data(find(indx~=0),:);
Dexp1 = Dexp.data(:,indx(find(indx~=0)));
genes = Dexp.textdata(2:end,1);
XtrainCS6 = XYZ;
XtrainCS6(:,1:3) = XtrainCS6(:,1:3) / 400;
XtestCS6 = D2(:,[1:3]);
XtestCS6(:,1:3) = XtestCS6(:,1:3)/400;

LabBinCS6 = NaN*zeros(length(CellIDs2),1);
LabBinCS6(find(strcmp(CellIDs2,'Am_CS6')==1)) = 0;
LabBinCS6(find(strcmp(CellIDs2,'EmDisc_CS6')==1)) = 1;
LabBinCS6(find(strcmp(CellIDs2,'VE_CS6')==1)) = 4;
LabBinCS6(find(strcmp(CellIDs2,'SYS_CS6')==1)) = 2;
LabBinCS6(find(strcmp(CellIDs2,'Tb_CS6')==1)) = 5;
LabBinCS6(find(strcmp(CellIDs2,'ExMes_CS6')==1)) = 3;
LabBinCS6(find(strcmp(CellIDs2,'PGC_CS6')==1)) = 6;

indCS60 = find(D2(:,4)==0);
indCS61 = find(D2(:,4)==1);
indCS62 = find(D2(:,4)==2);
indCS63 = find(D2(:,4)==3);
indCS64 = find(D2(:,4)==4);
indCS65 = find(D2(:,4)==5);
indCS66 = find(D2(:,4)==6);

%Carnnegie Stage 7
fid = fopen('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/SmoothSurfaces/CS7/Proj_shots_BLEND_high.csv', 'rt'); 
D0_ = textscan(fid,'%f %f %f %s %s','headerLines', 1,'Delimiter',',');
fclose(fid);
D1.data = cell2mat(D0_(:,1:3));
D1.textdata =D0_{:,4};

OBJ3=importdata('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/3D_Modellinng/Blender meshes/CS7.mat'); %read_wobj('Medium_resolution_CS5_version10-fordylan-310.obj')
%OBJ4=importdata('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/3D_Modellinng/Blender meshes/CS7_section_cut.mat'); %read_wobj('Medium_resolution_CS5_version10-fordylan-310.obj')

Xex2 = [OBJ3.vertices(OBJ3.objects(24).data.vertices,:);
        OBJ3.vertices(OBJ3.objects(26).data.vertices,:);
        OBJ3.vertices(OBJ3.objects(34).data.vertices,:)];
Xem2    = OBJ3.vertices(OBJ3.objects(20).data.vertices,:);
Xtroph2 = OBJ3.vertices(OBJ3.objects(4).data.vertices,:);
Xam2    = OBJ3.vertices(OBJ3.objects(16).data.vertices,:);
Xsys2   = [OBJ3.vertices(OBJ3.objects(8).data.vertices,:);
           OBJ3.vertices(OBJ3.objects(12).data.vertices,:)];
Xexst2  = OBJ3.vertices(OBJ3.objects(30).data.vertices,:);
D2 = [Xam2,zeros(size(Xam2,1),1) ; Xem2,ones(size(Xem2,1),1) ; Xexst2, 1*ones(size(Xexst2,1),1); Xsys2,2*ones(size(Xsys2,1),1) ; Xex2,3*ones(size(Xex2,1),1) ;  Xtroph2,5*ones(size(Xtroph2,1),1) ];

%idx = NaN*zeros(length(locations),1);
locations = D1.textdata(1:end,1);
for i = 1:length(locations)
    try
    idx(i,1) = find(strcmp(D3.textdata(2:end,5) , locations{i}));   
    end    
end
CellIDs = D3.textdata(idx+1,1);
D3subs = D3.data(idx,:);
CellIDsu = D3.textdata(idx+1,6);
CellIDs3 = D3.textdata(idx+1,9);
CellIDs2 = D3.textdata(idx+1,11);
CellIDs4 = D3.textdata(idx+1,5);
genes = Dexp.textdata(2:end,1);
headers = Dexp.textdata(1,2:end);
headers = strrep(headers,'"','');
clear indx
for i = 1:length(CellIDs)
    try
        indx(i,1) = find(strcmp(headers , CellIDs{i} ));
    catch
        try
            indx(i,1) = find(strcmp(headers , strrep(CellIDs{i}, '-','.' )));   
        catch
        end
    end
end
CellIDsCS7 = CellIDs(find(indx~=0));
CellIDs2 = CellIDs2(find(indx~=0));
CellIDs3 = CellIDs3(find(indx~=0));
CellIDsu = CellIDsu(find(indx~=0));
CellIDs4 = CellIDs4(find(indx~=0));
XYZ = D1.data(find(indx~=0),:);
Dexp1 = Dexp.data(:,indx(find(indx~=0)));
genes = Dexp.textdata(2:end,1);
XtrainCS7 = XYZ;
XtrainCS7(:,1:3) = XtrainCS7(:,1:3) / 800;
XtestCS7 = D2(:,[1:3]);
XtestCS7(:,1:3) =XtestCS7(:,1:3) / 800;
LabBinCS7 = NaN*zeros(length(CellIDs2),1);
LabBinCS7(find(strcmp(CellIDs2,'Am_CS7')==1)) = 0;
LabBinCS7(find(strcmp(CellIDs2,'EmDisc_CS7')==1)) = 1;
LabBinCS7(find(strcmp(CellIDs2,'SYS_CS7')==1)) = 2;
LabBinCS7(find(strcmp(CellIDs2,'ExMes_CS7')==1)) = 3;
LabBinCS7(find(strcmp(CellIDs2,'Tb_CS7')==1)) = 5;
indCS70 = find(D2(:,4)==0);
indCS71 = find(D2(:,4)==1);
indCS72 = find(D2(:,4)==2);
indCS73 = find(D2(:,4)==3);
indCS75 = find(D2(:,4)==5);


Da1 = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/Correlation_alltissues_pearson_Int.csv')
%Da1 = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/Correlation_pearson_Int.csv')
ID1 = Da1.textdata(2:end,1)
ID2 = Da1.textdata(1,2:end)

ID1 = strrep(ID1,'"','');
ID2 = strrep(ID2,'"','');

Da1 = Da1.data;

ind0_CS5 = find(LabBinCS5==0 | LabBinCS5==1 | LabBinCS5==4 | LabBinCS5==6);
ind0_CS6 = find(LabBinCS6==0 | LabBinCS6==1 | LabBinCS6==4 | LabBinCS6==6);
ind0_CS7 = find(LabBinCS7==0 | LabBinCS7==1 | LabBinCS7==4 | LabBinCS7==6);

XB0 = zeros(length(CellIDsCS5),size(Da1,2));
 for i = 1:length(CellIDsCS5)
     try
     IndsMatch(i,1) = find(strcmp([CellIDsCS5{i}],ID1)==1);
     XB0(i,:) = Da1(IndsMatch(i,1),:);
     catch
     end
 end

XB1 = zeros(length(CellIDsCS6),size(Da1,2));
for i = 1:length(CellIDsCS6)
     try
     IndsMatch(i,1) = find(strcmp(['X' CellIDsCS6{i}],ID1)==1);
     XB1(i,:) = Da1(IndsMatch(i,1),:);
     catch
         XB1(i,:) = 0;
     end
end

XB2 = zeros(length(CellIDsCS7),size(Da1,2));
 for i = 1:length(CellIDsCS7)
     try
     IndsMatch(i,1) = find(strcmp(['X' CellIDsCS7{i}],ID1)==1);
     XB2(i,:) = Da1(IndsMatch(i,1),:);
     catch         
         try
     IndsMatch(i,1) = find(strcmp([CellIDsCS7{i}],ID1)==1);
     XB2(i,:) = Da1(IndsMatch(i,1),:);
         catch        
             try
     IndsMatch(i,1) = find(strcmp(strrep(CellIDsCS7{i},'SLX-','SLX.'),ID1)==1);
     XB2(i,:) = Da1(IndsMatch(i,1),:);
             catch
             end
             end
     end
 end

set(0, 'DefaultTextInterpreter', 'none') 
MD{1} = XB0;
MD{2} = XB1;
MD{3} = XB2;

OBJ2=importdata('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/3D_Modellinng/Blender meshes/CS6.mat'); %read_wobj('Medium_resolution_CS5_version10-fordylan-310.obj')
OBJ2.vertices = Quaternion3(19.9,[0,0,1],OBJ2.vertices); 
OBJ2B=OBJ2;
OBJ2B.vertices = OBJ2B.vertices - repmat([100,0,0],size(OBJ2B.vertices,1),1); %Quaternion3(10,[0,1,0],OBJ2.vertices); 
OBJ2C=OBJ2;
OBJ2C.vertices = OBJ2C.vertices - repmat([200,0,0],size(OBJ2C.vertices,1),1); %Quaternion3(10,[0,1,0],OBJ2.vertices); 
OBJ2D=OBJ2;
OBJ2D.vertices = OBJ2D.vertices - repmat([250,0,0],size(OBJ2D.vertices,1),1); %Quaternion3(10,[0,1,0],OBJ2.vertices); 
OBJ2E=OBJ2;
OBJ2E.vertices = OBJ2E.vertices - repmat([80,0,0],size(OBJ2E.vertices,1),1); %Quaternion3(10,[0,1,0],OBJ2.vertices); 
Xtrain2 = Quaternion3(19.9,[0,0,1],XtrainCS6); 
set(0, 'DefaultTextInterpreter', 'none')
Arc22 = Quaternion3(19.9,[0,0,1],Arc); 


%GP priors
pg1 = {@priorGauss,0,1};  pg2 = {@priorGauss,0,1};  pg3 = {@priorGauss,log(0.5),1};  pc = {@priorClamped};
prior1.mean = {[]}; prior1.cov  = {pg1;pg2}; prior1.lik  = {pg3};
im1 = {@infPrior,@infExact,prior1};                % inference method


Idents1 = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/Idents_alltissues1_Int.csv')
Idents2 = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/Idents_alltissues2_Int.csv')

keyboard

uLab = [0,1,6];
%k = 4;
k = 4;
for i = 1:315
    
    Xtrainin1 = XtrainCS5;
    Xtrainin2 = Xtrain2;
    Xtrainin3 = XtrainCS7;
       
    Ytrainin1 = MD{1}(:,i);
    Ytrainin2 = MD{2}(:,i);
    Ytrainin3 = MD{3}(:,i);
    
    Ymean = ([Ytrainin1;Ytrainin2;Ytrainin3]);
    Ymean = mean(Ymean(find(Ymean~=0)));
    Ytrainin1 = 1 ./ ( 1 + exp(- k*(Ytrainin1-Ymean) ));
    Ytrainin2 = 1 ./ ( 1 + exp(- k*(Ytrainin2-Ymean) ));
    Ytrainin3 = 1 ./ ( 1 + exp(- k*(Ytrainin3-Ymean) ));

    %Ytrainin = Ytrainin/sum(Ytrainin);
    for j = 1:length(uLab)
        substrain = find(LabBinCS5==uLab(j));      
        hyp1.cov  = [log(2); log(0.1) ]; 
        hyp1.lik  = [log(0.1/2)];  
        hyp1.mean = mean(Ytrainin1(substrain,1));       
        par  = {@meanConst,@covSEiso,'likGauss',Xtrainin1(substrain,:), Ytrainin1(substrain,1)};
        hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});
        [m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrainin1(substrain,:), Ytrainin1(substrain,1), [OBJ1.vertices]/400);        
        ConvCS5{i,j} = m_1';
        ConvfCS5{i,j} = s_1';  
    end
    

    for j = 1:length(uLab)
        substrain = find(LabBinCS6==uLab(j));      
        hyp1.cov  = [log(2); log(0.1) ]; 
        hyp1.lik  = [log(0.1/2)];  
        hyp1.mean = mean(Ytrainin2(substrain,1));       
        par  = {@meanConst,@covSEiso,'likGauss',Xtrainin2(substrain,:), Ytrainin2(substrain,1)};
        hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});
        [m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrainin2(substrain,:), Ytrainin2(substrain,1), [OBJ2.vertices]/400);        
        ConvCS6{i,j} = m_1';
        ConvfCS6{i,j} = s_1';  
    end    

    for j = 1:length(uLab)
        substrain = find(LabBinCS7==uLab(j));      
        hyp1.cov  = [log(2); log(0.1) ]; 
        hyp1.lik  = [log(0.1/2)];  
        hyp1.mean = mean(Ytrainin3(substrain,1));       
        par  = {@meanConst,@covSEiso,'likGauss',Xtrainin3(substrain,:), Ytrainin3(substrain,1)};
        hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});
        [m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrainin3(substrain,:), Ytrainin3(substrain,1), [OBJ3.vertices]/800);        
        ConvCS7{i,j} = m_1';
        ConvfCS7{i,j} = s_1';  
    end     
    
end

%First get the most varied cells and 
for i = 1:315
    m1(i,1) = max([ConvCS6{i,1},ConvCS6{i,2}]) - min(max([ConvCS6{i,1};ConvCS6{i,2}]));
end


listType = find(strcmp(Idents2.textdata(2:end,2),"Amnioid_#2")==1);
m2 = m1(listType);
[aa bb] =  sort(m2,'descend')
h = figure(1)
for i = 1:min(length(listType),25)
subplot(5,5,i);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',ConvCS6{listType(bb(i)),1}','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',ConvCS6{listType(bb(i)),2}','FaceColor','interp','LineStyle','none');
%f5 = patch('Faces',OBJ2C.objects(24).data.vertices,'Vertices',OBJ2C.vertices,'FaceVertexCData',Conv{i,3}','FaceColor','interp','LineStyle','none');
%f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',ConvCS6{listType(i),3}','FaceColor','interp','LineStyle','none');
view([-239.7545,27.0685])
%view([-92.8171,-38.7904])
%colorbar
set(gca,'clim',[.2 .7] )
axis equal
axis off
camlight('headlight')
end
%savefig('/Volumes/GoogleDrive/Shared drives/Munger et al./Chris'' Sequencing Dungeon/Updated Figures/3DFigures/Am_reps.fig')
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 10])
%print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/EmD_All.png'])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 30])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/AmB_reps_k=4.png'])


clf
listType = find(strcmp(Idents2.textdata(2:end,2),"Amnioid_#1")==1);
m2 = m1(listType);
[aa bb] =  sort(m2,'descend')
h = figure(1)
for i = 1:min(length(listType),25)
subplot(5,5,i);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',ConvCS6{listType(bb(i)),1}','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',ConvCS6{listType(bb(i)),2}','FaceColor','interp','LineStyle','none');
%f5 = patch('Faces',OBJ2C.objects(24).data.vertices,'Vertices',OBJ2C.vertices,'FaceVertexCData',Conv{i,3}','FaceColor','interp','LineStyle','none');
%f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',ConvCS6{listType(i),3}','FaceColor','interp','LineStyle','none');
view([-239.7545,27.0685])
%view([-92.8171,-38.7904])
%colorbar
set(gca,'clim',[.2 .7] )
axis equal
axis off
camlight('headlight')
end
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 30])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/Am_reps_k=4.png'])


clf
listType = find(strcmp(Idents2.textdata(2:end,2),"EmDisc_#1")==1);
m2 = m1(listType);
[aa bb] =  sort(m2,'descend')
h = figure(1)
for i = 1:min(length(listType),25)
subplot(5,5,i);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',ConvCS6{listType(bb(i)),1}','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',ConvCS6{listType(bb(i)),2}','FaceColor','interp','LineStyle','none');
%f5 = patch('Faces',OBJ2C.objects(24).data.vertices,'Vertices',OBJ2C.vertices,'FaceVertexCData',Conv{i,3}','FaceColor','interp','LineStyle','none');
%f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',ConvCS6{listType(i),3}','FaceColor','interp','LineStyle','none');
%view([-239.7545,27.0685])
view([-92.8171,-38.7904])
%colorbar
set(gca,'clim',[.2 .7] )
axis equal
axis off
camlight('right')
material shiny 
material dull 

end
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 30])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/EmD_reps_k=4.png'])




material([1 1 1])


clf
listType = find(strcmp(Idents2.textdata(2:end,2),"EmDisc_#2")==1);
m2 = m1(listType);
[aa bb] =  sort(m2,'descend')
h = figure(1)
for i = 1:min(length(listType),25)
subplot(5,5,i);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',ConvCS6{listType(bb(i)),1}','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',ConvCS6{listType(bb(i)),2}','FaceColor','interp','LineStyle','none');
%f5 = patch('Faces',OBJ2C.objects(24).data.vertices,'Vertices',OBJ2C.vertices,'FaceVertexCData',Conv{i,3}','FaceColor','interp','LineStyle','none');
%f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',ConvCS6{listType(i),3}','FaceColor','interp','LineStyle','none');
%view([-239.7545,27.0685])
view([-92.8171,-38.7904])
%colorbar
set(gca,'clim',[.2 .7] )
axis equal
axis off
camlight('right')
end
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 30])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/EmDisc_reps_k=4.png'])



%for i = 1:315
%    Xtrainin = Xtrain2;
%    Ytrainin = MD{1}(:,i);%sum(MD{Dataset}(:,CL{i}),2) / sum(DL.data(CL{i}));    
%    for j = 1:length(uLab)
%        substrain = find(LabBin==uLab(j));      
%        hyp1.cov  = [log(2); log(0.1) ]; 
%        hyp1.lik  = [log(0.1/2)];  
%        hyp1.mean = mean(Ytrainin(substrain,1));       
%        par  = {@meanConst,@covSEiso,'likGauss',Xtrainin(substrain,:), Ytrainin(substrain,1)};
%        hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});
%        [m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrainin(substrain,:), Ytrainin(substrain,1), [OBJ2.vertices]/400);        
%        Conv{i,j} = m_1';
%        Convf{i,j} = s_1';  
%end
%end



%Averages.......
listType = find(strcmp(Idents2.textdata(2:end,2),"EmDisc_#1")==1);
X_EmD_CS5_1 = ConvCS5{listType(1),1};X_EmD_CS5_2 = ConvCS5{listType(1),2};
for i = 2:length(listType)
X_EmD_CS5_1 = X_EmD_CS5_1 + ConvCS5{listType(i),1};
X_EmD_CS5_2 = X_EmD_CS5_2 + ConvCS5{listType(i),2};
end
X_EmD_CS5_1 = X_EmD_CS5_1/i;
X_EmD_CS5_2 = X_EmD_CS5_2/i;
X_EmD_CS6_1 = ConvCS6{listType(1),1};X_EmD_CS6_2 = ConvCS6{listType(1),2};X_EmD_CS6_3 = ConvCS6{listType(1),3};
for i = 2:length(listType)
X_EmD_CS6_1 = X_EmD_CS6_1 + ConvCS6{listType(i),1};
X_EmD_CS6_2 = X_EmD_CS6_2 + ConvCS6{listType(i),2};
X_EmD_CS6_3 = X_EmD_CS6_3 + ConvCS6{listType(i),3};
end
X_EmD_CS6_1 = X_EmD_CS6_1/i;
X_EmD_CS6_2 = X_EmD_CS6_2/i;
X_EmD_CS6_3 = X_EmD_CS6_3/i;
X_EmD_CS7_1 = ConvCS7{listType(1),1};X_EmD_CS7_2 = ConvCS7{listType(1),2};
for i = 2:length(listType)
X_EmD_CS7_1 = X_EmD_CS7_1 + ConvCS7{listType(i),1};
X_EmD_CS7_2 = X_EmD_CS7_2 + ConvCS7{listType(i),2};
end
X_EmD_CS7_1 = X_EmD_CS7_1/i;
X_EmD_CS7_2 = X_EmD_CS7_2/i;

save('X_EmD_CS6_1_k=4.mat','X_EmD_CS6_1')
save('X_EmD_CS6_2_k=4.mat','X_EmD_CS6_2')
subplot(1,3,1);
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',X_EmD_CS5_2','FaceColor','interp','LineStyle','none');
view([-59.3039,41.1722])
%colorbar
set(gca,'clim',[.2 .7] )
axis equal
axis off
camlight('headlight')
subplot(1,3,2);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',X_EmD_CS6_1','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',X_EmD_CS6_2','FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',X_EmD_CS6_3','FaceColor','interp','LineStyle','none');
view([ 165.7876,21.7629])
%colorbar
set(gca,'clim',[.2 .7] )
axis equal
axis off
camlight('headlight')
subplot(1,3,3);
f1 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_EmD_CS7_2','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ3.objects(16).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_EmD_CS7_1','FaceColor','interp','LineStyle','none');
%colorbar
set(gca,'clim',[.2 .7] )
view([24.5784,1.2879])
axis equal
axis off
camlight('headlight')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 10])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/EmD_All.png'])

listType = find(strcmp(Idents2.textdata(2:end,2),"EmDisc_#2")==1);
X_EmDisc_CS5_1 = ConvCS5{listType(1),1};X_EmDisc_CS5_2 = ConvCS5{listType(1),2};
for i = 2:length(listType)
X_EmDisc_CS5_1 = X_EmDisc_CS5_1 + ConvCS5{listType(i),1};
X_EmDisc_CS5_2 = X_EmDisc_CS5_2 + ConvCS5{listType(i),2};
end
X_EmDisc_CS5_1 = X_EmDisc_CS5_1/i;
X_EmDisc_CS5_2 = X_EmDisc_CS5_2/i;
X_EmDisc_CS6_1 = ConvCS6{listType(1),1};X_EmDisc_CS6_2 = ConvCS6{listType(1),2};X_EmDisc_CS6_3 = ConvCS6{listType(1),3};
for i = 2:length(listType)
X_EmDisc_CS6_1 = X_EmDisc_CS6_1 + ConvCS6{listType(i),1};
X_EmDisc_CS6_2 = X_EmDisc_CS6_2 + ConvCS6{listType(i),2};
X_EmDisc_CS6_3 = X_EmDisc_CS6_3 + ConvCS6{listType(i),3};
end
X_EmDisc_CS6_1 = X_EmDisc_CS6_1/i;
X_EmDisc_CS6_2 = X_EmDisc_CS6_2/i;
X_EmDisc_CS6_3 = X_EmDisc_CS6_3/i;
X_EmDisc_CS7_1 = ConvCS7{listType(1),1};X_EmDisc_CS7_2 = ConvCS7{listType(1),2};
for i = 2:length(listType)
X_EmDisc_CS7_1 = X_EmDisc_CS7_1 + ConvCS7{listType(i),1};
X_EmDisc_CS7_2 = X_EmDisc_CS7_2 + ConvCS7{listType(i),2};
end
X_EmDisc_CS7_1 = X_EmDisc_CS7_1/i;
X_EmDisc_CS7_2 = X_EmDisc_CS7_2/i;

save('X_EmDisc_CS6_1_k=4.mat','X_EmDisc_CS6_1')
save('X_EmDisc_CS6_2_k=4.mat','X_EmDisc_CS6_2')

clf
subplot(1,3,1);
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',X_EmDisc_CS5_2','FaceColor','interp','LineStyle','none');
view([-59.3039,41.1722])
%colorbar
set(gca,'clim',[.2 .7] )
axis equal
axis off
camlight('headlight')
subplot(1,3,2);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',X_EmDisc_CS6_1','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',X_EmDisc_CS6_2','FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',X_EmDisc_CS6_3','FaceColor','interp','LineStyle','none');
view([ 165.7876,21.7629])
set(gca,'clim',[.2 .7] )
axis equal
axis off
camlight('headlight')
subplot(1,3,3);
f1 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_EmDisc_CS7_2','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ3.objects(16).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_EmDisc_CS7_1','FaceColor','interp','LineStyle','none');
%colorbar
set(gca,'clim',[.2 .7] )
view([24.5784,1.2879])
axis equal
axis off
camlight('headlight')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 10])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/EmDisc_All.png'])

listType = find(strcmp(Idents2.textdata(2:end,2),"Amnioid_#1")==1);
X_Am_CS5_1 = ConvCS5{listType(1),1};X_Am_CS5_2 = ConvCS5{listType(1),2};
for i = 2:length(listType)
X_Am_CS5_1 = X_Am_CS5_1 + ConvCS5{listType(i),1};
X_Am_CS5_2 = X_Am_CS5_2 + ConvCS5{listType(i),2};
end
X_Am_CS5_1 = X_Am_CS5_1/i;
X_Am_CS5_2 = X_Am_CS5_2/i;
X_Am_CS6_1 = ConvCS6{listType(1),1};X_Am_CS6_2 = ConvCS6{listType(1),2};X_Am_CS6_3 = ConvCS6{listType(1),3};
for i = 2:length(listType)
X_Am_CS6_1 = X_Am_CS6_1 + ConvCS6{listType(i),1};
X_Am_CS6_2 = X_Am_CS6_2 + ConvCS6{listType(i),2};
X_Am_CS6_3 = X_Am_CS6_3 + ConvCS6{listType(i),3};
end
X_Am_CS6_1 = X_Am_CS6_1/i;
X_Am_CS6_2 = X_Am_CS6_2/i;
X_Am_CS6_3 = X_Am_CS6_3/i;
X_Am_CS7_1 = ConvCS7{listType(1),1};X_Am_CS7_2 = ConvCS7{listType(1),2};
for i = 2:length(listType)
X_Am_CS7_1 = X_Am_CS7_1 + ConvCS7{listType(i),1};
X_Am_CS7_2 = X_Am_CS7_2 + ConvCS7{listType(i),2};
end
X_Am_CS7_1 = X_Am_CS7_1/i;
X_Am_CS7_2 = X_Am_CS7_2/i;
save('X_Am_CS6_1_k=4.mat','X_Am_CS6_1')
save('X_Am_CS6_2_k=4.mat','X_Am_CS6_2')
clf
subplot(1,3,1);
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',X_Am_CS5_2','FaceColor','interp','LineStyle','none');
view([-59.3039,41.1722])
%colorbar
set(gca,'clim',[.2 .7] )
axis equal
axis off
camlight('headlight')
subplot(1,3,2);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',X_Am_CS6_1','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',X_Am_CS6_2','FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',X_Am_CS6_3','FaceColor','interp','LineStyle','none');
view([ 165.7876,21.7629])
set(gca,'clim',[.2 .7] )
axis equal
axis off
camlight('headlight')
subplot(1,3,3);
f1 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_Am_CS7_2','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ3.objects(16).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_Am_CS7_1','FaceColor','interp','LineStyle','none');
%colorbar
set(gca,'clim',[.2 .7] )
view([24.5784,1.2879])
axis equal
axis off
camlight('headlight')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 10])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/Am_All.png'])

listType = find(strcmp(Idents2.textdata(2:end,2),"Amnioid_#2")==1);
X_AmB_CS5_1 = ConvCS5{listType(1),1};X_AmB_CS5_2 = ConvCS5{listType(1),2};
for i = 2:length(listType)
X_AmB_CS5_1 = X_AmB_CS5_1 + ConvCS5{listType(i),1};
X_AmB_CS5_2 = X_AmB_CS5_2 + ConvCS5{listType(i),2};
end
X_AmB_CS5_1 = X_AmB_CS5_1/i;
X_AmB_CS5_2 = X_AmB_CS5_2/i;
X_AmB_CS6_1 = ConvCS6{listType(1),1};X_AmB_CS6_2 = ConvCS6{listType(1),2};X_AmB_CS6_3 = ConvCS6{listType(1),3};
for i = 2:length(listType)
X_AmB_CS6_1 = X_AmB_CS6_1 + ConvCS6{listType(i),1};
X_AmB_CS6_2 = X_AmB_CS6_2 + ConvCS6{listType(i),2};
X_AmB_CS6_3 = X_AmB_CS6_3 + ConvCS6{listType(i),3};
end
X_AmB_CS6_1 = X_AmB_CS6_1/i;
X_AmB_CS6_2 = X_AmB_CS6_2/i;
X_AmB_CS6_3 = X_AmB_CS6_3/i;
X_AmB_CS7_1 = ConvCS7{listType(1),1};X_AmB_CS7_2 = ConvCS7{listType(1),2};
for i = 2:length(listType)
X_AmB_CS7_1 = X_AmB_CS7_1 + ConvCS7{listType(i),1};
X_AmB_CS7_2 = X_AmB_CS7_2 + ConvCS7{listType(i),2};
end
X_AmB_CS7_1 = X_AmB_CS7_1/i;
X_AmB_CS7_2 = X_AmB_CS7_2/i;
save('X_AmB_CS6_1_k=4.mat','X_AmB_CS6_1')
save('X_AmB_CS6_2_k=4.mat','X_AmB_CS6_2')
clf
subplot(1,3,1);
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',X_AmB_CS5_2','FaceColor','interp','LineStyle','none');
view([-59.3039,41.1722])
%colorbar
set(gca,'clim',[.2 .7] )
axis equal
axis off
camlight('headlight')
subplot(1,3,2);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',X_AmB_CS6_1','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',X_AmB_CS6_2','FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',X_AmB_CS6_3','FaceColor','interp','LineStyle','none');
view([ 165.7876,21.7629])
set(gca,'clim',[.2 .7] )
axis equal
axis off
camlight('headlight')
subplot(1,3,3);
f1 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_AmB_CS7_2','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ3.objects(16).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_AmB_CS7_1','FaceColor','interp','LineStyle','none');
%colorbar
set(gca,'clim',[.2 .7] )
view([24.5784,1.2879])
axis equal
axis off
camlight('headlight')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 10])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/Amnioid_bead_All.png'])




%Am lineagesFGF_noMEF
listType = find(strcmp(Idents2.textdata(2:end,2),"FGF_noMEF")==1);
X_FGFNo_CS5_1 = ConvCS5{listType(1),1};X_FGFNo_CS5_2 = ConvCS5{listType(1),2};
for i = 2:length(listType)
X_FGFNo_CS5_1 = X_FGFNo_CS5_1 + ConvCS5{listType(i),1};
X_FGFNo_CS5_2 = X_FGFNo_CS5_2 + ConvCS5{listType(i),2};
end
X_FGFNo_CS5_1 = X_FGFNo_CS5_1/i;
X_FGFNo_CS5_2 = X_FGFNo_CS5_2/i;
X_FGFNo_CS6_1 = ConvCS6{listType(1),1};X_FGFNo_CS6_2 = ConvCS6{listType(1),2};X_FGFNo_CS6_3 = ConvCS6{listType(1),3};
for i = 2:length(listType)
X_FGFNo_CS6_1 = X_FGFNo_CS6_1 + ConvCS6{listType(i),1};
X_FGFNo_CS6_2 = X_FGFNo_CS6_2 + ConvCS6{listType(i),2};
X_FGFNo_CS6_3 = X_FGFNo_CS6_3 + ConvCS6{listType(i),3};
end
X_FGFNo_CS6_1 = X_FGFNo_CS6_1/i;
X_FGFNo_CS6_2 = X_FGFNo_CS6_2/i;
X_FGFNo_CS6_3 = X_FGFNo_CS6_3/i;
X_FGFNo_CS7_1 = ConvCS7{listType(1),1};X_FGFNo_CS7_2 = ConvCS7{listType(1),2};
for i = 2:length(listType)
X_FGFNo_CS7_1 = X_FGFNo_CS7_1 + ConvCS7{listType(i),1};
X_FGFNo_CS7_2 = X_FGFNo_CS7_2 + ConvCS7{listType(i),2};
end
X_FGFNo_CS7_1 = X_FGFNo_CS7_1/i;
X_FGFNo_CS7_2 = X_FGFNo_CS7_2/i;
save('X_FGFNo_CS6_1_k=4.mat','X_FGFNo_CS6_1')
save('X_FGFNo_CS6_2_k=4.mat','X_FGFNo_CS6_2')

clf
subplot(1,3,1);
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',X_FGFNo_CS5_2' - X_AmB_CS5_2','FaceColor','interp','LineStyle','none');
view([-59.3039,41.1722])
%colorbar
set(gca,'clim',[-0.1 .1] )
axis equal
axis off
camlight('headlight')
subplot(1,3,2);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',X_FGFNo_CS6_1' - X_AmB_CS6_1','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',X_FGFNo_CS6_2' - X_AmB_CS6_2','FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',X_FGFNo_CS6_3' - X_AmB_CS6_3','FaceColor','interp','LineStyle','none');
view([ 165.7876,21.7629])
set(gca,'clim',[-0.1 .1] )
%colorbar
axis equal
axis off
camlight('headlight')
subplot(1,3,3);
f1 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_FGFNo_CS7_2' - X_AmB_CS7_2','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ3.objects(16).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_FGFNo_CS7_1' - X_AmB_CS7_1','FaceColor','interp','LineStyle','none');
%colorbar
set(gca,'clim',[-0.1 .1] )
view([24.5784,1.2879])
axis equal
axis off
camlight('headlight')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 10])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/FGF_vs_AmB_All.png'])



listType = find(strcmp(Idents2.textdata(2:end,2),"ActA_noMEF")==1);
X_ActANo_CS5_1 = ConvCS5{listType(1),1};X_ActANo_CS5_2 = ConvCS5{listType(1),2};
for i = 2:length(listType)
X_ActANo_CS5_1 = X_ActANo_CS5_1 + ConvCS5{listType(i),1};
X_ActANo_CS5_2 = X_ActANo_CS5_2 + ConvCS5{listType(i),2};
end
X_ActANo_CS5_1 = X_ActANo_CS5_1/i;
X_ActANo_CS5_2 = X_ActANo_CS5_2/i;
X_ActANo_CS6_1 = ConvCS6{listType(1),1};X_ActANo_CS6_2 = ConvCS6{listType(1),2};X_ActANo_CS6_3 = ConvCS6{listType(1),3};
for i = 2:length(listType)
X_ActANo_CS6_1 = X_ActANo_CS6_1 + ConvCS6{listType(i),1};
X_ActANo_CS6_2 = X_ActANo_CS6_2 + ConvCS6{listType(i),2};
X_ActANo_CS6_3 = X_ActANo_CS6_3 + ConvCS6{listType(i),3};
end
X_ActANo_CS6_1 = X_ActANo_CS6_1/i;
X_ActANo_CS6_2 = X_ActANo_CS6_2/i;
X_ActANo_CS6_3 = X_ActANo_CS6_3/i;
X_ActANo_CS7_1 = ConvCS7{listType(1),1};X_ActANo_CS7_2 = ConvCS7{listType(1),2};
for i = 2:length(listType)
X_ActANo_CS7_1 = X_ActANo_CS7_1 + ConvCS7{listType(i),1};
X_ActANo_CS7_2 = X_ActANo_CS7_2 + ConvCS7{listType(i),2};
end
X_ActANo_CS7_1 = X_ActANo_CS7_1/i;
X_ActANo_CS7_2 = X_ActANo_CS7_2/i;
save('X_ActANo_CS6_1_k=4.mat','X_ActANo_CS6_1')
save('X_ActANo_CS6_2_k=4.mat','X_ActANo_CS6_2')



clf
subplot(1,3,1);
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',X_ActANo_CS5_2' - X_AmB_CS5_2','FaceColor','interp','LineStyle','none');
view([-59.3039,41.1722])
%colorbar
set(gca,'clim',[-0.1 .1] )
axis equal
axis off
camlight('headlight')
subplot(1,3,2);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',X_ActANo_CS6_1' - X_AmB_CS6_1','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',X_ActANo_CS6_2' - X_AmB_CS6_2','FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',X_ActANo_CS6_3' - X_AmB_CS6_3','FaceColor','interp','LineStyle','none');
view([ 165.7876,21.7629])
set(gca,'clim',[-0.1 .1] )
%colorbar
axis equal
axis off
camlight('headlight')
subplot(1,3,3);
f1 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_ActANo_CS7_2' - X_AmB_CS7_2','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ3.objects(16).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_ActANo_CS7_1' - X_AmB_CS7_1','FaceColor','interp','LineStyle','none');
%colorbar
set(gca,'clim',[-0.1 .1] )
view([24.5784,1.2879])
axis equal
axis off
camlight('headlight')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 10])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/ActANo_vs_AmB_All.png'])



listType = find(strcmp(Idents2.textdata(2:end,2),"BMP_noMEF")==1);
X_BMPNo_CS5_1 = ConvCS5{listType(1),1};X_BMPNo_CS5_2 = ConvCS5{listType(1),2};
for i = 2:length(listType)
X_BMPNo_CS5_1 = X_BMPNo_CS5_1 + ConvCS5{listType(i),1};
X_BMPNo_CS5_2 = X_BMPNo_CS5_2 + ConvCS5{listType(i),2};
end
X_BMPNo_CS5_1 = X_BMPNo_CS5_1/i;
X_BMPNo_CS5_2 = X_BMPNo_CS5_2/i;
X_BMPNo_CS6_1 = ConvCS6{listType(1),1};X_BMPNo_CS6_2 = ConvCS6{listType(1),2};X_BMPNo_CS6_3 = ConvCS6{listType(1),3};
for i = 2:length(listType)
X_BMPNo_CS6_1 = X_BMPNo_CS6_1 + ConvCS6{listType(i),1};
X_BMPNo_CS6_2 = X_BMPNo_CS6_2 + ConvCS6{listType(i),2};
X_BMPNo_CS6_3 = X_BMPNo_CS6_3 + ConvCS6{listType(i),3};
end
X_BMPNo_CS6_1 = X_BMPNo_CS6_1/i;
X_BMPNo_CS6_2 = X_BMPNo_CS6_2/i;
X_BMPNo_CS6_3 = X_BMPNo_CS6_3/i;
X_BMPNo_CS7_1 = ConvCS7{listType(1),1};X_BMPNo_CS7_2 = ConvCS7{listType(1),2};
for i = 2:length(listType)
X_BMPNo_CS7_1 = X_BMPNo_CS7_1 + ConvCS7{listType(i),1};
X_BMPNo_CS7_2 = X_BMPNo_CS7_2 + ConvCS7{listType(i),2};
end
X_BMPNo_CS7_1 = X_BMPNo_CS7_1/i;
X_BMPNo_CS7_2 = X_BMPNo_CS7_2/i;

save('X_BMPNo_CS6_1_k=4.mat','X_BMPNo_CS6_1')
save('X_BMPNo_CS6_2_k=4.mat','X_BMPNo_CS6_2')

clf
subplot(1,3,1);
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',X_BMPNo_CS5_2' - X_AmB_CS5_2','FaceColor','interp','LineStyle','none');
view([-59.3039,41.1722])
%colorbar
set(gca,'clim',[-0.1 .1] )
axis equal
axis off
camlight('headlight')
subplot(1,3,2);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',X_BMPNo_CS6_1' - X_AmB_CS6_1','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',X_BMPNo_CS6_2' - X_AmB_CS6_2','FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',X_BMPNo_CS6_3' - X_AmB_CS6_3','FaceColor','interp','LineStyle','none');
view([ 165.7876,21.7629])
set(gca,'clim',[-0.1 .1] )
%colorbar
axis equal
axis off
camlight('headlight')
subplot(1,3,3);
f1 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_BMPNo_CS7_2' - X_AmB_CS7_2','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ3.objects(16).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_BMPNo_CS7_1' - X_AmB_CS7_1','FaceColor','interp','LineStyle','none');
%colorbar
set(gca,'clim',[-0.1 .1] )
view([24.5784,1.2879])
axis equal
axis off
camlight('headlight')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 10])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/BMP_vs_AmB_All.png'])



listType = find(strcmp(Idents2.textdata(2:end,2),"ActA_MEF")==1);
X_ActA_CS5_1 = ConvCS5{listType(1),1};X_ActA_CS5_2 = ConvCS5{listType(1),2};
for i = 2:length(listType)
X_ActA_CS5_1 = X_ActA_CS5_1 + ConvCS5{listType(i),1};
X_ActA_CS5_2 = X_ActA_CS5_2 + ConvCS5{listType(i),2};
end
X_ActA_CS5_1 = X_ActA_CS5_1/i;
X_ActA_CS5_2 = X_ActA_CS5_2/i;
X_ActA_CS6_1 = ConvCS6{listType(1),1};X_ActA_CS6_2 = ConvCS6{listType(1),2};X_ActA_CS6_3 = ConvCS6{listType(1),3};
for i = 2:length(listType)
X_ActA_CS6_1 = X_ActA_CS6_1 + ConvCS6{listType(i),1};
X_ActA_CS6_2 = X_ActA_CS6_2 + ConvCS6{listType(i),2};
X_ActA_CS6_3 = X_ActA_CS6_3 + ConvCS6{listType(i),3};
end
X_ActA_CS6_1 = X_ActA_CS6_1/i;
X_ActA_CS6_2 = X_ActA_CS6_2/i;
X_ActA_CS6_3 = X_ActA_CS6_3/i;
X_ActA_CS7_1 = ConvCS7{listType(1),1};X_ActA_CS7_2 = ConvCS7{listType(1),2};
for i = 2:length(listType)
X_ActA_CS7_1 = X_ActA_CS7_1 + ConvCS7{listType(i),1};
X_ActA_CS7_2 = X_ActA_CS7_2 + ConvCS7{listType(i),2};
end
X_ActA_CS7_1 = X_ActA_CS7_1/i;
X_ActA_CS7_2 = X_ActA_CS7_2/i;
save('X_ActA_CS6_1_k=4.mat','X_ActA_CS6_1')
save('X_ActA_CS6_2_k=4.mat','X_ActA_CS6_2')



%Emb lineages

%Emb lineages
listType = find(strcmp(Idents2.textdata(2:end,2),"BMP_MEF")==1);
X_BMP_CS5_1 = ConvCS5{listType(1),1};X_BMP_CS5_2 = ConvCS5{listType(1),2};
for i = 2:length(listType)
X_BMP_CS5_1 = X_BMP_CS5_1 + ConvCS5{listType(i),1};
X_BMP_CS5_2 = X_BMP_CS5_2 + ConvCS5{listType(i),2};
end
X_BMP_CS5_1 = X_BMP_CS5_1/i;
X_BMP_CS5_2 = X_BMP_CS5_2/i;
X_BMP_CS6_1 = ConvCS6{listType(1),1};X_BMP_CS6_2 = ConvCS6{listType(1),2};X_BMP_CS6_3 = ConvCS6{listType(1),3};
for i = 2:length(listType)
X_BMP_CS6_1 = X_BMP_CS6_1 + ConvCS6{listType(i),1};
X_BMP_CS6_2 = X_BMP_CS6_2 + ConvCS6{listType(i),2};
X_BMP_CS6_3 = X_BMP_CS6_3 + ConvCS6{listType(i),3};
end
X_BMP_CS6_1 = X_BMP_CS6_1/i;
X_BMP_CS6_2 = X_BMP_CS6_2/i;
X_BMP_CS6_3 = X_BMP_CS6_3/i;
X_BMP_CS7_1 = ConvCS7{listType(1),1};X_BMP_CS7_2 = ConvCS7{listType(1),2};
for i = 2:length(listType)
X_BMP_CS7_1 = X_BMP_CS7_1 + ConvCS7{listType(i),1};
X_BMP_CS7_2 = X_BMP_CS7_2 + ConvCS7{listType(i),2};
end
X_BMP_CS7_1 = X_BMP_CS7_1/i;
X_BMP_CS7_2 = X_BMP_CS7_2/i;
clf
save('X_BMP_CS6_1_k=4.mat','X_BMP_CS6_1')
save('X_BMP_CS6_2_k=4.mat','X_BMP_CS6_2')
subplot(1,3,1);
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',X_BMP_CS5_2' - X_EmDisc_CS5_2','FaceColor','interp','LineStyle','none');
view([-59.3039,41.1722])
%colorbar
set(gca,'clim',[-0.1 .1] )
axis equal
axis off
camlight('headlight')
subplot(1,3,2);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',X_BMP_CS6_1' - X_EmDisc_CS6_1','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',X_BMP_CS6_2' - X_EmDisc_CS6_2','FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',X_BMP_CS6_3' - X_EmDisc_CS6_3','FaceColor','interp','LineStyle','none');
view([ 165.7876,21.7629])
set(gca,'clim',[-0.1 .1] )
%colorbar
axis equal
axis off
camlight('headlight')
subplot(1,3,3);
f1 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_BMP_CS7_2' - X_EmDisc_CS7_2','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ3.objects(16).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_BMP_CS7_1' - X_EmDisc_CS7_1','FaceColor','interp','LineStyle','none');
%colorbar
set(gca,'clim',[-0.1 .1] )
view([24.5784,1.2879])
axis equal
axis off
camlight('headlight')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 10])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/BMP_vs_EmDisc_All.png'])




listType = find(strcmp(Idents2.textdata(2:end,2),"SB43_MEF")==1);
X_SB43_CS5_1 = ConvCS5{listType(1),1};X_SB43_CS5_2 = ConvCS5{listType(1),2};
for i = 2:length(listType)
X_SB43_CS5_1 = X_SB43_CS5_1 + ConvCS5{listType(i),1};
X_SB43_CS5_2 = X_SB43_CS5_2 + ConvCS5{listType(i),2};
end
X_SB43_CS5_1 = X_SB43_CS5_1/i;
X_SB43_CS5_2 = X_SB43_CS5_2/i;
X_SB43_CS6_1 = ConvCS6{listType(1),1};X_SB43_CS6_2 = ConvCS6{listType(1),2};X_SB43_CS6_3 = ConvCS6{listType(1),3};
for i = 2:length(listType)
X_SB43_CS6_1 = X_SB43_CS6_1 + ConvCS6{listType(i),1};
X_SB43_CS6_2 = X_SB43_CS6_2 + ConvCS6{listType(i),2};
X_SB43_CS6_3 = X_SB43_CS6_3 + ConvCS6{listType(i),3};
end
X_SB43_CS6_1 = X_SB43_CS6_1/i;
X_SB43_CS6_2 = X_SB43_CS6_2/i;
X_SB43_CS6_3 = X_SB43_CS6_3/i;
X_SB43_CS7_1 = ConvCS7{listType(1),1};X_SB43_CS7_2 = ConvCS7{listType(1),2};
for i = 2:length(listType)
X_SB43_CS7_1 = X_SB43_CS7_1 + ConvCS7{listType(i),1};
X_SB43_CS7_2 = X_SB43_CS7_2 + ConvCS7{listType(i),2};
end
X_SB43_CS7_1 = X_SB43_CS7_1/i;
X_SB43_CS7_2 = X_SB43_CS7_2/i;
save('X_SB43_CS6_1_k=4.mat','X_SB43_CS6_1')
save('X_SB43_CS6_2_k=4.mat','X_SB43_CS6_2')

clf
subplot(1,3,1);
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',X_SB43_CS5_2' - X_EmDisc_CS5_2','FaceColor','interp','LineStyle','none');
view([-59.3039,41.1722])
%colorbar
set(gca,'clim',[-0.1 .1] )
axis equal
axis off
camlight('headlight')
subplot(1,3,2);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',X_SB43_CS6_1' - X_EmDisc_CS6_1','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',X_SB43_CS6_2' - X_EmDisc_CS6_2','FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',X_SB43_CS6_3' - X_EmDisc_CS6_3','FaceColor','interp','LineStyle','none');
view([ 165.7876,21.7629])
set(gca,'clim',[-0.1 .1] )
%colorbar
axis equal
axis off
camlight('headlight')
subplot(1,3,3);
f1 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_SB43_CS7_2' - X_EmDisc_CS7_2','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ3.objects(16).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_SB43_CS7_1' - X_EmDisc_CS7_1','FaceColor','interp','LineStyle','none');
%colorbar
set(gca,'clim',[-0.1 .1] )
view([24.5784,1.2879])
axis equal
axis off
camlight('headlight')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 10])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/SB43_vs_EmDisc_All.png'])



listType = find(strcmp(Idents2.textdata(2:end,2),"CHIR_MEF")==1);
X_CHIR_CS5_1 = ConvCS5{listType(1),1};X_CHIR_CS5_2 = ConvCS5{listType(1),2};
for i = 2:length(listType)
X_CHIR_CS5_1 = X_CHIR_CS5_1 + ConvCS5{listType(i),1};
X_CHIR_CS5_2 = X_CHIR_CS5_2 + ConvCS5{listType(i),2};
end
X_CHIR_CS5_1 = X_CHIR_CS5_1/i;
X_CHIR_CS5_2 = X_CHIR_CS5_2/i;
X_CHIR_CS6_1 = ConvCS6{listType(1),1};X_CHIR_CS6_2 = ConvCS6{listType(1),2};X_CHIR_CS6_3 = ConvCS6{listType(1),3};
for i = 2:length(listType)
X_CHIR_CS6_1 = X_CHIR_CS6_1 + ConvCS6{listType(i),1};
X_CHIR_CS6_2 = X_CHIR_CS6_2 + ConvCS6{listType(i),2};
X_CHIR_CS6_3 = X_CHIR_CS6_3 + ConvCS6{listType(i),3};
end
X_CHIR_CS6_1 = X_CHIR_CS6_1/i;
X_CHIR_CS6_2 = X_CHIR_CS6_2/i;
X_CHIR_CS6_3 = X_CHIR_CS6_3/i;
X_CHIR_CS7_1 = ConvCS7{listType(1),1};X_CHIR_CS7_2 = ConvCS7{listType(1),2};
for i = 2:length(listType)
X_CHIR_CS7_1 = X_CHIR_CS7_1 + ConvCS7{listType(i),1};
X_CHIR_CS7_2 = X_CHIR_CS7_2 + ConvCS7{listType(i),2};
end
X_CHIR_CS7_1 = X_CHIR_CS7_1/i;
X_CHIR_CS7_2 = X_CHIR_CS7_2/i;
save('X_CHIR_CS6_1_k=4.mat','X_CHIR_CS6_1')
save('X_CHIR_CS6_2_k=4.mat','X_CHIR_CS6_2')

clf
subplot(1,3,1);
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',X_CHIR_CS5_2' - X_EmDisc_CS5_2','FaceColor','interp','LineStyle','none');
view([-59.3039,41.1722])
%colorbar
set(gca,'clim',[-0.1 .1] )
axis equal
axis off
camlight('headlight')
subplot(1,3,2);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',X_CHIR_CS6_1' - X_EmDisc_CS6_1','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',X_CHIR_CS6_2' - X_EmDisc_CS6_2','FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',X_CHIR_CS6_3' - X_EmDisc_CS6_3','FaceColor','interp','LineStyle','none');
view([ 165.7876,21.7629])
set(gca,'clim',[-0.1 .1] )
%colorbar
axis equal
axis off
camlight('headlight')
subplot(1,3,3);
f1 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_CHIR_CS7_2' - X_EmDisc_CS7_2','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ3.objects(16).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_CHIR_CS7_1' - X_EmDisc_CS7_1','FaceColor','interp','LineStyle','none');
%colorbar
set(gca,'clim',[-0.1 .1] )
view([24.5784,1.2879])
axis equal
axis off
camlight('headlight')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 10])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/CHIR_vs_EmDisc_All.png'])


return

D1 = importdata('/Users/christopherpenfold/Desktop/Thorsten/ESCs_3DPCAcoordinates.csv');
 Distancematrix = dist(D1.data');
 IDs = D1.textdata(2:end,1:3);
 esc = find(strcmp(IDs(:,2),'ESC_primed')==1);
 naive = find(strcmp(IDs(:,2),'ESC_naive')==1);

  n1 = find(strcmp(IDs(:,2),'EmDisc_CS5')==1);
  n2 = find(strcmp(IDs(:,2),'EmDisc_CS6')==1);


  
 for i = 1:length(CellIDs)
     IndsMatch(i,1) = find(strcmp( ['X' CellIDs{i}],D1.textdata(2:end,1))==1);
 end
 Dhat1 = ( max(max(Distancematrix)) - mean(Distancematrix(IndsMatch,esc),2) ) ./max(max(Distancematrix));
 Dhat2 = ( max(max(Distancematrix)) - mean(Distancematrix(IndsMatch,naive),2) ) ./max(max(Distancematrix));
 Dhat3 = ( max(max(Distancematrix)) - min(Distancematrix(IndsMatch,esc)')' ) ./max(max(Distancematrix));
 Dhat4 = ( max(max(Distancematrix)) - min(Distancematrix(IndsMatch,naive)')' ) ./max(max(Distancematrix)); 
 Dhat5 = ( max(max(Distancematrix)) - Distancematrix(IndsMatch,esc) ) ./max(max(Distancematrix));
 Dhat6 = ( max(max(Distancematrix)) - Distancematrix(IndsMatch,naive) ) ./max(max(Distancematrix)); 
 

% D1 = importdata('/Users/christopherpenfold/Desktop/Thorsten/ESCs_3DPCAcoordinates.csv');
% Distancematrix = dist(D1.data');
% IDs = D1.textdata(2:end,1:3);
% 
% esc = find(strcmp(IDs(:,2),'ESC_primed')==1);
% naive = find(strcmp(IDs(:,2),'ESC_naive')==1);
% 
% INVIVO = setdiff(1:1:length(IDs),[esc;naive]);
% 
% for i = 1:length(IDs)    
%     type = IDs(i,2);
%     %indexclude = find(strcmp(IDs(:,2),type)==1);    
%     %type = IDs(find(Distancematrix(i,:)<15),2);
%     Closest = find(Distancematrix(i,:)<10);
%     N1(i,1) = length(intersect(Closest,esc));    
%     N2(i,1) = length(intersect(Closest,naive));
%     N3(i,1) = length(Closest);    
% 
%     Closest = find(Distancematrix(i,:)<12);
%     N1(i,2) = length(intersect(Closest,esc));    
%     N2(i,2) = length(intersect(Closest,naive));
%     N3(i,2) = length(Closest);    
%   
%     Closest = find(Distancematrix(i,:)<15);
%     N1(i,3) = length(intersect(Closest,esc));   
%     N2(i,3) = length(intersect(Closest,naive));
%     N3(i,3) = length(Closest);    
%     
%     Closest = find(Distancematrix(i,:)<18);
%     N1(i,4) = length(intersect(Closest,esc));   
%     N2(i,4) = length(intersect(Closest,naive));
%     N3(i,4) = length(Closest);    
% 
%     Closest = find(Distancematrix(i,:)<20);
%     N1(i,5) = length(intersect(Closest,esc));   
%     N2(i,5) = length(intersect(Closest,naive));
%     N3(i,5) = length(Closest);        
% end
% 
% %imagesc([N1./N3,N2./N3])
% %set(gca,'YTick',1:1:size(IDs,1),'YTickLabel',IDs(:,2))
% 
% 
% %DMAP1 = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/CCAFigForESCs/DimRed/ESC_all.csv')
% %targ = DMAP1.textdata(1,2:end);
% %targ = strrep(targ,'"','');
% %targ = strrep(targ,'RNA.','');
% 
% %estarg = DMAP1.textdata(2:end,1);
% %estarg = strrep(estarg,'RNA.','');
% 
% %Targ = targ;
% %Targ = strrep(Targ,'EmDisc_CS5_','');
% %Targ = strrep(Targ,'ExMes_CS5_','');
% %Targ = strrep(Targ,'VE_CS5_','');
% %Targ = strrep(Targ,'Tb_CS5_','');
% %Targ = strrep(Targ,'Am_CS5_','');
% for i = 1:length(CellIDs)
%     try
%     P1(i,:) = N1(find(strcmp(IDs(:,1),[CellIDs{i}])==1),:)./N3(find(strcmp(IDs(:,1),[CellIDs{i}])==1),:);
%     P2(i,:) = N2(find(strcmp(IDs(:,1),[CellIDs{i}])==1),:)./N3(find(strcmp(IDs(:,1),[CellIDs{i}])==1),:);
%     catch
%             P1(i,:) = N1(find(strcmp(IDs(:,1),['X',CellIDs{i}])==1),:)./N3(find(strcmp(IDs(:,1),['X',CellIDs{i}])==1),:);
%     P2(i,:) = N2(find(strcmp(IDs(:,1),['X',CellIDs{i}])==1),:)./N3(find(strcmp(IDs(:,1),['X',CellIDs{i}])==1),:);
%     end
% end    
%     

set(0, 'DefaultTextInterpreter', 'none')
%for i = 1:size(P1,1)

    
%     Ytrain2 = P2(find(LabBin==1),3);
%     Ytrain = P1(find(LabBin==1),3);
%     %Ytrain = (Pt(find(isnan(Pt(:,i))~=1),i) );
%     Xtrain = Output.Xtrain(find(LabBin==1),:);
 % 
 
 
  pg1 = {@priorGauss,0,2};  
  pg2 = {@priorGauss,0,2};  
  pg3 = {@priorGauss,log(0.5),2};  
  pc = {@priorClamped};
  prior1.mean = {[]};  
  prior1.cov  = {pg1;pg2;pg2;pg2};
  prior1.lik  = {pg3};
  im1 = {@infPrior,@infExact,prior1};
%  pg1 = {@priorGauss,0,1};  pg2 = {@priorGauss,0,1};  pg3 = {@priorGauss,log(0.5),1};  pc = {@priorClamped};
%  prior1.mean = {[]}; prior1.cov  = {pg1;pg2}; prior1.lik  = {pg3};
%  im1 = {@infPrior,@infExact,prior1};                % inference method
%  
%  hyp1.cov  = [var(Ytrain); log(2)];
%  hyp1.lik  = [var(Ytrain)/2];    
%  hyp1.mean = mean(Ytrain);       
% 
% par  = {@meanConst,@covSEiso,'likGauss',Xtrain(:,1:3), Ytrain};
% hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
% 
% %[m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(:,1:3), Ytrain, Xtrain(:,1:3));
% [m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(:,1:3), Ytrain, Output.Xtest(Output.ind1,:));
%  hyp1.cov  = [var(Ytrain2); log(2)];
%  hyp1.lik  = [var(Ytrain2)/2];    
%  hyp1.mean = mean(Ytrain2);       
% 
% [m_2 s_2] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(:,1:3), Ytrain2, Output.Xtest(Output.ind1,:));


%scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,[0,0,0],'fill')
%hold on
%scatter3(Xtrain(:,1),Xtrain(:,2),Xtrain(:,3),250,Ytrain,'fill')


%h = figure(1)
% subplot(1,2,1);
% scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,m_1,'fill')
% hold on
% scatter3(Xtrain(:,1),Xtrain(:,2),Xtrain(:,3),250,Ytrain,'fill')
% %t%itle([estarg{i}])
% view([ 48, 18])
% alpha = 0.05;
% %set(h1,'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha)
% 
% yl = ylim;
% xl = xlim;
% zl = zlim;
% set(gca,'clim',[0 0.2] )
% 
% set(gca,'ZTick',linspace(zl(1),zl(2),5),'YTick',linspace(yl(1),yl(2),5),'XTick',linspace(xl(1),xl(2),5),'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))
% set(gca,'fontsize', 12)
% set(gca,'TickLength',[0 0])
% cb=colorbar;
% %cb.Position = cb.Position + [0.1 -0.02 0 0]
% %cy=get(cb,'YTick');
% %set(cb,'YTick',[0 0.1)])
% set(gca,'fontsize', 12)
% %axis equal
% 
% subplot(1,2,2);
% %h = figure(2)
% scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,m_2,'fill')
% hold on
% scatter3(Xtrain(:,1),Xtrain(:,2),Xtrain(:,3),250,Ytrain2,'fill')
% %title([estarg{i}])
% view([ 48, 18])
% alpha = 0.05;
% %set(h1,'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha)
% 
% yl = ylim;
% xl = xlim;
% zl = zlim;
% set(gca,'clim',[0 0.2] )
% 
% set(gca,'ZTick',linspace(zl(1),zl(2),5),'YTick',linspace(yl(1),yl(2),5),'XTick',linspace(xl(1),xl(2),5),'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))
% set(gca,'fontsize', 12)
% set(gca,'TickLength',[0 0])
% cb=colorbar;
% %cb.Position = cb.Position + [0.1 -0.02 0 0]
% %cy=get(cb,'YTick');
% %set(cb,'YTick',[0 0.1)])
% set(gca,'fontsize', 12)
% axis equal
% 
% clf

% 
% 
% Dcorr = importdata("~/Desktop/Thorsten/FINAL/CCAFigForPaper_wCS6_v2//DimRed/ESC_corr_wCS6.csv")
% 
% 
% targs = Dcorr.textdata(1,2:end);
% targs = strrep(targs,'"','');
% targs = strrep(targs,'RNA.','');
% escs = Dcorr.textdata(2:end,1);
% 
% CoC = zeros(length(escs),length(CellIDs4));
% for i = 1:length(CellIDs4)
%     idnss = find(contains(targs,CellIDs4{i}));
%     try
%     CoC(:,i) = Dcorr.data(:,idnss);
%     catch
%     CoC(:,1) =0;
%     end
% end
% 
% 
% 
% %   
% % 
% % for i = 1:length(escs)
% %     
% %     train = CoC(i,find(LabBin==1))';
% %     Xtrain = Output.Xtrain(find(LabBin==1),:);
% %     
% %     hyp1.cov  = [var(Ytrain); log(2)]; hyp1.lik  = [var(Ytrain)/2];  hyp1.mean = mean(Ytrain);       
% % 
% %     par  = {@meanConst,@covSEiso,'likGauss',Xtrain(:,1:3), Ytrain};
% %     hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
% % 
% %     [m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(:,1:3), Ytrain, Output.Xtest(Output.ind1,:));
% %     
% %     
% %     h = figure(1)
% %     scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,m_1,'fill')
% %     hold on
% %     scatter3(Xtrain(:,1),Xtrain(:,2),Xtrain(:,3),250,Ytrain,'fill','MarkerEdgeColor','k')
% %     title([escs{i}])
% %     view([ -78.8,14.9])
% %     alpha = 0.05;
% %     yl = ylim;
% %     xl = xlim;
% %     zl = zlim;
% %     set(gca,'clim',[0.2 0.9] )
% %     set(gca,'ZTick',linspace(zl(1),zl(2),5),'YTick',linspace(yl(1),yl(2),5),'XTick',linspace(xl(1),xl(2),5),'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))
% %     set(gca,'fontsize', 12)
% %     set(gca,'TickLength',[0 0])
% %     cb=colorbar;
% %     set(gca,'fontsize', 12)
% % 
% %    %pause 
% %     print(h,['~/Desktop/Sections/CS6_' escs{i} '.pdf'],'-dpdf','-painters');
% %     clf
% % end
% 
% %print(h, '-dpdf', ['~/Desktop/' estarg{i} '.pdf'],'-r500','-painters')
% %clf
% %end
% 
% %plot(m_1,Ytrain,'o')
% 
% %[L1 ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(:,1:3), Ytrain');
% %[L2 ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(:,1:3), Ytrain');
% 
% %for i = 1:25
% %subplot(5,5,i);scatter3(Output.Xtrain(:,1),Output.Xtrain(:,2),Output.Xtrain(:,3),50,Pt(:,i),'fill')
% %
% %end
% % h = figure(5)
% % scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,m_1,'fill')
% % hold on
% % scatter3(Output.Xtrain(find(Output.LabBin==5),1),Output.Xtrain(find(Output.LabBin==5),2),Output.Xtrain(find(Output.LabBin==5),3),95,'fill')
% % text(Output.Xtrain(find(Output.LabBin==5),1),Output.Xtrain(find(Output.LabBin==5),2),Output.Xtrain(find(Output.LabBin==5),3), CellIDs4(find(Output.LabBin==5)));
% % savefig(h,['~/Desktop/FIG/CS5_Tb'])
% % clf
% % 
% 
% %escs(1:45)
% %escs(46:61)
% 
% 
%   Ytrain1 = reshape(CoC([1:45],find(LabBin==1))',45*46,1);
%   Ytrain2 = reshape(CoC([46:61],find(LabBin==1))',16*46,1);  
%   
%   Xtrain1 = repmat(Output.Xtrain(find(LabBin==1),:),45,1);
%   Xtrain2 = repmat(Output.Xtrain(find(LabBin==1),:),16,1);   
%       
%   hyp1.cov  = [var(Ytrain1); log(2)]; hyp1.lik  = [var(Ytrain1)/2];  hyp1.mean = mean(Ytrain1);       
%   par  = {@meanConst,@covSEiso,'likGauss',Xtrain(:,1:3), Ytrain1};
%   hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
%   [m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain1(:,1:3), Ytrain1, [OBJ2.vertices]/400);
%         
%   Output.Conv = m_1;
% 
% %    h = figure(2)
% %    subplot(2,3,1);
% %    scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,m_1,'fill')
% %    hold on
% %    %scatter3(Xtrain(:,1),Xtrain(:,2),Xtrain(:,3),250,Ytrain,'fill','MarkerEdgeColor','k')
% %     title(['Conv'])
% %     view([-121.7479,14.1676])
% %     alpha = 0.05;
% %     yl = ylim;
% %     xl = xlim;
% %     zl = zlim;
% %     set(gca,'clim',[0.4 0.7] )
% %     set(gca,'ZTick',linspace(zl(1),zl(2),5),'YTick',linspace(yl(1),yl(2),5),'XTick',linspace(xl(1),xl(2),5),'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))
% %     set(gca,'fontsize', 12)
% %     set(gca,'TickLength',[0 0])
% %     cb=colorbar;
% %     set(gca,'fontsize', 12)
% %     subplot(2,3,4);
% %     scatter3(Xtrain1(:,1)+randn(2070,1)*0.05,Xtrain1(:,2)+randn(2070,1)*0.05,Xtrain1(:,3)+randn(2070,1)*0.05,250,Ytrain1,'fill','MarkerEdgeColor','k')
% %     view([-121.7479,14.1676])
% %     yl = ylim;
% %     xl = xlim;
% %     zl = zlim;
% %     set(gca,'clim',[0.4 0.7] )
% %     set(gca,'ZTick',linspace(zl(1),zl(2),5),'YTick',linspace(yl(1),yl(2),5),'XTick',linspace(xl(1),xl(2),5),'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))
% %     set(gca,'fontsize', 12)
% %     set(gca,'TickLength',[0 0])
% %     cb=colorbar;
% %     set(gca,'fontsize', 12)
%    %pause 
%     %print(h,['~/Desktop/Sections/CS6_' escs{i} '.pdf'],'-dpdf','-painters');
%     %clf
%     
%     hyp1.cov  = [var(Ytrain2); log(2)]; hyp1.lik  = [var(Ytrain2)/2];  hyp1.mean = mean(Ytrain2);       
%     par  = {@meanConst,@covSEiso,'likGauss',Xtrain(:,1:3), Ytrain2};
%     hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
%     [m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain2(:,1:3), Ytrain2, [OBJ2.vertices]/400);
%     
%     
%     %h = figure(1)
% %     subplot(2,3,2);
% %     scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,m_1,'fill')
% %     hold on
% %     %scatter3(Xtrain(:,1),Xtrain(:,2),Xtrain(:,3),250,Ytrain,'fill','MarkerEdgeColor','k')
% %     title(['Naive'])
% %     view([-121.7479,14.1676])
% %     alpha = 0.05;
% %     yl = ylim;
% %     xl = xlim;
% %     zl = zlim;
% %     set(gca,'clim',[0.4 0.7] )
% %     set(gca,'ZTick',linspace(zl(1),zl(2),5),'YTick',linspace(yl(1),yl(2),5),'XTick',linspace(xl(1),xl(2),5),'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))
% %     set(gca,'fontsize', 12)
% %     set(gca,'TickLength',[0 0])
% %     cb=colorbar;
% %     set(gca,'fontsize', 12)
% % 
% %     
% %     subplot(2,3,5); scatter3(Xtrain2(:,1)+randn(736,1)*0.05,Xtrain2(:,2)+randn(736,1)*0.05,Xtrain2(:,3)+randn(736,1)*0.05,250,Ytrain2,'fill','MarkerEdgeColor','k')
% %     view([-121.7479,14.1676])
% %     yl = ylim;
% %     xl = xlim;
% %     zl = zlim;
% %     set(gca,'clim',[0.4 0.7] )
% %     set(gca,'ZTick',linspace(zl(1),zl(2),5),'YTick',linspace(yl(1),yl(2),5),'XTick',linspace(xl(1),xl(2),5),'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))
% %     set(gca,'fontsize', 12)
% %     set(gca,'TickLength',[0 0])
% %     cb=colorbar;
% %     set(gca,'fontsize', 12)
%     
% Output.Naive = m_1;
% 
%     %print(h,['~/Desktop/Sections/CS6_' escs{i} '.pdf'],'-dpdf','-painters');
%     %clf    




set(0, 'DefaultTextInterpreter', 'none')

%pg1 = {@priorGauss,0,2};  
%pg2 = {@priorGauss,0,2};  
%pg3 = {@priorGauss,log(0.5),2};  
%pc = {@priorClamped};
%prior1.mean = {[]};  
%prior1.cov  = {pg1;pg2;pg2;pg2};
%prior1.lik  = {pg3};
%im1 = {@infPrior,@infExact,prior1};
pg1 = {@priorGauss,0,1};  pg2 = {@priorGauss,0,1};  pg3 = {@priorGauss,log(0.5),1};  pc = {@priorClamped};
prior1.mean = {[]}; prior1.cov  = {pg1;pg2}; prior1.lik  = {pg3};
im1 = {@infPrior,@infExact,prior1};                % inference method

Dcorr = importdata("~/Desktop/Thorsten/FINAL/CCAFigForPaper_wCS6_v2//DimRed/ESC_corr_wCS6.csv");


targs = Dcorr.textdata(1,2:end);
targs = strrep(targs,'"','');
targs = strrep(targs,'RNA.','');
escs = Dcorr.textdata(2:end,1);

CoC = zeros(length(escs),length(CellIDs4));
for i = 1:length(CellIDs4)
    idnss = find(contains(targs,CellIDs4{i}));
    try
    CoC(:,i) = Dcorr.data(:,idnss);
    catch
    CoC(:,1) =0;
    end
end

Dhat5 = Dhat5';
Dhat6 = Dhat6';

for i = 1:size(CoC)
    weightedDist(i,1:3) = sum(repmat(CoC(i,:)',1,3).*Xtrain,1);
end


%keyboard
uLab = unique(LabBin);
for i = 1:length(uLab)
    
    
  Ytrain1 = reshape(CoC([1:45],find(LabBin==uLab(i)))',45*length(find(LabBin==uLab(i))),1);
  Ytrain2 = reshape(CoC([46:61],find(LabBin==uLab(i)))',16*length(find(LabBin==uLab(i))),1);  
  
  Xtrain1 = repmat(Output.Xtrain(find(LabBin==uLab(i)),:),45,1);
  Xtrain1 = Xtrain1 + randn(size(Xtrain1,1),3)*0.05;
  Xtrain2 = repmat(Output.Xtrain(find(LabBin==uLab(i)),:),16,1);
  Xtrain2 = Xtrain2 + randn(size(Xtrain2,1),3)*0.05;   
      
  hyp1.cov  = [log(2); var(Ytrain1)]; hyp1.lik  = [var(Ytrain2)/2];  hyp1.mean = mean(Ytrain1);       
  par  = {@meanConst,@covSEiso,'likGauss',Xtrain1(:,1:3), Ytrain1};
  hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
  [m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain1(:,1:3), Ytrain1, [OBJ2.vertices]/400);
        
       
  H{i,1} = hyp2;
  Output.Conv{i,1} = m_1;

    hyp1.cov  = [log(2); var(Ytrain2)]; hyp1.lik  = [var(Ytrain2)/2];  hyp1.mean = mean(Ytrain2);       
    par  = {@meanConst,@covSEiso,'likGauss',Xtrain2(:,1:3), Ytrain2};
    hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
    [m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain2(:,1:3), Ytrain2, [OBJ2.vertices]/400);
    
    
    Output.Naive{i,1} = m_1;
    
      H{i,2} = hyp2;
      
      
  
 
Ytrain3 = Dhat1(find(LabBin==uLab(i)));
Ytrain4 = Dhat2(find(LabBin==uLab(i)));
Ytrain5 = Dhat3(find(LabBin==uLab(i)));
Ytrain6 = Dhat4(find(LabBin==uLab(i)));
Xtrain3 = Output.Xtrain(find(LabBin==uLab(i)),:);
Xtrain4 = Xtrain3;
Xtrain5 = Xtrain3;
Xtrain6 = Xtrain3;

hyp1.cov  = [log(2); var(Ytrain3)]; hyp1.lik  = [var(Ytrain3)/2];  hyp1.mean = mean(Ytrain3);       
par  = {@meanConst,@covSEiso,'likGauss',Xtrain1(:,1:3), Ytrain3};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
[m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain3(:,1:3), Ytrain3, [OBJ2.vertices]/400);

 Output.Conv{i,2} = m_1;
      

hyp1.cov  = [log(2); var(Ytrain4)]; hyp1.lik  = [var(Ytrain4)/2];  hyp1.mean = mean(Ytrain4);       
par  = {@meanConst,@covSEiso,'likGauss',Xtrain1(:,1:3), Ytrain4};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
[m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain4(:,1:3), Ytrain4, [OBJ2.vertices]/400);

Output.Naive{i,2} = m_1;

hyp1.cov  = [log(2); var(Ytrain5)]; hyp1.lik  = [var(Ytrain5)/2];  hyp1.mean = mean(Ytrain5);       
par  = {@meanConst,@covSEiso,'likGauss',Xtrain1(:,1:3), Ytrain5};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
[m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain5(:,1:3), Ytrain5, [OBJ2.vertices]/400);

Output.Conv{i,3} = m_1;


hyp1.cov  = [log(2); var(Ytrain6)]; hyp1.lik  = [var(Ytrain6)/2];  hyp1.mean = mean(Ytrain6);       
par  = {@meanConst,@covSEiso,'likGauss',Xtrain1(:,1:3), Ytrain6};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
[m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain6(:,1:3), Ytrain6, [OBJ2.vertices]/400);

Output.Naive{i,3} = m_1;
    

  Ytrain7 = reshape(Dhat5(:,find(LabBin==uLab(i)))',45*length(find(LabBin==uLab(i))),1);
  Ytrain8 = reshape(Dhat6(:,find(LabBin==uLab(i)))',16*length(find(LabBin==uLab(i))),1);  

  
  
hyp1.cov  = [log(2); var(Ytrain7)]; hyp1.lik  = [var(Ytrain7)/2];  hyp1.mean = mean(Ytrain5);       
par  = {@meanConst,@covSEiso,'likGauss',Xtrain1(:,1:3), Ytrain7};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
[m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain1(:,1:3), Ytrain7, [OBJ2.vertices]/400);

Output.Conv{i,4} = m_1;


hyp1.cov  = [log(2); var(Ytrain8)]; hyp1.lik  = [var(Ytrain8)/2];  hyp1.mean = mean(Ytrain8);       
par  = {@meanConst,@covSEiso,'likGauss',Xtrain1(:,1:3), Ytrain8};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
[m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain2(:,1:3), Ytrain8, [OBJ2.vertices]/400);
 
Output.Naive{i,4} = m_1;

      
      %f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',m_1,'FaceColor','interp','LineStyle','none');
%hold on

end

Output.weightedDist = weightedDist;
Output.CoC = CoC;
Output.Xtrain = Xtrain;