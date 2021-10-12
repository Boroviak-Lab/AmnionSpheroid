function [Dexp1,Arc,Xtest,Output] = loadCS6Mesh_ClaraData(D3,Dexp)

fid = fopen('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/SmoothSurfaces/CS6/Proj_shots_BLEND_high_update.csv', 'rt'); 
D0_ = textscan(fid,'%f %f %f %s %s','headerLines', 1,'Delimiter',',');
fclose(fid);
D1.data = cell2mat(D0_(:,1:3));
D1.textdata =D0_{:,4};

Arc = importdata('EmDiscCS6_Arc.mat');

OBJ2=importdata('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/3D_Modellinng/Blender meshes/CS6.mat')

%Ignore bubble
Xex2 = [OBJ2.vertices(OBJ2.objects(28).data.vertices,:)]; %OBJ2.vertices(OBJ2.objects(32).data.vertices,:)
Xem2 = [OBJ2.vertices(OBJ2.objects(8).data.vertices,:)];
Xtroph2 = [OBJ2.vertices(OBJ2.objects(20).data.vertices,:)];
Xam2 = [OBJ2.vertices(OBJ2.objects(4).data.vertices,:)];
Xsys2 = [OBJ2.vertices(OBJ2.objects(16).data.vertices,:)];
Xve2 = [OBJ2.vertices(OBJ2.objects(24).data.vertices,:)];
Xpgc2 = [OBJ2.vertices(OBJ2.objects(12).data.vertices,:)];

D2 = [Xam2,zeros(size(Xam2,1),1) ; Xem2,ones(size(Xem2,1),1) ; Xsys2,2*ones(size(Xsys2,1),1) ; Xex2,3*ones(size(Xex2,1),1) ; Xve2,4*ones(size(Xve2,1),1) ; Xtroph2,5*ones(size(Xtroph2,1),1); Xpgc2,6*ones(size(Xpgc2,1),1)  ];

locations = D1.textdata(1:end,1)
for i = 1:length(locations)
    try
    idx(i,1) = find(strcmp(D3.textdata(2:end,5) , ['E15C_',strrep(locations{i},'-','_')]));   
    end    
end
CellIDs = D3.textdata(idx+1,1);
D3subs = D3.data(idx,:);
CellIDsu = D3.textdata(idx+1,6);
CellIDs2 = D3.textdata(idx+1,11);
CellIDs4 = D3.textdata(idx+1,5);

genes = Dexp.textdata(2:end,1);
headers = Dexp.textdata(1,2:end);
headers = strrep(headers,'"','');
headers = strrep(headers,'X3536STDY','3536STDY');

for i = 1:length(CellIDs)
    try
    indx(i,1) = find(strcmp(headers , CellIDs{i} ));
    catch
    end
end

%Are there any that don't exist?
CellIDs = CellIDs(find(indx~=0));
CellIDs2 = CellIDs2(find(indx~=0));
CellIDs4 = CellIDs4(find(indx~=0));
CellIDsu = CellIDsu(find(indx~=0));

XYZ = D1.data(find(indx~=0),:);


Dexp1 = Dexp.data(:,indx(find(indx~=0)));
genes = Dexp.textdata(2:end,1);
Xtrain = XYZ;
Xtrain(:,1:3) = Xtrain(:,1:3) / 400;
Xtest = D2(:,[1:3]);
Xtest(:,1:3) =Xtest(:,1:3)/400;

%Training labels
LabBin(find(strcmp(CellIDs2,'Am_CS6')==1)) = 0;
LabBin(find(strcmp(CellIDs2,'EmDisc_CS6')==1)) = 1;
LabBin(find(strcmp(CellIDs2,'VE_CS6')==1)) = 4;
LabBin(find(strcmp(CellIDs2,'SYS_CS6')==1)) = 2;
LabBin(find(strcmp(CellIDs2,'Tb_CS6')==1)) = 5;
LabBin(find(strcmp(CellIDs2,'ExMes_CS6')==1)) = 3;
LabBin(find(strcmp(CellIDs2,'PGC_CS6')==1)) = 6;

alpha = 0.05;

ind0 = find(D2(:,4)==0);
ind1 = find(D2(:,4)==1);
ind2 = find(D2(:,4)==2);
ind3 = find(D2(:,4)==3);
ind4 = find(D2(:,4)==4);
ind5 = find(D2(:,4)==5);
ind6 = find(D2(:,4)==6);

Output.Xtest = Xtest;
Output.Xtrain = Xtrain;
Output.LabBin = LabBin;
Output.genes = genes;
Output.ind0 = ind0;
Output.ind1 = ind1;
Output.ind2 = ind2;
Output.ind3 = ind3;
Output.ind4 = ind4;
Output.ind5 = ind5;
Output.ind6 = ind6;


%Output.m = Conv;
%Output.s = Convf;


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