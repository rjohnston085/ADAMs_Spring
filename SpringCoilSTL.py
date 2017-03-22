
# this is the original Matlab script that i will transfer into python


## function [out tube] = SpringCoilSTL (spring_rad, wire_dia, pitch, n_seg, NCoil, tubeDia)
## %to create an STL file of a spring element to be used by a ADAMs spring
## %generator
## 
## % clear all
## % clc
## % clf
## %
## % spring_rad=2;
## % wire_dia=1;
## % pitch=1;
## % n_seg=3;
## % NCoil=10;
## % tubeDia=5;
## 
## out=[]; %Spring patch for later graphing
## %%%%%%% Create Segment Degree Increments %%%%%%%%%%
## seg=[]; offset=2.5; ang_inc=360/n_seg;
## for i=1:n_seg
##     if isempty(seg)==1
##         seg(end+1,:)=[offset ang_inc-offset];
##     else
##         seg(end+1,:)=[seg(end,1)+ang_inc seg(end,2)+ang_inc];
##     end
##     
## end
## 
## %%%%%%%% Create Each Secgment of Spring and Save STL
## for Scoil=1:size(seg,1)
##     
##     trace=[];
##     pitch_inc=linspace(pitch*seg(Scoil,1)/360,pitch*seg(Scoil,2)/360,15);
##     count=0;
##     for i=linspace(seg(Scoil,1),seg(Scoil,2),15)
##         count=count+1;
##         trace(:,end+1)=[spring_rad*cos(i*pi()/180) spring_rad*sin(i*pi()/180) pitch_inc(count)];
##         
##     end
##     
##     %     plot3(trace(1,:),trace(2,:),trace(3,:)); hold on
##     %     plot3(trace(1,:)*0,trace(2,:)*0,trace(3,:),'or')
##     %     axis equal
##     % axis([-2.5 2.5 -2.5 2.5 0 4])
##     for  i=1:size(trace,2)
##         vect=[trace(1,i); trace(2,i); trace(3,i)];
##         vect(:,end+1)=[0 0 trace(3,i)];
##         %         plot3(vect(1,:),vect(2,:),vect(3,:),'g')
##     end
##     
##     %Wire Circle
##     skel=[]; %[(spring segement) (segement pt) (xyz)
##     for k=1:size(trace,2)-1
##         a=[trace(1,k); trace(2,k); trace(3,k)];
##         b=[trace(1,k+1); trace(2,k+1); trace(3,k+1)];
##         c=[0; 0; trace(3,k)];
##         
##         u_rad=(c-a)/magf(a,c);
##         u_ab=(b-a)/magf(a,b);
##         u_per=cross(u_rad,u_ab);
##         ang=10;
##         s=[]; degs=[0:40:120 130:10:230 240:40:320];
##         for i=1:18
##             if i>=5 && i<=17
##                 r=(wire_dia/2)*u_rad.*cos(degs(i)*pi()/180);
##                 n=(wire_dia/2)*u_per.*sin(degs(i)*pi()/180);
##                 s(end+1,:)=a+r+n;
##             else
##                 r=(wire_dia/2)*u_rad.*cos(degs(i)*pi()/180);
##                 n=(wire_dia/2)*u_per.*sin(degs(i)*pi()/180);
##                 s(end+1,:)=a+r+n;
##             end
##             
## %             plot3(s(:,1),s(:,2),s(:,3),'m')
##         end
##         
##         skel(end+1,:,:)=s;
##     end
##     
##     % Create patch surface
##     x=[]; y=[]; z=[];
##     for i=1:size(skel,1)-1
##         for j=1:size(skel,2)
##             if j<size(skel,2)
##                 x(:,end+1)=[skel(i,j,1) skel(i,j+1,1) skel(i+1,j,1)];
##                 x(:,end+1)=[skel(i+1,j+1,1) skel(i,j+1,1) skel(i+1,j,1)];
##                 y(:,end+1)=[skel(i,j,2) skel(i,j+1,2) skel(i+1,j,2)];
##                 y(:,end+1)=[skel(i+1,j+1,2) skel(i,j+1,2) skel(i+1,j,2)];
##                 z(:,end+1)=[skel(i,j,3) skel(i,j+1,3) skel(i+1,j,3)];
##                 z(:,end+1)=[skel(i+1,j+1,3) skel(i,j+1,3) skel(i+1,j,3)];
##             else
##                 x(:,end+1)=[skel(i,j,1) skel(i,1,1) skel(i+1,j,1)];
##                 x(:,end+1)=[skel(i+1,1,1) skel(i,1,1) skel(i+1,j,1)];
##                 y(:,end+1)=[skel(i,j,2) skel(i,1,2) skel(i+1,j,2)];
##                 y(:,end+1)=[skel(i+1,1,2) skel(i,1,2) skel(i+1,j,2)];
##                 z(:,end+1)=[skel(i,j,3) skel(i,1,3) skel(i+1,j,3)];
##                 z(:,end+1)=[skel(i+1,1,3) skel(i,1,3) skel(i+1,j,3)];
##             end
##         end
##         
##     end
##     %     patch(x,y,z,'b'), alpha(0.5)
##     
##     %  Create Normal for each patch
##     norm=[];
##     for i=1:size(x,2)
##         % Highlight Patch
##         %         patch(x(:,i),y(:,i),z(:,i),'r'); alpha(0.5)
##         %     find incenter
##         incenter = incenterTri (x(:,i),y(:,i),z(:,i));
##         %         plot3(incenter(1),incenter(2),incenter(3),'*k')
##         
##         if isodd(i)==0
##             v1=[x(2,i)-x(1,i) y(2,i)-y(1,i) z(2,i)-z(1,i)]; %Vector 1
##             v1=v1/magf(v1); %normal Vector 1
##             v2=[x(3,i)-x(2,i) y(3,i)-y(2,i) z(3,i)-z(2,i)]; % Vector 2
##             v2=v2/magf(v2); %normal vector 2
##         else
##             v1=[x(3,i)-x(1,i) y(3,i)-y(1,i) z(3,i)-z(1,i)]; %Vector 1
##             v1=v1/magf(v1); %normal Vector 1
##             v2=[x(2,i)-x(3,i) y(2,i)-y(3,i) z(2,i)-z(3,i)]; % Vector 2
##             v2=v2/magf(v2); %normal vector 2
##         end
##         
##         v=cross(v1,v2);
##         n=incenter+0.5*wire_dia*v;
##         %         plot3([incenter(1) n(1)],[incenter(2) n(2)],[incenter(3) n(3)],'m')
##         norm(:,end+1)=v;
##     end
##     
##     %create ends of Spring
##     for j=[1 size(skel,1)]
##         sprEnd=skel(j,:,:);
##         xe=[]; ye=[]; ze=[]; ne=[];
##         for i=2:size(sprEnd,2)-1 %patch
##             xe(:,end+1)=[skel(j,1,1) skel(j,i,1) skel(j,i+1,1)];
##             ye(:,end+1)=[skel(j,1,2) skel(j,i,2) skel(j,i+1,2)];
##             ze(:,end+1)=[skel(j,1,3) skel(j,i,3) skel(j,i+1,3)];
##         end
##         %         patch(xe,ye,ze,'b'); alpha(0.5)
##         
##         for i=1:size(xe,2) %create Normals
##             incenter=incenterTri(xe(:,i),ye(:,i),ze(:,i));
##             %             plot3(incenter(1),incenter(2),incenter(3),'*k')
##             if j==1
##                 v1=[xe(2,i)-xe(1,i) ye(2,i)-ye(1,i) ze(2,i)-ze(1,i)]; %Vector 1
##                 v1=v1/magf(v1); %normal Vector 1
##                 v2=[xe(3,i)-xe(2,i) ye(3,i)-ye(2,i) ze(3,i)-ze(2,i)]; % Vector 2
##                 v2=v2/magf(v2); %normal vector 2
##             else
##                 v1=[xe(3,i)-xe(1,i) ye(3,i)-ye(1,i) ze(3,i)-ze(1,i)]; %Vector 1
##                 v1=v1/magf(v1); %normal Vector 1
##                 v2=[xe(2,i)-xe(3,i) ye(2,i)-ye(3,i) ze(2,i)-ze(3,i)]; % Vector 2
##                 v2=v2/magf(v2); %normal vector 2
##             end
##             
##             v=cross(v1,v2);
##             n=incenter+0.5*wire_dia*v;
##             %             plot3([incenter(1) n(1)],[incenter(2) n(2)],[incenter(3) n(3)],'m')
##             ne(:,end+1)=v;
##         end
##         x=[x xe];
##         y=[y ye];
##         z=[z ze];
##         norm=[norm ne];
##     end
##     
##     if isempty(out)==1
##         out=zeros(size(seg,1),3,3,size(x,2));
##     end
##     out(Scoil,1,:,:)=x;
##     out(Scoil,2,:,:)=y;
##     out(Scoil,3,:,:)=z;
##     
##     %%%%%%%%%%%%%% Write STL %%%%%%%%%%%%%%%%%%%%
##     content={};
##     
##     content(end+1,:)=mat2cell('solid ascii');
##     for i=1:size(x,2)
##         content(end+1,:)=mat2cell(['facet normal ' num2str(norm(1,i)) ' ' num2str(norm(2,i)) ' ' num2str(norm(3,i))]);
##         content(end+1,:)=mat2cell('outer loop');
##         content(end+1,:)=mat2cell(['vertex ' num2str(x(1,i)) ' ' num2str(y(1,i)) ' ' num2str(z(1,i))]);
##         content(end+1,:)=mat2cell(['vertex ' num2str(x(2,i)) ' ' num2str(y(2,i)) ' ' num2str(z(2,i))]);
##         content(end+1,:)=mat2cell(['vertex ' num2str(x(3,i)) ' ' num2str(y(3,i)) ' ' num2str(z(3,i))]);
##         content(end+1,:)=mat2cell('endloop');
##         content(end+1,:)=mat2cell('endfacet');
##     end
##     content(end+1,:)=mat2cell('endsolid');
##     
##     % content
##     content = cellfun(@ex_func,content,'UniformOutput',0);
##     size_content = cellfun(@length,content,'UniformOutput',0);
##     str_length = max(max(cell2mat(size_content)));
##     content = cellfun(@(x) ex_func2(x,str_length),content,'uniformoutput',0);
##     content = cell2mat(content);
##     
##     
##     fid = fopen(['Spring_Coil' num2str(Scoil) '.stl'],'wt');
##     for i = 1:size(content,1)
##         fprintf(fid,'%s\n',content(i,:));
##     end
##     fclose(fid)
## end
## clc
## 
## 
## 
## %%%%%%%% Create Tube %%%%%%%%%%%%%
## x=[]; y=[]; z=[]; norm=[];
## % clf
## %%% Inner Ring %%%
## tube=tubeDia/2; %Dia to Rad
## deg_inc=2; inner=zeros(2,4,1);
## iter=0:deg_inc:360-deg_inc;
## for i=1:length(iter)
##     inner(1,:,end)=[tube*cos(iter(i)*pi()/180) tube*sin(iter(i)*pi()/180) 0 iter(i)];
##     inner(2,:,end)=[tube*cos(iter(i)*pi()/180) tube*sin(iter(i)*pi()/180) pitch*NCoil iter(i)];
##     if i~=360-deg_inc
##         inner(:,:,end+1)=zeros(2,4,1);
##     end
## end
## %%% Patch inner Ring %%%
## for i=1:size(inner,3)-1 %increments aorund tube
##     if i<size(inner,3)-1
##         x(:,end+1)=[inner(1,1,i) inner(1,1,i+1) inner(2,1,i)];
##         x(:,end+1)=[inner(2,1,i+1) inner(2,1,i) inner(1,1,i+1)];
##         y(:,end+1)=[inner(1,2,i) inner(1,2,i+1) inner(2,2,i)];
##         y(:,end+1)=[inner(2,2,i+1) inner(2,2,i) inner(1,2,i+1)];
##         z(:,end+1)=[inner(1,3,i) inner(1,3,i+1) inner(2,3,i)];
##         z(:,end+1)=[inner(2,3,i+1) inner(2,3,i) inner(1,3,i+1)];
##     else
##         x(:,end+1)=[inner(1,1,i) inner(1,1,1) inner(2,1,i)];
##         x(:,end+1)=[inner(2,1,1) inner(2,1,i) inner(1,1,1)];
##         y(:,end+1)=[inner(1,2,i) inner(1,2,1) inner(2,2,i)];
##         y(:,end+1)=[inner(2,2,1) inner(2,2,i) inner(1,2,1)];
##         z(:,end+1)=[inner(1,3,i) inner(1,3,1) inner(2,3,i)];
##         z(:,end+1)=[inner(2,3,1) inner(2,3,i) inner(1,3,1)];
##     end
##     %     patch(x,y,z,'b'); hold on
## end
## %%% Normals for Inner Ring
## for i=1:size(x,2) %create Normals
##     incenter=incenterTri(x(:,i),y(:,i),z(:,i));
##     %     plot3(incenter(1),incenter(2),incenter(3),'*k')
##     v1=[x(3,i)-x(1,i) y(3,i)-y(1,i) z(3,i)-z(1,i)]; %Vector 1
##     v1=v1/magf(v1); %normal Vector 1
##     v2=[x(2,i)-x(3,i) y(2,i)-y(3,i) z(2,i)-z(3,i)]; % Vector 2
##     v2=v2/magf(v2); %normal vector 2
##     v=cross(v1,v2);
##     ne=incenter+5*wire_dia*v;
##     %     plot3([incenter(1) ne(1)],[incenter(2) ne(2)],[incenter(3) ne(3)],'m')
##     norm(:,end+1)=v;
## end
## CarryOn=size(x,2);
## %%% Determine Outer Segment Number %%%
## outDia=tube*1.25;
## inc=3; %Number of segments
## r=0;
## while r<tube
##     r=outDia*cos(360/(inc*2)*pi()/180);
##     inc=inc+1;
## end
## %%% Outer Ring %%%
## outer=zeros(2,4,1);
## deg_inc=360/inc;
## iter=0:deg_inc:360-deg_inc;
## for i=1:length(iter)
##     
##     outer(1,:,end)=[outDia*cos(iter(i)*pi()/180) outDia*sin(iter(i)*pi()/180) 0 iter(i)];
##     outer(2,:,end)=[outDia*cos(iter(i)*pi()/180) outDia*sin(iter(i)*pi()/180) pitch*NCoil iter(i)];
##     if iter(i)~=360
##         outer(:,:,end+1)=zeros(2,4,1);
##     end
## end
## %%% Patch Outer Rings %%%
## for i=1:size(outer,3)-1 %increments aorund tube
##     if i<size(outer,3)-1
##         x(:,end+1)=[outer(1,1,i) outer(2,1,i) outer(1,1,i+1)];
##         x(:,end+1)=[outer(2,1,i+1) outer(1,1,i+1) outer(2,1,i)];
##         y(:,end+1)=[outer(1,2,i) outer(2,2,i) outer(1,2,i+1)];
##         y(:,end+1)=[outer(2,2,i+1) outer(1,2,i+1) outer(2,2,i)];
##         z(:,end+1)=[outer(1,3,i) outer(2,3,i) outer(1,3,i+1)];
##         z(:,end+1)=[outer(2,3,i+1) outer(1,3,i+1) outer(2,3,i)];
##     else
##         x(:,end+1)=[outer(1,1,i) outer(2,1,i) outer(1,1,1)];
##         x(:,end+1)=[outer(2,1,1) outer(1,1,1) outer(2,1,i)];
##         y(:,end+1)=[outer(1,2,i) outer(2,2,i) outer(1,2,1)];
##         y(:,end+1)=[outer(2,2,1) outer(1,2,1) outer(2,2,i)];
##         z(:,end+1)=[outer(1,3,i) outer(2,3,i) outer(1,3,1)];
##         z(:,end+1)=[outer(2,3,1) outer(1,3,1) outer(2,3,i)];
##     end
##     %     patch(x,y,z,'b'); hold on
## end
## %%% Normals for Outer Ring
## for i=CarryOn+1:size(x,2) %create Normals
##     incenter=incenterTri(x(:,i),y(:,i),z(:,i));
##     %     plot3(incenter(1),incenter(2),incenter(3),'*k')
##     v1=[x(3,i)-x(1,i) y(3,i)-y(1,i) z(3,i)-z(1,i)]; %Vector 1
##     v1=v1/magf(v1); %normal Vector 1
##     v2=[x(2,i)-x(3,i) y(2,i)-y(3,i) z(2,i)-z(3,i)]; % Vector 2
##     v2=v2/magf(v2); %normal vector 2
##     v=cross(v1,v2);
##     ne=incenter+0.5*wire_dia*v;
##     %     plot3([incenter(1) ne(1)],[incenter(2) ne(2)],[incenter(3) ne(3)],'m')
##     norm(:,end+1)=v;
## end
## CarryOn=size(x,2);
## %%% Patch the Lids %%%
## % iter=0:deg_inc:360-deg_inc;
## inner(:,:,end)=[];
## for k=1:2 %upper and Lower Lids
##     p_last=[]; %Last point from previous section
##     for i=1:size(outer,3)-1 % outer ring
##         if isempty(p_last)==1 %carry over the last point
##             in_p=[];
##         else
##             in_p=p_last;
##         end
##         for j=1:size(inner,3)% inner ring
##             if i==size(outer,3)-1
##                 if inner(k,4,j)>=outer(k,4,i) && inner(k,4,j)<=360
##                     in_p(:,end+1)=inner(k,:,j);
##                 end
##                 if j==size(inner,3) % for last round, include first point to close the loop
##                     in_p(:,end+1)=inner(k,:,1);
##                 end
##             else
##                 if inner(k,4,j)>=outer(k,4,i) && inner(k,4,j)<=outer(k,4,i+1)
##                     in_p(:,end+1)=inner(k,:,j);
##                 end
##             end
##         end
##         if in_p(4,1)==in_p(4,2) %deletes the carry over point if common
##             in_p(:,1)=[];
##         end
##         for m=1:size(in_p,2)-1
##             if in_p(4,m)<=median(in_p(4,:)) %first section
##                 x(:,end+1)=[outer(k,1,i) in_p(1,m) in_p(1,m+1)];
##                 y(:,end+1)=[outer(k,2,i) in_p(2,m) in_p(2,m+1)];
##                 z(:,end+1)=[outer(k,3,i) in_p(3,m) in_p(3,m+1)];
##             elseif in_p(4,m)>median(in_p(4,:)) && i==size(outer,3)-1 % To close the loop
##                 x(:,end+1)=[outer(k,1,1) in_p(1,m) in_p(1,m+1)];
##                 y(:,end+1)=[outer(k,2,1) in_p(2,m) in_p(2,m+1)];
##                 z(:,end+1)=[outer(k,3,1) in_p(3,m) in_p(3,m+1)];
##                 
##             else %second section
##                 x(:,end+1)=[outer(k,1,i+1) in_p(1,m) in_p(1,m+1)];
##                 y(:,end+1)=[outer(k,2,i+1) in_p(2,m) in_p(2,m+1)];
##                 z(:,end+1)=[outer(k,3,i+1) in_p(3,m) in_p(3,m+1)];
##             end
##             %             patch(x,y,z,'b'); hold on
##             if in_p(4,m)==median(in_p(4,:))
##                 mid=in_p(:,m+1); %collect Median point
##             end
##         end
##         p_last=in_p(:,end);
##         % Outer Triangle
##         if i==size(outer,3)-1 % To close the loop
##             x(:,end+1)=[outer(k,1,i) mid(1) outer(k,1,1)];
##             y(:,end+1)=[outer(k,2,i) mid(2) outer(k,2,1)];
##             z(:,end+1)=[outer(k,3,i) mid(3) outer(k,3,1)];
##             p_last=[];
##             %             patch(x(:,end),y(:,end),z(:,end),'r'); alpha(0.5); hold on
##             if k==1 %stamp end of 1st lid
##                 CarryOn(end+1)=size(x,2);
##             end
##         else
##             x(:,end+1)=[outer(k,1,i) mid(1) outer(k,1,i+1)];
##             y(:,end+1)=[outer(k,2,i) mid(2) outer(k,2,i+1)];
##             z(:,end+1)=[outer(k,3,i) mid(3) outer(k,3,i+1)];
##             
##         end
## %         patch(x,y,z,'b'); alpha(0.5); hold on
##     end
## end
## % clf
## % patch(x,y,z,'b'); alpha(0.5); hold on
## %%% Create Normals on Lids %%%
## for i=CarryOn(1)+1:size(x,2) %create Normals
##     incenter=incenterTri(x(:,i),y(:,i),z(:,i));
##     %     plot3(incenter(1),incenter(2),incenter(3),'*r')
##     if i<CarryOn(2)
##         v1=[x(2,i)-x(1,i) y(2,i)-y(1,i) z(2,i)-z(1,i)]; %Vector 1
##         v1=v1/magf(v1); %normal Vector 1
##         v2=[x(3,i)-x(2,i) y(3,i)-y(2,i) z(3,i)-z(2,i)]; % Vector 2
##         v2=v2/magf(v2); %normal vector 2
##         v=cross(v1,v2);
##     else
##         v1=[x(3,i)-x(1,i) y(3,i)-y(1,i) z(3,i)-z(1,i)]; %Vector 1
##         v1=v1/magf(v1); %normal Vector 1
##         v2=[x(2,i)-x(3,i) y(2,i)-y(3,i) z(2,i)-z(3,i)]; % Vector 2
##         v2=v2/magf(v2); %normal vector 2
##         v=cross(v1,v2);
##     end
##     ne=incenter+2*(v/magf(v));
##     %     plot3([incenter(1) ne(1)],[incenter(2) ne(2)],[incenter(3) ne(3)],'m')
##     norm(:,end+1)=v;
## end
## 
## tube=zeros(3,3,size(x,2));
## tube(1,:,:)=x;
## tube(2,:,:)=y;
## tube(3,:,:)=z;
## 
## 
## %%%%%%%%%%%%%% Write STL %%%%%%%%%%%%%%%%%%%%
## content={};
## content(end+1,:)=mat2cell('solid ascii');
## for i=1:size(x,2)
##     content(end+1,:)=mat2cell(['facet normal ' num2str(norm(1,i)) ' ' num2str(norm(2,i)) ' ' num2str(norm(3,i))]);
##     content(end+1,:)=mat2cell('outer loop');
##     content(end+1,:)=mat2cell(['vertex ' num2str(x(1,i)) ' ' num2str(y(1,i)) ' ' num2str(z(1,i))]);
##     content(end+1,:)=mat2cell(['vertex ' num2str(x(2,i)) ' ' num2str(y(2,i)) ' ' num2str(z(2,i))]);
##     content(end+1,:)=mat2cell(['vertex ' num2str(x(3,i)) ' ' num2str(y(3,i)) ' ' num2str(z(3,i))]);
##     content(end+1,:)=mat2cell('endloop');
##     content(end+1,:)=mat2cell('endfacet');
## end
## content(end+1,:)=mat2cell('endsolid');
## 
## % content
## content = cellfun(@ex_func,content,'UniformOutput',0);
## size_content = cellfun(@length,content,'UniformOutput',0);
## str_length = max(max(cell2mat(size_content)));
## content = cellfun(@(x) ex_func2(x,str_length),content,'uniformoutput',0);
## content = cell2mat(content);
## 
## 
## fid = fopen('Tube.stl','wt');
## for i = 1:size(content,1)
##     fprintf(fid,'%s\n',content(i,:));
## end
## fclose(fid)
## 
## end
## function out = magf (in1, in2)
## if nargin < 2
##     out=(in1(1)^2+in1(2)^2+in1(3)^2)^0.5;
## else
##     out=((in1(1)-in2(1))^2+(in1(2)-in2(2))^2+(in1(3)-in2(3))^2)^0.5;
## end
## end
## 
## function x = isodd (number)
## % isodd(number)
## % returns 1 if the number is Odd, 0 if it is even.
## a = number/2;
## whole = floor(a);
## part = a-whole;
## if part > 0;
##     x = 1;
## else
##     x = 0;
## end
## end
## 
## function out = incenterTri (x,y,z)
## p1=[x(1) y(1) z(1)];
## p2=[x(2) y(2) z(2)];
## p3=[x(3) y(3) z(3)];
## a=magf(p2,p3);
## b=magf(p1,p3);
## c=magf(p1,p2);
## p=a+b+c;
## out=[(a*p1(1)+b*p2(1)+c*p3(1))/p ...
##     (a*p1(2)+b*p2(2)+c*p3(2))/p ...
##     (a*p1(3)+b*p2(3)+c*p3(3))/p];
## end
## 
## 
## 
## 

