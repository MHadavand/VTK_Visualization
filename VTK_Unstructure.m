% Author: Mostafa Hadavand, hadavand@ualberta.ca
% Copyright: Mostafa Hadavand
function VTK_Unstructure(Layers,Top_Col, Base_Col, GridFl,Model,Dir, Name_Grid)

% Note this code is based on onlap
%%
% Layers='./FinalLayers.dat';
% GridFl='./GRID_Spec_R.grd';
% Model='./Real_CL_1.out';
% Dir=1;
% Name_Grid='Grid.vtk';
% Top_Col=1; Base_Col=2;
%%
% Reading the surface
gsl_Surface = readgeo( Layers );
gsl_Top=gsl_Surface(:,Top_Col); gsl_Base=gsl_Surface(:,Base_Col);

grd_set = dlmread( GridFl );

% Reading grid information
Nx = grd_set(1,1); Xmn = grd_set(1,2); Xsiz = grd_set(1,3);
Ny = grd_set(2,1); Ymn = grd_set(2,2); Ysiz = grd_set(2,3);
Nz = grd_set(3,1); Zmn = grd_set(3,2); Zsiz = grd_set(3,3);

% Set the surface matrix
[x,y] = meshgrid (	Xmn-Xsiz/2:Xsiz:(Nx-1)*Xsiz+Xmn-Xsiz/2, ...
                    Ymn-Ysiz/2:Ysiz:(Ny-1)*Ysiz+Ymn-Ysiz/2 ...
					);
% Getting the indices to load surface gslib file to matlab matrix image
% format
mat_xy = [ x(:) y(:) ];
[gsl_xy,mat2gsl] = sortrows( mat_xy, [2 1] );
[~,gsl2mat] = sortrows( gsl_xy, [1 2] );               
                
Top_srf = zeros(size(x));  Base_srf = zeros(size(x)); 
Top_srf(:) = gsl_Top(gsl2mat); Base_srf(:) = gsl_Base(gsl2mat);

[x,y,z] = meshgrid (	Xmn-Xsiz/2:Xsiz:(Nx-1)*Xsiz+Xmn-Xsiz/2, ...
						Ymn-Ysiz/2:Ysiz:(Ny-1)*Ysiz+Ymn-Ysiz/2, ...
						Zmn-Zsiz/2:Zsiz:(Nz-1)*Zsiz+Zmn-Zsiz/2  ...
					);
[sx,sy,sz] = meshgrid (	Xmn-Xsiz/2:Xsiz:(Nx)*Xsiz+Xmn-Xsiz/2, ...
						Ymn-Ysiz/2:Ysiz:(Ny)*Ysiz+Ymn-Ysiz/2, ...
						Zmn-Zsiz/2:Zsiz:(Nz)*Zsiz+Zmn-Zsiz/2  ...
					);

str_z = z; sstr_z=sz;
for i = size(z,3):-1:1
    str_z(:,:,i) = Top_srf - Dir*(Nz-i)* Zsiz - Dir*(Zsiz/2);
end

sstr_z(1:end-1,1:end-1,1:end-1) = str_z; sstr_z(:,:,end)=sstr_z(:,:,end-1)-(sstr_z(:,:,end-2)-sstr_z(:,:,end-1)); 
sstr_z(:,end,:)=sstr_z(:,end-1,:)-(sstr_z(:,end-2,:)-sstr_z(:,end-1,:));
sstr_z(end,:,:)=sstr_z(end-1,:,:)-(sstr_z(end-2,:,:)-sstr_z(end-1,:,:));

%% Reading the model
[FM, Names] = readgeo( Model );
Nvar=size(FM,2);

%%
Nb=Nx*Ny*Nz;


if (Nb~=size(FM,1)) 
    msg = 'Size not consistent!';
    error(msg)
end

NPoints=(Nx+1)*(Ny+1)*(Nz+1);
XM=zeros(1,NPoints);
YM=zeros(1,NPoints);
ZM=zeros(1,NPoints);

count=0;
for k=1:Nz+1
    for j=1:Ny+1
        for i=1:Nx+1
            
                count=count+1;
                XM(count)= Xmn - (Xsiz/2) + (i-1)*Xsiz;
                YM(count)= Ymn - (Ysiz/2) + (j-1)*Ysiz;
    %             ZM(count)= Zmn - (Zsiz/2) + (k-1)*Zsiz;
                ZM(count) = sstr_z(j,i,k);
        end
    end
end

% Connection array
CM=zeros(Nb,8);
count=0;
for k=1:Nz
    for j=1:Ny
        for i=1:Nx

                count=count+1;
                CM(count,1) = i + (j-1)*(Nx+1) + (k-1)*(Nx+1)*(Ny+1) - 1; 
                CM(count,2) = CM(count,1) + 1;
                CM(count,3) = CM(count,1) + (Nx+1);
                CM(count,4) = CM(count,3) + 1;

                CM(count,5) = CM(count,1) + (Nx+1)*(Ny+1);
                CM(count,6) = CM(count,5) + 1;
                CM(count,7) = CM(count,5) + (Nx+1);
                CM(count,8) = CM(count,7) + 1;
        end
    end
end

%%
fid = fopen(Name_Grid,'w');
fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,'Example Format\n');
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid,'POINTS %5.0f float\n',NPoints);
fprintf(fid,'%12.6f %12.6f %12.6f\n',[XM(:) YM(:) ZM(:)]');
fprintf(fid,'\nCELLS %7.0f %8.0f\n',[Nb Nb*9]);
fprintf(fid,'%3.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f\n', ...
    [ ones(Nb,1)*8 CM(1:Nb,:) ]');
fprintf(fid,'\n\nCELL_TYPES %7.0f\n',Nb);
fprintf(fid,'%3.0f %3.0f %3.0f %3.0f %3.0f %3.0f\n',ones(Nb,1)*11);
fprintf(fid,'\n\nCELL_DATA %7.0f\n',Nb);
for i=1:Nvar
    string=sprintf('SCALARS %s float%s',Names{i});
    fprintf(fid,string);
    fprintf(fid,'\nLOOKUP_TABLE default\n');
    fprintf(fid,'%5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n',FM(1:Nb,i));
end
fprintf(fid,'\n\n\n');
fclose(fid);

end
