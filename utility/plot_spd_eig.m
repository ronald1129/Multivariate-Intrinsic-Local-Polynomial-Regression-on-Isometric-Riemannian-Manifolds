function [delta]=plot_spd_eig(S,delta)
% modified from https://www.mathworks.com/matlabcentral/fileexchange/27462-diffusion-tensor-field-dti-visualization

sz=size(S);
if length(sz)==2
    nx=1;
    ny=1;
    nz=1;
elseif length(sz)==3
    nx=sz(3);
    ny=1;
    nz=1;
elseif length(sz)==4
    nx=sz(3);
    ny=sz(4);
    nz=1;
elseif length(sz)==5
    nx=sz(3);
    ny=sz(4);
    nz=sz(5);
end
Nface=10;

if nargin<2|| isempty(delta)
    delta=1;
end


for i=1:nx
    for j=1:ny
        for k=1:nz
            [v,l]=eig(S(:,:,i,j));
            [X,Y,Z]=ellipsoid(0,0,0,l(1,1),l(2,2),l(3,3),Nface);
            sz=size(X);
            for x=1:sz(1)
                for y=1:sz(2)
                    A=[X(x,y) Y(x,y) Z(x,y)]';
                    A=v*A;
                    X(x,y)=A(1);Y(x,y)=A(2);Z(x,y)=A(3);
                end
            end
            X=X+(i-1)*delta*2;
            Y=Y+(j-1)*delta*2;
            Z=Z+(k-1)*delta*2;
            h(i)=surf(X,Y,Z);
            if i==1 && j==1 && k==1
                hold on
                if nargin<2
                    delta=1*max(l,[],"all");
                end
            end
        end
    end
end

axis equal
view([0 50]);
set(gca,'GridLineStyle','none')
% set(gca,'XTick',[])
% set(gca,'YTick',[])
% set(gca,'ZTick',[])
xlabel('X')
ylabel('Y')
zlabel('Z')
shading interp
colormap([0.8 0.8 0.8])
lighting phong
light('Position',[0 0 1],'Style','infinite','Color',[ 1.000 0.584 0.000]);
hold on




% fprintf(1,'\nIf you use plotDTI.m please cite the following work:\n');
% fprintf(1,'A. Barmpoutis, B. C. Vemuri, T. M. Shepherd, and J. R. Forder "Tensor splines for\n');
% fprintf(1,'interpolation and approximation of DT-MRI with applications to segmentation of\n');
% fprintf(1,'isolated rat hippocampi", IEEE TMI: Transactions on Medical Imaging, Vol. 26(11),\n');
% fprintf(1,'pp. 1537-1546\n');

