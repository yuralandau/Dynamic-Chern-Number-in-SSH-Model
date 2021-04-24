% SSH模型动力学陈数数值计算
tic
clear all
% i 初始， f 末态
delta_i=0.4;     delta_f=-0.2;      % ∆

n_grid=2^10; % 格子数

% SSH模型的不动点为pi
k_vector=linspace(0,pi,n_grid); delta_k=(pi-0)/n_grid; % 从0积分到第一个不动点pi
n_grid=n_grid+1; % 区别两个变量，后期可删去
t_vector=linspace(0,pi,n_grid); delta_t=(pi-0)/(n_grid); % 从0积分到pi

[k_mesh,t_mesh]=meshgrid(k_vector,t_vector); % k 横向 t 纵向

d_i(:,:,1)=(1+delta_i)+(1-delta_i)*cos(k_mesh); % d_ix
d_i(:,:,2)=(1-delta_i)*sin(k_mesh); % d_iy
d_i(:,:,3)=zeros(size(k_mesh)); % d_iz
% d_i=[0,-delta_i*sin(k_vector),-omega*cos(k_vector)-mu];

d_f(:,:,1)=(1+delta_f)+(1-delta_f)*cos(k_mesh); % d_fx
d_f(:,:,2)=(1-delta_f)*sin(k_mesh); % d_fy
d_f(:,:,3)=zeros(size(k_mesh)); % d_fz
% d_f=[0,-delta_f*sin(k_vector),-omega*cos(k_vector)-mu];

% dhat=d_i/sqrt(d_i*d_i')*cos(2*sqrt(d_f*d_f')*t_vector)+2*d_f/sqrt(d_f*d_f')*((d_i/sqrt(d_i*d_i'))*(d_f/sqrt(d_f*d_f'))')*sin(sqrt(d_f*d_f')*t_vector)^2+cross(d_i/sqrt(d_i*d_i'),d_f/sqrt(d_f*d_f'))*sin(2*sqrt(d_f*d_f')*t_vector);
dhat_i=Normalize3D(d_i); % dhat_i对确定的k与t是一个三维矢量。
dhat_f=Normalize3D(d_f);
NormOfd_f=VectorNorm3D(d_f); % dhat_i对确定的k与t是一个标量。
NormOfd_f=1; % t'=t/|d_f|
dhat=dhat_i.*cos(2*NormOfd_f.*t_mesh) + 2*dhat_f.*dot(dhat_i,dhat_f,3).*sin(NormOfd_f.*t_mesh).*sin(NormOfd_f.*t_mesh) + cross(dhat_i,dhat_f,3).*sin(2*NormOfd_f.*t_mesh);

partial_kdhat=diff(dhat,1,2)/delta_k; % 对k的三维向量函数导数（横向求差分），维度由n_grid减小1
partial_tdhat=diff(dhat,1,1)/delta_t; % 对t的三维向量函数导数（纵向求差分），维度由n_grid+1减小1

% Bdym=cross(dhat,diff(dhat,t,1))*diff(dhat,k,1)';
Bdym=dot(cross(dhat(1:end-1,1:end-1,:),partial_tdhat(:,1:end-1,:),3),partial_kdhat(1:end-1,:,:),3);

% Cdym=int(int(Bdym,t_vector,0,pi),k_vector,0,2*pi)/(4*pi)
Cdym=sum(Bdym*delta_k*delta_t,'all')/(4*pi); % 为什么不×权重sin(t)?
Cdym
% 输出结果

% FIG. 2

hold on
for ni=1:5
    nn=1+(ni-1)*fix(n_grid/5);
    dhatPlot3D=dhat(:,nn,:);
    plot3(dhatPlot3D(:,1),dhatPlot3D(:,2),dhatPlot3D(:,3),'-*');
end
axis([-1 1 -1 1 -1 1]);
pbaspect([1 1 1]);
view(45,45);
title("动力学陈数数值计算结果"+num2str(Cdym,'%.4f'));
hold off
toc


function xhat3d=Normalize3D(x3d)

% for ni=1:1:size(x3d,1)
%     for nj=1:1:size(x3d,2)
%         x3d(ni,nj,:)=x3d(ni,nj,:)/sqrt(sum(x3d(ni,nj,:).*x3d(ni,nj,:)));
%     end
% end
xhat3d=x3d./VectorNorm3D(x3d);

end


function y2d=VectorNorm3D(x3d)

y2d=sqrt(dot(x3d,x3d,3));

end

