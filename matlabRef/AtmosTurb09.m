function out=AtmosTurb09(in)
%% 2018/4/19
%发现phi_r1求导计算中的系数错误
%% 2018/4/18 修改
%修改了delta_r1/delta_r偏导计算的指数，由-1/2改为-1
%% 2018/4/24 根据文献A Ring-Vortex Downburst Model for Flight Simulations
% 作部分修改。。。
%主要修改部分包括
% 1,涡量的计算
% 2、

%% 2018/4/16 修改，添加涡环旋转程序
%与原来相比，其变化为,
% % |Wx|                 |Wxp|                   |Wxl|
% % |Wy|=L(\psi)L(\theta)|Wyp|+L(-\psi)L(-\theta)|Wyl|
% % |Wz|                 |Wzp|                   |Wxl|
% % Wxp指代上涡环，wxl指代下涡环速度.
%% 流场初始化
%%in=[-919 0 -610];

x0=in(1);
y0=in(2);
z0=in(3);

H=610;
R0=915;
r0=458;
Wz0=10;%原始数据
gamma=2*R0*Wz0/(1-1/((1+(2*H/R0)^2)^1.5)); %涡量
x_vor=0;
y_vor=0;
z_vor=-H;%%涡环中心位置
% Flag=0;   %% Flag=1 代表涡核外  Flag=2 代表涡核内
%% 旋转矩阵
% 默认值设为0
psi=0;
theta=pi/12;
Lphi=[1     0           0;
    0     cos(psi)    sin(psi);
    0     -sin(psi)   cos(psi)];
Ltheta=[cos(theta)     0        -sin(theta);
    0                1       0;
    sin(theta)       0       cos(theta)];
L_phi=Lphi';
L_theta=Ltheta';
%% 涡环微分计算
% %计算phi_zm 与phi_r
% syms rs zs
% r11=sqrt((rs-R0)^2+(zs-H)^2);%距离主涡核最近点
% r12=sqrt((rs+R0)^2+(zs-H)^2);%距离主涡核最远点
% r21=sqrt((rs-R0)^2+(zs+H)^2);%距离镜像涡核最近点
% r22=sqrt((rs+R0)^2+(zs+H)^2);%距离镜像涡核最远点
%
% k1=((r12-r11)/(r11+r12));
% k2=((r22-r21)/(r21+r22));
% phi=gamma/2/pi*(0.788*k1^2*(r11+r12)/(0.25+0.75*sqrt(1-k1^2))-0.788*k2^2*(r21+r22)/(0.25+0.75*sqrt(1-k2^2)));

% phi_z=diff(phi,zs);
% phi_r=diff(phi,rs);

%% 对上部分程序做修改，改成两个涡系分开计算;在下面的计算中，l指代下涡环，u指代上涡环
% % % %计算phi_zm 与phi_r
% % syms r1 r2
% % %********************原始公式
% % % %phi=ksquare*(r11+r12)/(0.25+0.75*sqrt(1-ksquare));   %%
% % %*************************
% % % %化简出来的式子有点非常的长...考虑先对其进行初步化简。。。
% % phi=(r2-r1)^2/(0.25*(r1+r2)+1.5*(r1*r2)^（1/2）);
% % % 然后又发现，其实在计算过程中，可以考虑，先把r1,r2算出来，然后再带入微分方程的公式
% % phi_r1=diff(phi,r1);    %%delta-phi/delta-r1;
% % phi_r2=diff(phi,r2);    %%delta-phi/delta-r2;
% %
% % % % 得到的结果
% % phi_r1=(2*r1 - 2*r2)/(r1/4 + r2/4 + (3*(r1*r2)^(1/2))/2)...
% %     - ((r1 - r2)^2*((3*r2)/(4*(r1*r2)^(1/2)) + 1/4))/(r1/4 + r2/4 + (3*(r1*r2)^(1/2))/2)^2
% %
% % phi_r2=- (2*r1 - 2*r2)/(r1/4 + r2/4 + (3*(r1*r2)^(1/2))/2)...
% %     - ((r1 - r2)^2*((3*r1)/(4*(r1*r2)^(1/2)) + 1/4))/(r1/4 + r2/4 + (3*(r1*r2)^(1/2))/2)^2
% %
% % % % 接着计算r1,r2对r,z的偏导数
% % %********************原始公式
% % % % r1=sqrt((r-R0)^2+(z-H)^2);%距离主涡核最近点
% % % % r2=sqrt((r+R0)^2+(z-H)^2);%距离主涡核最远点
% % % % 经过手动化简后得到的求导公式
% % r1_r=(r1)^(-1/2)*(r-R0);       %%  delta_r1/delta_r式中r1是有具体值的r1,即将rs,zs带入后得到的
% % r1_z=(r1)^(-1/2)*(z-H);        %%  delta_r1/delta_z;
% %
% % r2_r=(r2)^(-1/2)*(r+R0);       %%  delta_r1/delta_r式中r1是有具体值的r1,即将rs,zs带入后得到的
% % r2_z=(r2)^(-1/2)*(z-H);        %%  delta_r1/delta_z;
% %
% % % % 则最后的求导公式应为
% % phi_r=phi_r1*r1_r+phi_r2*r2_r;
% % phi_z=phi_r1*r1_z+phi_r2*r2_z;

% % 在程序中，为了防止出现与程序中的自变量冲突，又不想再搞一个函数，所以，考虑将该部分的r,z做变换，改为r_cal;z_cal;

%% 此部分为坐标系转化部分
% 考虑先从体轴系转换到涡轴系，然后再转回到体轴系下。
%% 微下暴击流流场速度计算程序
% 上涡环
P=Ltheta*Lphi*[x0-x_vor y0-y_vor z0-z_vor]';
x=P(1);y=P(2);z=P(3);
r=sqrt(x*x+y*y);  %样本点与中心轴的距离   %% 取中心点的坐标为（0,0）  即为文中r_m
r_wohe=sqrt((r-R0)^2+(z)^2);    %样本点到涡核的距离

if -1<=x&&x<=1&&-1<=y&&y<=1               %涡环轴中心速度
    Wx=0;
    Wy=0;
    Wz=gamma/((2*R0)*(1+((z)/R0)^2)^(3/2));
    Wu=[Wx,Wy,Wz]';
    Flag=1;
else if r_wohe>=r0          %样本点在涡核外
        Flag=1;
        r_cal=r;
        z_cal=z;
        [phi_ru,phi_zu]=STRMFN(r_cal,z_cal,R0);
        % %         Wxu=gamma/2/pi*0.788*x/r^2*phi_zu;  %% 乘以gamma/2/pi*0.788*
        % %         Wyu=gamma/2/pi*0.788*y/r^2*phi_zu;
        % %         Wzu=-1*gamma/2/pi*0.788/r*phi_ru;
        Wu =-gamma/2/pi*0.7889*[x/r^2*phi_zu y/r^2*phi_zu -1/r*phi_ru]';   
    
    else        %涡核内速度
         Flag=2;
        if r==R0&&z==0
            Wu=[0 0 0]';
            lamda=0;
        else
            alpha=atan(y/x);                     %% 原来的仅考虑了x轴为正的情形
            if x<0
                alpha=alpha+pi;
            end
            xo=R0*cos(alpha);
            yo=R0*sin(alpha);
            zo=0;
            lamda=sqrt((r-R0)^2+(z-zo)^2)/r0;
            
            % % 计算边缘点N的速度
            xn=(x+(lamda-1)*xo)/lamda;
            yn=(y+(lamda-1)*yo)/lamda;
            zn=(z+(lamda-1)*zo)/lamda;
            rn=sqrt(xn*xn+yn*yn);
            [phi_ru,phi_zu]=STRMFN(rn,zn,R0);   
            % % 将该点转至地轴系，方便下面对其的计算
            p=L_phi*L_theta*[xn,yn,zn]'+[x_vor,y_vor,z_vor]';
            x0=p(1);y0=p(2);z0=p(3);   
      
            Wu =-gamma/2/pi*0.7889*[x/r^2*phi_zu y/r^2*phi_zu -1/r*phi_ru]';
        end
    end
end


% 下涡环
P=L_theta*L_phi*[x0-x_vor y0-y_vor z0-(-z_vor)]';      %镜像涡坐标系z轴是反的
x=P(1);y=P(2);z=P(3);
r=sqrt(x*x+y*y);  %样本点与中心轴的距离   %% 取中心点的坐标为（0,0）  即为文中r_m
r_wohe=sqrt((r-R0)^2+(z)^2);    %样本点到涡核的距离
if -1<=x&&x<=1&&-1<=y&&y<=1         %涡环中心速度
    Wx=0;
    Wy=0;
    Wz=-gamma/((2*R0)*(1+((z)/R0)^2)^(3/2));
    Wl=[Wx,Wy,Wz]';
else
    if r_wohe>=r0          %样本点在涡核外
        
        [phi_r,phi_z]=STRMFN(r,z,R0);
        Wl =(-1)*-gamma/2/pi*0.7889*[x/r^2*phi_z y/r^2*phi_z -1/r*phi_r]';
    else
        fprintf('there is warning,because the point is in the core of mirror vortex');
        Wl=[0,0,0]';
    end
end
if Flag==1
    out=L_phi*L_theta*Wu+Lphi*Ltheta*Wl;
elseif Flag==2
    out=lamda*(L_phi*L_theta*Wu+Lphi*Ltheta*Wl);
else
    fprintf('Flag=0,未赋值');
end

function [phi_r,phi_z]=STRMFN(r_cal,z_cal,R0)
            
            r1=sqrt((r_cal-R0)^2+(z_cal)^2);%距离镜像涡核最近点
            r2=sqrt((r_cal+R0)^2+(z_cal)^2);%距离镜像涡核最远点
            r1_r=(r1)^(-1)*(r_cal-R0);       %%  delta_r1/delta_r
            r1_z=(r1)^(-1)*(z_cal);        %%  delta_r1/delta_z;
            r2_r=(r2)^(-1)*(r_cal+R0);       %%  delta_r2/delta_r
            r2_z=(r2)^(-1)*(z_cal);        %%  delta_r2/delta_z;
            phi_r1 =(2*r1 - 2*r2)/(r1/4 + r2/4 + 6*(r1*r2)^(1/2)/4)...
                - (((3*(r2/r1)^(1/2))/4 + 1/4)*(r1 - r2)^2)/(r1/4 + r2/4 + 6*(r1*r2)^(1/2)/4)^2;
            
            phi_r2=- (2*r1 - 2*r2)/(r1/4 + r2/4 + 6*(r1*r2)^(1/2)/4)...
                - (((3*(r1/r2)^(1/2))/4 + 1/4)*(r1 - r2)^2)/(r1/4 + r2/4 + 6*(r1*r2)^(1/2)/4)^2;
            phi_r=phi_r1*r1_r+phi_r2*r2_r;
            phi_z=phi_r1*r1_z+phi_r2*r2_z;
            
            
            