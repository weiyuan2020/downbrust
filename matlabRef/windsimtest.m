%% 算例验证 风场模拟
[x,z] = meshgrid(-1400:100:1400,0:30:610);
xlength=size(x,2);
zlength=size(x,1);
vx=zeros(zlength,xlength);
vz=zeros(zlength,xlength);
vy=zeros(zlength,xlength);
for i=1:zlength
    for j=1:xlength
        W=AtmosTurb09([x(i,j),100,-z(i,j)]);
        vx(i,j)=W(1);
        vz(i,j)=-W(3);
        vy(i,j)=W(2);
    end
end
quiver(x,z,vx,vz)