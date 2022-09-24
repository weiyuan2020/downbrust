A=importdata("../output/PlaneSpeedoutput.txt");
N = 21;
x1=reshape(A(:,1),N,[]);
z1=reshape(-A(:,3),N,[]);
vx1=reshape(A(:,4),N,[]);
vz1=reshape(A(:,6),N,[]);
quiver(x1,z1,vx1,vz1)
