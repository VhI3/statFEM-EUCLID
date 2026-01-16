%% Matlab mesh
%% quarterPlate_withHole_100, Created by Gmsh
%% ASCII
clear msh;
msh.nbNod = 5;
msh.POS = [
0 10 0;
10 0 0;
100 0 0;
100 100 0;
0 100 0;
];
msh.MAX = max(msh.POS);
msh.MIN = min(msh.POS);
msh.TRIANGLES =[
 1 4 5 4
 1 2 4 4
 2 3 4 4
];
