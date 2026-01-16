%% Matlab mesh
%% quarterPlate_withHole_100, Created by Gmsh
%% ASCII
clear msh;
msh.nbNod = 7;
msh.POS = [
0 10 0;
10 0 0;
100 0 0;
100 100 0;
0 100 0;
7.071067827788018 7.071067795942932 0;
36.17851130463134 36.17851129932382 0;
];
msh.MAX = max(msh.POS);
msh.MIN = min(msh.POS);
msh.LINES =[
 2 3 3
 3 4 1
 5 1 2
];
msh.TRIANGLES =[
 4 5 7 4
 3 4 7 4
 5 1 7 4
 1 6 7 4
 2 3 7 4
 6 2 7 4
];
