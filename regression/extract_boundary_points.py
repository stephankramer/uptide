import numpy
import sys

xy_file_name = sys.argv[1]
msh_file_name = sys.argv[2]
surface_ids = [int(x) for x in sys.argv[3:]]

f = open(msh_file_name, 'r')

while f.readline().strip() != '$Nodes':
    pass

nnodes = int(f.readline())
xy = numpy.zeros((nnodes,2))
for i in range(nnodes):
    line = f.readline().strip().split(' ')
    xy[int(line[0])-1,:] = [float(line[1]), float(line[2])]

f.readline() # $EndNodes
f.readline() # $Elements
nelements =  int(f.readline())
boundary_nodes = set()
for i in range(nelements):
    line = f.readline().strip().split(' ')
    elm_type = int(line[1])
    if int(line[1]) != 1:
        continue
    sid = int(line[3])
    if sid in surface_ids:
        boundary_nodes.add(int(line[4])-1)
        boundary_nodes.add(int(line[5])-1)

f.close()
numpy.savetxt(xy_file_name, xy[list(boundary_nodes),:])
