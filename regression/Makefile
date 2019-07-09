%.msh: %.geo
	gmsh -2 -format msh2 $<

orkney.xy: orkney.msh
	python extract_boundary_points.py $@ $< 2 3 4 5 6
