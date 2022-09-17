#
# A demo program for a closed surface of genus 3.
#

from dev_map import *

# Create the dual (oriented) graph of a genus 3 surface.
v0 = Vertex(); v1 = Vertex(); v2 = Vertex(); v3 = Vertex()
e0 = Edge(v0, 0, v1, 0); e1 = Edge(v0, 1, v1, 2); e2 = Edge(v0, 2, v2, 1)
e3 = Edge(v2, 0, v3, 0); e4 = Edge(v2, 2, v3, 1); e5 = Edge(v1, 1, v3, 2)

# Optional: Set colors
v0.colors = ["#ff9999", "#ffeeee" ]; v1.colors=["#99ff99", "#eeffee"];
v2.colors = ["#9999ff", "#eeeeff" ]; v3.colors=["#ff99ff", "#ffeeff"];

# Set the length parameters.
l0 = 2.0; l1 = 2.0; l2 = 2.0; l3 = 2.0; l4 = 2.05; l5 = 2.05
# Convert length parameters to eigenvalue parameters.
L0 = -exp(l0/2); L1 = -exp(l1/2); L2 = -exp(l2/2);
L3 = -exp(l3/2); L4 = -exp(l4/2); L5 = -exp(l5/2)
# Set eigenvalue parameters.
e0.set_eigenvalue(L0); e1.set_eigenvalue(L1); e2.set_eigenvalue(L2)
e3.set_eigenvalue(L3); e4.set_eigenvalue(L4); e5.set_eigenvalue(L5)

# Set the twist parameters.
t0 = 0.1*l0; t1 = 0.1*l1; t2 = 0.1*l2; t3 = 0.3*l3; t4 = 0.0*l4; t5 = 0.1*l5
# Set the bending parameters.
b0 = 0.49*pi; b1 = 0.48*pi; b2 = 0.51*pi; b3 = 0.55*pi; b4 = 0.49*pi; b5 = 0.43*pi
# Set the twist-bending parameters.
e0.set_twist_by_complexFN(complex(t0, b0))
e1.set_twist_by_complexFN(complex(t1, b1))
e2.set_twist_by_complexFN(complex(t2, b2))
e3.set_twist_by_complexFN(complex(t3, b3))
e4.set_twist_by_complexFN(complex(t4, b4))
e5.set_twist_by_complexFN(complex(t5, b5))

# Fix one ideal triangle on CP^1.
l = LiftedVertex(v0)
l.set_fixed_points(sqrt(3)/3*1J, -0.5-sqrt(3)/6*1J, 0.5-sqrt(3)/6*1J)

# Develope ideal triangles.
#
# CAUTION: 
#  This process consumes memory, in particular, if the representation is indiscrete.
#  Try first num_iteration = 3~8, then 50~100.
#
num_iteration = 50
l.develop_lifted_vertices(num_iteration)

# Draw ideal triangles to "demo_g3.eps"
drawToEPS(l, "demo_g3.eps", OX=300,OY=400,RATIO=100)

# Draw ideal triangles to "demo_g3_grafting.eps"
# This does not work if the eigenvalue parameters are not real.
drawGraftingToEPS(l, "demo_g3_grafting.eps", OX=300, OY=400, RATIO=100, draw_geodesics=False, draw_endpoints_of_geodesics=False)

# Create a povray file.
# This does not work if the eigenvalue parameters are not real.
writeToPovrayFile(l,"demo_g3.pov",RATIO=500)

# Create a povray file in the ball model.
# This does not work if the eigenvalue parameters are not real.
writeToBallModelPovrayFile(l,"demo_g3_ball.pov")
