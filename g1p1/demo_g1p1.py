#
# A demo program for a once-punctured torus.
#

from dev_map import *

# Create the dual (oriented) graph of a once-puctured torus.
v0 = Vertex()
e0 = Edge(v0, 0, v0, 2); e1 = Edge(v0, 1, None, -1)
# Optional: Set colors
v0.colors = ["#0000cd", "#1e90ff"]



# Set the length parameters.
l0 = 1.5; l1 = 0 # l1 corresponds to the cusp. 
# Convert length parameters to eigenvalue parameters.
L0 = -exp(l0/2); L1 = -exp(l1/2)
# Set eigenvalue parameters.
e0.set_eigenvalue(L0)
e1.set_eigenvalue(L1)

# Set the twist parameters.
t0 = l0*0.48
# Set the bending parameters.
b0 = 0.63*pi
# Set the twist-bending parameters.
e0.set_twist_by_complexFN(complex(t0, b0))

# Fix one ideal triangle on CP^1.
l = LiftedVertex(v0)
l.set_fixed_points(sqrt(3)/3*1J, -0.5-sqrt(3)/6*1J, 0.5-sqrt(3)/6*1J)

# Develope ideal triangles.
#
# CAUTION: 
#  This process consumes memory, in particular, if the representation is indiscrete.
#  Try first num_iteration = 3~8, then 50~100.
#
num_iteration = 100
l.develop_lifted_vertices(num_iteration)

# Draw ideal triangles to "demo_g1p1.eps"
drawToEPS(l, "demo_g1p1.eps", OX=300,OY=400,RATIO=100)

# Draw ideal triangles to "demo_g1p1_grafting.eps"
# This does not work if the eigenvalue parameters are not real.
drawGraftingToEPS(l, "demo_g1p1_grafting.eps", OX=300, OY=400, RATIO=100, draw_geodesics=False, draw_endpoints_of_geodesics=False)

# Create a povray file.
# This does not work if the eigenvalue parameters are not real.
writeToPovrayFile(l,"demo_g1p1.pov",RATIO=500)

# Create a povray file in the ball model.
# This does not work if the eigenvalue parameters are not real.
writeToBallModelPovrayFile(l,"demo_g1p1_ball.pov")
