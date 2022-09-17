#
# A demo program for a closed surface of genus 2.
#

from dev_map import *

# Create the dual (oriented) graph of a genus 2 surface.
v0 = Vertex()
v1 = Vertex()
e0 = Edge(v0, 0, v1, 0) 
e1 = Edge(v0, 1, v1, 2)
e2 = Edge(v0, 2, v1, 1)
# Optional: Set colors
v0.colors = ["#ff9999", "#ffeeee"]
v1.colors=["#99ff99", "#eeffee"]

# Set the length parameters.
l0 = 2.0; l1 = 2.0; l2 = 2.0
# Convert length parameters to eigenvalue parameters.
L0 = -exp(l0/2); L1 = -exp(l1/2); L2 = -exp(l2/2);
# Set eigenvalue parameters.
e0.set_eigenvalue(L0)
e1.set_eigenvalue(L1)
e2.set_eigenvalue(L2)

# Set the twist parameters.
t0 = l0*0.96; t1 = t0; t2 = t0;
# Set the bending parameters.
b0 = 0.5*pi; b1 = b0; b2 = b0;
# Set the twist-bending parameters.
e0.set_twist_by_complexFN(complex(t0, b0))
e1.set_twist_by_complexFN(complex(t1, b1))
e2.set_twist_by_complexFN(complex(t2, b2))

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

# Draw ideal triangles to "demo_g2.eps"
drawToEPS(l, "demo_g2.eps", OX=300,OY=400,RATIO=100)

# Draw ideal triangles to "demo_g2_grafting.eps"
# This does not work if the eigenvalue parameters are not real.
drawGraftingToEPS(l, "demo_g2_grafting.eps", OX=300, OY=400, RATIO=100, draw_geodesics=False, draw_endpoints_of_geodesics=False)

# Create a povray file.
# This does not work if the eigenvalue parameters are not real.
writeToPovrayFile(l,"demo_g2.pov",RATIO=500)

# Create a povray file in the ball model.
# This does not work if the eigenvalue parameters are not real.
writeToBallModelPovrayFile(l,"demo_g2_ball.pov")
