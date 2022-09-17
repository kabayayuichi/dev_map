# 
# dev_map.py
# Ver. 0.1
#
# A python program drawing the developing map of a surface group representation to PSL(2,C) 
# parametrized by Fenchel-Nielsen coordinates.
#
# Copyright: Yuichi Kabaya

from math import sqrt, cos, sin, pi, exp, tan, acos
import cmath

# We represent a compact oriented surface with a pants decomposition
# by its dual graph. The dual graph has trivalent or monovalent
# vertices. Each trivarlent vertex corresponds to a pair of pants,
# and each monovalent vertex corresponds to a boundary component.
# A trivalent vertex has a cyclic order coming from the orientation
# of the surface.
# We give an orientations to each edge for parametrizing PSL(2,C)-representations.
class Surface():
    def __init__(self):
        self.vertices = None
        self.edges = None
    
# Vertices of the dual graph
# The corresponding pair of pants is decomposed into two ideal triangles.
# We can set the colors of these two ideal triangles, which are used when
# we draw the developing map.
class Vertex():
    def __init__(self):
        self.edges = [None, None, None] 
        # edge[0], edges[1], edges[2] are the adjacent edges.
        
        self.is_outgoing = [True, True, True] 
        # True: outgoing, False: incoming?

        self.colors = ["#5555ff", "#eeeeff"] 
        # (Only for drawing.)
        # The colors of two ideal triangles.

    def set_i_th_edge(self, e, i):
        self.edges[i] = e       
    def opposite_vertex(self, i):
        if( self.edges[i] == None ):
            return None
        if( self.is_outgoing[i] ):
            return self.edges[i].terminal_vertex
        else:
            return self.edges[i].initial_vertex                    
    def order_at_opposite_vertex(self,i):
        if( self.opposite_vertex(i) == None ):
            return -1
        else:
            if( self.is_outgoing[i] ): # This lifted edge is the initial vertex of edges[i].
                return self.edges[i].order_at_terminal 
            else:
                return self.edges[i].order_at_initial
                
# For the edges of the dual graph.
class Edge():
    def __init__(self, initial_vertex, order_at_initial, terminal_vertex, order_at_terminal):
        self.initial_vertex = initial_vertex
        self.order_at_initial = order_at_initial
        if( self.initial_vertex != None ):
            self.initial_vertex.edges[self.order_at_initial] = self
            self.initial_vertex.is_outgoing[self.order_at_initial] = True 

        self.terminal_vertex = terminal_vertex
        self.order_at_terminal = order_at_terminal
        if( self.terminal_vertex != None ):
            self.terminal_vertex.edges[self.order_at_terminal] = self
            self.terminal_vertex.is_outgoing[self.order_at_terminal] = False

        # Default values
        self.eigenvalue = -2.0 # Eigenvalue parameter (a non zero complex number)
        self.twist = 0.0 # Twist parameter (a non zero complex number)

    def __repr__(self):
        return "eigen:" + str(self.eigenvalue) + "\t" + "twist:" + str(self.twist)

    def set_eigenvalue(self, eigenvalue):
        self.eigenvalue = eigenvalue
    def set_twist(self, twist):
        self.twist = twist

    # Eigenvalue parameters and twist parameters are always defined for 
    # any PSL(2,C)-representation which is not elemenetary when ristricted 
    # to each pair-of-pants subgroup.
    # But the compex Fenchel-Nielsen parameters are not always defined.
    def set_twist_by_complexFN(self, twist):
        L1 = self.eigenvalue
        if self.initial_vertex.is_outgoing[(self.order_at_initial+1)%3] == True:
            L2 = self.initial_vertex.edges[(self.order_at_initial+1)%3].eigenvalue
        else :
            L2 = 1/self.initial_vertex.edges[(self.order_at_initial+1)%3].eigenvalue

        if self.initial_vertex.is_outgoing[(self.order_at_initial+2)%3] == True:
            L3 = self.initial_vertex.edges[(self.order_at_initial+2)%3].eigenvalue
        else :
            L3 = 1/self.initial_vertex.edges[(self.order_at_initial+2)%3].eigenvalue

        if self.terminal_vertex.is_outgoing[(self.order_at_terminal+1)%3] == True:
            L4 = self.terminal_vertex.edges[(self.order_at_terminal+1)%3].eigenvalue
        else :
            L4 = 1/self.terminal_vertex.edges[(self.order_at_terminal+1)%3].eigenvalue

        if self.terminal_vertex.is_outgoing[(self.order_at_terminal+2)%3] == True:
            L5 = self.terminal_vertex.edges[(self.order_at_terminal+2)%3].eigenvalue
        else :
            L5 = 1/self.terminal_vertex.edges[(self.order_at_terminal+2)%3].eigenvalue
            
        self.twist = cmath.exp(twist)*cmath.sqrt( ((L1*L2-L3)/(L1*L3-L2))*((1-L1*L2*L3)/(L2*L3-L1))* \
                                                  ((L1*L5-L4)/(L1*L4-L5))*((1-L1*L4*L5)/(L4*L5-L1)) )

                
class LiftedVertex(Vertex):
    def __init__(self, vertex):
        Vertex.__init__(self)
        self.underlying_vertex = vertex
        self.generation = -1 # If the generation is not set, generation = -1.
        self.fixed_points = [None, None, None] 
        self.opposite_fixed_points = [None, None, None] 
        self.adjacent_triangles = [None,None,None]
        self.is_front = True

    def __repr__(self):
        string = ""
        for i in range(3):
            if( self.underlying_vertex.is_outgoing[i] ): 
                eigen = self.underlying_vertex.edges[i].eigenvalue 
            else:
                eigen = 1/self.underlying_vertex.edges[i].eigenvalue
            string = string + "eigen:" + str(eigen) + "\t" + "f_pt:" + str(self.fixed_points[i]) + ("\n" if i != 2 else "")
        return string
    def set_fixed_points(self, f0, f1, f2): # f1, f2, f3 are fixed points
        self.fixed_points = [f0, f1, f2]
        self.compute_opposite_fixed_points()
    def set_i_th_fixed_point(self, f, i):
        self.fixed_points[i] = f
    def compute_opposite_fixed_points(self):
        e = [1.0,1.0,1.0]
        for i in range(3):
            if( self.underlying_vertex.is_outgoing[i] ): # This lifted edge is the initial vertex of edges[i].
                e[i] = self.underlying_vertex.edges[i].eigenvalue
            else:
                e[i] = 1/self.underlying_vertex.edges[i].eigenvalue
        if(   self.is_front == True ): 
            for i in range(3):
                self.opposite_fixed_points[i] = ( \
                                                  e[i]*e[i]*e[(i+1)%3]*self.fixed_points[(i+1)%3]*(self.fixed_points[i]      -self.fixed_points[(i+2)%3]) \
                                                  +         e[(i+1)%3]*self.fixed_points[(i+2)%3]*(self.fixed_points[(i+1)%3]-self.fixed_points[i]) \
                                                  +    e[i]*e[(i+2)%3]*      self.fixed_points[i]*(self.fixed_points[(i+2)%3]-self.fixed_points[(i+1)%3]) \
                                              )/( e[i]*e[i]*e[(i+1)%3]*                           (self.fixed_points[i]      -self.fixed_points[(i+2)%3]) \
                                                  +         e[(i+1)%3]*                           (self.fixed_points[(i+1)%3]-self.fixed_points[i]) \
                                                  +    e[i]*e[(i+2)%3]*                           (self.fixed_points[(i+2)%3]-self.fixed_points[(i+1)%3]) \
                                              )
        else: 
            for i in range(3):
                self.opposite_fixed_points[i] = ( \
                                                  e[i]*e[i]*e[(i+2)%3]*self.fixed_points[(i+2)%3]*(self.fixed_points[i]      -self.fixed_points[(i+1)%3]) \
                                                  +         e[(i+2)%3]*self.fixed_points[(i+1)%3]*(self.fixed_points[(i+2)%3]-self.fixed_points[i]) \
                                                  +    e[i]*e[(i+1)%3]*      self.fixed_points[i]*(self.fixed_points[(i+1)%3]-self.fixed_points[(i+2)%3]) \
                                              )/( e[i]*e[i]*e[(i+2)%3]*                           (self.fixed_points[i]      -self.fixed_points[(i+1)%3]) \
                                                  +         e[(i+2)%3]*                           (self.fixed_points[(i+2)%3]-self.fixed_points[i]) \
                                                  +    e[i]*e[(i+1)%3]*                           (self.fixed_points[(i+1)%3]-self.fixed_points[(i+2)%3]) \
                                              )

    def area(self):
        return abs( (self.fixed_points[1]-self.fixed_points[0]).real*(self.fixed_points[2]-self.fixed_points[0]).imag \
                     -(self.fixed_points[1]-self.fixed_points[0]).imag*(self.fixed_points[2]-self.fixed_points[0]).real )

    def signed_area(self):
        return (self.fixed_points[1]-self.fixed_points[0]).real*(self.fixed_points[2]-self.fixed_points[0]).imag \
            -(self.fixed_points[1]-self.fixed_points[0]).imag*(self.fixed_points[2]-self.fixed_points[0]).real

    def create_adjacent_vertex(self,i):        
        if( self.underlying_vertex.opposite_vertex(i) != None): # The opposite vertex exists.
            if( self.underlying_vertex.is_outgoing[i] ): # This lifted edge is the initial vertex of edges[i].
                self.edges[i] = Edge(self,i, LiftedVertex(self.underlying_vertex.edges[i].terminal_vertex), self.underlying_vertex.edges[i].order_at_terminal)
                e1 = self.underlying_vertex.edges[i].eigenvalue
                j = self.edges[i].order_at_terminal # The order at the opposite edge.                   
            else: # This lifted edge is the terminal vertex of edges[i].
                self.edges[i] = Edge(LiftedVertex(self.underlying_vertex.edges[i].initial_vertex), self.underlying_vertex.edges[i].order_at_initial, self,i)
                e1 = 1/self.underlying_vertex.edges[i].eigenvalue
                j = self.edges[i].order_at_initial # The order at the opposite edge.
                    
            # Compute fixed points
            if( self.underlying_vertex.is_outgoing[(i+1)%3] ):
                e2 = self.underlying_vertex.edges[(i+1)%3].eigenvalue
            else:
                e2 = 1/self.underlying_vertex.edges[(i+1)%3].eigenvalue
            if( self.underlying_vertex.is_outgoing[(i+2)%3] ):
                e3 = self.underlying_vertex.edges[(i+2)%3].eigenvalue
            else:
                e3 = 1/self.underlying_vertex.edges[(i+2)%3].eigenvalue
            if (self.underlying_vertex.opposite_vertex(i).is_outgoing[(j+1)%3]):
                e4 = self.underlying_vertex.opposite_vertex(i).edges[(j+1)%3].eigenvalue
            else:
                e4 = 1/self.underlying_vertex.opposite_vertex(i).edges[(j+1)%3].eigenvalue
            if (self.underlying_vertex.opposite_vertex(i).is_outgoing[(j+2)%3]):
                e5 = self.underlying_vertex.opposite_vertex(i).edges[(j+2)%3].eigenvalue
            else:
                e5 = 1/self.underlying_vertex.opposite_vertex(i).edges[(j+2)%3].eigenvalue
                    
            if( self.is_front == True ):
                self.opposite_vertex(i).is_front = True

                if( self.underlying_vertex.is_outgoing[i] ): # This lifted edge is the initial vertex of edges[i].
                    t1 = self.underlying_vertex.edges[i].twist
                else: # This lifted edge is the terminal vertex of edges[i].
                    t1 = ((e1*e2-e3)/(e1*e3-e2))*((e1*e5-e4)/(e1*e4-e5))/self.underlying_vertex.edges[i].twist

                x1 = self.fixed_points[i] 
                x2 = self.fixed_points[(i+1)%3]
                x3 = self.fixed_points[(i+2)%3]

                x4 = ( e1*(-(e1*e3-e2)*(e1*e4-e5)*t1+e3*(e1*e5-e4))*x1*(x2-x3)+ e1*e1*e2*(e1*e5-e4)*x2*(x3-x1)+e2*(e1*e5-e4)*x3*(x1-x2) ) / \
                     ( e1*(-(e1*e3-e2)*(e1*e4-e5)*t1+e3*(e1*e5-e4))*(x2-x3)+e1*e1*e2*(e1*e5-e4)*(x3-x1)+e2*(e1*e5-e4)*(x1-x2) )
                x5 = ( ((e1*e3-e2)*t1+e1*e3)*x1*(x2-x3)+e1*e1*e2*x2*(x3-x1)+e2*x3*(x1-x2) ) / \
                     ( ((e1*e3-e2)*t1+e1*e3)*(x2-x3)+e1*e1*e2*(x3-x1)+e2*(x1-x2) )
                self.opposite_vertex(i).set_i_th_fixed_point(x1, j)
                self.opposite_vertex(i).set_i_th_fixed_point(x4, (j+1)%3)
                self.opposite_vertex(i).set_i_th_fixed_point(x5, (j+2)%3)
                self.opposite_vertex(i).compute_opposite_fixed_points()
            else:
                self.opposite_vertex(i).is_front = False

                if( self.underlying_vertex.is_outgoing[i] ): # This lifted edge is the initial vertex of edges[i].
#                   t1 = (( (1/e1)*e5-e4)/( (1/e1)*e4-e5))*(( (1/e1)*e2-e3)/( (1/e1)*e3-e2))*self.underlying_vertex.edges[i].twist  # This is equal to the below one.
                    t1 = ((e1*e3-e2)/(e1*e2-e3))*((e1*e4-e5)/(e1*e5-e4))*self.underlying_vertex.edges[i].twist
                else: # This lifted edge is the terminal vertex of edges[i].
                    t1 = 1/self.underlying_vertex.edges[i].twist

                x1 = self.fixed_points[i] 
                x3 = self.fixed_points[(i+1)%3]
                x2 = self.fixed_points[(i+2)%3]
                
                x4 = ( e1*(-(e1*e2-e3)*(e1*e5-e4)*t1+e2*(e1*e4-e5))*x1*(x2-x3)+ e1*e1*e3*(e1*e4-e5)*x2*(x3-x1)+e3*(e1*e4-e5)*x3*(x1-x2) ) / \
                     ( e1*(-(e1*e2-e3)*(e1*e5-e4)*t1+e2*(e1*e4-e5))*(x2-x3)+e1*e1*e3*(e1*e4-e5)*(x3-x1)+e3*(e1*e4-e5)*(x1-x2) )
                x5 = ( ((e1*e2-e3)*t1+e1*e2)*x1*(x2-x3)+e1*e1*e3*x2*(x3-x1)+e3*x3*(x1-x2) ) / \
                     ( ((e1*e2-e3)*t1+e1*e2)*(x2-x3)+e1*e1*e3*(x3-x1)+e3*(x1-x2) )
                self.opposite_vertex(i).set_i_th_fixed_point(x1, j)
                self.opposite_vertex(i).set_i_th_fixed_point(x5, (j+1)%3)
                self.opposite_vertex(i).set_i_th_fixed_point(x4, (j+2)%3)
                self.opposite_vertex(i).compute_opposite_fixed_points()
        else: # The opposite vertex does not exist.
            if( self.underlying_vertex.is_outgoing[i] ): # This lifted edge is the initial vertex of edges[i].                        
                self.edges[i] = Edge(self,i, None, -1)
            else:
                self.edege[i] = Edge(None, -1, self,i)
    def create_adjacent_vertex(self,i):        
        values_of_area = []
        index_of_max_value = 0
        if( self.underlying_vertex.opposite_vertex(i) != None): # The opposite vertex exists.
            if( self.underlying_vertex.is_outgoing[i] ): # This lifted edge is the initial vertex of edges[i].
                self.edges[i] = Edge(self,i, LiftedVertex(self.underlying_vertex.edges[i].terminal_vertex), self.underlying_vertex.edges[i].order_at_terminal)
                e1 = self.underlying_vertex.edges[i].eigenvalue
                j = self.edges[i].order_at_terminal # The order at the opposite edge.                   
            else: # This lifted edge is the terminal vertex of edges[i].
                self.edges[i] = Edge(LiftedVertex(self.underlying_vertex.edges[i].initial_vertex), self.underlying_vertex.edges[i].order_at_initial, self,i)
                e1 = 1/self.underlying_vertex.edges[i].eigenvalue
                j = self.edges[i].order_at_initial # The order at the opposite edge.
                    
            # Compute fixed points
            if( self.underlying_vertex.is_outgoing[(i+1)%3] ):
                e2 = self.underlying_vertex.edges[(i+1)%3].eigenvalue
            else:
                e2 = 1/self.underlying_vertex.edges[(i+1)%3].eigenvalue
            if( self.underlying_vertex.is_outgoing[(i+2)%3] ):
                e3 = self.underlying_vertex.edges[(i+2)%3].eigenvalue
            else:
                e3 = 1/self.underlying_vertex.edges[(i+2)%3].eigenvalue
            if (self.underlying_vertex.opposite_vertex(i).is_outgoing[(j+1)%3]):
                e4 = self.underlying_vertex.opposite_vertex(i).edges[(j+1)%3].eigenvalue
            else:
                e4 = 1/self.underlying_vertex.opposite_vertex(i).edges[(j+1)%3].eigenvalue
            if (self.underlying_vertex.opposite_vertex(i).is_outgoing[(j+2)%3]):
                e5 = self.underlying_vertex.opposite_vertex(i).edges[(j+2)%3].eigenvalue
            else:
                e5 = 1/self.underlying_vertex.opposite_vertex(i).edges[(j+2)%3].eigenvalue
                    
            if( self.is_front == True ):
                self.opposite_vertex(i).is_front = True

                if( self.underlying_vertex.is_outgoing[i] ): # This lifted edge is the initial vertex of edges[i].
                    t1 = self.underlying_vertex.edges[i].twist
                else: # This lifted edge is the terminal vertex of edges[i].
                    t1 = ((e1*e2-e3)/(e1*e3-e2))*((e1*e5-e4)/(e1*e4-e5))/self.underlying_vertex.edges[i].twist

                x1 = self.fixed_points[i] 
                x2 = self.fixed_points[(i+1)%3]
                x3 = self.fixed_points[(i+2)%3]

                values_of_area = []
                for k in range(-5,5):
                    T1 = pow(e1*e1,k)*t1
                    x4 = ( e1*(-(e1*e3-e2)*(e1*e4-e5)*T1+e3*(e1*e5-e4))*x1*(x2-x3)+ e1*e1*e2*(e1*e5-e4)*x2*(x3-x1)+e2*(e1*e5-e4)*x3*(x1-x2) ) / \
                         ( e1*(-(e1*e3-e2)*(e1*e4-e5)*T1+e3*(e1*e5-e4))*(x2-x3)+e1*e1*e2*(e1*e5-e4)*(x3-x1)+e2*(e1*e5-e4)*(x1-x2) )
                    x5 = ( ((e1*e3-e2)*T1+e1*e3)*x1*(x2-x3)+e1*e1*e2*x2*(x3-x1)+e2*x3*(x1-x2) ) / \
                         ( ((e1*e3-e2)*T1+e1*e3)*(x2-x3)+e1*e1*e2*(x3-x1)+e2*(x1-x2) )
                    self.opposite_vertex(i).set_i_th_fixed_point(x1, j)
                    self.opposite_vertex(i).set_i_th_fixed_point(x4, (j+1)%3)
                    self.opposite_vertex(i).set_i_th_fixed_point(x5, (j+2)%3)
                    values_of_area.append(self.opposite_vertex(i).area())
                index_of_max_value = 0
                for k in range(0,10):
                    if (values_of_area[k] >= values_of_area[index_of_max_value] ):
                        index_of_max_value = k
                t1 = pow(e1*e1,index_of_max_value-5)*t1
                x4 = ( e1*(-(e1*e3-e2)*(e1*e4-e5)*t1+e3*(e1*e5-e4))*x1*(x2-x3)+ e1*e1*e2*(e1*e5-e4)*x2*(x3-x1)+e2*(e1*e5-e4)*x3*(x1-x2) ) / \
                     ( e1*(-(e1*e3-e2)*(e1*e4-e5)*t1+e3*(e1*e5-e4))*(x2-x3)+e1*e1*e2*(e1*e5-e4)*(x3-x1)+e2*(e1*e5-e4)*(x1-x2) )
                x5 = ( ((e1*e3-e2)*t1+e1*e3)*x1*(x2-x3)+e1*e1*e2*x2*(x3-x1)+e2*x3*(x1-x2) ) / \
                     ( ((e1*e3-e2)*t1+e1*e3)*(x2-x3)+e1*e1*e2*(x3-x1)+e2*(x1-x2) )
                self.opposite_vertex(i).set_i_th_fixed_point(x1, j)
                self.opposite_vertex(i).set_i_th_fixed_point(x4, (j+1)%3)
                self.opposite_vertex(i).set_i_th_fixed_point(x5, (j+2)%3)
                self.opposite_vertex(i).compute_opposite_fixed_points()
            else:
                self.opposite_vertex(i).is_front = False

                if( self.underlying_vertex.is_outgoing[i] ): # This lifted edge is the initial vertex of edges[i].
#                   t1 = (( (1/e1)*e5-e4)/( (1/e1)*e4-e5))*(( (1/e1)*e2-e3)/( (1/e1)*e3-e2))*self.underlying_vertex.edges[i].twist  # This is equal to the below one.
                    t1 = ((e1*e3-e2)/(e1*e2-e3))*((e1*e4-e5)/(e1*e5-e4))*self.underlying_vertex.edges[i].twist
                else: # This lifted edge is the terminal vertex of edges[i].
                    t1 = 1/self.underlying_vertex.edges[i].twist

                x1 = self.fixed_points[i] 
                x3 = self.fixed_points[(i+1)%3]
                x2 = self.fixed_points[(i+2)%3]
                
                values_of_area = []
                for k in range(-5,5):
                    T1 = pow(e1*e1,k)*t1
                    x4 = ( e1*(-(e1*e2-e3)*(e1*e5-e4)*T1+e2*(e1*e4-e5))*x1*(x2-x3)+ e1*e1*e3*(e1*e4-e5)*x2*(x3-x1)+e3*(e1*e4-e5)*x3*(x1-x2) ) / \
                         ( e1*(-(e1*e2-e3)*(e1*e5-e4)*T1+e2*(e1*e4-e5))*(x2-x3)+e1*e1*e3*(e1*e4-e5)*(x3-x1)+e3*(e1*e4-e5)*(x1-x2) )
                    x5 = ( ((e1*e2-e3)*T1+e1*e2)*x1*(x2-x3)+e1*e1*e3*x2*(x3-x1)+e3*x3*(x1-x2) ) / \
                         ( ((e1*e2-e3)*T1+e1*e2)*(x2-x3)+e1*e1*e3*(x3-x1)+e3*(x1-x2) )
                    self.opposite_vertex(i).set_i_th_fixed_point(x1, j)
                    self.opposite_vertex(i).set_i_th_fixed_point(x5, (j+1)%3)
                    self.opposite_vertex(i).set_i_th_fixed_point(x4, (j+2)%3)
                    values_of_area.append(self.opposite_vertex(i).area())
                index_of_max_value = 0
                for k in range(0,10):
                    if (values_of_area[k] >= values_of_area[index_of_max_value] ):
                        index_of_max_value = k
                t1 = pow(e1*e1,index_of_max_value-5)*t1
                x4 = ( e1*(-(e1*e2-e3)*(e1*e5-e4)*t1+e2*(e1*e4-e5))*x1*(x2-x3)+ e1*e1*e3*(e1*e4-e5)*x2*(x3-x1)+e3*(e1*e4-e5)*x3*(x1-x2) ) / \
                     ( e1*(-(e1*e2-e3)*(e1*e5-e4)*t1+e2*(e1*e4-e5))*(x2-x3)+e1*e1*e3*(e1*e4-e5)*(x3-x1)+e3*(e1*e4-e5)*(x1-x2) )
                x5 = ( ((e1*e2-e3)*t1+e1*e2)*x1*(x2-x3)+e1*e1*e3*x2*(x3-x1)+e3*x3*(x1-x2) ) / \
                     ( ((e1*e2-e3)*t1+e1*e2)*(x2-x3)+e1*e1*e3*(x3-x1)+e3*(x1-x2) )
                self.opposite_vertex(i).set_i_th_fixed_point(x1, j)
                self.opposite_vertex(i).set_i_th_fixed_point(x5, (j+1)%3)
                self.opposite_vertex(i).set_i_th_fixed_point(x4, (j+2)%3)
                self.opposite_vertex(i).compute_opposite_fixed_points()
        else: # The opposite vertex does not exist.
            if( self.underlying_vertex.is_outgoing[i] ): # This lifted edge is the initial vertex of edges[i].                        
                self.edges[i] = Edge(self,i, None, -1)
            else:
                self.edege[i] = Edge(None, -1, self,i)

    def create_adjacent_triangle(self, i):
        e0 = self.underlying_vertex.edges[0].eigenvalue if self.underlying_vertex.is_outgoing[0] else 1/self.underlying_vertex.edges[0].eigenvalue
        e1 = self.underlying_vertex.edges[1].eigenvalue if self.underlying_vertex.is_outgoing[1] else 1/self.underlying_vertex.edges[1].eigenvalue
        e2 = self.underlying_vertex.edges[2].eigenvalue if self.underlying_vertex.is_outgoing[2] else 1/self.underlying_vertex.edges[2].eigenvalue
        x0 = self.fixed_points[0]
        x1 = self.fixed_points[1]
        x2 = self.fixed_points[2]
                
        if(   self.is_front == True ):
            if(   i==0 ):
                self.adjacent_triangles[0] = LiftedVertex(self.underlying_vertex)
                self.adjacent_triangles[0].is_front = False
                self.adjacent_triangles[0].set_fixed_points( ( e1*e2*x2*(x0-x1)+e0*x1*(x2-x0) )/( e1*e2*(x0-x1)+e0*(x2-x0) ), x1, x2 )
            elif( i==1 ):
                self.adjacent_triangles[1] = LiftedVertex(self.underlying_vertex)
                self.adjacent_triangles[1].is_front = False
                self.adjacent_triangles[1].set_fixed_points( x0, ( e2*e0*x0*(x1-x2)+e1*x2*(x0-x1) )/( e2*e0*(x1-x2)+e1*(x0-x1) ), x2 )
            elif( i==2 ):
                self.adjacent_triangles[2] = LiftedVertex(self.underlying_vertex)
                self.adjacent_triangles[2].is_front = False
                self.adjacent_triangles[2].set_fixed_points( x0, x1, ( e0*e1*x1*(x2-x0)+e2*x0*(x1-x2) )/( e0*e1*(x2-x0)+e2*(x1-x2) ) )
        elif( self.is_front == False ):
            if(   i==0 ):
                self.adjacent_triangles[0] = LiftedVertex(self.underlying_vertex)
                self.adjacent_triangles[0].is_front = True
                self.adjacent_triangles[0].set_fixed_points( ( e1*e2*x1*(x0-x2)+e0*x2*(x1-x0) )/( e1*e2*(x0-x2)+e0*(x1-x0) ), x1, x2 )
            elif( i==1 ):
                self.adjacent_triangles[1] = LiftedVertex(self.underlying_vertex)
                self.adjacent_triangles[1].is_front = True
                self.adjacent_triangles[1].set_fixed_points( x0, ( e2*e0*x2*(x1-x0)+e1*x0*(x2-x1) )/( e2*e0*(x1-x0)+e1*(x2-x1) ), x2 )
            elif( i==2 ):
                self.adjacent_triangles[2] = LiftedVertex(self.underlying_vertex)
                self.adjacent_triangles[2].is_front = True
                self.adjacent_triangles[2].set_fixed_points( x0, x1, ( e0*e1*x0*(x2-x1)+e2*x1*(x0-x2) )/( e0*e1*(x2-x1)+e2*(x0-x2) ) )
                
        self.adjacent_triangles[i].adjacent_triangles[i] = self


    def develop_lifted_vertices(self, generation):
        if (generation < 0):
            return False
        self.generation = generation

        for i in range(3):        
            self.create_adjacent_triangle(i)
            self.adjacent_triangles[i].from_adjacent_triangle(i, generation-1)

        for i in range(3):
            self.create_adjacent_vertex(i)
            if( self.opposite_vertex(i) != None ):
                self.opposite_vertex(i).from_another_POP(self.order_at_opposite_vertex(i), generation-1)
            
    def from_adjacent_triangle(self, i, generation): # INTERNALLY USED
        if (generation < 0):
            return False
        self.generation = generation

        self.create_adjacent_vertex(i)
        if( (self.opposite_vertex(i) != None) and (self.opposite_vertex(i).is_small() == False) ):
            self.opposite_vertex(i).from_another_POP(self.order_at_opposite_vertex(i), generation-1)
            
        for j in range(1,3):
            self.create_adjacent_triangle((i+j)%3)
            if( self.adjacent_triangles[(i+j)%3].is_small() == False ):
                self.adjacent_triangles[(i+j)%3].from_adjacent_triangle((i+j)%3, generation-1)
            
    def from_another_POP(self, i, generation): # INTERNALLY USED
        if (generation < 0):
            return False
        self.generation = generation

        for j in range(1,3):
            self.create_adjacent_vertex((i+j)%3)
            if( (self.opposite_vertex((i+j)%3) != None) and (self.opposite_vertex((i+j)%3).is_small() == False) ):
                self.opposite_vertex((i+j)%3).from_another_POP(self.order_at_opposite_vertex((i+j)%3), generation-1)

        for j in range(3):
            self.create_adjacent_triangle(j)
            if( self.adjacent_triangles[j].is_small() == False ):
                self.adjacent_triangles[j].from_adjacent_triangle(j, generation-1)


    def is_small(self):
        if(True): # True: spherical metric  False: Euclidean metric
            if( square_of_spherical_distance(self.fixed_points[0],self.fixed_points[1]) < 0.001 or
                square_of_spherical_distance(self.fixed_points[0],self.fixed_points[2]) < 0.001 or
                square_of_spherical_distance(self.fixed_points[1],self.fixed_points[2]) < 0.001 ):
                return True
            else:
                return False
        else:
            if( abs(self.fixed_points[0]-self.fixed_points[1]) < 0.3/RATIO or
                abs(self.fixed_points[0]-self.fixed_points[2]) < 0.3/RATIO or
                abs(self.fixed_points[1]-self.fixed_points[2]) < 0.3/RATIO ):
                return True
            else:
                return False

    def center(self):
        s = ((self.fixed_points[2].real-self.fixed_points[0].real)*(self.fixed_points[2].real - self.fixed_points[1].real) \
             +(self.fixed_points[2].imag-self.fixed_points[0].imag)*(self.fixed_points[2].imag - self.fixed_points[1].imag))/ \
            ((self.fixed_points[1].imag-self.fixed_points[0].imag)*(self.fixed_points[2].real-self.fixed_points[0].real) \
             -(self.fixed_points[1].real-self.fixed_points[0].real)*(self.fixed_points[2].imag-self.fixed_points[0].imag))
        return (self.fixed_points[0]+self.fixed_points[1])/2.0 + s*1J*(self.fixed_points[0]-self.fixed_points[1])/2.0

    def radius(self):
        return abs(self.center()-self.fixed_points[0]);

def square_of_spherical_distance(z1,z2):
    v1 = to_ball_model(z1)
    v2 = to_ball_model(z2)
    return acos(inner_product(v1,v2))

def to_ball_model(z):
    z = complex(z)
    r_sq = (z.real)**2+(z.imag)**2
    x1 = (2*z.real)/(1+r_sq)
    x2 = (2*z.imag)/(1+r_sq)
    x3 = (r_sq-1)/(r_sq+1)
    return [x1,x2,x3]

def inner_product(v1,v2):
    return sum([v1[i]*v2[i] for i in range(len(v1))])

def scalar_mult(c, v):
    return [c*x for x in v]

def normalized_vector(v):
    r = sqrt( sum([x**2 for x in v]) )
    if r==0.0:
        return [0 for i in range(len(v))]
    return [x/r for x in v]

def signed_area(z0,z1,z2):
    return ((z1-z0).real)*((z2-z0).imag) - ((z1-z0).imag)*((z2-z0).real)

#
# For EPS files
#
def drawToEPS(lifted_vertex, filename, width = 800, height=800, OX=400, OY=400, RATIO=200, draw_geodesics=False, draw_endpoints_of_geodesics=True):
    f = open(filename, "w")
    f.write("%!PS-Adobe-3.0 EPSF-3.0\n")
    f.write("%%%%BoundingBox: 0 0 %d %d\n" % (width, height))

    drawRelatedTrianglesToEPS(lifted_vertex, f, OX, OY, RATIO, draw_geodesics, draw_endpoints_of_geodesics)
    
    f.write("showpage\n")
    f.close()

def drawRelatedTrianglesToEPS(lifted_vertex, f, OX, OY, RATIO, draw_geodesics, draw_endpoints_of_geodesics):
    if(False): # Only draw positively (negatively?) oriented triangles
        if(lifted_vertex.signed_area() > 0 and lifted_vertex.is_front == True):
            drawTriangleToEPS(lifted_vertex, f, OX, OY, RATIO)
        if(lifted_vertex.signed_area() < 0 and lifted_vertex.is_front == False):
            drawTriangleToEPS(lifted_vertex, f, OX, OY, RATIO)
                
        for i in range(3):
            if( (lifted_vertex.adjacent_triangles[i] != None) and (lifted_vertex.adjacent_triangles[i].generation < lifted_vertex.generation) ): 
                drawRelatedTrianglesToEPS(lifted_vertex.adjacent_triangles[i],f, OX, OY, RATIO, draw_geodesics, draw_endpoints_of_geodesics)
        for i in range(3):
            if( (lifted_vertex.opposite_vertex(i) != None) and (lifted_vertex.opposite_vertex(i).generation < lifted_vertex.generation ) ):
                drawLineToEPS(lifted_vertex, i,f, OX, OY, RATIO, draw_geodesics, draw_endpoints_of_geodesics)
                drawRelatedTrianglesToEPS(lifted_vertex.opposite_vertex(i),f, OX, OY, RATIO, draw_geodesics, draw_endpoints_of_geodesics)
    else:
        drawTriangleToEPS(lifted_vertex, f, OX, OY, RATIO)
                
        for i in range(3):
            if( (lifted_vertex.adjacent_triangles[i] != None) and (lifted_vertex.adjacent_triangles[i].generation < lifted_vertex.generation) ): 
                drawRelatedTrianglesToEPS(lifted_vertex.adjacent_triangles[i],f, OX, OY, RATIO, draw_geodesics, draw_endpoints_of_geodesics)
        for i in range(3):
            if( (lifted_vertex.opposite_vertex(i) != None) and (lifted_vertex.opposite_vertex(i).generation < lifted_vertex.generation ) ):
                drawLineToEPS(lifted_vertex, i,f, OX, OY, RATIO, draw_geodesics, draw_endpoints_of_geodesics)
                drawRelatedTrianglesToEPS(lifted_vertex.opposite_vertex(i),f, OX, OY, RATIO, draw_geodesics, draw_endpoints_of_geodesics)
def drawTriangleToEPS(lifted_vertex, f, OX, OY, RATIO):
    if(True):
        f.write("0 setlinewidth\n")
        f.write("newpath\n")
        f.write("%.1f %.1f moveto\n" % (OX+lifted_vertex.fixed_points[0].real*RATIO, OY+lifted_vertex.fixed_points[0].imag*RATIO))
        f.write("%.1f %.1f lineto\n" % (OX+lifted_vertex.fixed_points[1].real*RATIO, OY+lifted_vertex.fixed_points[1].imag*RATIO))
        f.write("%.1f %.1f lineto\n" % (OX+lifted_vertex.fixed_points[2].real*RATIO, OY+lifted_vertex.fixed_points[2].imag*RATIO))
        f.write("closepath\n")
        if lifted_vertex.is_front:
            rgb = lifted_vertex.underlying_vertex.colors[0]
        else:
            rgb = lifted_vertex.underlying_vertex.colors[1]
        f.write("%f %f %f setrgbcolor\n" % (r_of_rgb(rgb), g_of_rgb(rgb), b_of_rgb(rgb)))
        f.write("fill\n\n")
def drawLineToEPS(lifted_vertex, i, f, OX, OY, RATIO, draw_geodesics, draw_endpoints_of_geodesics): 
    if(draw_geodesics): # If True, draw geodesics.
        f.write("0.5 setlinewidth\n")
        f.write("0 0 0 setrgbcolor\n")
        f.write("newpath\n")
        f.write("%f %f moveto\n" % (OX+lifted_vertex.fixed_points[i].real*RATIO, OY+lifted_vertex.fixed_points[i].imag*RATIO))
        f.write("%f %f lineto\n" % (OX+lifted_vertex.opposite_fixed_points[i].real*RATIO, OY+lifted_vertex.opposite_fixed_points[i].imag*RATIO))
        f.write("stroke\n\n")
    if(draw_endpoints_of_geodesics): # If True, draw fixed points.
        f.write("0 setlinewidth\n")
        f.write("newpath\n")
        f.write("%.1f %.1f moveto\n" % (OX+lifted_vertex.fixed_points[i].real*RATIO-0.2, OY+lifted_vertex.fixed_points[i].imag*RATIO-0.2))
        f.write("%.1f %.1f lineto\n" % (OX+lifted_vertex.fixed_points[i].real*RATIO-0.2, OY+lifted_vertex.fixed_points[i].imag*RATIO+0.2))
        f.write("%.1f %.1f lineto\n" % (OX+lifted_vertex.fixed_points[i].real*RATIO+0.2, OY+lifted_vertex.fixed_points[i].imag*RATIO+0.2))
        f.write("%.1f %.1f lineto\n" % (OX+lifted_vertex.fixed_points[i].real*RATIO+0.2, OY+lifted_vertex.fixed_points[i].imag*RATIO-0.2))
        f.write("closepath\n")
        f.write("1 0 0 setrgbcolor\n")
        f.write("fill\n\n")
        f.write("newpath\n")
        f.write("%.1f %.1f moveto\n" % (OX+lifted_vertex.opposite_fixed_points[i].real*RATIO-0.2, OY+lifted_vertex.opposite_fixed_points[i].imag*RATIO-0.2))
        f.write("%.1f %.1f lineto\n" % (OX+lifted_vertex.opposite_fixed_points[i].real*RATIO-0.2, OY+lifted_vertex.opposite_fixed_points[i].imag*RATIO+0.2))
        f.write("%.1f %.1f lineto\n" % (OX+lifted_vertex.opposite_fixed_points[i].real*RATIO+0.2, OY+lifted_vertex.opposite_fixed_points[i].imag*RATIO+0.2))
        f.write("%.1f %.1f lineto\n" % (OX+lifted_vertex.opposite_fixed_points[i].real*RATIO+0.2, OY+lifted_vertex.opposite_fixed_points[i].imag*RATIO-0.2))
        f.write("closepath\n")
        f.write("1 0 0 setrgbcolor\n")
        f.write("fill\n\n")
def r_of_rgb(rgb):
    r = int(rgb[1],16)*16+int(rgb[2],16)
    return r/256.0
def g_of_rgb(rgb):
    g = int(rgb[3],16)*16+int(rgb[4],16)
    return g/256.0
def b_of_rgb(rgb):
    b = int(rgb[5],16)*16+int(rgb[6],16)
    return b/256.0


# 
# For grafting
#
def drawGraftingToEPS(lifted_vertex, filename, width=800, height=800, OX=400, OY=400, RATIO=200, draw_geodesics=False, draw_endpoints_of_geodesics=False):
    f = open(filename, "w")
    f.write("%!PS-Adobe-3.0 EPSF-3.0\n")
    f.write("%%%%BoundingBox: 0 0 %d %d\n" % (width, height))

    drawRelatedIdealTrianglesToEPS(lifted_vertex, f, OX, OY, RATIO, draw_geodesics, draw_endpoints_of_geodesics)
    
    f.write("showpage\n")
    f.close()

def drawRelatedIdealTrianglesToEPS(lifted_vertex, f, OX, OY, RATIO, draw_geodesics, draw_endpoints_of_geodesics):
    drawLuneToEPS(lifted_vertex, f, OX, OY, RATIO)
    drawIdealTriangleToEPS(lifted_vertex, f, OX, OY, RATIO)
    if draw_endpoints_of_geodesics:
        for i in range(3):
            f.write("0 setlinewidth\n")
            f.write("newpath\n")
            f.write("%.1f %.1f moveto\n" % (OX+lifted_vertex.fixed_points[i].real*RATIO-0.2, OY+lifted_vertex.fixed_points[i].imag*RATIO-0.2))
            f.write("%.1f %.1f lineto\n" % (OX+lifted_vertex.fixed_points[i].real*RATIO-0.2, OY+lifted_vertex.fixed_points[i].imag*RATIO+0.2))
            f.write("%.1f %.1f lineto\n" % (OX+lifted_vertex.fixed_points[i].real*RATIO+0.2, OY+lifted_vertex.fixed_points[i].imag*RATIO+0.2))
            f.write("%.1f %.1f lineto\n" % (OX+lifted_vertex.fixed_points[i].real*RATIO+0.2, OY+lifted_vertex.fixed_points[i].imag*RATIO-0.2))
            f.write("closepath\n")
            f.write("1 0 0 setrgbcolor\n")
            f.write("fill\n\n")
            f.write("newpath\n")
            f.write("%.1f %.1f moveto\n" % (OX+lifted_vertex.opposite_fixed_points[i].real*RATIO-0.2, OY+lifted_vertex.opposite_fixed_points[i].imag*RATIO-0.2))
            f.write("%.1f %.1f lineto\n" % (OX+lifted_vertex.opposite_fixed_points[i].real*RATIO-0.2, OY+lifted_vertex.opposite_fixed_points[i].imag*RATIO+0.2))
            f.write("%.1f %.1f lineto\n" % (OX+lifted_vertex.opposite_fixed_points[i].real*RATIO+0.2, OY+lifted_vertex.opposite_fixed_points[i].imag*RATIO+0.2))
            f.write("%.1f %.1f lineto\n" % (OX+lifted_vertex.opposite_fixed_points[i].real*RATIO+0.2, OY+lifted_vertex.opposite_fixed_points[i].imag*RATIO-0.2))
            f.write("closepath\n")
            f.write("1 0 0 setrgbcolor\n")
            f.write("fill\n\n")
    if(draw_geodesics): # If True, draw geodesics.
        for i in range(3):
            f.write("0.5 setlinewidth\n")
            f.write("0 0 0 setrgbcolor\n")
            f.write("newpath\n")
            f.write("%f %f moveto\n" % (OX+lifted_vertex.fixed_points[i].real*RATIO, OY+lifted_vertex.fixed_points[i].imag*RATIO))
            f.write("%f %f lineto\n" % (OX+lifted_vertex.opposite_fixed_points[i].real*RATIO, OY+lifted_vertex.opposite_fixed_points[i].imag*RATIO))
            f.write("stroke\n\n")
        

    for i in range(3):
        if( (lifted_vertex.adjacent_triangles[i] != None) and (lifted_vertex.adjacent_triangles[i].generation < lifted_vertex.generation) ): 
            drawRelatedIdealTrianglesToEPS(lifted_vertex.adjacent_triangles[i],f, OX, OY, RATIO, draw_geodesics, draw_endpoints_of_geodesics)
    for i in range(3):
        if( (lifted_vertex.opposite_vertex(i) != None) and (lifted_vertex.opposite_vertex(i).generation < lifted_vertex.generation ) ):
            drawRelatedIdealTrianglesToEPS(lifted_vertex.opposite_vertex(i),f, OX, OY, RATIO, draw_geodesics, draw_endpoints_of_geodesics)

def drawIdealTriangleToEPS(lifted_vertex, f, OX, OY, RATIO):
    c = lifted_vertex.center()
    z0 = lifted_vertex.fixed_points[0]
    z1 = lifted_vertex.fixed_points[1]
    z2 = lifted_vertex.fixed_points[2]

    th0 = cmath.phase( (z2 - c)/(z1 - c) )
    th1 = cmath.phase( (z0 - c)/(z2 - c) )
    th2 = cmath.phase( (z1 - c)/(z0 - c) )

    l0 = abs(z0-c)*tan(th0/2.0)
    l1 = abs(z1-c)*tan(th1/2.0)
    l2 = abs(z2-c)*tan(th2/2.0)
    
    det = (1.0/(1+l0/l1))*(1.0/(1+l1/l2))*(1.0/(1+l2/l0)) + (1.0/(1+l1/l0))*(1.0/(1+l2/l1))*(1.0/(1+l0/l2))
    
    w0 = ( -(1.0/(1+l1/l0))*(1.0/(1+l2/l0))*z0 + (1.0/(1+l1/l0))*(1.0/(1+l2/l1))*z1 + (1.0/(1+l1/l2))*(1.0/(1+l2/l0))*z2)/det
    w1 = (  (1.0/(1+l0/l1))*(1.0/(1+l2/l0))*z0 - (1.0/(1+l2/l1))*(1.0/(1+l0/l1))*z1 + (1.0/(1+l0/l2))*(1.0/(1+l2/l1))*z2)/det
    w2 = (  (1.0/(1+l0/l2))*(1.0/(1+l1/l0))*z0 + (1.0/(1+l1/l2))*(1.0/(1+l0/l1))*z1 - (1.0/(1+l1/l2))*(1.0/(1+l0/l2))*z2)/det

    f.write("0 setlinewidth\n")
    f.write("newpath\n")

    if l0 > 0:
        f.write("%.1f %.1f %.1f %.1f %.1f arcn\n" % (OX+w0.real*RATIO, OY+w0.imag*RATIO, abs(l0)*RATIO, cmath.phase(z1-w0)*(180.0/pi), cmath.phase(z2-w0)*(180.0/pi)))
    else:
        f.write("%.1f %.1f %.1f %.1f %.1f arc\n" % (OX+w0.real*RATIO, OY+w0.imag*RATIO, abs(l0)*RATIO, cmath.phase(z1-w0)*(180.0/pi), cmath.phase(z2-w0)*(180.0/pi)))        
    if l1 > 0:
        f.write("%.1f %.1f %.1f %.1f %.1f arcn\n" % (OX+w1.real*RATIO, OY+w1.imag*RATIO, abs(l1)*RATIO, cmath.phase(z2-w1)*(180.0/pi), cmath.phase(z0-w1)*(180.0/pi)))
    else:
        f.write("%.1f %.1f %.1f %.1f %.1f arc\n" % (OX+w1.real*RATIO, OY+w1.imag*RATIO, abs(l1)*RATIO, cmath.phase(z2-w1)*(180.0/pi), cmath.phase(z0-w1)*(180.0/pi)))
    if l2 > 0:
        f.write("%.1f %.1f %.1f %.1f %.1f arcn\n" % (OX+w2.real*RATIO, OY+w2.imag*RATIO, abs(l2)*RATIO, cmath.phase(z0-w2)*(180.0/pi), cmath.phase(z1-w2)*(180.0/pi)))
    else:
        f.write("%.1f %.1f %.1f %.1f %.1f arc\n" % (OX+w2.real*RATIO, OY+w2.imag*RATIO, abs(l2)*RATIO, cmath.phase(z0-w2)*(180.0/pi), cmath.phase(z1-w2)*(180.0/pi)))

    if lifted_vertex.is_front:
        rgb = lifted_vertex.underlying_vertex.colors[0]
    else:
        rgb = lifted_vertex.underlying_vertex.colors[1]
    f.write("%f %f %f setrgbcolor\n" % (r_of_rgb(rgb), g_of_rgb(rgb), b_of_rgb(rgb)))
    f.write("fill\n\n")


def drawLuneToEPS(lifted_vertex, f, OX, OY, RATIO):
    c = lifted_vertex.center()
    for i in range(3):
        if( (lifted_vertex.opposite_vertex(i) != None) and (lifted_vertex.opposite_vertex(i).generation < lifted_vertex.generation ) ):
            z = lifted_vertex.fixed_points[i]
            w = lifted_vertex.opposite_fixed_points[i]

            oc = lifted_vertex.opposite_vertex(i).center()

            th = cmath.phase( (w - c)/(z - c) )
            oth = cmath.phase( (z - oc)/(w - oc) )
            
            if abs(th) > 0.0001 and abs(oth) > 0.0001: 
                l = abs(z-c)/cos(th/2.0)
                ol = abs(z-oc)/cos(oth/2.0)

                p = c + (l/abs(z+w-2*c))*((z-c) + (w-c))
                op = oc + (ol/abs(z+w-2*oc))*((z-oc) + (w-oc))

                f.write("1 setlinewidth\n")
                f.write("newpath\n")
                if signed_area(p,z,w) < 0:
                    f.write("%.1f %.1f %.1f %.1f %.1f arcn\n" % (OX+p.real*RATIO, OY+p.imag*RATIO, abs(z-p)*RATIO, cmath.phase(z-p)*(180.0/pi), cmath.phase(w-p)*(180.0/pi)))
                else:
                    f.write("%.1f %.1f %.1f %.1f %.1f arc\n" % (OX+p.real*RATIO, OY+p.imag*RATIO, abs(z-p)*RATIO, cmath.phase(z-p)*(180.0/pi), cmath.phase(w-p)*(180.0/pi)))
                if signed_area(op,z,w) > 0:
                    f.write("%.1f %.1f %.1f %.1f %.1f arcn\n" % (OX+op.real*RATIO, OY+op.imag*RATIO, abs(z-op)*RATIO, cmath.phase(w-op)*(180.0/pi), cmath.phase(z-op)*(180.0/pi)))
                else:
                    f.write("%.1f %.1f %.1f %.1f %.1f arc\n" % (OX+op.real*RATIO, OY+op.imag*RATIO, abs(z-op)*RATIO, cmath.phase(w-op)*(180.0/pi), cmath.phase(z-op)*(180.0/pi)))

                f.write("%f %f %f setrgbcolor\n" % (1, 1, 0.6))
                f.write("fill\n\n")

#
# For POV-Ray file
#
def writeToPovrayFile(lifted_vertex, filename, OX=400, OY=400, RATIO=200):
    color_list = []
    
    f = open(filename, "w")

    f.write("#include \"colors.inc\"\n\n")
    f.write("camera {\n")
    f.write(" location <100, 2500, 0>\n")
    f.write(" look_at <0, 0, 0>\n")
    f.write(" angle 80\n")
    f.write(" translate <clock,0,0>\n")
    f.write("}\n\n")

    f.write("sky_sphere {\n")
    f.write(" pigment {\n")
    f.write("  gradient y\n")
    f.write("  color_map {\n")
    f.write("   [0.0 rgb <0.6,0.7,1.0>]\n")
    f.write("   [0.7 rgb <0.0,0.1,0.8>]\n")
    f.write("  }\n")
    f.write(" }\n")
    f.write("}\n\n")

    f.write("light_source {\n")
    f.write(" <0, 300, -1000>\n")
    f.write(" color red 1 green 1 blue 1\n")
    f.write("}\n\n")

    f.write("light_source {\n")
    f.write(" <500, 300, 1000>\n")
    f.write(" color red 1 green 1 blue 1\n")
    f.write("}\n\n")

    f.write("light_source {\n")
    f.write(" <1000, 300, 0>\n")
    f.write(" color red 1 green 1 blue 1\n")
    f.write("}\n\n")

    f.write("plane {\n")
    f.write(" y, 0\n")
    f.write(" texture {\n")
    f.write("  pigment {\n")
    f.write("   color rgb <0, 0.6, 0.6>\n")
    f.write("  }\n")
    f.write(" }\n")
    f.write("}\n\n")

    rgb = lifted_vertex.underlying_vertex.colors[0]
    color_list.append(rgb[1:7])
    f.write("#declare c_%s=\n"%rgb[1:7])    
    f.write("texture {\n")
    f.write(" pigment {\n")
    f.write("  color rgb <%f,%f,%f>\n" % (r_of_rgb(rgb), g_of_rgb(rgb), b_of_rgb(rgb)))
    f.write(" }\n")
    f.write(" finish{\n")
    f.write("  diffuse 0.3\n")
    f.write("  ambient 0.0\n")
    f.write("  specular 0.6\n")
    f.write("  reflection {\n")
    f.write("   0.2\n") 
    f.write("   metallic\n")
    f.write("  }\n")
    f.write("  conserve_energy\n")
    f.write(" }\n")
    f.write("}\n\n")

    f.write("#declare i_%s=\n"%rgb[1:7])    
    f.write("interior {}\n\n")

    f.write("sphere {\n")
    f.write(" <%f, 0, %f>\n" % (lifted_vertex.center().real*RATIO, lifted_vertex.center().imag*RATIO))
    f.write(" %f\n" % (lifted_vertex.radius()*RATIO))
    f.write(" texture {c_%s}\n"%rgb[1:7])
    f.write(" interior {i_%s}\n"%rgb[1:7])
    f.write("}\n")

    drawRelatedCirclesToPovrayFile(lifted_vertex, f, OX, OY, RATIO, color_list)

    f.close()

def drawRelatedCirclesToPovrayFile(lifted_vertex, f, OX, OY, RATIO, color_list):
    for i in range(3):
        if( (lifted_vertex.adjacent_triangles[i] != None) and (lifted_vertex.adjacent_triangles[i].generation < lifted_vertex.generation) ): 
            drawRelatedCirclesToPovrayFile(lifted_vertex.adjacent_triangles[i],f, OX, OY, RATIO, color_list)
    for i in range(3):
        if( (lifted_vertex.opposite_vertex(i) != None) and (lifted_vertex.opposite_vertex(i).generation < lifted_vertex.generation ) ):
            drawCircleToPovrayFile(lifted_vertex.opposite_vertex(i),f, OX, OY, RATIO, color_list)
            drawRelatedCirclesToPovrayFile(lifted_vertex.opposite_vertex(i),f, OX, OY, RATIO, color_list)

def drawCircleToPovrayFile(lifted_vertex, f, OX, OY, RATIO, color_list):
    rgb = lifted_vertex.underlying_vertex.colors[0]
    if color_list.count(rgb[1:7])==0:
        f.write("\n#declare c_%s=\n"%rgb[1:7])    
        f.write("texture {\n")
        f.write(" pigment {\n")
        f.write("  color rgb <%f,%f,%f>\n" % (r_of_rgb(rgb), g_of_rgb(rgb), b_of_rgb(rgb)))
        f.write(" }\n")
        f.write(" finish{\n")
        f.write("  diffuse 0.3\n")
        f.write("  ambient 0.0\n")
        f.write("  specular 0.6\n")
        f.write("  reflection {\n")
        f.write("   0.2\n") 
        f.write("   metallic\n")
        f.write("  }\n")
        f.write("  conserve_energy\n")
        f.write(" }\n")
        f.write("}\n\n")
        color_list.append(rgb[1:7])

        f.write("#declare i_%s=\n"%rgb[1:7])    
        f.write("interior {}\n\n")

    f.write("sphere {\n")
    f.write("<%f, 0, %f>\n" % (lifted_vertex.center().real*RATIO, lifted_vertex.center().imag*RATIO))
    f.write("%f\n" % (lifted_vertex.radius()*RATIO))
    f.write("texture {c_%s}\n"%rgb[1:7])
    f.write("interior {i_%s}\n"%rgb[1:7])
    f.write("}\n")


#
# For POV-Ray file  -- Ball Model --
#
def writeToBallModelPovrayFile(lifted_vertex, filename):
    color_list = []
    
    f = open(filename, "w")

    f.write("#include \"colors.inc\"\n\n")

    f.write("camera {\n")
    f.write(" location  <0.5, -1.5, -1.5>\n")
    f.write(" look_at   <0, 0, 0>\n")
    f.write(" angle     80\n")
    f.write(" translate <clock,0,0>\n")
    f.write("}\n\n")

    f.write("sky_sphere {\n")
    f.write(" pigment {\n")
    f.write("  gradient y\n")
    f.write("  color_map {\n")
    f.write("   [0.0 rgb <0.6,0.7,1.0>]\n")
    f.write("   [0.7 rgb <0.0,0.1,0.8>]\n")
    f.write("  }\n")
    f.write(" }\n")
    f.write("}\n\n")

    f.write("light_source {\n")
    f.write(" <0, 300, -1000>\n")          
    f.write(" color red 1 green 1 blue 1\n")
    f.write("}\n\n")

    f.write("light_source {\n")
    f.write(" <500, 300, 1000>\n")           
    f.write(" color red 1 green 1 blue 1\n")
    f.write("}\n\n")

    f.write("light_source {\n")
    f.write(" <1000, 300, 0>\n")         
    f.write(" color red 1 green 1 blue 1\n")
    f.write("}\n\n")

    f.write("difference{sphere {<0.0,0.0,0.0>, 1.00 pigment { rgbf <0.8,0.8,1,0.9> }}")	

    rgb = lifted_vertex.underlying_vertex.colors[0]
    color_list.append(rgb[1:7])

    if abs(lifted_vertex.center()) != 0:
        max_pt = lifted_vertex.center() + \
                 abs(lifted_vertex.radius())/abs(lifted_vertex.center())*lifted_vertex.center()
        min_pt = lifted_vertex.center() - \
                 abs(lifted_vertex.radius())/abs(lifted_vertex.center())*lifted_vertex.center()
    else:
        max_pt = lifted_vertex.radius()
        min_pt = -lifted_vertex.radius()
        
    max_on_sp = to_ball_model(max_pt)
    min_on_sp = to_ball_model(min_pt)
    theta = acos(inner_product(max_on_sp, min_on_sp))

    radius = tan(theta/2)
    vect = [max_on_sp[0]+min_on_sp[0], max_on_sp[1]+min_on_sp[1], max_on_sp[2]+min_on_sp[2]]
    point = scalar_mult(1/(cos(theta/2)),normalized_vector(vect))

    f.write("sphere {\n")
    f.write("<%f, %f, %f>\n" % (point[0], point[1], point[2]))
    f.write("%f\n" % (radius))
    f.write("texture{ pigment{ color rgb <%f,%f,%f>} finish { phong 1 }}\n" % (r_of_rgb(rgb), g_of_rgb(rgb), b_of_rgb(rgb)))
    f.write("}\n")

    drawRelatedCirclesToBallModelPovrayFile(lifted_vertex, f, color_list)

    f.write("}")            
    f.close()

def drawRelatedCirclesToBallModelPovrayFile(lifted_vertex, f, color_list):
    for i in range(3):
        if( (lifted_vertex.adjacent_triangles[i] != None) and (lifted_vertex.adjacent_triangles[i].generation < lifted_vertex.generation) ): 
            drawRelatedCirclesToBallModelPovrayFile(lifted_vertex.adjacent_triangles[i],f, color_list)
    for i in range(3):
        if( (lifted_vertex.opposite_vertex(i) != None) and (lifted_vertex.opposite_vertex(i).generation < lifted_vertex.generation ) ):
            drawCircleToBallModelPovrayFile(lifted_vertex.opposite_vertex(i),f, color_list)
            drawRelatedCirclesToBallModelPovrayFile(lifted_vertex.opposite_vertex(i),f, color_list)

def drawCircleToBallModelPovrayFile(lifted_vertex, f, color_list):
    rgb = lifted_vertex.underlying_vertex.colors[0]
    if color_list.count(rgb[1:7])==0:
        color_list.append(rgb[1:7])

    if abs(lifted_vertex.center()) != 0:
        max_pt = (1+abs(lifted_vertex.radius())/abs(lifted_vertex.center()))*lifted_vertex.center()
        min_pt = (1-abs(lifted_vertex.radius())/abs(lifted_vertex.center()))*lifted_vertex.center()
    else:
        max_pt = lifted_vertex.radius()
        min_pt = -lifted_vertex.radius()

    max_on_sp = to_ball_model(max_pt)
    min_on_sp = to_ball_model(min_pt)
    theta = acos(inner_product(max_on_sp, min_on_sp))

    radius = tan(theta/2)
    vect = [max_on_sp[0]+min_on_sp[0], max_on_sp[1]+min_on_sp[1], max_on_sp[2]+min_on_sp[2]]
    point = scalar_mult(1/(cos(theta/2)),normalized_vector(vect))
    
    if radius > 0.00001:
        f.write("sphere {\n")
        f.write("<%f, %f, %f>\n" % (point[0], point[1], point[2]))
        f.write("%f\n" % (radius))
        f.write("texture{ pigment{ color rgb <%f,%f,%f>} finish { phong 1 }}\n" % (r_of_rgb(rgb), g_of_rgb(rgb), b_of_rgb(rgb)))
        f.write("}\n")
