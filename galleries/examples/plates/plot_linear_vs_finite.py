"""
===============================================================
A square patch made of one quadrilateral plate elements
===============================================================

Basic implementation test with applied loads.
    Testing the tangent stiffness computation.

.. code::

    free   free
     ^     ^
     |     |
     3-----2 -> free
     |     | >
     |  a  | > (w = 1.0)
     |     | >
     0-----1 -> free

    width:  10.
    height: 10.

    Material parameters: St. Venant-Kirchhoff, plane stress
        E  = 10.0
        nu =  0.30
        t  =  1.0

    Element loads:
        node 0: [ 0.0, 0.0]
        node 1: [ 5.0, 0.0]
        node 2: [ 5.0, 0.0]
        node 3: [ 0.0, 0.0]

Author: Peter Mackenzie-Helnwein
"""

from femedu.examples import Example
from femedu.domain import System, Node
from femedu.solver import NewtonRaphsonSolver
from femedu.elements.linear import Truss
from femedu.elements.wrapper import CorotQuad
from femedu.elements.linear import Quad as Quad_lin
from femedu.elements.finite import Quad as Quad_fin
from femedu.materials import PlaneStress, FiberMaterial


class ExampleLinearVsFinite(Example):

    # sphinx_gallery_start_ignore
    def docString(self):
        s = """
    ## A square patch made of one quadrilateral plate elements

    Basic implementation test with applied loads.
    Testing the tangent stiffness computation. 

    free   free
     ^     ^
     |     |
     3-----2 -> free
     |     | >
     |  a  | > (w = 1.0)
     |     | >
     0-----1 -> free
    
    width:  10.
    height: 10.
    
    Material parameters: St. Venant-Kirchhoff, plane stress
        E  = 10.0
        nu =  0.30
        t  =  1.0
        
    Author: Peter Mackenzie-Helnwein 
    """
        return s

    # sphinx_gallery_end_ignore
    def problem(self):

        params2D = dict(
            E  = 10., # Young's modulus
            nu = 0.3,   # Poisson's ratio
            t  = 1.0,   # thickness of the plate
            fy = 1.e30  # yield stress
        )

        params1D = {'E': 10., 'A': 1., 'nu': 0.0, 'fy': 1.e30}

        a = 10.     # length of the plate in the x-direction
        b = 10.     # length of the plate in the y-direction
        c = 1.5

        model = System()
        model.setSolver(NewtonRaphsonSolver())


        # linear plate nodes
        nd0 = Node( 0.0, 0.0)
        nd1 = Node(   a, 0.0)
        nd2 = Node(   a,   b)
        nd3 = Node( 0.0,   b)
        nd4 = Node(   a,  -c)

        nd0.fixDOF('ux', 'uy')
        nd4.fixDOF('ux', 'uy')

        # finite plate nodes
        x_offset = a + 10
        nd10 = Node( 0.0 + x_offset, 0.0)
        nd11 = Node(   a + x_offset, 0.0)
        nd12 = Node(   a + x_offset,   b)
        nd13 = Node( 0.0 + x_offset,   b)
        nd14 = Node(   a + x_offset,  -c)

        nd10.fixDOF('ux', 'uy')
        nd14.fixDOF('ux', 'uy')

        # corotational plate nodes
        x_offset *= 2
        nd20 = Node( 0.0 + x_offset, 0.0)
        nd21 = Node(   a + x_offset, 0.0)
        nd22 = Node(   a + x_offset,   b)
        nd23 = Node( 0.0 + x_offset,   b)
        nd24 = Node(   a + x_offset,  -c)

        nd20.fixDOF('ux', 'uy')
        nd24.fixDOF('ux', 'uy')

        model.addNode(nd0, nd1, nd2, nd3, nd4, nd10, nd11, nd12, nd13, nd14, nd20, nd21, nd22, nd23, nd24)


        # small strain model
        plateA = Quad_lin(nd0, nd1, nd2, nd3, PlaneStress(params2D))
        trussA = Truss(nd4, nd1, FiberMaterial(params1D))

        # large deformation model with total lagrangian formulation
        plateB = Quad_fin(nd10, nd11, nd12, nd13, PlaneStress(params2D))
        trussB = Truss(nd14, nd11, FiberMaterial(params1D))

        # large deformation model with corotational formulation
        plateC = Quad_lin(nd20, nd21, nd22, nd23, PlaneStress(params2D))
        plateC_corot = CorotQuad(plateC)
        trussC = Truss(nd24, nd21, FiberMaterial(params1D))

        model.addElement(plateA, trussA, plateB, trussB, plateC_corot, trussC)


        # add loads
        # .. load only the trusses
        nd1.setLoad([0.0, 1.0], ('ux', 'uy'))
        nd11.setLoad([0.0, 1.0], ('ux', 'uy'))
        nd21.setLoad([0.0, 1.0], ('ux', 'uy'))

        #plateA.setSurfaceLoad(face=3, pn=1.0)
        #plateB.setSurfaceLoad(face=1, pn=1.0)

        model.plot(factor=0.0, title="Undeformed system", filename="linear_vs_finite_undeformed.png", show_bc=1)

        model.setLoadFactor(50.0)
        model.solve()

        # %%
        # The output shows that we do have a quadratic rate of convergence.

        # %%
        # Let's finish off with a nice plot of the deformed system.

        model.plot(factor=1.0, filename="linear_vs_finite_deformed.png")

        #model.report()


# %%
# Run the example by creating an instance of the problem and executing it by calling :py:meth:`Example.run()`
#

if __name__ == "__main__":
    ex = ExampleLinearVsFinite()
    ex.run()


