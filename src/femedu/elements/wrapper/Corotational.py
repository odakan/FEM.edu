"""
==============================================
4-node quadrilateral - small displacement
==============================================
incomplete bi-linear interpolation

.. code::

         1
      x     y
        x*y


"""
import numpy as np
from copy import deepcopy

from ..Element import *
from ...domain.Node import *

class CorotQuad(Element):
    """
    class: representing a corotational finite deformation element wrapper

    This element wrapper works as an element and wraps an existing linear element. 
    The resulting element is a finite deformation transformation of the wrapped linear
    element using the corotational formulation.

    * For 2D plate behavior, define nodes as two-dimensional nodes
    * For 3D membrane behavior, define nodes as three-dimensional nodes
    """

    def __init__(self, element):
        self.wrapped_element = element
        super().__init__(element.nodes, element.material)
        self.element_type = element.element_type
        self.createFaces()

        if element.nodes[0].getPos().size == 3:
            dof_list = ('ux','uy','uz')
            ndof = 3
        elif element.nodes[0].getPos().size == 2:
            dof_list = ('ux','uy')
            ndof = 2
        else:
            raise TypeError("dimension of nodes must be 2 or 3")

        self._requestDofs(dof_list)

        self.force    = element.force
        self.Forces   = element.Forces
        self.Kt       = element.Kt
        self.ndof = element.ndof
    
    def __str__(self):
        return str(self.wrapped_element)
    
    def setSurfaceLoad(self, face, pn, ps=0):
        self.wrapped_element.setSurfaceLoad(self, face, pn, ps=0)

    def resetLoads(self):
        self.wrapped_element.resetLoads()

    def updateState(self):
        self.wrapped_element.updateState()
        self.Forces = self.wrapped_element.Forces
        self.Kt = self.wrapped_element.Kt

        

    def computeSurfaceLoads(self):
        self.wrapped_element.computeSurfaceLoads()
        self.Loads = self.wrapped_element.Loads

    def getStress(self):
        return self.wrapped_element.Stress
