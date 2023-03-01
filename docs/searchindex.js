Search.setIndex({"docnames": ["design", "design/element", "design/material", "design/node", "design/plotter", "design/solver", "design/system", "examples", "examples/beam_examples", "examples/continuum_examples", "examples/frame_examples", "examples/mixed_structures_examples", "examples/plate_examples", "examples/truss_examples", "implementation", "implementation/Beam2D_class", "implementation/DOF_class", "implementation/ElasticSection", "implementation/ElementPlotter3D_class", "implementation/ElementPlotter_class", "implementation/Element_class", "implementation/FiberMaterial", "implementation/FiberSection", "implementation/Frame2D_class", "implementation/LinearSolver_class", "implementation/LinearTriangle_class", "implementation/Material_class", "implementation/NewtonRaphsonSolver_class", "implementation/Node_class", "implementation/PlaneStrain", "implementation/PlaneStress", "implementation/PlateSection", "implementation/Plot_Support_Classes", "implementation/Plotter3D_class", "implementation/Plotter_class", "implementation/Sections", "implementation/Solver_class", "implementation/System_class", "implementation/Transformation_class", "implementation/Truss_class", "index"], "filenames": ["design.rst", "design/element.rst", "design/material.rst", "design/node.rst", "design/plotter.rst", "design/solver.rst", "design/system.rst", "examples.rst", "examples/beam_examples.rst", "examples/continuum_examples.rst", "examples/frame_examples.rst", "examples/mixed_structures_examples.rst", "examples/plate_examples.rst", "examples/truss_examples.rst", "implementation.rst", "implementation/Beam2D_class.rst", "implementation/DOF_class.rst", "implementation/ElasticSection.rst", "implementation/ElementPlotter3D_class.rst", "implementation/ElementPlotter_class.rst", "implementation/Element_class.rst", "implementation/FiberMaterial.rst", "implementation/FiberSection.rst", "implementation/Frame2D_class.rst", "implementation/LinearSolver_class.rst", "implementation/LinearTriangle_class.rst", "implementation/Material_class.rst", "implementation/NewtonRaphsonSolver_class.rst", "implementation/Node_class.rst", "implementation/PlaneStrain.rst", "implementation/PlaneStress.rst", "implementation/PlateSection.rst", "implementation/Plot_Support_Classes.rst", "implementation/Plotter3D_class.rst", "implementation/Plotter_class.rst", "implementation/Sections.rst", "implementation/Solver_class.rst", "implementation/System_class.rst", "implementation/Transformation_class.rst", "implementation/Truss_class.rst", "index.rst"], "titles": ["Program Design", "Element", "Material", "Node", "Plotter", "Solver", "System", "Example Problems", "Beam Examples", "Continuum Examples", "Frame Examples", "Mixed Structures Examples", "Plate Examples", "Truss Examples", "Implementation", "Beam2D class", "DoF class", "ElasticSection", "ElementPlotter3D class", "ElementPlotter class", "Element classes", "FiberMaterial material class", "FiberSection", "Frame2D class", "Linear Solver class", "LinearTriangle class", "Material class", "Newton-Raphson Solver class", "Node class", "PlaneStrain material class", "PlaneStress material class", "PlateSection", "Plot Support classes", "Plotter3D class", "Plotter class", "Section Material classes", "Solver class", "System class", "Transformation class", "Truss class", "Welcome to the FEM.edu documentation!"], "terms": {"The": [0, 4, 15, 17, 19, 20, 22, 23, 24, 27, 28, 31, 36, 37, 39], "goal": 0, "i": [0, 2, 3, 4, 5, 15, 17, 19, 20, 22, 23, 28, 31, 32, 33, 34, 35, 37, 38, 39], "creat": [0, 4, 5, 6, 19, 32, 33, 34, 37], "an": [0, 1, 5, 6, 20, 27, 32, 35, 36, 37, 38], "object": [0, 1, 4, 5, 6, 15, 23, 28, 33, 34, 37, 39], "orient": 0, "finit": [0, 19], "element": [0, 4, 5, 6, 14, 15, 19, 23, 24, 25, 27, 28, 32, 35, 36, 37, 38, 39, 40], "analysi": [0, 27, 36], "base": [0, 4, 19, 32, 35, 37], "method": [0, 14, 19, 24, 28, 32, 36, 37], "handl": [0, 6], "arbitrari": 0, "truss": [0, 1, 6, 7, 14, 20, 25, 40], "2d": [0, 15, 23, 29, 30, 39], "allow": 0, "load": [0, 3, 5, 6, 8, 9, 10, 11, 12, 13, 15, 20, 23, 24, 27, 28, 36, 37, 38], "can": [0, 20, 28, 32], "plot": [0, 4, 5, 6, 14, 18, 19, 20, 33, 34, 36, 37, 40], "its": [0, 1, 20], "undeform": [0, 4, 20, 37], "deform": [0, 3, 4, 15, 19, 20, 28, 32, 33, 34, 35, 37, 39], "shape": [0, 4, 20, 28, 36], "need": [0, 1, 3, 4, 19, 20, 27, 28, 36], "us": [0, 1, 5, 8, 9, 10, 11, 12, 13, 15, 19, 20, 23, 24, 27, 28, 32, 33, 34, 35, 36, 37, 38, 39], "follow": [0, 35, 39], "node": [0, 1, 4, 5, 6, 14, 15, 16, 19, 20, 23, 24, 27, 32, 36, 37, 38, 39, 40], "materi": [0, 1, 14, 15, 17, 20, 22, 23, 25, 31, 39, 40], "system": [0, 3, 4, 14, 19, 20, 24, 27, 28, 32, 33, 34, 36, 40], "solver": [0, 6, 14, 37, 40], "plotter": [0, 5, 6, 14, 18, 19, 32, 33, 40], "each": [1, 3, 4, 5], "instanc": [1, 3, 5, 6, 22, 31, 37], "repres": [1, 3, 4, 15, 16, 20, 21, 23, 25, 26, 28, 29, 30, 33, 34, 35, 37, 39], "one": [1, 3, 5, 6, 32, 37], "member": 1, "input": [1, 2, 3, 4, 5, 6, 13, 17, 22, 31], "return": [1, 2, 3, 4, 5, 6, 15, 20, 21, 23, 26, 28, 29, 30, 36, 37, 39], "descript": [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 28], "__init__": [1, 2, 3, 4, 5, 6], "nd0": 1, "nd1": 1, "two": [1, 3, 4], "constructor": [1, 2, 3, 4, 5, 6, 20], "getforc": [1, 5, 20], "list": [1, 3, 4, 5, 6, 15, 16, 19, 20, 23, 24, 27, 28, 32, 36, 37, 39], "1d": [1, 21], "np": [1, 3, 4, 5, 15, 23, 39], "arrai": [1, 3, 4, 5, 15, 19, 23, 28, 30, 32, 34, 37, 39], "A": [1, 2, 17, 19, 20, 21, 24, 26, 27, 28, 35, 36, 37, 38], "nodal": [1, 3, 5, 6, 15, 19, 20, 23, 28, 32, 33, 34, 38, 39], "forc": [1, 3, 4, 5, 6, 15, 17, 19, 20, 22, 23, 24, 27, 31, 32, 33, 34, 36, 39], "shall": [1, 4, 20, 28, 36, 38], "compon": [1, 3, 30, 31], "getstiff": [1, 2, 5, 20, 26], "2": [1, 2, 4, 5, 22, 30, 31, 39], "matrix": [1, 36], "contain": [1, 4, 5], "tangent": [1, 2, 15, 22, 23, 26, 31, 36, 39], "matric": 1, "note": [1, 3], "mai": [1, 15, 19, 23, 24, 35, 36, 38], "have": [1, 3, 32], "chang": [1, 15, 17, 22, 23, 31, 37], "state": [1, 2, 15, 19, 21, 26, 29, 30, 37, 39], "between": [1, 35, 39], "call": [1, 20, 27, 28], "so": [1, 3, 32], "you": 1, "recomput": [1, 5], "everi": [1, 20], "time": 1, "name": [1, 2, 3, 4, 5, 6], "type": [1, 2, 3, 4, 5, 6, 16, 20, 36], "either": [1, 35], "end": [1, 4, 30], "pointer": [1, 6, 15, 23, 24, 27, 28, 32, 36, 37, 39], "comput": [1, 5, 6, 15, 23, 31, 36, 39], "stress": [1, 2, 20, 26, 29, 30, 31], "modulu": [1, 2], "hold": [1, 2, 3, 5, 6], "p0": [1, 24, 27, 36], "p1": 1, "kt": [1, 15, 23, 39], "stiff": [1, 2, 5, 15, 17, 22, 23, 26, 31, 36, 39], "all": [1, 2, 5, 6, 8, 9, 10, 11, 12, 13, 20, 28, 32, 35, 36, 38], "equat": [1, 2, 3, 5, 37], "bf": [1, 3, 5, 24, 28, 36], "l": 1, "x": [1, 3, 4, 15, 23, 28], "_1": [1, 3], "_0": [1, 3], "ell": 1, "n": [1, 20, 21, 29, 30, 31], "frac": [1, 2, 30], "1": [1, 2, 3, 5, 6, 15, 19, 20, 21, 23, 26, 28, 29, 30, 32, 35, 36, 37], "strain": [1, 2, 17, 21, 22, 26, 29, 30, 31], "varepsilon": [1, 2, 17, 22, 31], "cdot": 1, "u": [1, 3, 4, 5, 6, 23, 24, 28, 39], "f": [1, 2, 5, 6, 17, 22, 28, 36], "sigma": [1, 2, 22, 26, 31], "setstrain": [1, 2, 26, 29], "ep": [1, 2, 21, 26, 29, 30], "getstress": [1, 2, 26], "vector": [1, 3, 4, 5, 15, 19, 20, 23, 24, 27, 28, 32, 34, 36, 37, 39], "p": [1, 5, 6, 23, 24, 30], "e": [1, 2, 8, 9, 10, 11, 12, 13, 17, 20, 21, 26, 28, 29, 30, 32, 35, 36, 37], "k": [1, 5, 24, 36], "e_t": [1, 2], "otim": 1, "find": [1, 37], "_": [1, 5, 31], "00": 1, "11": [1, 5], "01": 1, "10": [1, 2, 27], "thi": [2, 3, 4, 5, 20, 22, 24, 28, 31, 32, 36, 37], "provid": [2, 3, 15, 20, 21, 23, 26, 28, 29, 30, 32, 37], "demonstr": 2, "exampl": [2, 40], "paramet": [2, 15, 17, 19, 20, 21, 24, 26, 27, 28, 29, 30, 32, 33, 34, 36, 37], "0": [2, 3, 5, 6, 15, 19, 20, 21, 23, 26, 28, 29, 30, 31, 32, 35, 36, 37], "set": [2, 3, 4, 5, 20, 27, 28, 36, 37], "initi": [2, 3, 4, 24, 27, 28, 36, 37], "intern": [2, 3, 19, 20, 24, 27, 32, 33, 34, 36], "getarea": 2, "cross": [2, 22, 31, 35], "section": [2, 14, 15, 17, 22, 23, 26, 31], "area": 2, "from": [2, 5, 8, 9, 10, 11, 12, 13, 19, 20, 22, 24, 27, 31, 32, 35], "request": [2, 5, 6, 20, 26, 28, 36], "axial": [2, 10, 17, 21, 22, 23, 26, 29, 30, 31, 35, 39], "updat": [2, 5, 6, 21, 24, 26, 27, 29, 30, 36, 37], "user": [2, 4, 19, 20, 21, 26, 28, 29, 30, 38], "valu": [2, 4, 15, 19, 20, 21, 23, 26, 28, 29, 30, 32, 37, 39], "param": [2, 17, 21, 26, 29, 30, 35], "dict": 2, "default": [2, 3, 4, 19, 20, 23, 36, 37], "100": 2, "nu": [2, 21, 26, 29, 30, 35], "fy": [2, 21, 26, 29, 30, 35], "0e30": 2, "moe": 2, "poisson": 2, "": [2, 5, 15, 23, 24, 27, 28, 32, 36], "ratio": [2, 32], "yield": 2, "plastic_strain": 2, "float": [2, 3, 16, 39], "sig": 2, "current": [2, 3, 8, 9, 10, 11, 12, 13, 15, 16, 20, 24, 27, 28, 36, 37, 39], "et": [2, 26], "materil": 2, "elast": [2, 15, 17, 23, 35], "trial": 2, "varepsilon_p": 2, "check": 2, "f_y": 2, "IF": 2, "ge": 2, "3": 2, "delta": [2, 23], "text": [2, 37], "sign": 2, "y": [3, 4, 15, 23, 28], "coordin": [3, 4, 20, 28], "point": [3, 4, 8, 15, 20, 28], "posit": [3, 4, 5, 6, 15, 28, 39], "displac": [3, 4, 5, 6, 15, 16, 19, 23, 24, 27, 28, 32, 36, 37], "zero": [3, 5, 28, 36, 37], "fixdof": [3, 28], "idx": 3, "degre": [3, 16], "freedom": [3, 16], "dof": [3, 14, 15, 20, 23, 28, 37, 39, 40], "flag": 3, "accordingli": 3, "isfix": [3, 16, 28], "true": [3, 19, 27, 28, 32, 33, 34, 36, 37], "fals": [3, 4, 19, 20, 24, 27, 28, 32, 33, 34, 36, 37], "test": [3, 10], "function": [3, 19, 20, 27, 36, 37], "fix": [3, 5, 28], "otherwis": 3, "setdisp": [3, 5, 28], "v": [3, 4, 5, 15, 23, 39], "overwrit": 3, "getdisp": [3, 28], "getpo": [3, 28], "getdeformedpo": [3, 28], "factor": [3, 5, 6, 19, 20, 28, 32, 36, 37], "magnifi": 3, "would": 3, "good": 3, "none": [3, 4, 19, 20, 28, 32, 33, 34, 37], "given": [3, 4, 19, 20, 24, 28, 32, 33, 34, 37], "addload": [3, 28], "px": 3, "py": 3, "add": [3, 4, 5, 6, 19, 32, 33, 34, 37], "setload": [3, 28], "replac": [3, 4, 36], "getload": [3, 5, 20, 28], "po": 3, "index": [3, 5, 36, 37, 40], "int": [3, 16, 27], "addnod": [3, 5, 6, 37], "thisnod": 3, "team": 3, "disp": [3, 4, 5, 16, 32, 34], "fixiti": [3, 16, 38], "th": 3, "sensibl": 4, "setmesh": [4, 32, 34], "vert": [4, 34], "line": [4, 20, 34], "indic": [4, 5, 16, 19, 28, 32, 33, 34, 37], "self": [4, 15, 16, 20, 23, 28, 39], "vertic": 4, "inform": [4, 19, 27, 32], "setdisplac": [4, 32, 34], "setvalu": [4, 32, 34], "val": [4, 15, 23, 32, 34], "displacementplot": [4, 19, 32, 33, 34], "file": [4, 19, 32, 33, 34, 37], "string": [4, 8, 9, 10, 11, 12, 13, 19, 32, 37], "show": [4, 37], "black": 4, "model": [4, 5, 6, 19, 22, 23, 31, 37, 39], "red": 4, "If": [4, 5, 19, 20, 28, 32, 33, 34, 36, 37, 38], "save": [4, 16], "copi": 4, "valueplot": [4, 19, 32, 33, 34, 37], "color": [4, 19, 32, 33, 34, 37], "colormap": 4, "colorbar": 4, "legend": 4, "pair": 4, "start": [4, 20, 37], "respect": [4, 19, 36], "must": 4, "ident": 4, "entri": 4, "newnod": [5, 6, 37], "your": [5, 6], "addel": [5, 6, 37], "newelem": [5, 6], "solv": [5, 6, 24, 27, 36, 37], "assembl": [5, 6, 24, 27, 36], "k_t": [5, 6], "loop": [5, 6], "through": [5, 6, 19, 20, 31, 32], "unbalanc": [5, 6], "r": [5, 6, 19, 32, 34, 36], "collect": [5, 6], "info": [5, 6], "send": [5, 6, 28], "report": [5, 6, 28, 37], "print": [5, 6, 8, 9, 10, 11, 12, 13, 37], "summari": [5, 6, 37], "size": 5, "ha": [5, 28, 37], "node0": [5, 25], "elem": [5, 32], "node1": [5, 25], "j": [5, 15, 23, 39], "local": [5, 15, 20, 28], "d": [5, 22, 28, 31, 36], "o": [5, 28, 36], "belong": [5, 16], "global": [5, 20, 36], "4": 5, "m": [5, 17, 22, 31], "assembli": [5, 28], "sy": 5, "over": [5, 31], "should": [5, 20, 28], "0e20": 5, "step": 5, "5": 5, "6": 5, "do": [5, 28], "repeat": 5, "7": 5, "everyth": 5, "wa": 5, "done": 5, "correctli": 5, "support": [5, 14, 18, 19, 33, 34, 40], "reaction": 5, "free": 5, "numer": [5, 31], "beam": [7, 10, 15, 20, 32, 35, 37, 40], "frame": [7, 20, 23, 35, 40], "plate": [7, 8, 20, 31, 35, 40], "continuum": [7, 20, 40], "mix": [7, 36, 40], "structur": [7, 40], "ar": [8, 9, 10, 11, 12, 13, 15, 19, 20, 23, 28, 32, 35, 36], "packag": [8, 9, 10, 11, 12, 13], "To": [8, 9, 10, 11, 12, 13, 28], "run": [8, 9, 10, 11, 12, 13], "specif": [8, 9, 10, 11, 12, 13, 28, 32], "g": [8, 9, 10, 11, 12, 13, 28, 32], "femedu": [8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 21, 23, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37, 38, 39], "import": [8, 9, 10, 11, 12, 13], "beam01": 8, "ex": [8, 9, 10, 11, 12, 13], "examplebeam01": 8, "doc": [8, 9, 10, 11, 12, 13, 38], "actual": [8, 9, 10, 11, 12, 13, 32], "problem": [8, 9, 10, 11, 12, 13, 40], "singl": [8, 16, 20, 25, 27, 28, 39], "span": 8, "beam02": 8, "three": 8, "continu": 8, "solid": [9, 20], "solid01": 9, "examplesolid01": 9, "wip": [9, 11, 12], "frame01": 10, "exampleframe01": 10, "later": [10, 28, 37], "frame02": 10, "column": [10, 35], "rotat": [10, 15, 23, 28], "invari": 10, "frame03": 10, "build": [10, 36], "mixed01": 11, "examplemixed01": 11, "plate01": 12, "exampleplate01": 12, "truss01": 13, "exampletruss01": 13, "basic": 13, "triangl": [13, 20], "truss02": 13, "simpl": [13, 17, 20], "static": 13, "determin": 13, "bridg": 13, "truss03": 13, "explor": 13, "altern": 13, "format": [13, 19, 32, 33, 34, 37], "truss04": 13, "3d": [13, 32, 33, 39], "class": [14, 31, 40], "inherit": [14, 24, 27, 32], "deriv": 14, "beam2d": [14, 20], "frame2d": [14, 20], "lineartriangl": [14, 20], "fibermateri": [14, 22, 26], "planestress": [14, 26, 31], "planestrain": [14, 26], "plotter3d": [14, 32], "elementplott": [14, 32], "elementplotter3d": [14, 32], "linear": [14, 15, 17, 23, 27, 35, 36], "newton": [14, 36], "raphson": [14, 36], "transform": [14, 20, 28, 40], "nodei": [15, 23, 39], "nodej": [15, 23, 39], "assumpt": [15, 23, 24], "defin": [15, 20, 28, 31, 32, 38], "along": 15, "axi": [15, 19, 28, 32, 33, 34], "ignor": 15, "small": [15, 23], "navier": [15, 23], "bernoulli": [15, 23], "euler": [15, 23], "plane": [15, 23, 29, 30, 31, 35], "shear": [15, 23, 32, 37], "rigid": [15, 23], "futur": [15, 23, 36], "releas": [15, 23, 36], "ui": [15, 23, 28, 39], "theta": [15, 23], "rz": [15, 23, 28], "tupl": [15, 20, 23, 28, 37, 39], "_v_": [15, 23], "moment": [15, 17, 22, 23, 31, 32, 37], "_m_": [15, 23], "getinternalforc": [15, 23], "normal": [15, 23], "locat": [15, 23], "ndarrai": [15, 23, 28, 37], "which": [15, 19, 20, 23, 32, 33, 34, 38], "interv": [15, 23], "resetload": [15, 20, 23, 28, 37], "setdistload": 15, "w": [15, 39], "uniform": 15, "direct": [15, 28], "updatest": [15, 20, 21, 23, 25, 30, 39], "domain": [16, 28, 37, 38], "str": [16, 19, 32, 33, 34, 37], "nodeid": 16, "statu": 16, "bool": [16, 28], "record": [16, 37], "whether": [16, 24], "histori": [16, 37], "implement": [17, 19, 20, 22, 23, 28, 31, 32, 36, 40], "relat": 17, "ea": 17, "qquad": [17, 22, 31], "ei": 17, "phi": [17, 22], "where": [17, 20, 22, 30, 31, 36], "result": [17, 22, 24, 31], "output": [17, 22, 28, 31, 32], "curvatur": [17, 22, 31], "constant": 17, "flexur": [17, 22, 31, 35], "take": 19, "mesh": [19, 37], "potenti": 19, "shade": 19, "specifi": [19, 28, 38], "variabl": [19, 32, 37], "obtain": [19, 24, 31], "standard": 19, "version": 19, "those": [19, 20, 32], "overload": [19, 36], "addforc": [19, 32, 33, 34], "ax": [19, 32, 33, 34], "shown": [19, 32, 33, 34], "beamvalueplot": [19, 32, 37], "variable_nam": [19, 32], "filenam": [19, 32, 33, 34, 37], "show_arrow": 19, "store": [19, 28, 32, 33, 34, 37], "proper": [19, 32, 33, 34, 37], "extens": [19, 32, 33, 34, 37], "desir": [19, 32, 33, 34, 37], "png": [19, 32, 33, 34, 37], "pdf": [19, 32, 33, 34, 37], "code": [19, 20, 28, 32, 37], "scale": [19, 32], "modeshap": [19, 20, 28], "kwarg": [19, 20, 24, 27, 28, 34, 36, 37], "setreact": [19, 32, 34], "identifi": [19, 32, 33, 34], "magnitud": [19, 32, 33, 34], "particular": 20, "rather": 20, "common": [20, 28], "fem": 20, "edu": 20, "drawel": [20, 32], "declar": 20, "gener": [20, 26, 28, 36], "element_typ": 20, "mechan": 20, "work": [20, 27, 32], "unknown": 20, "mean": 20, "see": [20, 28, 36], "below": 20, "more": [20, 37], "detail": 20, "abstract": [20, 26, 32, 36], "addtransform": [20, 28], "t": [20, 28, 29, 30], "local_nod": 20, "attach": [20, 28, 38], "ani": [20, 28, 36, 37], "empti": 20, "hand": 20, "appli": [20, 28, 36, 37, 38], "non": 20, "go": 20, "number": [20, 27], "remov": 20, "assign": [20, 37], "getdof": 20, "driven": 20, "onli": [20, 27, 36], "apply_load_factor": [20, 28, 36], "like": 20, "surfac": 20, "bodi": 20, "reset": [20, 23, 36, 37], "setloadfactor": [20, 28, 36, 37], "lam": [20, 28, 36, 37], "target": [20, 28, 36, 37], "instead": [20, 28, 37], "enter": [20, 28, 36, 37], "pattern": [20, 28, 36, 37], "consid": [20, 28, 36, 37], "refer": [20, 24, 27, 28, 36, 37], "multipli": [20, 28, 36, 37], "scalar": [20, 28, 36, 37], "lambda": [20, 28, 36, 37], "explicitli": [20, 28, 36, 37], "assum": [20, 24, 28, 36, 37], "entir": [20, 28, 35, 36, 37], "full": [20, 27, 28, 36, 37], "draw": 20, "curv": 20, "nice": 20, "shell": [20, 35], "tetrahedron": 20, "quad": 20, "brick": 20, "seri": 20, "_x_": 20, "_y_": 20, "_z_": 20, "For": 20, "magnif": [20, 28, 32, 37], "drawbrick": 20, "drawcurv": 20, "drawlin": 20, "drawquad": 20, "drawtetrahedron": 20, "drawtriangl": 20, "1e": [21, 26, 29, 30, 35], "30": [21, 26, 29, 30, 35], "fiber": [21, 22], "tensor": [21, 26, 29, 30, 31], "yet": [22, 31], "avail": [22, 31], "complex": [22, 31], "built": [22, 31], "multipl": [22, 31], "z": [22, 28, 31, 32], "varepsilon_0": 22, "int_a": 22, "da": 22, "df": 22, "mathcal": [22, 31], "c": [22, 30, 31], "mathbb": [22, 31], "dm": [22, 31], "b": [22, 31], "distanc": [22, 31], "coupl": [22, 31, 35], "moder": 23, "ux": [23, 28, 39], "_f_": 23, "linearsolv": 24, "impli": 24, "force_onli": [24, 27, 36], "pushstat": [24, 27, 36], "push": [24, 27, 28, 36, 37], "data": [24, 27, 36, 37], "requir": [24, 27, 28, 32, 36], "pref": [24, 27, 36], "u1": [24, 27, 36], "converg": [24, 27, 36], "un": [24, 27, 36], "previou": [24, 27, 36], "lam1": [24, 27, 36], "level": [24, 27, 36, 37], "lamn": [24, 27, 36], "total": [24, 30], "verifi": 24, "correct": 24, "out": [24, 36], "equilibrium": 24, "experi": 24, "nonlinear": [24, 27], "behavior": 24, "under": 24, "node2": 25, "control": 27, "iter": 27, "newtonraphsonsolv": 27, "getresiduum": 27, "redesign": 27, "TO": [27, 28], "WITH": 27, "smart": 27, "max_step": 27, "verbos": [27, 36], "maximum": 27, "addit": [27, 32], "solvesinglestep": 27, "helper": 27, "perform": [27, 31, 36], "solut": [27, 32], "x0": 28, "y0": 28, "z0": 28, "prescrib": 28, "furthermor": 28, "veloc": 28, "acceler": 28, "arefix": 28, "routin": 28, "num_dof": 28, "restrain": 28, "also": 28, "caller": 28, "adjust": 28, "dure": 28, "nd": 28, "sequenc": 28, "same": 28, "order": 28, "dof_list": 28, "exist": 28, "forget": 28, "definit": 28, "kei": 28, "hasload": 28, "popu": [28, 37], "restor": [28, 37], "previous": [28, 37], "pushu": [28, 37], "individu": 28, "uz": [28, 39], "rx": 28, "about": 28, "ry": 28, "usual": 28, "sent": 28, "resetal": [28, 37], "resetdisp": [28, 37], "BE": 28, "mode": [28, 36], "left": 30, "begin": 30, "s_": 30, "xx": 30, "yy": 30, "xy": 30, "right": 30, "ccc": 30, "e_": 30, "2e_": 30, "ij": [30, 31], "plastic": 30, "varepsilon_": 31, "kl": 31, "0_": 31, "phi_": 31, "int_": 31, "h": 31, "dz": 31, "layer": 31, "thick": 31, "integr": [31, 35], "sum": 31, "further": 31, "dn": 31, "ijkl": 31, "membran": 31, "fundament": 32, "featur": 32, "In": 32, "abstractplott": 32, "tradit": [32, 37], "diagram": [32, 37], "link": 32, "get": [32, 37], "them": 32, "set_axes_equ": 32, "make": 32, "equal": 32, "sphere": 32, "appear": 32, "cube": 32, "etc": 32, "possibl": 32, "matplotlib": 32, "set_aspect": 32, "plt": 32, "gca": 32, "sourc": 32, "http": 32, "stackoverflow": 32, "com": 32, "question": 32, "13685386": 32, "unit": 32, "length": 32, "aspect": 32, "perpendicular": 35, "middl": 35, "That": 35, "directli": [35, 37], "without": 35, "sectionmateri": 35, "elasticsect": 35, "fibersect": 35, "platesect": 35, "interfac": 36, "describ": 36, "balanc": 36, "residuum": 36, "_t": 36, "most": 36, "special": 36, "while": 36, "awar": 36, "thei": 36, "ask": 36, "explic": 36, "pass": 36, "access": [36, 37], "residu": 36, "checkstabl": 36, "stabil": [36, 37], "mathop": 36, "det": 36, "less": 36, "than": 36, "25": 36, "min": 36, "lambda_i": 36, "eigenvalu": 36, "fetchstat": 36, "fetch": 36, "getbucklingmod": 36, "eigen": 36, "smallest": 36, "absolut": 36, "plotbucklingmod": [36, 37], "eigenmod": 36, "lambda_": 36, "mathtt": 36, "resetdisplac": 36, "resetforc": 36, "newel": 37, "fetchrecord": 37, "loadfactor": 37, "getsolv": 37, "give": 37, "instruct": 37, "initrecord": 37, "gather": 37, "been": 37, "mark": 37, "deprec": 37, "appropri": 37, "plotdof": 37, "plotsystem": 37, "recordthisstep": 37, "setsolv": 37, "upon": 37, "success": 37, "old": 37, "new": 37, "startrecord": 37, "stoprecord": 37, "stop": 37, "contour": 37, "select": 37, "when": 38, "dode": 38, "getaxialforc": 39, "tension": 39, "program": 40, "design": 40, "modul": 40, "search": 40, "page": 40}, "objects": {"femedu.domain": [[16, 0, 0, "-", "DoF"], [28, 0, 0, "-", "Node"], [37, 0, 0, "-", "System"], [38, 0, 0, "-", "Transformation"]], "femedu.domain.DoF": [[16, 1, 1, "", "DoF"]], "femedu.domain.Node": [[28, 1, 1, "", "Node"]], "femedu.domain.Node.Node": [[28, 2, 1, "", "addLoad"], [28, 2, 1, "", "addTransformation"], [28, 2, 1, "", "areFixed"], [28, 2, 1, "", "fixDOF"], [28, 2, 1, "", "getDeformedPos"], [28, 2, 1, "", "getDisp"], [28, 2, 1, "", "getLoad"], [28, 2, 1, "", "getPos"], [28, 2, 1, "", "hasLoad"], [28, 2, 1, "", "isFixed"], [28, 2, 1, "", "popU"], [28, 2, 1, "", "pushU"], [28, 2, 1, "", "request"], [28, 2, 1, "", "resetAll"], [28, 2, 1, "", "resetDisp"], [28, 2, 1, "", "resetLoad"], [28, 2, 1, "", "setDisp"], [28, 2, 1, "", "setLoad"], [28, 2, 1, "", "setLoadFactor"]], "femedu.domain.System": [[37, 1, 1, "", "System"]], "femedu.domain.System.System": [[37, 2, 1, "", "addElement"], [37, 2, 1, "", "addNode"], [37, 2, 1, "", "beamValuePlot"], [37, 2, 1, "", "fetchRecord"], [37, 2, 1, "", "getSolver"], [37, 2, 1, "", "initRecorder"], [37, 2, 1, "", "plot"], [37, 2, 1, "", "plotBucklingMode"], [37, 2, 1, "", "plotDOF"], [37, 2, 1, "", "plotSystem"], [37, 2, 1, "", "popU"], [37, 2, 1, "", "pushU"], [37, 2, 1, "", "recordThisStep"], [37, 2, 1, "", "report"], [37, 2, 1, "", "resetAll"], [37, 2, 1, "", "resetDisp"], [37, 2, 1, "", "resetLoad"], [37, 2, 1, "", "setLoadFactor"], [37, 2, 1, "", "setSolver"], [37, 2, 1, "", "solve"], [37, 2, 1, "", "startRecorder"], [37, 2, 1, "", "stopRecorder"], [37, 2, 1, "", "valuePlot"]], "femedu.domain.Transformation": [[38, 1, 1, "", "Transformation"]], "femedu.elements": [[15, 0, 0, "-", "Beam2D"], [20, 0, 0, "-", "DrawElement"], [20, 0, 0, "-", "Element"], [23, 0, 0, "-", "Frame2D"], [25, 0, 0, "-", "LinearTriangle"], [39, 0, 0, "-", "Truss"]], "femedu.elements.Beam2D": [[15, 1, 1, "", "Beam2D"]], "femedu.elements.Beam2D.Beam2D": [[15, 2, 1, "", "getInternalForce"], [15, 2, 1, "", "resetLoads"], [15, 2, 1, "", "setDistLoad"], [15, 2, 1, "", "updateState"]], "femedu.elements.DrawElement": [[20, 1, 1, "", "DrawElement"]], "femedu.elements.DrawElement.DrawElement": [[20, 2, 1, "", "draw"], [20, 2, 1, "", "drawBrick"], [20, 2, 1, "", "drawCurve"], [20, 2, 1, "", "drawLine"], [20, 2, 1, "", "drawQuad"], [20, 2, 1, "", "drawTetrahedron"], [20, 2, 1, "", "drawTriangle"]], "femedu.elements.Element": [[20, 1, 1, "", "Element"]], "femedu.elements.Element.Element": [[20, 2, 1, "", "addTransformation"], [20, 2, 1, "", "getDofs"], [20, 2, 1, "", "getForce"], [20, 2, 1, "", "getLoad"], [20, 2, 1, "", "getStiffness"], [20, 2, 1, "", "resetLoads"], [20, 2, 1, "", "setLoadFactor"], [20, 2, 1, "", "updateState"]], "femedu.elements.Frame2D": [[23, 1, 1, "", "Frame2D"]], "femedu.elements.Frame2D.Frame2D": [[23, 2, 1, "", "getInternalForce"], [23, 2, 1, "", "resetLoads"], [23, 2, 1, "", "updateState"]], "femedu.elements.LinearTriangle": [[25, 1, 1, "", "LinearTriangle"]], "femedu.elements.LinearTriangle.LinearTriangle": [[25, 2, 1, "", "updateState"]], "femedu.elements.Truss": [[39, 1, 1, "", "Truss"]], "femedu.elements.Truss.Truss": [[39, 2, 1, "", "getAxialForce"], [39, 2, 1, "", "updateState"]], "femedu.materials": [[17, 0, 0, "-", "ElasticSection"], [21, 0, 0, "-", "FiberMaterial"], [26, 0, 0, "-", "Material"], [29, 0, 0, "-", "PlaneStrain"], [30, 0, 0, "-", "PlaneStress"], [35, 0, 0, "-", "SectionMaterial"]], "femedu.materials.ElasticSection": [[17, 1, 1, "", "ElasticSection"]], "femedu.materials.FiberMaterial": [[21, 1, 1, "", "FiberMaterial"]], "femedu.materials.FiberMaterial.FiberMaterial": [[21, 2, 1, "", "updateState"]], "femedu.materials.Material": [[26, 1, 1, "", "Material"]], "femedu.materials.Material.Material": [[26, 2, 1, "", "getStiffness"], [26, 2, 1, "", "getStress"], [26, 2, 1, "", "setStrain"]], "femedu.materials.PlaneStrain": [[29, 1, 1, "", "PlaneStrain"]], "femedu.materials.PlaneStrain.PlaneStrain": [[29, 2, 1, "", "setStrain"]], "femedu.materials.PlaneStress": [[30, 1, 1, "", "PlaneStress"]], "femedu.materials.PlaneStress.PlaneStress": [[30, 2, 1, "", "updateState"]], "femedu.materials.SectionMaterial": [[35, 1, 1, "", "SectionMaterial"]], "femedu.plotter": [[32, 0, 0, "-", "AbstractPlotter"], [19, 0, 0, "-", "ElementPlotter"], [18, 0, 0, "-", "ElementPlotter3D"], [34, 0, 0, "-", "Plotter"], [33, 0, 0, "-", "Plotter3D"]], "femedu.plotter.AbstractPlotter": [[32, 1, 1, "", "AbstractPlotter"]], "femedu.plotter.AbstractPlotter.AbstractPlotter": [[32, 2, 1, "", "addForces"], [32, 2, 1, "", "beamValuePlot"], [32, 2, 1, "", "displacementPlot"], [32, 2, 1, "", "setDisplacements"], [32, 2, 1, "", "setMesh"], [32, 2, 1, "", "setReactions"], [32, 2, 1, "", "setValues"], [32, 2, 1, "", "set_axes_equal"], [32, 2, 1, "", "valuePlot"]], "femedu.plotter.ElementPlotter": [[19, 1, 1, "", "ElementPlotter"]], "femedu.plotter.ElementPlotter.ElementPlotter": [[19, 2, 1, "", "addForces"], [19, 2, 1, "", "beamValuePlot"], [19, 2, 1, "", "displacementPlot"], [19, 2, 1, "", "setReactions"], [19, 2, 1, "", "valuePlot"]], "femedu.plotter.ElementPlotter3D": [[18, 1, 1, "", "ElementPlotter3D"]], "femedu.plotter.Plotter": [[34, 1, 1, "", "Plotter"]], "femedu.plotter.Plotter.Plotter": [[34, 2, 1, "", "addForces"], [34, 2, 1, "", "displacementPlot"], [34, 2, 1, "", "setDisplacements"], [34, 2, 1, "", "setMesh"], [34, 2, 1, "", "setReactions"], [34, 2, 1, "", "setValues"], [34, 2, 1, "", "valuePlot"]], "femedu.plotter.Plotter3D": [[33, 1, 1, "", "Plotter3d"]], "femedu.plotter.Plotter3D.Plotter3d": [[33, 2, 1, "", "addForces"], [33, 2, 1, "", "displacementPlot"], [33, 2, 1, "", "valuePlot"]], "femedu.solver": [[24, 0, 0, "-", "LinearSolver"], [27, 0, 0, "-", "NewtonRaphsonSolver"], [36, 0, 0, "-", "Solver"]], "femedu.solver.LinearSolver": [[24, 1, 1, "", "LinearSolver"]], "femedu.solver.LinearSolver.LinearSolver": [[24, 2, 1, "", "assemble"], [24, 2, 1, "", "pushState"], [24, 2, 1, "", "solve"]], "femedu.solver.NewtonRaphsonSolver": [[27, 1, 1, "", "NewtonRaphsonSolver"]], "femedu.solver.NewtonRaphsonSolver.NewtonRaphsonSolver": [[27, 2, 1, "", "assemble"], [27, 2, 1, "", "getResiduum"], [27, 2, 1, "", "pushState"], [27, 2, 1, "", "solve"], [27, 2, 1, "", "solveSingleStep"]], "femedu.solver.Solver": [[36, 1, 1, "", "Solver"]], "femedu.solver.Solver.Solver": [[36, 2, 1, "", "assemble"], [36, 2, 1, "", "checkStability"], [36, 2, 1, "", "fetchState"], [36, 2, 1, "", "getBucklingMode"], [36, 2, 1, "", "initialize"], [36, 2, 1, "", "pushState"], [36, 2, 1, "", "reset"], [36, 2, 1, "", "resetDisplacements"], [36, 2, 1, "", "resetForces"], [36, 2, 1, "", "setLoadFactor"], [36, 2, 1, "", "solve"]]}, "objtypes": {"0": "py:module", "1": "py:class", "2": "py:method"}, "objnames": {"0": ["py", "module", "Python module"], "1": ["py", "class", "Python class"], "2": ["py", "method", "Python method"]}, "titleterms": {"program": 0, "design": 0, "element": [1, 2, 20], "class": [1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37, 38, 39], "method": [1, 2, 3, 4, 5, 6, 20], "variabl": [1, 2, 3, 4, 5, 6, 15, 23, 39], "materi": [2, 21, 26, 29, 30, 35], "node": [3, 28], "plotter": [4, 34], "solver": [5, 24, 27, 36], "system": [5, 6, 37], "exampl": [7, 8, 9, 10, 11, 12, 13], "problem": 7, "beam": 8, "avail": [8, 9, 10, 11, 12, 13], "continuum": 9, "frame": 10, "mix": 11, "structur": 11, "plate": 12, "truss": [13, 39], "implement": 14, "beam2d": 15, "parent": [15, 17, 18, 19, 21, 22, 23, 24, 25, 27, 29, 30, 33, 34, 35, 39], "doc": [15, 18, 19, 21, 23, 24, 25, 27, 29, 30, 33, 34, 39], "intern": [15, 23, 39], "dof": 16, "elasticsect": 17, "elementplotter3d": 18, "elementplott": 19, "inherit": 20, "deriv": [20, 26, 32, 35, 36], "fibermateri": 21, "fibersect": 22, "frame2d": 23, "linear": 24, "state": [24, 27, 36], "i": [24, 27, 36], "defin": [24, 27, 36], "dictionari": [24, 27, 36], "follow": [24, 27, 36], "content": [24, 27, 36, 40], "lineartriangl": 25, "newton": 27, "raphson": 27, "planestrain": 29, "planestress": 30, "platesect": 31, "plot": 32, "support": 32, "plotter3d": 33, "section": 35, "transform": 38, "welcom": 40, "fem": 40, "edu": 40, "document": 40, "indic": 40, "tabl": 40}, "envversion": {"sphinx.domains.c": 2, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 6, "sphinx.domains.index": 1, "sphinx.domains.javascript": 2, "sphinx.domains.math": 2, "sphinx.domains.python": 3, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx": 56}})