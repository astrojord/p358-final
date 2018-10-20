"""Implementation of Barnes-Hut tree algorithm (quadtree, octatree...) in n-dimensions.

Barnes-Hut algorithm is an approximation algorithm for performing an n-body simulation.
BHTree generation recursively divides n-dim space into cells, which contain 0 or 1 bodies.
This algorithm is used to approximate forces acting on a body. Group of bodies sufficently away
from queried body can be approximated to one center of mass.

See:
    https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation
    http://arborjs.org/docs/barnes-hut

See example implementation:
    https://codereview.stackexchange.com/questions/43749/barnes-hut-n-body-quadtree?newreg=e5c75739678d47d58bc1963b43ddd2e4

Terminology:
    1.  Node            - Basic element of BHTree structure. Sector can be either:
        a.  Empty       - Doesn't contain any bodies,
        b.  External    - Contains precisely one body,
        c.  Internal    - Contains 2**ndim children nodes.
    2.  Body            - Element fit into BHTree topology. It has position in n-dimensions and mass.
"""

import numpy as np


class Node:
    """Node is a basic element of a BHTree. It is described by its:
        a. position     - position of center of the cell in n-dimensional space,
        b. length       - cell extends by length/2 in every direction.
        c. type         - EMPTY, EXTERNAL, INTERNAL
        d. body         - for EXTERNAL nodes only.
        e. children     - for INTERNAL nodes only. Children are stored in an array with length of 2**ndim.
        f. mass,
        g. center of mass.
    """

    def __init__(self, pos, length):
        if isinstance(pos, list):
            pos = np.array(pos)
        if isinstance(pos, tuple):
            pos = np.array(list(pos))
        assert isinstance(pos, np.ndarray), "Position should be either a numpy.ndarray, list, or a tuple."

        assert (isinstance(length, float) or isinstance(length, int)), "Length should be either a float, or int."

        self.pos = pos
        self.length = length
        self.type = "EMPTY"
        self.body = None
        self.children = None
        self.com = pos
        self.mass = 0

    def fit(self, body):
        """Fits a body into the node.
        Recognized inputs:
            ndbh.Body                - body object      + list of multiple ndbh.Body objects
            [list, float]            - position, mass   + list of multiple [list, float] objects
            (list, float)            - position, mass   + list of multiple (list, float) objects
            [numpy.ndarray, float]   - position, mass   + list of multiple [numpy.ndarray, float] objects
            (numpy.ndarray, float)   - position, mass   + list of multiple (numpy.ndarray, float) objects
        """

        # input sanitation:
        if isinstance(body, Body):
            bodies = [body]
        elif isinstance(body, tuple):
            bodies = [Body(body[0], body[1])]
        elif isinstance(body, list):
            try:
                if isinstance(body[1], float):
                    bodies = [Body(body[0], body[1])]
                else:
                    bodies = []
                    for body in body:
                        if isinstance(body, Body):
                            bodies.append(body)
                        else:
                            bodies.append(Body(body[0], body[1]))
            except IndexError:
                raise AssertionError("Body format not recognized.")
        else:
            raise AssertionError("Body format not recognized.")
        # input sanitation END

        for body in bodies:

            assert len(body.pos) == len(self.pos), "Body and node dimensionality don't match."

            # first, check for out of bounds
            bounds_max = self.pos + self.length * 0.5
            bounds_min = self.pos - self.length * 0.5
            if any(body.pos > bounds_max) or any(body.pos < bounds_min):
                raise AssertionError("Body is out of bounds!")

            def child_node_index(body):
                """Returns an index of a child node from self.children to put body into"""

                # evaluate position of body relative to node's center
                ndim = len(self.pos)
                relative_pos = np.array(body.pos > self.pos, dtype=int)
                multiplier = np.array([2 ** (ndim - 1 - i) for i in range(ndim)])
                index = sum(relative_pos * multiplier)
                return index

            if self.type == "EMPTY":
                self.type = "EXTERNAL"
                self.body = body

            elif self.type == "EXTERNAL":

                # first check if new body has the same position as the occupant
                if np.array_equal(self.body.pos, body.pos):
                    self.body += body
                    return

                # DIVIDE SELF
                # calculate new centers
                ndim = len(self.pos)
                offset = self.length * 0.25
                centers = []
                for i in range(2 ** ndim):
                    s = np.binary_repr(i, width=ndim)  # creates strings like '000', '010', '111' (for ndim=3)
                    pos = self.pos + [(lambda c: offset if c == '1' else -offset)(c) for c in s]
                    centers.append(pos)

                self.children = [Node(i, self.length * 0.5) for i in centers]

                # find new place for occupant body
                idx = child_node_index(self.body)
                self.children[idx].fit(self.body)
                self.body = None
                self.type = "INTERNAL"

            if self.type == "INTERNAL":
                idx = child_node_index(body)
                try:
                    self.children[idx].fit(body)
                except RecursionError:
                    # just add to existing body
                    self.children[idx] += body

    def summary(self, include_empty=False, _final=True):
        """Returns node and all its children in a dictionary form. For debugging / un-black-boxing purposes."""

        return_dict = {'type': self.type, 'pos': str(self.pos.tolist())}

        if self.type != "EMPTY":
            return_dict['center_of_mass'] = str(self.com.tolist())
            return_dict['mass'] = self.mass
            return_dict['length'] = self.length

        if self.type == "INTERNAL":
            children = []
            for child in self.children:
                if child is None:
                    continue
                if (not include_empty) & (child.type == "EMPTY"):
                    continue
                children.append(child.summary(_final=False, include_empty=include_empty))
            return_dict['children'] = children

        if _final:
            import json
            return json.dumps(return_dict, indent=4)
        else:
            return return_dict

    def calculate_coms(self):
        """Calculates centers of mass for all nodes."""

        nodes = self._get_all_nodes()
        from operator import attrgetter
        sorted_nodes = sorted(nodes, key=attrgetter("length"))
        for node in sorted_nodes:
            node._calculate_center_of_mass()

    def _get_all_nodes(self):
        """Used for calculate_coms(). Returns node and all its children's nodes."""

        nodes = []
        if self.type == "INTERNAL":
            for child in self.children:
                nodes += child._get_all_nodes()

        nodes.append(self)
        return nodes

    def _calculate_center_of_mass(self):
        """Used for calculate_coms(). Calculates a center of mass of one node."""

        if self.type == "EMPTY":
            self.com = self.pos
            self.mass = 0
        elif self.type == "EXTERNAL":
            self.com = self.body.pos
            self.mass = self.body.mass
        else:
            sum_pos = np.zeros(len(self.pos))
            sum_mass = 0
            for child in self.children:
                if child.type == "EMPTY": continue
                if (child.mass == 0) & (child.type == "EXTERNAL"):
                    if child.occupant.mass != 0:
                        print("Error: Child seems to have wrongly calculated mass/center of mass. Recalculating.")
                        child._calculate_center_of_mass()
                sum_pos += child.com * child.mass
                sum_mass += child.mass
            self.com = sum_pos / sum_mass
            self.mass = sum_mass

    def neighbors(self, body, theta=0.75):
        """Returns a list of (position = numpy.ndarray, mass = float) tuples of bodies/nodes affecting a given body.
        Distance is controlled by theta. Lower theta = faster search = less accurate.
        Recognized inputs:
            ndbh.Body                - body object
            list                     - position
            numpy.ndarray            - position
            [list, float]            - position, mass
            (list, float]            - position, mass
            [numpy.ndarray, float]   - position, mass
            (numpy.ndarray, float]   - position, mass
        """

        # input sanitation:
        if (isinstance(body, list)) or (isinstance(body, tuple)):
            if isinstance(body[0], list) or (isinstance(body[0], np.ndarray)):
                assert len(body) == 2, "Body format not recognized."
                body = Body(body[0], body[1])
            else:
                body = Body(body, 0)
        elif isinstance(body, np.ndarray):
            body = Body(body, 0)
        assert isinstance(body, Body), "Body format not recognized."
        # input sanitation END

        neighbors = []
        if self.type == "EXTERNAL":
            if self.body == body:
                pass
            neighbors = [(self.com, float(self.mass))]
        elif self.type == "INTERNAL":
            dist = np.linalg.norm(body.pos - self.com)
            if self.length / dist < theta:
                neighbors = [(self.com, float(self.mass))]
            else:
                for child in self.children:
                    neighbors += child.neighbors(body=body, theta=theta)

        return neighbors

    def __repr__(self):
        return "<ndbh.Node: %s at %s, length: %d>" % (self.type, self.pos, self.length)


class Body:
    """Body is an object populating Nodes. It is described by its:
        a. position     - Position in n-dimensional space,
        b. mass.
        """

    def __init__(self, pos, mass):
        if isinstance(pos, list):
            pos = np.array(pos)
        if isinstance(pos, tuple):
            pos = np.array(list(pos))
        assert isinstance(pos, np.ndarray), "Position should be either a numpy.ndarray, list, or a tuple."
        assert (isinstance(mass, float) or isinstance(mass, int)), "Mass should be either a float, or int."

        self.pos = pos
        self.mass = mass

    def __eq__(self, other):
        return np.array_equal(self.pos, other.pos) and self.mass == other.mass

    def __add__(self, other):
        if isinstance(other, self.__class__):
            return Body(self.pos, self.mass + other.mass)
        else:
            raise TypeError("unsupported operand type(s) for +: '{}' and '{}'").format(self.__class__, type(other))

    def __repr__(self):
        return "<ndbh.Body: %s, mass: %d>" % (self.pos, self.mass)
