import sys
from collections import deque

#from multiprocessing import Pool
from threading import Thread

#from . import quad_edge_mesh
from .quad_edge_mesh import QEMesh

class AABBTree (object):
    """ An Axis-aligned bounding box tree.
    For speedy intersection queries!
    """

    def __init__(self, qemesh):
        """ Construct an AABB Tree for this mesh.
        qemesh is a quad_edge_mesh.QEMesh.
        """
        if not isinstance(qemesh, QEMesh):
            raise TypeError("qemesh must be of type QEMesh!")
        
        # Construct the tree
        self._tree = AABBNode(list(qemesh._faces))

    def collides_with(self, other_face):
        """ Search tree for collision """
        pt_max = other_face.max
        pt_min = other_face.min
        return self._tree.collides_with(pt_max, pt_min)

    def collides_with_orthogonal_plane(self, orientation, position):
        """ Check if this tree collides with orientation.
        Orientation must be 0, 1, or 2 corresponding to the axis
        parallel to the plane normal.
        position is the coordinate of this plane along orientation axis
        """
        if orientation not in range(0,3):
            raise TypeError("orientation must be 0, 1, or 2!")

        collision_pairs = []
        self._tree.collides_with_plane(orientation, position, collision_pairs)

        return collision_pairs

    def collides_with_tree(self, other_tree):
        """ Return a list of pairs of faces whose aabb's collide """
        if not isinstance(other_tree, AABBTree):
            raise TypeError("Can only collide with other AABBTree's")
        #return self._tree.collides_with_tree(other_tree._tree)
        collision_pairs = []
        self._tree.collides_with_tree(other_tree._tree, collision_pairs)

        return collision_pairs

    def collides_with_tree_mt(self, other_tree):
        collision_pairs = []
        threads = []
        results = list()
        threads.append(Thread(
                target=self._tree.left_node.left_node.collides_with_tree,
                args=(other_tree._tree, results)))
        threads.append(Thread(
                target=self._tree.left_node.right_node.collides_with_tree,
                args=(other_tree._tree, results)))
        threads.append(Thread(
                target=self._tree.right_node.left_node.collides_with_tree,
                args=(other_tree._tree, results)))
        threads.append(Thread(
                target=self._tree.right_node.right_node.collides_with_tree,
                args=(other_tree._tree, results)))

        for t in threads:
            t.start()
        for t in threads:
            t.join()

        ret = []
        while len(results) != 0:
            ret.extend(results.pop())
        return ret

    def update_bbs(self):
        """ Update this tree's AABB's. This will result
        in a higher degree of overlap.
        """
        self._tree.update_bbs()
        
    def update_bbs_mt(self):
        """ Try to use parallelization to speed up the bounding box search.
        """
        threads = []
        threads.append(Thread(target=self._tree.left_node.left_node.update_bbs,
                              args=()))
        threads.append(Thread(target=self._tree.left_node.right_node.update_bbs,
                              args=()))
        threads.append(Thread(target=self._tree.right_node.left_node.update_bbs,
                              args=()))
        threads.append(Thread(target=self._tree.right_node.right_node.update_bbs,
                              args=()))

        for t in threads:
            t.start()
        for t in threads:
            t.join()
        
        # Assemble final bb's
        self._tree.left_node.update_bb_from_children()
        self._tree.right_node.update_bb_from_children()
        self._tree.update_bb_from_children()
        

class AABBNode (object):
    """ A node of an Axis-aligned bounding box tree. """

    def __init__(self, faces):
        """ Initialize this AABB Node. Faces must be a list of type QEFace

        Find longest axis of the AABB and split points along the halfway point
        into two new nodes, left and right.

        Leaves hold a pointer to the face's bounding box values so it will break
        if the point is replaced
        """
        
        if len(faces) is 1:
            self.is_leaf = True
            self.leaf = faces[0]
            self.left_node = None
            self.right_node = None
            self.max_pt = faces[0].max_coord
            self.min_pt = faces[0].min_coord
            return
        
        self.is_leaf = False
        # First compute BB for this set of faces
        sysmin = sys.float_info.min
        sysmax = sys.float_info.max
        self.min_pt = [sysmax, sysmax, sysmax]
        self.max_pt = [sysmin, sysmin, sysmin]
        
        # It's actually faster to sort 3 times and use
        # as an ESTIMATE of the max/min than to iterate through
        # all faces.
        maxd = 0
        maxi = -1

        fs = [None, None, None]
        # faces.sort(key=lambda face:face.min[0])
        fs[0] = sorted(faces, key=lambda face:face.min_coord[0])
        d = faces[-1].max_coord[0] - faces[0].min_coord[0]
        if d > maxd:
            maxd = d
            maxi = 0
        # faces.sort(key=lambda face:face.min[1])
        fs[1] = sorted(faces, key=lambda face:face.min_coord[1])
        d = faces[-1].max_coord[1] - faces[0].min_coord[1]
        if d > maxd:
            maxd = d
            maxi = 1
        # faces.sort(key=lambda face:face.min[2])
        fs[2] = sorted(faces, key=lambda face:face.min_coord[2])

        d = faces[-1].max_coord[2] - faces[0].min_coord[2]
        if d > maxd:
            maxd = d
            maxi = 2

        # Sort according to minimum index. It is faster
        # to do in-place sort because faces might already be sorted
        # in this direction.
        # faces.sort(key=lambda face:face.min[maxi])
        # sorted_faces = faces
        sorted_faces = fs[maxi]
        
        split_index = len(sorted_faces)//2
        left_sorted_faces = sorted_faces[:split_index]
        right_sorted_faces = sorted_faces[split_index:]

        self.left_node = AABBNode(left_sorted_faces)
        self.right_node = AABBNode(right_sorted_faces)

        self.update_bb_from_children()
        
    # def is_leaf(self):
    #     return self.left_node is None and self.right_node is None

    def collides_with(self, pt_max, pt_min):
        if (self.min_pt[0] > pt_max[0] or self.max_pt[0] < pt_min[0]): 
            return False
        if (self.min_pt[1] > pt_max[1] or self.max_pt[1] < pt_min[1]): 
            return False
        if (self.min_pt[2] > pt_max[2] or self.max_pt[2] < pt_min[2]): 
            return False
            
        if (self.left_node == None and self.right_node == None):
            return True
        
        return (self.left_node.collides_with(pt_max, pt_min) or
                self.right_node.collides_with(pt_max, pt_min))

    def collides_with_plane(self, orientation, position, pairs):
        """ Orientation must be 0, 1, or 2 (unchecked).
        position is a float.
        Checks if this bounding box and its children are intersected
        by this plane.
        pairs is a list of QEFaces (the leaves of this tree)
        If yes, appends leaves to pairs
        """
        if (self.min_pt[orientation] > position or
            self.max_pt[orientation] < position):
            return

        if self.is_leaf:
            pairs.append(self.leaf)
            return

        self.left_node.collides_with_plane(orientation, position, pairs)
        self.right_node.collides_with_plane(orientation, position, pairs)
        return

    def collides_with_tree(self, other_tree, pairs):
        if (self.min_pt[0] > other_tree.max_pt[0] or
            self.max_pt[0] < other_tree.min_pt[0]): 
            return
        if (self.min_pt[1] > other_tree.max_pt[1] or
            self.max_pt[1] < other_tree.min_pt[1]): 
            return
        if (self.min_pt[2] > other_tree.max_pt[2] or
            self.max_pt[2] < other_tree.min_pt[2]): 
            return

        if self.is_leaf:
            if other_tree.is_leaf:
                pairs.append((self.leaf, other_tree.leaf))
                return
            else:
                self.collides_with_tree(other_tree.left_node, pairs)
                self.collides_with_tree(other_tree.right_node, pairs)
                return
        elif other_tree.is_leaf:
            self.left_node.collides_with_tree(other_tree, pairs)
            self.right_node.collides_with_tree(other_tree, pairs)
            return
        else:
            self.left_node.collides_with_tree(other_tree.left_node, pairs)
            self.left_node.collides_with_tree(other_tree.right_node, pairs)
            self.right_node.collides_with_tree(other_tree.left_node, pairs)
            self.right_node.collides_with_tree(other_tree.right_node, pairs)
            return
        
    def collides_with_tree_append(self, other_tree, results=None):
        if (self.min_pt[0] > other_tree.max_pt[0] or
            self.max_pt[0] < other_tree.min_pt[0]): 
            return []
        if (self.min_pt[1] > other_tree.max_pt[1] or
            self.max_pt[1] < other_tree.min_pt[1]): 
            return []
        if (self.min_pt[2] > other_tree.max_pt[2] or
            self.max_pt[2] < other_tree.min_pt[2]): 
            return []

        result = []
        if self.is_leaf:
            if other_tree.is_leaf:
                result = [(self.leaf, other_tree.leaf)]
            else:
                result = self.collides_with_tree_append(other_tree.left_node) + \
                    self.collides_with_tree_append(other_tree.right_node)
        elif other_tree.is_leaf:
            result = self.left_node.collides_with_tree_append(other_tree) + \
                self.right_node.collides_with_tree_append(other_tree)
        else:
            result = self.left_node.collides_with_tree_append(other_tree.left_node) + \
                self.left_node.collides_with_tree_append(other_tree.right_node) + \
                self.right_node.collides_with_tree_append(other_tree.left_node) + \
                self.right_node.collides_with_tree_append(other_tree.right_node)
            
        if results is not None:
            results.extend(result)

        return result

    def update_bbs(self):
        """ Recursively update this bounding box and recurisvely call
        on children. """

        if self.is_leaf:
            # self.max_pt = self.leaf.max_coord
            # self.min_pt = self.leaf.min_coord

            ### DEBUG
            # print ("maxpos: (" + str(self.max_pt) +")")
            # print ("minpos: (" + str(self.min_pt) +")")
            ### /DEBUG
            return self.leaf.bb_updated

        left_updated = self.left_node.update_bbs()
        right_updated = self.right_node.update_bbs()
        
        if (not left_updated) and (not right_updated):
            return False
        
        self.update_bb_from_children()

        return True
        
    def update_bb_from_children(self):
        """ Update the bounding box of this node without descending tree.
        """
        self.min_pt[0] = min(self.left_node.min_pt[0],
                             self.right_node.min_pt[0])
        self.min_pt[1] = min(self.left_node.min_pt[1],
                             self.right_node.min_pt[1])
        self.min_pt[2] = min(self.left_node.min_pt[2],
                             self.right_node.min_pt[2])
        self.max_pt[0] = max(self.left_node.max_pt[0],
                             self.right_node.max_pt[0])
        self.max_pt[1] = max(self.left_node.max_pt[1],
                             self.right_node.max_pt[1])
        self.max_pt[2] = max(self.left_node.max_pt[2],
                             self.right_node.max_pt[2])
        
        # Rolled up loop below
        # for i in range(0,3):
        #     self.min_pt[i] = min(self.left_node.min_pt[i],
        #                          self.right_node.min_pt[i])
        #     self.max_pt[i] = max(self.left_node.max_pt[i],
        #                          self.right_node.max_pt[i])

