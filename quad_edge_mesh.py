from time import time
import sys

from . import obj_reader

class QEMeshBuilder (object):
    """ Defalut mesh builder.
    With a list of verts and faces, construct a Quad Edge Mesh
    """
    @classmethod
    def construct(cls, verts, faces):
        """ Constructs a new QEMesh with verts and faces.

        verts is a list of vertex tuples (x,y,z) and
        faces is a list of lists of indices which point to vertices.
        faces does not need to be all triangles, or all uniform.
        """
        qem = QEMesh()
        
        for i in range(0, len(verts)):
            v = verts[i]
            qev = QEVertex(qem, i)
            qev.pos[0] = v[0]
            qev.pos[1] = v[1]
            qev.pos[2] = v[2]
            
            qem.add_vertex(qev)

        eidx = 0 # an index for edges
        for fidx in range(0, len(faces)):
            f = faces[fidx]
            qef = QEFace(qem, fidx)

            nvs = len(f)
            for ptidx in range(0, nvs):
                vidx = f[ptidx]-1
                if ptidx+1 is nvs:
                    vidx_1 = f[0]-1
                else:
                    vidx_1 = f[ptidx+1]-1
                qef.verts.append(qem.get_vertex(vidx))

                qee = qem.get_edge_by_verts(vidx, vidx_1)
                if qee is None:
                    qee = QEEdge(qem, eidx)
                    eidx = eidx + 1
                    qee.b_vert = qem.get_vertex(vidx)
                    qee.t_vert = qem.get_vertex(vidx_1)
                    qee.l_face = qef
                    qem.add_edge(qee)
                else:
                    if qee.r_face is not None:
                        raise ValueError("This edge already has two faces.")
                    qee.r_face = qef

                qef.edges.append(qee)
                
            for face_eidx in range(0, len(qef.edges)):
                # Update the rest of the edges internal pointers 
                qee = qef.edges[face_eidx]

                if qee.r_face is None:
                    if face_eidx+1 is len(qef.edges):
                        qee.tl_edge = qef.edges[0]
                        qee.bl_edge = qef.edges[0]
                    else:
                        qee.tl_edge = qef.edges[face_eidx+1]
                        qee.bl_edge = qef.edges[face_eidx+1]
                else:
                    # else this is the right face
                    if face_eidx+1 is len(qef.edges):
                        qee.br_edge  = qef.edges[0]
                    else:
                        qee.br_edge  = qef.edges[face_eidx+1]
                    if face_eidx == 0:
                        qee.tr_edge = qef.edges[-1]
                    else:
                        qee.tr_edge = qef.edges[face_eidx-1]
            if len(qef.edges) is not len(qef.verts):
                raise ValueError("Uh oh")
            qem.add_face(qef)
            
        return qem

class QEMesh (object):
    """ A Quad-edge mesh.
    """

    def __init__(self):
        self._edges = []
        self._vertices = []
        self._faces = []
        self._existing_edges = {}

    def update_bounding_boxes(self):
        for face in self._faces:
            face.update_bounding_box()

    def get_vertex(self, vidx):
        """ return a QEVertex by index.
        """
        return self._vertices[vidx]

 
    def get_edge_by_idx(self, eidx):
        """ Get a QEEdge from this QEMesh by the index in the internal
        list of QEEdges.
        """
        return self._edges[eidx]

    def get_edge_by_verts(self, vidx1, vidx2):
        """ Get a QEEdge from this QEMesh by the two QEVertex Indices
        which the edge connects, or None if not existent.

        
        """
        if (vidx1, vidx2) in self._existing_edges:
            return self._existing_edges[(vidx1, vidx2)]
        elif (vidx2, vidx1) in self._existing_edges:
            return self._existing_edges[(vidx2, vidx1)]
        else:
            return None
        # try:
        #     return self._existing_edges[(vidx1, vidx2)]
        # except KeyError:
        #     try:
        #         return self._existing_edges[(vidx2, vidx1)]
        #     except KeyError:
        #         return None

    def get_face(self, fidx):
        """ Get a QEFace from this mesh by index.
        """
        return self._faces[fidx]

    def add_face(self, face):
        """ Add a QEFace to this QEMesh.

        Raises TypeError if face is not QEFace.
        """
        if not isinstance(face, QEFace):
            raise TypeError("face must be QEFace")
        
        self._faces.append(face)

    def add_edge(self, edge):
        """ Add a QEEdge to this QEMesh.

        Raises TypeError if not QEEdge.
        Raises ValueError if edge already exists in this mesh.

        Internally, edge is appended to _edges list, and added
        to dictionary _existing_edges with key 2-tuple (edge.t_vert, edge.b_vert).
        """
        if not isinstance(edge, QEEdge):
            raise TypeError("edge must be QEEdge")
        
        if (edge.b_vert.index, edge.t_vert.index) in self._existing_edges:
            print("edge.b_vert.index: %d edge.t_vert.index: %d" %
                  (edge.b_vert.index, edge.t_vert.index))
            raise ValueError("edge already exists!")
        elif (edge.t_vert.index, edge.b_vert.index) in self._existing_edges:
            print("edge.b_vert.index: %d edge.t_vert.index: %d" %
                  (edge.b_vert.index, edge.t_vert.index))
            raise ValueError("edge already exists!")
        
        self._edges.append(edge)
        self._existing_edges[(edge.t_vert.index, edge.b_vert.index)] = edge

    def add_vertex(self, vert):
        """ Add a QEVertex to this QEMesh.
        
        Raises TypeError if vert isn't QEVertex.
        Raises ValueError if vert already exists.
        """
        if not isinstance(vert, QEVertex):
            raise TypeError("vert must be QEVertex!")

        if vert in self._vertices:
            raise ValueError("vert already exists in me!")
        
        self._vertices.append(vert)

    # def add_vertex_by_coords(self, coords):
    #     """ NOT IMPLEMENTED
    #     Not sure if this will be useful or not.
    #     """
    #     raise NotImplementedError("This method not yet implemented!")

    def is_valid(self):
        """ Performs a check (for testing) whether this mesh is valid or not.

        The following tests are performed:
          Every face's edges are attached to it,
          Every face has at least 3 edges and 3 vertices
          Every edge has one or more faces attached to it,
          Every edge's t and b verts are contained in this mesh
          Every edge is unique
          For any edge, following the edge.next or edge.prev chain
          will lead back to itself, and these edges will be contained exactly
          in the face you're iterating about.
        """

        # All edges in face are attached back to it
        for f in self._faces:
            for e in f.edges:
                if e.l_face is not f and e.r_face is not f:
                    raise ValueError("Edge does not contain face!")

        # Every face has at least 3 verts and 3 edges
        for f in self._faces:
            if len(f.edges) < 3 or len(f.verts) < 3:
                raise ValueError("Face should have at least 3 edges and 3 verts!")

        # Every edge has one or more faces attached to it
        for e in self._edges:
            if e.l_face is None and e.r_face is None:
                raise ValueError("Edge does not have any faces attached!")
            if e.l_face is None and e.r_face is not None:
                raise ValueError("Edge has right face but not left!")

        # Every edge's vertices are contained in mesh
        for e in self._edges:
            if (e.b_vert.index >= len(self._vertices) or
                e.t_vert.index >= len(self._vertices)):
                print("e.b_vert: %d, e.t_vert: %d" %
                      (e.b_vert.index, e.t_vert.index))
                raise ValueError("Edge's vertices aren't contained in me!")

        # Every edge has two faces (only for closed mesh)
        for e in self._edges:
            if e.l_face is None or e.r_face is None:
                print("This mesh is not closed!")
                break

        # Every edge is unique
        for e in self._edges:
            for e_other in self._edges:
                if e is e_other:
                    continue
                if (e.t_vert == e_other.t_vert and
                    e.b_vert == e_other.b_vert):
                    print("e.index: %d, e_other.index: %d" % (
                            e.index, e_other.index))
                    print("e.b_vert: %d, e.t_vert %d" % (
                            e.b_vert.index, e.t_vert.index))
                    print("other:")
                    print("e.b_vert: %d, e.t_vert %d" % (
                            e_other.b_vert.index, e_other.t_vert.index))
                    raise ValueError('Edge is duplicated!')
                if (e.t_vert == e_other.b_vert and
                    e.b_vert == e_other.t_vert):
                    print("e.index: %d, e_other.index: %d" % (
                            e.index, e_other.index))
                    print("e.b_vert: %d, e.t_vert %d" % (
                            e.b_vert.index, e.t_vert.index))
                    print("other:")
                    print("e.b_vert: %d, e.t_vert %d" % (
                            e_other.b_vert.index, e_other.t_vert.index))
                    raise ValueError('Edge is duplicated!')

        # Successive calls to edge.next_edge will wrap around a face
        for f in self._faces:
            e0 = f.edges[0]
            next_edge = e0.next_edge(f)
            counter = 1
            while next_edge is not e0:
                if next_edge not in f.edges:
                    print("currently at edge: %d in face %d" % (counter, f.index))
                    raise ValueError("Edge doesn't live in f.edges!")
                counter = counter+1
                next_edge = next_edge.next_edge(f)
            if counter is not len(f.edges):
                raise ValueError("Number of edges to go around face is not the " +
                                 "same as f.edges length! walk around: " +
                                 "%d, len(f.edges): %d" % counter, len(f.edges))

class QEVertex (object):
    """ A Quad-edge vertex.

    Holds pointers to parent (mesh), unique index in parent mesh (index),
    and a unique index in the parent mesh (index)
    """
    def __init__(self, parent_mesh, index):
        self.mesh = parent_mesh
        self.index = index
        self.pos = [0, 0, 0]

    def get_pos(self):
        return self.pos


class QEEdge (object):
    """ A Quad-edge edge.

    Holds pointers to parent (mesh), unique index in parent mesh (index),
    left and right faces (l_face, r_face),
    top and bottom vertices (t_vert, b_vert), and four other quad-edges
    (tl_vert, tr_vert, bl_vert, br_vert).

    The edge is considered to be directed (or point) from bottom to top.
    """
    def __init__(self, parent_mesh, index):
        self.mesh = parent_mesh
        self.index = index
        self.l_face = None
        self.r_face = None
        self.t_vert = None
        self.b_vert = None
        self.tl_edge = None
        self.tr_edge = None
        self.bl_edge = None
        self.br_edge = None

    def next_edge(self, face):
        """ Given a face, give the next QEEdge in a CCW direction around it.
        """
        if self.l_face == face:
            return self.tl_edge
        elif self.r_face == face:
            return self.br_edge
        else:
            return None

    def prev_edge(self, face):
        """ Given a face, give the previous QEEdge in a CCW direction around it.
        (that is, give the next edge in a CW direction)
        """
        if self.l_face == face:
            return self.bl_edge
        elif self.r_face == face:
            return self.tr_edge
        else:
            return None

class QEFace (object):
    """ A Quad-edge face.

    Holds pointers to parent (mesh), unique index in parent mesh (index),
    and a list of QEVerts and QEEdges which are adjacent to this face.

    """
    def __init__(self, parent_mesh, index):
        self.mesh = parent_mesh
        self.index = index
        self.verts = []
        self.edges = []
        self.min_coord = [0, 0, 0]
        self.max_coord = [0, 0, 0]

    def update_bounding_box(self):
        """ Update my bounding box defined by self.min_coord, self.max_coord.
        min_coord and max_coord must not be assigned to another list in order for
        AABB to work properly!!!
        """

        # self.min_coord[0] = sys.float_info.max
        # self.min_coord[1] = sys.float_info.max
        # self.min_coord[2] = sys.float_info.max
        # self.max_coord[0] = sys.float_info.min
        # self.max_coord[1] = sys.float_info.min
        # self.max_coord[2] = sys.float_info.min
        # for vert in self.verts:
        #     self.min_coord[0] = min(self.min_coord[0], vert.pos[0])
        #     self.min_coord[1] = min(self.min_coord[1], vert.pos[1])
        #     self.min_coord[2] = min(self.min_coord[2], vert.pos[2])
        #     self.max_coord[0] = max(self.max_coord[0], vert.pos[0])
        #     self.max_coord[1] = max(self.max_coord[1], vert.pos[1])
        #     self.max_coord[2] = max(self.max_coord[2], vert.pos[2])

        # TODO: this only works for tris. Above code is broken somehow
        self.min_coord[0] = min(self.verts[0].pos[0],
                                self.verts[1].pos[0],
                                self.verts[2].pos[0])
        self.min_coord[1] = min(self.verts[0].pos[1],
                                self.verts[1].pos[1],
                                self.verts[2].pos[1])
        self.min_coord[2] = min(self.verts[0].pos[2],
                                self.verts[1].pos[2],
                                self.verts[2].pos[2])
        self.max_coord[0] = max(self.verts[0].pos[0],
                                self.verts[1].pos[0],
                                self.verts[2].pos[0])
        self.max_coord[1] = max(self.verts[0].pos[1],
                                self.verts[1].pos[1],
                                self.verts[2].pos[1])
        self.max_coord[2] = max(self.verts[0].pos[2],
                                self.verts[1].pos[2],
                                self.verts[2].pos[2])
        
        # v0 = self.verts[0].get_pos()
        # v1 = self.verts[1].get_pos()
        # v2 = self.verts[2].get_pos()
        # self.min_coord[0] = min(v0[0], v1[0], v2[0])
        # self.min_coord[1] = min(v0[1], v1[1], v2[1])
        # self.min_coord[2] = min(v0[2], v1[2], v2[2])
        # self.max_coord[0] = max(v0[0], v1[0], v2[0])
        # self.max_coord[1] = max(v0[1], v1[1], v2[1])
        # self.max_coord[2] = max(v0[2], v1[2], v2[2])

def run_qemesh_test (fname):
    face_list = obj_reader.read_faces_to_list(fname)
    vert_list = obj_reader.read_verts_to_list(fname)
        
    print("Running test on %s..." % fname)
    print("  Constructing...")
    start = time()
    qem = QEMeshBuilder.construct(vert_list, face_list)
    print("  took %1.5f seconds" % (time()-start))
    
    print("  length of verts is %d" % len(vert_list))
    print("  length of faces is %d" % len(face_list))
    print("  length of edges is %d" % len(qem._edges)) #pylint: disable=W0212
    print("  Validating...")
    start = time()
    qem.is_valid()
    print("  took %1.5f seconds" % (time()-start))
    print("done!")

if __name__ == "__main__":
    run_qemesh_test('quad_edge_mesh/testMeshes/box.obj')
    run_qemesh_test('quad_edge_mesh/testMeshes/bbox.obj')

