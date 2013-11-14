quad_edge_mesh
==============

A Python class for a Quad Edge mesh, an AABB tree class for quickly computing collisions between QEMs.

quad_edge_mesh.py contains classes QEMesh, QEVertex, QEEdge, and QEFace. Includes test method run_qemesh_test which will read an obj file name, create a QEMesh, and check if the mesh is valid with QEMesh.is_valid(). Not an extremely thorough test.

aabb_collision contains an AABB collision tree class.

Includes obj_reader for reading *very* basic wavefront .obj files.

Includes two box meshes in a testMeshes subdirectory, one box with square faces, one box with tri faces.
