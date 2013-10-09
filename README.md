quad_edge_mesh
==============

A Python file for a Quad Edge mesh

quad_edge_mesh.py contains classes QEMesh, QEVertex, QEEdge, and QEFace. Includes test method run_qemesh_test which will read an obj file name, create a QEMesh, and check if the mesh is valid with QEMesh.is_valid(). Not an extremely thorough test.

Includes obj_reader for reading *very* basic wavefront .obj files.

Includes two box meshes in a testMeshes subdirectory, one box with square faces, one box with tri faces.
