def read_faces_to_list(fname):
    """ Read faces from an obj file and return as a list
    of lists of indices.
    """
    f = open(fname, 'r')
    faces = []
    line = f.readline()
    while (line != ''):
        tokens = line.split()
        if not tokens:
            line = f.readline()
            continue
  
        if tokens[0] is '#' or line[0] is '#':
            line = f.readline()
            continue
        if tokens[0] is 'f':
            indices = []
            for index in tokens[1:]:
                indices.append (int(index))
  
            faces.append(indices)
      
        line = f.readline()

    f.close()
        
    return faces

def read_verts_to_list(fname):
    """ Read vertices from an obj file and return as a list
    of 3-lists of vertex coordinates
    """
    
    f = open(fname, 'r')
    vertices = []
    line = f.readline()
    while (line != ''):
        tokens = line.split()
        if not tokens:
            line = f.readline()
            continue
  
        if tokens[0] is '#' or line[0] is '#':
            line = f.readline()
            continue
        if tokens[0] is 'v':
            coords = []
            for coord in tokens[1:]:
                coords.append (float(coord))
  
            vertices.append(coords)
  
        line = f.readline()

    f.close()
    
    return vertices
