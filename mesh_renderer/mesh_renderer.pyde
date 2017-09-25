rotate_flag = True    # automatic rotation of model?
time = 0   # keep track of passing time, for automatic rotation

geometry_table = [] # polyhedral data structure
vertex_table = []
opposite_table = []
rand_face_colors_list = []
vertex_normal_list = []

num_vertices = 0
num_faces = 0

per_vertex_shading = 0
rand_color_faces = 0
color_faces_white = 0
dual_iteration = 0


# initalize
def setup():
    size (600, 600, OPENGL)
    noStroke()
    

# draw the current mesh
def draw():
    global time
    
    background(0)    # clear screen to black

    perspective (PI*0.333, 1.0, 0.01, 1000.0)
    camera (0, 0, 5, 0, 0, 0, 0, 1, 0)    # place the camera in the scene
    scale (1, -1, 1)    # change to right-handed coordinate system
    
    # create an ambient light source
    ambientLight (102, 102, 102)
  
    # create two directional light sources
    lightSpecular (204, 204, 204)
    directionalLight (102, 102, 102, -0.7, -0.7, -1)
    directionalLight (152, 152, 152, 0, 0, -1)
    
    pushMatrix()
    if (color_faces_white):
        fill(255, 255, 255)
    elif (rand_color_faces):
        pass
    else:
        fill (50, 50, 200)            # set polygon color
    ambient (200, 200, 200)
    specular (0, 0, 0)            # no specular highlights
    shininess (1.0)
  
    rotate (time, 1.0, 0.0, 0.0)  
       
    draw_mesh()
    
    popMatrix()
    
    # step forward in time (for object rotation)
    if rotate_flag:
        time += 0.02
        
def reset():
    global geometry_table, vertex_table, opposite_table, rand_face_colors_list
    geometry_table = []
    vertex_table = []
    rand_face_colors_list = []
    opposite_table = []
    
def dual():
    global vertex_table, geometry_table
    centroid_table = []
    for i in range(0,len(vertex_table),3):
        vert1 = geometry_table[vertex_table[i]]
        vert2 = geometry_table[vertex_table[i+1]]
        vert3 = geometry_table[vertex_table[i+2]]
        vertex_xyz = triangle_centroid(vert1, vert2, vert3)
        centroid_table.append(myVertex(vertex_xyz[0], vertex_xyz[1], vertex_xyz[2], i/3))
    
    new_vertex_table = []
    visited = []
    for i in range (len(vertex_table)):
        visited.append(False)
        
    for i in range(len(vertex_table)):
        si = i
        neighbor_table = []
        while(True):
            
            if (visited[si] == True):
                break

            visited[si] = True
            si = swing(si)

            vert1 = geometry_table[vertex_table[(si/3) * 3]]
            vert2 = geometry_table[vertex_table[((si/3) * 3) + 1]]
            vert3 = geometry_table[vertex_table[((si/3) * 3) + 2]]
            vec_cent = triangle_centroid(vert1, vert2, vert3)
            vert_cent = myVertex(vec_cent[0], vec_cent[1], vec_cent[2], 0)
            for j in range(len(centroid_table)):
                if (vertices_equal(centroid_table[j], vert_cent)):
                    neighbor_table.append(j)
            
        if (len(neighbor_table) > 3):
            poly_centroid = polygon_centroid(neighbor_table, centroid_table)
            centroid_table.append(myVertex(poly_centroid[0], poly_centroid[1], poly_centroid[2], 0))
            new_neighbor_table = []
            for k in range(0,len(neighbor_table)):
                new_neighbor_table.append(len(centroid_table)-1)
                new_neighbor_table.append(neighbor_table[k])
                new_neighbor_table.append(neighbor_table[(k+1)%len(neighbor_table)])
                
            neighbor_table = new_neighbor_table
                    
        new_vertex_table = new_vertex_table + neighbor_table
            
    geometry_table = centroid_table
    vertex_table = new_vertex_table
    calculate_opposite_table()
    
def polygon_centroid(vertList, centroid_table):
    global geometry_table
    sum = [0, 0, 0]

    for i in range(len(vertList)):
        v = centroid_table[vertList[i]]
        vec = [v.x, v.y, v.z]
        sum = vector_add(sum, vec)
    sum = [sum[0]/len(vertList), sum[1]/len(vertList), sum[2]/len(vertList)]
    return sum
    
def triangle_centroid(vert1, vert2, vert3):
    v1 = [vert1.x, vert1.y, vert1.z]
    v2 = [vert2.x, vert2.y, vert2.z]
    v3 = [vert3.x, vert3.y, vert3.z]
    sum = vector_add(v1, v2)
    sum = vector_add(sum, v3)
    sum = [sum[0]/3.0, sum[1]/3.0, sum[2]/3.0]
    return sum
    
def vertices_equal(vert1, vert2):
    if (vert1.x != vert2.x):
        return False
    if (vert1.y != vert2.y):
        return False
    if (vert1.z != vert2.z):
        return False
    return True 
    
    
def calculate_rand_face_colors():
    global rand_face_colors_list, vertex_table
    for i in range(0, len(vertex_table), 3):
        rand_face_colors_list.append([random(0,255), random(0,255), random(0,255)])
    

# process key presses
def keyPressed():
    global rotate_flag, per_vertex_shading, rand_color_faces, color_faces_white, dual_iteration, vertex_normal_list
    if key == ' ':
        rotate_flag = not rotate_flag
    elif key == '1':
        read_mesh ('tetra.ply')
        dual_iteration = 0
    elif key == '2':
        read_mesh ('octa.ply')
        dual_iteration = 0
    elif key == '3':
        read_mesh ('icos.ply')
        dual_iteration = 0
    elif key == '4':
        read_mesh ('star.ply')
        dual_iteration = 0
    elif key == '5':
        read_mesh ('torus.ply')
        dual_iteration = 0
    elif key == 'n':
        # toggle per-vertex shading
        if (per_vertex_shading == 1):
            per_vertex_shading = 0
        else:
            per_vertex_shading = 1
    elif key == 'r':
        # randomly color faces
        rand_color_faces = 1
        color_faces_white = 0
    elif key == 'w':
        # color faces white
        color_faces_white = 1
        rand_color_faces = 0
    elif key == 'd':
        # calculate the dual mesh
        dual()
    elif key == 'q':
        exit()
        
        
    calculate_rand_face_colors()

# read in a mesh file
def read_mesh(filename):
    reset()
    global geometry_table, vertex_table, num_vertices, num_faces
    
    fname = "data/" + filename
    # read in the lines of a file
    with open(fname) as f:
        lines = f.readlines()
        
    # determine number of vertices (on first line)
    words = lines[0].split()
    num_vertices = int(words[1])
    print "number of vertices =", num_vertices

    # determine number of faces (on first second)
    words = lines[1].split()
    num_faces = int(words[1])
    print "number of faces =", num_faces

    # read in the vertices
    for i in range(num_vertices):
        words = lines[i+2].split()
        x = float(words[0])
        y = float(words[1])
        z = float(words[2])
        print "vertex = ", x, y, z
        v = myVertex(x, y, z, i)
        geometry_table.append(v)
        
    
    # read in the faces
    for i in range(num_faces):
        j = i + num_vertices + 2
        words = lines[j].split()
        nverts = int(words[0])
        if nverts != 3:
            print "error: this face is not a triangle"
            exit()
        
        index1 = int(words[1])
        index2 = int(words[2])
        index3 = int(words[3])
        print "face =", index1, index2, index3
        
        vertex_table.append(index1)
        vertex_table.append(index2)
        vertex_table.append(index3)
        
    calculate_opposite_table()
    
    
    
        
def calculate_opposite_table():
    global opposite_table, vertex_table
    opposite_table = [None for i in range(len(vertex_table))]
    for i in range(0, len(vertex_table)):
        for j in range(0, len(vertex_table)):
            a_n_v = vertex_table[next(i)]
            b_p_v = vertex_table[prev(j)]
            a_p_v = vertex_table[prev(i)]
            b_n_v = vertex_table[next(j)]
            
            if (a_n_v == b_p_v and a_p_v == b_n_v):
                opposite_table[i] = j
                opposite_table[j] = i
                
            
def tri(i): 
    t = i/3
    return t

def next(i): #gets index of next vertex    
    n = 3*tri(i) + (i+1)%3
    return n

def prev(i): # gets index of prev vertex
    p = next(next(i))
    return p

def swing(i):
    global opposite_table
    c_n = next(i)
    c_n_o = opposite_table[c_n]
    c_n_o_n = next(c_n_o)
    return c_n_o_n

def calculate_vertex_normal(i):
    global vertex_table, geometry_table
    corner = vertex_table[i]
    num_sides = 0
    for j in range(len(vertex_table)):
        if (vertex_table[j] == corner):
            num_sides = num_sides + 1
            
    neighbors = []
    for j in range(num_sides):
        neighbors.append(next(swing(i)))
        i = swing(i)
        
    sum = [0,0,0]
    for j in range(len(neighbors)):
        face_num = neighbors[j]/3
        starting_vertex_num = face_num * 3
        vert1 = geometry_table[vertex_table[starting_vertex_num]]
        vert2 = geometry_table[vertex_table[starting_vertex_num+1]]
        vert3 = geometry_table[vertex_table[starting_vertex_num+2]]
        
        face_normal = face_normal_calculation(vert1, vert2, vert3)
        sum = vector_add(sum, face_normal)
    
    sum = [sum[0]/len(neighbors), sum[1]/len(neighbors), sum[2]/len(neighbors)]
    
    return sum
            

def vertex_sub(v1, v2):
    return [ v1.x-v2.x , v1.y-v2.y , v1.z-v2.z ]

def vector_add(v1, v2):
    return [ v1[0]+v2[0] , v1[1]+v2[1] , v1[2]+v2[2] ]

def vector_sub(v1, v2):
    return [v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]]

def vector_cross_product(v2, v1):
    a = v1[0]
    b = v1[1]
    c = v1[2]
    d = v2[0]
    e = v2[1]
    f = v2[2]
    
    return [ b*f - c*e , c*d - a*f , a*e - b*d ]

def normalize_vector(v):
    d = sqrt( sq(v[0]) + sq(v[1]) + sq(v[2]) )
    return [ v[0]/d , v[1]/d , v[2]/d ]
    
            
class myVertex():
    
    def __init__(self, x, y, z, i):
        self.x = x
        self.y = y
        self.z = z
        self.c = i/3
        
    def getX(self):
        return self.x
    
    def getX(self):
        return self.y
    
    def getX(self):
        return self.z
    

    
    
class myFace():
    
    def __init__(self, i1, i2, i3):
        self.i1 = i1
        self.i2 = i2
        self.i3 = i3
        
        
def draw_mesh():
    global vertex_table, geometry_table, per_vertex_shading, rand_color_faces, color_faces_white, dual_iteration, opposite_table, rand_face_colors_list
    
    #print(model_data)
    if (per_vertex_shading):
       for i in range(0, len(vertex_table), 3):
            if (rand_color_faces):
               c = rand_face_colors_list[i/3]
               fill(c[0], c[1], c[2])
               
            beginShape(TRIANGLES)
            
            v1 = geometry_table[vertex_table[i]]
            n1 = calculate_vertex_normal(i)
            v2 = geometry_table[vertex_table[i+1]]
            n2 = calculate_vertex_normal(i+1)
            v3 = geometry_table[vertex_table[i+2]]
            n3 = calculate_vertex_normal(i+2)
            
            
            normal(n1[0],n1[1],n1[2])
            vertex(v1.x, v1.y, v1.z)
            normal(n2[0],n2[1],n2[2])
            vertex(v2.x, v2.y, v2.z)
            normal(n3[0],n3[1],n3[2])
            vertex(v3.x, v3.y, v3.z)
            
            endShape()
            
    else:
        for i in range(0, len(vertex_table), 3):
            
            
            
            if (rand_color_faces):
                c = rand_face_colors_list[i/3]
                fill(c[0], c[1], c[2])
            
            beginShape(TRIANGLES)
            
            v1 = geometry_table[vertex_table[i]]
            v2 = geometry_table[vertex_table[i+1]]
            v3 = geometry_table[vertex_table[i+2]]
            
            
            
            face_normal = face_normal_calculation(v1, v2, v3)
            
            normal(face_normal[0], face_normal[1], face_normal[2])
            vertex(v1.x, v1.y, v1.z)
            vertex(v2.x, v2.y, v2.z)
            vertex(v3.x, v3.y, v3.z)
            
            endShape()
            
                
def face_normal_calculation(vert1, vert2, vert3):
    
    v1 = [vert1.x, vert1.y, vert1.z]
    v2 = [vert2.x, vert2.y, vert2.z]
    v3 = [vert3.x, vert3.y, vert3.z]
    
    v1_v2 = vector_sub(v1, v2)
    v1_v3 = vector_sub(v1, v3)
    
    vcp = vector_cross_product(v1_v2, v1_v3)
    return normalize_vector(vcp)
    
    
    
    
    
    
    
    
    
    
    
    
        