import ezdxf

def read_dxf(file_path):
    doc = ezdxf.readfile(file_path)
    return extract_polygons(doc)

def extract_polygons(doc):
    polygons = []
    min_x, min_y = float('inf'), float('inf')
    max_x, max_y = float('-inf'), float('-inf')

    # Iterate over all entities in the modelspace
    msp = doc.modelspace()

    #Determine drawings units
    for entity in msp:
        if entity.dxftype() == 'POLYLINE':
            for vertex in entity.vertices:
                x = vertex.dxf.location[0]
                y = vertex.dxf.location[1]
                min_x = min(min_x, x)
                min_y = min(min_y, y)
                max_x = max(max_x, x)
                max_y = max(max_y, y)
    if(abs(max_x)-abs(min_x) > 15 or abs(max_y)-abs(min_y)):
        SCALE_FACTOR = 25.4
    else:
        SCALE_FACTOR = 1


    for entity in msp:
        if entity.dxftype() == 'POLYLINE':
            try:
                entity.scale(SCALE_FACTOR,SCALE_FACTOR,1)
            except AttributeError:
                print(f"Entity {entity.dxftype()} does not support scaling and was skipped.")

            vertices = [(vertex.dxf.location[0], vertex.dxf.location[1]) for vertex in entity.vertices]
            layer_name = entity.dxf.layer
            polygons.append((vertices, layer_name))
    
    return polygons
