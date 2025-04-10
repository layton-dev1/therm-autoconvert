import ezdxf
from config import settings

def read_dxf(file_path):
    #Open the dxf file in ezdxf and create a modelspace
    doc = ezdxf.readfile(file_path)
    msp = doc.modelspace()
    
    #Make sure drawing units are always in metric
    set_units(msp)
    polygons = extract_polygons(msp)
    
    return polygons

def set_units(msp):
    '''
    Finds the min and max x and y points to determine how large the model is. Anything above 15" is probably in Metric
    '''
    min_x, min_y = float('inf'), float('inf')
    max_x, max_y = float('-inf'), float('-inf')

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
    if(abs(max_x)-abs(min_x) > 15 or abs(max_y)-abs(min_y) > 15):
        settings.scale_factor = 1
    else:
        settings.scale_factor = 25.4
        print(f"Detecing imperial units, converting to metric")

def extract_polygons(msp):
    '''
    Loop through all entities in drawing and scale them appropriately
    '''
    polygons = []

    for entity in msp:
        #Only interested in polylines
        if entity.dxftype() == 'POLYLINE':
            try:
                entity.scale(settings.scale_factor, settings.scale_factor,1)
            except AttributeError:
                #Should probably just remove it
                print(f"Entity {entity.dxftype()} does not support scaling and was skipped.")

            #Might want to change this later to store as shapely polygons directly (instead of list of vertices)
            vertices = [(vertex.dxf.location[0], vertex.dxf.location[1]) for vertex in entity.vertices]
            layer_name = entity.dxf.layer
            polygons.append((vertices, layer_name))
    
    return polygons
