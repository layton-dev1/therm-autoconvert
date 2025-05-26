import xml.etree.ElementTree as ET
from config import config

def create_spacer(materials_element, polygons_element, origin, desired_size=12.7):
    file_path = config.data["spacer_file_path"] 
    tree = ET.parse(file_path)
    root = tree.getroot()
    
    xmlns={"therm": "http://windows.lbl.gov"}
    
    polygons = root.findall(".//therm:Polygon", xmlns)
    
    min_x = float('inf')
    min_y = float('inf')
    max_x = float('-inf')
    bottom_left_corner = (0, 0)
    
    for poly in polygons:
        points = poly.findall(".//therm:Point", xmlns)
        for point in points:
            x = float(point.attrib["x"])
            y = float(point.attrib["y"])
            if (y < min_y) or (y == min_y and x < min_x):
                min_x = x
                min_y = y
                bottom_left_corner = [x, y]
            elif x > max_x:
                max_x = x
    
    print(f"min x: {min_x}, min y: {min_y}")
    print(f"blc: {bottom_left_corner[0]}, {bottom_left_corner[1]}")
    print(f"origin x: {origin[0]} origin y: {origin[1]}")
    spacer_size = max_x-min_x
    scale = spacer_size/12.7
    desired_size = desired_size*scale
    to_add = desired_size-spacer_size
    halfway_point = min_x + spacer_size/2
    
    diff_from_origin = (float(bottom_left_corner[0])-origin[0], float(bottom_left_corner[1])-origin[1])
    
    print(f"Diff: {diff_from_origin}")
    for poly in polygons:
        points = poly.findall(".//therm:Point", xmlns)
        for point in points:
            x_value = float(point.attrib["x"])
            if(x_value > halfway_point):
                point.attrib["x"] = str(x_value + to_add)
                
    for poly in polygons:
        points = poly.findall(".//therm:Point", xmlns)
        polygon = []
        for point in points:
            point.attrib["x"] = str(float(point.attrib["x"])-diff_from_origin[0])
            point.attrib["y"] = str(float(point.attrib["y"])-diff_from_origin[1])
            polygon.append((point.attrib["x"], point.attrib["y"]))
        material = get_material_properties_by_name(root, poly.attrib["Material"], xmlns)
        create_spacer_element(materials_element, polygons_element, polygon, material)
    
    # min_x = float('inf')
    # min_y = float('inf')
    # for poly in polygons:
    #     points = poly.findall(".//therm:Point", xmlns)
    #     for point in points:
    #         x_value = float(point.attrib["x"])
    #         y_value = float(point.attrib["y"])
    #         if x_value < min_x:
    #             min_x = x_value
    #             if y_value < min_y:
    #                 bottom_left_corner = point
    #                 min_y = y_value
    #         elif x_value > max_x:
    #             max_x = x_value
    
    # print(f"new min x: {min_x}, new min y: {min_y}")
    #ET.register_namespace('', 'http://windows.lbl.gov')
    #tree.write("lib/data/Test_spacer.thmx", encoding="utf-8", xml_declaration=True)
    
def create_spacer_element(materials_element, polygons_element, polygon, material):
    if material == None:
        material = {
                'Name': "Spacer",
                'Conductivity': str(1),
                'Emissivity': "0.900000",
                'Color': "0xC0C0C0",
                'Type': "0"
            }
        
    last_id = materials_element[-1].attrib.get("Index")
    # Create the Material element
    material_element = ET.SubElement(
        materials_element,
        "Material",
        Name=material['Name'],
        Index=str(int(last_id)+1),
        Type=material['Type'],
        Conductivity=material['Conductivity'],
        Tir="0.000000",
        EmissivityFront=material['Emissivity'],
        EmissivityBack=material['Emissivity'],
        RGBColor=material['Color']
    )
    
    for side in ["Front", "Back"]:
                for range_type in ["Visible", "Solar"]:
                    for specularity in ["Direct", "Diffuse"]:
                        ET.SubElement(
                            material_element,
                            "Property",
                            Side=side,
                            Range=range_type,
                            Specularity=specularity,
                            T="0.000000",
                            R="0.000000"
                        )
    
    last_id = polygons_element[-1].attrib.get("ID")
    polygon_element = ET.SubElement(
        polygons_element,
        "Polygon",
        ID=str(int(last_id)+1),
        Material=material['Name'],
        NSides=str(len(polygon)),
        Type="1",
        units="mm"
    )
    
    # Add each point to the polygon
    for j, (x, y) in enumerate(polygon):
        ET.SubElement(
            polygon_element,
            "Point",
            index=str(j),
            x=x,#f"{x:.6f}",
            y=y#f"{y:.6f}"
        )
        
def get_material_properties_by_name(root, name, xmlns):
    """
    Search for a material by name in the XML tree and return its properties.

    Parameters:
        root (Element): The root of the ElementTree.
        name (str): The material name to search for.

    Returns:
        dict: A dictionary of the material's properties if found, else None.
    """
    materials_node = root.find(".//therm:Materials", xmlns)
    if materials_node is None:
        print("Found no materials node")
        return None
    
    for material in materials_node.findall(".//therm:Material", xmlns):
        if material.attrib.get("Name") == name:
            return {
                'Name': name,
                'Conductivity': material.attrib.get("Conductivity", "0.0"),
                'Emissivity': material.attrib.get("EmissivityFront", "0.9"),
                'Color': material.attrib.get("RGBColor", "0xC0C0C0"),
                'Type': material.attrib.get("Type", "0")
            }
    
    return None  # Material not found
