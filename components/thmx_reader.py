import xml.etree.ElementTree as ET
from config import config

def create_spacer(materials_element, polygons_element, origin, desired_size=12.7):
    file_path = config.data["spacer_file_path"]
    tree = ET.parse(file_path)
    root = tree.getroot()
    
    xmlns = {"therm": "http://windows.lbl.gov"}
    polygons = root.findall(".//therm:Polygon", xmlns)

    # Determine original bounding box
    min_x = float('inf')
    min_y = float('inf')
    max_x = float('-inf')

    for poly in polygons:
        for point in poly.findall(".//therm:Point", xmlns):
            x = float(point.attrib["x"])
            y = float(point.attrib["y"])
            min_x = min(min_x, x)
            min_y = min(min_y, y)
            max_x = max(max_x, x)

    original_width = max_x - min_x

    # Compute X scale only
    scale_x = desired_size / original_width
    x0, y0 = origin
    translate_x = x0 - min_x * scale_x
    translate_y = y0 - min_y  # no y-scale

    # Apply X-scale and full translation to each point
    for poly in polygons:
        polygon = []
        for point in poly.findall(".//therm:Point", xmlns):
            x = float(point.attrib["x"])
            y = float(point.attrib["y"])
            new_x = round(x * scale_x + translate_x, 6)
            new_y = round(y + translate_y, 6)  # Y is just translated
            point.attrib["x"] = str(new_x)
            point.attrib["y"] = str(new_y)
            polygon.append((str(new_x), str(new_y)))

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
