import xml.etree.ElementTree as ET
from datetime import datetime
from utils import utils
import pyodbc
from config import config

def create_xml(polygons, materials, output_file):
    # Create the root element
    root = ET.Element("THERM-XML", xmlns="http://windows.lbl.gov")
    
    # Add metadata elements
    ET.SubElement(root, "ThermVersion").text = "Version 7.8.74.0"
    ET.SubElement(root, "FileVersion").text = "1"
    ET.SubElement(root, "SaveDate").text = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    ET.SubElement(root, "Title")
    ET.SubElement(root, "CreatedBy")
    ET.SubElement(root, "Company")
    ET.SubElement(root, "Client")
    ET.SubElement(root, "CrossSectionType").text = "Jamb"
    ET.SubElement(root, "Notes")
    ET.SubElement(root, "Units").text = "SI"
    
    # Add MeshControl element
    mesh_control = ET.SubElement(root, "MeshControl", MeshLevel="9", ErrorCheckFlag="1", ErrorLimit="10.000000", MaxIterations="5", CMAflag="0")
    
    # Add RadianceModeBC element
    radiance_mode_bc = ET.SubElement(root, "RadianceModeBC")
    ET.SubElement(radiance_mode_bc, "RadianceModeInBCName").text = "ShadeInE"
    ET.SubElement(radiance_mode_bc, "RadianceModeInBCTag").text = "ShadeInETag"
    ET.SubElement(radiance_mode_bc, "RadianceModeOutBCName").text = "ShadeOutE"
    ET.SubElement(radiance_mode_bc, "RadianceModeOutBCTag").text = "ShadeOutETag"

    materials_element = ET.SubElement(root, "Materials")
    boundaryconditions_element = ET.SubElement(root, "BoundaryConditions")
    polygons_element = ET.SubElement(root, "Polygons")
    boundaries_element = ET.SubElement(root, "Boundaries")
    
    root, bc_id = polygon_xml_section(root, polygons, materials, materials_element, polygons_element)
    
    mdb_file_location = r"C:\Users\tyler.henderson\Documents\DXFtoTHERM\Glazing\DB.mdb"
    mdb_glzsys_name = "G1"

    #root = import_glass_from_mdb_to_xml(root, polygons_element, mdb_file_location, mdb_glzsys_name, materials_element)
    
    root = boundary_xml_section(root, materials_element, polygons_element, boundaryconditions_element, boundaries_element, bc_id + 1)

    # Write the XML to a file
    tree = ET.ElementTree(root)
    tree.write(output_file, encoding="utf-8", xml_declaration=True)

def polygon_xml_section(root, polygons, materials, materials_element, polygons_element):
    unique_materials = {}
    material_index = 1
    bc_id = 0
    cavity_count = 2

    for i, (polygon, layer_name) in enumerate(polygons):
        # Find the closest matching material for the layer name
        if "frame" in layer_name.lower().strip() and "cavity" in layer_name.lower().strip():
            #Define the material
            material = {
                'Name': f"Frame Cavity NFRC 100" + "_Cavity_" + str(cavity_count),
                'Conductivity': "-1",
                'Emissivity': "0.900000",
                'Color': "0x00FF00",
                'Type': "1",
                'CavityModel': "4"
            }

            #Create material Element
            material_element = ET.SubElement(
                    materials_element,
                    "Material",
                    Name=material['Name'],
                    Index=str(material_index),
                    Type=material['Type'],
                    Conductivity=material['Conductivity'],
                    Tir="1.000000",
                    EmissivityFront=material['Emissivity'],
                    EmissivityBack=material['Emissivity'],
                    RGBColor=material['Color'],
                    CavityModel=material['CavityModel']
                )
            cavity_count = cavity_count + 1
        else:
            #Define the material
            material = utils.find_closest_material(layer_name, materials)

            #If none found add in placeholder
            if not material:
                    print(f"Warning: No matching material found for layer '{layer_name}'. Using default values.")
                    material = {
                        'Name': f"{layer_name}, Appendix A",
                        'Conductivity': "0.170000",
                        'Emissivity': "0.900000",
                        'Color': "0x000000",
                        'Type': "0"
                    }
            
            # Add material to the Materials section if not already added
            if layer_name not in unique_materials:
                material_index = len(unique_materials) + 1
                
                # Create the Material element
                material_element = ET.SubElement(
                    materials_element,
                    "Material",
                    Name=material['Name'],
                    Index=str(material_index),
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
                
                # Add Property elements for Front and Back sides
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
                
                unique_materials[layer_name] = material_index
        
        # Add the polygon
        polygon_element = ET.SubElement(
            polygons_element,
            "Polygon",
            ID=str(i + 1),
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
                x=f"{x:.6f}",
                y=f"{y:.6f}"
            )
        bc_id = i + 1

    return root, bc_id

def boundary_xml_section(root, materials_element, polygons_element, boundaryconditions_element, boundaries_element, bc_id):
    boundary_conditions = [
        {
            "Name": "Adiabatic",
            "Type": "0",
            "H": "0.000000",
            "HeatFLux": "0.000000",
            "Temperature": "0.000000",
            "RGBColor": "0x000000",
            "Tr": "0.000000",
            "Hr": "0.000000",
            "Ei": "1.000000",
            "Viewfactor": "1.000000",
            "RadiationModel": "0",
            "ConvectionFlag": "0",
            "FluxFlag": "1",
            "RadiationFlag": "0",
            "ConstantTemperatureFlag": "0",
            "EmisModifier": "1.000000"
        },
        {
            "Name": "Interior Thermally Broken Frame (convection only)",
            "Type": "1",
            "H": "3.000000",
            "HeatFLux": "0.000000",
            "Temperature": "21.000000",
            "RGBColor": "0xFF8080",
            "Tr": "21.000000",
            "Hr": "-431602080.000000",
            "Ei": "1.000000",
            "Viewfactor": "1.000000",
            "RadiationModel": "3",
            "ConvectionFlag": "1",
            "FluxFlag": "0",
            "RadiationFlag": "1",
            "ConstantTemperatureFlag": "0",
            "EmisModifier": "1.000000"
        },
        {
            "Name": "NFRC 100-2010 Exterior",
            "Type": "1",
            "H": "26.000000",
            "HeatFLux": "0.000000",
            "Temperature": "-18.000000",
            "RGBColor": "0x0080C0",
            "Tr": "-18.000000",
            "Hr": "-431602080.000000",
            "Ei": "1.000000",
            "Viewfactor": "1.000000",
            "RadiationModel": "1",
            "ConvectionFlag": "1",
            "FluxFlag": "0",
            "RadiationFlag": "1",
            "ConstantTemperatureFlag": "0",
            "EmisModifier": "1.000000"
        },
        {
            "Name": "Frame Cavity Surface",
            "Type": "7",
            "H": "0.000000",
            "HeatFLux": "0.000000",
            "Temperature": "0.000000",
            "RGBColor": "0xFF0000",
            "Tr": "0.000000",
            "Hr": "0.000000",
            "Ei": "1.000000",
            "Viewfactor": "1.000000",
            "RadiationModel": "0",
            "ConvectionFlag": "0",
            "FluxFlag": "0",
            "RadiationFlag": "0",
            "ConstantTemperatureFlag": "0",
            "EmisModifier": "1.000000"
        }
    ]
    
    # Add each boundary condition to the XML
    for condition in boundary_conditions:
        boundary_condition_element = ET.SubElement(
            boundaryconditions_element,
            "BoundaryCondition",
            Name=condition["Name"],
            Type=condition["Type"],
            H=condition["H"],
            HeatFLux=condition["HeatFLux"],
            Temperature=condition["Temperature"],
            RGBColor=condition["RGBColor"],
            Tr=condition["Tr"],
            Hr=condition["Hr"],
            Ei=condition["Ei"],
            Viewfactor=condition["Viewfactor"],
            RadiationModel=condition["RadiationModel"],
            ConvectionFlag=condition["ConvectionFlag"],
            FluxFlag=condition["FluxFlag"],
            RadiationFlag=condition["RadiationFlag"],
            ConstantTemperatureFlag=condition["ConstantTemperatureFlag"],
            EmisModifier=condition["EmisModifier"]
        )

    materials = []

    for material in materials_element:
        materials.append({
            "Name": material.attrib.get("Name"),
            "Emissivity": material.attrib.get("EmissivityFront")
            })

    # Add BCPolygon for each polygon
    for polygon_element in polygons_element:
        # Extract the ID and MaterialName of each polygon from the XML
        polygon_id = int(polygon_element.attrib.get('ID'))
        material_name = polygon_element.attrib.get('Material')
        emissivity = ""
        for material in materials_element:
            if material.attrib.get("Name") == material_name:
                emissivity = material.attrib.get("EmissivityFront")

        # Extract the points of the polygon
        polygon_points = []
        for point_element in polygon_element.findall("Point"):
            x = float(point_element.attrib["x"])
            y = float(point_element.attrib["y"])
            polygon_points.append((x, y))

        if not utils.is_counterclockwise(polygon_points):
            polygon_points.reverse()
        
        touching_edges = utils.find_touching_edges(polygons_element, polygon_id, polygon_points)

        #This is a debug function to see which edges are being marked as touching
        #print(bc_id)
        #print(polygon_id)
        #plot_polygon_debug(polygon_points, touching_edges)

        bc_id = create_solid_boundary(polygon_points, touching_edges, boundaries_element, polygon_id, material_name, bc_id, emissivity)

        #Frame cavities need an extra boundary condition called Frame Cavity Surface 
        if "Frame Cavity NFRC 100" in material_name:
            polygon_points.reverse()
        
            touching_edges = utils.find_touching_edges(polygons_element, polygon_id, polygon_points)
            bc_id = create_frame_cavity_boundary(polygon_points, materials_element, touching_edges, boundaries_element, polygon_id, bc_id, emissivity)

    return root



def create_frame_cavity_boundary(polygon_points, materials_element, touching_edges, boundaries_element, polygon_id, bc_id, emissivity):
    for side in range(len(polygon_points)):  # Loop through each side of the polygon
        p1 = polygon_points[side]
        p2 = polygon_points[(side + 1) % len(polygon_points)]
        edge = (p1, p2)

        bc_polygon_name = "Frame Cavity Surface"
        material_name = ""

        if edge in touching_edges:  # If the edge is touching another polygon, the material name should be set to that adjacent polygon's material
            material_name = touching_edges[edge]
            for material in materials_element:
                if material.attrib.get("Name") == material_name:
                    emissivity = material.attrib.get("EmissivityFront")
        else:
            emissivity = "0.9"
        
        create_bc_polygon(boundaries_element, bc_id, bc_polygon_name, material_name, polygon_id, p1, p2, emissivity, "1")
                
        bc_id = bc_id + 1

    return bc_id

def create_solid_boundary(polygon_points, touching_edges, boundaries_element, polygon_id, material_name, bc_id, emissivity):
    for side in range(len(polygon_points)):  # Loop through each side of the polygon
        p1 = polygon_points[side]
        p2 = polygon_points[(side + 1) % len(polygon_points)]
        edge = (p1, p2)

        if not edge in touching_edges: #Boundaries for solid materials should only appear when edge it not touching another polygon
            # Check the boundary's orientation to determine if it's exterior or interior
            if p1[0] > p2[0]:  # If the y-coordinate of p1 is less than p2, it's facing the exterior (upward)
                bc_polygon_name = "NFRC 100-2010 Exterior"
                ufactortag = "SHGC Exterior"
            elif p1[0] < p2[0]:  # If the y-coordinate of p1 is greater than p2, it's facing the interior (downward)
                bc_polygon_name = "Interior Thermally Broken Frame (convection only)"
                ufactortag = "Frame"
            elif p1[1] == p2[1]:
                bc_polygon_name = "Adiabatic"
                ufactortag = ""
            else:  # Handle vertical boundaries (when x-coordinates are the same)
                if p1[1] > p2[1]:  # If the y-coordinate of p1 is less than p2, it's facing the exterior (upward)
                    bc_polygon_name = "NFRC 100-2010 Exterior"
                    ufactortag = "SHGC Exterior"
                elif p1[1] < p2[1]:  # If the y-coordinate of p1 is greater than p2, it's facing the interior (downward)
                    bc_polygon_name = "Interior Thermally Broken Frame (convection only)"
                    ufactortag = "Frame"
                else:  # Handle the case where the points are exactly the same
                    bc_polygon_name = "Adiabatic"
                    ufactortag = ""
            #bc_polygon_name = "Adiabatic"
            create_bc_polygon(boundaries_element, bc_id, bc_polygon_name, material_name, polygon_id, p1, p2, emissivity, ufactortag=ufactortag)
            bc_id = bc_id + 1

    return bc_id

def create_bc_polygon(boundaries_element, bc_id, bc_polygon_name, material_name, polygon_id, p1, p2, emissivity, enclosure="0", ufactortag=""):
    bc_polygon = ET.SubElement(
        boundaries_element,
        "BCPolygon",
        ID=str(bc_id),
        BC=bc_polygon_name,
        units="mm",
        MaterialName=material_name,
        PolygonID=str(polygon_id),
        EnclosureID=enclosure,
        UFactorTag=ufactortag,
        Emissivity=emissivity,
        MaterialSide="Front",  # Could be dynamic depending on side
        IlluminatedSurface="FALSE"
    )

    # Add points for the BCPolygon
    for idx, (x, y) in enumerate([p1, p2]):
        ET.SubElement(
            bc_polygon,
            "Point",
            index=str(idx),
            x=f"{x:.6f}",
            y=f"{y:.6f}"
        )

def import_glass_from_mdb_to_xml(root, polygons_element, mdb_file_location, mdb_glzsys_name, materials_element):
    # Step 1: Open the MDB file and access the "GlzSys" table to find the "GlzSys" ID.
    conn_str = r'DRIVER={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=' + mdb_file_location
    conn = pyodbc.connect(conn_str)
    cursor = conn.cursor()

    # Search for the rectangle on the Glass1 layer
    glass_layer = "Glass1"
    bottom_left_corner = ""
    for polygon_element in polygons_element:
        material_name = polygon_element.attrib.get('Material')

        if glass_layer.strip().lower() in material_name.strip().lower():
            bottom_left_corner = (float(polygon_element.findall("Point")[0].attrib["x"]), 
                                  float(polygon_element.findall("Point")[0].attrib["y"]))

            for point_element in polygon_element.findall("Point"):
                x = float(point_element.attrib["x"])
                y = float(point_element.attrib["y"])
                
                if x < bottom_left_corner[0] and y < bottom_left_corner[1]:
                    bottom_left_corner = (x, y)
                        
            polygons_element.remove(polygon_element)

    # Query to get the GlzSys ID by Name
    cursor.execute("SELECT ID FROM GlzSys WHERE Name = ?", (mdb_glzsys_name,))
    glzsys_id = cursor.fetchone()

    if not glzsys_id:
        print(f"No GlzSys entry found with Name = {mdb_glzsys_name}")
        return
    glzsys_id = glzsys_id[0]

    # Step 2: Retrieve child elements from the GlassList table
    cursor.execute("SELECT ID FROM GlassList WHERE ParentID = ?", (glzsys_id,))
    glass_childs_ids = cursor.fetchall()
    
    cursor.execute("SELECT ID FROM GapList WHERE ParentID = ?", (glzsys_id,))
    gap_child_ids = cursor.fetchall()
    gap_child_ids.insert(0, [0,])

    # Step 4: For each child ID, find the corresponding rows in the "Glass" table and import them.
    origin = bottom_left_corner
    
    for glass_count, child_id in enumerate(glass_childs_ids):
        origin = insert_glass_polygon(cursor, gap_child_ids, glass_count, materials_element, polygons_element, child_id, origin)
    
    # Commit changes and close the connection
    conn.close()
    return root

def insert_glass_polygon(cursor, gap_child_ids, glass_count, materials_element, polygons_element, child_id, origin):
    child_id = child_id[0]
    # Query to find the corresponding glass elements
    cursor.execute("SELECT ID, Thickness FROM Glass WHERE ID = ?", (child_id,))
    glass_row = cursor.fetchone()
    
    glass_id, thickness = glass_row
    thickness = thickness*config["units"]["scale_factor"]
    print(origin)

    # Step 5: Define the polygon (rectangle) based on the thickness

    box_polygon = [origin, (origin[0]+thickness, origin[1]), (origin[0]+thickness, origin[1]+150), 
                (origin[0], origin[1]+150)]

    polygon = []
    for point in box_polygon:
        polygon.append((point[0] + gap_child_ids[glass_count][0]*config["units"]["scale_factor"], point[1]))

    origin = polygon[-3]

    material = {
                'Name': "G1_Glass_" + str(glass_count),
                'Conductivity': "1",
                'Emissivity': "0.900000",
                'Color': "0x00ffff",
                'Type': "0"
            }
        
    # Create the Material element
    material_element = ET.SubElement(
        materials_element,
        "Material",
        Name=material['Name'],
        Index=str(100+glass_count),
        Type=material['Type'],
        Conductivity=material['Conductivity'],
        Tir="0.000000",
        EmissivityFront=material['Emissivity'],
        EmissivityBack=material['Emissivity'],
        RGBColor=material['Color']
    )

    polygon_element = ET.SubElement(
        polygons_element,
        "Polygon",
        ID=str(glass_count + 200),
        Material="G1_Glass_" + str(glass_count),
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
            x=f"{x:.6f}",
            y=f"{y:.6f}"
        )

    polygon_element.set("id", str(glass_id))
    polygon_element.set("thickness", str(thickness))

    return origin
