import xml.etree.ElementTree as ET
from datetime import datetime
from utils import utils
from xml_writer import igu_importer

def create_xml(polygons, output_file):
    # Create the root element
    root = create_headings()
    materials_element = ET.SubElement(root, "Materials")
    boundaryconditions_element = ET.SubElement(root, "BoundaryConditions")
    polygons_element = ET.SubElement(root, "Polygons")
    boundaries_element = ET.SubElement(root, "Boundaries")
    
    # Create polygon, igu, and boundary elements
    polygon_xml_section(root, polygons, materials_element, polygons_element)

    glass_left, glass_right = igu_importer.import_glass_from_mdb_to_xml(root, polygons_element, materials_element)
    
    boundary_xml_section(root, materials_element, polygons_element, boundaryconditions_element, boundaries_element, glass_left, glass_right)

    # Write the XML to a file
    tree = ET.ElementTree(root)
    tree.write(output_file, encoding="utf-8", xml_declaration=True)

def create_headings():
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
    ET.SubElement(root, "MeshControl", MeshLevel="9", ErrorCheckFlag="1", ErrorLimit="10.000000", MaxIterations="5", CMAflag="0")
    
    # Add RadianceModeBC element
    radiance_mode_bc = ET.SubElement(root, "RadianceModeBC")
    ET.SubElement(radiance_mode_bc, "RadianceModeInBCName").text = "ShadeInE"
    ET.SubElement(radiance_mode_bc, "RadianceModeInBCTag").text = "ShadeInETag"
    ET.SubElement(radiance_mode_bc, "RadianceModeOutBCName").text = "ShadeOutE"
    ET.SubElement(radiance_mode_bc, "RadianceModeOutBCTag").text = "ShadeOutETag"
    
    return root
    
def create_material(materials_element, material, is_cavity=False):
    if(len(materials_element) > 0):
        last_id = int(materials_element[-1].attrib.get("Index"))
    else:
        last_id = 1
    
    if(is_cavity):
        # Create frame material element
        material_element = ET.SubElement(
            materials_element,
            "Material",
            Name=material['Name'],
            Index=str(last_id+1),
            Type=material['Type'],
            Conductivity=material['Conductivity'],
            Tir="1.000000",
            EmissivityFront=material['Emissivity'],
            EmissivityBack=material['Emissivity'],
            RGBColor=material['Color'],
            CavityModel=material['CavityModel']
        ) 
    else:
        # Create material element
        material_element = ET.SubElement(
            materials_element,
            "Material",
            Name=material['Name'],
            Index=str(last_id + 1),
            Type=material['Type'],
            Conductivity=material['Conductivity'],
            Tir="0.000000",
            EmissivityFront=material['Emissivity'],
            EmissivityBack=material['Emissivity'],
            RGBColor=material['Color']
        ) 
        
        # Material elements need property subelements for Front and Back sides
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
    
    return material_element   

def create_polygon_element(polygons_element, material, polygon):
    if(len(polygons_element) > 0):
        last_id = int(polygons_element[-1].attrib.get("ID"))
    else:
        last_id = 1
        
    # Add the polygon
    polygon_element = ET.SubElement(
        polygons_element,
        "Polygon",
        ID=str(last_id + 1),
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

def polygon_xml_section(root, polygons, materials_element, polygons_element):
    unique_materials = []
    cavity_count = 2

    for i, (polygon, layer_name) in enumerate(polygons):
        # Find the closest matching material for the layer name
        if "frame" in layer_name.lower().strip() and "cavity" in layer_name.lower().strip():
            #Define the frame cavity material
            material = {
                'Name': f"Frame Cavity NFRC 100" + "_Cavity_" + str(cavity_count),
                'Conductivity': "-1",
                'Emissivity': "0.900000",
                'Color': "0x00FF00",
                'Type': "1",
                'CavityModel': "4"
            }

            #Create material Element
            create_material(materials_element, material, is_cavity=True)
            cavity_count = cavity_count + 1
        else:
            #Define the material
            material = utils.find_closest_material(layer_name)

            #If none found add in placeholder
            if not material:
                    print(f"Warning: No matching material found for layer '{layer_name}'. Using default values.")
                    material = {
                        'Name': f"{layer_name}",
                        'Conductivity': "0.170000",
                        'Emissivity': "0.900000",
                        'Color': "0x000000",
                        'Type': "0"
                    }
            
            # Add material to the materials section if not already existing
            if layer_name not in unique_materials:
                # Create the Material element
                create_material(materials_element, material)

                # Add to unique materials list so we don't add it twice
                unique_materials.append(layer_name)
        
        create_polygon_element(polygons_element, material, polygon)

    return root

def boundary_xml_section(root, materials_element, polygons_element, boundaryconditions_element, boundaries_element, glass_left, glass_right):
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
        #plot_polygon_debug(polygon_points, touching_edges)

        create_solid_boundary(polygon_points, touching_edges, boundaries_element, polygon_id, material_name, emissivity, polygons_element, glass_left, glass_right)

        #Frame cavities need an extra boundary condition called Frame Cavity Surface 
        if "Frame Cavity NFRC 100" in material_name:
            polygon_points.reverse()
        
            touching_edges = utils.find_touching_edges(polygons_element, polygon_id, polygon_points)
            create_frame_cavity_boundary(polygon_points, materials_element, touching_edges, boundaries_element, polygon_id, emissivity, polygons_element)

    return root

def create_frame_cavity_boundary(polygon_points, materials_element, touching_edges, boundaries_element, polygon_id, emissivity, polygons_element):
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
        
        create_bc_polygon(boundaries_element, bc_polygon_name, material_name, polygon_id, p1, p2, emissivity, polygons_element, "1")

def create_solid_boundary(polygon_points, touching_edges, boundaries_element, polygon_id, material_name, emissivity, polygons_element, glass_left, glass_right):
    for side in range(len(polygon_points)):  # Loop through each side of the polygon
        p1 = polygon_points[side]
        p2 = polygon_points[(side + 1) % len(polygon_points)]
        edge = (p1, p2)

        if not edge in touching_edges: #Boundaries for solid materials should only appear when edge it not touching another polygon
            # Check the boundary's orientation to determine if it's exterior or interior
            if p1[0] < glass_left or p2[0] < glass_left:  # If the y-coordinate of p1 is less than p2, it's facing the exterior (upward)
                bc_polygon_name = "NFRC 100-2010 Exterior"
                ufactortag = "SHGC Exterior"
            elif p1[0] > glass_right or p2[0] > glass_right:  # If the y-coordinate of p1 is greater than p2, it's facing the interior (downward)
                bc_polygon_name = "Interior Thermally Broken Frame (convection only)"
                ufactortag = "Frame"
            elif p1[1] == p2[1]:
                bc_polygon_name = "Adiabatic"
                ufactortag = ""
            else:  # Handle vertical boundaries (when x-coordinates are the same)
                if p1[0] < glass_left:  # If the y-coordinate of p1 is less than p2, it's facing the exterior (upward)
                    bc_polygon_name = "NFRC 100-2010 Exterior"
                    ufactortag = "SHGC Exterior"
                elif p1[0] > glass_right:  # If the y-coordinate of p1 is greater than p2, it's facing the interior (downward)
                    bc_polygon_name = "Interior Thermally Broken Frame (convection only)"
                    ufactortag = "Frame"
                else:  # Handle the case where the points are exactly the same
                    bc_polygon_name = "Adiabatic"
                    ufactortag = ""
            #bc_polygon_name = "Adiabatic"
            create_bc_polygon(boundaries_element, bc_polygon_name, material_name, polygon_id, p1, p2, emissivity, polygons_element, ufactortag=ufactortag)

def create_bc_polygon(boundaries_element, bc_polygon_name, material_name, polygon_id, p1, p2, emissivity, polygons_element, enclosure="0", ufactortag=""):
    if(len(boundaries_element) > 0):
        last_id = int(boundaries_element[-1].attrib.get("ID"))
    else:
        #First boundary condition polygon ID should start after regular polygon IDs
        last_id = int(polygons_element[-1].attrib.get("ID"))
        
    # Create the bondary condition polygon element
    bc_polygon = ET.SubElement(
        boundaries_element,
        "BCPolygon",
        ID=str(last_id + 1),
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
