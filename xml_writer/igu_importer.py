import pyodbc
from config import config
import xml.etree.ElementTree as ET
from components import thmx_reader

def import_glass_from_mdb_to_xml(root, polygons_element, materials_element):
    mdb_file_location = r"C:\Users\tyler.henderson\Documents\DXFtoTHERM\therm-autoconvert\lib\data\therm-autoconvert - DB.mdb"
    mdb_glzsys_id = 2
    spacer_file = config.data["spacer_file_path"]
    
    # Step 1: Open the MDB file and access the "GlzSys" table to find the "GlzSys" ID.
    conn_str = r'DRIVER={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=' + mdb_file_location
    conn = pyodbc.connect(conn_str)
    cursor = conn.cursor()

    # Search for the rectangle on the Glass1 layer
    glass_layer = "GLASS1"
    bottom_left_corner = ""
    for polygon_element in polygons_element:
        material_name = polygon_element.attrib.get('Material')

        if glass_layer.strip().lower() in material_name.strip().lower():
            bottom_left_corner = (float(polygon_element.findall("Point")[0].attrib["x"]), 
                                  float(polygon_element.findall("Point")[0].attrib["y"]))

            for point_element in polygon_element.findall("Point"):
                x = float(point_element.attrib["x"])
                y = float(point_element.attrib["y"])
                
                if x < bottom_left_corner[0] or y < bottom_left_corner[1]:
                    bottom_left_corner = (x, y)
                        
            polygons_element.remove(polygon_element)

    glass_left = bottom_left_corner[0]

    igu = {
        "glass": [],
        "gap": []
    }

    sql_query = """
    SELECT 
        GlzSys.ID AS GlzSysID,
        GlassList.ID AS GlassListIndex,
        GlassList.Flip,
        Glass.ID AS GlassID,
        Glass.Name AS GlassName,
        Glass.Thickness AS GlassThickness,
        Glass.emis1,
        Glass.emis2,
        Glass.Conductivity
    FROM ((GlzSys
        LEFT JOIN GlassList ON GlassList.ParentID = GlzSys.ID)
        LEFT JOIN Glass ON Glass.ID = GlassList.ID)
    WHERE 
        GlzSys.ID = ?
    """
    
    cursor.execute(sql_query, mdb_glzsys_id)
    rows = cursor.fetchall()

    # Print each row
    for row in rows:
        # GlzSys
        glzsys_id = row.GlzSysID

        # GlassList
        glasslist_index = row.GlassListIndex
        flip = row.Flip
        
        if flip == 1:
            emis1, emis2 = emis2, emis1

        # Glass
        glass_id = row.GlassID
        glass_name = row.GlassName
        glass_thickness = row.GlassThickness
        emis1 = row.emis1
        emis2 = row.emis2
        conductivity = row.Conductivity
        
        # TODO make rounding according to ini file?
        igu["glass"].append({
            "id": glass_id,
            "name": glass_name,
            "index": glasslist_index,
            "thickness": round(glass_thickness, 0),
            "emis1": emis1,
            "emis2": emis2,
            "conductivity": conductivity,
        })
        
    sql_query = """
    SELECT 
        GlzSys.ID AS GlzSysID,
        GapList.Index AS GapListIndex,
        GapList.Thickness AS GapThickness,
        GapList.Keff,
        Gap.ID AS GapID,
        Gap.Name AS GapName
    FROM ((GlzSys
        LEFT JOIN GapList ON GapList.ParentID = GlzSys.ID)
        LEFT JOIN Gap ON Gap.ID = GapList.ID)
    WHERE 
        GlzSys.ID = ?
    """
    
    cursor.execute(sql_query, mdb_glzsys_id)
    rows = cursor.fetchall()
    
    for row in rows:
        # GapList
        gaplist_index = row.GapListIndex
        gap_thickness = row.GapThickness
        keff = row.Keff

        # Gap
        gap_id = row.GapID
        gap_name = row.GapName
        
        igu["gap"].append({
            "id": gap_id,
            "index": gaplist_index,
            "name": gap_name,
            "thickness": gap_thickness,
            "keff": keff,
        })
        
    highest_point_y = float('-inf')
    for polygon in polygons_element:
        for point in polygon.findall("Point"):
            y_value = float(point.attrib["y"])
            if(y_value>highest_point_y):
                highest_point_y = y_value
    
    glass_right = insert_glass_polygon(igu, materials_element, polygons_element, bottom_left_corner, 9.7, float(highest_point_y))
    
    # Commit changes and close the connection
    conn.close()
    return glass_left, glass_right[0]

def insert_glass_polygon(igu, materials_element, polygons_element, origin, spacer_height, highest_point_y):
    # Loop through glass
    polygon = []
    if(spacer_height + origin[1] > highest_point_y):
        print("Spacer is taller than the highest point")
        highest_point_y = spacer_height + origin[1]
    igu_height = (highest_point_y-origin[1])+150
    edge_height = (highest_point_y-origin[1])+63.5
    
    print(f"Highest Point: {highest_point_y}")
    print(f"IGU Height: {igu_height}")
    print(f"Edge Height: {edge_height}")
    print(f"Origin: {origin[0]}, {origin[1]}")
    
    for i, glass in enumerate(igu["glass"]):
        # Create a rectangle for class starting at origin
        glass_polygon = [origin, (origin[0]+glass["thickness"], origin[1]), (origin[0]+glass["thickness"], origin[1]+igu_height), 
                (origin[0], origin[1]+igu_height)]
        
        if(i==len(igu["glass"])-1):
            glass_polygon.insert(2, (origin[0]+glass["thickness"], origin[1]+edge_height))
            
        origin = (origin[0]+glass["thickness"], origin[1])

        create_glass_element(materials_element, polygons_element, glass_polygon, glass, i)
        if not i == len(igu["glass"]) - 1:
            gap = igu["gap"][i]
            gap_polygon = [(origin[0], origin[1]+spacer_height), (origin[0]+gap["thickness"], origin[1]+spacer_height), (origin[0]+gap["thickness"], origin[1]+igu_height), 
                (origin[0], origin[1]+igu_height)]
            create_gap_element(materials_element, polygons_element, gap_polygon, gap, i)
            thmx_reader.create_spacer(materials_element, polygons_element, origin, gap["thickness"])
            origin = (origin[0]+gap["thickness"], origin[1])

    return origin

def create_gap_element(materials_element, polygons_element, gap_polygon, gap, i):
    material = {
                'Name': "G1_Gap_" + str(i),
                'Conductivity': str(gap["keff"]),
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
        NSides=str(len(gap_polygon)),
        Type="1",
        units="mm"
    )
    
    # Add each point to the polygon
    for j, (x, y) in enumerate(gap_polygon):
        ET.SubElement(
            polygon_element,
            "Point",
            index=str(j),
            x=f"{x:.6f}",
            y=f"{y:.6f}"
        )

def create_glass_element(materials_element, polygons_element, glass_polygon, glass, i):
    material = {
                'Name': "G1_Glass_" + str(i),
                'Conductivity': str(glass["conductivity"]),
                'Emissivity': "0.900000",
                'Color': "0x00ffff",
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
        NSides=str(len(glass_polygon)),
        Type="1",
        units="mm"
    )
    
    # Add each point to the polygon
    for j, (x, y) in enumerate(glass_polygon):
        ET.SubElement(
            polygon_element,
            "Point",
            index=str(j),
            x=f"{x:.6f}",
            y=f"{y:.6f}"
        )