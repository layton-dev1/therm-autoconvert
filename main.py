from dxf_reader import reader
from utils import material_loader, utils
from xml_writer import writer
from config import config

def main():
    # Path to your DXF file
    data = config.data
    input_file_path = data["input_file_path"]
    output_file_path = data["output_file_path"]
    csv_file_path = data["materials_csv"]
    files = data["dxfs"]

    for file_name in files:
        print(f"DXF file '{file_name}' is being processed")
    
        input_dxf_file = input_file_path + file_name + ".dxf"
        output_xml_file = output_file_path + file_name + ".thmx"

        #Read and extract the drawing details
        polygons = reader.read_dxf(input_dxf_file)

        #This adds vertices necessary to properly define boundary conditions
        polygons = utils.add_vertices_to_edges(polygons)
        materials = material_loader.from_csv(csv_file_path)

        writer.create_xml(polygons, materials, output_xml_file)
        print(f"XML file '{output_xml_file}' has been created successfully.")
        
if __name__ == "__main__":
    main()