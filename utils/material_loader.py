import csv

def from_csv(csv_file):
    # Load materials from the CSV file
    materials = []
    with open(csv_file, mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            materials.append(row)
    return materials