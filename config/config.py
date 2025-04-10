import yaml

with open("config/config.yaml", "r") as file:
    config = yaml.safe_load(file)
    
data = config["data"]
units = config["units"]