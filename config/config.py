import yaml

with open("config/config.yaml", "r") as file:
    settings = yaml.safe_load(file)