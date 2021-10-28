import matplotlib.pyplot as plt
import networkx as nx


if __name__ == "__main__":
    print("hello world")
    
    group_id = snakemake.params.group_ids
    
    file_path = f'images/{group_id}.txt'
    f = open(file_path, "a")
    f.close()