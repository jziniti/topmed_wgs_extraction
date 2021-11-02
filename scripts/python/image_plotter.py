import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
from itertools import combinations
import math

class GenerateGraph:
    # initialize object with manifest data and cdata
    def __init__(self, manifest_df, cdata_df):
        self.df = manifest_df
        self.dfc = cdata_df
        self.node_mappings = { 
            "ObservedGender": { "F": "o", "M": "s", "U": "d", "nan": "v" },
            "SampleLabel": { "Number": "s", "NWD": "d", "TOE": "v" },
            "SampleAlias": { "Number": "s", "NWD": "d", "TOE": "v" },
            "SampleTypeCode": { "d": "s", "m": "p", "r": "d" },
            "dataset_id": { "tmp0": "s", "tmp1": "d", "tmp2": "v", "tmp3": "o" },
            "S_SAMPLEID": { "tmp0": "s", "tmp1": "d", "tmp2": "v", "tmp3": "o" }
        }
        self.mapping_funcs = {
            "ObservedGender": lambda x: x,
            "SampleLabel": self.samplelabel_func,
            "SampleAlias": self.samplelabel_func,
            "SampleTypeCode": lambda x: x,
            "dataset_id": self.dataset_id_func,
            "S_SAMPLEID": lambda x: x[0:3] + x[4]
        }
    
    # method to call to generate graph
    def graph(self, group_id, shape_by):
        # filtering by group_id
        query_groupid = f'group_id == {group_id}'
        df_filter_groupid = self.df.query(query_groupid)
        # get nodes
        nodes = df_filter_groupid["SampleLabel"].values

        # group nodes by shape mappings
        shape_nodes = self.map_nodes(nodes, shape_by)

        # getting all combinations of samplelabels
        edges, labels = self.generate_edges(nodes)

        self.generate_graph(nodes, shape_nodes, edges, labels, group_id)
        
    def generate_graph(self, nodes, shape_nodes, edges, labels, groupid):
        graph_path = "images/" + "group" + str(groupid) + ".dot"

        fig = plt.figure(figsize=(12,12))
        ax = plt.subplot(111)
        ax.set_title('Output Graph', fontsize=10)

        G = nx.Graph()
        G.add_nodes_from(nodes)
        G.add_edges_from(edges)
        pos = nx.spring_layout(G, k=0.6)

        for shape in shape_nodes.keys():
            nx.draw_networkx_nodes(G, pos, nodelist=shape_nodes[shape]["nodes"], node_shape=shape, label=shape_nodes[shape]["attribute"])
            nx.draw_networkx_labels(G, pos, labels={node:node.split("_")[0] for node in G.nodes()})
    
        nx.draw_networkx_edges(G, pos, edgelist=edges)
        nx.draw_networkx_edge_labels(G,pos,edge_labels=labels,font_color='red',label_pos=0.4)
        nx.drawing.nx_pydot.write_dot(G, graph_path)
        #plt.legend()
        #plt.savefig(graph_path, format="PNG")
        #plt.close(fig)
        
    # returns a list of edges and a list of edge labels
    # given a list of nodes
    def generate_edges(self, n):
        all_edges = list(combinations(n, 2))
        filtered_edges = []
        labels = {}
        
        for edge in all_edges:
            query_string = f'ID1 == "{edge[0]}" and ID2 == "{edge[1]}"'
            concord = self.dfc.query(query_string)["Concord"].values[0]
            
            # only create edge for pairs that have concord values
            if not math.isnan(concord):
                labels[edge] = concord
                filtered_edges.append(edge)
                
        return filtered_edges, labels
    
    # returns a dictionary with nodes grouped by shape such as { "shape": { "nodes": [nodes], "attribute": "attr" }... }
    # given a list of nodes
    def map_nodes(self, n, by):
        shape_nodes = {}
        
        try:
            curr_map = self.node_mappings[by]
        except Exception:
            raise Exception("shape_by is not a valid option!")
        
        for node in n:
            # get attribute for node
            query_string = f'SampleLabel == "{node}"'
            # cast to String for nan case
            attr = str(self.df.query(query_string)[by].values[0])

            # get final attribute String
            # will not throw exception because of try except above
            final_attr = self.mapping_funcs[by](attr)

            # get mapped shape
            shape = curr_map[final_attr]
            
            if shape in shape_nodes:
                shape_nodes[shape]["nodes"].append(node)
            else:
                shape_nodes[shape] = { "nodes": [node], "attribute": final_attr }
                
        return shape_nodes
    
    # allows a client to add shape by functionality by providing the attribute name exactly as it appears in the 
    # data, the mapping of an attribute to shapes as a dictionary { "final attribute": "shape"... }, 
    # and a function(attribute) -> final attribute
    def add_shape_by(self, attr_name, attr_mapping, attr_func):
        if attr_name in self.node_mappings:
            raise Exception("Shape by attribute is already supported")
        
        self.node_mappings[attr_name] = attr_mapping
        self.mapping_funcs[attr_name] = attr_func

    def samplelabel_func(self, attr):
        prefix = attr[0:3]
        
        if prefix == "NWD":
            return "NWD"
        elif prefix == "TOE":
            return "TOE"
        else:
            return "Number"
        
    def dataset_id_func(self, attr):
        num = int(attr)

        if num == 0:
            return "tmp0"
        elif num == 1:
            return "tmp1"
        elif num == 2:
            return "tmp2"
        elif num == 3:
            return "tmp3"

# concord values must be floats, or nan(some c_data files have strings in concord column)
# assuming that c_data has concord value for both combinations of ID1 and ID2 (line 74)
# S_SAMPLEID must be in format "tmp/#..."

if __name__ == "__main__":
    if 'snakemake' in locals():
        group_ids = snakemake.params.group_ids
        shape_by = snakemake.params.shape_by
        manifest_path = snakemake.input.annotated_manifest
        c_data_path = snakemake.input.concordance_data
    else:
        manifest_path = "/proj/regeps/regep00/studies/TopMed/data/dna/whole_genome/TopMed/data/freezes/freeze.10.cdnm/plate103/ANNOTATED_MANIFEST.csv"
        c_data_path = "/proj/regeps/regep00/studies/TopMed/data/dna/whole_genome/TopMed/data/freezes/freeze.10.cdnm/plate103/c_data.build_sample_groups"
        shape_by = "dataset_id"
        group_ids = [1009, 1010, 102, 1081, 1088]
    
    df = pd.read_csv(manifest_path)
    dfc = pd.read_csv(c_data_path)
    
    myGraph = GenerateGraph(df, dfc)
    
    if group_ids[0] == "all":
        # satisfy requirement for outputing all.txt file
        f = open("images/groupall.dot", "a")
        f.write(str(len(group_ids)))
        f.close()
        
        # generate graphs for all group ids
        all_groups = df["group_id"].unique()
        for group in all_groups:
            myGraph.graph(group, shape_by)
        
    else:
        for group_id in group_ids:
            myGraph.graph(group_id, shape_by)
