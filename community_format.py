file_path = "./dataset/polbooks/"
file_name = "community.dat"
import pandas as pd

def read_community_file(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            node_id, community_id = line.strip().split()
            community_id = int(community_id)
            if community_id not in community_dict:
                community_dict[community_id] = []
            community_dict[community_id].append(node_id)
def save_to_file(df, output_file):
    with open(output_file, 'w') as file:
        for community_id, nodes in community_dict.items():
            file.write('\t'.join(nodes) + '\n')


community_dict = {}
output_file = 'community_origin.dat'  # Output file name
df = read_community_file(file_path+file_name)

# Save the DataFrame to a file
save_to_file(df, file_path+output_file)
print("File saved successfully!")
