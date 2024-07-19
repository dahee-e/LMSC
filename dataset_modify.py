file_path = f"./dataset/polblog/"
file_name = "community_origin.dat"

#rewrite the file
#original file consist of

import pandas as pd

# def read_community_file(file_path):
#     # Initialize an empty list to hold the rows of the dataframe
#     data = []
#
#     # Open the file and read line by line
#     with open(file_path, 'r') as file:
#         community_id = 1  # Start with community ID 1
#         for line in file:
#             # Split the line into node IDs
#             node_ids = line.strip().split()
#             # Create a row for each node ID with its community ID
#             for node_id in node_ids:
#                 data.append({'Node ID': int(node_id), 'Community ID': community_id})
#             # Increment the community ID for the next line
#             community_id += 1
#
#     # Create a DataFrame from the list of rows
#     df = pd.DataFrame(data)
#     return df
#
#
# def save_to_file(df, output_file):
#     df_sorted = df.sort_values(by='Node ID')
#     df_sorted.to_csv(output_file, index=False, header=False, sep='\t')



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