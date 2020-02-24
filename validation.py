from collections import defaultdict
import pandas as pd

elements = pd.read_csv('validation_data/element_coords.csv',comment="#",header=0)
nodes = pd.read_csv('validation_data/node_coords.csv',comment="#",header=0)
#bendingStress = pd.read_csv('validation_data/bending.csv',comment="#")
#print(bendingStress)


dict = defaultdict(list)
for index,row in elements.iterrows():
    dict[row['element']].append(row['node1'])
    dict[row['element']].append(row['node2'])
    dict[row['element']].append(row['node3'])
    dict[row['element']].append(row['node4'])
print('Elements added to dict')
print(dict)

dict2 = defaultdict(list)
for index,row in nodes.iterrows():
    dict2[row['node']].append(row['x'])
    dict2[row['node']].append(row['y'])
    dict2[row['node']].append(row['z'])
print('Nodes added to dict')
print(dict2)





dfObj = pd.DataFrame(columns=['element', 'node', 'x','y','z'])

for el in dict:
    for nod in dict[el]:
        dfObj = dfObj.append({'element': el, 'node': nod, 'x': dict2[nod][0], 'y': dict2[nod][1], 'z': dict2[nod][2]}, ignore_index=True)

print(dfObj)