from collections import defaultdict
import pandas as pd

def processData(case):
    elements = pd.read_csv('validation_data/element_coords.csv', comment="#", header=0)
    nodes = pd.read_csv('validation_data/node_coords.csv', comment="#", header=0)
    if case=="bending":
        stress = pd.read_csv('validation_data/bending.csv',comment="#",header=0)
        rf = pd.read_csv('validation_data/bending_rf.csv',comment="#",header=0)
        u = pd.read_csv('validation_data/bending_u.csv', comment="#", header=0)
    elif case=="jam_bent":
        stress = pd.read_csv('validation_data/jam_bent.csv', comment="#", header=0)
        rf = pd.read_csv('validation_data/jam_bent_rf.csv', comment="#", header=0)
        u = pd.read_csv('validation_data/jam_bent_u.csv', comment="#", header=0)

    elif case=="jam_straight":
        stress = pd.read_csv('validation_data/jam_straight.csv', comment="#", header=0)
        rf = pd.read_csv('validation_data/jam_straight_rf.csv', comment="#", header=0)
        u = pd.read_csv('validation_data/jam_straight_u.csv', comment="#", header=0)
    else:
        print("Case not found")
        return

    dictElement = defaultdict(list)
    for index,row in elements.iterrows():
        dictElement[row['element']].append(row['node1'])
        dictElement[row['element']].append(row['node2'])
        dictElement[row['element']].append(row['node3'])
        dictElement[row['element']].append(row['node4'])
    print('Elements added to dict')
    print(dictElement)

    dictNode = defaultdict(list)
    for index,row in nodes.iterrows():
        dictNode[row['node']].append(row['x'])
        dictNode[row['node']].append(row['y'])
        dictNode[row['node']].append(row['z'])
    print('Nodes added to dict')
    print(dictNode)

    dictStress = defaultdict(list)
    for index, row in stress.iterrows():
        dictStress[row['element']].append(row['sMisesLoc1'])
        dictStress[row['element']].append(row['sMisesLoc2'])
        dictStress[row['element']].append(row['s12Loc1'])
        dictStress[row['element']].append(row['s12Loc2'])
    print('Stress added to dict')
    print(dictStress)

    dictRF = defaultdict(list)
    for index, row in rf.iterrows():
        dictRF[row['node']].append(row['rfMagnitudeLoc1'])
        dictRF[row['node']].append(row['rf1Loc1'])
        dictRF[row['node']].append(row['rf2Loc1'])
        dictRF[row['node']].append(row['rf3Loc1'])
    print('RF added to dict')
    print(dictRF)

    dictU = defaultdict(list)
    for index, row in u.iterrows():
        dictU[row['node']].append(row['uMagnitudeLoc1'])
        dictU[row['node']].append(row['u1Loc1'])
        dictU[row['node']].append(row['u2Loc1'])
        dictU[row['node']].append(row['u3Loc1'])
    print('U added to dict')
    print(dictU)

    dfObj = pd.DataFrame(columns=['element', 'node', 'x','y','z'])

    for el in dictElement:
        for nod in dictElement[el]:
            dfObj = dfObj.append({'element': el, 'node': nod, 'x': dictNode[nod][0], 'y': dictNode[nod][1], 'z': dictNode[nod][2],
                                  'sMisesLoc1': dictStress[el][0], 'sMisesLoc2': dictStress[el][1], 's12Loc1': dictStress[el][2], 's12Loc2': dictStress[el][3],
                                  'uMagnitudeLoc1': dictU[nod][0], 'u1Loc1': dictU[nod][1], 'u2Loc1': dictU[nod][2], 'u3Loc1': dictU[nod][3],
                                  'rfMagnitudeLoc1': dictRF[nod][0], 'rf1Loc1': dictRF[nod][1], 'rf2Loc1': dictRF[nod][2], 'rf3Loc1': dictRF[nod][3]
                                  }, ignore_index=True)
    print(dfObj)

    if case == "bending":
        dfObj.to_csv(r'validation_processed/bending_processed.csv')
        print("Data exported")
    elif case == "jam_bent":
        dfObj.to_csv(r'validation_processed/jam_bent_processed.csv')
        print("Data exported")
    elif case == "jam_straight":
        dfObj.to_csv(r'validation_processed/jam_straight_processed.csv')
        print("Data exported")
    return

#processData("jam_straight")
