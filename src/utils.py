import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class PTN:
    '''
    Class to initialize the public transport network which is build from
    nodes and verticses.
    V: set (1,2,3,4)
    E: dict: {(i, j): {'fmin': int, 'fmax': int}}
    '''
    def __init__(self, V: set, E: dict) -> None:
        if V is not None:
            self.V = V
        if E is not None:
            self.E = E

    def get_unique_edges(self) -> set:
        return set(self.E.keys())


class EAN:
    '''
    Class for initializing a event-activity network from a public transport 
    network to solve for example the PESP.

    Example of the data structure:
    Events: {1: (1, 'departure', 'S1_f'), ...} 
    Activities: {((1, 'departure', 'S1_f'), (3, 'arrival', 'S1_f')): {'type': 'drive', 'lb': 5, 'ub': 7, 'weight': 130}, ...}
    '''
    def __init__(self, Events: dict, Activities: dict) -> None:
        self.Events = Events
        self.Activities = Activities

    def get_unique_events(self) -> set:
        return set(self.Events.values())
    

def preprocess_data(BASE_PATH: str) -> (PTN, EAN, dict):
    '''
    Preprocesses the data of LinTim - Open Source Project into the desired
        data structure of a ptn, ean and linepool
    '''

    # BASE_PATH = '../example_ring/'
    # PATH_TO_DATASET = '../LinTim/datasets/toy/basis/'
    PATH_TO_LOAD = f'{BASE_PATH}Load_new.giv'
    PATH_TO_EDGE = f'{BASE_PATH}Edge.giv'
    PATH_TO_POOL = f'{BASE_PATH}Pool.giv'
    PATH_TO_POOL_COST = f'{BASE_PATH}Pool-Cost.giv'
    PATH_TO_EVENTS = f'{BASE_PATH}Events-periodic.giv'
    PATH_TO_ACTIVITIES = f'{BASE_PATH}Activities-periodic.giv'

    data_load = pd.read_csv(PATH_TO_LOAD, sep=';')
    data_edge = pd.read_csv(PATH_TO_EDGE, sep=';')
    data_pool = pd.read_csv(PATH_TO_POOL, sep=';')
    data_pool_cost = pd.read_csv(PATH_TO_POOL_COST, sep=';')
    ### Creating Line Pool
    line_pool_undirected = {}
    grouped = data_pool.groupby('# line-id')[' edge'].apply(list)
    line_dict = grouped.to_dict()

    # sequence contains only the edge ID's, but not the actual stops
    for key, value in line_dict.items():
        line_pool_undirected[f'S{key}'] = {
            'stations': value,
            'cost': data_pool_cost.loc[data_pool_cost['# line-id'] == key, ' cost'].values[0]
        }
    # now rewrite edgeIDs into stations
    for key, value in line_pool_undirected.items():
        stations = []
        for index in value['stations']:
            row = data_edge.iloc[index-1]
            stations.append(row['left-stop'])
            stations.append(row['right-stop'])
        line_pool_undirected[key]['stations'] = stations

    for key, value in line_pool_undirected.items():
        a = value['stations']
        for i in range(0, len(a)-2, 2):
            j, u, v = i+1, i+2, i+3
            if a[i] == a[v]:
                a[i], a[j] = a[j], a[i]
                a[u], a[v] = a[v], a[u]
            elif a[i] == a[u]:
                a[i], a[i+1] = a[i+1], a[i]
            elif a[j] == a[v]:
                a[u], a[v] = a[v], a[u]
            
        line_pool_undirected[key]['stations'] = remove_consecutive_numbers(a)
    # Now we have to create a forward and backward direction
    line_pool = {}
    for key, value in line_pool_undirected.items():
        line_pool[f'{key}_f'] = {
            'stations': list(map(int, value['stations'])),
            'cost': value['cost']
        }
        line_pool[f'{key}_b'] = {
            'stations': list(map(int, value['stations'][::-1])),
            'cost': value['cost']
        }
    linepool = line_pool

    ### Creating PTN (edges first, then vertices)
    data_edge.rename(columns={'#edge-ID': '# link_index'}, inplace=True)
    df_ptn_edges = pd.merge(left=data_edge, right=data_load, on='# link_index')
    ptn_edges_undirected = {}
    for index, row in df_ptn_edges.iterrows():
        i, j = row['left-stop'], row['right-stop']
        fmin, fmax = row['1'], row[' max_freq']
        ptn_edges_undirected[(int(i), int(j))] = {'fmin': fmin, 'fmax': fmax}

    # again we have to create a forward and backward direction
    ptn_edges = {}
    for key, value in ptn_edges_undirected.items():
        ptn_edges[(key[0], key[1])] = value
        ptn_edges[(key[1], key[0])] = value

    # PTN vertices
    ptn_vertices = set()
    for (i, j) in ptn_edges:
        ptn_vertices |= {i, j}
    ptn = PTN(ptn_vertices, ptn_edges)

    ### Creating EAN
    data_events = pd.read_csv(PATH_TO_EVENTS, sep=';')

    # data cleaning
    data_events[' type'] = data_events[' type'].apply(lambda x: x.replace('"', ''))
    data_events[' type'] = data_events[' type'].apply(lambda x: x.replace(' ', ''))

    ean_events = {}
    for index, row in data_events.iterrows():
        station, event_type = row[' stop-id'], row[' type']
        if row[' line-direction'] == ' >':
            line = f"S{row[' line-id']}_f"
        else: 
            line = f"S{row[' line-id']}_b"
        ean_events[index+1] = (station, event_type, line)

    data_activities = pd.read_csv(PATH_TO_ACTIVITIES, sep=';')

    # data cleaning
    data_activities[' type'] = data_activities[' type'].apply(lambda x: x.replace('"', ''))
    data_activities[' type'] = data_activities[' type'].apply(lambda x: x.replace(' ', ''))

    ean_activities = {}

    for index, row in data_activities.iterrows():
        i = ean_events[row[' from_event']] 
        j = ean_events[row[' to_event']]
        ean_activities[(i,j)] = {
            'type': row[' type'],
            'lb': row[' lower_bound'],
            'ub': row[' upper_bound'],
            'weight': row[' passengers']
        }
        
    ean = EAN(Events=ean_events, Activities=ean_activities)
    return ptn, ean, line_pool


def remove_consecutive_numbers(numbers):
    result = []
    for num in numbers:
        if not result or num != result[-1]:  # Add number if it's not the same as the last in result
            result.append(num)
    return result


def pesp_part_to_betas(
        ean: dict,
        schedule: dict,
        linepool: dict,
        y_connection_param: dict | None,
        z_modulo_param: dict,
        T=60) -> np.array:
    '''
    Calculates the optimal beta values for the integrated model with T as
    time priod. It is done by iterating over all activities and addding
    its weight times duration in the corresponding entry in a |L|x|L|
    matrix.
    '''
    look_up = {line: index for index, line in enumerate(linepool)}
    betas = np.zeros((len(linepool), len(linepool)), dtype=int)
    y, z = y_connection_param, z_modulo_param

    for activity in ean.Activities:
        i, j = activity[0], activity[1]
        line1, line2 = i[2], j[2]
        weight = ean.Activities[activity]['weight']

        x_ij = 0 # if we solved the IM and line has been chosen
        if (y is not None) and (y[activity]-1 == 0):
            x_ij = (schedule[j] - schedule[i] + z[activity]*T)
        elif (y is None): # if we only have solved the PESP in the init step
            x_ij = (schedule[j] - schedule[i] + z[activity]*T)
        betas[look_up[line1], look_up[line2]] += weight * x_ij

    return betas


def pesp_weights_to_betas(  # This has to be added in utils
        ean: dict,
        linepool: dict,
        T=60) -> np.array:
    '''
    Calculates initial beta values based on the given weights for 
     the PESP.
    '''
    look_up = {line: index for index, line in enumerate(linepool)}
    betas = np.zeros((len(linepool), len(linepool)), dtype=int)

    for activity in ean.Activities:
        i, j = activity[0], activity[1]
        line1, line2 = i[2], j[2]
        weight = ean.Activities[activity]['weight']
        betas[look_up[line1], look_up[line2]] += weight

    return betas


def is_positive_definite(matrix: np.array) -> bool:
    '''Checks if all eigenvalues of the input matrix are positive'''
    eigenvalues = np.linalg.eigvals(matrix)
    return np.all(eigenvalues > 0)


def modify_matrix_efficiently(input_array, matrix):

    matrix[input_array, :] = 0
    matrix[:, input_array] = 0

    for i in range(len(matrix)):
        if matrix[i, i] == 0:
            matrix[i, i] = 1

    return matrix



