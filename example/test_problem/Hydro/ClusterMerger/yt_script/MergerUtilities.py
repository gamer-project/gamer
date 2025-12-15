#====================================================================================================
# Import
#====================================================================================================
import numpy as np
import h5py



#====================================================================================================
# Function
#====================================================================================================
def getRecordData(filename, target_cols=[]):
    '''
    Get data from Record__ClusterCenter
    '''
    data = [[] for _ in target_cols]
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line[0] == '#': continue
            temp = line.split()
            N = len(temp)
            for i, c in enumerate(target_cols):
               if c > N-1: continue
               data[i].append(float(temp[c]))

    return data


def getPassageTime(filename):
    time_col = 0
    c0_x_col = 2
    c0_y_col = 3
    c0_z_col = 4
    c1_x_col = 42
    c1_y_col = 43
    c1_z_col = 44
    cols = [time_col, c0_x_col, c0_y_col, c0_z_col, c1_x_col, c1_y_col, c1_z_col]
    data = getRecordData(filename, cols)
    N0 = len(data[1])
    N1 = len(data[4])

    prev = float('inf')
    for i in range(N1):
        dist2 = (data[1][i]-data[4][i])**2 + (data[2][i]-data[5][i])**2 + (data[3][i]-data[6][i])**2
        if dist2 > prev: return data[0][i]
        prev = dist2

    return -1.0


def getMergerTime(filename):
    '''
    Get BH merger time from Record__ClusterCenter
    '''
    time_col = 0
    c1_x_col = 42
    cols = [time_col, c1_x_col]
    data = getRecordData(filename, cols)
    N = len(data[1])
    return data[0][N-1]


def getClusterCen(data_path, cluster_id, start, end, diff=1):
    '''
    Get cluster center from Data_*
    '''
    cenX, cenY, cenZ = [], [], []
    idx_x = "ClusterCen_%d_0"%cluster_id
    idx_y = "ClusterCen_%d_1"%cluster_id
    idx_z = "ClusterCen_%d_2"%cluster_id
    for i in range(start, end+1, diff):
        f = h5py.File( data_path+"/Data_%06d"%i, "r" )
        if idx_x not in f["User"]["UserPara"].dtype.names:
            print("Can not find %s in Data_%06d"%(idx_x, i))
            break
        if idx_y not in f["User"]["UserPara"].dtype.names:
            print("Can not find %s in Data_%06d"%(idx_y, i))
            break
        if idx_z not in f["User"]["UserPara"].dtype.names:
            print("Can not find %s in Data_%06d"%(idx_z, i))
            break

        cenX.append( f["User"]["UserPara"][idx_x] )
        cenY.append( f["User"]["UserPara"][idx_y] )
        cenZ.append( f["User"]["UserPara"][idx_z] )
    return cenX, cenY, cenZ



#====================================================================================================
# Main (test)
#====================================================================================================
if __name__ == "__main__":
    print("Testing: getClusterCen")
    data_mcf5 = "/projectU/chunyenc/Develop/gamer_cool_core_new/bin/mcf5/"
    x, y, z = getClusterCen(data_mcf5, 1, 0, 100, 1)
    print(x)
    print(y)
    print(z)
    print("Test Done!!")

    print("Testing: getPassageTime")
    record_mcf5 = "/projectU/chunyenc/Develop/gamer_cool_core_new/bin/mcf5/Record__ClusterCenter"
    time = getPassageTime(record_mcf5)
    print(time)
    print("Test Done!!")

    print("Testing: getMergerTime")
    record_mcf5 = "/projectU/chunyenc/Develop/gamer_cool_core_new/bin/mcf5/Record__ClusterCenter"
    time = getMergerTime(record_mcf5)
    print(time)
    print("Test Done!!")
