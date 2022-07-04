'''
Created on Jul 8, 2021

@author: nassim
'''
import time
from math import floor, isinf, sqrt
from utilities import indexToCoordinates, coordinatesToIndex,indexToCoordinates3D, d_


def computeGradientDescent(s_idx, g_idx, nodes, d_, dim):
    print '--> OFPSearch processing...'
    cur_point, cur_coord, step, idx, gradient, path = [0.0 for _ in range(dim)], indexToCoordinates3D(s_idx, dim),\
                                                            1.0, s_idx, [0.0 for _ in range(dim)], []
    
    print cur_point, len(cur_point), dim
    #time.sleep(1000)
    
    for i in range(dim):
        cur_point[i] = float(cur_coord[i])
    print(cur_point)    
    path.append([cur_point[0],cur_point[1],cur_point[2]])
    
    while(idx != g_idx):
        
        gradient[0] = -nodes[idx-1].cost /2 + nodes[idx+1].cost /2 
        
        if (isinf(gradient[0]) and gradient[0]<0):
            gradient[0] = -1
        
        elif (isinf(gradient[0]) and gradient[0]>0):
            gradient[0] = 1
        
        max_grad = gradient[0]
        
        for i in range(1, dim):

            if (isinf(gradient[i]) and gradient[i]<0):
                gradient[i] = -1
            
            elif (isinf(gradient[i]) and gradient[i]>0):
                gradient[i] = 1
            
            if abs(max_grad)<abs(gradient[i]):
                max_grad = gradient[i]
                
        for i in range(dim):
            cur_point[i] =((cur_point[i] - step*gradient[i]/abs(max_grad)))
            cur_coord[i] = int(cur_point[i]+0.5)
            
        print(cur_point)
        path.append([cur_point[0],cur_point[1],cur_point[2]])
        idx = coordinatesToIndex(cur_coord, dim)    
        
    cur_coord=indexToCoordinates3D(idx, dim)
    
    for i in range(dim):
        cur_point[i] = float(cur_coord[i])
       
    path.append([cur_point[0],cur_point[1],cur_point[2]])    
    print 'path found !',path
    return path

def Gradient(START, GOAL, nodes, d_):


    temps = []
    path = []
    current_point=[0.0,0.0,0.0]
    idx = START   
    grads = [0.0,0.0,0.0]
    current_coord = indexToCoordinates3D(idx, 3)
    print current_coord, nodes[START].indice,nodes[START].z,nodes[START].y,nodes[START].x, START
    time.sleep(2)
    current_point[0] = float(current_coord[0])
    current_point[1] = float(current_coord[1])
    current_point[2] = float(current_coord[2])
    
    path.append([current_point[0],current_point[1],current_point[2]])
    temps.append(nodes[idx].cost)
    
   
    while (idx !=GOAL):
        idx_x_l = idx - 1
        idx_x_u = idx+1
        
        index1 = (max(idx_x_l,0))
        index2 = (max(idx_x_u,0))

        grad = 0.0
        cnt = 9
    
        for i in range(-1,2):
            for j in range(-1,2):
                for t in range(-1,2):
                    
                    index_n_l = index1 + i+j*int(d_[0])+t*int(d_[1])
                    
                    index_s_l = min((max(index_n_l,0)),d_[2]-1)
                    value_l = nodes[index_s_l].cost
                    
                    index_n_u = index2 + i+j*int(d_[0])+t*int(d_[1])
                    index_s_u = min((max(index_n_u,0)),d_[2]-1)
                    value_u = nodes[index_s_u].cost

                    
                    if isinf(value_l) or isinf(value_u):
                        cnt-=1
                
                    else : 
                        grad +=(value_u-value_l)/2.0
        grad/=cnt
        grads[0] = grad      
    
    
        ######for Y dimension
        idx_y_l = idx - d_[0]
        idx_y_u = idx + d_[0]    
        
        index1 =(max(idx_y_l,0))
        index2 =(max(idx_y_u,0))
        grad = 0.0
        cnt = 9
        
        for i in range(-1,2):
            for j in range(-1,2):
                for t in range(-1,2):
                    
                    index_n_l = index1 + i+ j*d_[0]+t**d_[1]
                    index_s_l = min((max(index_n_l,0)),d_[2]-1)
                    
                    value_l  = nodes[index_s_l].cost
                    
                    
                    index_n_u = index2 + i + j*d_[0]+t*d_[1]
                    index_s_u = min((max(index_n_u,0)),d_[2]-1)
                    
                    value_u  = nodes[index_s_u].cost
        
        
                    if isinf(value_l) or isinf(value_u):
                        cnt-=1
                
                    else : 
                        grad +=(value_u-value_l)/2.0
        grad/=cnt
        grads[1] = grad
    
        max_dif = -1000.0
        best_idx = idx
        
    
        ######for Z dimension
        idx_z_l = idx - d_[1]
        idx_z_u = idx + d_[1]    
        
        indx1 =(max(idx_z_l,0))
        indx2 =(max(idx_z_u,0))
        grad = 0.0
        cnt = 9
        
        for i in range(-1,2):
            for j in range(-1,2):
                for t in range(-1,2):
                    
                    index_n_l = indx1 + i+ j*d_[0]+t*d_[1]
                    index_s_l = min((max(index_n_l,0)),d_[2]-1)
                    
                    value_l  = nodes[index_s_l].cost
                    
                    
                    index_n_u = indx2 + i*1 + j*d_[0]+t*d_[1]
                    index_s_u = min((max(index_n_u,0)),d_[2]-1)
                    
                    value_u  = nodes[index_s_u].cost
        
        
                    if isinf(value_l) or isinf(value_u):
                        cnt-=1
                
                    else : 
                        grad +=(value_u-value_l)/2.0
        grad/=cnt
        grads[2] = grad
    
        max_dif = -1000.0
        best_idx = idx

        
        
        if isinf(grads[0]) or isinf(grads[1]) or isinf(grads[2]) or grads[0]==0.0 or grads[1]==0.0 or grads[2]==0.0:
            
            for i in range(-1,2):
                for j in range(-1,2):
                    for t in range(-1,2):
                       
                        if(i==0 and j==0 and t==0):
                            continue
                        """
                        idx_delta = i+d_[0]*j
                        idx_n = idx+idx_delta
                        
                        idx_tmp = (idx_n)
                        """  
                        idx_tmp = min((max(idx+i+d_[0]*j+d_[1]*t,0)),d_[2]-1)
                        
                        if isinf(nodes[idx_tmp].cost):
                            continue
                        dif = (-nodes[idx_tmp].cost + nodes[idx].cost)/sqrt(i*i+j*j+t*t)
                        if dif > max_dif:
                            max_dif = dif
                            best_idx = idx_tmp
                
            best_idx = max(best_idx,0)
            
            current_coord = indexToCoordinates3D(best_idx, 3)
            current_point[0] = float(current_coord[0])
            current_point[1] = float(current_coord[1])
            current_point[2] = float(current_coord[2])
        else:
            
            grad_norm = sqrt(grads[0]*grads[0]+grads[1]*grads[1]+grads[2]*grads[2])
            
            for i in range(3):
                current_point[i] = current_point[i] -grads[i]/grad_norm
                current_coord[i] = int(current_point[i]+0.5)
        
        
        print(best_idx, [current_point[0],current_point[1],current_point[2]]) 
        path.append([current_point[0],current_point[1],current_point[2]])
        
        temps.append(nodes[idx].cost)
                  
        idx = coordinatesToIndex(current_coord, 3)
    
    indexToCoordinates3D(idx, 3)
    current_point[0] = current_coord[0]
    current_point[1] = current_coord[1]
    current_point[2] = current_coord[2]
    path.append([current_point[0],current_point[1],current_point[2]])
    
    return path 

def Gradient_3D(START, GOAL,nodeslist,nodes, dim):
    
    current_point =  [0.0,0.0,0.0]
    path = []
    step = 1.0
    current = START
    idx= current
    current_coord= indexToCoordinates3D(idx, dim)
    
    current_point[0] = current_coord[0]
    current_point[1] = current_coord[1]
    current_point[2] = current_coord[2]
    path_velocity = []
    path.append([current_point[0],current_point[1],current_point[2]])
    
    path_velocity.append(nodeslist[nodes[idx].block][nodes[idx].idx])
    grads = [0.0,0.0,0.0]
    
    print 'OFPS Execution ...', current_coord, current_point, idx, nodes[idx].TAG
        
    while (idx !=GOAL):
        block = nodes[idx].block    
        grads[0] = -nodeslist[nodes[idx-1].block][nodes[idx-1].idx].cost/2+nodeslist[nodes[idx+1].block][nodes[idx+1].idx].cost/2 

        if (isinf(grads[0]) and grads[0]<0):
            grads[0] = -1
        elif (isinf(grads[0]) and grads[0]>0):
            grads[0] = 1
        
        max_grad = grads[0]
        
        for i in range(1,dim):
            grads[i]= (nodeslist[nodes[idx+d_[i-1]].block][nodes[idx+d_[i-1]].idx].cost/2)-(nodeslist[nodes[idx-d_[i-1]].block][nodes[idx-d_[i-1]].idx].cost/2)  
            
            if (isinf(grads[i]) and grads[i]<0):
                grads[i] = -1
            elif (isinf(grads[i]) and grads[i]>0):
                grads[i] = 1
            
            if abs(max_grad)<abs(grads[i]):
                max_grad = grads[i]
            print max_grad, grads[i]
        
        for i in range(dim):
            current_point[i] =((current_point[i] - step*grads[i]/abs(max_grad)))
            current_coord[i] = int(current_point[i]+0.5)
        
        path.append([round(current_point[0],3),round(current_point[1],3),round(current_point[2],3)])
        path_velocity.append(nodeslist[nodes[idx].block][nodes[idx].idx].v)
        idx = coordinatesToIndex(current_coord, dim)        
        print idx, current_coord, current_point
       
    current_coord=indexToCoordinates3D(idx, dim)
    current_point[0] = current_coord[0]
    current_point[1] = current_coord[1]
    current_point[2] = current_coord[2]
    path.append([round(current_point[0],3),round(current_point[1],3),round(current_point[2],3)])
    path_velocity.append(nodeslist[nodes[idx].idx].v)
    print 'Done !'
    
    return path