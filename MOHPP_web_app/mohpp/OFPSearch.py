'''
Created on Jul 8, 2021

@author: nassim
'''
from math import floor, isinf, sqrt
from utilities import indexToCoordinates, coordinatesToIndex


def computeGradientDescent(s_idx, g_idx, nodes, d_):
    print '--> OFPSearch processing...'
    nbrdim, cur_point, cur_coord, step, idx, gradient, path = 2,[0.0, 0.0], indexToCoordinates(s_idx, d_),\
                                                            1.0, s_idx, [0.0, 0.0], []
    
    cur_point[0],cur_point[1] = float(cur_coord[0]), float(cur_coord[1])
    path.append([cur_point[0],cur_point[1]])
    
    while(nodes[idx].cost != 0):
        
        gradient[0] = -nodes[idx-1].cost /2 + nodes[idx+1].cost /2 
        
        if (isinf(gradient[0]) and gradient[0]<0):
            gradient[0] = -1
        
        elif (isinf(gradient[0]) and gradient[0]>0):
            gradient[0] = 1
        
        max_grad = gradient[0]
        
        for i in range(1, nbrdim):
            
            gradient[i] = -nodes[idx-d_[i-1]].cost /2 + nodes[idx+d_[i-1]].cost /2   
            
            if (isinf(gradient[i]) and gradient[i]<0):
                gradient[i] = -1
            
            elif (isinf(gradient[i]) and gradient[i]>0):
                gradient[i] = 1
            
            if abs(max_grad)<abs(gradient[i]):
                max_grad = gradient[i]
                
        for i in range(nbrdim):
            cur_point[i] =((cur_point[i] - step*gradient[i]/abs(max_grad)))
            cur_coord[i] = int(cur_point[i]+0.5)
        
        path.append([round(cur_point[0],5),round(cur_point[1],5)])
        idx = coordinatesToIndex(cur_coord, d_)            
    
    cur_coord=indexToCoordinates(idx, d_)
    cur_point[0] = float(cur_coord[0])
    cur_point[1] = float(cur_coord[1])
   
    path.append([round(cur_point[0],5),round(cur_point[1],5)])    
    print '                    path found !'
    return path


def Gradient(START, GOAL, nodes, d_):


    temps = []
    path = []
    current_point=[0.0,0.0]
    idx = START   
    grads = [0.0,0.0]
    current_coord = indexToCoordinates(idx, d_)
    current_point[0] = float(current_coord[0])
    current_point[1] = float(current_coord[1])
    
    path.append([current_point[0],current_point[1]])
    temps.append(nodes[idx].cost)
    
   
    while (nodes[idx].cost !=0):
        idx_x_l = idx - 1
        idx_x_u = idx+1
    
        index1 = (max(idx_x_l,0))
        index2 = (max(idx_x_u,0))
    
        grad = 0.0
        cnt = 9
    
        for i in range(-1,2):
        
            for j in range(-1,2):
                index_n_l = index1 + i*int(d_[0])+j*int(d_[1])
                index_s_l = min((max(index_n_l,0)),d_[1]-1)
                value_l = nodes[index_s_l].cost
                
                index_n_u = index2 + i*int(d_[0])+j*int(d_[1])
                index_s_u = min((max(index_n_u,0)),d_[1]-1)
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
                
                index_n_l = index1 + i*1 + j*d_[1]
                index_s_l = min((max(index_n_l,0)),d_[1]-1)
                
                value_l  = nodes[index_s_l].cost
                
                
                index_n_u = index2 + i*1 + j*d_[1]
                index_s_u = min((max(index_n_u,0)),d_[1]-1)
                
                value_u  = nodes[index_s_u].cost
    
    
                if isinf(value_l) or isinf(value_u):
                    cnt-=1
            
                else : 
                    grad +=(value_u-value_l)/2.0
        grad/=cnt
        grads[1] = grad
    
        max_dif = -1000.0
        best_idx = idx
        
        if isinf(grads[0]) or isinf(grads[1]) or grads[0]==0.0 or grads[1]==0.0:
            
            for i in range(-1,2):
                for j in range(-1,2):
                   
                    if(i==0 and j==0):
                        continue
                    """
                    idx_delta = i+d_[0]*j
                    idx_n = idx+idx_delta
                    
                    idx_tmp = (idx_n)
                    """  
                    idx_tmp = min((max(idx+i+d_[0]*j,0)),d_[1]-1)
                    
                    if isinf(nodes[idx_tmp].cost):
                        continue
                    dif = (-nodes[idx_tmp].cost + nodes[idx].cost)/sqrt(i*i+j*j)
                    if dif > max_dif:
                        max_dif = dif
                        best_idx = idx_tmp
                
            best_idx = max(best_idx,0)
            current_coord = indexToCoordinates(best_idx, d_)
            current_point[0] = float(current_coord[0])
            current_point[1] = float(current_coord[1])
        else:
            
            grad_norm = sqrt(grads[0]*grads[0]+grads[1]*grads[1])
            
            for i in range(2):
                current_point[i] = current_point[i] -grads[i]/grad_norm
                current_coord[i] = int(current_point[i]+0.5)
        
        
       
        path.append([current_point[0],current_point[1]])
        
        temps.append(nodes[idx].cost)
                  
        idx = coordinatesToIndex(current_coord, d_)
    
    indexToCoordinates(idx, d_)
    current_point[0] = current_coord[0]
    current_point[1] = current_coord[1]
    path.append([current_point[0],current_point[1]])
    
    return path 