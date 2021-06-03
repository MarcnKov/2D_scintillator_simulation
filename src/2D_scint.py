import numpy as np
import math
from random import random
import matplotlib.pyplot as plt

def count_exit_photons(exitX):
    gt0, lt0 = 0,0
    for x in exitX:
        if x >= 0:
            gt0 += 1
        elif x < 0:
            lt0 += 1
    return gt0, lt0

def draw_angle():

    theta = np.random.uniform(0, 2*np.pi)
    if (theta == 0):
        return np.random.uniform(0, 2*np.pi)
    elif (theta == np.pi/2):
        return np.random.uniform(0, 2*np.pi)
    elif (theta == np.pi):
        return np.random.uniform(0, 2*np.pi)
    elif (theta == 1.5*np.pi):
        return np.random.uniform(0, 2*np.pi)
    else:
        return theta

def draw_len():
    
    l0 = 150
    return np.random.exponential(l0)

def draw_lambertian_angle():
    
    u = np.random.uniform(-1,1)
    return np.arcsin(u)

def compute_dist(XY1, XY2):
    return np.sqrt( (XY2[0]-XY1[0])**2 \
                +   (XY2[1]-XY1[1])**2 )

def snells_descartes_law(incident_angle):
    
    n = 1.357 #n1/n2 = 1.9/1.4
    if (np.abs(n*np.sin(incident_angle)) <= 1): 
        return True
    else:
        return False

def compute_next_collision(XY1,velXY,segments):
    '''
    La fonction calcule le temps minimum de
    collision
    '''
    coll_t_seg = {} 
    t_tmp = 1

    #parcourir tous les segments
    for seg in segments:
        XY2 = seg.coordinates()
        
        #calculer le temps selon X ou Y
        if (XY2[0] == 0):
            t = (XY2[1]-XY1[1])/velXY[1]
        else:
            t = (XY2[0]-XY1[0])/velXY[0]
        
        #choisir le temps minimum parmi quatres segments 
        if (t > 0 and t < t_tmp and t > 1e-20):
            coll_t_seg[0], coll_t_seg[1] = t, seg
            t_tmp = t

    #renvoyer le temps minimum et le segments correspondant à temps minimum 
    return coll_t_seg[0], coll_t_seg[1]

def calculate_velocity(theta):
    v0 = 3e11
    velXY = [0,0]
    velXY = [v0*np.cos(theta),v0*np.sin(theta)]

    return velXY

def compute_new_velocity(posXY,velXY,seg_ID,theta):
    '''
    La fonction calcule la nouvelle vitesse
    '''
    velXY_new = calculate_velocity(theta)
    
    if (velXY[0] < 0 and velXY_new[0] > 0):
        velXY[0] = -velXY_new[0]
    elif (velXY[1] < 0 and velXY_new[1] > 0):
        velXY[1] = -velXY_new[1]
    
    #photon se trouve sur les bords
    if (np.abs(round(posXY[0])) == 25.0 and np.abs(round(posXY[1])) == 5.0):
        velXY[0] = -velXY[0]
        velXY[1] = -velXY[1]
    #changer la composante (du vecteur vitesse) orthogonale à la face de collision
    elif (seg_ID == "AB"):
        velXY[0] = -velXY[0]
    elif (seg_ID == "BC"):
        velXY[1] = -velXY[1]
    elif (seg_ID == "CD"):
        velXY[0] = -velXY[0]
    elif (seg_ID == "DA"):
        velXY[1] = -velXY[1]
    
    return velXY[0], velXY[1] 

def compute_new_position(posXY,velXY,t):
    
    posXY_new = [0,0]
    posXY_new[0] = posXY[0] + velXY[0]*t  
    posXY_new[1] = posXY[1] + velXY[1]*t
    
    return posXY_new[0], posXY_new[1]

def compute_collision_type(coll_seg,posXY,velXY): 
    
    #calculer l'angle de collision selon x
    theta = np.arctan2(velXY[1],velXY[0])
    
    isAlive = True
    isAbsorbed = False
    
    if (coll_seg.absorbing() == True):
        isAlive = False
        isAbsorbed = True

    elif(coll_seg.exit_seg() == True):
        
        #projeter selon y
        if theta < 0:
            theta = np.pi/2 + theta
        else:
            theta = np.pi/2 - theta 
        
        #vérifier si la transmission est possible 
        theta_trans = snells_descartes_law(theta)
        if(theta_trans == True):
            isAlive = False
            isAbsorbed = False
    
    elif(coll_seg.lambertian() == True):

        theta = draw_lambertian_angle()
        if (theta > 0):
            theta = theta + np.pi/2
        else:
            theta = theta - np.pi/2
    
    return theta,isAlive,isAbsorbed

class Segment:
    '''
    Objet segment dont les faces ont les différentes
    propriétés
    '''
    coord = []
    def __init__(self,name,absorb,lamb,exit):

        self.name = name
        self.absorb = absorb
        self.lamb = lamb
        self.exit = exit
        
        if self.name == 'AB':
            self.coord = [25,0]
        elif self.name == 'BC':
            self.coord = [0,-5]
        elif self.name == 'CD':
            self.coord = [-25,0]
        elif self.name == 'DA':
            self.coord = [0,5]
        
    def coordinates(self):
        return self.coord
    def seg_name(self):
        return self.name
    def absorbing(self):
        return self.absorb
    def lambertian(self):
        return self.lamb
    def exit_seg(self):
        return self.exit
    

def begin_simulation(N,AB_abs,CD_abs,DA_lamb):
    
    #Initializer les segments
    AB = Segment('AB',AB_abs,False,False)
    BC = Segment('BC',False,False,True)
    CD = Segment('CD',CD_abs,False,False)
    DA = Segment('DA',False,DA_lamb,False)
    
    segments = [AB,BC,CD,DA]
    
    exit_x_coord = []
    for k in range(N):
        
        theta = draw_angle() 
        dist_max = draw_len()
        posXY = [0,0]
        posXY_new = [0,0]
        velXY = calculate_velocity(theta)
            
        dist_traveled = 0
        while (dist_traveled < dist_max):
            
            t_new, coll_seg = compute_next_collision(posXY,velXY,segments)
            posXY_new[0], posXY_new[1] = compute_new_position(posXY,velXY,t_new)
            theta, isAlive, isAbsorbed = compute_collision_type(coll_seg,posXY_new,velXY)  
            
            if(not(isAlive) and isAbsorbed):
                break
            elif(not(isAlive) and not(isAbsorbed)):
                exit_x_coord.append(posXY_new[0])
                break

            dist_traveled += compute_dist(posXY,posXY_new)
            velXY[0], velXY[1] = compute_new_velocity(posXY_new,velXY,coll_seg.seg_name(),theta)
            posXY = posXY_new.copy()
    
    return exit_x_coord


'''
          D  LAMB/SPEC   A
            ------------
ABS/SPEC    |          |    ABS/SPEC
            |          | 
            ------------            
          C    SORTIE    B        
'''


#L'histograme pour la densité linéique de lumière, modifier les différents faces

exitX_per_simul = []
num_photons = 10000 
AB_abs = False
CD_abs = False
DA_lamb = False
exitX_per_simul.append(begin_simulation(num_photons,AB_abs,CD_abs,DA_lamb))

gt0, lt0 = count_exit_photons(exitX_per_simul[0])
print("% exit : ", len(exitX_per_simul[0])/num_photons)
print("[0,25] = ", gt0/num_photons)
print("[-25,0] = ",lt0/num_photons) 

kwargs = dict(alpha=0.5, bins=60, density=True, stacked=True)
plt.hist(exitX_per_simul,**kwargs)
plt.title('Densité linéique de lumière', fontsize = 20)
plt.ylabel('Densité de lumière', fontsize = 20)
plt.xlabel('Coordonnées (x, -5) en (mm) sur la face de sortie', fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.show()

#L'histograme pour la quantité de lumière
#choisir num_simulations = N, modifier les differentes faces

'''
num_simulations = 1000
params = [[False,False,False],[True,False,False],[True,True,False],[False,False,True]]

for par in params:

    num_photons_exit = []
    exitX_per_simul = []
    for i in range(num_simulations):
            
        num_photons = 10000
        exitX_per_simul.append(begin_simulation(num_photons,par[0],par[1],par[2]))
        
        num_photons_exit.append(len(exitX_per_simul[i]))
    
    print(str(par[0])+ str(par[1])+ str(par[2]))
    print("Mean = ", np.mean(num_photons_exit))
    print("Std = ", np.std(num_photons_exit))
    print("Median = ", np.median(num_photons_exit))
    print("Max-Min = ",np.max(num_photons_exit)-np.min(num_photons_exit))
 

    w = 6.5
    num_photons_exit = np.array(num_photons_exit)
    n = math.ceil((num_photons_exit.max()-num_photons_exit.min())/w)
    kwargs = dict(alpha=0.5, bins=n, density=True, stacked=True)
    plt.hist(num_photons_exit,**kwargs)
    plt.title('Quantité de lumière', fontsize = 20)
    plt.ylabel('Quantité de lumière', fontsize = 20)
    plt.xlabel('Nombre des photons sortis', fontsize = 20)
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.show()
'''
