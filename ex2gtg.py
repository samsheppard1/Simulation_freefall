import numpy as np
import matplotlib . pyplot as plt


def y_val_drag_func(t):     # Function for equation for height with drag
            
    val_y = y_0 - (m/k)*(np.log( np.cosh( np.sqrt((k*g/m)) * t)))
    return val_y

def vy_val_drag_func(t):   # Equation for velocity with drag
            
    val_v = - (np.sqrt(m*g/k)) * (np.tanh(np.sqrt((k*g/m)) * tvals[i]))
    return val_v

# Defining variables

g = 9.81 # m/s^2, acceleration 
MyInput = 0
C_d = 1.0 # Drag coefficient
A = 1 # m^2, cross section area estimation for human 
k = C_d * 1.2 * A / 2 # Drag constant where the 1.2 is p_0
m = 200 # Mass for free fall

# Creation of arrays for equations of motion using drag


# User navigation menu

while MyInput != 'q':
    print('')   
    print('Choice "a" is for free-falling object under constant gravity and constant drag factor')
    print('Choice "b" is using eulers method to calculate the motion of an object in free-falling')
    print('Choice "c" is using a variable altitude dependent atmospheric pressure in the free-fall')
    print('Choice "d" is finding the maximum mach number ')
    print('')
    MyInput = input ( 'Enter a choice, "a", "b", "c", "d" or "q" to quit: ')
    print('') 

    print('You entered the choice: ', MyInput)
    print('') 
    if MyInput == 'a' :
        print ('You have chosen part (a)' )

        NumPoints = 300 # Number of values in a array for equations of motion and time values
        
        # The creation of an array for equations of motion with drag
        
        yvals_drag = np.zeros(NumPoints)  
        vyvals_drag = np.zeros(NumPoints)
        y_0 = 1000 # m, initial height 
        
        tmin = 0.0 # seconds
        tmax = 20.0 # seconds
        tvals = np.linspace(tmin,tmax,NumPoints) 

        for i in range(NumPoints) : # iteration to add values to arrays for both equations of motion              
            
            yvals_drag[i] = y_val_drag_func(tvals[i]) 
           
            vyvals_drag[i] = vy_val_drag_func(tvals[i])

        fig , (ax3, ax4) = plt .subplots(1, 2, figsize=(12,4)) # Creation of both graphs
        fig . suptitle ('Height and velocity data for freefall with drag')
        ax3.set(xlabel = 'Time (s)', ylabel = 'Height (m)' , title= 'Altitude')
        ax4.set(xlabel='Time (s)', ylabel='y-velocity (m/s)', title='Vertical speed')
        ax3.plot(tvals ,yvals_drag,'tab:orange')
        ax4.plot(tvals ,vyvals_drag,'tab:blue')
        plt.show()
        
    elif MyInput == 'b':
        print ('You have chosen part (b)' )
        
        delta_t = 0.05 # Value for delta t 
        y_eul = [0]
        v_eul = [0]
        t_eul = [0]
        y_eul[0] = 1000
        l = len(y_eul) # Length of array 
         
        while y_eul[l-1] > 0: # Iterations to add values to arrays for free falll until it hits the ground
            l = len(y_eul)
            y_eul.append(float(0))
            v_eul.append(float(0))
            t_eul.append(float(0))
             
            t_eul[l] = t_eul[l-1] + delta_t
            v_eul[l] = v_eul[l-1] - delta_t*(g + (k/m)*v_eul[l-1]*abs(v_eul[l-1]))
            y_eul[l] = y_eul[l-1] + delta_t*v_eul[l-1]

        t_eul.pop(l-1)
        y_eul.pop(l-1)
        v_eul.pop(l-1)

        fig , (ax5, ax6) = plt .subplots(1, 2, figsize=(12,4)) # Creation of graph 
        fig . suptitle ('Height and velocity data fusing Eulers method')
        ax5.set(xlabel = 'Time (s)', ylabel = 'Height (m)' , title= 'Altitude')
        ax6.set(xlabel='Time (s)', ylabel='y-velocity (m/s)', title='Vertical speed')
        ax5.plot(t_eul ,y_eul,'tab:green')
        ax6.plot(t_eul ,v_eul,'tab:purple')
        plt.show()
    
    elif MyInput == 'c':
        print ('You have chosen part (b)' )
 
        y_0 = 1000 # m, initial height 
        NumPoints = 300 # Number of values in arrays for equation of motions and time
        tmin = 0.0 # seconds
        tmax = 20 # seconds
        tvals = np.linspace(tmin,tmax,NumPoints) 
        
        p_0 = 1.2
        h = 7640 #

        #e Creation of arrays for equation of motions 
        yvals_drag = np.zeros(NumPoints) 
        vyvals_drag = np.zeros(NumPoints)
        
        yvals_drag[0] = 1000 # First array value for height
        
        p_val = np.zeros(NumPoints) # Creation of array for air density values
        
        
        p_val[0] = 1.2 # First array value for density
        
        k_val = np.zeros(NumPoints) # Creation of array for drag constants
        k_val[0] = (C_d*1.2*A)/2 # First array value for drag constants

        for i in range(1,NumPoints) : # Iteration to add values to array of motions and air density and drag constants
             
            yvals_drag[i] = y_0 - (m/k_val[i-1])*(np.log( np.cosh( np.sqrt((k_val[i-1]*g/m)) * tvals[i])))
            
            vyvals_drag[i] = - (np.sqrt(m*g/k_val[i-1])) * (np.tanh(np.sqrt((k_val[i-1]*g/m)) * tvals[i]))
            
            p_val[i] = p_0 * np.exp(yvals_drag[i]/h)
            
            k_val[i] = (C_d*p_val[i]*A)/2
            
        fig , (ax7, ax8) = plt .subplots(1, 2, figsize=(12,4)) # Creation of graphs
        fig . suptitle ('Height and velocity data for using air density equation')
        ax7.set(xlabel = 'Time (s)', ylabel = 'Height (m)' , title= 'Altitude')
        ax8.set(xlabel='Time (s)', ylabel='y-velocity (m/s)', title='Vertical speed')
        ax7.plot(tvals ,yvals_drag,'tab:green')
        ax8.plot(tvals ,vyvals_drag,'tab:purple')
        plt.show()
        
    elif MyInput == 'd': 
        
        
        NumPoints = 400 # Number of values in a array for equations of motion and time values
        
        # Creation of arrays for equation of motion
        
        yvals_drag = np.zeros(NumPoints)  
        vyvals_drag = np.zeros(NumPoints)
        
        tmin = 0.0 # seconds
        tmax = 50.0 # seconds
        tvals = np.linspace(tmin,tmax,NumPoints) 
        y_0 = 1000 # Starting fall height, m
    
        def T(H) : # funciton for absolute temperature, kelvin
            if H <= 11000:
                T_output = 288 - 0.0065*H
            elif 11000 < H <= 25100:
                T_output = 216.5
            elif H > 25100:
                T_output = 141.3 + 0.003*H
            return T_output
        
        print ('You have chosen part (d)' )
        
        mol_mass = 0.0289645 # Molar mass, kg/mol
        R = 8.31 # Molar gas constant 
        v_s = np.zeros(NumPoints) # creation of array for speed of sound in air 
        
        Bau_max = np.zeros(NumPoints) # creation of array for the ratio of velocity and speed of sound in air 
        
        for i in range(NumPoints) : # Iteration to add values to each array 
           
           yvals_drag[i] = y_val_drag_func(tvals[i])
           vyvals_drag[i] = vy_val_drag_func(tvals[i])
           v_s[i] = np.sqrt((1.4*R*T(yvals_drag[i]))/mol_mass)
           Bau_max[i] = (vyvals_drag[i] / v_s[i])
        
        plt.plot(vyvals_drag, Bau_max) # Graph creation 
        plt.xlabel('y-velocity (m/s)')
        plt.ylabel('Baumgartners Mach Number')
        plt.title('Baumgartners Mach Number against velocity')
        plt.show()
           
           
        
        
        
    elif MyInput != 'q':
        print('This is not a valid choice')
        print('You have chosen to finish - goodbye.')