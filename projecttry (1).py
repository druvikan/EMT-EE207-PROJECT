import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
while (True):
    choice = int(input(
        "1.Simulation of cyclotron\n2.Variation with particle.\n3.Variation with Magnetic Field\n4.Exit\nEnter your choice number: "))
    if choice == 1:
        a = int(input("Choose a particle\n1. Deutron\n2. Proton\n3. alpha particle\n"))
        if a == 1:
            # Deutron
            q = 1.6e-19
            m = 2*1.67e-27
        elif a == 2:
            # Proton
            q = 1.6e-19
            m = 1.67e-27
        elif a == 3:
            # Alpha Particle
            q = 2*1.6e-19
            m = 4*1.67e-27
        else:
            print("Enter correct choice")

        V = float(input("Enter Voltage: "))  # Voltage between plates
        d = 90e-6  # separation between the plates
        E_0 = V/(d)  # Electric Field

        b = float(input("Enter magnetic field in positive z direction: "))
        # Set magnetic field to bT in the +Z direction
        B = np.array([0.0, 0.0, b])
        r_cyclotron = .05  # set the radius of the D's to 5cm
        particlepos = np.array([0.0, 0.0, 0.0])
        # Set the initial particle speed to 0
        particlev = np.array([0.0, 0.0, 0.0])

        # create an array for the x position of the particle
        particleposx = [particlepos[0]]
        # create an array for the y position of the particle
        particleposy = [particlepos[1]]

        w = q*np.linalg.norm(B)/m  # define the cyclotron frequency

        t = 0  # initialize time to 0
        dt = 5e-12  # Set timestep to 5 picoseconds

        # loop while the magnitude of the proton's position remains within the cyclotron radius
        while (np.linalg.norm(particlepos) < r_cyclotron):

            # create a vector for the net force on the particle.
            Fnet = np.array([0.0, 0.0, 0.0])

            # if the particle is between the two D's calculate the electric force
            if np.abs(particlepos[0]) < d/2:
                Fnet[0] = q*E_0*np.cos(w*t)
            else:  # if the particle is not, calculate the magnetic force
                Fnet = q*np.cross(particlev, B)

            particlev = particlev + Fnet*dt/m  # Update the velocity of the particle
            # Use velocity to update the position of the particle
            particlepos = particlepos + particlev*dt

            # append the x position to the x-position list
            particleposx = np.append(particleposx, particlepos[0])
            # append the y position to the y-position list
            particleposy = np.append(particleposy, particlepos[1])
            t = t + dt  # update the timestep

        print("The final speed of the particle is",
              np.linalg.norm(particlev), "m/s")

        plt.figure(figsize=(12, 12))  # create the figure
        plt.plot(particleposx, particleposy)  # create the plot
           
        plt.show()

    elif choice == 2:
        V = float(input("Enter Voltage: "))
        d = 90e-6
        E_0 = V/(d)
        b = float(input("Enter magnetic field in positive z direction: "))
        # Set magnetic field to bT in the +Z direction
        B = np.array([0.0, 0.0, b])
        r_cyclotron = .05  # set the radius of the D's to 5cm
        particlepos = np.array([0.0, 0.0, 0.0])
        # Set the initial particle speed to 0
        particlev = np.array([0.0, 0.0, 0.0])

        # Store the final velocities
        final_vel = np.array([0.0, 0.0, 0.0])

        # Store the value of q and m in a contiguous fashion
        qs = np.array([1.6e-19, 1.6e-19, 2*1.6e-19])
        ms = np.array([2*1.67e-27, 1.6e-27, 4*1.6e-27])
        t = 0  # initialize time to 0
        dt = 5e-12  # Set timestep to 5 picoseconds

        # Iterate the simulation for each particle by varying q and m
        for i in range(0, 3):
            q = qs[i]
            m = ms[i]
            w = q*np.linalg.norm(B)/m  # define the cyclotron frequency
            # loop while the magnitude of the proton's position remains within the cyclotron radius
            while (np.linalg.norm(particlepos) < r_cyclotron):

                # create a vector for the net force on the particle.
                Fnet = np.array([0.0, 0.0, 0.0])

                # if the particle is between the two D's calculate the electric force
                if np.abs(particlepos[0]) < d/2:
                    Fnet[0] = q*E_0*np.cos(w*t)
                else:  # if the particle is not, calculate the magnetic force
                    Fnet = q*np.cross(particlev, B)

                particlev = particlev + Fnet*dt/m  # Update the velocity of the particle
                # Use velocity to update the position of the particle
                particlepos = particlepos + particlev*dt
                t = t + dt  # update the timestep

            # Update the final vel and reset position and velocity
            final_vel[i] = np.linalg.norm(particlev)
            particlepos = np.array([0.0, 0.0, 0.0])
            particlev = np.array([0.0, 0.0, 0.0])

        print("Final speeds are: ")
        print("Particle \t Final Speed")
        print("-------------------------------------------------------")
        print("Deutron \t", final_vel[0])
        print("Proton \t\t", final_vel[1])
        print("Alpha \t\t", final_vel[2])

    elif choice == 3:
        a = int(input("Choose a particle\n1. Deutron\n2. Proton\n3. alpha particle\n"))
        if a == 1:
            q = 1.6e-19  # Set the charge of the particle to the charge of a proton
            m = 2*1.67e-27  # Set mass of the particle to the mass of a proton
        elif a == 2:
            q = 1.6e-19  # Set the charge of the particle to the charge of a proton
            m = 1.67e-27  # Set mass of the particle to the mass of a proton
        elif a == 3:
            q = 2*1.6e-19  # Set the charge of the particle to the charge of a proton
            m = 4*1.67e-27  # Set mass of the particle to the mass of a proton
        else:
            print("Enter correct choice")
        V = float(input("Enter Voltage: "))
        d = 90e-6
        E_0 = V/(d)
        r_cyclotron = .05
        final_vel = np.array([0.0, 0.0, 0.0])
        # Stores the magnetic fields
        magneticfields = np.array([0.0, 0.0, 0.0])

        # Input the magentic fields:
        for i in range(0, 3):
            magneticfields[i] = float(
                input("Enter magnetic field in positive z direction: "))

        # Iterate the simulation for each of the magnetic fields.
        for i in range(0, 3):
            B = np.array([0.0, 0.0, magneticfields[i]])
            particlepos = np.array([0.0, 0.0, 0.0])
            # Set the initial particle speed to 0
            particlev = np.array([0.0, 0.0, 0.0])
            w = q*np.linalg.norm(B)/m  # define the cyclotron frequency
            t = 0  # initialize time to 0
            dt = 5e-12  # Set timestep to 5 picoseconds

            # loop while the magnitude of the proton's position remains within the cyclotron radius
            while (np.linalg.norm(particlepos) < r_cyclotron):

                # create a vector for the net force on the particle.
                Fnet = np.array([0.0, 0.0, 0.0])

                # if the particle is between the two D's calculate the electric force
                if np.abs(particlepos[0]) < d/2:
                    Fnet[0] = q*E_0*np.cos(w*t)
                else:  # if the particle is not, calculate the magnetic force
                    Fnet = q*np.cross(particlev, B)

                particlev = particlev + Fnet*dt/m  # Update the velocity of the particle
                # Use velocity to update the position of the particle
                particlepos = particlepos + particlev*dt
                t = t + dt  # update the timestep
            final_vel[i] = np.linalg.norm(particlev)

        print("Final speeds are: ")
        print("Particle \t Final Speed")
        print("-------------------------------------------------------")
        print(magneticfields[0], "\t", final_vel[0])
        print(magneticfields[1], "\t", final_vel[1])
        print(magneticfields[2], "\t", final_vel[2])

    elif choice == 4:
        print("......................................................program terminated.....................................................")
        break
    else:
        print("Enter a valid choice!")
