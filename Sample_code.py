## Sample code of calculating grating coupler performance. The rest of the application is merely standard Genetic algorithm
## and Tkinter for GUI

from camfr import *
from numpy import *

angl = cos(theta_in*3.141592654/180)-sin(theta_in*3.141592654/180)

if GC_polarization == 0: # input from Tkinter to determine the chosen polarization
    Tem = TE
elif GC_polarization == 1:
    Tem = TM
if GC_structure == 0: # input from Tkinter to determine straight or apodized grating coupler
    # define parameters
    set_lambda(lambda0) # wavelength
    set_N(N_mode)      #number of modes
    set_polarisation(Tem)
    set_chunk_tracing(0)
    set_degenerate(0)
    set_orthogonal(False)

    # define material
    substrate=Material(sub_n)
    guiding=Material(dev_n)
    iml=Material(sup_n)

    # create output file
    outfile = open(save_file + ".txt",'w')

    # define our own parameters
    ff = fill
    d = Air_thickness
    pml = PML_strength
    set_lower_PML(-pml)
    set_upper_PML(-pml)
    guide_thickness = dev_t
    dclad = t_sub
    period = Period
    groove_depth = etch

    # sweep the wavelength
    for P in arange(lambda1*1000,lambda2*1000+1,2):
        set_lambda(P/1000.0)
        # define slabs
        waveguide = Slab(guiding(d-dclad)+substrate(dclad)+guiding(guide_thickness)+
                         iml(d))
        etched = Slab(guiding(d-dclad)+substrate(dclad)+
                      guiding(guide_thickness-groove_depth)+iml(groove_depth+d))
    
        # define stack
        stack = Stack(waveguide(1.0) + 20*(etched(period*(1.0-ff)) + waveguide(period*ff)) + waveguide(2.0))

        # find the guided mode
        waveguide.calc()
        guided = 0
        niguided = 1
        for t in range(0,60):
            if abs(waveguide.mode(t).n_eff().imag) < niguided:
                guided = t
                niguided = abs(waveguide.mode(t).n_eff().imag)

        # set input for calculating the fields
        inc = zeros(N())
        inc[guided] = 1
        stack.set_inc_field(inc)

        stack.calc()

        R = abs(stack.R12(guided,guided))
        T = abs(stack.T12(guided,guided))
        up = stack.lateral_S_flux(d+1.5)
        down = stack.lateral_S_flux(d-1.5)
        x = d+1.5
        xd = d-1.5

        #calculate the coupling efficiency to fibre
        powerup = 0.0 + 0.0*1j
        powerdown = 0.0 + 0.0*1j
        soverlapint = zeros(100,complex)  # 8 degrees
        align = zeros(100)  # to find optimal fiber position
        pfib = 0.0

        #normalization gausian profile
        Zzero = 377
        nZ = sup_n/Zzero
        for z in arange(-10.0,10.0,0.01):
            pfib+=0.01*exp(-((z/5.2)**2))*nZ*exp(-((z/5.2)**2))

        #calculate different fibre position simultaniously
        for counter in range(100):
            align[counter]=-1.0-period*(counter/5)
        for z in arange(0.01, 23*period , 0.01):
            if Tem == TE:
                powerup+=(0.01*stack.field(Coord(x,0,z)).E2()*
                            conjugate(stack.field(Coord(x,0,z)).Hz()))
                powerdown+=(0.01*stack.field(Coord(xd,0,z)).E2()*
                            conjugate(stack.field(Coord(xd,0,z)).Hz()))
                soverlapint+=(0.01*stack.field(Coord(x,0,z)).E2()*
                                nZ*exp(-(((align+z)/5.2)**2))*exp(angl*1j*z))
            elif Tem == TM :
                powerup+=(0.01*stack.field(Coord(x,0,z)).H2()*
                            conjugate(stack.field(Coord(x,0,z)).Ez()))
                powerdown+=(0.01*stack.field(Coord(xd,0,z)).H2()*
                            conjugate(stack.field(Coord(xd,0,z)).Ez()))
                soverlapint+=(0.01*stack.field(Coord(x,0,z)).H2()*
                                exp(-(((align+z)/5.2)**2))*exp(angl*1j*z))
            
        coupling = ((abs(overlapint))**2)/(pfib*powerup.real)
        scoupling = ((abs(soverlapint))**2)/(pfib*powerup.real)

        if GC_type == 1:
            Radius = 3.141592654*width0*waveguide.mode(t).n_eff().real*5.2/lambda0
            sectionAngle = 2*atan(5.2/Radius)
        elif GC_type == 0:
            Radius = 'infinity'
            sectionAngle = 'zero'

        print >> outfile, P,ff,groove_depth,dclad,R*R,T*T,powerup.real, powerdown.real,powerup.real*max(scoupling), Radius, sectionAngle
        outfile.flush()
        print 'outpt file row sequence: 1-Period 2-fill factor 3-etch depth 4-substrte thickness\
5-Reflection 6-Transmission 7-Power up 8-Power down 9-Coupling 10-Curved coupler radius 11- Curved coupler section angle.'
        free_tmps()
    outfile.close()
