from setup          import *
from cp0_analytic   import *

class Aggregate:
    def __init__(self,num,L,b,cm0,cp0):
        # {'Micellar':3,'Hexagonal':2,'Lamellar':1}
        self.Dim    = num
        self.L      = L                                      # size of aqueous region,         [A]^3
        self.b      = b                                      # radius of surfactant aggregate charge surface,  [A]^3
        self.cm0    = cm0                                    # minus ion concentration at cell boundary, [mM]
        self.cp0    = cp0                                    # plus  ion concentration at cell boundary, [mM]

        self.v1     = v1                                     # monomer volume  v1 = V_tail+V_head
        self.v0     = v0                                     # Tanford hydrophobic tail volume
        self.M_amph = M_amph                                 # amphiphile molecule mass, including counterion K
        self.l0_exd = l_tail_exd                             # Tanford hydrophobic tail length extend

    def Dependants(self):
        d           = self.Dim                               # dimension of aggregate
        self.B      = self.b+self.L                          # size of cell
        self.x      = np.linspace(self.b/self.B,1.0,nx)      # x is for the PBE_Solver, Dimensionless

        self.l0     = (self.v0/self.v1)**(1/d)*self.b        # tail length in the hydrocarbon core

        self.ba     = self.b                                 # aggregate   size
        self.bc     = self.b                                 # charge      size
        self.bt     = self.b                                 # tension     size
        self.bh     = self.l0                                # hydrophobic size

        self.Aa     = 2**(d-1)*pi/1*self.ba**(d-1)*B0**(3-d) # aggregate   surface area
        self.Ac     = 2**(d-1)*pi/1*self.bc**(d-1)*B0**(3-d) # charge      surface area
        self.At     = 2**(d-1)*pi/1*self.bt**(d-1)*B0**(3-d) # tention     surface area
        self.Ah     = 2**(d-1)*pi/1*self.bh**(d-1)*B0**(3-d) # hydrophobic surface area

        self.Vagga  = 2**(d-1)*pi/d*self.ba**(d-0)*B0**(3-d) # aggregate   volume
        self.Vchag  = 2**(d-1)*pi/d*self.bc**(d-0)*B0**(3-d) # charge      volume
        self.Vtens  = 2**(d-1)*pi/d*self.bt**(d-0)*B0**(3-d) # tension     volume
        self.Vhpbc  = 2**(d-1)*pi/d*self.bh**(d-0)*B0**(3-d) # hydrophobic volume

        self.Vcell  = 2**(d-1)*pi/d*self.B**(d-0)*B0**(3-d)  # cell        volume
        self.Vaque  = self.Vcell-self.Vagga                  # Aqueous

        self.n      = self.Vagga/self.v1                     # aggregation number

        self.sigma  = -1.0*e*self.n/self.Ac*1e20             # [sigma] = C/m^2
        self.Cm     = 1.0/self.Vcell*1e30/NA                 # Micell, aggregate concentration, mM

    def Calculation(self,out_opt=0):
        """
        Set up for  PB equaiton solving
        generate the optimized b, cm0, cp0
        equilibrate the amphiphile ion in aqueous and in aggregates
        """
        # update the variables,and localized into the Calculation function

        d       = self.Dim
        L       = self.L
        b       = self.b
        cm0     = self.cm0                      # minus ion concentration at cell boundary, [mM]
        cp0     = self.cp0                      # plus  ion concentration at cell boundary, [mM]

        self.Dependants()

        bt      = self.bt
        B       = self.B

        Aa      = self.Aa
        Ac      = self.Ac
        At      = self.At
        Ah      = self.Ah

        v1      = self.v1
        v0      = self.v0
        l0      = self.l0
        l0_exd  = self.l0_exd
        
        n       = self.n
        sigma   = self.sigma
        x       = self.x

        PBsolution  = self.PBE_Solver()
        Ui          = self.Ui_Cal(x)

        y       = PBsolution[:,0]               # reduced potential y = e*Phi/kT,  dimensionless
        y_prime = PBsolution[:,1]               # derivative of reduced potential, dimensionless
        yp      = y_prime[-1]

        print('y[-1]  = %8.4f y[0]  = %8.4f'%(y[-1],y[0]))
        print('Ui[-1] = %8.4f Ui[0] = %8.4f'%(Ui[-1],Ui[0]))
        Im      = (Ac*b)*cumtrapz(exp( y)/x**(1+d),   x,initial=0)
        Ip      = (Ac*b)*cumtrapz(exp(-y-Ui)/x**(1+d),x,initial=0)
        Nm      = self.cm0*NA*1e-30*Im[-1]
        Np      = self.cp0*NA*1e-30*Ip[-1]
        Eel     = 0.5*D*eps0*(Ac/b)*cumtrapz(x**(3-d)*(y_prime*kT/e)**2,x,initial=0)*1e-10  # [J]

        mu10 =  @mu10
        g1   =  @g1

        mu20 =  @mu20
        g0   =  @g0
        H0   =  @H0
        kc   =  (g1-g0)/(2*H0**2)

        mu30 = @mu30
        kg   = kc*(H0*bt-2)

        if d == 1:
            c1          = 0.0
            c2          = 0.0
            H           = 0.5*(c1+c2)
            KG          = c1*c2
            gamma       = g1
            Gob_curv    = 0.0
            mu0_correct = mu10
        
        if d == 2:
            c1          = 1/bt
            c2          = 0.0
            H           = 0.5*(c1+c2)
            KG          = c1*c2
            gamma       = g0+2*kc*(H-H0)**2
            Gob_curv    = At*(2*kc*(H-H0)*(-1)/bt)*1e-23
            mu0_correct = mu20

        if d == 3:
            c1          = 1/bt
            c2          = 1/bt
            H           = 0.5*(c1+c2)
            KG          = c1*c2
            gamma       = g0+2*kc*(H-H0)**2+kg*KG
            Gob_curv    = At*(2*kc*(H-H0)*(-2)/bt+kg*(-2)/bt**2)*1e-23
            mu0_correct = mu30

        m0           = (mu0_transfer+mu0_correct)*kT
        m1           = -1.0*y[-1]*kT
        m2           = -1.0*At*gamma*1e-23
        m3           = Eel[-1]
        m4           = kT*(Np+Nm)
        m5           = -kT*NA*self.Vaque*(cp0+cm0)*1e-30
        Mu_Ama       = m0+m1-1/n*(m2+m3+m4+m5)                      # chemical potential of Amphiphile in aggregate, [0]
        Mu_Amw       = kT*log(cm0/55500.)                           # [J]
        MuEx_Ion_m   = kT*log(cm0/1000.)-RT*self.cp0*1e-30*v1       # v1 [A]^3 monomer volume
        MuEx_Ion_p   = kT*log(cp0/1000.)-RT*self.cm0*1e-30*v2       # v2 [A]^3 ion volume defined in setup.py
        MuEx_Am      = (MuEx_Ion_m + MuEx_Ion_p)
        MuEx_H2O     = -RT*(self.cp0+self.cm0)*1e-30*VH2O

        PhiDb        = -1.0*PBsolution[-1,1]*kT/(e*self.b)          # [J]/([C]*[A]) # [A], Angstrom 
        BCV          = -1.0*sigma/(eps0*D)*1e-10                    # [J]/([C]*[A])
        Gob          = (2*m3+m2+Gob_curv)/b

        if d==3:
            Mu_Ama      = Mu_Ama   + kT*log(self.Cm/55500.)/n
            MuEx_H2O    = MuEx_H2O - RT*self.Cm*1e-30*VH2O
            Gob         = Gob      - 3*kT*log(self.Cm/55500.)/(b+L)    

#--------------------------#
#    objective function    #
#--------------------------#
        Diff0   = (Mu_Ama-Mu_Amw)/(Mu_Ama)
        Diff1   = (PhiDb-BCV)/BCV
        Diff2   = Gob/kT

        if out_opt == 2:
            output  = np.array([Diff0,Diff1])
            print(' '.join('%-8.4f '%it for it in output)+'\n')
            return output
        if out_opt == 3:
            output  = np.array([Diff0,Diff1,Diff2])
            print(' '.join('%-8.4f '%it for it in output)+'\n')
            return output

        print('')
        print('L    = %8.4f   b     = %-+8.4f'%(L,b))
        print('l0  = %-+8.4f, bt    = %-+8.4f'%(l0,bt))
        print('g0  = %8.4f,   H0    = %8.4f,  kc    = %8.4f'%(g0,H0,kc))
        print('mu0_transfer/kT      = %-+8.4e'%(mu0_transfer))
        print('m0/kT                = %-+8.4e'%(m0/kT))
        print('m1/kT                = %-+8.4e'%(m1/kT))
        print('m2/kT                = %-+8.4e, -gamma*At'%(m2/kT))
        print('m3/kT                = %-+8.4e, Eel'%(m3/kT))
        print('m4/kT                = %-+8.4e'%(m4/kT))
        print('m5/kT                = %-+8.4e'%(m5/kT))
        print('MuEx_Ion_m/kT        = %-+8.4e\nMuEx_Ion_p/kT        = %-+8.4e'%(MuEx_Ion_m/kT,MuEx_Ion_p/kT))
        print('')
        print('Mu_Ama/kT            = %-+8.4g\nMu_Amw/kT            = %-+8.4g'%(Mu_Ama/kT,Mu_Amw/kT))
        print('PhiDb                = %-+8.4g\nBCV                  = %-+8.4g'%(PhiDb,BCV))
        print('Gob/kT               = %-+8.4g'%(Gob/kT))
        print('Gob_curv/kT          = %-+8.4f'%(Gob_curv/kT))
        print('2*m3+m2              = %-+8.4g'%((2*m3+m2)/kT))
        print('')
        print('Np                   = %-+8.4g\nNm                   = %-+8.4g'%(Np,Nm))
        print('Np-Nm                = %-+8.4g\nn                    = %-+8.4g'%(Np-Nm,n))
        print('')
        print('Diff0                = %-+8.4g'%Diff0)
        print('Diff1                = %-+8.4g'%Diff1)
        print('Diff2                = %-+8.4g'%Diff2)

        R_MuEx_Am   = MuEx_Am /kT                     # Reduced Excess chem. potential of neutral Amphiphile
        R_MuEx_H2O  = MuEx_H2O/kT*1e3                 # Reduced Excess chem. potential of water
        CA          = Np/NA/self.Vcell*1e30
        CamW        = Nm/NA/self.Vaque*1e30
        CamA        = (Np-Nm)/NA/self.Vaque*1e30
        MassF       = (Np*self.M_amph)/(Np*self.M_amph+(self.Vaque-Np*v2-Nm*v1)/VH2O*MH2O)*100
        output      = np.array((L,b,cm0,cp0,R_MuEx_Am,R_MuEx_H2O,MassF,Diff0,Diff1,Diff2))
        return output


    def DE_BCMU(self, guess = None):
        """
        set up 2 non-linear equations to be estimated by DE
        """
        if guess is not None:
            self.cm0, self.cp0 = guess
        try:
            Calculated = self.Calculation(out_opt = 2)
        except (Warning,RuntimeError,RuntimeWarning,ZeroDivisionError):
            return 1000.00
        sq = 0.0
        for item in Calculated:
            sq = sq + abs(item)
        return(sq)

    def DE_BCMUGOB(self, guess = None):
        """
        set up 3 non-linear equations to be estimated by DE
        """
        if guess is not None:
            self.b, self.cm0, self.cp0 = guess
        try:
            Calculated = self.Calculation(out_opt = 3)
        except (Warning,RuntimeError,RuntimeWarning,ZeroDivisionError):
            return 1000.00
        sq = 0.0
        for item in Calculated:
            sq = sq + abs(item)
        return(sq)

    def Ui_Cal(self,x):
        Bi  = @Bi*1e-20              # J*[A]^3
        pis = np.pi**0.5
        b   = self.b
        L   = self.L
        if type(x) == float:
            if x == 1:
                r1  = 1e-6
                output1 = 16*Bi/(3*pis*ai**3)/kT
            else:
                r1  = b/x-b
                r12 = r1**2;    ai2 = ai**2
                r14 = r1**4;    ai4 = ai**4
                hi1 = 1+2*r1/(pis*ai)*(2*r12/ai2-1)*exp(-r12/ai2)-(1+4*r14/ai4)*erfc(r1/ai)
                output1 = Bi/r1**3*hi1/kT     # [J]
            r2  = 2*L-r1
            r22 = r2**2; ai2 = ai**2
            r24 = r2**4; ai4 = ai**4
            hi2 = 1+2*r2/(pis*ai)*(2*r22/ai2-1)*exp(-r22/ai2)-(1+4*r24/ai4)*erfc(r2/ai)
            output2 = Bi/r2**3*hi2/kT     # [J]
            output  = output1 - output2

        if type(x) == np.ndarray:
            r1 = b/x-b
            r1[-1] = 1e-6
            r12  = r1**2; ai2 = ai**2
            r14  = r1**4; ai4 = ai**4
            hi1  = 1+2*r1/(pis*ai)*(2*r12/ai2-1)*exp(-r12/ai2)-(1+4*r14/ai4)*erfc(r1/ai)
            output1 = Bi/r1**3*hi1/kT     # [J]
            output1[-1] = 16*Bi/(3*pis*ai**3)/kT
            r2   = 2*L-r1
            r22  = r2**2; ai2 = ai**2
            r24  = r2**4; ai4 = ai**4
            hi2  = 1+2*r2/(pis*ai)*(2*r22/ai2-1)*exp(-r22/ai2)-(1+4*r24/ai4)*erfc(r2/ai)
            output2 = Bi/r2**3*hi2/kT     # [J]
            output  = output1 - output2

        return output

    def PB_Equation(self,initC,x):
        """
        Poisson-Boltzmann equation
        y and y_prime, reduced potential and electrofiled, y = e*phi/kT
        initC,  the vector of initial Condition for y and y_prime
        x,      the variable, x = b/r
        return the y_prime and y_pp
        """
        cm0        = self.cm0*1.0e-30       # mol/[A]^3 
        cp0        = self.cp0*1.0e-30       # mol/[A]^3
        K1         = NA*(e*self.b)**2/(eps0*D*kT)*1.0e+10          # [A]^3; Josson1987, Equ.20
        Dim        = 3 - self.Dim
        y,y_prime  = initC                   
        Ui         = self.Ui_Cal(x)
        dydx       = [y_prime, -Dim/x*y_prime - K1/x**4*(cp0*exp(-y-Ui)-cm0*exp(y))]
        return dydx

    def PBE_Solver(self):
        """
        solving Poisson-Boltzmann equation by scipy odeint
        """
        y0 = [0.0,0.0]
        x  = np.linspace(self.b/self.B,1.0,nx)  # x is for the PBE_Solver
        with stdout_redirected():
            pb_sol = odeint(self.PB_Equation,y0,x)
        return pb_sol

def computing(mode,num):

    # mode, 0 bopt;     1 bfix;
    # num,  1 lamellar; 2 hexagonal; 3 micelle;

    out_data    = 'out'+str(mode)+str(num)
    filename    = 'summary'+str(mode)+str(num)

    Gotten      = []

    L0          = np.arange(1.5,15.6,0.1)
    lmax        = l_tail_exd*(v1/v0)**(1/num)

    for index in range(len(L0)):
        L   = L0[index]
        print('\nL = %8.4f Start off\n'%L)
        guess,bound = compute_guess(num,mode,Gotten,L,lmax)
        b,cm0,cp0   = guess
        b1,b2       = bound[0]
        cm01,cm02   = bound[1]
        cp01,cp02   = bound[2]

        with open(filename,'a') as dump:
            dump.write('L   = %8.4f\n'%(L))
            dump.write('Guess:  b = %8.4f cm0 = %8.4f cp0 = %8.4f\n'%(b,cm0,cp0))
            dump.write('Bounds    = [(%8.4f,%8.4f),(%8.4f,%8.4f),(%8.4f,%8.4f)]\n'%(b1,b2,cm01,cm02,cp01,cp02))

        A0 = Aggregate(num,L,b,cm0,cp0) 

        fail_times = 0
        for seed1 in range(1,6):
            try:
                with stdout_redirected():
                    if mode == 1:
                        Bounds2 = [(cm01,cm02),(cp01,cp02)]
                        result  = differential_evolution(A0.DE_BCMU,   Bounds2,seed=seed1)
                    if mode == 0:
                        Bounds3 = [(b1,b2),(cm01,cm02),(cp01,cp02)]
                        result  = differential_evolution(A0.DE_BCMUGOB,Bounds3,seed=seed1)
            except (Warning,RuntimeWarning):
                continue
            
            with open(filename,'a') as dump:
                dump.write('seed1 = %2d result = %12.8f\n'%(seed1,result.fun))

            if result.fun > Accuracy:
                with open(filename,'a') as dump:
                    dump.write('Fails, computation continue\n')
                print('Fails, computation continue')
                fail_times  = fail_times+1
                continue
            else:
                if mode == 1:
                    cm0,cp0 = result.x
                    with open(filename,'a') as dump:
                        dump.write('Succeed b = %8.4f cm0 = %8.4f cp0 = %8.4f\n\n'%(b,cm0,cp0))
                    print('compute DE2 Succeed  b = %8.4f cm0 = %8.4f cp0 = %8.4f'%(b,cm0,cp0))
                if mode == 0:
                    b,cm0,cp0 = result.x
                    with open(filename,'a') as dump:
                        dump.write('Succeed b = %8.4f cm0 = %8.4f cp0 = %8.4f\n\n'%(b,cm0,cp0))
                    print('compute DE3 Succeed  b = %8.4f cm0 = %8.4f cp0 = %8.4f'%(b,cm0,cp0))
                Gotten.append((L,b,cm0,cp0))    
                break

        if fail_times < 5:
            computed = 1
            A1       = Aggregate(num,L,b,cm0,cp0)
            output   = A1.Calculation()
            with open(out_data,'a') as log_dump:
                log_dump.write(' '.join('%-12.8f'%it for it in output)+'\n')
            print('L = %8.4f Done\n'%L)
        else:
            with open(filename,'a') as dump:
                dump.write('cm0 cp0 DE fails 5 times, do it manually\n')
            computed = 0
            print('cm0 cp0 DE fails 5 times, do it manually\n')
            break

    return computed




