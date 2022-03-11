from setup import *

def L_K_Equ(K,L,sigma):
    y1 = -1.0*e*sigma/(D*eps0*2*kT)*1e-10/K
    y2 = np.tan(K*L)
    return abs(y1**2-y2**2)


def L_K_Cal(L,sigma):
    bounds = [(1e-3,np.pi*0.5/L)]
    for i in range(21):
        try:
            result = differential_evolution(L_K_Equ,bounds,args=(L,sigma),seed=i)
        except (Warning,RuntimeWarning):
            print('first DE Failure')
        if result.fun < 1e-10:
            K        = result.x[0]
            accuracy = result.fun
            cp0      = (K*1e+10)**2/(e**2/(2*D*eps0*kT))/NA
            break
    return cp0

    
def H_K_Equ(K,L,b,sigma):
    coeff   = 2*D*eps0*RT/(sigma*Fc*b)*1e+10
    y1      = K/(1+coeff*(1+K**2))
    y2      = np.tan(K*np.log((L+b)/b))
    return abs(y1**2-y2**2)


def H_K_Cal(L,b,sigma,bounds):
    got_it = 0
    for i in range(21):
        try:
            result = differential_evolution(H_K_Equ,bounds,args=(L,b,sigma),seed=i)
        except (Warning,RuntimeWarning):
            print('first DE Failure')
        if result.fun < 1e-6:
            print(bounds)
            K        = result.x[0]
            accuracy = result.fun
            coeff    = 2*D*eps0*RT/(sigma*Fc*b)*1e+10
            y1       = K/(1+coeff*(1+K**2))
            y2       = np.tan(K*np.log((L+b)/b))
            print('y1 = %8.4f y2 = %8.4f'%(y1,y2))
            print('accuracy = %8.4g '%accuracy)
            print('K        = %8.4g [1/A]'%K)
            cp0      = (1+K**2)*2*RT*D*eps0/(Fc*(b+L))**2*1e+20
            print('cp0      = %8.4g [mM]'%cp0)
            got_it = 1
            break
    if not got_it:
        cp0 = None
    return cp0


def compute_cp0(num,L,b):
    L       = float(L)
    b       = float(b)
    sa      = num*VAMP/b
    sigma   = -1.0*e/sa*1e+20
    got_it  = 0
    if num == 1:
        cp0_analytic = L_K_Cal(L,sigma)
        if cp0_analytic is not None:
            cp0 = cp0_analytic
            got_it = 1
    if num == 2:
        b_up    = pi*0.5/np.log((b+L)/b)
        for scalor in np.linspace(0.95,0.15,17):
            b_low   = b_up*scalor
            bounds  = [(b_low,b_up),]
            cp0_analytic = H_K_Cal(L,b,sigma,bounds)
            if cp0_analytic is not None:
                cp0 = cp0_analytic
                got_it = 1
                break
            else:
                continue
    if num == 3:
        cp0 = 0.7e+4
        got_it = 1

    if got_it:
        return(cp0)
    else:
        return (None)


def compute_guess(num,mode,Gotten,L,b):

    #mode 0 bopt, 1 bfix
    lmax = (v1/v0)**(1/num)*l_tail_exd

    No_Gotten = len(Gotten)

    if No_Gotten == 0:
        cm0     = 1e-6
        cp0     = compute_cp0(num,L,b)
 
    if No_Gotten > 0:
        GA      = np.array(Gotten[-4:])
        L_set   = GA[:,0]
        b_set   = GA[:,1]
        cm0_set = GA[:,2]
        cp0_set = GA[:,3]

        L       = L_set[-1]+0.5
        if mode == 0:
            b   = b_set[-1]
        if mode == 1:
            b   = lmax
        if num == 3:
            cm0 = cm0_set[-1]
            cp0 = cp0_set[-1]
        else:
            cm0 = cm0_set[-1]
            cp0 = compute_cp0(num,L,b)
            if cp0 is None:
                cp0 = cp0_set[-1]


    b1   = b*0.125
    b2   = b*8.000
    cm01 = 1.0e-9
    cm02 = cp0*0.250
    cp01 = cp0*0.125
    cp02 = cp0*8.000

    guess = (b,cm0,cp0)
    bound = [(b1,b2),(cm01,cm02),(cp01,cp02)]
    return(guess,bound)

