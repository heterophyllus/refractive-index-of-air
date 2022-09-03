import numpy as np

class Air:
    def __init__(self) -> None:
        pass

    def refractive_index_ciddor(self, lambdamicron:float, t:float, p:float, RH:float, xCO2:float) -> float:
        
        w0 = 295.235
        w1 = 2.6422
        w2 = -0.03238
        w3 = 0.004028

        k0 = 238.0185
        k1 = 5792105
        k2 = 57.362
        k3 = 167917

        p_R1 = 101325
        T_R1 = 288.15

        Za = 0.9995922115

        rho_vs = 0.00985938

        R = 8.314472
        Mv = 0.018015

        S = 1.0/(lambdamicron**2)

        r_as = (10**(-8))*( k1/(k0-S) + k3/(k2-S) )
        r_vs = (1.022e-8)*(w0 + w1*S + w2*S**2 + w3*S**3)

        Ma = 0.0289635 + (1.2011e-8)*(xCO2-400)

        r_axs = r_as*( 1.0 + (5.34e-7)*(xCO2-450) )

        T = t + 273.15
        xv = self.humidity_for_ciddor(t,p,RH)
        Zm = self.compressibility(t,p,RH)
        rho_axs = p_R1*Ma/(Za*R*T_R1)
        rho_v   = xv*p*Mv/(Zm*R*T)
        rho_a   = (1-xv)*p*Ma/(Zm*R*T)

        n = 1.0 + (rho_a/rho_axs)*r_axs + (rho_v/rho_vs)*r_vs

        return n

    def compressibility(self, t:float, p:float, RH:float) -> float:
        T = t + 273.15

        a0 = 1.58123e-6
        a1 = -2.9331e-8
        a2 = 1.1043e-10

        b0 = 5.707e-6
        b1 = -2.051e-8

        c0 = 1.9898e-4
        c1 = -2.376e-6

        d = 1.83e-11
        e = -0.765e-8

        xv = self.humidity_for_ciddor(t,p,RH)
        tmp = a0 + a1*t + a2*t**2 + (b0+b1*t)*xv + (c0+c1*t)*xv**2
        Zm = 1.0 - (p/T)*tmp + (d + e*xv**2)*(p/T)**2

        return Zm

    def humidity_for_ciddor(self, t, p, RH):
        alpha = 1.00062
        beta  = 3.14e-8
        gamma = 5.60e-7

        f = alpha + beta*p + gamma*(t**2)

        # we are given relative humidity
        psv = self.saturation_vapor_pressure(t)
        xv = (RH/100)*f*psv/p

        return xv

    def saturation_vapor_pressure(self, t):
            K1 = 1.16705214528e+3
            K2 = -7.24213167032e+5
            K3 = -1.70738469401e+1
            K4 = 1.20208247025e+4
            K5 = -3.23255503223e+6
            K6 = 1.49151086135e+1
            K7 = -4.82326573616e+3
            K8 = 4.05113405421e+5
            K9 = -2.38555575678e-1
            K10 = 6.50175348448e+2
            
            T = t + 273.15
            omega = T + K9/(T-K10)
            A = omega**2 + K1*omega + K2
            B = K3*omega**2 + K4*omega + K5
            C = K6*omega**2 + K7*omega + K8
            X = -B + np.sqrt(B**2 - 4*A*C)

            psv = (10**6)*(2*C/X)**4
            
            return psv


