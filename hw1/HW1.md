
EEE 545 home work 1
===================
Pozar problems 4.14 and 4.18.

Due 1/21/15


    # Stuff usefull for all problems
    import numpy as np


    def dB(value):
        ''' converts value to dB scale
        
        Parameter
        ---------
        value:
            value to be conerted to dB
        
        return
        ------
        dB value of value
        '''
        return 20*np.log10(np.abs(value))
    
    def phase(value,deg=True):
        ''' Converts value to phase
        
        Parameter
        ---------
        value:
            Value to be converted to phase
        
        deg:
            degree flag, defualt is True
        
        return
        ------
        Phase of the value
        '''
        return np.angle(value,deg=deg)


    def latc(mag, angle):
        ''' latc or Linear magnitude and Angle To Complex
        is a function that converts a tuple of Magitude 
        and angle in degrees to a complex number
        
        Parameters
        ----------
        mag: 
            magitude
        angle:
            angle in degrees
        
        return
        ------
        complex number form of the magitude and degree form number
        '''
        return mag * np.exp( 1j * np.pi * angle / 180.)


    def S_ABCD(S, Zo=50.):
        ''' Converts a 2x2 S matrix into a 2x2 ABCD matrix
        
        Parameters
        ----------
        S: numpy.ndarray
            S matrix to be converted
        
        Zo:
            Impedance
        
        return
        ------
        ABCD matrix
        '''
        S11, S12, S21, S22 = S[0][0], S[0][1], S[1][0], S[1][1]
        A, B, C, D = (( 1 + S11) * ( 1 - S22) + S12 * S21)/( 2 * S21),\
            Zo * (( 1 + S11) * ( 1 + S22) - S12 * S21) / ( 2 * S21),\
            1/Zo * ((1 - S11) * (1 - S22) - S12 * S21) / ( 2 * S21),\
            ((1 - S11) * ( 1 + S22) + S12 * S21)/( 2 * S21)
        return (np.array([A, B, C, D],dtype=S.dtype)).reshape(S.shape)
    
    def ABCD_S(ABCD, Zo=50.):
        ''' Converts ABCD matrix to S amtrix
        
        Parameter
        ---------
        ABCD: numpy.ndarray
            ABCD matrix to be converted
            
        Zo:
            Impedance
            
        return
        ------
        S matrix
        '''
        A, B, C, D = ABCD[0][0], ABCD[0][1], ABCD[1][0], ABCD[1][1]
        S11, S12, S21, S22 = \
            (A + (B / Zo) - C * Zo - D) / (A + (B / Zo) + C * Zo + D),\
            2 * ( A * D - B * C ) / (A + (B / Zo) + C * Zo + D),\
            2 / (A + (B / Zo) + C * Zo + D),\
            (-1*A + (B / Zo) - C * Zo + D) / (A + (B / Zo) + C * Zo + D)
        return ( np.array([S11,S12, S21, S22], dtype=ABCD.dtype)).reshape(ABCD.shape)


    def Casecade( networks, Zo=50.):
        ''' Casecades python list of S matrixes
        
        Parameter
        ---------
        networks: Python list
            List of the S matrixes in order
        
        Zo:
            Impedance
        
        return
        ------
        Combined network
        '''
        abcd = S_ABCD(networks.pop(0), Zo=Zo)
        for n in networks:
            abcd = np.dot(abcd, S_ABCD(n))
        return ABCD_S( abcd )

Problem 1, Pozar 4.14
---------------------
A four-port network has the scattering matrix shown as follows.

a: Is this network lossless?

b: Is this network reciprocal?

c: What is the return loss at port 1 when all other ports are terminated with
matched loads?

d: What is the insertion loss and phase delay between ports 2 and 4 when all
other ports are terminated with matched loads?

e: What is the reflection coefficient seen at port 1 if a short circuit is
placed at the terminal plane of port 3 and all other ports are terminated with
matched loads?

$$ \begin{equation*} \mathbf{\mathbf{[S]}} = \begin{vmatrix} 0.178 \angle 90 &
0.6 \angle 45 & 0.4 \angle 45 & 0 \\
0.6 \angle 45 & 0 & 0 & 0.3 \angle -45\\
0.4 \angle 45 & 0 & 0 & 0.5 \angle -45\\
0 & 0.3 \angle -45 & 0.5 \angle -45 & 0 \end{vmatrix} \end{equation*}$$


    # Create the S matrix from the provide information
    S_abs = np.array([[ 0.178, 0.6  , 0.4  , 0    ],
                      [ 0.6  , 0    , 0    , 0.3  ],
                      [ 0.4  , 0    , 0    , 0.5  ],
                      [ 0    , 0.3  , 0.5  , 0    ]],dtype=float)
    S_angle = np.array([[ 90, 45, 45,  0],
                        [ 45,  0,  0,-45],
                        [ 45,  0,  0,-45],
                        [  0,-45,-45,  0]],dtype=float)
    S = latc(S_abs, S_angle)

###Problem 1, 4.14 a)
Is this network lossless?


    def is_lossless(network, max_loss=1.e-6):
        '''is_lossless 
        is a function to test if a network is within the bounds of being lossless
        
        Parameters
        ----------
        network: ndarray
            the scatering network for a device
        max_loss: float
            the maximum loss for a network to be lossless
        
        return
        ------
        True or False, conjurgate transpose of the network times the network'''
        SH_S = np.dot(network.conj().T,network)
        if np.all(np.abs(SH_S - np.eye(SH_S.shape[0])) <= max_loss):
            return True, SH_S
        else:
            return False, SH_S
        


    lossless, SH_S = is_lossless(S)
    if lossless is True:
        print('The network is lossless')
    else:
        print('The network is not lossless, unless %f is within the bounds of loss' % np.max(np.abs(SH_S)))

    The network is not lossless, unless 0.551684 is within the bounds of loss


###Problem 1, 4.14 b)
Is this network reciprocal?


    def is_reciprocal(network, tolorance=1.e-6):
        ''' is_reciprocal
        functioin to test of the network is reciprocal
        
        Parameters
        ----------
        network: ndarray
            the scattering parameter matrix for the network
        tolorance: float
            The tolorance value for if a network being reciprocal
        
        return
        ------
        True or False, the difference from being reciprocalible
        '''
        dif_rcp_net = network.T-network
        if np.all( np.abs( dif_rcp_net) <= tolorance):
            return True, dif_rcp_net
        else:
            return False, dif_rcp_net


    reciprocal, difference = is_reciprocal(S)
    if reciprocal is True:
        print('The network is reciprocal')
    else:
        print('The network is not reciprocal, unless %f is within the bounds of less' % np.max(np.abs(difference)))

    The network is reciprocal


###Problem 1, 4.14 c)
What is the return loss at port 1 when all other ports are terminated with
matched loads?


    V_in = np.array([1,0,0,0],dtype=complex).reshape((4,1))


    V_out = np.dot( S, V_in)


    print('Return loss for port 1 is %.5f linear mag or %.3f dB, and %.2f phase' % 
          (np.abs(V_out[0]),dB(V_out[0]),np.angle(V_out[0],deg=True)))

    Return loss for port 1 is 0.17800 linear mag or -14.992 dB, and 90.00 phase


###Problem 1, 4.14 d)
What is the insertion loss and phase delay between ports 2 and 4 when all other
ports are terminated with matched loads?


    V_in = np.array([0,1,0,0],dtype=complex).reshape(4,1)


    V_out = np.dot( S, V_in)


    print('Return loss from port 2 to  is %.3f linear mag or %.3f dB, and %.2f phase'%
          (np.abs(V_out[1]), 20*np.log10(np.abs(V_out[1])), np.angle(V_out[1],deg=True)))
    print('Insertion Loss and Phase delay from 2 to 4 is %.3f linear mag or %.3f dB, and %.2f phase'%
          (np.abs(V_out[3]),20*np.log10(np.abs(V_out[3])), np.angle(V_out[3],deg=True)))

    Return loss from port 2 to  is 0.000 linear mag or -inf dB, and 0.00 phase
    Insertion Loss and Phase delay from 2 to 4 is 0.300 linear mag or -10.458 dB, and -45.00 phase


    -c:2: RuntimeWarning: divide by zero encountered in log10


###Problem 1, 4.14 e)
What is the reflection coefficient seen at port 1 if a short circuit is placed
at the terminal plane of port 3 and all other ports are terminated with matched
loads?


    V_in = np.array([1, 0, 0, 0],dtype=complex).reshape(4,1)
    gamma = np.array([ 0, 0, -1, 0], dtype=complex).reshape(4,1)


    V_out = np.dot( S, V_in)
    V_out = np.dot( S, gamma*V_out+V_in)


    print('Return loss for port 1 is %f linear mag, and %f phase' % (np.abs(V_out[0]),np.angle(V_out[0],deg=True)))

    Return loss for port 1 is 0.018000 linear mag, and 90.000000 phase


Problem 2, Pozar 4.18
---------------------

A four-port network has the scattering matrix shown as follows.
If ports 3 and 4 are connected with a lossless matched transmission line with an
electrical length of 45,
find the resulting insertion loss and phase delay between ports 1 and 2.

$$ \begin{equation*} \mathbf{\mathbf{[S]}} = \begin{vmatrix}
0.2 \angle 50  & 0              & 0              & 0.4 \angle -45\\
0              & 0.6 \angle 45  & 0.7 \angle -45 & 0 \\
0              & 0.7 \angle -45 & 0.6 \angle 45  & 0 \\
0.4 \angle -45 & 0              & 0              & 0.3 \angle 45 \end{vmatrix}
\end{equation*}$$


    # Create the S matrix from the provide information
    S_abs = np.array([[ 0.2  , 0    , 0    , 0.4  ],
                      [ 0    , 0.6  , 0.7  , 0    ],
                      [ 0    , 0.7  , 0.6  , 0    ],
                      [ 0.4  , 0    , 0    , 0.3  ]],dtype=float)
    S_angle = np.array([[ 50,  0,  0,-45],
                        [  0, 45,-45,  0],
                        [  0,-45, 45,  0],
                        [-45,  0,  0, 45]],dtype=float)
    S = latc(S_abs, S_angle)


    S_line_abs = np.array([[ 0, 1],
                           [ 1, 0]], dtype=float)
    S_line_angle = np.array([[ 0, -45],
                             [-45,  0]], dtype=float)
    S_line = latc(S_line_abs, S_line_angle)

###Simple method
This is a simple method of solving this problem. Using

  1. Starting point

        <img src="files/images/4_18_1.png"/>

  2. First step, as the network can be separated into 2 separate networks.

        <img src="files/images/4_18_2a.png"/>
  3. Signal Flow Diagram of the network.

        <img src="files/images/4_18_3a.png"/>
  4. The simplified signal flow graph.

        <img src="files/images/4_18_4.png"/>
  5. From here the values for insertion loss and phase delay can be found.
      $$ T = 1 \angle -45 = {e}^{-1j* {\pi} \over {4}} $$

      $$ IL = {{{S}_{41} * {S}_{23} } \over {1 - {S}_{33} * {S}_{44} * T } }$$

      $$ IL = {{{0.4 \angle -45} * {0.7 \angle -45} } \over {1 - {0.6 \angle 45}
* {0.3 \angle 45} * 1 \angle -45 } }$$

      $$ IL = {{0.28 \angle -135} \over {1 - 0.18 \angle 0}} $$

      $$ IL = {{0.28 \angle -135} \over {0.82 \angle 0}} $$

      $$ IL = {0.34145 \angle -135.  }$$

      $$ IL = 9.333 dB \angle -135.$$


    IL = ((S[0][3]*S[1][2]*np.exp(-11j*np.pi/4))\
        /(1-S[2][2]*S[3][3]*(np.exp(-1j*np.pi/4))**2))
    print("Insertion Loss is %.3f dB and Phase delay is %.2f degrees" % (-dB(IL), phase(IL)))

    Insertion Loss is 9.333 dB and Phase delay is 135.00 degrees


### More advanced method
This is another more advanced method using convertion to ABCD matrixes and back
to S amtrixes.


    Sp14, Sp23 = (S[[0,3]][:,[0,3]]), (S[[1,2]][:,[1,2]])
    networks = [ Sp14, S_line, Sp23]
    eq_net = Casecade(networks)


    IL = eq_net[0,1]
    print("Insertion loss from Port 1 to 2 %.3f dB and %.2f degrees phase" % (-dB(IL), phase(IL)))

    Insertion loss from Port 1 to 2 9.333 dB and -135.00 degrees phase

