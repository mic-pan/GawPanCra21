""" Parameter identification using BondGraphTools """
## Quadratic programming stuff.
import quadprog

## Copy structures
import copy

## NumPy
import numpy as np

## Function from https://scaron.info/blog/quadratic-programming-in-python.html
def quadprog_solve_qp(P, q, G=None, h=None, A=None, b=None):
    qp_G = .5 * (P + P.T)   # make sure P is symmetric
    qp_a = -q
    if A is not None:
        qp_C = -np.vstack([A, G]).T
        qp_b = -np.hstack([b, h])
        meq = A.shape[0]
    else:  # no equality constraint
        qp_C = -G.T
        qp_b = -h
        meq = 0
    return quadprog.solve_qp(qp_G, qp_a, qp_C, qp_b, meq)[0]

## Function to compute phi from Phi subject to Phi>positive number
## NN Reduced N corresponding to known Phi
def quadsolve_phi(N0,N1,Phi0,A_ec=None,b_ec=None,Phi_min=1e-3,mu=1e-10):
    """Solve for species potentials  phi from known reaction potentials phi0"""
    ## Deduce number of species n_X and unknown reactions n_V
    (n_X,n_V) = N1.shape
    P = 1.0*N0@(N0.T) + mu*np.eye(n_X)
    q = (N0@Phi0).T
    G = 1.0*N1.T
    h = -Phi_min*np.ones((n_V))
    phi = quadprog_solve_qp(P, q, G=G, h=h, A=A_ec, b=b_ec)
    
    return phi

def Phi2phi(s,Phi_known,Ignore=[],Phi_min=1e-3,A_ec=None,b_ec=None,mu=1e-10,quiet = False):

    ## Decompose stoichiometric matrix N
    ## N0 for known and N1 for unknown Phi
    N = copy.copy(s['N'])
    N_0 = None
    N_1 = None
    Phi_0 = []
    for j,reac in enumerate(s['reaction']):
        if (reac in Phi_known.keys()) and not np.isnan(Phi_known[reac]):
            Phi_0.append(Phi_known[reac])
            if N_0 is None:
                N_0 = N[:,j]
            else:
                N_0 = np.vstack((N_0,N[:,j]))
        else:
            if reac in Ignore:
                if not quiet:
                        print(f'Ignoring {reac}')
            else:
                if N_1 is None:
                    N_1 = N[:,j]
                else:
                    N_1 = np.vstack((N_1,N[:,j]))

    ## Convert to numpy array          
    Phi_0 = np.array(Phi_0)

    ## Transpose
    N_0 = N_0.T
    N_1 = N_1.T
    
    n_X,n_V = N_0.shape
    print(f'Extracting {n_X} values of phi from {n_V} values of Phi')

    ## QP solution
    phi_est = quadsolve_phi(N_0,N_1,Phi_0,A_ec=A_ec,b_ec=b_ec,Phi_min=Phi_min,mu=mu)
    Phi_est = -N.T@phi_est
    Phi_0_est = -N_0.T@phi_est

    ## Error
    if not quiet:
        err = (np.linalg.norm(Phi_0_est-Phi_0))/len(Phi_0)
        print(f'Phi error = {err:.2}')

    ## Create dict versions
    phi_est_d = {}
    Phi_est_d = {}

    for i,spec in enumerate(s['species']):
        phi_est_d[spec] = phi_est[i]

    for i,reac in enumerate(s['reaction']):
        Phi_est_d[reac] = Phi_est[i]
    
    
    return phi_est,Phi_est,phi_est_d,Phi_est_d

    
