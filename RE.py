import BondGraphTools as bgt
import sympy as sp

def model():
    """ Acausal bond graph ER_abg.py
    Created by svgBondGraph at Wed Jul 24 08:53:50 2019 from ER_abg.svg

    Usage:
    import ER_abg; model = ER_abg.model()
    """

    model = bgt.new(name="ER")

    ## Junction 0:BGT1
    BGT1 = bgt.new('0')
    model.add(BGT1)

    ## Junction 0:BGT2
    BGT2 = bgt.new('0')
    model.add(BGT2)

    ## Junction 0:BGT3
    BGT3 = bgt.new('0')
    model.add(BGT3)

    ## Junction 0:BGT4
    BGT4 = bgt.new('0')
    model.add(BGT4)

    ## Junction 1:BGT0
    BGT0 = bgt.new('1')
    model.add(BGT0)

    ## Junction 1:BGT5
    BGT5 = bgt.new('1')
    model.add(BGT5)

    ## Component Ce:A
    K_A =  sp.symbols('K_A')
    RT = sp.symbols('RT')
    A = bgt.new('Ce',name='A',value={'k':K_A,'R':RT,'T':1},library='BioChem')
    model.add(A)

    ## Component Ce:AE
    K_AE =  sp.symbols('K_AE')
    RT = sp.symbols('RT')
    AE = bgt.new('Ce',name='AE',value={'k':K_AE,'R':RT,'T':1},library='BioChem')
    model.add(AE)

    ## Component Ce:B
    K_B =  sp.symbols('K_B')
    RT = sp.symbols('RT')
    B = bgt.new('Ce',name='B',value={'k':K_B,'R':RT,'T':1},library='BioChem')
    model.add(B)

    ## Component Ce:E
    K_E =  sp.symbols('K_E')
    RT = sp.symbols('RT')
    E = bgt.new('Ce',name='E',value={'k':K_E,'R':RT,'T':1},library='BioChem')
    model.add(E)

    ## Component Re:r1
    kappa_r1 =  sp.symbols('kappa_r1')
    RT = sp.symbols('RT')
    r1 = bgt.new('Re',name='r1',value={'r':kappa_r1,'R':RT,'T':1},library='BioChem')
    model.add(r1)

    ## Component Re:r2
    kappa_r2 =  sp.symbols('kappa_r2')
    RT = sp.symbols('RT')
    r2 = bgt.new('Re',name='r2',value={'r':kappa_r2,'R':RT,'T':1},library='BioChem')
    model.add(r2)

    ## Bonds
    bgt.connect(BGT5,BGT2)
    bgt.connect(BGT2,BGT0)
    bgt.connect(BGT1,BGT0)
    bgt.connect(BGT5,BGT4)
    bgt.connect(BGT1,A)
    bgt.connect(BGT2,E)
    bgt.connect(BGT3,AE)
    bgt.connect(BGT4,B)
    bgt.connect((r1,1),BGT3)
    bgt.connect(BGT3,(r2,0))
    bgt.connect((r2,1),BGT5)
    bgt.connect(BGT0,(r1,0))

    return model
