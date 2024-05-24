with open('selected.csv') as f:
    lines = f.read().strip().split('\n')

for l in lines[1:]:
    tokens = l.split(',')
    #parse 
    #param	m	k	w	tau	k0	k1	t0	t1	pk	sig	sig_em	ell
    #128s	1296	737	216	11	12	11	7	4	86	3653	4566	616
    param = tokens[0]
    m = int(tokens[1])
    k = int(tokens[2])
    w = int(tokens[3])
    tau = int(tokens[4])
    k0 = int(tokens[5])
    k1 = int(tokens[6])
    t0 = int(tokens[7])
    t1 = int(tokens[8])
    pk = int(tokens[9])
    sig = int(tokens[10])
    sig_em = int(tokens[11])
    ell = int(tokens[12])

    #output format is
    """
    param_128s = configuration_data()
    param_128s.set('PARAM', '128S')
    param_128s.set('PARAM_L', '128s')
    param_128s.set('LAMBDA', 128)
    param_128s.set('Nwd', 4)
    param_128s.set('L', 696)
    param_128s.set('TAU', 15)
    param_128s.set('K0', 9)
    param_128s.set('K1', 8)
    param_128s.set('T0', 8)
    param_128s.set('T1', 7)
    param_128s.set('SIG_SIZE', 4375)
    param_128s.set('PK_SIZE', 95) 
    param_128s.set('SK_SIZE', 32)
    param_128s.set('M', 1464)
    param_128s.set('N', 1464 - 833)
    param_128s.set('W', 244)
    param_128s.set('D', 6)
    parameter_sets += param_128s
    """
    
    print("param_%s = configuration_data()" % param)
    print("param_%s.set('PARAM', '%s')" % (param, param.upper()))
    print("param_%s.set('PARAM_L', '%s')" % (param, param))
    print("param_%s.set('LAMBDA', %s)" % (param, param[:3]))
    print("param_%s.set('Nwd', %d)" % (param, int(param[:3])/32))
    print("param_%s.set('L', %d)" % (param, ell))
    print("param_%s.set('TAU', %d)" % (param, tau))
    print("param_%s.set('K0', %d)" % (param, k0))
    print("param_%s.set('K1', %d)" % (param, k1))
    print("param_%s.set('T0', %d)" % (param, t0))
    print("param_%s.set('T1', %d)" % (param, t1))
    print("param_%s.set('SIG_SIZE', %d)" % (param, sig))
    print("param_%s.set('PK_SIZE', %d)" % (param, pk))
    print("param_%s.set('SK_SIZE', %d)" % (param, int(param[:3])/4))
    print("param_%s.set('M', %d)" % (param, m))
    print("param_%s.set('N', %d)" % (param, m - k))
    print("param_%s.set('W', %d)" % (param, w))
    print("param_%s.set('D', %d)" % (param, 6))
    print("parameter_sets += param_%s" % param)
    print('')




