def assure_tree(T, tmp_dir):
    if T is None:
        from treetime.utils import tree_inference
        T = os.path.basename(params.aln)+'.nwk'
        print("No tree given: inferring tree")
        tree_inference(params.aln, T, tmp_dir = tmp_dir)
        if os.path.isdir(tmp_dir):
            shutil.rmtree(tmp_dir)
        return 0
    elif not os.path.isfile(T):
        print("Input tree file does not exist:", params.tree)
        return 1


def create_gtr(params):
    from treetime import GTR
    model = params.gtr
    gtr_params = params["gtr_params"]
    if model == 'infer':
        gtr = GTR.standard('jc')
        infer_gtr = True
    else:
        infer_gtr = False
        try:
            kwargs = {}
            if gtr_params is not None:
                for param in gtr_params:
                    keyval = param.split('=')
                    if len(keyval)!=2: continue
                    if keyval[0] in ['pis', 'pi', 'Pi', 'Pis']:
                        keyval[1] = map(float, keyval[1].split(','))
                    elif keyval[0] not in ['alphabet']:
                        keyval[1] = float(keyval[1])
                    kwargs[keyval[0]] = keyval[1]
            else:
                print ("GTR params are not specified. Creating GTR model with default parameters")


            gtr = GTR.standard(model, **kwargs)
        except:
            print ("Could not create GTR model from input arguments. Using default (Jukes-Cantor 1969)")
            gtr = GTR.standard('jc')
    return gtr
