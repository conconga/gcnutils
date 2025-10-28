##@@##@@##@@ ##@@##@@##@@ ##@@##@@##@@ ###@@##@@##@ @#@@##@@##@@ ##@@##@@##@@ ##@@##@@##@@
#                                                                                        @
#             SoSode: System of Systems of Ordinary Differential Equations               @
#                                                                                        @
##@@##@@##@@ ##@@##@@##@@ ##@@##@@##@@ ###@@##@@##@ @#@@##@@##@@ ##@@##@@##@@ ##@@##@@##@@
#    \author Luciano A. Kruk         @ @                                                 @
#    \date   2020.03.13              @ @                                                 @
##@@##@@##@@ ##@@##@@##@@ ##@@##@@##@@ ###@@##@@##@ @#@@##@@##@@ ##@@##@@##@@ ##@@##@@##@@
#                                                                                        @
#     This object concatenate differential equations and parameter equations from        @
#     different sources into a unique global function to be called by the numerical      @
#     integration solver. The source of equations can also be other SoSode objects.      @
#                                                                                        @
#     The goal is to generate a function as:                                             @
#                                                                                        @
#         d_y/dt = f(y,t)                                                                @
#                                                                                        @
#     where f(y,t) is callable by the integrator, and y is the global state              @
#     [y1...yL].  Internally, hidden to the Integrator, f(y,t) can be a complicated      @
#     mess of ODE's and parameters.                                                      @
#                                                                                        @
#     The constituent functions to generate f(y,t) are f_i(y_i,t,P_i), where y_i is a    @
#     subset of y, and P_i is a subset of all available parameters. The function f_i     @
#     is registered with four descriptors:                                               @
#                                                                                        @
#       i_state: list with labels of input states of f_i                                 @
#       i_param: list with labels of input params of f_i                                 @
#       o_state: list with labels of y states by which derivatives are calculated by f_i @
#       o_param: list with labels of output params calculated by f_i                     @
#                                                                                        @
#     The order of parameters and returns is relevant when creating the lists. For       @
#     example, the function [a,b] = f([c,d],t,[e,f]) needs i_state=['c','d'] with        @
#     this order, as much as i_param=['e','f']. The return with also be                  @
#     o_state=['a','b'].                                                                 @
#                                                                                        @
#     A sanity check is provided to certify that all necessary inputs and params are     @
#     available to calculate the global state y.                                         @
#                                                                                        @
##@@##@@##@@ ##@@##@@##@@ ##@@##@@##@@ ###@@##@@##@ @#@@##@@##@@ ##@@##@@##@@ ##@@##@@##@@
##>>                                                                                  <<##
##@@##@@##@@ ##@@##@@##@@ ##@@##@@##@@ ###@@##@@##@ @#@@##@@##@@ ##@@##@@##@@ ##@@##@@##@@
class kSosodeFunction:
    def __init__(self, handler):
        self.handler = handler

        self.i_state = []
        self.i_param = []
        self.o_state = []
        self.o_param = []

    def set_i_state(self, i_state):
        self.i_state = i_state

    def set_i_param(self, i_param):
        self.i_param = i_param

    def set_o_state(self, o_state):
        self.o_state = o_state

    def set_o_param(self, o_param):
        self.o_param = o_param

    def __call__(self, *args):
        return self.handler(*args)

##@@##@@##@@ ##@@##@@##@@ ##@@##@@##@@ ###@@##@@##@ @#@@##@@##@@ ##@@##@@##@@ ##@@##@@##@@
##>>                                                                                  <<##
##@@##@@##@@ ##@@##@@##@@ ##@@##@@##@@ ###@@##@@##@ @#@@##@@##@@ ##@@##@@##@@ ##@@##@@##@@
class kSosode:

    def __init__(self, *fn_objs, reverse=False, order_states=None):
        """
        Use:
            reverse      :  True if integrator-input order will be (y,t) instead of (t,y)
            order_states :  list with the order of states in the state vector.
        """

        self.list_fn      = []
        self.reverse      = reverse
        self.order_states = order_states

        for i in fn_objs:
            self.register(i)


    def register(self, function):
        self.list_fn.append(function)


    def showregisteredfunctions(self):
        print()
        for i,j in enumerate(self.list_fn):
            print("-> function #{:d}:".format(i))
            print("     input states: {:s}".format(j.i_state.__str__()))
            print("     input params: {:s}".format(j.i_param.__str__()))
            print("    output states: {:s}".format(j.o_state.__str__()))
            print("    output params: {:s}".format(j.o_param.__str__()))


    def create_nets(self):

        ret = self._check_sanity()

        if ret:
            self._create_all_sets()
            self._create_net_params()
            self._create_net_states()
        else:
            raise(NameError("create_net: the net will not be created."))

        return ret


    def __call__(self, *args):
        """
        Calculates d(state)/dt : derivative of the states (same order as self.list_of_states)

        args:
            t        :  time
            y        :  state vector with same order as self.list_of_states
        """

        if self.reverse:
            t      = args[1]
            state  = args[0]
        else:
            t      = args[0]
            state  = args[1]

        if False:
            print('states(t={:f}) ='.format(t))
            print(state)

        params = self._calc_all_parameters(t)
        return self._calc_all_ddtstates(t, state, params)


    def _create_all_sets(self):
        self.list_of_outputs = list( set( sum( [i.o_state + i.o_param for i in self.list_fn ], [] ) ) )
        self.list_of_inputs  = list( set( sum( [i.i_state + i.i_param for i in self.list_fn ], [] ) ) )
        self.list_of_params  = list( set( sum( [i.i_param + i.o_param for i in self.list_fn ], [] ) ) )
        self.list_of_states  = list( set( sum( [i.i_state + i.o_state for i in self.list_fn ], [] ) ) )

        if self.order_states is not None:
            if set(self.list_of_states) == set(self.order_states):
                # set new order:
                self.list_of_states = self.order_states
            else:
                # not all states are to be updated?
                raise(NameError("are you sure not all states are to be updated?"))


    def _who_calcs_param(self, param):
        """
        Identifies which function calculates one parameter.

        return:
            index of the function calculating that particular parameter
            -1 when no functions calculates 'param'.
        """

        ret = -1
        for i,j in enumerate(self.list_fn):
            if param in j.o_param:
                # i,j calculate the parameter.
                ret = i
                break

        return ret


    def _sequence_calc_parameter(self, param, _cur_seq=[]):
        """
        Calculates a list of functions necessary to calculate the parameter 'param', in cascade.

        return:
            [...]   list of functions, from the first to the last in the recursion.
            -1      if there is a infinity loop to calculate 'param'.
        """

        if False:
            level  = len(_cur_seq)
            spaces = "".join([' ' for i in range(4*level)])
            print("{:s}starting _sequence_calc_parameter({:s}, {:s})".format(spaces, param.__str__(), _cur_seq.__str__()))

        ret = []
        i   = self._who_calcs_param(param)

        # no function calculates 'param'
        # or
        # infinity loop:
        if (i < 0) or (i in _cur_seq): return -1

        ret.append(i)

        # recursively look for the next functions:
        #print("{:s}  {:s} depends on {:s}...".format(spaces, param, self.list_fn[i].i_param.__str__()))
        set_fn_this_step = []
        for j in self.list_fn[i].i_param:
            k = self._sequence_calc_parameter(j, _cur_seq+ret)
            #print("{:s}{:s} = self._sequence_calc_parameter({:s}, {:s})".format(spaces, k.__str__(), j, (_cur_seq+ret).__str__()))

            if not hasattr(k, '__iter__'):
                if k < 0: return -1

            set_fn_this_step += list(k)

        # remove redundant items:
        set_fn_this_step = set(set_fn_this_step)

        if len(set_fn_this_step) > 0:
            ret += list(set_fn_this_step)

        #print("{:s}Return: {:s}".format(spaces, ret.__str__()))
        return ret


    def _check_sanity(self):
        self._create_all_sets()
        ret = True

        # 1) each output is calculated by a unique function.
        for i in self.list_of_outputs:
            c = 0
            for j in self.list_fn:
                c += 1 if i in j.o_param + j.o_state else 0

            if c != 1:
                print("sanity check: output '{:s}' is calculated by more than one function.".format(i))
                ret = False

        # 2) each function calculates at least one output.
        for i,j in enumerate(self.list_fn):
            if len(j.o_param + j.o_state) == 0:
                print("sanity check: function #{:d} do not calculate any output.".format(i))
                ret = False

        # 3) the necessary inputs for the functions are in the list_of_outputs.
        for i in self.list_of_inputs:
            if i not in self.list_of_outputs:
                print("sanity check: input '{:s}' is not calculated by any function.".format(i))
                ret = False

        # 4) derivative functions cannot calculate parameters AND parameter
        #    functions cannot calculate derivatives.
        for i,j in enumerate(self.list_fn):
            if (len(j.o_state) > 0) and (len(j.o_param) > 0):
                print("sanity check: funtion #{:d} calculates derivatives AND parameters.".format(i))
                ret = False

        # 5) at least one state to integrate.
        if len(self.list_of_states) == 0:
            print("sanity check: there is no derivative of states to integrate.")
            ret = False

        return ret


    def _create_net_params(self):
        """
        Creates a constant net to calculate the parameters. It describes the flow necessary to
        calculate all parameters inside the list self.net_params.

        It only need to be called once, as soon as the equations are defines.
        """

        net            = []
        nb_param       = len(self.list_of_params)
        list_of_params = list(self.list_of_params)

        # True when a parameter is calculated:
        flag_param_calculated = [0 for i in range(nb_param)]

        # True when a function is called:
        flag_called_functions = [0 for i in range(len(self.list_fn))]

        # Results for the calculated parameters:
        #values_params         = [0 for i in range(nb_param)]

        for idx_p, p in enumerate(list_of_params):

            # already calculated?
            if flag_param_calculated[idx_p]: continue

            # with the reverse(), the first item does not have dependences other than time:
            seq_fn = self._sequence_calc_parameter(p)
            if len(seq_fn) > 0: seq_fn.reverse()

            for i in seq_fn:

                # already called?
                if flag_called_functions[i]: continue

                # create a list/tuple with the arguments to call the function[i]:
                net.append({ "cmd":"new arg", "count": len(self.list_fn[i].i_param) })

                #args = []
                for ip in self.list_fn[i].i_param:
                    idx = list_of_params.index(ip)
                    #args.append(values_params[idx])
                    net.append({ "cmd":"arg", "fn":i, "input_param":ip, "idx_param":idx })

                net.append({ "cmd":"callfn", "fn":i })
                #fn_out = self.list_fn[i](t, *args)
                flag_called_functions[i] = True
                net.append({ "cmd":"save result", "count": len(self.list_fn[i].o_param) })

                for k,j in enumerate(self.list_fn[i].o_param):
                    idx                        = list_of_params.index(j)
                    if len(self.list_fn[i].o_param) > 1:
                        net.append({ "cmd":"value vector", "idx_param":idx, "idx_fnout":k })
                        #values_params[idx]     = fn_out[k]
                    else:
                        net.append({ "cmd":"value scalar", "idx_param":idx })
                        #values_params[idx]     = fn_out
                    flag_param_calculated[idx] = True

                if all(flag_param_calculated): break

        self.net_params = net


    def _calc_all_parameters(self, t):
        """
        Call self._create_net_params() before.
        """

        idx = 0
        net = self.net_params
        val = [0 for i in range(len(self.list_of_params))]

        while idx < len(net):
            if net[idx]['cmd'] == "new arg":
                args  = []
                count = net[idx]['count']
                for i in range(count):
                    args.append(val[ net[idx+i+1]['idx_param'] ])

                idx += net[idx]['count'] + 1

            elif net[idx]['cmd'] == "callfn":
                fn_out = self.list_fn[net[idx]['fn']](t, *args)
                idx += 1

            elif net[idx]['cmd'] == "save result":
                count = net[idx]['count']
                if count == 1:
                    val[ net[idx+1]['idx_param'] ] = fn_out
                else:
                    for i in range(count):
                        val[ net[idx+i+1]['idx_param'] ] = fn_out[i]

                idx += net[idx]['count'] + 1

            else:
                print("what do you mean here??")

        return val


    def _calc_all_parameters_full(self, t):
        """
        Calculates all parameters in the set 'self.list_of_params', in the same order.
        Do not use this one: call self._calc_all_parameters() instead.
        """

        nb_param       = len(self.list_of_params)
        list_of_params = self.list_of_params

        # True when a parameter is calculated:
        flag_param_calculated = [0 for i in range(nb_param)]

        # True when a function is called:
        flag_called_functions = [0 for i in range(len(self.list_fn))]

        # Results for the calculated parameters:
        values_params         = [0 for i in range(nb_param)]

        for idx_p, p in enumerate(list_of_params):

            # already calculated?
            if flag_param_calculated[idx_p]: continue

            # with the reverse(), the first item does not have dependences other than time:
            seq_fn = self._sequence_calc_parameter(p)
            if len(seq_fn) > 0: seq_fn.reverse()

            for i in seq_fn:

                # already called?
                if flag_called_functions[i]: continue

                # create a list/tuple with the arguments to call the function[i]:
                args = []
                for ip in self.list_fn[i].i_param:
                    idx = list_of_params.index(ip)
                    args.append(values_params[idx])

                fn_out = self.list_fn[i](t, *args)
                flag_called_functions[i] = True

                for k,j in enumerate(self.list_fn[i].o_param):
                    idx                        = list_of_params.index(j)
                    if hasattr(fn_out, '__iter__'):
                        values_params[idx]     = fn_out[k]
                    else:
                        values_params[idx]     = fn_out
                    flag_param_calculated[idx] = True

                if all(flag_param_calculated): break

        return values_params


    def _create_net_states(self):
        """
        Creates a constant net to calculate the states. It describes the flow necessary to
        calculate all states inside the list self.net_states.

        It only need to be called once, as soon as the equations are defines.
        """

        net            = []
        nb_states      = len(self.list_of_states)
        list_of_states = self.list_of_states

        for i in range(len(self.list_fn)): # all state functions are going to be called

            # only state functions:
            if len(self.list_fn[i].o_state) == 0: continue

            # prepare the arguments for function [i]:
            net.append({ "cmd":"new params", "count":len(self.list_fn[i].i_param) })
            for j in self.list_fn[i].i_param:
                idx = self.list_of_params.index(j)
                net.append({ "cmd":"arg", "idx_param":idx })
                #params.append(val_params[idx])

            net.append({ "cmd":"new state", "count":len(self.list_fn[i].i_state) })
            for j in self.list_fn[i].i_state:
                idx = self.list_of_states.index(j)
                net.append({ "cmd":"state", "idx_param":idx })
                #state.append(y[idx])

            # call the function:
            net.append({ "cmd":"callfn", "fn":i })
            #res = self.list_fn[i](t, state, *params)
            net.append({ "cmd":"save result", "count":len(self.list_fn[i].o_state) })

            # spread over the output list:
            for j,k in enumerate(self.list_fn[i].o_state):
                idx = list_of_states.index(k)
                net.append({ "cmd":"value state", "idx_param":idx })
                #val_ddtstates[idx] = res[j]

        self.net_states = net


    def _calc_all_ddtstates(self, t, y, val_params):
        """
        Call self._create_net_states() before.
        """

        idx = 0
        net = self.net_states
        ret = [0 for i in range(len(self.list_of_states))]

        while idx < len(net):
            if net[idx]['cmd'] == "new params":
                args  = []
                count = net[idx]['count']
                for i in range(count):
                    args.append(val_params[ net[idx+i+1]['idx_param'] ])

                idx += net[idx]['count'] + 1

            elif net[idx]['cmd'] == 'new state':
                state = []
                count = net[idx]['count']
                for i in range(count):
                    state.append( y[ net[idx+i+1]['idx_param'] ])

                idx += net[idx]['count'] + 1

            elif net[idx]['cmd'] == 'callfn':
                fn_out = self.list_fn[ net[idx]['fn'] ](t, state, *args)
                idx += 1

            elif net[idx]['cmd'] == 'save result':
                count = net[idx]['count']
                for i in range(count):
                    if hasattr(fn_out, '__iter__'):
                        ret[ net[idx+i+1]['idx_param'] ] = fn_out[i]
                    else:
                        ret[ net[idx+i+1]['idx_param'] ] = fn_out


                idx += net[idx]['count'] + 1

            else:
                print("you did something wrong, maaann!")

        return ret


    def _calc_all_ddtstates_full(self, t, y, val_params):
        """
        Calculates the derivatives of all states in self.list_of_states.
        """

        nb_states      = len(self.list_of_states)
        list_of_states = self.list_of_states
        val_ddtstates  = [0 for i in range(nb_states)]

        for i in range(len(self.list_fn)): # all state functions are going to be called

            # only state functions:
            if len(self.list_fn[i].o_state) == 0: continue

            # prepare the arguments for function [i]:
            params = []
            for j in self.list_fn[i].i_param:
                idx = self.list_of_params.index(j)
                params.append(val_params[idx])

            state = []
            for j in self.list_fn[i].i_state:
                idx = self.list_of_states.index(j)
                state.append(y[idx])

            # call the function:
            res = self.list_fn[i](t, state, *params)

            # spread over the output list:
            for j,k in enumerate(self.list_fn[i].o_state):
                idx = list_of_states.index(k)
                val_ddtstates[idx] = res[j]

        return val_ddtstates


##@@##@@##@@ ##@@##@@##@@ ##@@##@@##@@ ###@@##@@##@ @#@@##@@##@@ ##@@##@@##@@ ##@@##@@##@@
##>>                                                                                  <<##
##@@##@@##@@ ##@@##@@##@@ ##@@##@@##@@ ###@@##@@##@ @#@@##@@##@@ ##@@##@@##@@ ##@@##@@##@@

class kSosodeTests:
    def do_tests(self):
        import numpy            as np
        import scipy.integrate  as Int
        import matplotlib.pylab as plt

        def f1(y, t, PARAM):
            tta1 = PARAM[0]
            tta2 = PARAM[1]

            yp = np.asarray([
                (   4*y[0]) - (tta1*y[0]*y[1]),
                (tta2*y[2]) - (  2.*y[2]*y[3]),
            ])

            return yp

        ##@@##@@##@@ ##@@##@@##@@

        def f2(y, t):

            yp = np.asarray([
                (  -3*y[1]) + (  3.*y[0]*y[1]),
                (  -3*y[3]) + (  3.*y[2]*y[3]),
            ])

            return yp

        ##@@##@@##@@ ##@@##@@##@@

        def g1(y, t):
            return (2,4) if t < 5 else (1,5)

        ##@@##@@##@@ ##@@##@@##@@

        fn1 = kSosodeFunction(f1)
        fn1.set_i_state([ "y 1", "y 2", "y 3", "y 4" ])
        fn1.set_i_param([ "theta 1", "theta 2" ])
        fn1.set_o_state([ "y 1", "y 3" ])

        fn2 = kSosodeFunction(f2)
        fn2.set_i_state([ "y 1", "y 2", "y 3", "y 4" ])
        fn2.set_o_state([ "y 2", "y 4" ])

        fn3 = kSosodeFunction(g1)
        fn3.set_o_param([ "theta 1", "theta 2" ])

        a = fn1( [1,2,3,4], 0, [1,2] )
        
        ##@@##@@##@@ ##@@##@@##@@
        # the main object:
        print()
        print("--  REGULAR TESTS...  --")
        b = kSosode(fn1, fn2)
        b.register(fn3)
        b._check_sanity()


        ##@@##@@##@@ ##@@##@@##@@
        print()
        print("--  SANITY CHECKS...  --")

        # erro: 1
        fn1.set_o_state([ "y 1", "y 3", "y 2" ])
        b._check_sanity()

        # erro: 2 and 3
        fn1.set_o_state([])
        fn1.set_o_param([])
        b._check_sanity()

        # erro: 4
        fn1.set_o_state([ "y 1", "y 3" ])
        fn1.set_o_param([ "y 1", "y 3" ])
        b._check_sanity()

        ##@@##@@##@@ ##@@##@@##@@
        print()
        print("--  SEQUENCE CALC PARAMS  --")

        fn1 = kSosodeFunction(None)
        fn1.set_o_param(['t1','t2'])
        fn1.set_i_param(['t3','t4'])

        fn2 = kSosodeFunction(None)
        fn2.set_o_param(['t3','t4'])
        fn2.set_i_param(['t5','t6'])

        fn3 = kSosodeFunction(None)
        fn3.set_o_param(['t5'])
        #fn3.set_i_param(['t1'])

        fn4 = kSosodeFunction(None)
        fn4.set_o_param(['t6'])

        b = kSosode( fn1, fn2, fn3, fn4 )
        if b._who_calcs_param('t1') != 0: print("ERRO 1")
        if b._who_calcs_param('t4') != 1: print("ERRO 2")

        for i in range(6):
            j = "t{:d}".format(i+1)
            k = b._sequence_calc_parameter(j)

            if k != [ [0, 1, 2, 3], [0, 1, 2, 3], [1, 2, 3], [1, 2, 3], [2], [3] ][i]:
                print("ERROR while calculating {:s}".format(j))
            else:
                print("for '{:s}', these functions are called: {:s}".format(j, k.__str__()))

        # this line will create and infinity loop:
        print("adding an infinity loop...")
        fn3.set_i_param(['t1'])

        for i in range(6):
            j = "t{:d}".format(i+1)
            k = b._sequence_calc_parameter(j)

            if k != [ -1, -1, -1, -1, -1, [3] ][i]:
                print("ERROR while calculating {:s}".format(j))
            else:
                print("for '{:s}', these functions are called: {:s}".format(j, k.__str__()))

        # test calculation of parameters:
        print()
        print("--  TEST CALC of PARAMETERS  --")

        def g0(t, *args): return (args[0]+args[1])
        def g1(t, *args): return (3 if t < 2 else 4)
        def g2(t, *args): return (args[0] + 1.)*args[1]
        def g3(t, *args): return [ (-2 if t < 2.5 else 0), 2. ]

        gn0 = kSosodeFunction(g0)
        gn0.set_o_param(['t1'])
        gn0.set_i_param(['t2', 't3'])

        gn1 = kSosodeFunction(g1)
        gn1.set_o_param(['t2'])

        gn2 = kSosodeFunction(g2)
        gn2.set_o_param(['t3'])
        gn2.set_i_param(['t4', 't5'])

        gn3 = kSosodeFunction(g3)
        gn3.set_o_param(['t4', 't5'])

        b = kSosode( gn0, gn1, gn2, gn3 )
        b._create_all_sets()

        try:
            ok = False
            b.create_nets()
        except:
            ok = True
        assert ok

        b._create_net_params()
        for t in np.arange(0,5,0.1):
            i = b._calc_all_parameters_full(t)
            j = b._calc_all_parameters(t) # needs _create_net_params()

            if i != j:
                print("ERROR calculating the arguments since {:s} != {:s}".format(i.__str__(), j.__str__()))

        # test calculation of parameters:
        print()
        print("--  TEST CALC of ddt_STATES  --")

        def f0(t, y, *args): 
            tta1 = args[0]
            tta2 = args[1]
            yp = np.asarray([
                (   4*y[0]) - (tta1*y[0]*y[1]),
                (tta2*y[2]) - (  2.*y[2]*y[3]),
            ])
            return yp

        def f1(t, y):
            yp = np.asarray([
                (  -3*y[1]) + (  3.*y[0]*y[1]),
                (  -3*y[3]) + (  3.*y[2]*y[3]),
            ])
            return yp

        fn0 = kSosodeFunction(f0)
        fn0.set_i_state([ 'y0', 'y1', 'y2', 'y3' ])
        fn0.set_i_param([ 't1', 't2' ])
        fn0.set_o_state([ 'y0', 'y2' ])

        fn1 = kSosodeFunction(f1)
        fn1.set_i_state([ 'y0', 'y1', 'y2', 'y3' ])
        fn1.set_o_state([ 'y1', 'y3' ])

        b.register(fn0)
        b.register(fn1)
        b.showregisteredfunctions()
        b.create_nets()

        for t in np.arange(0,5,0.1):
            i = b._calc_all_parameters(t)
            m = np.random.randn(4)
            j = b(t, m)
            k = b._calc_all_ddtstates_full(t, m, i)

            if j != k:
                print("ERROR calculating the ddt_states since {:s} != {:s}".format(j.__str__(), k.__str__()))

        # test calculation of parameters:
        print()
        print("--  TEST FULL SYSTEM  --")

        # from:
        #   http://support.ptc.com/help/mathcad/en/index.html#page/PTC_Mathcad_Help/example_solving_a_first_order_system_of_ODEs.html

        def tta(t):
            return 8 if t < 10 else -2

        def eq0(t, y, *args):
            tta = args[0]
            return (-8.*y[0]) + (tta*y[1])

        def eq1(t, y):
            return (30*y[0]) + y[1] - (y[0]*y[2])

        def eq2(t, y):
            return (y[0]*y[1]) - (8.*y[2]/3)

        fn0 = kSosodeFunction(eq0)
        fn0.set_i_state([ 'y0', 'y1' ])
        fn0.set_i_param([ 'tta' ])
        fn0.set_o_state([ 'y0' ])

        fn1 = kSosodeFunction(eq1)
        fn1.set_i_state([ 'y0', 'y1', 'y2' ])
        fn1.set_o_state([ 'y1' ])

        fn2 = kSosodeFunction(eq2)
        fn2.set_i_state([ 'y0', 'y1', 'y2' ])
        fn2.set_o_state([ 'y2' ])

        g0 = kSosodeFunction(tta)
        g0.set_o_param([ 'tta' ])

        b = kSosode( fn0, fn1, fn2, g0, reverse=True, order_states= [ 'y0', 'y1', 'y2' ] )
        b.create_nets()
        b.showregisteredfunctions()

        T = np.linspace(0,12.0,1000)
        R = Int.odeint(b, [-1, 0, 1], T, args=() )

        print()
        print("DONE.")

        # figure:
        fig = plt.figure(0).clf()
        ax  = plt.subplot(1,1,1)
        ax.plot(T, R)
        ax.grid(True)
        plt.show(block=False)

##@@##@@##@@ ##@@##@@##@@ ##@@##@@##@@ ###@@##@@##@ @#@@##@@##@@ ##@@##@@##@@ ##@@##@@##@@
