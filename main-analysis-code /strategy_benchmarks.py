import numpy as np
import pandas as pd
import random
import support_funcs as sf


#----------------------------------------------------------------------------------
# strategy benchmarks for both interval (simulate_strategies) and probabilistic (simulate_random_prob) schedules
#----------------------------------------------------------------------------------

def simulate_strategies(check_times,bmeta,ses_n,stagger=True,strategy='random'):
    
    """
    Simulate port-checking behavior under different sampling strategies.

    This function replays a sequence of port-check times and assigns choices
    according to a specified strategy, then evaluates whether each choice would
    have been rewarded based on the reward schedule for a given session.

    Available strategies are:
        - 'random'  : choose ports randomly
        - 'io'      : ideal observer; choose a currently baited port when possible
        - 'winstay' : stay on the previous port after reward, otherwise switch

    Reward availability is updated across discrete timesteps using the session's
    programmed reward intervals. Rewards can be generated either with staggered
    timing offsets (to approximate the real task structure) or with simple
    fixed periodic schedules.

    Parameters
    ----------
    check_times : array-like
        Times of port checks in seconds, regardless of which port was checked.
        These times define when sampling occurs in the simulation.
    bmeta : dict
        Behavioral metadata dictionary. Must contain session reward interval
        information under:
            bmeta['intervals'][ses_n - 1]
    ses_n : int
        Session number to simulate. Used to select the corresponding set of
        reward intervals from `bmeta`.
    stagger : bool, default=True
        If True, simulate staggered reward availability using offset reward
        times for each port. If False, use simple periodic reward schedules
        based only on each port's reward interval.
    strategy : {'random', 'io', 'winstay'}, default='random'
        Strategy used to assign port choices at each check time.

    Returns
    -------
    strategy_df : pandas.DataFrame
        DataFrame with one row per simulated check, containing:
            - 'outcomes'    : int
                Whether the sampled port was rewarded (1) or not (0)
            - 'check_times' : float
                Check times in seconds
            - 'rewards_col' : list or ndarray
                Running total of rewards collected at each port after each check
            - 'port_list'   : int
                Zero-indexed port chosen on each simulated check

    Notes
    -----
    - Ports are represented internally as indices 0-5.
    - For the 'random' strategy, port identities are assigned to check times
      before the simulation begins.
    - For the 'io' strategy, if multiple ports are currently baited, the port
      with the longest reward interval is selected.
    - If no ports are baited under 'io', a random port is chosen.
    - For the 'winstay' strategy, the first choice is random.
    - If the number of simulated outcomes is shorter than `check_times`,
      the final check time is dropped to keep array lengths aligned.
    """
    
    
    rand_ports = [random.choice([1,2,3,4,5,6]) for nchecks in range(0,len(check_times))]

    longform_checks = sf.longform(check_times) 
    if len(longform_checks[longform_checks>0]) == len(rand_ports):
        longform_checks[longform_checks>0] = rand_ports
    else:
        rand_ports = rand_ports[:-1]
        longform_checks[longform_checks>0] = rand_ports
    max_steps = len(longform_checks)
    
    # get intervals - and (in the case of within session changes) when they change
    ris_all     = bmeta['intervals'][ses_n-1]
    avail_r     = np.zeros(6)  
    outcomes    = []
    RC          = np.zeros(6)  
    RC_t        = []
    port_list   = []
    
    # editing intervals to be true (staggered) 
    jitter         = [(2*i) for i, ris in enumerate(ris_all)]
    stag_rew_times = [ri*np.arange(1,max_steps/ri)+jitter[i] for i,ri in enumerate(ris_all)]

    for ts in range(0,max_steps): 
        
        for i,ri in enumerate(ris_all):
            if ts == 0:
                avail_r[i] = 0
            
            # using staggered reward times
            if stagger:
                if ts in stag_rew_times[i]: # rewarded timestep 
                    avail_r[i] = 1 
            # using simple reward times 
            else:
                if (ts+1)%ri ==0: 
                    avail_r[i] = 1 # rewarded timestep

        if longform_checks[ts]: # if this is a timestep that the mouse sampled a port 
            if strategy=='random':
                port = longform_checks[ts]-1
                
            elif strategy=='io':
                port = [i for i,avail in enumerate(avail_r) if avail == 1]
                if len(port) == 0: #if empty (no ports baited)
                    port =  int(random.choice([0,1,2,3,4,5]))
                elif len(port)>1: # if more than one port baited, pick the worst one 
                    port = port[np.argmax(np.array(ris_all)[port])] # index of highest number (interval) = worst port 
                else:
                    port = port
                    
            elif strategy=='winstay':
                if len(outcomes) == 0: # random first choice
                    port = int(random.choice([0,1,2,3,4,5]))
                else:
                    if outcomes[-1] == 1: # if the last outcome was a reward
                        port = int(port_list[-1]) # choose the same port again
                    else:
                        prev_port = port_list[-1]
                        while True: 
                            port = int(random.choice([0,1,2,3,4,5]))
                            if port != prev_port:
                                break # get out the loop 
                            else:
                                pass
                            
            outcome = int(avail_r[port])

            outcomes.append(outcome)
            RC[port] += outcome
            RC_t.append(RC.copy())
            port_list.append(port)
            
            if outcome == 1:
                avail_r[port] = 0
                
    if len(check_times)>len(outcomes):
        check_times = check_times[:-1]
        
    # create a handy dataframe
    strategy_df = pd.DataFrame()
    strategy_df['outcomes']    = outcomes
    strategy_df['check_times'] = check_times
    strategy_df['rewards_col'] = RC_t
    strategy_df['port_list']   = port_list

    return strategy_df


def simulate_random_prob(check_times,bmeta,ses_n,strategy='random'):
    """
    Simulate random port sampling behavior in a probabilistic reward task.

    At each check time, a port is selected uniformly at random (1–6), and
    reward delivery is determined probabilistically based on the port-specific
    reward probabilities for the given session.

    Parameters
    ----------
    check_times : array-like
        Times of port checks in seconds. Defines when sampling occurs.
    bmeta : dict
        Behavioral metadata dictionary. Must contain reward probabilities under:
            bmeta['intervals'][ses_n]
        where each entry is a list or array of length 6 (one per port).
    ses_n : int
        Session index used to select the corresponding reward probabilities.
        (Note: indexing here is assumed to be 0-based, unlike some other functions.)
    strategy : str, optional
        Included for API consistency with other simulation functions.
        Currently unused; only 'random' behavior is implemented.

    Returns
    -------
    strategy_df : pandas.DataFrame
        DataFrame with one row per check, containing:
            - 'outcomes'    : int
                Binary reward outcome (1 = rewarded, 0 = not rewarded)
            - 'check_times' : float
                Check times in seconds
            - 'rewards_col' : int
                Rewards collected at each check (identical to 'outcomes' here)
            - 'port_list'   : int
                Port chosen at each check (1–6 indexing)

    Notes
    -----
    - Port selection is uniform random across all 6 ports.
    - Reward outcomes are sampled independently using a Bernoulli distribution:
          outcome ~ Binomial(n=1, p=port_probability)
    - Unlike interval-based tasks, there is no reward accumulation or state;
      each check is independent.
    - `rewards_col` is redundant here (same as `outcomes`) but included for
      consistency with other simulation outputs.
    """
    rand_ports = [random.choice([1,2,3,4,5,6]) for nchecks in range(0,len(check_times))]
    probs      = bmeta['intervals'][ses_n]
    
    outcomes = []
    for choice in rand_ports:
        chosen_prob = probs[choice-1]
        outcomes.append(np.random.binomial(n=1, p=chosen_prob))
        
    # create a handy dataframe
    strategy_df = pd.DataFrame()
    strategy_df['outcomes']    = outcomes
    strategy_df['check_times'] = check_times
    strategy_df['rewards_col'] = outcomes 
    strategy_df['port_list']   = rand_ports
    
    return strategy_df
