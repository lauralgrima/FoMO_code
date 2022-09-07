import numpy as np
import lick_df as ldf
# reward learning model for concurrent 6 port task

# things to add
# 1. right now every timepoint involves a check (which is obviously not reasonable). Instead use intervals from the mouse behaviour 
# 2. hazard function 
# 3. two second jitter between updating rewards becoming available at each port 
# 4. how to initialise probability at each port? 
# 5. how to calculate hazard rate
# 6. add ideal observer 
# 7. add space as a dimension. if the probability of reward at multiple ports is the same, then choose based on which port the agent is currently closest to
# softmax decision rule 

# akrami



# variables
n_ports   = 6
rrs       = [30, 60, 240, 240, 1200, 2400] # reward schedules in seconds 
policy    = 'win_stay'
decay     = 0.8                            # decay factor for win stay, between 0 and 1
max_steps = 11000                          # max no. timesteps in seconds 






# state machine
def state_machine(raw_data_folder,policy,rrs,decay,max_steps,n_ports,check_threshold,subject_ID,region):
    
    # get matrix of real data - port visits and times 
    
    # initial parameters
    pr      = np.zeros(n_ports)                # agent's estimated probability of reward at each port 
    avail_r = np.zeros(n_ports)                # available rewards 
    #g       = []                               # current goal state
    #ss      = list(range(0,n_ports+1))         # state space. 0 is when agent is not at any port 
    #cs      = []                               # current state 
    #s_hist  = []                               # state history
    RC      = np.zeros(n_ports)                # rewards collected at each port
        
    
    RC_t         = []
    pr_t         = []
    port_check_t = []
    avail_r_t    = []
    
    for t in list(range(0,max_steps)): 
        
        # assume that initially estimated prob of reward at each port is the same (maybe in future assume 0?)
        pr = np.array([1/6 for prob in avail_r]) 
    
        # update reward availability at this timepoint
        avail_r = update_avail_rr(avail_r,rrs,t)
        avail_r_t.append(avail_r.copy())
        
        # at every timepoint, choose whether to sample or not 
        # independent from where to check 
        
        # choose which port to check (here, whichever has maximum probability of giving reward. If more than one port, then choose randomly)
        port_check = np.random.choice(np.flatnonzero(pr == pr.max()))
        port_check_t.append(port_check.copy())
        
        # IMPLEMENT SOFTMAX RULE HERE
        
        
        
        # call check_outcome to give outcome of check. Update RC array
        outcome = check_outcome(avail_r,port_check)
        RC[port_check] += outcome
        RC_t.append(RC.copy())
        
        # if reward was delivered and port checked, collect reward 
        if outcome:
            avail_r[port_check] = 0
        
        # update pr based on check outcome 
        pr = update_pr(policy,pr,outcome,port_check,t,RC,decay)
        pr_t.append(pr.copy())

    return(pr_t,RC_t,port_check_t,avail_r_t)



# support functions 
def check_outcome(avail_r,port):
    'Determine whether check is rewarded'
    outcome = avail_r[port]
    return(outcome)
    

def update_avail_rr(avail_r,rrs,t):
    'Update reward availability at each port'
    # check if t is a multiple of the given rr. If it is, make reward available at that port 
    for i,rr in enumerate(rrs):
        if (t/rr).is_integer():
            avail_r[i] = 1
        else:
            avail_r[i] = 0
    return(avail_r)


def update_pr(policy,pr,outcome,port,t,RC,decay_factor):
    'Update pr based on policy'
    updated_pr = pr
    
    if policy == 'random':
        # random policy: pr does not update based on reward
        updated_pr = pr
        
    elif policy == 'win_stay':
        # win stay: reward increases probability of staying at port to 1, then decays at some rate to random 
        # normalise across all ports. gamma * distance. Epsilon for exploration (separate)
        if outcome == 1:
            updated_pr[port] = 1
        elif outcome == 0:
            decayed_pr = pr[port]*decay_factor   
            if decayed_pr < 1/6:
                decayed_pr = 1/6
            updated_pr[port] = decayed_pr

    elif policy == 'matching':
        # naive rate estimate (total rewards at that port/total time)
        # local rate estimate vs. global rate estimate. Exponentially weight more recent rewards, one free parameter that can be adjusted 
        tot_rew          = RC[port]
        updated_pr[port] = tot_rew/t
        
   # elif policy == 'hazard':
        # as a function of time, likelihood of reward occurring 
        # cdf/ pdf
        # 
    
        
# Pr(t) = Hazard(t) (assuming an agent knows/correctly estimates the true hazard for each port) is an optimal strategy that matches estimates to the likelihood that a reward is available at a given port at time t.    
#cumulant and pdf (plus offset to avoid it going to zero/infinity)
        
        
    return(updated_pr)


# ADD THIS 
# def softmax(pr):
    
    
#     def softmax(x):
#     """Compute softmax values for each sets of scores in x."""
#         e_x = np.exp(x - np.max(x))
#         return e_x / e_x.sum()



def get_check_matrix(raw_data_folder,check_threshold,subject_ID,region,save_vid_csv=False):


    multi_lick_df = ldf.gen_multises_lick_df(raw_data_folder=raw_data_folder, on_only=True, check_threshold=check_threshold,
                     subject_ID=subject_ID, region=region, task ='conc', save_vid_csv=save_vid_csv, save_multises_csv=False)
    




# ordered list of visits 







# 1. An moment-to-moment state estimate of the probability of reward at all six ports. This is Pr(t) where Pr is a vector with as many entries as reward ports. 
# 2. A decision process that evaluates Pr(t) to render a current decision about which port is the current goal. So some process mapping Pr(t) -> G(t) where G is a scalar. 
# 3. A navigational process that takes the agent from its current location to G(t).


# Pr(t) = Hazard(t) (assuming an agent knows/correctly estimates the true hazard for each port) is an optimal strategy that matches estimates to the likelihood that a reward is available at a given port at time t.
# Pr(t,x,y) - could imagine formulations in which Pr is not only the hazard but also scaled by the distance of the agent from the reward port. 
# If Pr(t) is total rewards / time + distance this gets very close to MVT because it is the current rate penalized by some scaled version of travel cost (distance to ports from current location).

# ideal observer: noiseless, just know what the intervals are - but still have sampling rate of checking    

# what to plot from this model? cumulative rewards collected over time. Can fit 6 dimensional curve  
    
    
# agent's probability of checking is matched to the total number of checks that a mouse does - can query whether the agent wants to check at some frequency, or some probabilistic model that generates checks - can use actual distribution of checks
# can take histogram of all intercheck intervals, sample from histogram 

# integrating hazard functions over time, cross some threshold to check a port 
# compare 15s to 30s in terms of check rate in mouse data - if they're different then it's dependent on the environment, if not then there's some intrinsic rate'

# use montecarlo method to resample distribution - run simulation a bunch of times from empircial samples
# can split histogram into just rewarded versus just unrewarded (because after reward maybe on average they take more or less time to check another port) - conditioned distributions
# one back reward
# can do a GLM (montee carlos)    
    
# We can use empirical estimates of the inter-portcheck interval for each mouse as a baseline model “state” for how often the **P(t,id)** is evaluated and a decision to check **G(t)=softmax(P(t,id))** is rendered. 
# This can be done as a monte carlo method in which we re-sample the empirical distribution (with repetition) to generate probabilistic predictions about the cumulative checks per port over time as a function of different models of how **P(t,id)** is computed. 
# These empirical distributions and simulations can be computed from the matrix described above.



# Note: given that we have already observed some dependence upon spatial location of ports we could calculate probability normalized to distance to port more akin to a utility **U(t,id)=P(t,id)/Cost(id)** and then **G(t)=softmax(U(t,id))**.