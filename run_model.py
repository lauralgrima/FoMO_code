from model import state_machine


# variables
n_ports   = 6
rrs       = [30, 60, 240, 240, 1200, 2400] # reward schedules in seconds 
policy    = 'win_stay'
decay     = 0.8                            # decay factor for win stay, between 0 and 1
max_steps = 11000                          # max no. timesteps in seconds 


[port_probabilities,rewards_collected,ports_checked,rewards_available] = state_machine(policy,rrs,decay,max_steps,n_ports)