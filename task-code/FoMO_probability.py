from pyControl.utility import *
import hardware_definition as hw

# parameters
v.session_duration = 3*hour 
v.port_probs = [0.9,0.8,0.5,0.5,0.2,0.1]
v.reward_durations = [50, 50, 50, 50, 50, 50]
v.vis_time_thresh = 2*second
v.pwp = [2, 3, 4, 5, 6, 7] # pulse width percent 

# variables
v.n_rew_consumed  = [0, 0, 0, 0, 0, 0]
v.n_visits = [0, 0, 0, 0, 0, 0]
v.tot_rew = 0
v.active_port  = -1 
v.active_visit = 1 # toggle this when it's been long enough between licks to count as a new visit
v.prev_port = -1
v.port_sols = [hw.port1.SOL, hw.port2.SOL, hw.port3.SOL, hw.port4.SOL, hw.port5.SOL, hw.port6.SOL]

states = ['task_available',
		  'iti',
		  'reward']

events = ['rsync',
		  'session_timer',
		  'ili_timer', # inter lick interval timer
		  'lick1','lick_off1',
		  'lick2','lick_off2',
		  'lick3','lick_off3',
		  'lick4','lick_off4',
		  'lick5','lick_off5',
		  'lick6','lick_off6']

initial_state = 'task_available'

def run_start():  
	set_timer('session_timer', v.session_duration)  
	print(v.port_probs)

def run_end():
    # Turn off all hardware outputs.  
    print(v.tot_rew)
    hw.off()

def task_available(event):
	if event in ('lick1','lick2','lick3','lick4','lick5','lick6'):
		v.active_port = int(event[4])-1 
		reset_timer('ili_timer',v.vis_time_thresh) # start timer for inter lick interval
		if v.active_visit == 1 or v.prev_port!=v.active_port: # can count as a visit 
			v.n_visits[v.active_port] += 1
			outcome = withprob(v.port_probs[v.active_port]) # whether visit is rewarded 
			v.active_visit = 0 # reset to 0 
			v.prev_port = v.active_port
			if outcome == 1:
				v.n_rew_consumed[v.prev_port] += 1 # because prev port has been assigned to active port value 
				v.tot_rew += 1
				print('NV: {} RC: {} TR: {}'.format(v.n_visits,v.n_rew_consumed, v.tot_rew))
				goto_state('reward')
			else:
				print('NV: {} RC: {} TR: {}'.format(v.n_visits,v.n_rew_consumed, v.tot_rew))

def reward(event):
	if event == 'entry':
		v.port_sols[v.prev_port].on()
		timed_goto_state('iti', v.reward_durations[v.prev_port])

def iti(event):
	if event == 'entry':
		v.port_sols[v.prev_port].off()
		timed_goto_state('task_available',200*ms)

def all_states(event):
	if event == 'ili_timer': # if inter lick interval timer goes off (i.e. if the ili has passed) the visit can be counted as a new one
		v.active_visit = 1 
	elif event == 'session_timer':
		stop_framework()
