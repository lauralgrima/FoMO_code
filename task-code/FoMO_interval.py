from pyControl.utility import *
import hardware_definition as hw


# parameters
v.session_duration = 3*hour 
v.port_delay = [30*second,60*second,240*second,240*second,1200*second,2400*second]
v.reward_durations = [50, 50, 50, 50, 50, 50]

# variables
v.n_rew_available = [0, 0, 0, 0, 0, 0]
v.n_rew_consumed  = [0, 0, 0, 0, 0, 0]
v.tot_rew = 0
v.active_port  = -1 
v.active_timer = -1
v.port_sols = [hw.port1.SOL, hw.port2.SOL, hw.port3.SOL, hw.port4.SOL, hw.port5.SOL, hw.port6.SOL]

states = ['task_available',
		  'iti',
		  'reward']

events = ['rsync',
		  'session_timer',
		  'lick1','lick_off1','port1_timer',
		  'lick2','lick_off2','port2_timer',
		  'lick3','lick_off3','port3_timer',
		  'lick4','lick_off4','port4_timer',
		  'lick5','lick_off5','port5_timer',
		  'lick6','lick_off6','port6_timer']


initial_state = 'task_available'

def run_start():  
	# Start off all timers
	set_timer('port1_timer',v.port_delay[0])
	set_timer('port2_timer',(v.port_delay[1]+2*second))
	set_timer('port3_timer',(v.port_delay[2]+4*second))
	set_timer('port4_timer',(v.port_delay[3]+6*second))
	set_timer('port5_timer',(v.port_delay[4]+8*second))
	set_timer('port6_timer',(v.port_delay[5]+10*second))
	set_timer('session_timer', v.session_duration)  
	print(v.port_delay)

def run_end():
    # Turn off all hardware outputs.  
    print(v.tot_rew)
    hw.off()

def task_available(event):
	if event in ('lick1','lick2','lick3','lick4','lick5','lick6'):
		v.active_port = int(event[4])-1 
		if v.n_rew_available[v.active_port] > 0: # if reward is available 
			v.n_rew_available[v.active_port] = 0 
			v.n_rew_consumed[v.active_port] += 1
			v.tot_rew += 1
			goto_state('reward')

def reward(event):
	if event == 'entry':
		print('RA: {} RC: {} TR: {}'.format(v.n_rew_available, v.n_rew_consumed, v.tot_rew))
		v.port_sols[v.active_port].on()
		timed_goto_state('iti', v.reward_durations[v.active_port])

def iti(event):
	if event == 'entry':
		v.port_sols[v.active_port].off()
		timed_goto_state('task_available',200*ms)

def all_states(event):
	if event in ('port1_timer','port2_timer','port3_timer','port4_timer','port5_timer','port6_timer'):
		v.active_timer = int(event[4])-1
		v.n_rew_available[v.active_timer] += 1
		print('RA: {}'.format(v.n_rew_available))
		timer_to_be_set = 'port'+str(v.active_timer+1)+'_timer'
		set_timer(timer_to_be_set,v.port_delay[v.active_timer])
	elif event == 'session_timer':
		stop_framework()
