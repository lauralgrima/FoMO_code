from pyControl.utility import *
import hardware_definition as hw

# parameters
v.session_duration = 3*hour 
v.port_delay = [[15000, 3000, 1000, 12000, 60000, 12000], [1200000,120000,600000,15000,120000,30000], [120000,600000,30000,15000,120000,1200000],[30000,600000,15000,1200000,120000,120000]]
v.within_session_switch = 5 # after X number of rewards, switch delays across ports 
v.reward_durations = [50, 50, 50, 50, 50, 50]
v.pwp = [2, 3, 4, 5, 6, 7] # pulse width percent 

# variables
v.n_rew_available = [0, 0, 0, 0, 0, 0]
v.n_rew_consumed  = [0, 0, 0, 0, 0, 0]
v.tot_rew = 0
v.n_switches = 0
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
	set_timer('session_timer', v.session_duration)  
	set_timer('port1_timer',v.port_delay[v.n_switches][0])
	set_timer('port2_timer',(v.port_delay[v.n_switches][1]+2*second))
	set_timer('port3_timer',(v.port_delay[v.n_switches][2]+4*second))
	set_timer('port4_timer',(v.port_delay[v.n_switches][3]+6*second))
	set_timer('port5_timer',(v.port_delay[v.n_switches][4]+8*second))
	set_timer('port6_timer',(v.port_delay[v.n_switches][5]+10*second))  
	hw.PWM_channel.pulse_width_percent(1)
	print(v.port_delay)

def run_end():
    # Turn off all hardware outputs.  
    print(v.tot_rew)
    hw.off()

def task_available(event):
	if event in ('lick1','lick2','lick3','lick4','lick5','lick6'):
		v.active_port = int(event[4])-1 
		hw.PWM_channel.pulse_width_percent(v.pwp[v.active_port])
		if v.n_rew_available[v.active_port] > 0: # if reward is available 
			v.n_rew_available[v.active_port] = 0 
			v.n_rew_consumed[v.active_port] += 1
			v.tot_rew += 1
			goto_state('reward')
	elif event in ('lick_off1','lick_off2','lick_off3','lick_off4','lick_off5','lick_off6'):
		hw.PWM_channel.pulse_width_percent(1)

def reward(event):
	if event == 'entry':
		hw.PWM_channel.pulse_width_percent(1)
		print('RA: {} RC: {} TR: {} S: {}'.format(v.n_rew_available, v.n_rew_consumed, v.tot_rew, v.n_switches)) 
		v.port_sols[v.active_port].on()
		timed_goto_state('iti', v.reward_durations[v.active_port])

def iti(event):
	if event == 'entry':
		v.port_sols[v.active_port].off()
		if v.within_session_switch > 0:
			switch_threshold = v.within_session_switch * (v.n_switches + 1)
			if v.tot_rew == switch_threshold:
				v.n_switches +=1
				print('switch')
				print(v.port_delay[v.n_switches])
				reset_timer('port1_timer',v.port_delay[v.n_switches][0])
				reset_timer('port2_timer',(v.port_delay[v.n_switches][1]+2*second))
				reset_timer('port3_timer',(v.port_delay[v.n_switches][2]+4*second))
				reset_timer('port4_timer',(v.port_delay[v.n_switches][3]+6*second))
				reset_timer('port5_timer',(v.port_delay[v.n_switches][4]+8*second))
				reset_timer('port6_timer',(v.port_delay[v.n_switches][5]+10*second))  
				timed_goto_state('task_available',200*ms)
			else: 
				timed_goto_state('task_available',200*ms)
		else:
			timed_goto_state('task_available',200*ms)

def all_states(event):
	if event in ('port1_timer','port2_timer','port3_timer','port4_timer','port5_timer','port6_timer'):
		v.active_timer = int(event[4])-1
		v.n_rew_available[v.active_timer] += 1
		print('RA: {} PD: {}'.format(v.n_rew_available,v.port_delay[v.n_switches]))
		timer_to_be_set = 'port'+str(v.active_timer+1)+'_timer'
		set_timer(timer_to_be_set,v.port_delay[v.n_switches][v.active_timer])
	elif event == 'session_timer':
		stop_framework()
