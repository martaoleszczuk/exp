

from psychopy import visual, event, core
import random
import numpy as np
# Experiment instructions




def main_instructions_ieeg(win):
       
    inst = visual.TextStim(win, pos = [0,3])
    inst.wrapWidth = 20
    inst.text = "In this experiment your task is to remember the position of a stimulus that will briefly appear around the fixation point. After a while it will disappear, and by moving the mouse or looking in the right direction you will try to adjust the position of a red mark to match the remembered position of the stimulus. When you are more or less sure of the selected position, click with the mouse and you will see by how many degrees of angle your answer differs."
    
    inst.height = 1.0
    nextt = visual.TextStim(win, pos = [0,-8])
    nextt.wrapWidth = 20
    nextt.height = 0.7
    nextt.color = 'black'
    nextt.text = "Press space to continue "
    inst.draw()
    nextt.draw()
    win.flip()
    event.waitKeys(keyList = ["space"])

    inst.text = "IT IS VERY IMPORTANT TO TRY TO KEEP YOUR EYES ON THE FIXATION POINT DURING EACH TRIAL! IF THE EYE-TRACKER DETECTS THAT YOU ARE NOT LOOKING AT THE CENTRE, THE TRIAL WILL BE REPEATED. At the end of each block you can rest as much as you need. Also if you are tired you can move your eyes a little at the end of each trial. The experiment has a duration of 10 blocks of approximately 6 minutes each. Take it easy and thank you for participating!"
    inst.draw()
    nextt.draw()
    win.flip()
    event.waitKeys(keyList = ["space"])
    
    return



def main_instructions(win):
       
    inst = visual.TextStim(win, pos = [0,3])
    inst.wrapWidth = 20
    inst.text = "In this experiment your task is to remember the position of a stimulus that will briefly appear around the fixation point. After a while it will disappear, and by moving the mouse or looking in the right direction you will try to adjust the position of a red mark to match the remembered position of the stimulus. When you are more or less sure of the selected position, click with the mouse and you will see by how many degrees of angle your answer differs."
    
    inst.height = 1.0
    nextt = visual.TextStim(win, pos = [0,-8])
    nextt.wrapWidth = 20
    nextt.height = 0.7
    nextt.color = 'black'
    nextt.text = "Press space to continue "
    inst.draw()
    nextt.draw()
    win.flip()
    event.waitKeys(keyList = ["space"])

    inst.text = "Practice: \
    As well as remembering the position, you will need to pay attention to the fixation point because it will invariably blink subtly.\
    IF YOU DETECT THAT THE FIXING POINT HAS FLASHED, instead of reporting the orientation of the stimulus, IN THE RESPONSE PHASE YOU SHOULD PRESS THE SPACEBAR.\
    When you report your answer, you will receive feedback on whether you have made a mistake."
    inst.height = 0.7
    nextt = visual.TextStim(win, pos = [0,-8])
    nextt.wrapWidth = 20
    nextt.height = 0.7
    nextt.color = 'black'
    nextt.text = "Press space to continue "
    inst.draw()
    nextt.draw()
    win.flip()
    event.waitKeys(keyList = ["space"])

    
    inst.text = "Practice: \
    IT IS VERY IMPORTANT TO TRY TO KEEP YOUR EYES ON THE FIXATION POINT DURING EACH TRIAL! IF THE EYE-TRACKER DETECTS THAT YOU ARE NOT LOOKING AT THE CENTRE, THE TRIAL WILL BE REPEATED. \
    At the end of each block you can rest as much as you need. Also if you are tired you can move your eyes a little at the end of each trial. \
        \
    The experiment has a duration of 10 blocks of approximately 4 minutes each. Take it easy and thank you for participating!"
    inst.draw()
    nextt.draw()
    win.flip()
    event.waitKeys(keyList = ["space"])
    
    return

def block_start(win, nblock, blocklen, cond_resps):
    inst = visual.TextStim(win, pos = [0,-4])
    inst.wrapWidth = 20
    inst.height = 0.7
    inst.color = 'black'
    nextt = visual.TextStim(win, pos = [0,-8])
    nextt.height = 0.7
    nextt.color = "black"
    
    mod = visual.TextStim(win, pos = [0, 5])
    mod.height = 1.2
    mod.color = "white"
    
    inst.text = 'Starting the block ' + str(nblock+1) + ' / ' + str(blocklen)
    if cond_resps == 'drag':
        mod.text = 'In the following block you will be responding with the mouse. Look at the fixation point until the red dot appears - then respond by moving the mouse to the remembered position and clicking when you are ready.'
    else:
        mod.text = 'In the following block you will be responding with the eyes. Look at the fixation point until the red dot appears - then respond by looking in the remembered position and clicking the mouse when you are ready.'
        
    nextt.text = "When you are ready, press space to start "
    inst.draw()
    nextt.draw()
    mod.draw()
    win.flip()
    event.waitKeys( keyList=['space'])
    core.wait(1)
    return

    
def new_trial(win):
    inst = visual.TextStim(win, pos = [0,0])
    inst.wrapWidth = 20
    inst.height = 0.7
    inst.color = 'white'
    inst.text = "New trial"
    inst.draw()
    win.flip()
    core.wait(0.75)
    return



    

def block_ID(win, rep, autodraw):
    inst = visual.TextStim(win, pos = [0,8])
    inst.wrapWidth = 20
    inst.height = 0.9
    inst.color = 'red'
    if rep == 'repeat':
        inst.text = "Repeated sequences"
    else:
        inst.text = "Similar sequences"
    inst.autoDraw = autodraw   
    inst.draw()  
    
    
    


def end_experiment(win):
    inst = visual.TextStim(win, pos=[0,0], height = 1.2)
    inst.text = "End of the experiment! Notify the researcher"     
    inst.draw()
    win.flip()
    event.waitKeys(keyList = ["space"])
    return


def end_experiment_lot(win, corr_lotery, nblocks):
    inst = visual.TextStim(win, pos=[0,0], height = 1.2)
    inst.text = 'End of this part of the experiment!! The participant has earned ' + str(sum(corr_lotery)) + ' points out of ' + str(nblocks) + '. Notify the researcher'     
    inst.draw()
    win.flip()
    event.waitKeys(keyList = ["space"])
    return

    
    

def lotery(win, block, ifi):
    inst = visual.TextStim(win, pos=[0,0], height = 1.2)
    inst.text = "Loteria de trials. Pulsa -espacio- para comenzar el sorteo."     
    inst.draw()
    win.flip()
    event.waitKeys(keyList = ["space"])
    allKeys = []
    number_text = visual.TextStim(win, pos=[0,4], height = 2)
    inst.text = "Loteria de trials. \
    Pulsa -espacio- para parar el sorteo."     
    n_trials = len(block['data'])
    
    time_number = round(100/ifi) # time in frames that each lotery number is displayed on the screen
    
    
    allKeys = []
    while allKeys != ['space']:
        #angle = np.random.uniform(0,1,1)
        trial_ix =  random.randint(0,n_trials-1)
        corr = block['data'].loc[trial_ix,'correct']
        col = 'red' if corr == -1 else 'green'
        number_text.text = trial_ix
        number_text.color = col
        
        for iframe in  range(time_number):
             inst.draw()
             number_text.draw()
             win.flip()
        
        allKeys = event.getKeys()
        
    allKeys = []
    
    inst.text = "Genial, has ganado! Pulsa -espacio- para continuar"
    if corr == -1:
        inst.text = "Oh! Qué mala suerte. Pulsa -espacio- para continuar"
        

    inst.draw()
    number_text.draw()
    win.flip()
        
    event.waitKeys(keyList = ["space"]) 
    return corr 
        
        
        