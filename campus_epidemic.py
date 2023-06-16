# Code accompnaying paper by Best & Singh (Unpublished), "Comparing intervention measures in a model of a disease outbreak on a university campus"


### IMPORT LIBRARIES
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from random import *

# INITIALISE LISTS TO STORE OUTPUT
bigstore=[] # Used to store peaks
bigstore2=[] # Used to store totals
bigstorer0=[] # Used to store R0s

### GLOBAL CONSTANTS
POPSIZE=1000         # Total individuals in population
NOHOUSES=100         # No. of households (mean household size = POPSIZE / NOHOUSES)
REPS=100              # No. of replicate runs
I0=10                # Initial no. of infected individuals (mean)
MAX_TIME=180         # Max time
CLASS_START = 0.4    # Time in day when move to class 
CLASS_END =  0.48333 # Time in day when move to house

# EPIDEMIOLOGICAL CONSTANTS
OMEGA=1/7    # Latent rate
GAMMA=1/7    # Recovery rate
WANE=1/120   # Waning immunity rate

BETA=np.zeros(2)
BETA[0]=0.01 # Transmission coefficient at home
BETA[1]=0.01 # Transmission coefficient at class

# DEFAULT PARAMETERS TO VARY
MIXD=0.04    # Amount of extra-household mixing (keep low)
TRACINGD=5   # Max no. contact traced from a positive case
TESTNOD=200 # Proportion weekly testing
CLASSD=100    # Class sizes
ISOLATED=10   # Days isolating

# -----------------------------------#
### FUNCTIONS

# Function to count numbers in each compartment
def count_type(type):
    return np.count_nonzero(PERSON[:,3] == type)

# Function to check the current scale
def findscale(home):
    # Check location
    loc = NOHOUSES if home==1 else NOCLASSES
    counts = np.zeros((loc, 4))
    icounts = np.zeros((loc, 4))
    icountss = np.zeros((loc,4))
    for j in range(0, loc):
        for c in range(1, 4):
            counts[j, c] = np.count_nonzero(np.logical_and(PERSON[:,home]==j, PERSON[:,3]==c)) 
        for c in [0,2]:
            icounts[j, c] = np.count_nonzero(np.logical_and(np.logical_and(PERSON[0:1000,home]==j, PERSON[0:1000,3]==c), PERSON[0:1000,4]==0))
            #icountss[j, c] = np.count_nonzero(np.logical_and(np.logical_and(PERSON[0:200,home]==j, PERSON[0:200,3]), PERSON[0:200,4]==0))
    # Sum of all rates
    #scale = np.sum(OMEGA*counts[:,1] + GAMMA*counts[:,2] + WANE*counts[:,3]) + (1-MIX)*BETA[home-1]*np.sum((icounts[:,0]+icountss[:,0])*icounts[:,2])+(1-MIX)*BETASS*np.sum((icounts[:,0]+icountss[:,0])*icountss[:,2])+MIX*(BETA[home-1]*np.sum(icounts[:,2])+BETASS*np.sum(icountss[:,2]))*np.sum(icounts[:,0]+icountss[:,0])
    scale = np.sum(OMEGA*counts[:,1] + GAMMA*counts[:,2] + WANE*counts[:,3]) + (1-MIX)*BETA[home-1]*np.sum((icounts[:,0])*icounts[:,2])+MIX*(BETA[home-1]*np.sum(icounts[:,2]))*np.sum(icounts[:,0])

    return scale

def test(ct):
    # print(ct)
    # Only non-isolators test
    findx=np.where(PERSON[:,4]==0)
    # Adjust test number if too many isolating
    testnum = TESTNO if len(findx[0])>TESTNO else len(findx[0])
    # Randomise list
    shuffle(findx[0])
    # Isolate Exposed and Infected individuals for 10 days
    for tests in range(testnum):
        if PERSON[findx[0][tests],3]==1 or PERSON[findx[0][tests],3]==2:
            # False negative / ignore rate of 5%
            if np.random.uniform()<0.95:
                PERSON[findx[0][tests],4]=ISOLATE
                # contact trace
                num_traced = np.random.randint(TRACING+1)
                case_contacts=np.where((PERSON[:,1]==PERSON[findx[0][tests],1]) | (PERSON[:,2]==PERSON[findx[0][tests],2]))
                shuffle(case_contacts[0])
                #print(num_traced,case_contacts[0].size)
                for trace in range(num_traced):
                    if PERSON[case_contacts[0][trace],4]==0:
                        PERSON[case_contacts[0][trace],4]=ISOLATE
        else:
            # False positive rate of 1%
            if np.random.uniform()>0.99:
                PERSON[findx[0][tests],4]=ISOLATE
                # contact trace
                num_traced = np.random.randint(TRACING+1)
                case_contacts=np.where((PERSON[:,1]==PERSON[findx[0][tests],1]) | (PERSON[:,2]==PERSON[findx[0][tests],2]))
                shuffle(case_contacts[0])
                for trace in range(num_traced):
                    if PERSON[case_contacts[0][trace],4]==0:
                        PERSON[case_contacts[0][trace],4]=ISOLATE

solutionstore1=np.zeros((REPS,MAX_TIME+1))
solutionstore2=np.zeros((REPS,MAX_TIME+1))

# MAIN FUNCTION
for case in range(6):
    
    # Initialise stores
    peaks=[] # For storing peak no. of infecteds
    tots=[]  # For storing total no. infected
    R0store=np.zeros(REPS)
    
    # Set all but one to default, vary the other
    # The defaults are in lines 31-35
    MIX=MIXD          
    CLASS=CLASSD
    TESTNO=TESTNOD
    TRACING=TRACINGD
    ISOLATE=ISOLATED
    
    # Use these lines to define your cases
    # If varying one variable it can be done simply as, e.g. MIX=0.01*case
    TESTNO=200*case 
        
    classes=CLASS
    NOCLASSES=int(np.ceil(POPSIZE/(classes))) 

    # Main section of code
    for reps in range(0,REPS):
        print(reps)
        DTstore=0
        DTcheck=0
        shuff1=np.arange(POPSIZE)
        shuff2=np.arange(POPSIZE)

        shuffle(shuff1)
        shuffle(shuff2)

        PERSON=np.zeros((POPSIZE,5))
        # 0th col: ID - First 200 are superspreaders (if used)
        # 1st col: house
        # 2nd col: class
        # 3rd col: SIR status
        # 4th col: Isolation status
        infs=np.sort(np.random.choice(POPSIZE,I0))
        for i in range(0,POPSIZE):
            PERSON[i][0] = i
            if i in infs:
                PERSON[i][3] = 2 # Initially everyone susceptible except I0 individuals
                    
        for index in range(NOHOUSES):
            for i in shuff1[index::NOHOUSES]:
                PERSON[i][1]=index

        for index in range(NOCLASSES):
            for i in shuff2[index::NOCLASSES]:
                PERSON[i][2]=index 
                
        # Some local constants / lists
        tsteps=[0]
        infecteds=[count_type(2)]
        susceptibles=[count_type(0)]
        exposed=[count_type(1)]
        current_t=0  
        home=2              # Everyone starts at home
        total_infections=10 # Since immunity wanes, count every infection event
        
        test(current_t) # Run test at t=0 or R0 values are skewed

        # Main run
        while current_t < MAX_TIME and np.any(PERSON[:,3] == 2):

            # Find proposed time to next event
            scale = findscale(home)
            dt = -np.log(np.random.uniform()) / scale
            proposed_t = tsteps[-1] + dt
            
            # If new day, change isolation status and test, then start again
            if int(proposed_t) > int(current_t): # Has to come first or misses days
                # Update isolation status
                for ppl in range(POPSIZE):
                    if PERSON[ppl,4]>=1:
                        PERSON[ppl,4]-=1
                
                # Jump is only one new day
                current_t = int(current_t)+1
                
                # Set solutionstore for plotting
                if case==0:
                    solutionstore1[reps,current_t]=infecteds[-1]
                elif case==5:
                    solutionstore2[reps,current_t]=infecteds[-1]
                
                # Weekly Test
                if int(current_t)%7==0:
                    test(current_t)

            elif home == 1 and proposed_t > int(proposed_t) + CLASS_START and proposed_t < int(proposed_t) + CLASS_END:
                # If students are home and proposed time of next event is later 
                # than class starts, then no event occurs before class starts
                current_t = int(proposed_t)+CLASS_START
                home = 2
            elif home == 2 and proposed_t > int(proposed_t) + CLASS_END:
                # If students are in class and proposed time of next event is 
                # later than class ends, then no event occurs before class ends
                current_t = int(proposed_t)+CLASS_END
                home = 1 
            else:
                # Next event occurs before class starts/ends
                current_t = proposed_t
                                        
                # Find event
                eventcheck = np.random.uniform()
                # Need to find non-isolating infecteds, superspreaders and susceptibles
                countSS=np.count_nonzero(np.logical_and(PERSON[0:1000,3]==2, PERSON[0:1000,4]==0)) 
                #countNSS=np.count_nonzero(np.logical_and(PERSON[200:1000,3]==2, PERSON[200:1000,4]==0)) 
                countsus=np.count_nonzero(np.logical_and(PERSON[0:1000,3]==0, PERSON[0:1000,4]==0))
                
                if eventcheck < GAMMA*infecteds[-1]/scale: #Event is recovery  
                    
                    # If there are any infected people, randomly choose one to recover
                    infected_indices = np.where(PERSON[:,3] == 2)
                    if infected_indices:
                        choice=np.random.choice(infected_indices[0])
                        PERSON[choice, 3] = 3
                        if choice in infs:
                            infs=infs[infs!=choice]

                elif eventcheck < (GAMMA*infecteds[-1] + OMEGA*exposed[-1])/scale: # Event is latent->infected

                    # If there are any latents, randomly choose one to become infected
                    latent_indices = np.where(PERSON[:,3] == 1)
                    if latent_indices:
                        PERSON[np.random.choice(latent_indices[0]), 3] = 2                                                                                
                    
                elif eventcheck < (GAMMA*infecteds[-1] + OMEGA*exposed[-1]+WANE*count_type(3))/scale: 
                    # Event is waned immunity
                    # If there are any immunes, randomly choose one to become susceptible
                    immune_indices = np.where(PERSON[:,3] == 3)
                    if immune_indices:
                        PERSON[np.random.choice(immune_indices[0]), 3] = 0                                 
                
                elif eventcheck <(GAMMA*infecteds[-1] + OMEGA*exposed[-1]+WANE*count_type(3) + MIX*(BETA[home-1]*countSS)*countsus)/scale:
  
                    # Event is mixed transmission
                    
                    sus_indices = np.where((PERSON[:,3] == 0) & (PERSON[:,4]==0))
                    if sus_indices[0].size>0:
                        PERSON[np.random.choice(sus_indices[0]), 3] = 1
                        total_infections+=1
                        # Check if founder made infection
                        currentinfs=np.where((PERSON[infs,3]==2) & (PERSON[infs,4]==0))
                        totalinfs=np.where((PERSON[:,3]==2) & (PERSON[:,4]==0))
                        if totalinfs[0].size>0:
                            if np.random.rand() < currentinfs[0].size / totalinfs[0].size:
                                R0store[reps]+=1/I0
                                #print(current_t)
                
                else: #Event is Transmission by househould contact
                    findx=np.where((PERSON[:,3]==0) & (PERSON[:,4]==0)) # Find susceptible hosts who are not isolating
                    shuffle(findx[0]) # randomly shuffle
                    for tryx in findx[0]:
                        loc=PERSON[tryx,home]
                        contacts=np.where(np.logical_and(np.logical_and(PERSON[:,home]==loc, PERSON[:,3]==2), PERSON[:,4]==0))
                        #contacts=np.where(PERSON[:,home]==loc)
                        #infcontacts=0
                        #for c in contacts[0]:
                        #    # Find infectious non-isolating contacts
                        #    if PERSON[c,4]==0 and PERSON[c,3]==2:
                        #        infcontacts+=1
                        if contacts[0].size>0:
                            PERSON[tryx,3]=1
                            total_infections+=1
                            # Check if founder made infection, if so add to R0
                            currentinfs=np.where((PERSON[infs,3]==2) & (PERSON[infs,4]==0))
                            currentinfs_ind=(infs[currentinfs[0]])
                            if currentinfs[0].size>0:
                                idx = np.searchsorted(currentinfs_ind,contacts[0])
                                idx[idx==len(currentinfs_ind)] = 0
                                mask = currentinfs_ind[idx]==contacts[0]
                                if np.random.uniform()<sum(mask)/contacts[0].size:
                                    R0store[reps]+=1/I0
                                    #print(current_t)
                            break  
            
            # Update lists
            tsteps.append(current_t)
            infecteds.append(count_type(2))
            exposed.append(count_type(1))
            susceptibles.append(count_type(0))
            

            # Stop if infections has finished
            if infecteds[-1]==0 & exposed[-1]==0:
                break

        # Find peak no. infected
        peaks.append(max(infecteds)/10)
        tots.append(total_infections)

    print("The median peak is", np.median(peaks))
    print("The median R0 is", np.median(R0store))
    print("The median total is", np.median(tots))
    
    bigstore.append(peaks)
    bigstore2.append(tots)
    bigstorer0.append(R0store)

totsq=np.zeros((6,6))
peaksq=np.zeros((6,6))
peaksr0=np.zeros((6,6))
for i in range(6):
    totsq[i,0]=i
    peaksq[i,0]=i
    peaksr0[i,0]=i
    for j in range(5):
        totsq[i,j+1]=np.quantile(bigstore2[:][i],(j)*0.25)
        peaksq[i,j+1]=np.quantile(bigstore[:][i],(j)*0.25)
        peaksr0[i,j+1]=np.quantile(bigstorer0[:][i],(j)*0.25)

print(peaksq)
print(totsq)
print(peaksr0)