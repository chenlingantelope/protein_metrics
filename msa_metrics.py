import numpy as np
# HMM likelihood of sequence to HMM alignment
def hmm_likelihood(seq, hmm):
    j = -1
    ll = 0
    state = 'e'
    # if x0 is an emission state
    for i in range(len(seq)):
        x = seq[i]
        last_state = state
        if i == 0:
            step = hmm.start_step
        elif last_state !='i':
            step = hmm.steps[j]
        if x.isupper():
            state = 'e'
            ll += np.log10(step.p_emission_char[x])
            if last_state == 'e':
                ll += np.log10(step.p_emission_to_emission)
            elif last_state == 'd':
                ll += np.log10(step.p_deletion_to_emission)
            elif last_state == 'i':
                ll += np.log10(step.p_insertion_to_emission)
            j += 1
        elif x == '-':
            state = 'd'
            if last_state == 'e':
                ll += np.log10(step.p_emission_to_deletion)
            elif last_state == 'd':
                ll += np.log10(step.p_deletion_to_deletion)
            j += 1
        elif x.islower():
            state = 'i'
            x = x.upper()
            if j < len(hmm.steps):
                ll += np.log10(step.p_insertion_char[x])
                if last_state == 'e':
                    ll += np.log10(step.p_emission_to_insertion)
                elif last_state == 'i':
                    ll += np.log10(step.p_insertion_to_insertion)
        else:
            print("Invalid character in sequence")
    return ll
