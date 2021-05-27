import random
import math
import copy
import matplotlib.pyplot as plt
'''
Import metadata
'''
def alphabet():
    with open('data/alphabet.csv') as f:
        return f.read().rstrip().split(",")
alphabet = alphabet()
def transitions():
    matrix = {}
    with open('data/letter_transition_matrix.csv') as f:
        lines = f.readlines()
        for i in range(len(lines)):
            matrix[alphabet[i]] = {}
        for i in range(len(lines)):
            row = lines[i].rstrip().split(",")
            for t in range(len(row)):
                if(float(row[t])!=0):
                    matrix[alphabet[t]][alphabet[i]] = math.log(float(row[t]))
                else:
                    matrix[alphabet[t]][alphabet[i]] = -20
        return matrix
def initial_probs():
    probs = {}
    with open('data/letter_probabilities.csv') as f:
        line = f.read().rstrip().split(",")
        for p in range(len(line)):
            probs[alphabet[p]] = math.log(float(line[p]))
        return probs

#transitions[i][j] gives the probability of moving to j from state i.
transitions = transitions()
probabilities = initial_probs()
with open('data/sample/ciphertext.txt') as f:
    sample_text = f.read()
with open('data/sample/ciphertext_breakpoint.txt') as f:
    sample_bp = f.read()
with open('data/sample/plaintext.txt') as f:
    plain_text = f.read()
'''
Given ciphertext y and cipher function f, compute and return the likelihood,
p(y|f).
'''
def get_log_likelihood(f, y):
    inverse = get_inverse(f)
    current_letter = inverse[y[0]]
    likelihood = probabilities[current_letter]
    for i in range(1, len(y)):
        next_letter = inverse[y[i]]
        likelihood += transitions[current_letter][next_letter]
        current_letter = next_letter
    return likelihood
'''
Given a cipher f mapping a plaintext string x to a ciphertext string
y, return the inverse of f mapping a ciphertext y to plaintext x.
'''
def get_inverse(f):
    inverse = {}
    for x in f:
        inverse[f[x]] = x
    return inverse

'''
Given a distribution f, randomly select two elements in f and swap their
mappings to generate f'.
'''
def get_random_next(f):
    f_prime = copy.copy(f)
    to_swap = random.sample(alphabet, 2)
    temp = f_prime[to_swap[0]]
    f_prime[to_swap[0]] = f_prime[to_swap[1]]
    f_prime[to_swap[1]] = temp
    return f_prime

'''
Apply cipher f to invert the ciphertext sequence y back to plaintext x.
'''
def apply_cipher(f, y):
    inv = get_inverse(f)
    ans = []
    for i in y:
        ans.append(inv[i])
    return ''.join(ans)

'''
Using cipher f to decode ciphertext, compute the accuracy of f.
'''
def decode_accuracy(ciphertext, plaintext, f):
    correct = 0
    total = len(ciphertext)
    inv = get_inverse(f)
    for i in range(total):
        if(inv[ciphertext[i]] == plaintext[i]):
            correct =correct +1
    return correct/total
'''
Given ciphertext, iteratively perform M-H algorithm to decode it.
'''
def decode(ciphertext: str, has_breakpoint: bool) -> str:
    if(has_breakpoint):
        n = len(ciphertext)
        best_likelihood = -float("inf")
        best_b = int(n/2)
        best_f1 = {}
        best_f2 ={}
        for j in range(12):
            b = int(n/2)
            permute = copy.copy(alphabet)
            random.shuffle(permute)
            f1 = dict(zip(alphabet, permute))
            random.shuffle(permute)
            f2 = dict(zip(alphabet, permute))
            #pick breakpoints, evaluate log likelihood after doing M-H to find
            # maximum log likelihood for some ciphers.
            for i in range(10000):
                #First perform an M-H step to update f1/f2 based on current b.
                f1_prime = get_random_next(f1)
                f2_prime = get_random_next(f2)
                current_ll1 = get_log_likelihood(f1, ciphertext[0:b])
                new_ll1 = get_log_likelihood(f1_prime, ciphertext[0:b])
                a = min(0, new_ll1-current_ll1)
                x = random.uniform(0, 1)
                if(x<math.exp(a)):
                    f1 = f1_prime
                current_ll2 = get_log_likelihood(f2, ciphertext[b:n])
                new_ll2 = get_log_likelihood(f2_prime, ciphertext[b:n])
                a = min(0, new_ll2-current_ll2)
                x = random.uniform(0, 1)
                if(x<math.exp(a)):
                    f2 = f2_prime
                #Then perform an M-H step to update b based on current f1/f2.
                b_prime = (int(random.uniform(-1, 1)*15)+b)%n
                #catch edge case where b_prime is not a real breakpoint
                while(b_prime==0 or b_prime ==n-1):
                    b_prime = (int(random.uniform(-1, 1)*15)+b)%n
                current_ll = get_log_likelihood(f1, ciphertext[0:b])+get_log_likelihood(f2, ciphertext[b:n])
                new_ll = get_log_likelihood(f1, ciphertext[0:b_prime])+get_log_likelihood(f2, ciphertext[b_prime:n])
                a = min(0, new_ll-current_ll)
                x = random.uniform(0, 1)
                if(x<math.exp(a)):
                    b = b_prime
                    current_ll = new_ll
                if(current_ll>best_likelihood):
                    best_b = b
                    best_f1 = f1
                    best_f2 = f2
                    best_likelihood = current_ll
        plaintext = apply_cipher(best_f1, ciphertext[0:best_b])+apply_cipher(best_f2, ciphertext[best_b:n])
        return plaintext
    else:
        best_f = {}
        best_ll = -float('inf')
        for i in range(12):
            permute = copy.copy(alphabet)
            random.shuffle(permute)
            f = dict(zip(alphabet, permute))
            n = len(ciphertext)
            for i in range(10000): # use 10000 iterations for now
                f_prime = get_random_next(f)
                testLikelihood = get_log_likelihood(f_prime, ciphertext)
                currentLikelihood = get_log_likelihood(f, ciphertext)
                a = min(0, testLikelihood-currentLikelihood)
                x = random.uniform(0, 1)
                if(x<math.exp(a)):
                    f = f_prime
                    currentLikelihood = testLikelihood
                if(currentLikelihood>best_ll):
                    best_ll = currentLikelihood
                    best_f = f
        plaintext = apply_cipher(best_f, ciphertext)
        return plaintext
    # Code for plotting graphs
    #likelihood_array = []
    #acceptance_array = []
    #accuracy_array = []
    #likelihood_per = []
    #        likelihood_array.append(testLikelihood)
    #        likelihood_per.append(testLikelihood/n)
    #        acceptance_array.append(1) # 1 in array means accepted
    #    else:
    #        likelihood_array.append(currentLikelihood)
    #        likelihood_per.append(testLikelihood/n)
    #        acceptance_array.append(0) # 0 in array means rejected.
    #    accuracy_array.append(decode_accuracy(ciphertext, plain_text, f))
    #plt.plot(likelihood_array)
    #plt.ylabel('Log likelihood')
    #plt.xlabel('Iteration number')
    #plt.show()
    #acceptance = []
    #acceptance_count = 0
    #for i in range(len(acceptance_array)):
    #    if(i<20):
    #        acceptance_count = acceptance_count + acceptance_array[i]
    #    else:
    #        acceptance.append(acceptance_count/20)
    #        acceptance_count +=acceptance_array[i] - acceptance_array[i-20]
    #plt.plot(acceptance)
    #plt.ylabel("Acceptance rate for T=20")
    #plt.xlabel("Iteration number")
    #plt.show()
    #plt.plot(accuracy_array)
    #plt.ylabel("Decoding Accuracy")
    #plt.xlabel("Iteration number")
    #plt.show()
    #plt.plot(likelihood_per)
    #plt.ylabel("Log likelihood per symbol")
    #plt.xlabel("Iteration number")
    #plt.show()

