// markov-localization.cpp
// a discrete bayes filter for markov localization
// madeline lebreton, learning from 4AX3 notes Dr. Mohrenschildt

// define the map of green and red
PG = [0 1/3 1/3 0 0 0 1/3 0 0 0] // likelihood of seeing green at each location
PR = [1/7 0 0 1/7 1/7 1/7 0 1/7 1/7 1/7]

// state X is a 1x10 vectors rep. belief of where robot is
X0 = [.1 .1 .1 .1 .1 .1 .1 .1 .1 .1] // uniform prior, start with zero knowledge

// measurement update function using bayes' rule
function P = update(P1, P2)
    T = times(P1, P2); // element-wise multiplication of current belief by sensor observation likelihood
    P=T/sum(T); // divide by total sum to normalize
end 

// first measurement: sense green
X1 = update(X0, PG) // -> [0 1/3 1/3 0 0 0 1/3 0 0 0] after normalization

// motion update: shift the probability dis. array to account for movement
X2 = circshift(X1, [0,1]) // shift to the right by 1 index (arguments [row, column])
//X2 is the estimate now

// see green again:
X2 = update(X2, PG)
// now we know exactly where we are, because the probability is 1 at x2
// [0 0 1 0 0 0 0 0 0 0]
