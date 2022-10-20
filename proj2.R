# ------- Group member names----------------------------------------------------
# Julia Kaczmarczyk - s1977294 (Julia K)
# Akira Ishiyama - s2445245
# Julia Jose - s2421229 (Julia J)

# ------- Github repo ----------------------------------------------------------
# https://github.com/JuliaJ99/StatisticalProgramming

# -------- Individual Contributions --------------------------------------------
# Pone (Julia J, Julia K, Akira)
# Pall (Julia J, Julia K, Akira)
# dloop (Julia J, Akira)
# Task 6 - probability (Julia J)
# Task 6 - plot (Julia K)

# Rough proportions: Julia K - 1/3, Akira - 1/3, Julia J - 1/3

# ------- Outline of project----------------------------------------------------

#    We simulated the prisoner-box problem to estimate probabilities for 3 strategies.
# Prisoners open a sequence of boxes of length 2*n to find their number inside. Under
# strategy 1 they start the sequence by opening a box with their number,
# under strategy 2 they start with a random box, in both cases they continue 
# the search by opening boxes according to cards they found in the previous one.
# Under strategy 3 they pick all boxes in random order.
#
#    We constructed a function Pone to estimate the probability of one prisoner
# finding a card with their pre-assigned number inside shuffled boxes.
# We then built a function Pall which estimates the probability of all prisoners finding
# the cards with their pre-assigned number in shuffled boxes.
#
#    We also constructed a function dloop that estimates probabilities of each loop length
# ranging from 1 to 2*n occurring at least once in a given number of simulations.
# A plot is provided to visualize the probabilities. We find the probability that 
# there is no loop longer than 50 using the previously written dloop function for n=50.

# ------- Code -----------------------------------------------------------------

Pone <- function(n,k,strategy=1,nreps=10000){
  
  # The function estimates the probability of a single prisoner succeeding in 
  # finding their card by opening a sequence of numbered boxes.
  # The function takes arguments n (number of picks a prisoner have), 
  # k (the prisoner’s number), strategy (1, 2 or 3) and nreps. nreps represents 
  # the number of simulations run in order to estimate the probabilities. 
  
  success <- 0
  prisoners <- 2*n
  person <- k
  if (strategy == 3){                   # Check if strategy 3
    for (a in 1:nreps){
      picks <- sample(prisoners,n)      # Sample all picks randomly
      if (person %in% picks){           # We count the success if prisoner number appears
        success = success +1
      }
    }
  } else {for (a in 1:nreps){
    boxes <- sample(prisoners)          # We simulate numbers inside boxes 
    if (strategy == 2){             
      k <- sample(prisoners,1)          # If strategy 2 chosen, assign first box at random
    }
    card <- boxes[k]                    # First card
    for (i in 1:n-1){
      if (person == card){
        success <- success + 1          # We count the success if prisoner number appears
        break                           # We break the loop if successful
      }
      card <- boxes[card]               # Next cards
    }
  }
  }
  success/nreps                         # Probability estimate
}


Pall <- function(n,strategy=1,nreps=10000){
  
  # The function estimates the probability of all prisoners succeeding in 
  # finding their cards by opening a sequence of numbered boxes.
  # The function takes arguments n (number of picks each prisoner have), 
  # strategy (1, 2 or 3) and nreps. nreps represents 
  # the number of simulations run in order to estimate the probabilities. 
  
  success <- 0
  prisoners <- 2*n
  if (strategy == 3){                  # Check if strategy 3
    for (a in 1:nreps){
      p_success <- 0                   # Store success for each simulation
      for (j in 1:prisoners){
        picks <- sample(prisoners,n)   # Sample all picks randomly
        if (!(j %in% picks)){          # Break the loop if a prisoner is unsuccessful
          break
        }
        else{
          p_success <- p_success + 1   # Count individual successes
        }
      }
      if (p_success == prisoners){     # Count simulations with all prisoners successfull
        success <- success + 1
      }
    }
  }
  else {
    for (a in 1:nreps){
      boxes <- sample(prisoners)       # We simulate numbers inside boxes
      p_success <- 0                   # Store success for each simulation
      for (j in 1:prisoners){
        first <- j                     # Assign prisoner number to first box, strategy 1 default
        if (strategy == 2){           
          first <- sample(prisoners,1) # If strategy 2 chosen, reassign a random number to first box
        }
        card <- boxes[first]           # Find first card
        for (i in 1:n-1){
          if (j == card){          
            p_success <- p_success + 1 # If prisoner number found, record individual success
            break
          }
          card <- boxes[card]          # Find the sequence of cards
        }
        if(p_success == prisoners){    # Count simulations with all prisoners successful
          success <- success + 1
        }
      }
    }
    
  }
  success/nreps                        # Probability estimate
}

# Strategy 1, n=5 and n=50
Pone(5, k=1)
Pone(50, k=1)
Pall(5)
Pall(50)

# Strategy 2:
Pone(5, k=1, strategy=2)
Pone(50,k=1,strategy=2)
Pall(5, strategy=2)
Pall(50, strategy=2)

# Strategy 3
Pone(5, k=1, strategy=3)
Pone(50, k=1, strategy=3)
Pall(5, strategy=3)
Pall(50, strategy=3)

# Observations:
 #   What surprises us of the estimated probabilities is that the probability
 # of all prisoners (Pall) successfully finding their cards asymptotically 
 # converges to 0.31 as n goes to infinity for strategy 1, while for strategies 2 
 # and 3 the probabilities converges to 0 as n grows larger.
 #   
 #   The results we observed for Pall is a stark contrast compared to the 
 # probabilities we observed for Pone. For strategies 1 and 3 in Pone, the 
 # estimated probabilities converge to 0.5, and 0.375 for strategy 2. Hence, we
 # expected Pall to follow the same pattern but it does not.
 #   
 #   While all three strategies were somewhat successful for Pone, only strategy
 # 1 proved to be useful for Pall.


dloop <- function(n,nreps){
  
  # The function estimates the probability of the each loop length between 1 and 2*n
  # occuring at least once in nreps simulations. It returns a vector of probabilities.
  
  prisoners <- 2*n
  all_lengths_prob <- rep(0,prisoners) # We initialize the counts of each loop length appearing
  for (a in 1:nreps){
    boxes <- sample(prisoners)         # We simulate numbers inside boxes
    all_lengths <- rep(0,prisoners)    # We record which loop length has occurred
    for (p in 1:prisoners){
      person <- p 
      length <- 1                      # We record the loop length
      card <- boxes[p]
      for (k in 1:prisoners){
        if (card == person){
          all_lengths[length] <- 1
          length <- 0
          break
        }
        card <- boxes[card]           # We find the sequence of cards
        length <- length + 1          # We increase the loop length
      }
    } 
    all_lengths_prob <- all_lengths_prob + all_lengths # We sum the loop lengths that appeared across the simulations
  }
  all_lengths_prob/nreps              # Probability estimate
}


estimates <- dloop(50,10000)           # We find the probabilities for n=50
p_no_over_fifty <- 1-sum(estimates[51:100]) # We find the probability of no loop longer than 50
p_no_over_fifty                        # We expect this estimated probability to be close to
                                       # the estimate we obtained for Pall where n is 50 and 
                                       # strategy is 1.
plot(1:100, estimates, main = "Probabilities of each loop length appearing at least once for n=50 and nreps=10000",
     pch = 20,xlab="Length", ylab="Probability")
                                       # We plot the probabilities for each loop length
