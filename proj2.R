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

# Rough proportions: Julia K - %, Akira - %, Julia J - %

# ------- Outline of project----------------------------------------------------

# We simulated the prisoner-box problem to estimate probabilities for 3 strategies.
# We constructed a function PONE to estimate the probability of one prisoner
# finding a card with their pre-assigned number inside shuffled boxes.
# We then built a function PALL which estimates the probability of all prisoners finding
# the cards with their pre-assigned number in shuffled boxes.
# We also constructed a function DLOOP that estimates probabilities of each loop length
# ranging from 1 to 2*n occurring at least once in a given number of simulations.
# A plot is provided to visualise the probabilities. We find the probability that 
# there is no loop longer than 50 using the previously written DLOOP function for n=50.

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
        if(p_success == prisoners){    # Count simulations with all prisoners successfull
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


dloop <- function(n,nreps){
  
  # The function estimates the probability of the each loop length between 1 and 2*n
  # occuring at least once in nreps simulations. It returns a vector of probabilities.
  
  prisoners <- 2*n
  all_lengths_prob <- rep(0,prisoners) # We initialise the counts of each loop length appearing
  for (a in 1:nreps){
    boxes <- sample(prisoners)         # We simulate numbers inside boxes
    all_lengths <- rep(0,prisoners)    # We record which loop length has occurred
    for (p in 1:prisoners){
      person <- p 
      length <- 1                      # We record the loop length
      for (k in 1:prisoners){
        if (second == person){
          all_lengths[length] <- 1
          length <- 0
          break
        }
        second <- boxes[second]       # We find the sequence of cards
        length <- length + 1          # We increase the loop length
      }
    } 
    all_lengths_prob <- all_lengths_prob + all_lengths # We sum the loop lengths that appeared across the simulations
  }
  all_lengths_prob/nreps              # Probability estimate
}


estiamtes <- dloop(50)                # We find the probabilities for n=50
p_no_over_fifty <- 1-sum(estimates[51:100]) # We find the probability of no loop longer than 50
plot(1:100, estimates, main = "Probabilities of each loop length appearing at least once for n=50 and nreps=10000",
     pch = 20,xlab="Length", ylab="Probability")
                                      # We plot the probabilities for each loop length
