# ------- Group member names----------------------------------------------------
# Julia Kaczmarczyk - s1977294 (Julia K)
# Akira Ishiyama - s2445245
# Julia Jose - s2421229 (Julia J)

# -------- Individual Contributions --------------------------------------------
# Creation of repo (Julia J)
# Split_punct function (Julia K)
# Vector of most common unique words (Julia J and Akira)
# Creating matrices T, A and S (Julia J and Akira)
# Building model to simulate text (Julia K)
# Capitalizing words (Julia K and Akira)

# Rough proportions: Julia K - 37%, Akira - 30%, Julia J - 33%

# ------- Outline of project----------------------------------------------------

# We simulate 50 words of text from the Bible using 2nd order Markov Text model.
# We read the text into a vector, separate the words and punctuation marks, and 
# choose the 500 most common words. We calculate frequencies for unigrams,
# bigrams and trigrams and use it as probabilities to generate our simulated text. 
# We compare our 2nd order Markov model with one that uses only the unigrams.

setwd("C:\\Users\\kaklu\\StatisticalProgramming") # User needs to change the directory
a <- scan("pg10.txt",what="character",skip=104)
n <- length(a)
a <- a[-((n-2886):n)]
a <- a[-grep("[0123456789]:[0123456789]",a)]


# -------- Task 4: Creating the split_punct function ---------------------------

split_punct <- function(x) {
  # The function takes a vector x of words and punctuation joined together and 
  # outputs a vector x with words and punctuation separated. Firstly, it puts 
  # a £ sign between the word and the puctuation sign (,.;!:?) with sub() function, 
  # then it attaches £ sign at the end of each word, collapses all the words 
  # together and splits the string by £ signs.
  x <- sub('(.)([,.;!:?])','\\1£\\2',x) |>
    paste('£',sep="") |> 
    paste(collapse="") |>
    strsplit('£')
  x <- x[[1]] # We want the first element of the 2d array
}

# ------- Task 5: Using the function -------------------------------------------

a <- split_punct(a)

# ------- Task 6:Creating a vector of common unique words ----------------------

# (a) # Determine unique words
unique_words <- unique(tolower(a))

# (b) # Index each word w.r.t unique_words
index_vector <- match(tolower(a),unique_words)

# (c) # Find frequency of each word
tab <- tabulate(index_vector)

# (d) # Most frequent words
last <- tail(head(sort(tab,decreasing=TRUE),n=500),1)

# (e) # List of top 500 words
b <- (unique_words[tab>=last]) # On visual inspection of the output, it is exactly 500 words

# -------- Task 10: Capitalising the words -------------------------------------

new_index_vector <- match(a,unique_words) # Index of words with NA as uppercase words
new_tab <- tabulate(new_index_vector) # Count of the lowercase words
upper_probab <-(tab-new_tab)/tab # The proportion of the uppercase words overall
upper <- unique_words[upper_probab>0.5] # We choose the words that have the proportion of uppercase greater than 0.5
b_upper <- b[b%in%upper] # We compare the words with vector b
substr(b_upper, 1, 1) <- toupper(substr(b_upper, 1, 1)) # We capitalise the words
b[b%in%upper] <- b_upper # We substitute the lowercase words with uppercase

# -------- Task 7: Generating matrix T -----------------------------------------

# (a) # Index each word w.r.t 500 most common unique words
index_vector_freq <- match(tolower(a),tolower(b))

# (b) # The index vector shifted by one place, two places and matching the length by appending NA
vec2 <- index_vector_freq[2:length(index_vector_freq)] |> append(NA)
vec3 <- index_vector_freq[3:length(index_vector_freq)] |> append(NA) |> append(NA)
col_matrix <- cbind(index_vector_freq,vec2,vec3)

# (c) # Pulling out rows with no NA
col_matrix_T <- col_matrix[is.na(rowSums(col_matrix))== FALSE,1:3]

# (d) # Matrix T (word triplets/trigrams)
T_matrix <- array(0,dim=c(500, 500,500))

# Common word triplets, adding a 1 every time the kth word follows the pair i,j
for (row in 1:length(col_matrix_T[,1])){
  current_row <- col_matrix_T[row,]
  i <- current_row[1]
  j <- current_row[2]
  k <- current_row[3]
  T_matrix[i,j,k] <- T_matrix[i,j,k]+1
}

# (f) # Matrix A (words pairs/bigrams); following the same procedure

col_matrix_A <- cbind(index_vector_freq,vec2)
col_matrix_A <- col_matrix_A[is.na(rowSums(col_matrix_A))== FALSE,]
A_matrix <- array(0,dim=c(500,500))

for (row in 1:length(col_matrix_A[,1])){
  current_row <- col_matrix_A[row,]
  i <- current_row[1]
  j <- current_row[2]
  A_matrix[i,j] <- A_matrix[i,j]+1
}

# Vector S (single words/Unigrams)
vector_S <- tabulate(index_vector_freq)

# ------------ Task 8: Model ---------------------------------------------------
index_text <- rep(0,50)
first <- sample(500, size=1, prob=vector_S) # Finding first word using vector S
index_text[1] <- first

# Finding second word using matrix A, falling back on S when needed
if (sum(A_matrix[first,])>0){
  second <- sample(500, size=1, prob=A_matrix[first,])
} else {
  second <- sample(500, size=1, prob=vector_S)
}
index_text[2] <- second

# Finding consecutive 48 words;Using T first, if no frequency in T, falling back 
# on A and S as needed

for (i in 3:50){
  if (sum(T_matrix[first,second,])>0){
    index_text[i] <- sample(500, size=1, prob=T_matrix[first,second,])
  } else{
    if (sum(A_matrix[second,])>0){
      index_text[i] <- sample(500, size=1, prob=A_matrix[second,])
    } else{
      index_text[i] <- sample(500, size=1, prob=vector_S)
    }
  }
  first <- second
  second <- index_text[i]
}

# Output using matrix T,A and S
cat(b[index_text])

# ------------ Task 9: Comparing with simpler model using only vector S --------

cat(b[sample(500, size=50, prob=vector_S)])
