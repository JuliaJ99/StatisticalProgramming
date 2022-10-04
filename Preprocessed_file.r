setwd("C:\\Users\\julia\\Desktop\\Julia_Git_repo\\StatisticalProgramming")
  
start_time= Sys.time()
a <- scan("pg10.txt",what="character",skip=104)
n <- length(a)
a <- a[-((n-2886):n)]
a <- a[-grep("[0123456789]:[0123456789]",a)]

# Separate out punctuation
split_punct <- function(x) {
  x <- sub('(.)([,.;!:?])','\\1£\\2',x) |> 
    paste('£',sep="") |> 
    paste(collapse="") |>
    strsplit('£')
  x <- x[[1]]
}

a <- split_punct(a)
a[0:30]

# ------- #

# (a) # Determine unique words

unique_words <- unique(tolower(a))
unique_words[0:10]

#(b) # The index of each word in a w.r.t unique_words
index_vector <- match(tolower(a),unique_words)
index_vector

#(c) # Frequency of each word
tab <- tabulate(index_vector)
tab

#(d) # Most frequent words
tail(head(sort(tab,decreasing=TRUE),n=500))
###the word that occurs 163 times is the least of the 500 most common words.

#(e) # List of top 500 words
b <- (unique_words[tab>=163])
b 

# -------- Task 10: Capitalising the words --------------
new_index_vector <- match(a,unique_words) #NA are uppercase words
new_tab <- tabulate(new_index_vector) #count of the lowercase words
upper_probab <-(tab-new_tab)/tab #the proportion of the uppercase words overall
upper <- unique_words[upper_probab>0.5] #we choose the words that have the proportion of uppercase greater than 0.5
b_upper <- b[b%in%upper] #we compre the words with wector b
substr(b_upper, 1, 1) <- toupper(substr(b_upper, 1, 1)) #we capitalise the words
b[b%in%upper] <- b_upper #we substitute the lowercase words with uppercase

# -------- Task 7: Generating matrix T -------------------

# (a)
index_vector_freq <- match(tolower(a),tolower(b))
index_vector_freq[0:20]

vec2 <- index_vector_freq[2:length(index_vector_freq)] |> append(NA)
vec3 <- index_vector_freq[3:length(index_vector_freq)] |> append(NA) |> append(NA)

col_matrix <- cbind(index_vector_freq,vec2,vec3)

col_matrix_T <- col_matrix[is.na(rowSums(col_matrix))== FALSE,1:3]

# b may contain words dropped after NA, just using it for now
T_matrix <- array(0,dim=c(500, 500,500))
T_matrix

for (row in 1:length(col_matrix_T[,1])){
  current_row <- col_matrix_T[row,]
  i <- current_row[1]
  j <- current_row[2]
  k <- current_row[3]
  T_matrix[i,j,k] <- T_matrix[i,j,k]+1
}

T_matrix[105,2,107]

# --------- Creating matrix A (for pairs of words)---------

col_matrix_A <- cbind(index_vector_freq,vec2)
col_matrix_A

col_matrix_A <- col_matrix_A[is.na(rowSums(col_matrix_A))== FALSE,]

A_matrix <- array(0,dim=c(500,500)) #Should it be 500 or something else?
A_matrix

for (row in 1:length(col_matrix_A[,1])){
  current_row <- col_matrix_A[row,]
  i <- current_row[1]
  j <- current_row[2]
  A_matrix[i,j] <- A_matrix[i,j]+1
}

# ------------ Creating vector S
vector_S <- tabulate(index_vector_freq)

# ------------ Task 8: Model
index_text <- vector( "numeric" , 50 )
first <- sample(500, size=1, prob=vector_S)
index_text[1] <- first
if (sum(A_matrix[first,])>0){
  second <- sample(500, size=1, prob=A_matrix[first,])
} else {
  second <- sample(500, size=1, prob=vector_S)
}
index_text[2] <- second
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

cat(b[index_text])

# ------------ Task 9: Model simpler

cat(b[sample(1:500, size=50, prob=vector_S)])


# --------------Simulation
end_time <- Sys.time()
totaltime <- end_time-start_time

totaltime
